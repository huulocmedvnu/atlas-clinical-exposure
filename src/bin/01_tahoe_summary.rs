//! Step 1. Reads the Tahoe-100M pseudobulk differential expression table shard by shard from Hugging Face and writes,
//! per shard, the number of DE genes of each contrast at five |log2FC| thresholds, the DE rows of a list of compounds
//! and the full log2FC profiles of a second list.
//! Usage: tahoe_summary <out dir> <first shard> <last shard, exclusive> <compounds with DE rows> <compounds with full profiles>
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::{ArrowWriter, ProjectionMask};
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

const BASE: &str = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/refs%2Fconvert%2Fparquet/pseudobulk_differential_expression/train/";
const CUTS: [f32; 5] = [0.25, 0.5, 0.75, 1.0, 1.5];
const COLS: [&str; 9] = ["gene_name", "baseMean", "log2FoldChange", "stat", "padj", "plate", "drug", "concentration", "Cell_Name_Vevo"];

#[derive(Default)]
struct Rows { drug: Vec<String>, conc: Vec<f32>, line: Vec<String>, plate: Vec<String>, gene: Vec<String>, lfc: Vec<f32>, stat: Vec<f32>, padj: Vec<f32> }

fn write(path: &str, b: RecordBatch) -> Result<()> {
    let tmp = format!("{path}.tmp");
    let props = WriterProperties::builder().set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?)).build();
    let mut w = ArrowWriter::try_new(std::fs::File::create(&tmp)?, b.schema(), Some(props))?;
    w.write(&b)?; w.close()?; std::fs::rename(tmp, path)?; Ok(())
}

fn rows_batch(r: Rows, with_stat: bool) -> Result<RecordBatch> {
    let mut f = vec![Field::new("drug", DataType::Utf8, false), Field::new("concentration", DataType::Float32, false),
        Field::new("Cell_Name_Vevo", DataType::Utf8, false), Field::new("plate", DataType::Utf8, false),
        Field::new("gene_name", DataType::Utf8, false), Field::new("log2FoldChange", DataType::Float32, true)];
    let mut c: Vec<ArrayRef> = vec![Arc::new(StringArray::from(r.drug)), Arc::new(Float32Array::from(r.conc)),
        Arc::new(StringArray::from(r.line)), Arc::new(StringArray::from(r.plate)), Arc::new(StringArray::from(r.gene)),
        Arc::new(Float32Array::from(r.lfc))];
    if with_stat { f.push(Field::new("stat", DataType::Float32, true)); c.push(Arc::new(Float32Array::from(r.stat))); }
    f.push(Field::new("padj", DataType::Float32, true)); c.push(Arc::new(Float32Array::from(r.padj)));
    Ok(RecordBatch::try_new(Arc::new(Schema::new(f)), c)?)
}

fn s<'a>(b: &'a RecordBatch, n: &str) -> Result<&'a StringArray> {
    let a = b.column_by_name(n).ok_or(anyhow!("no {n}"))?;
    a.as_any().downcast_ref::<StringArray>().ok_or(anyhow!("{n} not utf8 {:?}", a.data_type()))
}
fn fl<'a>(b: &'a RecordBatch, n: &str) -> Result<&'a Float32Array> {
    let a = b.column_by_name(n).ok_or(anyhow!("no {n}"))?;
    a.as_any().downcast_ref::<Float32Array>().ok_or(anyhow!("{n} not f32 {:?}", a.data_type()))
}

fn shard(i: usize, out: &str, keep: &HashSet<String>, full: &HashSet<String>, cli: &reqwest::blocking::Client) -> Result<usize> {
    let dc = format!("{out}/c{i:04}.parquet");
    if std::path::Path::new(&dc).exists() { return Ok(0); }
    let bytes = cli.get(format!("{BASE}{i:04}.parquet")).send()?.error_for_status()?.bytes()?;
    let bld = ParquetRecordBatchReaderBuilder::try_new(bytes)?;
    let sch = bld.parquet_schema();
    let idx: Vec<usize> = COLS.iter().map(|c| sch.columns().iter().position(|x| x.name() == *c).unwrap()).collect();
    let mask = ProjectionMask::leaves(sch, idx);
    let rd = bld.with_projection(mask).with_batch_size(65536).build()?;
    let mut cnt: HashMap<(String, u32, String, String), [i64; 6]> = HashMap::new();
    let (mut g, mut f) = (Rows::default(), Rows::default());
    let mut nrow = 0usize;
    for b in rd {
        let b = b?; nrow += b.num_rows();
        let (gene, bm, lfc, st, padj, plate, drug, conc, line) = (s(&b, "gene_name")?, fl(&b, "baseMean")?, fl(&b, "log2FoldChange")?,
            fl(&b, "stat")?, fl(&b, "padj")?, s(&b, "plate")?, s(&b, "drug")?, fl(&b, "concentration")?, s(&b, "Cell_Name_Vevo")?);
        for r in 0..b.num_rows() {
            let d = drug.value(r);
            let key = (d.to_string(), conc.value(r).to_bits(), line.value(r).to_string(), plate.value(r).to_string());
            let e = cnt.entry(key).or_insert([0; 6]);
            let p = if padj.is_null(r) { f32::NAN } else { padj.value(r) };
            let l = if lfc.is_null(r) { f32::NAN } else { lfc.value(r) };
            if !p.is_nan() { e[5] += 1; }
            let sig = p < 0.05; let al = l.abs();
            for (k, c) in CUTS.iter().enumerate() { if sig && al >= *c { e[k] += 1; } }
            let push = |rw: &mut Rows, stv: bool| {
                rw.drug.push(d.to_string()); rw.conc.push(conc.value(r)); rw.line.push(line.value(r).to_string());
                rw.plate.push(plate.value(r).to_string()); rw.gene.push(gene.value(r).to_string()); rw.lfc.push(l); rw.padj.push(p);
                if stv { rw.stat.push(if st.is_null(r) { f32::NAN } else { st.value(r) }); }
            };
            if sig && al >= 0.25 && keep.contains(d) { push(&mut g, false); }
            if full.contains(d) && !bm.is_null(r) && bm.value(r) >= 5.0 { push(&mut f, true); }
        }
    }
    write(&format!("{out}/g{i:04}.parquet"), rows_batch(g, false)?)?;
    write(&format!("{out}/f{i:04}.parquet"), rows_batch(f, true)?)?;
    let mut ks: Vec<_> = cnt.into_iter().collect(); ks.sort_by(|a, b| a.0.cmp(&b.0));
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(ks.iter().map(|x| x.0 .0.clone()).collect::<Vec<_>>())),
        Arc::new(Float32Array::from(ks.iter().map(|x| f32::from_bits(x.0 .1)).collect::<Vec<_>>())),
        Arc::new(StringArray::from(ks.iter().map(|x| x.0 .2.clone()).collect::<Vec<_>>())),
        Arc::new(StringArray::from(ks.iter().map(|x| x.0 .3.clone()).collect::<Vec<_>>()))];
    let mut fields = vec![Field::new("drug", DataType::Utf8, false), Field::new("concentration", DataType::Float32, false),
        Field::new("Cell_Name_Vevo", DataType::Utf8, false), Field::new("plate", DataType::Utf8, false)];
    for (k, c) in CUTS.iter().enumerate() {
        fields.push(Field::new(format!("n_de_{}", if *c == 1.0 { "1.0".to_string() } else if *c == 1.5 { "1.5".to_string() } else { c.to_string() }), DataType::Int64, false));
        cols.push(Arc::new(Int64Array::from(ks.iter().map(|x| x.1[k]).collect::<Vec<_>>())));
    }
    fields.push(Field::new("n_tested", DataType::Int64, false)); cols.push(Arc::new(Int64Array::from(ks.iter().map(|x| x.1[5]).collect::<Vec<_>>())));
    write(&dc, RecordBatch::try_new(Arc::new(Schema::new(fields)), cols)?)?;
    Ok(nrow)
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (out, lo, hi) = (a[1].clone(), a[2].parse::<usize>()?, a[3].parse::<usize>()?);
    let rd = |p: &str| -> Result<HashSet<String>> { Ok(std::fs::read_to_string(p)?.lines().filter(|l| !l.is_empty()).map(String::from).collect()) };
    let (keep, full) = (rd(&a[4])?, rd(&a[5])?);
    std::fs::create_dir_all(&out)?;
    let cli = reqwest::blocking::Client::builder().timeout(std::time::Duration::from_secs(300)).build()?;
    let t0 = std::time::Instant::now();
    (lo..hi).into_par_iter().for_each(|i| {
        let mut last = String::new();
        for att in 0..6 {
            match shard(i, &out, &keep, &full, &cli) {
                Ok(n) => { eprintln!("{i} {n} {:.0}s", t0.elapsed().as_secs_f64()); return; }
                Err(e) => { last = e.to_string(); std::thread::sleep(std::time::Duration::from_secs(5 + 10 * att)); }
            }
        }
        eprintln!("{i} FAIL {last}");
    });
    eprintln!("done {:.0}s", t0.elapsed().as_secs_f64());
    Ok(())
}
