//! The port's side: gkl-deflate's level 5 on the same blocks, one thread, the same digest.
use sha2::{Digest, Sha256};
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let data = std::fs::read(args.next().expect("payload")).expect("payload");
    let reps: usize = args.next().expect("reps").parse().expect("reps");
    // The port's BGZF writer with the GKL deflater, serial, as htsjdk's stream writes it.
    let mut bgzf = Vec::new();
    for rep in 0..=reps {
        let start = Instant::now();
        let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), 5, htsjdk_bgzf::Deflater::Gkl)
            .parallel(false);
        std::io::Write::write_all(&mut writer, &data).expect("in memory");
        bgzf = writer.into_inner().expect("in memory");
        let seconds = start.elapsed().as_secs_f64();
        if rep > 0 {
            println!("port_bgzf_rate {:.1} MB/s", data.len() as f64 / 1e6 / seconds);
        }
    }
    let hex: String = Sha256::digest(&bgzf).iter().map(|b| format!("{b:02x}")).collect();
    println!("port_bgzf sha256 {hex}");
    for rep in 0..=reps {
        let mut digest = Sha256::new();
        let mut compressed = 0usize;
        let start = Instant::now();
        for block in data.chunks(65498) {
            let out = gkl_deflate::deflate_gkl(block, 5);
            digest.update(&out);
            compressed += out.len();
        }
        let seconds = start.elapsed().as_secs_f64();
        let hex: String = digest.finalize().iter().map(|b| format!("{b:02x}")).collect();
        if rep > 0 {
            println!(
                "port {:.1} MB/s ratio {:.4} sha256 {hex}",
                data.len() as f64 / 1e6 / seconds,
                compressed as f64 / data.len() as f64
            );
        }
    }
}
