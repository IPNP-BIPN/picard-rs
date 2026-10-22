# 0014. SortSam reads and sorts on every core, and the deflater is the floor

**Status:** accepted
**Date:** 2026-10-22

## Why this exists

Someone asked how picard-rs compares with sambamba. The first measurement answered for SortSam:
2.06 s against sambamba's 0.60 s for a coordinate sort, because SortSam decoded, validated and
re-encoded every record on one thread and sorted on one thread, as
[decision 0006](0006-per-record-transforms-parallelize-without-touching-the-bytes.md) had left it.
This records what moved it, and why none of it moves a byte.

## What changed

For a clean BAM input, in coordinate or queryname order:

* **The read is parallel.** The inflated stream is cut into batches of whole records on a
  thread of its own, and each batch is decoded, checked and encoded on the rayon pool. The batches
  are added in file order, and a spill falls on exactly the record where adding them one by one
  spills. Anything the serial read answers for (a damaged block, a record that does not decode or
  encode, a record that does not pass validation in silence) sends the file to the serial read,
  which reproduces the message. Nothing is written before every record is read, so falling back is
  only a second read.
* **The BGZF structure check runs beside the read**, not before it, for the same reason.
* **The coordinate sort is parallel.** `SAMRecordCoordinateComparator` is a total preorder, so
  every stable sort orders a batch as the serial one does; rayon's `par_sort_by` is stable. The
  queryname comparator is not transitive (a record flagged both first and second of pair is less
  and greater than its mate), so its batches keep the serial sort.
* **A spilled batch is sorted and written on a thread of its own** while the next one fills, through
  a 1 MiB buffer instead of 8 KiB, and the last batch is merged from RAM, as the last run it would
  have been. At most two batches are in RAM at once.
* **The read and the write each run on every core.** They never overlap, so the inflating workers'
  half of the cores is free once the input is read.

## The measurement

linux/arm64 in Docker, 16 cores, the 2,000,000-read benchmark BAM grouped by name, the defaults
(`MAX_RECORDS_IN_RAM=500000`, three spills), runs interleaved, median of nine:

| | median | peak RSS |
|---|---:|---:|
| SortSam before | 2.06 s (best of 3) | 242 MB |
| SortSam after, with the deflater of htsjdk-rs#273 | 0.67 s | 427 MB |
| sambamba 1.0.1 `sort -t 16 -l 5` | 0.65 s | 588 MB |

The output is the same 181 MB as before, byte for byte, with and without spills. Against the
pinned htsjdk-rs, without the faster deflater, the same build ran at 0.68 to 0.70 s.

## Where the floor is

sambamba spends more CPU than picard-rs (8.4 core-seconds against 7.7) and wins on how busy it
keeps the cores. About 60% of picard-rs's CPU is the deflater, which must write GKL's bytes:
the port of Intel's zlib runs at about 100 MB/s per core where libdeflate runs at 177 MB/s.
[IPNP-BIPN/htsjdk-rs#273](https://github.com/IPNP-BIPN/htsjdk-rs/pull/273) took it from 88 to
101 MB/s. The rest of the gap is that deflater, and it does not close without a faster one that
writes the same bytes.
