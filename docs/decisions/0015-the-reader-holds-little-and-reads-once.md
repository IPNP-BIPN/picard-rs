# 0015. The reader holds little and reads each file once

**Status:** accepted
**Date:** 2026-10-31

## Why this exists

Almost every tool sat at about 190 MB of peak RSS on the 2,000,000-read benchmark BAM, whatever it
did with the records. That memory was the lookahead of the decoding thread, not the tool's work:
it grew with the rayon pool, from 28 MB at one thread to 345 MB at sixteen. And each BAM was read
twice, once to check its blocks and once to inflate them. This records what replaced that, and what
was tried and dropped. Nothing here moves a byte: the outputs were compared with those of the
previous build, tool by tool.

## What changed

* **The lookahead is a fixed four batches of 1 MiB** rather than 32 batches of 256 KiB per round on
  16 cores. The tool takes the records on one thread, so more only holds more.
* **Three threads of the process's own decode the batches**, not rayon's pool. A first version kept
  rayon and shrank its rounds: RSS fell, but total CPU across the tools rose 7%, because the pool
  woke more often. These threads sleep on a channel, are shared by all of a tool's inputs, and at
  1 MiB wake a quarter as often as at 256 KiB.
* **A stream a tool takes record by record has a quarter of the cores inflating it**, not half.
  With two inputs, MergeSamFiles held 16 inflating workers, and its system time doubled once rayon's
  decoding no longer hid them. SortSam's parallel read keeps the half, because its read is limited by
  the stream.
* **Each input read at once makes the batches bigger**, so the threads wake no more often with two
  inputs than with one, and **decoding starts at the first record**, so a tool that opens an input
  for its header and opens it again does not decode ahead for nothing.
* **RevertSam reverts on the reading thread**: the work per record was too small to pay for waking
  rayon's pool, and that alone took 6% off its CPU.
* **Reading no longer zeroes 256 KiB per block read**, and a batch's record vector is sized from the
  lengths its records declare.
* **The block check reads only the block headers**, by position, so the blocks' bodies are read from
  disk once. A pipe, which has no positions, is still walked.
* **SortSam's merge reuses each run's record buffer**, and a run spends 5 bytes per record on what
  it adds rather than 13. An in-place sift replaces the pop and push under the coordinate order only:
  the queryname comparator is not transitive (decision 0014), so the walk of the heap can change its
  answer. A run also leaves out the part of each read name the record before it has, which is 12% of
  the temporary disk on Illumina names and 4% on `read123`, for a compare and a copy per record.
* **CollectAlignmentSummaryMetrics counts its small ids in tables** added to the tree histograms at
  the end, where a bin's value is the same sum of ones. **CollectQualityYieldMetrics** no longer
  allocates per read. The read-name maps of CollectWgsMetrics' slow collector, SamToFastq and
  EstimateLibraryComplexity use foldhash; none is iterated into an output.

## What was tried and dropped

LZ4 on SortSam's runs took about a third off the temporary disk (three runs of 72 MB against an
estimated 113 MB each) for 7% more CPU and 13% more peak RSS over the sort. No tool may pay more CPU
for this work, so the runs stay uncompressed.

mimalloc's purge settings (`MIMALLOC_PURGE_DELAY`, `MIMALLOC_PURGE_DECOMMITS`) moved neither CPU nor
RSS beyond noise. The samples in `swtch_pri`, rayon's yield, look like a large waste in a profile, but
the CPU time of a run did not change with the pool's size: they are a thread's idle time more than a
power cost, and were not the target they appear to be.

## The measurement

See the table in the pull request that introduced this, taken on one machine, an Apple M4 Max running
other work at the same time. Wall time there is noise, which is why instructions retired, CPU time
and RSS are what was compared, with runs interleaved. These are not x86-64 figures: the throughput job in CI is.
