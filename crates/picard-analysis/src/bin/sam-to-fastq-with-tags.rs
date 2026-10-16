//! `SamToFastqWithTags` as a runnable binary: the covering array's port side.
//!
//! `SamToFastq` (see `sam-to-fastq`), plus one more FASTQ per `SEQUENCE_TAG_GROUP` whose sequence
//! is the group's tag values joined by the group's separator and whose qualities are the
//! `QUALITY_TAG_GROUP` values, or `~` throughout when none is given. The tag files sit in
//! `OUTPUT_DIR`, or beside `FASTQ`, prefixed with the read group under `OUTPUT_PER_RG`.

fn main() {
    picard_analysis::command_line::check("SamToFastqWithTags");
    picard_analysis::sam_to_fastq_full::main_with(true)
}
