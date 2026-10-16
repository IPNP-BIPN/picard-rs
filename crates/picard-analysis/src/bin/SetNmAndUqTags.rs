//! `SetNmAndUqTags` as a runnable binary. Picard keeps it as a `@Deprecated` subclass of
//! `SetNmMdAndUqTags` that overrides nothing, so this binary is that tool's, under its own name.

#[path = "SetNmMdAndUqTags.rs"]
mod set_nm_md_and_uq_tags;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    set_nm_md_and_uq_tags::main()
}
