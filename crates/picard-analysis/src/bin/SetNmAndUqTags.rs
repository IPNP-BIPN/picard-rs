//! `SetNmAndUqTags` as a runnable binary. Picard keeps it as a `@Deprecated` subclass of
//! `SetNmMdAndUqTags` that overrides nothing, so this binary is that tool's, under its own name.

#[path = "SetNmMdAndUqTags.rs"]
mod set_nm_md_and_uq_tags;

/// The same declarations as the tool it subclasses, checked here too so that this binary, like
/// every other, refuses before it does anything (the included `main` checks again, harmlessly).
fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SetNmMdAndUqTags");
    set_nm_md_and_uq_tags::main()
}
