//! `JAVASCRIPT_FILE`: htsjdk's `AbstractJavascriptFilter`, which Picard runs on Nashorn.
//!
//! Ported from `htsjdk.samtools.filter.AbstractJavascriptFilter`,
//! `htsjdk.samtools.filter.JavascriptSamRecordFilter` and
//! `htsjdk.variant.variantcontext.filter.JavascriptVariantFilter` at the versions Picard 3.4.0
//! ships. The script is compiled once and evaluated once per record, in one global scope that
//! holds `header` and the record (`record` for a read, `variant` for a variant), so what a script
//! declares survives from one record to the next, as it does with Nashorn's bindings.
//!
//! The record is given to the script as an object whose methods are the Java getters a filter
//! script calls (`record.getMappingQuality()`, `variant.getStart()`), each answering what the Java
//! one answers for that record. The script's last value decides, as `accept` reads it: a boolean is
//! itself, a number is accepted when its `intValue()` is one, and anything else is a rejection.

use boa_engine::object::ObjectInitializer;
use boa_engine::property::Attribute;
use boa_engine::{js_string, Context, JsObject, JsValue, NativeFunction, Script, Source};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::TagValue;
use htsjdk_vcf::variant::{Value, VariantContext};

/// A compiled filter script and the engine it runs in.
pub struct JavascriptFilter {
    context: Context,
    script: Script,
}

/// A getter's answer, before it becomes a JavaScript value.
#[derive(Clone)]
enum Answer {
    Null,
    Bool(bool),
    Int(i64),
    Double(f64),
    Str(String),
}

impl Answer {
    fn value(&self) -> JsValue {
        match self {
            Answer::Null => JsValue::null(),
            Answer::Bool(b) => JsValue::from(*b),
            Answer::Int(i) => JsValue::from(*i as f64),
            Answer::Double(d) => JsValue::from(*d),
            Answer::Str(s) => JsValue::from(js_string!(s.as_str())),
        }
    }
}

impl JavascriptFilter {
    /// `new AbstractJavascriptFilter(scriptFile, header)`: the script is parsed here, and a syntax
    /// error is the reference's `RuntimeScriptException("Script error in input", ...)`.
    pub fn new(script: &str) -> Result<JavascriptFilter, String> {
        let mut context = Context::default();
        let script =
            Script::parse(Source::from_bytes(script), None, &mut context).map_err(|_| {
                "htsjdk.samtools.util.RuntimeScriptException: Script error in input".to_string()
            })?;
        // `bindings.put(DEFAULT_HEADER_KEY, header)`: the header is there, as an object.
        let header = ObjectInitializer::new(&mut context).build();
        context
            .register_global_property(js_string!("header"), header, Attribute::all())
            .map_err(|e| e.to_string())?;
        Ok(JavascriptFilter { context, script })
    }

    /// `accept(record)` with the record bound under `key` as an object of `getters`.
    fn accept(
        &mut self,
        key: &str,
        getters: Vec<(&'static str, Answer)>,
        maps: Vec<(&'static str, Vec<(String, Answer)>)>,
    ) -> Result<bool, String> {
        let object = {
            let mut init = ObjectInitializer::new(&mut self.context);
            for (name, answer) in getters {
                let captured = answer.value();
                init.function(
                    NativeFunction::from_copy_closure_with_captures(
                        |_, _, value, _| Ok(value.clone()),
                        captured,
                    ),
                    js_string!(name),
                    0,
                );
            }
            init.build()
        };
        // A getter that takes a key: the answers by key, and `null` for one the record lacks.
        for (name, entries) in maps {
            let table = ObjectInitializer::new(&mut self.context).build();
            for (k, answer) in entries {
                table
                    .set(
                        js_string!(k.as_str()),
                        answer.value(),
                        false,
                        &mut self.context,
                    )
                    .map_err(|e| e.to_string())?;
            }
            let getter = NativeFunction::from_copy_closure_with_captures(
                |_, args, table: &JsObject, context| {
                    let key = args
                        .first()
                        .cloned()
                        .unwrap_or_default()
                        .to_string(context)?;
                    if table.has_own_property(key.clone(), context)? {
                        table.get(key, context)
                    } else {
                        Ok(JsValue::null())
                    }
                },
                table,
            );
            let function =
                boa_engine::object::FunctionObjectBuilder::new(self.context.realm(), getter)
                    .name(js_string!(name))
                    .length(1)
                    .build();
            object
                .set(js_string!(name), function, false, &mut self.context)
                .map_err(|e| e.to_string())?;
        }
        self.context
            .global_object()
            .set(js_string!(key), object, false, &mut self.context)
            .map_err(|e| e.to_string())?;
        let result = self.script.evaluate(&mut self.context).map_err(|e| {
            format!("java.lang.RuntimeException: javax.script.ScriptException: {e}")
        })?;
        // `Boolean.TRUE.equals(result)`, or `((Number) result).intValue() == 1`.
        Ok(if let Some(b) = result.as_boolean() {
            b
        } else if let Some(d) = result.as_number() {
            d.is_finite() && d.trunc() == 1.0
        } else {
            false
        })
    }

    /// `JavascriptSamRecordFilter.filterOut(record) == !accept(record)`, with the read as `record`.
    pub fn accept_read(&mut self, record: &BamRecord, header: &SamHeader) -> Result<bool, String> {
        let flag = |bit: u16| Answer::Bool(record.flags & bit != 0);
        let contig = |index: i32| match header.sequences.get(index.max(0) as usize) {
            Some(sequence) if index >= 0 => Answer::Str(sequence.name.clone()),
            _ => Answer::Str("*".to_string()),
        };
        let unmapped = record.flags & 0x4 != 0;
        let reference_length: i64 = record
            .cigar
            .elements
            .iter()
            .filter(|e| e.op.consumes_reference_bases())
            .map(|e| e.length as i64)
            .sum();
        let tags: Vec<(String, Answer)> = record
            .tags
            .iter()
            .map(|(tag, value)| {
                let answer = match value {
                    TagValue::Str(s) => Answer::Str(s.clone()),
                    TagValue::Int(i) => Answer::Int(*i),
                    TagValue::Float(f) => Answer::Double(*f as f64),
                    TagValue::Char(c) => Answer::Str((*c as char).to_string()),
                    _ => Answer::Null,
                };
                (tag.to_string(), answer)
            })
            .collect();
        let bases = String::from_utf8_lossy(&record.read_bases).into_owned();
        let getters = vec![
            ("getReadName", Answer::Str(record.read_name.clone())),
            ("getReadLength", Answer::Int(record.read_bases.len() as i64)),
            (
                "getReadString",
                Answer::Str(if bases.is_empty() {
                    "*".to_string()
                } else {
                    bases
                }),
            ),
            ("getFlags", Answer::Int(record.flags as i64)),
            (
                "getMappingQuality",
                Answer::Int(record.mapping_quality as i64),
            ),
            (
                "getAlignmentStart",
                Answer::Int(record.alignment_start as i64),
            ),
            (
                "getAlignmentEnd",
                Answer::Int(if unmapped {
                    0
                } else {
                    record.alignment_start as i64 + reference_length - 1
                }),
            ),
            ("getReferenceName", contig(record.reference_index)),
            (
                "getContig",
                if unmapped {
                    Answer::Null
                } else {
                    contig(record.reference_index)
                },
            ),
            (
                "getReferenceIndex",
                Answer::Int(record.reference_index as i64),
            ),
            ("getMateReferenceName", contig(record.mate_reference_index)),
            (
                "getMateAlignmentStart",
                Answer::Int(record.mate_alignment_start as i64),
            ),
            (
                "getInferredInsertSize",
                Answer::Int(record.inferred_insert_size as i64),
            ),
            ("getCigarString", Answer::Str(record.cigar.to_text())),
            ("getReadPairedFlag", flag(0x1)),
            ("getProperPairFlag", flag(0x2)),
            ("getReadUnmappedFlag", flag(0x4)),
            ("getMateUnmappedFlag", flag(0x8)),
            ("getReadNegativeStrandFlag", flag(0x10)),
            ("getMateNegativeStrandFlag", flag(0x20)),
            ("getFirstOfPairFlag", flag(0x40)),
            ("getSecondOfPairFlag", flag(0x80)),
            ("getNotPrimaryAlignmentFlag", flag(0x100)),
            ("isSecondaryAlignment", flag(0x100)),
            ("getReadFailsVendorQualityCheckFlag", flag(0x200)),
            ("getDuplicateReadFlag", flag(0x400)),
            ("getSupplementaryAlignmentFlag", flag(0x800)),
        ];
        self.accept("record", getters, vec![("getAttribute", tags)])
    }

    /// `JavascriptVariantFilter.accept(variant)`, with the variant as `variant`.
    pub fn accept_variant(&mut self, variant: &VariantContext) -> Result<bool, String> {
        let attribute = |value: &Value| -> Answer {
            match value {
                Value::Int(i) => Answer::Int(*i),
                Value::Double(d) => Answer::Double(*d),
                Value::Bool(b) => Answer::Bool(*b),
                Value::Str(s) => Answer::Str(s.clone()),
                Value::Missing => Answer::Null,
                Value::List(_) => value.format().map_or(Answer::Null, Answer::Str),
            }
        };
        let attributes: Vec<(String, Answer)> = variant
            .attributes
            .iter()
            .map(|(k, v)| (k.clone(), attribute(v)))
            .collect();
        let has: Vec<(String, Answer)> = variant
            .attributes
            .iter()
            .map(|(k, _)| (k.clone(), Answer::Bool(true)))
            .collect();
        let snp =
            variant.alleles.iter().all(|a| a.base_string().len() == 1) && variant.alleles.len() > 1;
        let indel = variant.alleles.len() > 1
            && variant.alleles.iter().skip(1).all(|a| {
                let b = a.base_string();
                !b.is_empty() && b.chars().all(|c| "ACGTNacgtn".contains(c))
            })
            && variant
                .alleles
                .iter()
                .any(|a| a.base_string().len() != variant.alleles[0].base_string().len());
        let getters = vec![
            ("getContig", Answer::Str(variant.contig.clone())),
            ("getChr", Answer::Str(variant.contig.clone())),
            ("getStart", Answer::Int(variant.start)),
            ("getEnd", Answer::Int(variant.stop)),
            ("getID", Answer::Str(variant.id.clone())),
            (
                "getPhredScaledQual",
                Answer::Double(variant.phred_scaled_qual()),
            ),
            ("isFiltered", Answer::Bool(variant.is_filtered())),
            ("isVariant", Answer::Bool(variant.is_variant())),
            ("isSNP", Answer::Bool(snp)),
            ("isIndel", Answer::Bool(indel)),
            ("getNAlleles", Answer::Int(variant.alleles.len() as i64)),
            ("getNSamples", Answer::Int(variant.genotypes.len() as i64)),
        ];
        self.accept(
            "variant",
            getters,
            vec![("getAttribute", attributes), ("hasAttribute", has)],
        )
    }
}
