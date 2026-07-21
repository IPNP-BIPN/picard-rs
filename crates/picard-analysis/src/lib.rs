//! Port of Picard's `picard.analysis` tools.
//!
//! Ported from Picard 3.4.0, symbol by symbol, from the pinned clone in `picard/`.

pub mod quality_yield;

pub use quality_yield::{Options, QualityYieldMetrics, QualityYieldMetricsCollector};
