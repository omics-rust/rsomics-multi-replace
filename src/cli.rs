use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_multi_replace::run;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-multi-replace", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Composition table TSV (samples × components, corner cell ignored); reads stdin when "-" or omitted.
    #[arg(long, default_value = "-")]
    input: PathBuf,

    /// Replacement value for zeros; the default (1/D)^2 (D = number of components) matches skbio.
    #[arg(long)]
    delta: Option<f64>,

    /// Parse and write comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let delim = if self.csv { ',' } else { '\t' };

        let reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let mut out: Box<dyn Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        run(reader, &mut out, delim, self.delta)?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Multiplicative zero-replacement of a compositional matrix.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.stats.composition.multi_replace",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1023/A:1023866030544"),
    }),
    usage_lines: &["--input table.tsv [--delta F] [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("path"),
                required: false,
                default: Some("-"),
                description: "Composition table (- for stdin).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "delta",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("(1/D)^2"),
                description: "Value substituted for each zero.",
                why_default: Some("(1/D)^2 is skbio's default"),
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse and write comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Replace zeros with the default delta",
            command: "rsomics-multi-replace --input table.tsv -o replaced.tsv",
        },
        Example {
            description: "Use a custom delta before a log-ratio transform",
            command: "rsomics-multi-replace --input table.tsv --delta 1e-6 | rsomics-clr",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
