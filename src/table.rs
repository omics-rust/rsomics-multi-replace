use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A samples × components compositional matrix with row (sample) and column
/// (component) IDs, stored row-major.
pub struct Table {
    pub samples: Vec<String>,
    pub components: Vec<String>,
    pub data: Vec<f64>,
}

impl Table {
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// Header line is component IDs (corner cell ignored); each following line
    /// is a sample ID followed by its component values.
    pub fn parse(reader: impl BufRead, delim: char) -> Result<Self> {
        let mut lines = reader.lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty table".into()))?
            .map_err(RsomicsError::Io)?;
        let components: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if components.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "table header has no component columns".into(),
            ));
        }
        let m = components.len();

        let mut samples = Vec::new();
        let mut data = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() {
                continue;
            }
            let mut cells = line.split(delim);
            let id = cells
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("table row has no cells".into()))?;
            samples.push(id.trim().to_string());
            let before = data.len();
            for cell in cells {
                let v: f64 = cell.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric table value: '{cell}'"))
                })?;
                data.push(v);
            }
            if data.len() - before != m {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample '{}' has {} values, expected {m}",
                    samples.last().unwrap(),
                    data.len() - before
                )));
            }
        }
        if samples.is_empty() {
            return Err(RsomicsError::InvalidInput("table has no samples".into()));
        }
        Ok(Self {
            samples,
            components,
            data,
        })
    }
}
