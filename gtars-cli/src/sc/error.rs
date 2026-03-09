use std::fmt;
use std::io::Write;

use serde::Serialize;

/// Exit code categories per the CLI design plan.
pub const EXIT_INPUT: i32 = 1; // File not found, invalid matrix format
pub const EXIT_PARAM: i32 = 2; // Invalid config, incompatible flags
pub const EXIT_RUNTIME: i32 = 3; // All cells filtered out, PCA failed

/// Categorized CLI error with structured output support.
#[derive(Debug)]
pub struct ScError {
    pub message: String,
    pub category: ErrorCategory,
    pub details: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Copy)]
pub enum ErrorCategory {
    Input,
    Parameter,
    Runtime,
}

impl ErrorCategory {
    pub fn exit_code(self) -> i32 {
        match self {
            ErrorCategory::Input => EXIT_INPUT,
            ErrorCategory::Parameter => EXIT_PARAM,
            ErrorCategory::Runtime => EXIT_RUNTIME,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            ErrorCategory::Input => "input_error",
            ErrorCategory::Parameter => "parameter_error",
            ErrorCategory::Runtime => "runtime_error",
        }
    }
}

impl fmt::Display for ScError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for ScError {}

impl ScError {
    pub fn input(msg: impl Into<String>) -> Self {
        Self {
            message: msg.into(),
            category: ErrorCategory::Input,
            details: None,
        }
    }

    pub fn param(msg: impl Into<String>) -> Self {
        Self {
            message: msg.into(),
            category: ErrorCategory::Parameter,
            details: None,
        }
    }

    pub fn runtime(msg: impl Into<String>) -> Self {
        Self {
            message: msg.into(),
            category: ErrorCategory::Runtime,
            details: None,
        }
    }

    pub fn with_details(mut self, details: serde_json::Value) -> Self {
        self.details = Some(details);
        self
    }
}

#[derive(Serialize)]
struct JsonError {
    error: String,
    category: String,
    exit_code: i32,
    #[serde(skip_serializing_if = "Option::is_none")]
    details: Option<serde_json::Value>,
}

/// Write an error to stderr, using JSON if the format is JSON-like.
pub fn report_error(err: &anyhow::Error, json_mode: bool) {
    // Try to extract our structured error
    if let Some(sc_err) = err.downcast_ref::<ScError>() {
        if json_mode {
            let json_err = JsonError {
                error: sc_err.message.clone(),
                category: sc_err.category.as_str().to_string(),
                exit_code: sc_err.category.exit_code(),
                details: sc_err.details.clone(),
            };
            if let Ok(json) = serde_json::to_string(&json_err) {
                let _ = writeln!(std::io::stderr(), "{}", json);
            }
        } else {
            let _ = writeln!(std::io::stderr(), "Error: {}", sc_err.message);
        }
        std::process::exit(sc_err.category.exit_code());
    }

    // Fall back: categorize anyhow errors by message heuristics
    let msg = format!("{:#}", err);
    let category = categorize_error(&msg);

    if json_mode {
        let json_err = JsonError {
            error: msg,
            category: category.as_str().to_string(),
            exit_code: category.exit_code(),
            details: None,
        };
        if let Ok(json) = serde_json::to_string(&json_err) {
            let _ = writeln!(std::io::stderr(), "{}", json);
        }
    } else {
        let _ = writeln!(std::io::stderr(), "Error: {:#}", err);
    }
    std::process::exit(category.exit_code());
}

/// Best-effort categorization of anyhow errors by message content.
fn categorize_error(msg: &str) -> ErrorCategory {
    let lower = msg.to_lowercase();
    if lower.contains("not found")
        || lower.contains("no such file")
        || lower.contains("reading")
        || lower.contains("opening")
        || lower.contains("none of")
    {
        ErrorCategory::Input
    } else if lower.contains("must be")
        || lower.contains("invalid")
        || lower.contains("parsing")
        || lower.contains("config")
    {
        ErrorCategory::Parameter
    } else {
        ErrorCategory::Runtime
    }
}
