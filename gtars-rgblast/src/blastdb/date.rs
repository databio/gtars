//! The database date: makeblastdb's text format in `.nin`/`.nal`, ISO 8601 in `.njs`.

use crate::error::Error;

const MON: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

/// One instant (UTC, minute precision), formatted two ways.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastDate {
    /// Seconds since the Unix epoch, truncated to the minute.
    secs: i64,
}

fn civil_from_days(days: i64) -> (i64, i64, i64) {
    // H. Hinnant's algorithm
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + if m <= 2 { 1 } else { 0 };
    (y, m, d)
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = if m > 2 { m - 3 } else { m + 9 };
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

impl BlastDate {
    pub fn from_epoch(secs: i64) -> Self {
        BlastDate {
            secs: secs - secs.rem_euclid(60),
        }
    }

    pub fn now() -> Self {
        let secs = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0);
        Self::from_epoch(secs)
    }

    /// Parse makeblastdb's date text, e.g. `Nov 8, 2026  12:51 PM`.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let bad = || {
            Error::InvalidOptions(format!(
                "date '{text}' is not in makeblastdb's format (e.g. 'Nov 8, 2026  12:51 PM')"
            ))
        };
        let f: Vec<&str> = text.split_whitespace().collect();
        if f.len() != 5 {
            return Err(bad());
        }
        let m = MON.iter().position(|&x| x == f[0]).ok_or_else(bad)? as i64 + 1;
        let d: i64 = f[1]
            .strip_suffix(',')
            .ok_or_else(bad)?
            .parse()
            .map_err(|_| bad())?;
        let y: i64 = f[2].parse().map_err(|_| bad())?;
        let (h, mi) = f[3].split_once(':').ok_or_else(bad)?;
        let (h, mi): (i64, i64) = (
            h.parse().map_err(|_| bad())?,
            mi.parse().map_err(|_| bad())?,
        );
        if !(1..=31).contains(&d) || !(1..=12).contains(&h) || !(0..60).contains(&mi) {
            return Err(bad());
        }
        let h24 = match f[4] {
            "AM" => h % 12,
            "PM" => h % 12 + 12,
            _ => return Err(bad()),
        };
        Ok(BlastDate {
            secs: days_from_civil(y, m, d) * 86400 + h24 * 3600 + mi * 60,
        })
    }

    /// `explicit` (makeblastdb text) if given, else `SOURCE_DATE_EPOCH` if set, else now.
    pub fn resolve(explicit: Option<&str>) -> Result<Self, Error> {
        Self::resolve_with(explicit, std::env::var("SOURCE_DATE_EPOCH").ok())
    }

    fn resolve_with(explicit: Option<&str>, epoch: Option<String>) -> Result<Self, Error> {
        if let Some(t) = explicit {
            return Self::parse(t);
        }
        if let Some(e) = epoch {
            let s: i64 = e.trim().parse().map_err(|_| {
                Error::InvalidOptions(format!(
                    "SOURCE_DATE_EPOCH '{e}' is not a number of seconds"
                ))
            })?;
            return Ok(Self::from_epoch(s));
        }
        Ok(Self::now())
    }

    fn parts(&self) -> (i64, i64, i64, i64, i64) {
        let (y, m, d) = civil_from_days(self.secs.div_euclid(86400));
        let sod = self.secs.rem_euclid(86400);
        (y, m, d, sod / 3600, (sod % 3600) / 60)
    }

    /// makeblastdb's format, e.g. `Nov 8, 2026  12:51 PM`.
    pub fn blast(&self) -> String {
        let (y, m, d, h, mi) = self.parts();
        let (h12, ap) = (
            if h % 12 == 0 { 12 } else { h % 12 },
            if h < 12 { "AM" } else { "PM" },
        );
        format!(
            "{} {}, {}  {}:{:02} {}",
            MON[(m - 1) as usize],
            d,
            y,
            h12,
            mi,
            ap
        )
    }

    /// ISO 8601 with seconds, as `.njs` `last-updated`: `2026-11-08T12:51:00`.
    pub fn iso(&self) -> String {
        let (y, m, d, h, mi) = self.parts();
        format!("{y:04}-{m:02}-{d:02}T{h:02}:{mi:02}:00")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blast_format() {
        // 2026-11-08 12:51:30 UTC
        let t = days_from_civil(2026, 11, 8) * 86400 + 12 * 3600 + 51 * 60 + 30;
        let d = BlastDate::from_epoch(t);
        assert_eq!(d.blast(), "Nov 8, 2026  12:51 PM");
        assert_eq!(
            BlastDate::from_epoch(days_from_civil(2026, 1, 15) * 86400 + 5 * 60).blast(),
            "Jan 15, 2026  12:05 AM"
        );
        assert_eq!(BlastDate::parse("Nov 8, 2026  12:51 PM").unwrap(), d);
        assert_eq!(
            BlastDate::parse("Jan 15, 2026  12:05 AM").unwrap().blast(),
            "Jan 15, 2026  12:05 AM"
        );
        assert!(BlastDate::parse("2026-11-08").is_err());
    }

    #[test]
    fn iso_format() {
        let d = BlastDate::parse("Nov 8, 2026  2:18 PM").unwrap();
        assert_eq!(d.iso(), "2026-11-08T14:18:00");
        assert_eq!(
            BlastDate::parse("Feb 29, 2024  11:59 PM").unwrap().iso(),
            "2024-02-29T23:59:00"
        );
    }

    #[test]
    fn source_date_epoch() {
        let d = BlastDate::resolve_with(None, Some("1700000000".into())).unwrap();
        assert_eq!(d.iso(), "2023-11-14T22:13:00");
        assert_eq!(d.blast(), "Nov 14, 2023  10:13 PM");
        // an explicit date wins over the environment
        let e = BlastDate::resolve_with(Some("Nov 8, 2026  12:51 PM"), Some("1700000000".into()))
            .unwrap();
        assert_eq!(e.iso(), "2026-11-08T12:51:00");
        assert!(BlastDate::resolve_with(None, Some("soon".into())).is_err());
    }
}
