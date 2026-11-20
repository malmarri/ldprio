use std::fs::File;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

pub struct ReportData {
    pub version: String,
    pub out: String,
    pub command: String,
    pub runtime_secs: f64,
    pub weighting: String,
    pub r2: f64,
    pub window: usize,
    pub n_input: usize,
    pub n_snps: usize,
    pub n_fam: usize,
    pub n_ld: usize,
    pub pool: usize,
    pub prio_final: usize,
    pub greedy: usize,
    pub passes: Vec<(usize, usize, usize, usize)>, // (pass_idx, pairs, removed, keep)
    pub final_n: usize,
    pub converged: bool,
    pub samples: Vec<(String, usize, usize)>, // (label, avail, kept)
    pub tv_share: Option<(usize, usize)>,
    pub messages: Vec<String>,
}

const REPORT_CSS: &str = r#"
:root{--bg:#fff;--card:#f6f8fa;--text:#1f2328;--muted:#59636e;--line:#d1d9e0;
--track:#e6eaee;--accent:#0969da;--good:#1a7f37;--bad:#cf222e;--warn:#9a6700}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#0d1117;
--card:#151b23;--text:#e6edf3;--muted:#9198a1;--line:#3d444d;--track:#262c36;
--accent:#4493f8;--good:#3fb950;--bad:#f85149;--warn:#d29922}}
:root[data-theme="dark"]{--bg:#0d1117;--card:#151b23;--text:#e6edf3;--muted:#9198a1;
--line:#3d444d;--track:#262c36;--accent:#4493f8;--good:#3fb950;--bad:#f85149;--warn:#d29922}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);
font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:920px;margin:0 auto;padding:28px 16px 48px}
h1{font-size:24px;margin:0 0 4px}
h2{font-size:17px;margin:34px 0 4px}
.sub,.hint{color:var(--muted);font-size:13px;margin:0 0 12px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin-top:18px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.tile .k{color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.06em}
.tile .v{font-size:26px;font-weight:600;font-variant-numeric:tabular-nums}
.tile .s{color:var(--muted);font-size:13px}
.good{color:var(--good)}.bad{color:var(--bad)}
.row{display:grid;grid-template-columns:minmax(80px,150px) 1fr minmax(120px,auto);
gap:10px;align-items:center;padding:3px 0}
.row.wide{grid-template-columns:minmax(120px,250px) 1fr minmax(70px,auto)}
.lab{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:13px}
.track{height:14px;background:var(--track);border-radius:7px;overflow:hidden}
.fill{height:100%;background:var(--accent);border-radius:7px;min-width:2px}
.fill.zero{min-width:0}
.num{font-size:13px;text-align:right;color:var(--muted);font-variant-numeric:tabular-nums}
.hist{display:flex;gap:6px;align-items:flex-end;height:170px;margin:12px 0 4px}
.bin{flex:1;display:flex;flex-direction:column;justify-content:flex-end;align-items:center;height:100%}
.hb{width:100%;background:var(--accent);border-radius:4px 4px 0 0;min-height:2px}
.hb.zero{min-height:0}
.hn{font-size:12px;color:var(--muted);margin-bottom:2px;font-variant-numeric:tabular-nums}
.hl{font-size:11px;color:var(--muted);margin-top:4px}
table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:500}td.r,th.r{text-align:right;font-variant-numeric:tabular-nums}
dl{display:grid;grid-template-columns:max-content 1fr;gap:4px 18px;margin:8px 0;font-size:14px}
dt{color:var(--muted)}dd{margin:0;word-break:break-word}
code{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}
ul{padding-left:20px;margin:8px 0}li{margin:5px 0;font-size:14px}
li.w{color:var(--warn)}
footer{color:var(--muted);font-size:12px;margin-top:40px}
"#;

fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#x27;"),
            _ => out.push(c),
        }
    }
    out
}

pub fn fmt_int(x: usize) -> String {
    let s = x.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(b as char);
    }
    out
}

fn pct(a: usize, b: usize) -> String {
    if b == 0 {
        "n/a".to_string()
    } else {
        format!("{:.1}%", 100.0 * (a as f64) / (b as f64))
    }
}

/// UTC YYYY-MM-DD HH:MM using Howard Hinnant's civil-from-days algorithm (zero deps).
fn utc_timestamp() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    let hour = rem / 3600;
    let min = (rem % 3600) / 60;

    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{:04}-{:02}-{:02} {:02}:{:02} UTC", y, m, d, hour, min)
}

pub fn write_report(path: &str, r: &ReportData) -> Result<(), String> {
    let mut samples = r.samples.clone();
    samples.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let usable: Vec<&(String, usize, usize)> = samples.iter().filter(|t| t.1 > 0).collect();
    let mean_share = if usable.is_empty() {
        0.0
    } else {
        usable
            .iter()
            .map(|t| (t.2 as f64) / (t.1 as f64))
            .sum::<f64>()
            / (usable.len() as f64)
    };

    let tile = |k: &str, v: &str, sub: &str, cls: &str| {
        format!(
            "<div class=\"tile\"><div class=\"k\">{}</div><div class=\"v {}\">{}</div><div class=\"s\">{}</div></div>",
            k, cls, v, sub
        )
    };

    let ld_sub = if r.converged {
        format!("no pair above r&sup2; {} within {} variants", r.r2, r.window)
    } else {
        "raise --max-cleanup".to_string()
    };

    let tiles = [
        tile(
            "Final panel",
            &fmt_int(r.final_n),
            &format!("of {} autosomal SNPs", fmt_int(r.n_snps)),
            "",
        ),
        tile(
            "Priority samples",
            &fmt_int(samples.len()),
            &format!("{} weighting", esc(&r.weighting)),
            "",
        ),
        tile(
            "Average share kept",
            &format!("{:.0}%", 100.0 * mean_share),
            "of each priority sample's own SNPs",
            "",
        ),
        tile(
            "LD check",
            if r.converged { "Passed" } else { "Not verified" },
            &ld_sub,
            if r.converged { "good" } else { "bad" },
        ),
    ]
    .concat();

    let (per_sample_hint, per_sample) = if samples.len() <= 40 {
        let mut rows = String::new();
        for (iid, avail, kept) in &samples {
            let share = if *avail > 0 {
                (*kept as f64) / (*avail as f64)
            } else {
                0.0
            };
            let cls = if *kept == 0 { "fill zero" } else { "fill" };
            rows.push_str(&format!(
                "<div class=\"row\"><div class=\"lab\" title=\"{0}\">{0}</div><div class=\"track\"><div class=\"{1}\" style=\"width:{2:.1}%\"></div></div><div class=\"num\">{3} of {4} ({5})</div></div>",
                esc(iid),
                cls,
                100.0 * share,
                fmt_int(*kept),
                fmt_int(*avail),
                pct(*kept, *avail)
            ));
        }
        (
            "<p class=\"hint\">Bar = share of the sample's own SNPs (those it has calls at in the input panel) that survived pruning.</p>",
            rows,
        )
    } else {
        let mut bins = [0usize; 10];
        for &&(_, a, k) in &usable {
            let b = ((10.0 * (k as f64) / (a as f64)) as usize).min(9);
            bins[b] += 1;
        }
        let top = (*bins.iter().max().unwrap_or(&0)).max(1);
        let mut hist = String::new();
        for (i, &c) in bins.iter().enumerate() {
            let cls = if c == 0 { "hb zero" } else { "hb" };
            hist.push_str(&format!(
                "<div class=\"bin\"><div class=\"hn\">{}</div><div class=\"{}\" style=\"height:{:.1}%\"></div><div class=\"hl\">{}&ndash;{}%</div></div>",
                c,
                cls,
                100.0 * (c as f64) / (top as f64) * 0.78,
                10 * i,
                10 * i + 10
            ));
        }
        let mut worst = samples.clone();
        worst.sort_by(|a, b| a.2.cmp(&b.2).then_with(|| a.0.cmp(&b.0)));
        worst.truncate(15);
        let mut trs = String::new();
        for (i, a, k) in worst {
            trs.push_str(&format!(
                "<tr><td>{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td></tr>",
                esc(&i),
                fmt_int(a),
                fmt_int(k),
                pct(k, a)
            ));
        }
        (
            "<p class=\"hint\">Number of priority samples by share of their own SNPs (those they have calls at in the input panel) that survived pruning.</p>",
            format!(
                "<div class=\"hist\">{}</div><h3 style=\"font-size:14px;margin:22px 0 0\">15 samples with the fewest SNPs kept</h3><table><tr><th>Sample</th><th class=\"r\">SNPs with calls</th><th class=\"r\">SNPs kept</th><th class=\"r\">Share kept</th></tr>{}</table>",
                hist, trs
            ),
        )
    };

    let mut steps: Vec<(String, usize, String)> = Vec::new();
    steps.push(("Autosomal SNPs in input".to_string(), r.n_snps, String::new()));
    steps.push(("After first LD pass".to_string(), r.greedy, String::new()));
    for &(i, pairs, removed, keep) in &r.passes {
        steps.push((
            format!("Cleanup pass {}", i),
            keep,
            format!("{} residual pairs, removed {}", fmt_int(pairs), fmt_int(removed)),
        ));
    }
    steps.push(("Final panel".to_string(), r.final_n, String::new()));

    let mut funnel = String::new();
    for (lab, v, extra) in steps {
        funnel.push_str(&format!(
            "<div class=\"row wide\"><div class=\"lab\">{}</div><div class=\"track\"><div class=\"fill\" style=\"width:{:.1}%\"></div></div><div class=\"num\">{}</div></div>",
            esc(&lab),
            100.0 * (v as f64) / (r.n_snps as f64),
            fmt_int(v)
        ));
        if !extra.is_empty() {
            funnel.push_str(&format!(
                "<div class=\"hint\" style=\"margin:-2px 0 6px 0\">{}</div>",
                extra
            ));
        }
    }

    let mut settings: Vec<(&str, String)> = vec![
        ("Weighting", esc(&r.weighting)),
        ("LD samples", fmt_int(r.n_ld)),
        (
            "r&sup2; threshold / window",
            format!("{} / {} variants", r.r2, r.window),
        ),
        (
            "Priority pool",
            format!("{} SNPs ({} of panel)", fmt_int(r.pool), pct(r.pool, r.n_snps)),
        ),
        (
            "Pool SNPs retained",
            format!("{} ({})", fmt_int(r.prio_final), pct(r.prio_final, r.pool)),
        ),
        (
            "Non-autosomal SNPs dropped",
            fmt_int(r.n_input.saturating_sub(r.n_snps)),
        ),
        ("Samples in input", fmt_int(r.n_fam)),
    ];
    if let Some((n_tv, n_fin)) = r.tv_share {
        settings.push(("Transversions in output", pct(n_tv, n_fin)));
    }
    settings.push((
        "Run time",
        if r.runtime_secs < 10.0 {
            format!("{:.2} s", r.runtime_secs)
        } else {
            format!("{:.0} s", r.runtime_secs)
        },
    ));
    let mut dl = String::new();
    for (k, v) in settings {
        dl.push_str(&format!("<dt>{}</dt><dd>{}</dd>", k, v));
    }

    let msgs = if r.messages.is_empty() {
        "<li class=\"good\">No warnings.</li>".to_string()
    } else {
        r.messages
            .iter()
            .map(|m| {
                if let Some(rest) = m.strip_prefix("WARNING:") {
                    format!("<li class=\"w\">{}</li>", esc(rest.trim()))
                } else if let Some(rest) = m.strip_prefix("note:") {
                    format!("<li>{}</li>", esc(rest.trim()))
                } else {
                    format!("<li>{}</li>", esc(m))
                }
            })
            .collect::<String>()
    };

    let when = utc_timestamp();
    let doc = format!(
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">\n\
         <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n\
         <title>ldprio report</title><style>{css}</style></head><body><main>\n\
         <h1>ldprio report</h1>\n\
         <p class=\"sub\">v{ver} &middot; {when} &middot; <code>{out}</code></p>\n\
         <div class=\"tiles\">{tiles}</div>\n\
         <h2>What each priority sample kept</h2>\n\
         {per_sample_hint}{per_sample}\n\
         <h2>Panel size through the run</h2>\n\
         <p class=\"hint\">Bars are relative to the autosomal SNPs in the input.</p>\n\
         {funnel}\n\
         <h2>Run details</h2>\n\
         <dl>{dl}</dl>\n\
         <h2>Warnings and notes</h2>\n\
         <ul>{msgs}</ul>\n\
         <h2>Command</h2>\n\
         <p><code>{cmd}</code></p>\n\
         <footer>Generated by ldprio {ver}. Self-contained file: no scripts or external resources.</footer>\n\
         </main></body></html>\n",
        css = REPORT_CSS,
        ver = esc(&r.version),
        when = esc(&when),
        out = esc(&r.out),
        tiles = tiles,
        per_sample_hint = per_sample_hint,
        per_sample = per_sample,
        funnel = funnel,
        dl = dl,
        msgs = msgs,
        cmd = esc(&r.command),
    );

    let mut fh = File::create(path)
        .map_err(|e| format!("ERROR: failed to create {}: {}", path, e))?;
    fh.write_all(doc.as_bytes())
        .map_err(|e| format!("ERROR: failed to write {}: {}", path, e))?;
    Ok(())
}
