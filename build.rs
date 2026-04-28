fn main() {
    // No C++ compilation needed — all functionality is pure Rust.
    println!("cargo:rerun-if-changed=data/hyperloglogplus-bias.h");
    generate_hll_bias_tables();
}

fn generate_hll_bias_tables() {
    let header = std::fs::read_to_string("data/hyperloglogplus-bias.h")
        .expect("read hyperloglogplus-bias.h");
    let threshold =
        extract_braced_after(&header, "threshold[] =").expect("extract threshold table");
    let raw =
        extract_braced_after(&header, "rawEstimateData[] =").expect("extract raw estimate table");
    let bias = extract_braced_after(&header, "biasData[] =").expect("extract bias table");

    let mut out = String::new();
    out.push_str("// Generated from data/hyperloglogplus-bias.h by build.rs.\n");
    out.push_str("const THRESHOLD: &[u32] = &[\n");
    out.push_str(&threshold);
    out.push_str("\n];\n\n");
    out.push_str("const RAW_ESTIMATE_DATA: &[&[f64]] = &[\n");
    out.push_str(&format_nested_f64_arrays(&raw));
    out.push_str("];\n\n");
    out.push_str("const BIAS_DATA: &[&[f64]] = &[\n");
    out.push_str(&format_nested_f64_arrays(&bias));
    out.push_str("];\n");

    let out_dir = std::env::var("OUT_DIR").expect("OUT_DIR");
    std::fs::write(std::path::Path::new(&out_dir).join("hll_bias.rs"), out)
        .expect("write generated hll_bias.rs");
}

fn extract_braced_after(input: &str, marker: &str) -> Option<String> {
    let start = input.find(marker)?;
    let brace_start = input[start..].find('{')? + start;
    let mut depth = 0usize;
    let mut end = None;
    for (idx, ch) in input[brace_start..].char_indices() {
        match ch {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    end = Some(brace_start + idx);
                    break;
                }
            }
            _ => {}
        }
    }
    Some(input[brace_start + 1..end?].to_string())
}

fn format_nested_f64_arrays(input: &str) -> String {
    let mut arrays = Vec::new();
    let mut depth = 0usize;
    let mut start = None;
    for (idx, ch) in input.char_indices() {
        match ch {
            '{' => {
                if depth == 0 {
                    start = Some(idx + 1);
                }
                depth += 1;
            }
            '}' => {
                depth -= 1;
                if depth == 0 {
                    let body = &input[start.expect("array start")..idx];
                    arrays.push(body.trim());
                }
            }
            _ => {}
        }
    }

    let mut out = String::new();
    for arr in arrays {
        out.push_str("    &[\n");
        out.push_str(&format_f64_values(arr));
        out.push_str("\n    ],\n");
    }
    out
}

fn format_f64_values(input: &str) -> String {
    let mut out = String::new();
    for token in input.split(',') {
        let value = token.trim();
        if value.is_empty() {
            continue;
        }
        if !out.is_empty() {
            out.push_str(", ");
        }
        out.push_str(value);
        if !value.contains('.') && !value.contains('e') && !value.contains('E') {
            out.push_str(".0");
        }
    }
    out.push(',');
    out
}
