/// Mask low-complexity regions in a DNA sequence using SDust algorithm.
/// Port of C++ k2mask.cc logic.
/// Masked positions are set to lowercase (or a replacement character).
pub fn mask_low_complexity(seq: &mut [u8], threshold: usize) {
    if seq.len() < 3 {
        return;
    }

    let len = seq.len();
    let window_size = 64.min(len);

    // Trinucleotide counts within a sliding window
    let mut counts = [0u32; 64]; // 4^3 = 64 possible trinucleotides

    // Compute trinucleotide code
    let encode = |a: u8, b: u8, c: u8| -> Option<usize> {
        let enc = |ch: u8| -> Option<u8> {
            match ch.to_ascii_uppercase() {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'T' => Some(3),
                _ => None,
            }
        };
        Some((enc(a)? as usize) * 16 + (enc(b)? as usize) * 4 + (enc(c)? as usize))
    };

    // Simple windowed dust scoring
    let mut masked = vec![false; len];

    for start in (0..len.saturating_sub(window_size)).step_by(window_size / 2) {
        let end = (start + window_size).min(len);
        if end - start < 3 {
            continue;
        }

        // Reset counts
        counts.fill(0);
        let mut valid_trinucs = 0u32;

        // Count trinucleotides in this window
        for i in start..end - 2 {
            if let Some(code) = encode(seq[i], seq[i + 1], seq[i + 2]) {
                counts[code] += 1;
                valid_trinucs += 1;
            }
        }

        if valid_trinucs == 0 {
            continue;
        }

        // Compute score: sum of c*(c-1)/2 for each trinucleotide count c
        let score: u32 = counts.iter()
            .map(|&c| c * c.saturating_sub(1) / 2)
            .sum();

        // If score exceeds threshold, mask the window
        let dust_score = (score as f64 * 10.0 / valid_trinucs as f64) as usize;
        if dust_score > threshold {
            for i in start..end {
                masked[i] = true;
            }
        }
    }

    // Apply masking
    for i in 0..len {
        if masked[i] {
            seq[i] = seq[i].to_ascii_lowercase();
        }
    }
}
