/// Translation map: AGCT ordering (not ACGT).
/// Index = 6-bit codon (2 bits per nucleotide, A=0 G=1 C=2 T=3).
const TRANSLATION_MAP: &[u8; 64] = b"KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

/// Forward lookup: maps ASCII character to 2-bit AGCT code for codon position.
/// A/a -> 0x00, G/g -> 0x01, C/c -> 0x02, T/t -> 0x03, else -> 0xFF
fn fwd_lookup(ch: u8) -> u8 {
    match ch {
        b'A' | b'a' => 0x00,
        b'G' | b'g' => 0x01,
        b'C' | b'c' => 0x02,
        b'T' | b't' => 0x03,
        _ => 0xFF,
    }
}

/// Reverse lookup: maps ASCII character to 2-bit complement code, shifted for prepending.
/// A/a -> 0x30, G/g -> 0x20, C/c -> 0x10, T/t -> 0x00, else -> 0xFF
fn rev_lookup(ch: u8) -> u8 {
    match ch {
        b'A' | b'a' => 0x30,
        b'G' | b'g' => 0x20,
        b'C' | b'c' => 0x10,
        b'T' | b't' => 0x00,
        _ => 0xFF,
    }
}

/// Translate a DNA sequence to all 6 reading frames (3 forward + 3 reverse).
/// Exact port of C++ `TranslateToAllFrames()` from `aa_translate.cc`.
///
/// Returns a Vec of 6 strings. Forward frames are indices 0-2, reverse frames are 3-5.
/// Reverse frames are built by prepending from the end of a max-sized buffer.
pub fn translate_to_all_frames(dna_seq: &str) -> Vec<String> {
    let dna = dna_seq.as_bytes();
    let max_size = dna.len() / 3 + 1;
    let mut aa_seqs: Vec<Vec<u8>> = (0..6).map(|_| vec![b' '; max_size]).collect();

    if dna.len() < 3 {
        return aa_seqs.into_iter().map(|v| String::from_utf8(v).unwrap()).collect();
    }

    let mut fwd_codon: u8 = 0;
    let mut rev_codon: u8 = 0;
    let mut ambig_nt_countdown: i32 = 0;
    let mut frame_len = [0usize; 6];

    for i in 0..dna.len() {
        let frame = i % 3;
        fwd_codon <<= 2;
        fwd_codon &= 0x3f;
        rev_codon >>= 2;

        if ambig_nt_countdown > 0 {
            ambig_nt_countdown -= 1;
        }

        let fwd_code = fwd_lookup(dna[i]);
        let rev_code = rev_lookup(dna[i]);

        if fwd_code == 0xFF {
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_code;
            rev_codon |= rev_code;
        }

        if i >= 2 {
            // Forward frame
            let ch = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[fwd_codon as usize]
            } else {
                b'X'
            };
            aa_seqs[frame][frame_len[frame]] = ch;
            frame_len[frame] += 1;

            // Reverse frame (prepend from end)
            let ch = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[rev_codon as usize]
            } else {
                b'X'
            };
            aa_seqs[frame + 3][max_size - 1 - frame_len[frame + 3]] = ch;
            frame_len[frame + 3] += 1;
        }
    }

    // Resize forward frames
    for i in 0..3 {
        aa_seqs[i].truncate(frame_len[i]);
    }
    // Trim leading padding from reverse frames
    for i in 3..6 {
        let start = aa_seqs[i].len() - frame_len[i];
        aa_seqs[i] = aa_seqs[i][start..].to_vec();
    }

    aa_seqs.into_iter().map(|v| String::from_utf8(v).unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_known_values() {
        // "AAAGGGCCCTTT" with AGCT ordering
        let result = translate_to_all_frames("AAAGGGCCCTTT");
        assert_eq!(result.len(), 6);
        // Frame 0 gets codons starting at i=2 where i%3==0 → i=3,6,9
        // The frame indexing: at i=2 frame=2 gets first codon, i=3 frame=0, etc.
        // Frame 0: codons from positions 1-3, 4-6, 7-9 → KGP (3 codons)
        assert_eq!(result[0], "KGP");
        // Forward frames should have protein chars
        for i in 0..3 {
            assert!(!result[i].is_empty(), "Forward frame {} empty", i);
        }
        // Reverse frames should also have content
        for i in 3..6 {
            assert!(!result[i].is_empty(), "Reverse frame {} empty", i);
        }
    }

    #[test]
    fn test_translate_short_sequences() {
        // Sequences shorter than 3 should return space-filled frames
        let r1 = translate_to_all_frames("A");
        assert!(r1.iter().all(|f| f.is_empty() || f.chars().all(|c| c == ' ')));

        let r2 = translate_to_all_frames("AT");
        assert!(r2.iter().all(|f| f.is_empty() || f.chars().all(|c| c == ' ')));

        // Exactly 3 bases: frame 2 gets first codon (i=2, frame=2%3=2)
        let r3 = translate_to_all_frames("ATG");
        assert_eq!(r3.len(), 6);
        // Frame 2 should have 1 amino acid, frames 0 and 1 should be empty
        assert!(r3[0].is_empty());
        assert!(r3[1].is_empty());
        assert_eq!(r3[2].len(), 1);
    }

    #[test]
    fn test_translate_ambiguous() {
        // N in sequence should produce X
        let result = translate_to_all_frames("NNNAAAGGGCCC");
        // First codon NNN -> X, then normal
        assert!(result[0].starts_with('X'));
    }

    #[test]
    fn test_translate_roundtrip_consistency() {
        // Verify same input always gives same output
        let seq = "ATGCGATCGATCGATCGATCG";
        let r1 = translate_to_all_frames(seq);
        let r2 = translate_to_all_frames(seq);
        assert_eq!(r1, r2);
        assert_eq!(r1.len(), 6);
        // All frames should have content for a sequence this long
        for (i, frame) in r1.iter().enumerate() {
            assert!(!frame.is_empty(), "Frame {} is empty for seq len {}", i, seq.len());
        }
    }
}
