//! BLAST database to FASTA converter.
//! Port of C `blast_to_fasta.c`, `blast_defline.c`, `blast_utils.c`.
//!
//! Converts BLAST binary database files (.nin/.pin, .nhr/.phr, .nsq/.psq)
//! into standard FASTA format with optional taxid annotations.

use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// 2-bit nucleotide decoding table (BLAST uses NCBI2na encoding).
const NCBI2NA_DECODE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// 4-bit nucleotide decoding table (NCBI4na encoding for ambiguities).
#[allow(dead_code)]
const NCBI4NA_DECODE: [u8; 16] = [
    b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];
const MASK2DNA: [u8; 17] = [
    b'?', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
    b'N',
];

#[derive(Debug, Clone, Default)]
struct BlastString {
    string: Vec<u8>,
    cap: u32,
    len: u32,
}

struct BlastIdxData {
    idx_file: File,
    fmt_version: u32,
    db_seq_type: u32,
    volume: u32,
    title: BlastString,
    lmdb_file: BlastString,
    date: BlastString,
    num_oids: u32,
    vol_len: u64,
    max_seq_len: u32,
    hdr_arr: Vec<u8>,
    seq_arr: Vec<u8>,
    amb_arr: Vec<u8>,
}

struct BlastHdrData {
    hdr_file: Option<File>,
    fasta_hdr: BlastString,
    deflines: Vec<BlastDefline>,
    asn1: Asn1Data,
}

struct BlastSeqData {
    seq_file: File,
    buffer: BlastString,
    seq: BlastString,
    amb_data: Vec<Amb>,
    curr_pos: u32,
    amb_data_cap: u32,
}

type Integer = u32;
type VisibleString = Vec<u8>;

const SEQ_LOCAL: i32 = 0;
const SEQ_GIBBSQ: i32 = 1;
const SEQ_GIBBMT: i32 = 2;
const SEQ_GIIM: i32 = 3;
const SEQ_GENBANK: i32 = 4;
const SEQ_EMBL: i32 = 5;
const SEQ_PIR: i32 = 6;
const SEQ_SWISSPROT: i32 = 7;
const SEQ_PATENT: i32 = 8;
const SEQ_OTHER: i32 = 9;
const SEQ_GENERAL: i32 = 10;
const SEQ_GI: i32 = 11;
const SEQ_DDBJ: i32 = 12;
const SEQ_PRF: i32 = 13;
const SEQ_PDB: i32 = 14;
const SEQ_TPG: i32 = 15;
const SEQ_TPE: i32 = 16;
const SEQ_TPD: i32 = 17;
const SEQ_GPIPE: i32 = 18;
const SEQ_NAMED_ANNOT_TRACK: i32 = 19;
const SEQ_NONE: i32 = 20;

#[derive(Debug)]
struct Asn1Data {
    asn1_file: File,
    buffer: Vec<u8>,
    buf_len: usize,
    buf_cap: usize,
    curr_pos: usize,
    curr_oid: u32,
    block_lens: Vec<u32>,
    num_oids: u32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct DateStd {
    year: Integer,
    month: Integer,
    day: Integer,
    season: VisibleString,
    hour: Integer,
    minute: Integer,
    second: Integer,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct Date {
    str_value: VisibleString,
    std: DateStd,
    date_type: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct TextSeqId {
    name: VisibleString,
    acc: VisibleString,
    rel: VisibleString,
    ver: Integer,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct GiImportId {
    id: Integer,
    db: VisibleString,
    release: VisibleString,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct IdPat {
    country: VisibleString,
    number: VisibleString,
    app_number: VisibleString,
    doc_number_type: i32,
    doc_type: VisibleString,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct PatentSeqId {
    seqid: Integer,
    cit: IdPat,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct ObjectId {
    id: Integer,
    str_value: VisibleString,
    id_type: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct DbTag {
    db: VisibleString,
    tag: ObjectId,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct PdbSeqId {
    mol: VisibleString,
    chain: Integer,
    rel: Date,
    chain_id: VisibleString,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct SeqId {
    obj_id: ObjectId,
    int_id: Integer,
    db_tag_id: DbTag,
    pat_id: PatentSeqId,
    pdb_id: PdbSeqId,
    text_id: TextSeqId,
    giim_id: GiImportId,
    seq_id_type: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct BlastDefline {
    title: VisibleString,
    seq_ids: Vec<SeqId>,
    taxid: Integer,
    memberships: Vec<Integer>,
    links: Vec<Integer>,
    other_info: Vec<Integer>,
}

fn nlz(mut x: u32) -> u32 {
    if x == 0 {
        return 32;
    }

    let mut n = 1;
    if (x >> 16) == 0 {
        n += 16;
        x <<= 16;
    }
    if (x >> 24) == 0 {
        n += 8;
        x <<= 8;
    }
    if (x >> 28) == 0 {
        n += 4;
        x <<= 4;
    }
    if (x >> 30) == 0 {
        n += 2;
        x <<= 2;
    }
    n - (x >> 31)
}

fn next_power_of_2(n: u32) -> u32 {
    1u32 << (32 - nlz(n.wrapping_sub(1)))
}

fn alloc_memory(data: &mut Vec<u8>, element_size: u32, old_size: u32, new_size: u32) {
    if new_size > old_size {
        data.resize((element_size * new_size) as usize, 0);
    }
}

fn read_into_buffer<R: Read>(
    reader: &mut R,
    buffer: &mut [u8],
    element_size: u32,
    buffer_len: u32,
) -> io::Result<u32> {
    let bytes_to_read = (element_size * buffer_len) as usize;
    let bytes_read = reader.read(&mut buffer[..bytes_to_read])?;
    Ok((bytes_read / element_size as usize) as u32)
}

fn init_string() -> BlastString {
    BlastString::default()
}

fn string_length(s: &BlastString) -> u32 {
    s.len
}

fn string_capacity(s: &BlastString) -> u32 {
    s.cap
}

fn string_clear(s: &mut BlastString) {
    s.len = 0;
}

fn string_data(s: &BlastString) -> &[u8] {
    &s.string[..s.len as usize]
}

fn string_data_mut(s: &mut BlastString) -> &mut [u8] {
    let len = s.len as usize;
    &mut s.string[..len]
}

fn to_c_string(s: &mut BlastString) -> &[u8] {
    let len = s.len as usize;
    if s.string.len() <= len {
        s.string.resize(len + 1, 0);
        s.cap = s.cap.max((len + 1) as u32);
    }
    s.string[len] = 0;
    &s.string[..=len]
}

fn string_append_char(s: &mut BlastString, c: u8) {
    if s.len + 1 >= s.cap {
        let new_cap = next_power_of_2(s.cap + 1);
        alloc_memory(&mut s.string, 1, s.cap, new_cap);
        s.cap = new_cap;
    }
    s.string[s.len as usize] = c;
    s.len += 1;
}

fn string_append_int(s: &mut BlastString, mut num: u32) {
    if num == 0 {
        string_append_char(s, b'0');
        return;
    }

    let mut divisor = 1;
    while divisor <= num / 10 {
        divisor *= 10;
    }
    while divisor > 0 {
        let digit = (num / divisor) as u8 + b'0';
        num %= divisor;
        divisor /= 10;
        string_append_char(s, digit);
    }
}

fn string_append_str(s: &mut BlastString, string: &str) {
    for &byte in string.as_bytes() {
        string_append_char(s, byte);
    }
}

fn string_copy(
    s1: &mut BlastString,
    s2: &BlastString,
    s1_start: u32,
    s2_start: u32,
    mut len: u32,
) -> u32 {
    if s1_start > s1.len || s2_start > s2.len {
        return 0;
    }

    len = len.min(s2.len.saturating_sub(s2_start));
    let required_len = s1_start + len;
    if s1.cap < required_len + 1 {
        let new_len = next_power_of_2(required_len + 1);
        alloc_memory(&mut s1.string, 1, s1.cap, new_len);
        s1.cap = new_len;
    }

    let src = s2.string[s2_start as usize..(s2_start + len) as usize].to_vec();
    s1.string[s1_start as usize..required_len as usize].copy_from_slice(&src);
    s1.len = required_len;
    len
}

fn string_append(s1: &mut BlastString, s2: &BlastString, s2_start: u32, len: u32) -> u32 {
    string_copy(s1, s2, s1.len, s2_start, len)
}

fn string_reserve(s: &mut BlastString, size: u32) {
    alloc_memory(&mut s.string, 1, s.cap, size + 1);
    if size > s.cap {
        s.cap = size + 1;
    }
}

fn string_read_from_file<R: Read>(
    string: &mut BlastString,
    reader: &mut R,
    offset: u32,
    amt: u32,
) -> io::Result<()> {
    string_reserve(string, amt + offset);
    let start = offset as usize;
    let end = start + amt as usize;
    let bytes_read = read_into_buffer(reader, &mut string.string[start..end], 1, amt)?;
    string.len = bytes_read + offset;
    Ok(())
}

fn free_string(s: &mut BlastString) {
    s.string.clear();
    s.string.shrink_to_fit();
    s.cap = 0;
    s.len = 0;
}

fn ntoh_64(val: u64) -> u64 {
    u64::from_be(val)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Amb32 {
    offset: u32,
    length: u8,
    value: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Amb64 {
    offset: u32,
    unused: u16,
    length: u16,
    value: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Amb {
    Amb32(Amb32),
    Amb64(Amb64),
}

fn u32_to_amb32(n: u32) -> Amb32 {
    Amb32 {
        offset: n & 0x00ff_ffff,
        length: ((n >> 24) & 0x0f) as u8,
        value: ((n >> 28) & 0x0f) as u8,
    }
}

fn u64_to_amb64(n: u64) -> Amb64 {
    Amb64 {
        offset: (n & 0xffff_ffff) as u32,
        unused: ((n >> 32) & 0xffff) as u16,
        length: ((n >> 48) & 0x0fff) as u16,
        value: ((n >> 60) & 0x0f) as u8,
    }
}

fn read_int(file: &mut File) -> io::Result<u32> {
    read_be_u32(file)
}

fn read_long(file: &mut File) -> io::Result<u32> {
    let mut buf = [0u8; 8];
    file.read_exact(&mut buf)?;
    Ok(ntoh_64(u64::from_ne_bytes(buf)) as u32)
}

fn read_array(file: &mut File, element_size: u32, length: u32) -> io::Result<Vec<u8>> {
    let mut array = vec![0u8; (element_size * (length + 1)) as usize];
    let read_len = (element_size * length) as usize;
    file.read_exact(&mut array[..read_len])?;
    Ok(array)
}

fn open_file(filename: &str) -> io::Result<File> {
    File::open(filename)
}

fn read_string(file: &mut File) -> io::Result<BlastString> {
    let len = read_int(file)?;
    let mut string = init_string();
    string_reserve(&mut string, len + 1);
    file.read_exact(&mut string.string[..len as usize])?;
    string.len = len;
    Ok(string)
}

fn init_idx_data(filename: &str) -> io::Result<BlastIdxData> {
    Ok(BlastIdxData {
        idx_file: open_file(filename)?,
        fmt_version: 0,
        db_seq_type: 0,
        volume: 0,
        title: init_string(),
        lmdb_file: init_string(),
        date: init_string(),
        num_oids: 0,
        vol_len: 0,
        max_seq_len: 0,
        hdr_arr: Vec::new(),
        seq_arr: Vec::new(),
        amb_arr: Vec::new(),
    })
}

fn read_idx_data(idx_data: &mut BlastIdxData) -> io::Result<()> {
    idx_data.fmt_version = read_int(&mut idx_data.idx_file)?;
    idx_data.db_seq_type = read_int(&mut idx_data.idx_file)?;
    idx_data.volume = read_int(&mut idx_data.idx_file)?;
    idx_data.title = read_string(&mut idx_data.idx_file)?;
    idx_data.lmdb_file = read_string(&mut idx_data.idx_file)?;
    idx_data.date = read_string(&mut idx_data.idx_file)?;
    idx_data.num_oids = read_int(&mut idx_data.idx_file)?;
    idx_data.vol_len = read_long(&mut idx_data.idx_file)? as u64;
    idx_data.max_seq_len = read_int(&mut idx_data.idx_file)?;

    let array_bytes = 4 * (idx_data.num_oids + 1);
    idx_data.hdr_arr = read_array(&mut idx_data.idx_file, 1, array_bytes)?;
    idx_data.seq_arr = read_array(&mut idx_data.idx_file, 1, array_bytes)?;
    idx_data.amb_arr = read_array(&mut idx_data.idx_file, 1, array_bytes)?;
    Ok(())
}

fn free_idx_data(idx_data: &mut BlastIdxData) {
    free_string(&mut idx_data.title);
    free_string(&mut idx_data.lmdb_file);
    free_string(&mut idx_data.date);
    idx_data.hdr_arr.clear();
    idx_data.seq_arr.clear();
    idx_data.amb_arr.clear();
}

fn offset_at(offsets: &[u8], index: u32) -> u32 {
    let start = (index as usize) * 4;
    u32::from_be_bytes(offsets[start..start + 4].try_into().unwrap())
}

fn set_offset(offsets: &mut [u8], index: u32, value: u32) {
    let start = (index as usize) * 4;
    offsets[start..start + 4].copy_from_slice(&value.to_be_bytes());
}

fn init_asn1(filename: &str, block_lens: &[u32], num_oids: u32) -> io::Result<Asn1Data> {
    Ok(Asn1Data {
        asn1_file: open_file(filename)?,
        buffer: Vec::new(),
        buf_len: 0,
        buf_cap: 0,
        curr_pos: 0,
        curr_oid: 0,
        block_lens: block_lens.to_vec(),
        num_oids,
    })
}

fn free_asn1(asn1: &mut Asn1Data) {
    asn1.buffer.clear();
    asn1.buffer.shrink_to_fit();
    asn1.buf_len = 0;
    asn1.buf_cap = 0;
    asn1.curr_pos = 0;
}

fn load_more_data(asn1_data: &mut Asn1Data) -> io::Result<usize> {
    if asn1_data.curr_oid == asn1_data.num_oids {
        return Ok(0);
    }

    let block_len = asn1_data.block_lens[asn1_data.curr_oid as usize] as usize;
    asn1_data.curr_oid += 1;
    if (asn1_data.buf_cap - asn1_data.buf_len) < block_len {
        let new_cap = next_power_of_2((asn1_data.buf_len + block_len) as u32) as usize;
        alloc_memory(
            &mut asn1_data.buffer,
            1,
            asn1_data.buf_cap as u32,
            new_cap as u32,
        );
        asn1_data.buf_cap = new_cap;
    }
    read_into_buffer(
        &mut asn1_data.asn1_file,
        &mut asn1_data.buffer[asn1_data.buf_len..asn1_data.buf_len + block_len],
        1,
        block_len as u32,
    )?;
    asn1_data.buf_len += block_len;

    Ok(block_len)
}

fn asn1_get_byte(asn1_data: &mut Asn1Data) -> io::Result<u8> {
    if asn1_data.curr_pos == asn1_data.buf_len {
        load_more_data(asn1_data)?;
    }

    let byte = asn1_data.buffer[asn1_data.curr_pos];
    asn1_data.curr_pos += 1;
    Ok(byte)
}

fn asn1_backtrack(asn1_data: &mut Asn1Data) {
    asn1_data.curr_pos -= 1;
}

fn asn1_peek_byte(asn1_data: &mut Asn1Data) -> io::Result<u8> {
    if asn1_data.curr_pos == asn1_data.buf_len {
        load_more_data(asn1_data)?;
    }

    Ok(asn1_data.buffer[asn1_data.curr_pos])
}

fn asn1_sequence_start(asn1_data: &mut Asn1Data) -> io::Result<bool> {
    if asn1_peek_byte(asn1_data)? != 0x30 {
        return Ok(false);
    }

    let _tag = asn1_get_byte(asn1_data)?;
    let len = asn1_get_byte(asn1_data)?;
    Ok(len == 0x80)
}

fn asn1_indefinite_tag_end(asn1_data: &mut Asn1Data) -> io::Result<bool> {
    let null_octet1 = asn1_get_byte(asn1_data)?;
    let null_octet2 = asn1_get_byte(asn1_data)?;
    Ok(null_octet1 == 0 && null_octet2 == 0)
}

fn asn1_is_end_of_sequence(asn1_data: &mut Asn1Data) -> io::Result<bool> {
    let end = asn1_indefinite_tag_end(asn1_data)?;
    asn1_backtrack(asn1_data);
    asn1_backtrack(asn1_data);
    Ok(end)
}

fn asn1_get_visible_string(
    asn1_data: &mut Asn1Data,
    string: &mut VisibleString,
) -> io::Result<u32> {
    let mut str_len = 0u32;
    let tag = asn1_get_byte(asn1_data)?;
    debug_assert_eq!(tag, 0x1A);
    let mut octets = u32::from(asn1_get_byte(asn1_data)?);

    if (octets & 0x80) == 0x80 {
        let num_octets = octets & (0x80 - 1);
        for j in 0..num_octets {
            octets = u32::from(asn1_get_byte(asn1_data)?);
            str_len = (str_len << (j * 8)) | octets;
        }
    } else {
        str_len = octets;
    }
    if (asn1_data.curr_pos + str_len as usize) > asn1_data.buf_len {
        load_more_data(asn1_data)?;
    }

    let start = asn1_data.curr_pos;
    let end = start + str_len as usize;
    *string = asn1_data.buffer[start..end].to_vec();
    asn1_data.curr_pos = end;

    Ok(str_len)
}

fn asn1_get_integer(asn1_data: &mut Asn1Data, num: &mut Integer) -> io::Result<i32> {
    let tag = asn1_get_byte(asn1_data)?;
    debug_assert_eq!(tag, 0x02);

    let int_len = u32::from(asn1_get_byte(asn1_data)?);
    let mut value = u32::from(asn1_get_byte(asn1_data)?);
    for _ in 1..int_len {
        value = (value << 8) + u32::from(asn1_get_byte(asn1_data)?);
    }

    *num = value;
    Ok(1)
}

fn get_explicit_tag(asn1_data: &mut Asn1Data) -> io::Result<u8> {
    let field_no = asn1_get_byte(asn1_data)?;
    let _ = asn1_get_byte(asn1_data)?;
    Ok(field_no)
}

fn asn1_get_optional_visible_string_field(
    asn1_data: &mut Asn1Data,
    field_tag: u8,
    destination: &mut VisibleString,
) -> io::Result<()> {
    if asn1_peek_byte(asn1_data)? == field_tag {
        let _ = get_explicit_tag(asn1_data)?;
        let _ = asn1_get_visible_string(asn1_data, destination)?;
        let _ = asn1_indefinite_tag_end(asn1_data)?;
    }
    Ok(())
}

fn asn1_get_optional_integer_field(
    asn1_data: &mut Asn1Data,
    field_tag: u8,
    destination: &mut Integer,
) -> io::Result<()> {
    if asn1_peek_byte(asn1_data)? == field_tag {
        let _ = get_explicit_tag(asn1_data)?;
        let _ = asn1_get_integer(asn1_data, destination)?;
        let _ = asn1_indefinite_tag_end(asn1_data)?;
    }
    Ok(())
}

fn asn1_get_mandatory_integer_field(
    asn1_data: &mut Asn1Data,
    destination: &mut Integer,
) -> io::Result<()> {
    let _ = get_explicit_tag(asn1_data)?;
    let _ = asn1_get_integer(asn1_data, destination)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;
    Ok(())
}

fn asn1_get_mandatory_visible_string_field(
    asn1_data: &mut Asn1Data,
    destination: &mut VisibleString,
) -> io::Result<()> {
    let _ = get_explicit_tag(asn1_data)?;
    let _ = asn1_get_visible_string(asn1_data, destination)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;
    Ok(())
}

fn get_date(asn1_data: &mut Asn1Data, d: &mut Date) -> io::Result<i32> {
    let field_no = asn1_peek_byte(asn1_data)?;
    if field_no == 0xA0 {
        asn1_get_mandatory_visible_string_field(asn1_data, &mut d.str_value)?;
        d.date_type = 1;
    } else {
        let _ = get_explicit_tag(asn1_data)?;
        let _ = asn1_sequence_start(asn1_data)?;
        asn1_get_mandatory_integer_field(asn1_data, &mut d.std.year)?;
        asn1_get_optional_integer_field(asn1_data, 0xA1, &mut d.std.month)?;
        asn1_get_optional_integer_field(asn1_data, 0xA2, &mut d.std.day)?;
        asn1_get_optional_visible_string_field(asn1_data, 0xA3, &mut d.std.season)?;
        asn1_get_optional_integer_field(asn1_data, 0xA4, &mut d.std.hour)?;
        asn1_get_optional_integer_field(asn1_data, 0xA5, &mut d.std.minute)?;
        asn1_get_optional_integer_field(asn1_data, 0xA6, &mut d.std.second)?;
        let _ = asn1_indefinite_tag_end(asn1_data)?;
        let _ = asn1_indefinite_tag_end(asn1_data)?;
        d.date_type = 2;
    }

    Ok(1)
}

fn get_db_tag(asn1_data: &mut Asn1Data, tag: &mut DbTag) -> io::Result<i32> {
    *tag = DbTag::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_mandatory_visible_string_field(asn1_data, &mut tag.db)?;
    let _ = get_object_id(asn1_data, &mut tag.tag)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_pdb_seq_id(asn1_data: &mut Asn1Data, seqid: &mut PdbSeqId) -> io::Result<i32> {
    *seqid = PdbSeqId::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_mandatory_visible_string_field(asn1_data, &mut seqid.mol)?;
    asn1_get_optional_integer_field(asn1_data, 0xA1, &mut seqid.chain)?;
    if asn1_peek_byte(asn1_data)? == 0xA2 {
        let _ = get_explicit_tag(asn1_data)?;
        let _ = get_date(asn1_data, &mut seqid.rel)?;
        let _ = asn1_indefinite_tag_end(asn1_data)?;
    }
    asn1_get_optional_visible_string_field(asn1_data, 0xA3, &mut seqid.chain_id)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_gi_import_id(asn1_data: &mut Asn1Data, seqid: &mut GiImportId) -> io::Result<i32> {
    *seqid = GiImportId::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_mandatory_integer_field(asn1_data, &mut seqid.id)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA1, &mut seqid.db)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA2, &mut seqid.release)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_text_seq_id(asn1_data: &mut Asn1Data, seqid: &mut TextSeqId) -> io::Result<i32> {
    *seqid = TextSeqId::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA0, &mut seqid.name)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA1, &mut seqid.acc)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA2, &mut seqid.rel)?;
    asn1_get_optional_integer_field(asn1_data, 0xA3, &mut seqid.ver)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_patent_seq_id(asn1_data: &mut Asn1Data, seqid: &mut PatentSeqId) -> io::Result<i32> {
    *seqid = PatentSeqId::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_mandatory_integer_field(asn1_data, &mut seqid.seqid)?;
    let _ = get_id_pat(asn1_data, &mut seqid.cit)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_object_id(asn1_data: &mut Asn1Data, id: &mut ObjectId) -> io::Result<i32> {
    *id = ObjectId::default();

    let _ = get_explicit_tag(asn1_data)?;
    if asn1_peek_byte(asn1_data)? == 0xA0 {
        asn1_get_mandatory_integer_field(asn1_data, &mut id.id)?;
        id.id_type = 1;
    } else {
        asn1_get_mandatory_visible_string_field(asn1_data, &mut id.str_value)?;
        id.id_type = 2;
    }
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_id_pat(asn1_data: &mut Asn1Data, id: &mut IdPat) -> io::Result<i32> {
    *id = IdPat::default();

    let _ = asn1_sequence_start(asn1_data)?;
    asn1_get_optional_visible_string_field(asn1_data, 0xA0, &mut id.country)?;
    if asn1_peek_byte(asn1_data)? == 0xA1 {
        asn1_get_mandatory_visible_string_field(asn1_data, &mut id.number)?;
        id.doc_number_type = 1;
    } else {
        asn1_get_mandatory_visible_string_field(asn1_data, &mut id.app_number)?;
        id.doc_number_type = 2;
    }
    asn1_get_optional_visible_string_field(asn1_data, 0xA3, &mut id.doc_type)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_seq_id(asn1_data: &mut Asn1Data, seqid: &mut SeqId) -> io::Result<i32> {
    let field_no = asn1_peek_byte(asn1_data)?;
    match field_no {
        0xA0 => {
            let _ = get_object_id(asn1_data, &mut seqid.obj_id)?;
            seqid.seq_id_type = SEQ_LOCAL;
        }
        0xA1 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = asn1_get_integer(asn1_data, &mut seqid.int_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_GIBBSQ;
        }
        0xA2 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = asn1_get_integer(asn1_data, &mut seqid.int_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_GIBBMT;
        }
        0xA3 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_gi_import_id(asn1_data, &mut seqid.giim_id)?;
            seqid.seq_id_type = SEQ_GIIM;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA4 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            seqid.seq_id_type = SEQ_GENBANK;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA5 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            seqid.seq_id_type = SEQ_EMBL;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA6 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            seqid.seq_id_type = SEQ_PIR;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA7 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            seqid.seq_id_type = SEQ_SWISSPROT;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA8 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_patent_seq_id(asn1_data, &mut seqid.pat_id)?;
            seqid.seq_id_type = SEQ_PATENT;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xA9 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            seqid.seq_id_type = SEQ_OTHER;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xAA => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_db_tag(asn1_data, &mut seqid.db_tag_id)?;
            seqid.seq_id_type = SEQ_GENERAL;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
        }
        0xAB => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = asn1_get_integer(asn1_data, &mut seqid.int_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_GI;
        }
        0xAC => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_DDBJ;
        }
        0xAD => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_PRF;
        }
        0xAE => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_pdb_seq_id(asn1_data, &mut seqid.pdb_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_PDB;
        }
        0xAF => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_TPG;
        }
        0xB0 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_TPE;
        }
        0xB1 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_TPD;
        }
        0xB2 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_GPIPE;
        }
        0xB3 => {
            let _ = get_explicit_tag(asn1_data)?;
            let _ = get_text_seq_id(asn1_data, &mut seqid.text_id)?;
            let _ = asn1_indefinite_tag_end(asn1_data)?;
            seqid.seq_id_type = SEQ_NAMED_ANNOT_TRACK;
        }
        _ => {
            seqid.seq_id_type = SEQ_NONE;
        }
    }

    Ok(seqid.seq_id_type)
}

fn get_ints(asn1_data: &mut Asn1Data, vec: &mut Vec<Integer>) -> io::Result<i32> {
    let mut value = 0;
    let _ = asn1_sequence_start(asn1_data)?;
    loop {
        if asn1_peek_byte(asn1_data)? == 0x02 {
            let _ = asn1_get_integer(asn1_data, &mut value)?;
            vec.push(value);
        } else {
            break;
        }
    }
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_memberships(asn1_data: &mut Asn1Data, memberships: &mut Vec<Integer>) -> io::Result<i32> {
    if asn1_peek_byte(asn1_data)? != 0xA3 {
        return Ok(0);
    }

    let _ = get_explicit_tag(asn1_data)?;
    let _ = get_ints(asn1_data, memberships)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_links(asn1_data: &mut Asn1Data, links: &mut Vec<Integer>) -> io::Result<i32> {
    if asn1_peek_byte(asn1_data)? != 0xA4 {
        return Ok(0);
    }

    let _ = get_explicit_tag(asn1_data)?;
    let _ = get_ints(asn1_data, links)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_other_info(asn1_data: &mut Asn1Data, other_info: &mut Vec<Integer>) -> io::Result<i32> {
    if asn1_peek_byte(asn1_data)? != 0xA5 {
        return Ok(0);
    }

    let _ = get_explicit_tag(asn1_data)?;
    let _ = get_ints(asn1_data, other_info)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn get_blast_defline(asn1_data: &mut Asn1Data, defline: &mut BlastDefline) -> io::Result<i32> {
    if !asn1_sequence_start(asn1_data)? {
        return Ok(0);
    }
    asn1_get_optional_visible_string_field(asn1_data, 0xA0, &mut defline.title)?;
    let _ = get_explicit_tag(asn1_data)?;
    let _ = asn1_sequence_start(asn1_data)?;
    loop {
        let mut id = SeqId::default();
        let seq_id_type = get_seq_id(asn1_data, &mut id)?;
        if seq_id_type == SEQ_NONE {
            break;
        }
        defline.seq_ids.push(id);
    }
    let _ = asn1_indefinite_tag_end(asn1_data)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;
    asn1_get_optional_integer_field(asn1_data, 0xA2, &mut defline.taxid)?;
    let _ = get_memberships(asn1_data, &mut defline.memberships)?;
    let _ = get_links(asn1_data, &mut defline.links)?;
    let _ = get_other_info(asn1_data, &mut defline.other_info)?;
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(1)
}

fn init_blast_defline() -> BlastDefline {
    BlastDefline::default()
}

fn reset_blast_defline(defline: &mut BlastDefline) {
    defline.title.clear();
    defline.seq_ids.clear();
    defline.memberships.clear();
    defline.links.clear();
    defline.other_info.clear();
    defline.taxid = 0;
}

fn destroy_blast_defline(defline: &mut BlastDefline) {
    reset_blast_defline(defline);
    defline.seq_ids.shrink_to_fit();
    defline.memberships.shrink_to_fit();
    defline.links.shrink_to_fit();
    defline.other_info.shrink_to_fit();
}

fn get_blast_deflines(
    asn1_data: &mut Asn1Data,
    deflines: &mut Vec<BlastDefline>,
) -> io::Result<i32> {
    let mut i = 0usize;

    asn1_data.buf_len = 0;
    asn1_data.curr_pos = 0;

    if load_more_data(asn1_data)? == 0 {
        return Ok(0);
    }
    let _ = asn1_sequence_start(asn1_data)?;
    loop {
        if i == deflines.len() {
            deflines.push(init_blast_defline());
        }
        reset_blast_defline(&mut deflines[i]);
        let retval = get_blast_defline(asn1_data, &mut deflines[i])?;
        if retval == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "error parsing defline",
            ));
        }
        i += 1;
        if asn1_is_end_of_sequence(asn1_data)? {
            break;
        }
    }
    let _ = asn1_indefinite_tag_end(asn1_data)?;

    Ok(i as i32)
}

fn init_hdr_data(
    filename: &str,
    hdr_offsets: &mut [u8],
    num_oids: u32,
) -> io::Result<BlastHdrData> {
    for i in 1..(num_oids + 1) {
        let delta = offset_at(hdr_offsets, i) - offset_at(hdr_offsets, i - 1);
        set_offset(hdr_offsets, i - 1, delta);
    }
    let mut block_lens = Vec::with_capacity(num_oids as usize);
    for i in 0..num_oids {
        block_lens.push(offset_at(hdr_offsets, i));
    }
    set_offset(hdr_offsets, num_oids, 0);

    Ok(BlastHdrData {
        hdr_file: None,
        fasta_hdr: init_string(),
        deflines: Vec::new(),
        asn1: init_asn1(filename, &block_lens, num_oids)?,
    })
}

fn free_hdr_data(hdr: &mut BlastHdrData) {
    for defline in &mut hdr.deflines {
        destroy_blast_defline(defline);
    }
    hdr.deflines.clear();
    free_string(&mut hdr.fasta_hdr);
    free_asn1(&mut hdr.asn1);
    hdr.hdr_file = None;
}

fn db_tag_to_string(s: &mut BlastString, tag: &DbTag) {
    string_append_bytes(s, &tag.db);
    string_append_char(s, b'_');
    if tag.tag.id_type == 1 {
        string_append_int(s, tag.tag.id);
    } else {
        string_append_bytes(s, &tag.tag.str_value);
    }
}

fn text_seq_id_to_string(s: &mut BlastString, seqid: &TextSeqId) {
    string_append_bytes(s, &seqid.acc);
    if seqid.ver > 0 {
        string_append_char(s, b'.');
        string_append_int(s, seqid.ver);
    }
}

fn pdb_seq_id_to_string(s: &mut BlastString, seqid: &PdbSeqId) {
    string_append_bytes(s, &seqid.mol);
    if !seqid.chain_id.is_empty() {
        string_append_char(s, b'_');
        string_append_bytes(s, &seqid.chain_id);
    } else if seqid.chain != 0 {
        string_append_char(s, b'_');
        string_append_char(s, seqid.chain as u8);
    }
}

fn defline_to_header(
    s: &mut BlastString,
    deflines: &[BlastDefline],
    include_taxid: bool,
    curr_defline: usize,
) {
    string_append_char(s, b'>');
    let defline = &deflines[curr_defline];

    if include_taxid && defline.taxid != 0 {
        string_append_str(s, "kraken:taxid|");
        string_append_int(s, defline.taxid);
        string_append_char(s, b'|');
    }

    for seq_id in &defline.seq_ids {
        match seq_id.seq_id_type {
            SEQ_GI => {}
            SEQ_PDB => pdb_seq_id_to_string(s, &seq_id.pdb_id),
            SEQ_GENERAL => db_tag_to_string(s, &seq_id.db_tag_id),
            _ => text_seq_id_to_string(s, &seq_id.text_id),
        }
    }

    string_append_char(s, b' ');
    string_append_bytes(s, &defline.title);
}

fn deflines_to_header(
    s: &mut BlastString,
    deflines: &[BlastDefline],
    include_taxid: bool,
    num_deflines: i32,
) {
    for i in 0..(num_deflines as usize) {
        if string_length(s) > 0 {
            string_append_char(s, b' ');
        }
        defline_to_header(s, deflines, include_taxid, i);
    }
}

fn get_deflines(hdr_data: &mut BlastHdrData) -> io::Result<i32> {
    get_blast_deflines(&mut hdr_data.asn1, &mut hdr_data.deflines)
}

fn init_seq_data(filename: &str, max_seq_len: u32, max_buf_size: u32) -> io::Result<BlastSeqData> {
    let mut seq_file = open_file(filename)?;
    let mut first_byte = [0u8; 1];
    seq_file.read_exact(&mut first_byte)?;

    let mut seq = init_string();
    string_reserve(&mut seq, max_seq_len + 1);
    let mut buffer = init_string();
    string_reserve(&mut buffer, max_buf_size + 1);

    Ok(BlastSeqData {
        seq_file,
        buffer,
        seq,
        amb_data: Vec::new(),
        curr_pos: 0,
        amb_data_cap: 0,
    })
}

fn max_block_size(offsets: &[u8], length: u32) -> u32 {
    let mut max_length = 0;
    for i in 1..length {
        let delta = offset_at(offsets, i) - offset_at(offsets, i - 1);
        if delta > max_length {
            max_length = delta;
        }
    }
    max_length
}

fn has_ambiguous_data(block_end: u32, amb_start: u32) -> bool {
    block_end != amb_start
}

fn get_nucleotide_length(data: &[u8], data_length: usize) -> u32 {
    let remainder = data[data_length] as u32;
    ((data_length - 1) as u32) * 4 + remainder
}

fn read_seq_block(seq: &mut BlastSeqData, block_length: u32) -> io::Result<()> {
    string_clear(&mut seq.buffer);
    string_reserve(&mut seq.buffer, block_length + 1);
    seq.buffer.len = read_into_buffer(
        &mut seq.seq_file,
        &mut seq.buffer.string[..block_length as usize],
        1,
        block_length,
    )?;
    Ok(())
}

fn read_ambiguous_data(
    input: &[u8],
    output: &mut Vec<Amb>,
    output_length: u32,
    format: &mut u32,
) -> u32 {
    let mut amb_data_len = u32::from_be_bytes(input[..4].try_into().unwrap());
    *format = u32::from((amb_data_len & 0x8000_0000) == 0x8000_0000);
    amb_data_len &= 0x7fff_ffff;
    if amb_data_len > output_length {
        output.resize(
            amb_data_len as usize,
            Amb::Amb32(Amb32 {
                offset: 0,
                length: 0,
                value: 0,
            }),
        );
    }

    let mut pos = 4usize;
    if *format == 1 {
        let count = amb_data_len / 2;
        for i in 0..count as usize {
            let segment = u64::from_ne_bytes(input[pos..pos + 8].try_into().unwrap());
            output[i] = Amb::Amb64(u64_to_amb64(ntoh_64(segment)));
            pos += 8;
        }
        count
    } else {
        for i in 0..amb_data_len as usize {
            let segment = u32::from_be_bytes(input[pos..pos + 4].try_into().unwrap());
            output[i] = Amb::Amb32(u32_to_amb32(segment));
            pos += 4;
        }
        amb_data_len
    }
}

fn reconstruct_sequence(
    seq: &mut BlastSeqData,
    mut unamb_data_len: u32,
    amb_data_len: u32,
    format: u32,
) {
    let buffer = string_data(&seq.buffer);
    let final_byte = buffer[(unamb_data_len - 1) as usize];
    unamb_data_len -= 1;
    let mut nucs_in_final_byte = u32::from(final_byte & 0x3);
    let nuc_len = unamb_data_len * 4 + nucs_in_final_byte;

    string_reserve(&mut seq.seq, nuc_len + 1);
    let seq_buf = &mut seq.seq.string;
    let mut j = 0usize;
    for i in 0..(unamb_data_len as usize) {
        let byte = buffer[i];
        seq_buf[j] = NCBI2NA_DECODE[((byte >> 6) & 0x3) as usize];
        seq_buf[j + 1] = NCBI2NA_DECODE[((byte >> 4) & 0x3) as usize];
        seq_buf[j + 2] = NCBI2NA_DECODE[((byte >> 2) & 0x3) as usize];
        seq_buf[j + 3] = NCBI2NA_DECODE[(byte & 0x3) as usize];
        j += 4;
    }

    let mut k = 3i32;
    while final_byte != 0 && nucs_in_final_byte > 0 {
        seq_buf[j] = NCBI2NA_DECODE[((final_byte >> ((k << 1) as u8)) & 0x3) as usize];
        j += 1;
        nucs_in_final_byte -= 1;
        k -= 1;
    }

    if format == 1 {
        for amb in seq.amb_data.iter().take(amb_data_len as usize) {
            if let Amb::Amb64(a) = amb {
                for pos in a.offset..=(a.offset + u32::from(a.length)) {
                    seq_buf[pos as usize] = MASK2DNA[a.value as usize];
                }
            }
        }
    } else {
        for amb in seq.amb_data.iter().take(amb_data_len as usize) {
            if let Amb::Amb32(a) = amb {
                for pos in a.offset..=(a.offset + u32::from(a.length)) {
                    seq_buf[pos as usize] = MASK2DNA[a.value as usize];
                }
            }
        }
    }

    seq_buf[nuc_len as usize] = 0;
    seq.seq.len = nuc_len;
}

fn next_sequence(seq: &mut BlastSeqData, idx_data: &BlastIdxData) -> io::Result<i32> {
    let i = seq.curr_pos;
    let mut format = 0;
    let mut amb_data_len = 0;
    let seq_start = offset_at(&idx_data.seq_arr, i);
    let seq_end = offset_at(&idx_data.seq_arr, i + 1);
    let amb_start = offset_at(&idx_data.amb_arr, i);

    read_seq_block(seq, seq_end - seq_start)?;
    if has_ambiguous_data(seq_end, amb_start) {
        let buffer = string_data(&seq.buffer);
        amb_data_len = read_ambiguous_data(
            &buffer[(amb_start - seq_start) as usize..],
            &mut seq.amb_data,
            seq.amb_data_cap,
            &mut format,
        );
        if amb_data_len > seq.amb_data_cap {
            seq.amb_data_cap = amb_data_len;
        }
    }

    let unamb_data_len = amb_start - seq_start;
    reconstruct_sequence(seq, unamb_data_len, amb_data_len, format);
    seq.curr_pos += 1;

    Ok(0)
}

fn write_sequence(seq: &BlastString, seq_width: usize, out: &mut dyn Write) -> io::Result<()> {
    let raw_seq = string_data(seq);
    let full_width_lines = (seq.len as usize) / seq_width;
    let remainder = (seq.len as usize) % seq_width;

    for i in 0..full_width_lines {
        out.write_all(&raw_seq[i * seq_width..(i + 1) * seq_width])?;
        out.write_all(b"\n")?;
    }
    if remainder > 0 {
        out.write_all(
            &raw_seq[full_width_lines * seq_width..full_width_lines * seq_width + remainder],
        )?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn free_seq_data(seq_data: &mut BlastSeqData) {
    free_string(&mut seq_data.buffer);
    free_string(&mut seq_data.seq);
    seq_data.amb_data.clear();
    seq_data.amb_data.shrink_to_fit();
}

fn string_append_bytes(s: &mut BlastString, bytes: &[u8]) {
    for &byte in bytes {
        string_append_char(s, byte);
    }
}

fn c_string_bytes(s: &mut BlastString) -> &[u8] {
    let data = to_c_string(s);
    &data[..data.len().saturating_sub(1)]
}

/// BLAST database index header.
struct BlastIndex {
    fmt_version: u32,
    db_seq_type: u32, // 0 = nucleotide, 1 = protein
    num_oids: u32,
    title: String,
    #[allow(dead_code)]
    date: String,
    hdr_offsets: Vec<u32>,
    seq_offsets: Vec<u32>,
    #[allow(dead_code)]
    amb_offsets: Vec<u32>,
}

/// Read a 4-byte big-endian u32.
fn read_be_u32(file: &mut File) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    file.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

/// Read a 4-byte big-endian u32 as length-prefixed string.
fn read_length_prefixed_string(file: &mut File) -> io::Result<String> {
    let len = read_be_u32(file)? as usize;
    let mut buf = vec![0u8; len];
    file.read_exact(&mut buf)?;
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Parse a BLAST index file (.nin or .pin).
fn parse_blast_index(idx_path: &str) -> io::Result<BlastIndex> {
    let mut file = File::open(idx_path)?;

    let fmt_version = read_be_u32(&mut file)?;
    let db_seq_type = read_be_u32(&mut file)?;
    let _volume = read_be_u32(&mut file)?;
    let title = read_length_prefixed_string(&mut file)?;
    let date = read_length_prefixed_string(&mut file)?;
    let num_oids = read_be_u32(&mut file)?;
    let _vol_len_hi = read_be_u32(&mut file)?;
    let _vol_len_lo = read_be_u32(&mut file)?;
    let _max_seq_len = read_be_u32(&mut file)?;

    // Read offset arrays
    let mut hdr_offsets = Vec::with_capacity(num_oids as usize + 1);
    for _ in 0..=num_oids {
        hdr_offsets.push(read_be_u32(&mut file)?);
    }

    let mut seq_offsets = Vec::with_capacity(num_oids as usize + 1);
    for _ in 0..=num_oids {
        seq_offsets.push(read_be_u32(&mut file)?);
    }

    let mut amb_offsets = Vec::with_capacity(num_oids as usize + 1);
    if db_seq_type == 0 {
        // Nucleotide databases have ambiguity offsets
        for _ in 0..=num_oids {
            amb_offsets.push(read_be_u32(&mut file)?);
        }
    }

    Ok(BlastIndex {
        fmt_version,
        db_seq_type,
        num_oids,
        title,
        date,
        hdr_offsets,
        seq_offsets,
        amb_offsets,
    })
}

/// Decode a 2-bit packed nucleotide sequence.
fn decode_ncbi2na(packed: &[u8], length: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(length);
    for &byte in packed {
        for shift in (0..4).rev() {
            if seq.len() >= length {
                break;
            }
            let code = (byte >> (shift * 2)) & 0x03;
            seq.push(NCBI2NA_DECODE[code as usize]);
        }
    }
    seq.truncate(length);
    seq
}

/// Convert a BLAST database to FASTA format.
///
/// `db_prefix` is the path without extension (e.g., "mydb" for mydb.nin, mydb.nhr, mydb.nsq).
/// `output_path` is the FASTA output file path.
/// `include_taxid` adds kraken:taxid| prefix to headers if true.
pub fn blast_to_fasta(db_prefix: &str, output_path: &str, _include_taxid: bool) -> io::Result<()> {
    // Determine file extensions based on database type
    let idx_nuc = format!("{db_prefix}.nin");
    let idx_pro = format!("{db_prefix}.pin");

    let (idx_path, hdr_ext, seq_ext) = if Path::new(&idx_nuc).exists() {
        (idx_nuc, "nhr", "nsq")
    } else if Path::new(&idx_pro).exists() {
        (idx_pro, "phr", "psq")
    } else {
        return Err(io::Error::other(format!(
            "Cannot find BLAST index at {db_prefix}.nin or {db_prefix}.pin"
        )));
    };

    let idx = parse_blast_index(&idx_path)?;
    eprintln!("BLAST database: {}", idx.title);
    eprintln!(
        "  Type: {}",
        if idx.db_seq_type == 0 {
            "nucleotide"
        } else {
            "protein"
        }
    );
    eprintln!("  Sequences: {}", idx.num_oids);
    eprintln!("  Format version: {:#x}", idx.fmt_version);

    if idx.db_seq_type != 0 {
        let hdr_path = format!("{db_prefix}.{hdr_ext}");
        let seq_path = format!("{db_prefix}.{seq_ext}");
        let mut hdr_file = File::open(&hdr_path)?;
        let mut seq_file = File::open(&seq_path)?;
        let mut output = BufWriter::new(File::create(output_path)?);

        for oid in 0..idx.num_oids as usize {
            let hdr_start = idx.hdr_offsets[oid] as u64;
            let hdr_end = idx.hdr_offsets[oid + 1] as u64;
            let hdr_len = (hdr_end - hdr_start) as usize;
            hdr_file.seek(SeekFrom::Start(hdr_start))?;
            let mut hdr_buf = vec![0u8; hdr_len];
            hdr_file.read_exact(&mut hdr_buf)?;
            let header_text = extract_readable_header(&hdr_buf);

            let seq_start = idx.seq_offsets[oid] as u64;
            let seq_end = idx.seq_offsets[oid + 1] as u64;
            let seq_len = (seq_end - seq_start) as usize;
            seq_file.seek(SeekFrom::Start(seq_start))?;
            let mut seq_buf = vec![0u8; seq_len];
            seq_file.read_exact(&mut seq_buf)?;
            let seq: Vec<u8> = seq_buf.iter().map(|&b| ncbi_aa_decode(b)).collect();
            writeln!(output, ">{header_text}")?;
            for chunk in seq.chunks(70) {
                output.write_all(chunk)?;
                output.write_all(b"\n")?;
            }
        }

        eprintln!("Wrote {} sequences to {output_path}", idx.num_oids);
        return Ok(());
    }

    let idx_path = format!("{db_prefix}.nin");
    let hdr_path = format!("{db_prefix}.nhr");
    let seq_path = format!("{db_prefix}.nsq");

    let mut idx_data = init_idx_data(&idx_path)?;
    read_idx_data(&mut idx_data)?;
    let num_oids = idx_data.num_oids;
    let max_seq_len = idx_data.max_seq_len;
    let max_seq_block_size = next_power_of_2(max_block_size(&idx_data.seq_arr, idx_data.num_oids));
    let mut hdr_data = init_hdr_data(&hdr_path, &mut idx_data.hdr_arr, idx_data.num_oids)?;
    let mut seq_data = init_seq_data(&seq_path, max_seq_len, max_seq_block_size)?;
    let mut out_file = BufWriter::new(File::create(output_path)?);

    eprintln!("The version of this database is {:X}", idx_data.fmt_version);
    eprintln!(
        "Sequence type is: {}",
        if idx_data.db_seq_type == 0 {
            "nucleotide"
        } else {
            "protein"
        }
    );
    eprintln!("Volume number: {}", idx_data.volume);
    eprintln!(
        "Title of volume: {}",
        String::from_utf8_lossy(c_string_bytes(&mut idx_data.title))
    );
    eprintln!(
        "Date created: {}",
        String::from_utf8_lossy(c_string_bytes(&mut idx_data.date))
    );
    eprintln!("Number of OIDs: {}", idx_data.num_oids);
    eprintln!("Maximum sequence length: {}", idx_data.max_seq_len);

    for _ in 0..(idx_data.num_oids + 1) {
        let num_headers = get_deflines(&mut hdr_data)?;
        if num_headers == 0 {
            break;
        }
        let _ = next_sequence(&mut seq_data, &idx_data)?;
        string_clear(&mut hdr_data.fasta_hdr);
        deflines_to_header(
            &mut hdr_data.fasta_hdr,
            &hdr_data.deflines,
            _include_taxid,
            num_headers,
        );
        out_file.write_all(to_c_string(&mut hdr_data.fasta_hdr))?;
        out_file.write_all(b"\n")?;
        write_sequence(&seq_data.seq, 80, &mut out_file)?;
    }

    free_idx_data(&mut idx_data);
    free_hdr_data(&mut hdr_data);
    free_seq_data(&mut seq_data);
    eprintln!("Wrote {} sequences to {output_path}", num_oids);
    Ok(())
}

fn blast_to_fasta_usage(prog: &str) -> String {
    format!("{prog} [-hst] [-o out_file] [-w width] blast_volume")
}

fn blast_to_fasta_help(prog: &str) -> String {
    [
        blast_to_fasta_usage(prog),
        "  blast_volume: the base name of the blast volume e.g. core_nt.00".to_string(),
        "  -h: print this help message and exit".to_string(),
        "  -o: The filename that the output gets saved to (default: <blast_volume>.fna)"
            .to_string(),
        "  -s: BLAST merges the headers of FASTA entries with identical sequences, this option outputs a complete FASTA record for every such header".to_string(),
        "  -t: Prepend the tax ID to the FASTA header with format kraken:taxid|12345|"
            .to_string(),
        "  -w: The width of the FASTA sequences (default: 80)".to_string(),
    ]
    .join("\n")
}

pub fn blast_to_fasta_main(args: &[String]) -> io::Result<()> {
    let prog = args
        .first()
        .cloned()
        .unwrap_or_else(|| "blast_to_fasta".to_string());
    let mut seq_width = 80usize;
    let mut split_hdr = false;
    let mut include_taxid = false;
    let mut out_filename: Option<String> = None;
    let mut volume: Option<String> = None;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" => return Err(io::Error::new(io::ErrorKind::InvalidInput, blast_to_fasta_help(&prog))),
            "-o" => {
                i += 1;
                out_filename = Some(args.get(i).cloned().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, blast_to_fasta_usage(&prog))
                })?);
            }
            "-s" => split_hdr = true,
            "-t" => include_taxid = true,
            "-w" => {
                i += 1;
                seq_width = args
                    .get(i)
                    .and_then(|s| s.parse::<usize>().ok())
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "-w has to be at least 1"))?;
                if seq_width < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "-w has to be at least 1",
                    ));
                }
            }
            other if other.starts_with('-') => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    blast_to_fasta_usage(&prog),
                ));
            }
            other => {
                if volume.is_some() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        blast_to_fasta_usage(&prog),
                    ));
                }
                volume = Some(other.to_string());
            }
        }
        i += 1;
    }

    let volume = volume.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Missing basename of BLAST volume\n{}", blast_to_fasta_usage(&prog)),
        )
    })?;
    let out_filename = out_filename.unwrap_or_else(|| format!("{volume}.fna"));

    let idx_path = format!("{volume}.nin");
    let hdr_path = format!("{volume}.nhr");
    let seq_path = format!("{volume}.nsq");

    let mut idx_data = init_idx_data(&idx_path)?;
    read_idx_data(&mut idx_data)?;
    let max_seq_block_size = next_power_of_2(max_block_size(&idx_data.seq_arr, idx_data.num_oids));
    let mut hdr_data = init_hdr_data(&hdr_path, &mut idx_data.hdr_arr, idx_data.num_oids)?;
    let mut seq_data = init_seq_data(&seq_path, idx_data.max_seq_len, max_seq_block_size)?;

    eprintln!("The version of this database is {:X}", idx_data.fmt_version);
    eprintln!(
        "Sequence type is: {}",
        if idx_data.db_seq_type == 0 {
            "nucleotide"
        } else {
            "protein"
        }
    );
    eprintln!("Volume number: {}", idx_data.volume);
    eprintln!(
        "Title of volume: {}",
        String::from_utf8_lossy(c_string_bytes(&mut idx_data.title))
    );
    eprintln!(
        "Date created: {}",
        String::from_utf8_lossy(c_string_bytes(&mut idx_data.date))
    );
    eprintln!("Number of OIDs: {}", idx_data.num_oids);
    eprintln!("Maximum sequence length: {}", idx_data.max_seq_len);

    let mut out_file = BufWriter::new(File::create(&out_filename)?);
    for _ in 0..(idx_data.num_oids + 1) {
        let num_headers = get_deflines(&mut hdr_data)?;
        if num_headers == 0 {
            break;
        }
        let _ = next_sequence(&mut seq_data, &idx_data)?;
        if split_hdr {
            for header_idx in 0..(num_headers as usize) {
                string_clear(&mut hdr_data.fasta_hdr);
                defline_to_header(
                    &mut hdr_data.fasta_hdr,
                    &hdr_data.deflines,
                    include_taxid,
                    header_idx,
                );
                out_file.write_all(to_c_string(&mut hdr_data.fasta_hdr))?;
                out_file.write_all(b"\n")?;
                write_sequence(&seq_data.seq, seq_width, &mut out_file)?;
            }
        } else {
            string_clear(&mut hdr_data.fasta_hdr);
            deflines_to_header(
                &mut hdr_data.fasta_hdr,
                &hdr_data.deflines,
                include_taxid,
                num_headers,
            );
            out_file.write_all(to_c_string(&mut hdr_data.fasta_hdr))?;
            out_file.write_all(b"\n")?;
            write_sequence(&seq_data.seq, seq_width, &mut out_file)?;
        }
    }

    free_idx_data(&mut idx_data);
    free_hdr_data(&mut hdr_data);
    free_seq_data(&mut seq_data);
    Ok(())
}

/// Extract readable text from ASN.1 encoded header data.
/// This is a simplified extraction that finds printable ASCII strings.
fn extract_readable_header(data: &[u8]) -> String {
    // Find the longest printable ASCII run
    let mut best_start = 0;
    let mut best_len = 0;
    let mut curr_start = 0;
    let mut curr_len = 0;

    for (i, &b) in data.iter().enumerate() {
        if (0x20..0x7F).contains(&b) {
            if curr_len == 0 {
                curr_start = i;
            }
            curr_len += 1;
            if curr_len > best_len {
                best_start = curr_start;
                best_len = curr_len;
            }
        } else {
            curr_len = 0;
        }
    }

    if best_len > 0 {
        String::from_utf8_lossy(&data[best_start..best_start + best_len]).into_owned()
    } else {
        format!("sequence_{}", data.len())
    }
}

/// Decode NCBI amino acid encoding to ASCII.
fn ncbi_aa_decode(code: u8) -> u8 {
    const AA_TABLE: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    if (code as usize) < AA_TABLE.len() {
        AA_TABLE[code as usize]
    } else {
        b'X'
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn write_temp_file(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("{name}_{}.bin", std::process::id()));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    #[test]
    fn test_blast_utils_nlz_and_next_power_of_2() {
        assert_eq!(nlz(0), 32);
        assert_eq!(nlz(1), 31);
        assert_eq!(nlz(0x8000_0000), 0);
        assert_eq!(next_power_of_2(1), 1);
        assert_eq!(next_power_of_2(2), 2);
        assert_eq!(next_power_of_2(3), 4);
        assert_eq!(next_power_of_2(17), 32);
    }

    #[test]
    fn test_blast_string_helpers() {
        let mut s = init_string();
        string_append_str(&mut s, "abc");
        string_append_char(&mut s, b'_');
        string_append_int(&mut s, 42);
        assert_eq!(string_length(&s), 6);
        assert_eq!(string_data(&s), b"abc_42");

        let mut t = init_string();
        let copied = string_copy(&mut t, &s, 0, 4, 2);
        assert_eq!(copied, 2);
        assert_eq!(string_data(&t), b"42");

        let appended = string_append(&mut t, &s, 0, 3);
        assert_eq!(appended, 3);
        assert_eq!(string_data(&t), b"42abc");

        let c_string = to_c_string(&mut t);
        assert_eq!(&c_string[..5], b"42abc");
        assert_eq!(c_string[5], 0);

        string_clear(&mut t);
        assert_eq!(string_length(&t), 0);
        free_string(&mut t);
        assert_eq!(string_capacity(&t), 0);
    }

    #[test]
    fn test_blast_read_helpers_and_ambiguity_unpack() {
        let mut cursor = Cursor::new(vec![0x12, 0x34, 0x56, 0x78, b'h', b'i']);
        let mut buf = [0u8; 4];
        let count = read_into_buffer(&mut cursor, &mut buf, 1, 4).unwrap();
        assert_eq!(count, 4);
        assert_eq!(&buf, &[0x12, 0x34, 0x56, 0x78]);

        let amb32 = u32_to_amb32(0xA5BC_DEFA);
        assert_eq!(amb32.offset, 0x00BC_DEFA);
        assert_eq!(amb32.length, 0x5);
        assert_eq!(amb32.value, 0xA);

        let amb64 = u64_to_amb64(0xB123_4567_89AB_CDEF);
        assert_eq!(amb64.offset, 0x89AB_CDEF);
        assert_eq!(amb64.unused, 0x4567);
        assert_eq!(amb64.length, 0x123);
        assert_eq!(amb64.value, 0xB);
    }

    #[test]
    fn test_blast_index_helpers() {
        let path = write_temp_file(
            "blast_idx_test",
            &[
                0, 0, 0, 1, // fmt_version
                0, 0, 0, 2, // db_seq_type
                0, 0, 0, 3, // volume
                0, 0, 0, 1, b'T', // title
                0, 0, 0, 1, b'L', // lmdb
                0, 0, 0, 1, b'D', // date
                0, 0, 0, 1, // num_oids
                0, 0, 0, 0, 0, 0, 0, 7, // vol_len
                0, 0, 0, 9, // max_seq_len
                0, 0, 0, 0, 0, 0, 0, 5, // hdr arr
                0, 0, 0, 0, 0, 0, 0, 6, // seq arr
                0, 0, 0, 0, 0, 0, 0, 7, // amb arr
            ],
        );

        let mut idx = init_idx_data(path.to_str().unwrap()).unwrap();
        read_idx_data(&mut idx).unwrap();
        assert_eq!(idx.fmt_version, 1);
        assert_eq!(idx.db_seq_type, 2);
        assert_eq!(idx.volume, 3);
        assert_eq!(string_data(&idx.title), b"T");
        assert_eq!(idx.num_oids, 1);
        assert_eq!(idx.vol_len, 7);
        assert_eq!(idx.max_seq_len, 9);
        assert_eq!(idx.hdr_arr.len(), 9);
        free_idx_data(&mut idx);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_asn1_cursor_helpers() {
        let path = write_temp_file(
            "asn1_cursor_test",
            &[
                0x30, 0x80, // sequence start
                0x1A, 0x03, b'a', b'b', b'c', // visible string
                0x02, 0x02, 0x01, 0x2C, // integer 300
                0x00, 0x00, // end tag
            ],
        );
        let mut asn1 = init_asn1(path.to_str().unwrap(), &[2, 5, 4, 2], 4).unwrap();

        assert!(asn1_sequence_start(&mut asn1).unwrap());
        let mut s = VisibleString::new();
        assert_eq!(asn1_get_visible_string(&mut asn1, &mut s).unwrap(), 3);
        assert_eq!(s, b"abc");
        let mut value = 0;
        assert_eq!(asn1_get_integer(&mut asn1, &mut value).unwrap(), 1);
        assert_eq!(value, 300);
        assert!(asn1_is_end_of_sequence(&mut asn1).unwrap());
        assert_eq!(asn1.curr_pos, 11);
        assert!(asn1_indefinite_tag_end(&mut asn1).unwrap());

        free_asn1(&mut asn1);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_asn1_explicit_fields_and_date() {
        let path = write_temp_file(
            "asn1_fields_test",
            &[
                0xA0, 0x80, 0x1A, 0x03, b'x', b'y', b'z', 0x00, 0x00, // optional string
                0xA1, 0x80, 0x02, 0x01, 0x07, 0x00, 0x00, // optional int
                0xA0, 0x80, 0x1A, 0x04, b'2', b'0', b'2', b'6', 0x00, 0x00, // date string
            ],
        );
        let mut asn1 = init_asn1(path.to_str().unwrap(), &[26], 1).unwrap();

        let mut s = VisibleString::new();
        asn1_get_optional_visible_string_field(&mut asn1, 0xA0, &mut s).unwrap();
        assert_eq!(s, b"xyz");

        let mut value = 0;
        asn1_get_optional_integer_field(&mut asn1, 0xA1, &mut value).unwrap();
        assert_eq!(value, 7);

        let mut d = Date::default();
        assert_eq!(get_date(&mut asn1, &mut d).unwrap(), 1);
        assert_eq!(d.date_type, 1);
        assert_eq!(d.str_value, b"2026");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_blast_defline_readers() {
        let path = write_temp_file(
            "blast_defline_test",
            &[
                0x30, 0x80, // deflines sequence
                0x30, 0x80, // defline sequence
                0xA0, 0x80, 0x1A, 0x05, b't', b'i', b't', b'l', b'e', 0x00, 0x00, // title
                0xA1, 0x80, 0x30, 0x80, // seq-id list wrapper
                0xAB, 0x80, 0x02, 0x01, 0x2A, 0x00, 0x00, // gi seq-id
                0x00, 0x00, // seq-id list sequence end
                0x00, 0x00, // seq-id wrapper end
                0xA2, 0x80, 0x02, 0x02, 0x03, 0xE8, 0x00, 0x00, // taxid 1000
                0xA3, 0x80, 0x30, 0x80, // memberships wrapper/list
                0x02, 0x01, 0x01, 0x02, 0x01, 0x02, // memberships
                0x00, 0x00, // memberships sequence end
                0x00, 0x00, // memberships wrapper end
                0x00, 0x00, // defline sequence end
                0x00, 0x00, // deflines sequence end
            ],
        );
        let mut asn1 = init_asn1(path.to_str().unwrap(), &[56], 1).unwrap();
        let mut deflines = Vec::new();

        assert_eq!(get_blast_deflines(&mut asn1, &mut deflines).unwrap(), 1);
        assert_eq!(deflines.len(), 1);
        assert_eq!(deflines[0].title, b"title");
        assert_eq!(deflines[0].taxid, 1000);
        assert_eq!(deflines[0].memberships, vec![1, 2]);
        assert_eq!(deflines[0].seq_ids.len(), 1);
        assert_eq!(deflines[0].seq_ids[0].seq_id_type, SEQ_GI);
        assert_eq!(deflines[0].seq_ids[0].int_id, 42);

        destroy_blast_defline(&mut deflines[0]);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_blast_header_and_sequence_helpers() {
        let mut hdr_offsets = vec![];
        hdr_offsets.extend_from_slice(&0u32.to_be_bytes());
        hdr_offsets.extend_from_slice(&10u32.to_be_bytes());
        hdr_offsets.extend_from_slice(&18u32.to_be_bytes());
        let hdr_path = write_temp_file(
            "blast_hdr_blocks",
            &[
                0x30, 0x80, 0x30, 0x80, 0xA0, 0x80, 0x1A, 0x01, b'a', 0x00, 0x00, 0xA1, 0x80, 0x30,
                0x80, 0xAB, 0x80, 0x02, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xA2, 0x80,
                0x02, 0x01, 0x07, 0x00, 0x00, 0x00, 0x00,
            ],
        );
        let mut hdr = init_hdr_data(hdr_path.to_str().unwrap(), &mut hdr_offsets, 2).unwrap();
        assert_eq!(get_deflines(&mut hdr).unwrap(), 1);
        free_hdr_data(&mut hdr);
        std::fs::remove_file(hdr_path).unwrap();

        let mut header = init_string();
        let mut defline = init_blast_defline();
        defline.title = b"a".to_vec();
        defline.taxid = 7;
        defline.seq_ids.push(SeqId {
            text_id: TextSeqId {
                acc: b"ACC".to_vec(),
                ver: 2,
                ..Default::default()
            },
            seq_id_type: SEQ_GENBANK,
            ..Default::default()
        });
        deflines_to_header(&mut header, &[defline], true, 1);
        assert_eq!(c_string_bytes(&mut header), b">kraken:taxid|7|ACC.2 a");

        let mut seq_offsets = vec![];
        seq_offsets.extend_from_slice(&0u32.to_be_bytes());
        seq_offsets.extend_from_slice(&2u32.to_be_bytes());
        let mut amb_offsets = vec![];
        amb_offsets.extend_from_slice(&2u32.to_be_bytes());
        amb_offsets.extend_from_slice(&2u32.to_be_bytes());
        let seq_path = write_temp_file("blast_seq_blocks", &[0, 0b00_01_10_11, 0b00_00_00_00]);
        let mut seq = init_seq_data(seq_path.to_str().unwrap(), 8, 8).unwrap();
        let idx = BlastIdxData {
            idx_file: File::open(&seq_path).unwrap(),
            fmt_version: 0,
            db_seq_type: 0,
            volume: 0,
            title: init_string(),
            lmdb_file: init_string(),
            date: init_string(),
            num_oids: 1,
            vol_len: 0,
            max_seq_len: 8,
            hdr_arr: vec![],
            seq_arr: seq_offsets,
            amb_arr: amb_offsets,
        };
        assert_eq!(max_block_size(&idx.seq_arr, 2), 2);
        assert!(!has_ambiguous_data(2, 2));
        assert_eq!(next_sequence(&mut seq, &idx).unwrap(), 0);
        assert_eq!(string_data(&seq.seq), b"ACGT");
        assert_eq!(get_nucleotide_length(&[0, 0, 0, 2], 3), 10);

        let mut out = Vec::new();
        write_sequence(&seq.seq, 3, &mut out).unwrap();
        assert_eq!(&out, b"ACG\nT\n");

        free_seq_data(&mut seq);
        std::fs::remove_file(seq_path).unwrap();
    }

    #[test]
    fn test_extract_readable_header_and_ncbi_decode() {
        assert_eq!(extract_readable_header(b"\x00foo bar\x01z"), "foo bar");
        assert_eq!(decode_ncbi2na(&[0b00_01_10_11], 4), b"ACGT");
        assert_eq!(ncbi_aa_decode(1), b'A');
        assert_eq!(ncbi_aa_decode(255), b'X');
    }

    #[test]
    fn test_blast_to_fasta_usage_and_help() {
        assert_eq!(
            blast_to_fasta_usage("blast_to_fasta"),
            "blast_to_fasta [-hst] [-o out_file] [-w width] blast_volume"
        );
        assert!(blast_to_fasta_help("blast_to_fasta").contains("kraken:taxid|12345|"));
    }

    #[test]
    fn test_blast_to_fasta_main_requires_volume() {
        let err = blast_to_fasta_main(&["blast_to_fasta".to_string()]).unwrap_err();
        assert!(err.to_string().contains("Missing basename of BLAST volume"));
    }
}
