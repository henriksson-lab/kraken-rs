use std::collections::{BTreeMap, BTreeSet, VecDeque};
use ahash::AHashMap as HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use crate::mmap_file::MMapFile;
use crate::types::TaxonomyNode;

/// Magic header for the binary taxonomy file format.
/// C++ writes strlen("K2TAXDAT") = 8 bytes (no null terminator on disk).
const FILE_MAGIC: &[u8] = b"K2TAXDAT";

/// NCBI taxonomy — parses nodes.dmp and names.dmp, converts to Kraken 2 format.
/// Uses BTreeMap/BTreeSet to match the iteration order of C++ std::map/std::set.
pub struct NCBITaxonomy {
    parent_map: BTreeMap<u64, u64>,
    name_map: BTreeMap<u64, String>,
    rank_map: BTreeMap<u64, String>,
    child_map: BTreeMap<u64, BTreeSet<u64>>,
    marked_nodes: BTreeSet<u64>,
    known_ranks: BTreeSet<String>,
}

impl NCBITaxonomy {
    /// Parse NCBI nodes.dmp and names.dmp files.
    /// Exact port of C++ `NCBITaxonomy::NCBITaxonomy()`.
    pub fn new(nodes_filename: &str, names_filename: &str) -> io::Result<Self> {
        let mut tax = NCBITaxonomy {
            parent_map: BTreeMap::new(),
            name_map: BTreeMap::new(),
            rank_map: BTreeMap::new(),
            child_map: BTreeMap::new(),
            marked_nodes: BTreeSet::new(),
            known_ranks: BTreeSet::new(),
        };

        let delim = "\t|\t";

        // Parse nodes.dmp
        let nodes_file = BufReader::new(File::open(nodes_filename)?);
        for line in nodes_file.lines() {
            let mut line = line?;
            // Discard trailing "\t|"
            if line.ends_with("\t|") {
                line.truncate(line.len() - 2);
            }

            let mut node_id: u64 = 0;
            let mut parent_id: u64 = 0;
            let mut rank = String::new();

            let mut pos1 = 0;
            let mut field_ct = 0;

            loop {
                field_ct += 1;
                if field_ct > 10 {
                    break;
                }

                let (token, next_pos) = if let Some(pos2) = line[pos1..].find(delim) {
                    (line[pos1..pos1 + pos2].to_string(), Some(pos1 + pos2 + delim.len()))
                } else {
                    (line[pos1..].to_string(), None)
                };

                match field_ct {
                    1 => {
                        node_id = token.parse::<u64>().unwrap_or(0);
                        assert!(node_id != 0, "attempt to create taxonomy w/ node ID == 0");
                    }
                    2 => {
                        parent_id = token.parse::<u64>().unwrap_or(0);
                    }
                    3 => {
                        rank = token;
                        break; // finished = true
                    }
                    _ => {}
                }

                match next_pos {
                    Some(np) => pos1 = np,
                    None => break,
                }
            }

            if node_id == 1 {
                parent_id = 0;
            }
            tax.parent_map.insert(node_id, parent_id);
            tax.child_map.entry(parent_id).or_default().insert(node_id);
            tax.known_ranks.insert(rank.clone());
            tax.rank_map.insert(node_id, rank);
        }

        // Parse names.dmp
        let names_file = BufReader::new(File::open(names_filename)?);
        for line in names_file.lines() {
            let mut line = line?;
            if line.ends_with("\t|") {
                line.truncate(line.len() - 2);
            }

            let mut node_id: u64 = 0;
            let mut name = String::new();

            let mut pos1 = 0;
            let mut field_ct = 0;

            loop {
                field_ct += 1;
                if field_ct > 10 {
                    break;
                }

                let (token, next_pos) = if let Some(pos2) = line[pos1..].find(delim) {
                    (line[pos1..pos1 + pos2].to_string(), Some(pos1 + pos2 + delim.len()))
                } else {
                    (line[pos1..].to_string(), None)
                };

                match field_ct {
                    1 => {
                        node_id = token.parse::<u64>().unwrap_or(0);
                        assert!(node_id != 0, "attempt to create taxonomy w/ node ID == 0");
                    }
                    2 => {
                        name = token;
                    }
                    4 => {
                        if token == "scientific name" {
                            tax.name_map.insert(node_id, name.clone());
                        }
                        break; // finished = true
                    }
                    _ => {}
                }

                match next_pos {
                    Some(np) => pos1 = np,
                    None => break,
                }
            }
        }

        tax.marked_nodes.insert(1); // mark root node
        Ok(tax)
    }

    /// Mark a node and all its unmarked ancestors.
    /// Exact port of C++ `NCBITaxonomy::MarkNode()`.
    pub fn mark_node(&mut self, mut taxid: u64) {
        while !self.marked_nodes.contains(&taxid) {
            self.marked_nodes.insert(taxid);
            taxid = *self.parent_map.get(&taxid).unwrap_or(&0);
        }
    }

    /// Convert the marked NCBI taxonomy to Kraken 2 binary format and write to disk.
    /// Exact port of C++ `NCBITaxonomy::ConvertToKrakenTaxonomy()`.
    pub fn convert_to_kraken_taxonomy(&self, filename: &str) -> io::Result<()> {
        let node_count = self.marked_nodes.len() + 1; // +1 because 0 is illegal value
        let mut nodes = vec![TaxonomyNode::default(); node_count];

        let mut name_data = Vec::<u8>::new();
        let mut rank_data = Vec::<u8>::new();

        // Build rank offset map — deduplicated rank strings
        let mut rank_offsets: BTreeMap<String, u64> = BTreeMap::new();
        for rank in &self.known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.extend_from_slice(rank.as_bytes());
            rank_data.push(0); // null terminator
        }

        let mut internal_node_id: u64 = 0;
        let mut external_id_map: BTreeMap<u64, u64> = BTreeMap::new();
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1);

        // BFS through NCBI taxonomy
        let mut bfs_queue: VecDeque<u64> = VecDeque::new();
        bfs_queue.push_back(1);

        while let Some(external_node_id) = bfs_queue.pop_front() {
            internal_node_id += 1;
            external_id_map.insert(external_node_id, internal_node_id);

            let parent_ext = self.parent_map.get(&external_node_id).copied().unwrap_or(0);
            let parent_int = external_id_map.get(&parent_ext).copied().unwrap_or(0);

            let rank = self.rank_map.get(&external_node_id).map(|s| s.as_str()).unwrap_or("");
            let rank_offset = rank_offsets.get(rank).copied().unwrap_or(0);

            let name_offset = name_data.len() as u64;

            let first_child = internal_node_id + 1 + bfs_queue.len() as u64;
            let mut child_count: u64 = 0;

            if let Some(children) = self.child_map.get(&external_node_id) {
                for &child_node in children {
                    if self.marked_nodes.contains(&child_node) {
                        bfs_queue.push_back(child_node);
                        child_count += 1;
                    }
                }
            }

            nodes[internal_node_id as usize] = TaxonomyNode {
                parent_id: parent_int,
                first_child,
                child_count,
                name_offset,
                rank_offset,
                external_id: external_node_id,
                godparent_id: 0,
            };

            let name = self.name_map.get(&external_node_id).map(|s| s.as_str()).unwrap_or("");
            name_data.extend_from_slice(name.as_bytes());
            name_data.push(0); // null terminator
        }

        // Write to disk
        let taxonomy = Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map: HashMap::new(),
            _mmap: None,
        };
        taxonomy.write_to_disk(filename)
    }
}

/// Kraken 2 taxonomy — loaded from binary file, supports LCA queries.
pub struct Taxonomy {
    nodes: Vec<TaxonomyNode>,
    name_data: Vec<u8>,
    rank_data: Vec<u8>,
    external_to_internal_id_map: HashMap<u64, u64>,
    _mmap: Option<MMapFile>,
}

impl Taxonomy {
    /// Load taxonomy from a binary file.
    /// Exact port of C++ `Taxonomy::Init()`.
    pub fn from_file(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        if memory_mapping {
            Self::from_file_mmap(filename)
        } else {
            Self::from_file_read(filename)
        }
    }

    fn from_file_read(filename: &str) -> io::Result<Self> {
        let mut file = File::open(filename)?;

        // Read and verify magic
        let mut magic = vec![0u8; FILE_MAGIC.len()];
        file.read_exact(&mut magic)?;
        if magic != FILE_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("malformed taxonomy file {}", filename)));
        }

        // Read header
        let mut buf8 = [0u8; 8];
        file.read_exact(&mut buf8)?;
        let node_count = usize::from_le_bytes(buf8);
        file.read_exact(&mut buf8)?;
        let name_data_len = usize::from_le_bytes(buf8);
        file.read_exact(&mut buf8)?;
        let rank_data_len = usize::from_le_bytes(buf8);

        // Read nodes
        let mut node_bytes = vec![0u8; node_count * std::mem::size_of::<TaxonomyNode>()];
        file.read_exact(&mut node_bytes)?;
        let nodes = unsafe {
            let mut nodes = vec![TaxonomyNode::default(); node_count];
            std::ptr::copy_nonoverlapping(
                node_bytes.as_ptr(),
                nodes.as_mut_ptr() as *mut u8,
                node_bytes.len(),
            );
            nodes
        };

        // Read name and rank data
        let mut name_data = vec![0u8; name_data_len];
        file.read_exact(&mut name_data)?;
        let mut rank_data = vec![0u8; rank_data_len];
        file.read_exact(&mut rank_data)?;

        Ok(Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map: HashMap::new(),
            _mmap: None,
        })
    }

    fn from_file_mmap(filename: &str) -> io::Result<Self> {
        let mmap = MMapFile::open_read_only(filename)?;
        let data = mmap.as_slice();

        if data.len() < FILE_MAGIC.len() || &data[..FILE_MAGIC.len()] != FILE_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("malformed taxonomy file {}", filename)));
        }

        let mut offset = FILE_MAGIC.len();

        let node_count = usize::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        offset += 8;
        let name_data_len = usize::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        offset += 8;
        let rank_data_len = usize::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        offset += 8;

        let node_size = std::mem::size_of::<TaxonomyNode>();
        let nodes_bytes = &data[offset..offset + node_count * node_size];
        let nodes = unsafe {
            let mut nodes = vec![TaxonomyNode::default(); node_count];
            std::ptr::copy_nonoverlapping(
                nodes_bytes.as_ptr(),
                nodes.as_mut_ptr() as *mut u8,
                nodes_bytes.len(),
            );
            nodes
        };
        offset += node_count * node_size;

        let name_data = data[offset..offset + name_data_len].to_vec();
        offset += name_data_len;
        let rank_data = data[offset..offset + rank_data_len].to_vec();

        Ok(Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map: HashMap::new(),
            _mmap: Some(mmap),
        })
    }

    /// Create an empty taxonomy.
    pub fn new() -> Self {
        Taxonomy {
            nodes: Vec::new(),
            name_data: Vec::new(),
            rank_data: Vec::new(),
            external_to_internal_id_map: HashMap::new(),
            _mmap: None,
        }
    }

    /// Write taxonomy to binary file.
    /// Exact port of C++ `Taxonomy::WriteToDisk()`.
    pub fn write_to_disk(&self, filename: &str) -> io::Result<()> {
        let mut file = File::create(filename)?;

        file.write_all(FILE_MAGIC)?;

        let node_count = self.nodes.len() as u64;
        file.write_all(&node_count.to_le_bytes())?;
        let name_data_len = self.name_data.len() as u64;
        file.write_all(&name_data_len.to_le_bytes())?;
        let rank_data_len = self.rank_data.len() as u64;
        file.write_all(&rank_data_len.to_le_bytes())?;

        // Write nodes as raw bytes
        let node_bytes = unsafe {
            std::slice::from_raw_parts(
                self.nodes.as_ptr() as *const u8,
                self.nodes.len() * std::mem::size_of::<TaxonomyNode>(),
            )
        };
        file.write_all(node_bytes)?;

        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        Ok(())
    }

    /// Build the external-to-internal ID map.
    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);
        for i in 1..self.nodes.len() {
            self.external_to_internal_id_map.insert(self.nodes[i].external_id, i as u64);
        }
    }

    /// Look up an internal ID from an external (NCBI) ID.
    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        self.external_to_internal_id_map.get(&external_id).copied().unwrap_or(0)
    }

    /// Check if node A is an ancestor of node B (using internal IDs).
    /// Exact port of C++ `Taxonomy::IsAAncestorOfB()`.
    pub fn is_a_ancestor_of_b(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }
        let mut b = b;
        while b > a {
            b = self.nodes[b as usize].parent_id;
        }
        b == a
    }

    /// Compute the Lowest Common Ancestor of two nodes (using internal IDs).
    /// Exact port of C++ `Taxonomy::LowestCommonAncestor()`.
    pub fn lowest_common_ancestor(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }
        let mut a = a;
        let mut b = b;
        while a != b {
            if a > b {
                a = self.nodes[a as usize].parent_id;
            } else {
                b = self.nodes[b as usize].parent_id;
            }
        }
        a
    }

    /// Get the number of nodes (including the reserved node 0).
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Ensure taxonomy data is in heap memory (not mmap-backed).
    /// Port of C++ `Taxonomy::MoveToMemory()`.
    /// No-op if already heap-backed (Rust always copies to Vec on load).
    pub fn move_to_memory(&mut self) {
        // In the Rust implementation, data is always owned in Vecs,
        // so this is a no-op. The _mmap field is dropped when replaced.
        self._mmap = None;
    }

    /// Access node by internal ID.
    pub fn node(&self, internal_id: u64) -> &TaxonomyNode {
        &self.nodes[internal_id as usize]
    }

    /// Access the nodes slice.
    pub fn nodes(&self) -> &[TaxonomyNode] {
        &self.nodes
    }

    /// Get a null-terminated name string starting at the given offset.
    pub fn name_at_offset(&self, offset: u64) -> &str {
        let start = offset as usize;
        if start >= self.name_data.len() {
            return "";
        }
        let end = self.name_data[start..].iter().position(|&b| b == 0)
            .map(|p| p + start)
            .unwrap_or(self.name_data.len());
        std::str::from_utf8(&self.name_data[start..end]).unwrap_or("")
    }

    /// Get a null-terminated rank string starting at the given offset.
    pub fn rank_at_offset(&self, offset: u64) -> &str {
        let start = offset as usize;
        if start >= self.rank_data.len() {
            return "";
        }
        let end = self.rank_data[start..].iter().position(|&b| b == 0)
            .map(|p| p + start)
            .unwrap_or(self.rank_data.len());
        std::str::from_utf8(&self.rank_data[start..end]).unwrap_or("")
    }

    /// Get the raw name data bytes.
    pub fn name_data(&self) -> &[u8] {
        &self.name_data
    }

    /// Get the raw rank data bytes.
    pub fn rank_data(&self) -> &[u8] {
        &self.rank_data
    }
}

impl Default for Taxonomy {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn test_data_dir() -> String {
        format!("{}/kraken2/data", env!("CARGO_MANIFEST_DIR"))
    }

    /// Helper: create a seqid2taxid map from test FASTA headers.
    fn create_test_seqid2taxid(path: &str) -> io::Result<()> {
        // The test FASTA files have headers like ">NC_045512.2 ..." but no explicit taxid.
        // We need to create a map. Look at names.dmp for available taxids.
        // For the test, we just mark a few known taxids from names.dmp.
        let mut file = File::create(path)?;
        // These are taxids present in the test data names.dmp/nodes.dmp:
        // 2697049 = SARS-CoV-2, 694009 = SARS, 1335626 = MERS, etc.
        writeln!(file, "seq1\t2697049")?;
        writeln!(file, "seq2\t694009")?;
        writeln!(file, "seq3\t1335626")?;
        writeln!(file, "seq4\t11320")?; // Influenza A
        writeln!(file, "seq5\t11520")?; // Influenza B
        writeln!(file, "seq6\t11676")?; // HIV-1
        writeln!(file, "seq7\t11709")?; // HIV-2
        writeln!(file, "seq8\t10710")?; // Lambda phage
        Ok(())
    }

    #[test]
    fn test_taxonomy_roundtrip_matches_reference() {
        // Compare against C++-generated reference file
        let ref_path = format!("{}/tests/reference/taxo.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            eprintln!("Skipping: reference data not available");
            return;
        }

        let data_dir = test_data_dir();
        let nodes_file = format!("{}/nodes.dmp", data_dir);
        let names_file = format!("{}/names.dmp", data_dir);

        let tmp_dir = tempfile::tempdir().unwrap();
        let seqid2taxid = tmp_dir.path().join("seqid2taxid.map");
        let rust_taxo_file = tmp_dir.path().join("rust_taxo.k2d");

        // Create same seqid2taxid map used for reference generation
        let ref_map = format!("{}/tests/reference/seqid2taxid.map", env!("CARGO_MANIFEST_DIR"));
        std::fs::copy(&ref_map, &seqid2taxid).unwrap();

        // Generate via Rust
        let mut ncbi_tax = NCBITaxonomy::new(&nodes_file, &names_file).unwrap();
        let map_file = BufReader::new(File::open(&seqid2taxid).unwrap());
        for line in map_file.lines() {
            let line = line.unwrap();
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let Ok(taxid) = parts[1].parse::<u64>() {
                    if taxid != 0 {
                        ncbi_tax.mark_node(taxid);
                    }
                }
            }
        }
        ncbi_tax.convert_to_kraken_taxonomy(rust_taxo_file.to_str().unwrap()).unwrap();

        // Compare against C++ reference
        let ref_bytes = std::fs::read(&ref_path).unwrap();
        let rust_bytes = std::fs::read(&rust_taxo_file).unwrap();
        assert_eq!(ref_bytes.len(), rust_bytes.len(),
            "Taxonomy size mismatch: ref={}, rust={}", ref_bytes.len(), rust_bytes.len());
        assert_eq!(ref_bytes, rust_bytes, "Taxonomy files differ from C++ reference");
    }

    #[test]
    fn test_taxonomy_lca() {
        // Load reference taxonomy and test LCA properties
        let ref_path = format!("{}/tests/reference/taxo.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            eprintln!("Skipping: reference data not available");
            return;
        }

        let mut tax = Taxonomy::from_file(&ref_path, false).unwrap();
        tax.generate_external_to_internal_id_map();

        // LCA(x, 0) = x, LCA(0, x) = x
        assert_eq!(tax.lowest_common_ancestor(1, 0), 1);
        assert_eq!(tax.lowest_common_ancestor(0, 1), 1);
        assert_eq!(tax.lowest_common_ancestor(0, 0), 0);

        // LCA(x, x) = x
        for i in 1..tax.node_count().min(20) as u64 {
            assert_eq!(tax.lowest_common_ancestor(i, i), i);
        }

        // LCA(root, anything) = root
        for i in 1..tax.node_count().min(20) as u64 {
            assert_eq!(tax.lowest_common_ancestor(1, i), 1);
        }

        // LCA(child, parent) = parent
        for i in 2..tax.node_count().min(20) as u64 {
            let parent = tax.node(i).parent_id;
            assert_eq!(tax.lowest_common_ancestor(i, parent), parent);
        }
    }
}
