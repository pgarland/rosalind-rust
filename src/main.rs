// Copyright 2014 Phillip Garland
// Licensed under the GNU General Public License, Version 2, or later.


use std::os;
use std::io::File;
use std::io::BufferedReader;
use std::collections::HashMap;

#[test]
use std::num::abs;


fn main () {
    let args = os::args();   
    
    let problem = args[1].as_slice();
    match problem {
        "dna" => { dna(args[2].as_slice()) },

        "rna" => { rna(args[2].as_slice())},

        "revc" => { revc(args[2].as_slice()) },

        "iprb" => {
            if args.len() != 5 {
                panic!("Usage: rosalind iprb DOMINANT HET RECESSIVE");
            }
            let homo_d = from_str(args[2].as_slice()).unwrap();
            let het = from_str(args[3].as_slice()).unwrap();
            let homo_r = from_str(args[4].as_slice()).unwrap();
            do_iprb(homo_d, het, homo_r);
        }

        "fib" => {
            let n = from_str(args[2].as_slice()).unwrap();
            let k = from_str(args[3].as_slice()).unwrap();
            println!("{}", fib(n, k));
        }

        "gc" => {
            let (max_id, max_gc) = gc(args[2].as_slice());
            println!("{}\n{}", max_id, max_gc);
        }

        "prot" => {
            let codon_file = args[2].as_slice();
            let rna = read_nuc_file(args[3].as_slice());
            let codon_table = build_codon_table(codon_file);
            println!("{}", prot(rna.as_slice(), codon_table));
        }

        "subs" => {
            let genome = args[2].as_slice();
            let probe = args[3].as_slice();
            let positions = subs(genome, probe);
            match positions {
                Some(pos_list) => print!("{}", pos_list),
                None => print!("{}", "No matches found"),
            }
        }

        "hamm" => {
            let s = args[2].as_slice();
            let t = args[3].as_slice();
            let distance = hamming_distance(s,t);
            println!("{}", distance);
        }
        _ => { panic!("{} is not a valid problem.", problem) },
    }
}


fn read_nuc_file(fname: &str) -> String {
    let path = Path::new(fname);
    let nuc_result = File::open(&path).read_to_string();

    match nuc_result {
        Ok(nucs) => {return nucs},
        Err(e) => panic!("Opening and reading failed: {}", e),
    }

}


fn dna(fname: &str) -> () {
    let nucs = read_nuc_file(fname);
    let nuc_count = count_nucs(nucs.as_slice());
    
    match nuc_count {
        NucCount {a, c, g, t} => println!("{} {} {} {}", a, c, g, t)
    };
}


struct NucCount {
    a: u32,
    c: u32,
    g: u32,
    t: u32,
}


// Return a struct containing the number of occurances of each
// nucleotide in the nucs string.
fn count_nucs(nucs: &str) -> NucCount {
    let mut nuc_count = NucCount { a: 0, t: 0, c: 0, g: 0 };
            
    for nuc in nucs.chars() {
        match nuc {
            'A' => {nuc_count.a += 1}
            'C' => {nuc_count.c += 1}
            'G' => {nuc_count.g += 1}
            'T' => {nuc_count.t += 1}
            nuc @ _ => panic!("Unhandled nucleotide!: {}", nuc)
        }
    }

    return nuc_count;
}


fn rna(fname: &str) -> () {
    let dna = read_nuc_file(fname);
    let rna = transcribe(dna.as_slice());

    println!("{}", rna);
}


// Given a DNA string, return the expected RNA string when the DNA
// sequence is transcribed in a cell.
fn transcribe (dna: &str) -> String {
    let mut rna = String::new();
    
    for nuc in dna.chars() {
        // FIXME: learn a way to factor out the common functionality in transcribe/complement/revcomp
        match nuc  {
            'A' => {rna.push_str("A")}
            'C' => {rna.push_str("C")}
            'G' => {rna.push_str("G")}
            'T' => {rna.push_str("U")}
            nuc @ _ => panic!("Unhandled nucleotide!: {}", nuc)
        }
    }

    rna
}


// Given a DNA sequence, return the sequence on the opposite DNA
// strand
fn revcomp(dna: &str) -> String {
    let it = dna.chars();

    let mut cdna = String::new();

    for base in it.rev() {
        match base  {
            'A' => {cdna.push_str("T")}
            'C' => {cdna.push_str("G")}
            'G' => {cdna.push_str("C")}
            'T' => {cdna.push_str("A")}
            base @ _ => panic!("Unhandled nucleotide!: {}", base)
        }
    }

    cdna
}


fn revc(dna: &str) -> () {
    let cdna = revcomp(dna);
    println!("{}", cdna);
}


fn iprb(homo_d: uint, het: uint, homo_r: uint) -> f32 {
    let f_homo_d = homo_d as f32;
    let f_het = het as f32;
    let f_homo_r = homo_r as f32;
    let indv: f32  = f_homo_d + f_het + f_homo_r;

    let mut prob = 0.0;
    
    // (prob choosing each genotype) * (prob choosing the genotype,
    // with the first individual removed) * (prob of that pairing
    // resulting in a dominant phenotype
    prob += (f_homo_d/indv) * ((f_homo_d-1.0)/(indv-1.0)) * 1.0;
    prob += (f_homo_d/indv) * (f_het/(indv-1.0)) * 1.0;
    prob += (f_homo_d/indv) * (f_homo_r/(indv-1.0)) * 1.0;
    prob += (f_het/indv) * (f_homo_d/(indv-1.0));
    prob += (f_het/indv) * ((f_het-1.0)/(indv-1.0)) * 0.75;
    prob += (f_het/indv) * (f_homo_r/(indv-1.0)) * 0.5;
    prob += (f_homo_r/indv) * (f_homo_d/(indv-1.0)) * 1.0;
    prob += (f_homo_r/indv) * (f_het/(indv-1.0)) * 0.5;
    prob += (f_homo_r/indv) * ((f_homo_r-1.0)/(indv-1.0)) * 0.0;

    prob
}


fn do_iprb(homo_d: uint, het: uint, homo_r: uint) -> () {
    let prob = iprb(homo_d, het, homo_r);
    println!("{}", prob);
}


#[test]
fn irpb_test() {
    let prob = iprb(2, 2, 2);
    let epsilon: f32 = 5e-6;
    assert!( abs(prob - 0.78333) < epsilon);
}


fn fib( n: uint, k: uint) -> uint {
    if n == 1 || n == 2 { 1 }
    else { k* fib(n-2, k) + fib(n-1, k) } 
}


#[test]
fn fib_test() {
    let rabbits = fib(5, 3);
    assert!(rabbits == 19);
}

// Read in a FASTA file and return a hashmap where each key is an id
// and each corresponding value is the molecular sequence for that id
fn read_fasta(fname: &str) -> HashMap<String, String> {
    let path = Path::new(fname);
    let mut file = BufferedReader::new(File::open(&path));

    let mut entries = HashMap::new();
    let mut id = String::new();
    let mut seq = String::new();
    for line in file.lines() {
        let l = line.clone().unwrap();
        if !is_id_line(l.as_slice()) {
            if id.as_slice() == "" { panic!("Invalid FASTA file: sequence does not have id") }
            else {seq.push_str(l.as_slice().trim()) }
        }

        if is_id_line(l.as_slice()) {
            if id.as_slice() != "" {
                entries.insert(id, seq.clone());
                seq.clear();
            }
            id = extract_id(l.as_slice().trim());
            println!("{}", id)
        }
    }
    // push the last sequence
    entries.insert(id, seq.clone());
    
    entries
}


fn is_id_line(line: &str) -> bool {
    line.starts_with(">")
}

    
fn extract_id(line: &str) -> String {
    if !is_id_line(line) {panic!("{} is not a valid FASTA id line", line)}

    line.trim_left_chars('>').to_string()
}


// Given a DNA string, return the percentage of bases that are either
// guanine or cytosine
fn gc_content(dna: &str) -> f32 {
    let nuc_counts = count_nucs(dna);
    let len = dna.len() as f32;
    let gc = (nuc_counts.g + nuc_counts.c) as f32;

    gc/len
}


fn gc(fname: &str) -> (String, f32) {
    
    let mut max_id = String::new();
    let mut max_gc = 0.0;

    let seqs = read_fasta(fname);
    let mut this_gc = 0.0;
    for (id, seq) in seqs.iter() {
        this_gc = gc_content(seq.as_slice());

        if this_gc >= max_gc {
            max_gc = this_gc;
            max_id = id.clone();
        }
    }

    (max_id, max_gc * 100.0)
}


#[test]
fn gc_test() {
    let test_file = "data/gc.fasta";
    let epsilon = 1e-6;

    let (id, gc) = gc(test_file);
    assert!(id.as_slice() == "Rosalind_0808");
    assert!( abs(gc - 60.919540) < epsilon);
}


// Given a file where each line is mapping between RNA codons and
// amino acids, return a hash map where each key is an RNA codon and
// each value is the corresponding amino acid
fn build_codon_table(fname: &str) -> HashMap<String, String> {
    let path = Path::new(fname);
    let mut file = BufferedReader::new(File::open(&path));

    let mut codon_table = HashMap::new();
   
    for line in file.lines() {
        let l = line.clone().unwrap();
        let pair: Vec<&str> = l.as_slice().trim_right().split(' ').collect();
        codon_table.insert(pair[0].to_string(), pair[1].to_string());
    }

    return codon_table
}


// Given an RNA sequence and a codon table return the peptide expected
// when the RNA is translated.
fn translate <'a> (rna: &'a str, codon_table: HashMap<String, String>) -> String {
    let mut peptide = String::new();

    let mut it = range(0, rna.len()).filter(|idx| *idx % 3 == 0);
    for idx in it {

        let codon = rna.slice(idx, idx + 3).to_string();

        let aa = codon_table.get(&codon).unwrap();
        if aa.as_slice() == "(Stop)" {
            return peptide
        }

        peptide.push_str((*aa).as_slice());
    }

    peptide
}


fn prot(rna: &str, codon_table: HashMap<String, String>) -> String {
    translate(rna, codon_table)
}


#[test]
fn prot_test() {
    let rna = read_nuc_file("data/prot.fasta");
    let codon_table = build_codon_table("data/aa.txt");
    assert!(prot(rna.as_slice(), codon_table).as_slice() == "MAMAPRTEINSTRING");
}


// Return a hashmap of all substrings of genome of length kmer_len
fn hash_genome(genome: &str, kmer_len: uint) -> HashMap<String, Vec<uint>> {

    let mut genome_hash : HashMap<String, Vec<uint>> = HashMap::new();

    for pos in range(0, genome.len() - kmer_len) {
        let kmer = genome.slice(pos, pos + kmer_len).to_string();
  
        if !genome_hash.contains_key(&kmer) {
            genome_hash.insert(kmer, vec![pos + 1]);
        } else {
            let positions = genome_hash.get_mut(&kmer).unwrap();
            (*positions).push(pos + 1);
        }
    }
    genome_hash
}


// Return the starting position for all matches of probe in genome
fn subs(genome: &str, probe: &str) -> Option<Vec<uint>> {
    let gh = hash_genome(genome, probe.len());

    gh.find_copy(&probe.to_string())
}


#[test]
fn subs_test() {
    let genome = "GATATATGCATATACTT";
    let probe = "ATAT";
    assert_eq!(subs(genome, probe).unwrap(), vec![2, 4, 10]);
}


// return the number of positions at which the corresponding
// chararcters in s and t differ
fn hamming_distance(s: &str, t: &str) -> uint {
    let mut distance = 0u;

    // iterate over both strings simultaneously
    let mut it = s.chars().zip(t.chars());
    for (x, y) in it {
        if x != y { distance += 1 }
    }
    distance
}


#[test]
fn hamming_test() {
    let s = "GAGCCTACTAACGGGAT";
    let t = "CATCGTAATGACGGCCT";

    assert_eq!(7, hamming(s,t));
}

