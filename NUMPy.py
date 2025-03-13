#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NUMPy: NCBI, UniProt & MUSCLE Phylogenetic Tree Builder
with organism name appended to accession labels in the Newick file.
Spaces in organism names are replaced with underscores before saving to Excel.

Workflow:
  1) Load a protein sequence from a FASTA file.
  2) Perform BLAST searches (NCBI nr and UniProt via EBI).
  3) Save or reuse the BLAST XML results.
  4) Parse each "Hit Title" to extract the organism name (using either square brackets or the OS=...OX= pattern),
     then replace spaces with underscores. The extracted name is stored in the "Organism" field.
  5) Write an Excel file that includes an "Organism" column.
  6) Build a unique FASTA file for tree construction. For each hit, the label is built as:
         ACCESSION-ORGANISM
     (if an organism name is available) with colons replaced by underscores.
  7) The query sequence is added with a distinct label ("QUERY_" + simple accession).
  8) Force regeneration of the Newick file so that the updated labels are used.
  9) Build a phylogenetic tree from the alignment and save the Newick file.
 10) Render a PDF tree via ETE3 with the query leaf highlighted in red and thicker branch lines.

All output files use the base name of the input FASTA as a prefix.

Usage:
    python NUMPy.py <fasta_file>

Dependencies:
    - Biopython
    - pandas
    - requests
    - matplotlib
    - ETE3 (install via pip)
    - MUSCLE 5.3 installed and accessible as "muscle" (with options: -align, -output)
"""

import sys
import os
import time
import re
from io import StringIO
import subprocess
import requests
import pandas as pd
import xml.etree.ElementTree as ET

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')


def extract_organism_from_title(title: str) -> str:
    """
    Extract the organism's scientific name from a BLAST 'Hit Title'.

    Strategies:
      1) Look for text in square brackets, e.g. [Vitis vinifera]
      2) Otherwise, check for text between "OS=" and "OX=".
      
    Returns the extracted name or an empty string if not found.
    """
    if not isinstance(title, str):
        return ""
    bracket_match = re.search(r'\[([^\]]+)\]', title)
    if bracket_match:
        return bracket_match.group(1).strip()
    if "OS=" in title and "OX=" in title:
        try:
            return title.split("OS=")[1].split("OX=")[0].strip()
        except Exception:
            pass
    return ""


class ProteinBlastSearch:
    """
    Perform BLAST searches, parse results, and store them.
    
    Attributes:
      fasta_file: Path to the FASTA file.
      email: Email address for BLAST requests.
      query_sequence: Protein sequence read from FASTA.
      query_length: Length of the protein sequence.
      query_id: Identifier from the FASTA header.
      results: List of dictionaries with parsed BLAST results.
      prefix: Base name (without extension) of the FASTA file.
    """
    def __init__(self, fasta_file, email):
        self.fasta_file = fasta_file
        self.email = email
        self.query_sequence = ""
        self.query_length = 0
        self.query_id = ""
        self.results = []
        self.prefix = os.path.splitext(os.path.basename(self.fasta_file))[0]
    
    def load_sequence(self):
        records = list(SeqIO.parse(self.fasta_file, "fasta"))
        if not records:
            logging.error("No sequences found in the provided FASTA file.")
            sys.exit(1)
        self.query_sequence = str(records[0].seq)
        self.query_length = len(self.query_sequence)
        self.query_id = records[0].id
        logging.info(f"Loaded protein sequence: {self.query_id}")
    
    def run_ncbi_blast(self):
        xml_filename = f"{self.prefix}_ncbi_blast_result.xml"
        if os.path.exists(xml_filename):
            logging.info(f"NCBI BLAST result file {xml_filename} found. Reusing it...")
            with open(xml_filename, "r") as f:
                result_xml = f.read()
        else:
            logging.info("Running BLAST search on NCBI (nr database)...")
            try:
                result_handle = NCBIWWW.qblast("blastp", "nr", self.query_sequence)
            except Exception as e:
                logging.error(f"Error running NCBI BLAST: {e}")
                return []
            result_xml = result_handle.read()
            with open(xml_filename, "w") as f:
                f.write(result_xml)
        try:
            blast_record = NCBIXML.read(StringIO(result_xml))
        except Exception as e:
            logging.error(f"Error parsing NCBI BLAST XML: {e}")
            return []
        ncbi_results = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                coverage = round(100 * hsp.align_length / self.query_length, 2)
                ncbi_results.append({
                    "Database": "NCBI nr",
                    "Hit Title": alignment.title,
                    "Hit ID": "",
                    "Hit Length": alignment.length,
                    "E-value": hsp.expect,
                    "Score": hsp.score,
                    "Query Start": hsp.query_start,
                    "Query End": hsp.query_end,
                    "Hit Sequence": hsp.sbjct,
                    "Coverage (%)": coverage
                })
        return ncbi_results
    
    def run_uniprot_blast(self):
        xml_filename = f"{self.prefix}_uniprot_blast_result.xml"
        if os.path.exists(xml_filename):
            logging.info(f"UniProt BLAST result file {xml_filename} found. Reusing it...")
            with open(xml_filename, "r") as f:
                result_xml = f.read()
        else:
            logging.info("Running BLAST search on UniProt (via EBI BLAST service)...")
            url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
            params = {
                "email": self.email,
                "program": "blastp",
                "database": "uniprotkb",
                "sequence": self.query_sequence,
                "stype": "protein"
            }
            try:
                response = requests.post(url, data=params)
            except Exception as e:
                logging.error(f"Error sending UniProt BLAST request: {e}")
                return []
            logging.info(f"UniProt BLAST response status: {response.status_code}")
            logging.info(f"Response content: {response.text}")
            if response.status_code != 200:
                logging.error("Error: UniProt BLAST search request failed")
                return []
            job_id = response.text.strip()
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
            start_time = time.time()
            while True:
                try:
                    status_response = requests.get(status_url)
                except Exception as e:
                    logging.error(f"Error polling UniProt BLAST status: {e}")
                    return []
                if status_response.text.strip() == "FINISHED":
                    print()
                    break
                elapsed = time.time() - start_time
                minutes = int(elapsed // 60)
                seconds = int(elapsed % 60)
                print(f"Waiting for UniProt BLAST job to finish... Elapsed time: {minutes:02d}:{seconds:02d}",
                      end="\r", flush=True)
                time.sleep(5)
            result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/xml"
            try:
                result_response = requests.get(result_url)
            except Exception as e:
                logging.error(f"Error retrieving UniProt BLAST results: {e}")
                return []
            if result_response.status_code != 200:
                logging.error("Error retrieving UniProt BLAST results")
                return []
            result_xml = result_response.text
            with open(xml_filename, "w") as f:
                f.write(result_xml)
        try:
            blast_record = NCBIXML.read(StringIO(result_xml))
            uni_results = []
            for alignment in blast_record.alignments:
                hit_id = getattr(alignment, "hit_id", "")
                for hsp in alignment.hsps:
                    coverage = round(100 * hsp.align_length / self.query_length, 2)
                    uni_results.append({
                        "Database": "UniProt (EBI BLAST)",
                        "Hit Title": alignment.title,
                        "Hit ID": hit_id,
                        "Hit Length": alignment.length,
                        "E-value": hsp.expect,
                        "Score": hsp.score,
                        "Query Start": hsp.query_start,
                        "Query End": hsp.query_end,
                        "Hit Sequence": hsp.sbjct,
                        "Coverage (%)": coverage
                    })
            return uni_results
        except Exception as e:
            logging.error(f"Parsing UniProt BLAST XML with NCBIXML failed: {e}")
            logging.info("Attempting fallback parsing with ElementTree...")
            return self.parse_ebi_blast_xml(result_xml)
    
    def parse_ebi_blast_xml(self, result_xml):
        results = []
        try:
            root = ET.fromstring(result_xml)
            ns = {"ns": "http://www.ebi.ac.uk/schema"}
            hits_elem = root.find(".//ns:SequenceSimilaritySearchResult/ns:hits", ns)
            if hits_elem is None:
                logging.info("No hits found in EBI BLAST XML.")
                return results
            for hit in hits_elem.findall("ns:hit", ns):
                hit_id = hit.attrib.get("id", "")
                hit_title = hit.attrib.get("description", "")
                hit_length = hit.attrib.get("length", "")
                alignment = hit.find("ns:alignments/ns:alignment", ns)
                if alignment is None:
                    continue
                score_elem = alignment.find("ns:score", ns)
                expect_elem = alignment.find("ns:expectation", ns)
                query_seq_elem = alignment.find("ns:querySeq", ns)
                match_seq_elem = alignment.find("ns:matchSeq", ns)
                if query_seq_elem is not None:
                    try:
                        qs = int(query_seq_elem.attrib.get("start", "0"))
                        qe = int(query_seq_elem.attrib.get("end", "0"))
                        coverage = round(100 * (qe - qs + 1) / self.query_length, 2)
                    except Exception:
                        coverage = ""
                else:
                    coverage = ""
                results.append({
                    "Database": "UniProt (EBI BLAST)",
                    "Hit Title": hit_title,
                    "Hit ID": hit_id,
                    "Hit Length": hit_length,
                    "E-value": expect_elem.text if expect_elem is not None else "",
                    "Score": score_elem.text if score_elem is not None else "",
                    "Query Start": query_seq_elem.attrib.get("start", "") if query_seq_elem is not None else "",
                    "Query End": query_seq_elem.attrib.get("end", "") if query_seq_elem is not None else "",
                    "Hit Sequence": match_seq_elem.text if match_seq_elem is not None else "",
                    "Coverage (%)": coverage
                })
        except Exception as e:
            logging.error(f"Error in fallback XML parsing: {e}")
        return results
    
    def compile_results(self):
        self.load_sequence()
        results = []
        results.extend(self.run_ncbi_blast())
        results.extend(self.run_uniprot_blast())
        self.results = results
        return results
    
    def parse_organism_in_memory(self):
        for r in self.results:
            title = r.get("Hit Title", "")
            org = extract_organism_from_title(title)
            # Replace spaces with underscores in organism name
            if org:
                org = org.replace(" ", "_")
            r["Organism"] = org
    
    def save_to_excel(self):
        if not self.results:
            logging.info("No results to compile.")
            return
        df = pd.DataFrame(self.results)
        excel_filename = f"{self.prefix}_protein_search_results.xlsx"
        try:
            df.to_excel(excel_filename, index=False)
            logging.info(f"Results (including Organism) successfully saved to {excel_filename}")
        except Exception as e:
            logging.error(f"Error saving Excel file: {e}")


class PhylogeneticTreeBuilder:
    """
    Generate a phylogenetic tree from unique BLAST hit sequences.
    Leaf labels include the accession and organism name in the format:
         ACCESSION-ORGANISM
    """
    def __init__(self, results, prefix, query_id, query_sequence):
        self.results = results
        self.prefix = prefix
        self.seq_dict = {}
        self.query_id = query_id
        self.query_sequence = query_sequence
    
    def extract_unique_sequences(self):
        for r in self.results:
            seq = r.get("Hit Sequence", "").strip()
            if seq:
                if r.get("Database") == "UniProt (EBI BLAST)" and r.get("Hit ID"):
                    base = r.get("Hit ID")
                else:
                    base = r.get("Hit Title", "unknown").split()[0]
                org = r.get("Organism", "").strip()
                if org:
                    label = f"{base}-{org}"
                else:
                    label = base
                label = label.replace(":", "_")
                id_ = label
                counter = 1
                while id_ in self.seq_dict:
                    id_ = f"{label}_{counter}"
                    counter += 1
                self.seq_dict[id_] = seq
        # Add the query sequence with a distinct label using "QUERY_" prefix.
        simple_query = self.query_id.split("|")[-1].replace(":", "_")
        query_label = "QUERY_" + simple_query
        if query_label in self.seq_dict:
            counter = 1
            while query_label in self.seq_dict:
                query_label = "QUERY_" + simple_query + f"_{counter}"
                counter += 1
        self.seq_dict[query_label] = self.query_sequence
        return self.seq_dict
    
    def write_fasta(self):
        fasta_filename = f"{self.prefix}_hits.fasta"
        with open(fasta_filename, "w") as f:
            for id_, seq in self.seq_dict.items():
                f.write(f">{id_}\n{seq}\n")
        return fasta_filename
    
    def run_alignment(self, fasta_filename):
        aligned_filename = f"{self.prefix}_aligned_hits.fasta"
        try:
            subprocess.run(["muscle", "-align", fasta_filename, "-output", aligned_filename], check=True)
        except Exception as e:
            logging.error(f"Error running MUSCLE: {e}")
            return None
        return aligned_filename
    
    def build_tree(self, aligned_filename):
        alignment = AlignIO.read(aligned_filename, "fasta")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        return tree
    
    def generate_tree(self):
        newick_filename = f"{self.prefix}_phylogenetic_tree.nwk"
        pdf_filename = f"{self.prefix}_phylogenetic_tree.pdf"
        from ete3 import Tree, TreeStyle, TextFace
        
        # Force regeneration of Newick file by removing any existing file.
        if os.path.exists(newick_filename):
            os.remove(newick_filename)
        
        self.extract_unique_sequences()
        if not self.seq_dict:
            return
        fasta_filename = self.write_fasta()
        aligned_filename = self.run_alignment(fasta_filename)
        if not aligned_filename:
            return
        tree = self.build_tree(aligned_filename)
        from Bio import Phylo
        Phylo.write(tree, newick_filename, "newick")
        logging.info(f"Phylogenetic tree saved in Newick format to {newick_filename}")
        logging.info("Phylogenetic Tree (ASCII):")
        Phylo.draw_ascii(tree)
        try:
            # Load the freshly generated Newick file; use quoted_node_names=True
            ete_tree = Tree(newick_filename, format=1, quoted_node_names=True)
        except Exception as e:
            logging.error(f"Error loading Newick tree with format=1: {e}. Trying format=0.")
            ete_tree = Tree(newick_filename, format=0, quoted_node_names=True)
        
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.title.add_face(TextFace("Phylogenetic Tree", fsize=100), column=0)
        
        simple_query = self.query_id.split("|")[-1].replace(":", "_")
        for leaf in ete_tree.iter_leaves():
            if simple_query in leaf.name:
                face = TextFace(leaf.name, fgcolor="red")
            else:
                face = TextFace(leaf.name)
            leaf.add_face(face, column=0, position="branch-right")
        
        for node in ete_tree.traverse():
            nstyle = node.img_style
            nstyle["hz_line_width"] = 2
            nstyle["vt_line_width"] = 2
        
        n_seq = len(ete_tree.get_leaves())
        img_width = max(2000, int(n_seq * 200))
        try:
            ete_tree.render(pdf_filename, w=img_width, tree_style=ts)
            logging.info(f"Phylogenetic tree image saved to {pdf_filename}")
        except Exception as e:
            logging.error(f"Could not generate graphical tree: {e}")


def main():
    if len(sys.argv) < 2:
        logging.error("Usage: python NUMPy.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    email = "your.email@example.com"  # Replace with your actual email address.
    
    blast_searcher = ProteinBlastSearch(fasta_file, email)
    blast_searcher.compile_results()
    blast_searcher.parse_organism_in_memory()
    blast_searcher.save_to_excel()
    
    tree_builder = PhylogeneticTreeBuilder(
        blast_searcher.results,
        blast_searcher.prefix,
        blast_searcher.query_id,
        blast_searcher.query_sequence
    )
    tree_builder.generate_tree()


if __name__ == "__main__":
    main()
