#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein BLAST Search and Phylogenetic Tree Builder

This module performs the following tasks:
  1. Reads a protein sequence from an input FASTA file.
  2. Executes BLAST searches against:
       - NCBI’s nr database using BLASTP.
       - UniProt via the EBI BLAST service.
  3. Saves (or reuses) the raw BLAST XML results using the FASTA file’s base name as a prefix.
  4. Parses BLAST hit data into an Excel file that includes the following columns:
       - Database, Hit Title, Hit ID (for UniProt hits), Hit Length, E-value, Score,
         Query Start, Query End, Hit Sequence, and Coverage (%).
  5. Builds a phylogenetic tree from unique hit sequences. For UniProt hits, the "Hit ID"
     is used when available; otherwise, the first token of "Hit Title" is used.
     The tree is saved in Newick format, printed as an ASCII tree, and a graphical PDF image is generated.
     
All generated files are prefixed with the input FASTA file’s base name.

Usage:
    python prot_Search.py <fasta_file>

Dependencies:
    - Biopython
    - pandas
    - requests
    - matplotlib (for generating the tree image)
    - MUSCLE 5.3 installed and accessible via the command "muscle" (with options -align and -output)
"""

import sys
import os
import time
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

# Configure logging for detailed output.
logging.basicConfig(level=logging.INFO, format='%(message)s')


class ProteinBlastSearch:
    """
    Class to perform BLAST searches and compile results.

    Attributes:
        fasta_file (str): Path to the input FASTA file.
        email (str): Email address used for BLAST API requests.
        query_sequence (str): Protein sequence extracted from the FASTA file.
        query_length (int): Length of the protein sequence.
        results (list): List of dictionaries containing parsed BLAST hit data.
        prefix (str): Base name of the input FASTA file (used as a prefix for all generated files).
    """

    def __init__(self, fasta_file, email):
        """
        Initialize the ProteinBlastSearch instance.

        Args:
            fasta_file (str): Path to the input FASTA file.
            email (str): Email address to use for BLAST requests.
        """
        self.fasta_file = fasta_file
        self.email = email
        self.query_sequence = ""
        self.query_length = 0
        self.results = []
        # Set the prefix to the FASTA file's base name (without extension).
        self.prefix = os.path.splitext(os.path.basename(self.fasta_file))[0]

    def load_sequence(self):
        """
        Load the protein sequence from the FASTA file.

        Sets:
            self.query_sequence (str): The protein sequence.
            self.query_length (int): The length of the protein sequence.

        Raises:
            SystemExit: If no sequences are found.
        """
        records = list(SeqIO.parse(self.fasta_file, "fasta"))
        if not records:
            logging.error("No sequences found in the provided FASTA file.")
            sys.exit(1)
        self.query_sequence = str(records[0].seq)
        self.query_length = len(self.query_sequence)
        logging.info("Loaded protein sequence: {}".format(records[0].id))

    def run_ncbi_blast(self):
        """
        Run BLAST search on NCBI's nr database using BLASTP.

        Checks if an XML file (<prefix>_ncbi_blast_result.xml) exists. If so, it reuses it;
        otherwise, it performs the search and saves the XML.

        Returns:
            list: List of dictionaries with BLAST hit information. Each dictionary includes:
                  "Database", "Hit Title", "Hit ID" (empty for NCBI), "Hit Length", "E-value", 
                  "Score", "Query Start", "Query End", "Hit Sequence", "Coverage (%)".
        """
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
                # Calculate coverage as (alignment length / query length) * 100.
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
        """
        Run BLAST search on UniProt via EBI's BLAST service.

        Checks if an XML file (<prefix>_uniprot_blast_result.xml) exists. If so, it reuses it;
        otherwise, it performs the search and saves the XML.

        Returns:
            list: List of dictionaries with BLAST hit information. For UniProt hits, includes "Hit ID".
        """
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
            logging.info("UniProt BLAST response status: {}".format(response.status_code))
            logging.info("Response content: {}".format(response.text))
            if response.status_code != 200:
                logging.error("Error: UniProt BLAST search request failed")
                return []
            job_id = response.text.strip()
            status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
            
            # Poll for job completion with a dynamic elapsed time display.
            start_time = time.time()
            while True:
                try:
                    status_response = requests.get(status_url)
                except Exception as e:
                    logging.error(f"Error polling UniProt BLAST status: {e}")
                    return []
                if status_response.text.strip() == "FINISHED":
                    print()  # End the dynamic line.
                    break
                elapsed = time.time() - start_time
                minutes = int(elapsed // 60)
                seconds = int(elapsed % 60)
                print(f"Waiting for UniProt BLAST job to finish... Elapsed time: {minutes:02d}:{seconds:02d}", end="\r", flush=True)
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
        
        # First attempt to parse using NCBIXML.
        try:
            blast_record = NCBIXML.read(StringIO(result_xml))
            uni_results = []
            for alignment in blast_record.alignments:
                # Try to extract the hit ID from the alignment object if available.
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
            logging.error("Parsing UniProt BLAST XML with NCBIXML failed: {}".format(e))
            logging.info("Attempting fallback parsing with ElementTree...")
            return self.parse_ebi_blast_xml(result_xml)
    
    def parse_ebi_blast_xml(self, result_xml):
        """
        Fallback parser for EBI BLAST XML using ElementTree.

        Args:
            result_xml (str): The XML content as a string.

        Returns:
            list: A list of dictionaries with parsed BLAST hit data. For UniProt hits,
                  includes the "Hit ID" from the XML.
        """
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
        """
        Load the protein sequence, run both BLAST searches, and compile the results.

        Returns:
            list: A combined list of BLAST result dictionaries.
        """
        self.load_sequence()
        results = []
        results.extend(self.run_ncbi_blast())
        results.extend(self.run_uniprot_blast())
        self.results = results
        return results

    def save_to_excel(self):
        """
        Save the compiled BLAST results to an Excel file.

        The Excel file is named <prefix>_protein_search_results.xlsx.
        """
        if not self.results:
            logging.info("No results to compile.")
            return
        df = pd.DataFrame(self.results)
        excel_filename = f"{self.prefix}_protein_search_results.xlsx"
        try:
            df.to_excel(excel_filename, index=False)
            logging.info(f"Results successfully saved to {excel_filename}")
        except Exception as e:
            logging.error(f"Error saving Excel file: {e}")


class PhylogeneticTreeBuilder:
    """
    Class to generate and visualize a phylogenetic tree from BLAST hit sequences.

    Attributes:
        results (list): List of BLAST result dictionaries.
        prefix (str): File name prefix from the input FASTA.
        seq_dict (dict): Mapping of unique sequence identifiers to hit sequences.
    """
    
    def __init__(self, results, prefix):
        """
        Initialize the PhylogeneticTreeBuilder instance.

        Args:
            results (list): BLAST results from which to extract sequences.
            prefix (str): The file name prefix from the input FASTA.
        """
        self.results = results
        self.prefix = prefix
        self.seq_dict = {}

    def extract_unique_sequences(self):
        """
        Extract unique hit sequences from the BLAST results.

        For UniProt hits, the "Hit ID" is used as the base identifier if available;
        otherwise, the first token of "Hit Title" is used.

        Returns:
            dict: A dictionary mapping unique identifiers to sequence strings.
        """
        for r in self.results:
            seq = r.get("Hit Sequence", "").strip()
            if seq:
                if r.get("Database") == "UniProt (EBI BLAST)" and r.get("Hit ID"):
                    id_base = r.get("Hit ID")
                else:
                    id_base = r.get("Hit Title", "unknown").split()[0]
                id_ = id_base
                counter = 1
                while id_ in self.seq_dict:
                    id_ = f"{id_base}_{counter}"
                    counter += 1
                self.seq_dict[id_] = seq
        if not self.seq_dict:
            logging.info("No hit sequences available for phylogenetic analysis.")
        return self.seq_dict

    def write_fasta(self):
        """
        Write the unique hit sequences to a FASTA file.

        Returns:
            str: The FASTA filename (<prefix>_hits.fasta).
        """
        fasta_filename = f"{self.prefix}_hits.fasta"
        with open(fasta_filename, "w") as f:
            for id_, seq in self.seq_dict.items():
                f.write(f">{id_}\n{seq}\n")
        return fasta_filename

    def run_alignment(self, fasta_filename):
        """
        Run MUSCLE 5.3 to perform a multiple sequence alignment.

        Uses MUSCLE 5.3 options: -align <input> -output <output>.

        Args:
            fasta_filename (str): Input FASTA file with unique sequences.

        Returns:
            str or None: The aligned FASTA filename (<prefix>_aligned_hits.fasta) or None if MUSCLE fails.
        """
        aligned_filename = f"{self.prefix}_aligned_hits.fasta"
        try:
            subprocess.run(["muscle", "-align", fasta_filename, "-output", aligned_filename], check=True)
        except Exception as e:
            logging.error("Error running MUSCLE: {}".format(e))
            return None
        return aligned_filename

    def build_tree(self, aligned_filename):
        """
        Construct a phylogenetic tree from the multiple sequence alignment.

        Args:
            aligned_filename (str): Path to the aligned FASTA file.

        Returns:
            Bio.Phylo.BaseTree.Tree: The generated phylogenetic tree.
        """
        alignment = AlignIO.read(aligned_filename, "fasta")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        return tree

    def generate_tree(self):
        """
        Generate the phylogenetic tree and save outputs.

        The tree is saved in Newick format (<prefix>_phylogenetic_tree.nwk),
        printed as an ASCII tree, and a graphical PDF image is saved (<prefix>_phylogenetic_tree.pdf).
        Adjustments with plt.tight_layout() are applied for better readability.
        """
        self.extract_unique_sequences()
        if not self.seq_dict:
            return
        fasta_filename = self.write_fasta()
        aligned_filename = self.run_alignment(fasta_filename)
        if not aligned_filename:
            return
        tree = self.build_tree(aligned_filename)
        newick_filename = f"{self.prefix}_phylogenetic_tree.nwk"
        pdf_filename = f"{self.prefix}_phylogenetic_tree.pdf"
        # Save the tree in Newick format.
        Phylo.write(tree, newick_filename, "newick")
        logging.info(f"Phylogenetic tree saved in Newick format to {newick_filename}")
        # Print an ASCII representation of the tree.
        logging.info("Phylogenetic Tree (ASCII):")
        Phylo.draw_ascii(tree)
        # Generate a graphical tree image with improved layout.
        try:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(12, 10))
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, do_show=False, axes=axes)
            plt.tight_layout()  # Adjust layout to avoid overlapping text.
            # Optionally, adjust margins manually:
            # plt.subplots_adjust(left=0.3, right=0.95, top=0.95, bottom=0.1)
            plt.savefig(pdf_filename)
            logging.info(f"Phylogenetic tree image saved to {pdf_filename}")
        except Exception as e:
            logging.error("Could not generate graphical tree: {}".format(e))


def main():
    """
    Main pipeline function:

    1. Validates command-line input.
    2. Creates a ProteinBlastSearch instance to run BLAST searches.
    3. Saves BLAST results to an Excel file.
    4. Creates a PhylogeneticTreeBuilder instance using the BLAST results and generates the tree.
    """
    if len(sys.argv) < 2:
        logging.error("Usage: python prot_Search.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    email = "your.email@example.com"  # Replace with your actual email address.
    
    # Create the BLAST search object and run the searches.
    blast_searcher = ProteinBlastSearch(fasta_file, email)
    blast_searcher.compile_results()
    blast_searcher.save_to_excel()
    
    # Create the tree builder and generate the phylogenetic tree.
    tree_builder = PhylogeneticTreeBuilder(blast_searcher.results, blast_searcher.prefix)
    tree_builder.generate_tree()


if __name__ == "__main__":
    main()
