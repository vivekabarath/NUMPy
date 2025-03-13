Below is a README that documents the project, including installation instructions, usage examples, and detailed explanations of each functionality.

---

# NUMPy: NCBI, UniProt & MUSCLE Phylogenetic Tree Builder

**NUMPy** is a command-line tool that performs protein BLAST searches (via NCBI and UniProt), extracts organism names from the BLAST hit titles, compiles results into an Excel file, and builds a phylogenetic tree from the unique hit sequences. The tree is rendered both as a Newick file and as a graphical PDF with the query sequence highlighted.

> **Note:** This project is unrelated to the popular Python numerical library **NumPy**.

## Features

- **BLAST Searches:**
  - Runs BLASTP against NCBI's nr database.
  - Runs BLASTP on UniProt via the EBI BLAST service.
  - Reuses existing XML files if available.

- **Result Parsing and Organism Extraction:**
  - Parses BLAST hit titles and extracts organism names using two methods:
    - **Square brackets:** e.g. `[Vitis vinifera]`
    - **OS=...OX= pattern:** e.g. `OS=Camellia fraterna OX=542725`
  - Spaces in organism names are replaced with underscores.
  
- **Excel Output:**
  - Generates an Excel file with detailed columns:  
    Database, Hit Title, Hit ID, Hit Length, E-value, Score, Query Start, Query End, Hit Sequence, Coverage (%), and Organism.

- **Phylogenetic Tree Construction:**
  - Builds a unique FASTA file for tree construction with labels in the format:
  
    ```
    ACCESSION-ORGANISM
    ```
    
    (with no spaces around the hyphen and colons replaced with underscores).
  - The query sequence is added with a distinct `"QUERY_"` prefix.
  - Uses MUSCLE 5.3 for multiple sequence alignment.
  - Constructs the tree using the neighbor-joining method.

- **Graphical Tree Output:**
  - Saves the tree in Newick format.
  - Renders a PDF tree using ETE3 with the query leaf highlighted in red and thicker branch lines for clarity.

## Installation

### Prerequisites

- **Python 3.6 or later**
- **MUSCLE 5.3**  
  Download and install MUSCLE 5.3 from [Robert Edgar's MUSCLE website](https://drive5.com/muscle/). Ensure the executable is in your system's PATH.
- **ETE3**  
  Install via pip:
  ```bash
  pip install ete3
  ```
- **Other Required Packages:**  
  Install via pip:
  ```bash
  pip install biopython pandas requests matplotlib openpyxl
  ```

### Clone the Repository

```bash
git clone https://github.com/yourusername/NUMPy.git
cd NUMPy
```

## Usage

Run the script from the command line by passing a FASTA file as an argument:

```bash
python NUMPy.py <fasta_file>
```

### Example

```bash
python NUMPy.py CYP97C27.fasta
```

### Expected Outputs

After running the script, several output files will be generated (all prefixed with the base name of your input FASTA file):

- **BLAST XML Files:**
  - `CYP97C27_ncbi_blast_result.xml`
  - `CYP97C27_uniprot_blast_result.xml`

- **Excel File:**
  - `CYP97C27_protein_search_results.xlsx`  
    Contains detailed BLAST results along with an "Organism" column where spaces are replaced by underscores.

- **FASTA Files:**
  - `CYP97C27_hits.fasta` (unique sequences with labels formatted as `ACCESSION-ORGANISM`)
  - `CYP97C27_aligned_hits.fasta` (alignment from MUSCLE)

- **Phylogenetic Tree Files:**
  - `CYP97C27_phylogenetic_tree.nwk` (Newick format with updated labels)
  - `CYP97C27_phylogenetic_tree.pdf` (graphical tree with the query leaf highlighted in red)

## How It Works

1. **Sequence Loading:**  
   The script reads the protein sequence from the provided FASTA file.

2. **BLAST Searches:**  
   BLASTP is executed on NCBI’s nr database and on UniProt via the EBI BLAST service. Existing XML files are reused if present.

3. **Result Parsing and Organism Extraction:**  
   For each BLAST hit, the script parses the "Hit Title" to extract the organism name using:
   - A square-bracket pattern (e.g. `[Vitis_vinifera]`)
   - An OS=...OX= pattern (e.g. `OS=Camellia_fraterna`)
   
   Spaces in the extracted organism names are replaced with underscores.

4. **Excel File Generation:**  
   The complete set of results (including the new "Organism" column) is saved into an Excel file.

5. **Phylogenetic Tree Construction:**  
   Unique hit sequences are extracted and labeled as `ACCESSION-ORGANISM` (with no spaces around the hyphen and colons sanitized). The query sequence is added with a "QUERY_" prefix. MUSCLE 5.3 is used to align these sequences, and a neighbor-joining tree is built from the alignment.

6. **Tree Visualization:**  
   The tree is saved in Newick format and rendered as a PDF using ETE3. The query sequence is highlighted in red for easy identification, and branch line widths are increased for better visualization.

## Troubleshooting

- **Duplicate Names:**  
  The script ensures unique labels by appending numeric suffixes if duplicate names occur. The query is given a distinct label to avoid conflicts.

- **Newick Format Issues:**  
  Colons in labels are replaced with underscores to prevent conflicts with the Newick format. If you experience issues loading the Newick file, ensure that the labels are sanitized as described.

- **Dependency Problems:**  
  Verify that MUSCLE 5.3 and ETE3 are installed and accessible in your PATH. Also, confirm that all required Python packages are installed.

## Contributing

Contributions are welcome! Please fork this repository, make your changes, and submit a pull request for review.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

This README provides complete documentation for the current script, including installation instructions, detailed functionality, usage examples, and troubleshooting tips.
