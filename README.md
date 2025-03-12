Below is a sample complete GitHub README documentation for the project:

---

# NUMPy: NCBI, UniProt & MUSCLE Phylogenetic Tree Builder

NUMPy is a command-line tool for protein sequence analysis that integrates three key components:
- **NCBI BLAST**: Search NCBI's non-redundant (nr) protein database.
- **UniProt BLAST**: Retrieve BLAST results from UniProt via the EBI BLAST service.
- **MUSCLE Alignment & Phylogenetics**: Align sequences using MUSCLE 5.3 and build a phylogenetic tree.

All output files are consistently prefixed with the input FASTA file's base name, ensuring organized and traceable results.

## Features

- **BLAST Searches**  
  - Performs BLASTP searches on NCBI's nr database.
  - Executes BLAST searches on UniProt via the EBI BLAST service.
  - Reuses previously generated XML files (if available) to save time.

- **Comprehensive Reporting**  
  - Parses BLAST hit data into an Excel file with the following columns:
    - **Database**: Source of the hit (NCBI or UniProt).
    - **Hit Title**: Descriptive title of the hit.
    - **Hit ID**: Unique identifier (for UniProt hits).
    - **Hit Length**: Length of the hit sequence.
    - **E-value**: Statistical significance of the alignment.
    - **Score**: Alignment score.
    - **Query Start / End**: Alignment positions on the query.
    - **Hit Sequence**: Aligned segment from the database hit.
    - **Coverage (%)**: Percentage of the query sequence covered by the alignment.

- **Phylogenetic Analysis**  
  - Extracts unique hit sequences.
  - Performs multiple sequence alignment using MUSCLE 5.3.
  - Constructs a neighbor-joining phylogenetic tree.
  - Outputs the tree in:
    - **Newick format** for further analysis.
    - **ASCII** representation in the terminal.
    - **PDF image** with dynamic sizing based on the number of sequences.

- **Optimized Workflow**  
  - If all files up to the Newick tree are already generated, the tool will load the existing tree and immediately generate the PDF and print the ASCII tree.

## Installation

### Prerequisites

- **Python 3.x**  
  NUMPy is written in Python. Ensure you have Python 3 installed.

- **Python Dependencies**  
  Install the required Python libraries:
  ```bash
  pip install biopython pandas requests matplotlib
  ```
  Alternatively, use a virtual environment:
  ```bash
  python3 -m venv venv
  source venv/bin/activate
  pip install biopython pandas requests matplotlib
  ```

- **MUSCLE 5.3**  
  NUMPy uses MUSCLE 5.3 for multiple sequence alignment.  
  **On Ubuntu**, you can install MUSCLE via:
  ```bash
  sudo apt update
  sudo apt install muscle
  ```
  If the repository version is outdated, download and compile MUSCLE 5.3 from the [official MUSCLE website](https://www.drive5.com/muscle/). Make sure the executable is accessible as `muscle` in your PATH.

## Usage

Run the script with the input FASTA file:
```bash
python NUMPy.py <input_fasta_file>
```
For example, if your file is `CYP97C27.fasta`, execute:
```bash
python NUMPy.py CYP97C27.fasta
```

## Output Files

All output files will use the input FASTA file’s base name as a prefix. For example, for `CYP97C27.fasta`, you will get:

- `CYP97C27_ncbi_blast_result.xml` – NCBI BLAST XML result.
- `CYP97C27_uniprot_blast_result.xml` – UniProt BLAST XML result.
- `CYP97C27_protein_search_results.xlsx` – Excel file containing detailed BLAST hit data.
- `CYP97C27_hits.fasta` – FASTA file with unique hit sequences.
- `CYP97C27_aligned_hits.fasta` – Multiple sequence alignment of the unique hits.
- `CYP97C27_phylogenetic_tree.nwk` – Newick formatted phylogenetic tree.
- `CYP97C27_phylogenetic_tree.pdf` – Graphical PDF image of the phylogenetic tree (with dynamic sizing).

If the Newick tree file already exists, NUMPy will load it directly to generate the PDF and print the ASCII tree without re-running the alignment and tree-building pipeline.

## How It Works

1. **Sequence Loading**  
   The script reads the first protein sequence from the provided FASTA file.

2. **BLAST Searches**  
   - It performs a BLAST search against NCBI’s nr database.
   - It performs a BLAST search against UniProt via the EBI BLAST service.
   - For both services, if the XML results exist from a previous run (using the file prefix), those are reused.

3. **Result Compilation**  
   The parsed BLAST results are compiled into an Excel file with detailed information, including the unique "Hit ID" for UniProt results.

4. **Phylogenetic Tree Construction**  
   - Unique hit sequences are extracted.
   - MUSCLE 5.3 is invoked to align these sequences.
   - A neighbor-joining tree is constructed from the alignment.
   - The tree is saved in Newick format, printed as an ASCII tree, and a PDF image is generated. The PDF image size is dynamically adjusted based on the number of sequences to ensure readability.

## Example

Assuming you have a file `CYP97C27.fasta`, run:
```bash
python NUMPy.py CYP97C27.fasta
```
This will produce all the output files mentioned above, with the phylogenetic tree displayed in your terminal and saved as a PDF.

## Contributing

Contributions, issues, and feature requests are welcome!  
Feel free to check the [issues page](https://github.com/vivekabarath/NUMPy/issues).

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

Feel free to adjust the URLs (e.g., for issues and license) to match your repository. This documentation provides a comprehensive overview of the project's functionality, installation, usage, and output details.
