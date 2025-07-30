# PubMed Paper Fetcher

This command-line tool fetches research papers from PubMed based on a search query. It extracts key metadata such as PubMed ID, Title, Publication Date, Non-academic Authors, Company Affiliations, and Corresponding Author Email, and exports the results to a CSV file.

---

## Project Structure

pubmed-paper-fetcher/
├── pubmed_paper_fetcher/
│ ├── init.py
│ ├── cli.py # Command-line interface logic
│ └── fetcher.py # Helper functions for fetching and parsing data
├── pyproject.toml # Poetry configuration file
├── README.md # Project documentation
└── output.csv # Example output file (if saved using --file)


## Dependencies & Installation

This project uses **Python 3.10+** and **Poetry** for dependency management.

### Steps to set up:

1. **Clone the repository:**

```bash
git clone https://github.com//pubmed-paper-fetcher.git
cd pubmed-paper-fetcher


Install Poetry (if not already installed):
pip install poetry

Install dependencies using Poetry:
poetry install

To fetch papers using a search query:
poetry run get-papers-list "your search query"
Available Options:
-h, --help
Show usage instructions.

-d, --debug
Print debug information during execution.

-f, --file <filename>
Save the output to a CSV file. If not specified, the output is printed to the console.

Example:
poetry run get-papers-list "covid vaccine" --file covid.csv --debug
Output Fields
The CSV file contains the following columns:

PubmedID – Unique identifier for the paper

Title – Title of the paper

Publication Date – Date the paper was published

Non-academic Author(s) – Names of authors affiliated with non-academic institutions

Company Affiliation(s) – Names of pharmaceutical/biotech companies

Corresponding Author Email – Email of the corresponding author

Note: Some papers may not contain all fields (e.g., email or company affiliation), depending on the metadata available in PubMed.

Tools and Libraries Used
Biopython – For accessing PubMed via Entrez API

argparse – For command-line argument parsing

csv – For CSV output generation

Poetry – For packaging and dependency management