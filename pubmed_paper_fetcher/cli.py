# cli.py

import argparse
import csv
from .fetcher import fetch_pubmed_papers

def save_to_csv(papers, filename):
    if not papers:
        print("No papers found.")
        return

    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "PubmedID",
            "Title",
            "Publication Date",
            "Non-academic Author(s)",
            "Company Affiliation(s)",
            "Corresponding Author Email"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for paper in papers:
            writer.writerow(paper)

    print(f"Results saved to {filename}")

def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers based on a search query.")
    parser.add_argument("query", help="Search query for PubMed")
    parser.add_argument("--max", type=int, default=10, help="Maximum number of results to fetch")
    parser.add_argument("--file", default="results.csv", help="Output CSV file name")
    parser.add_argument("--debug", action="store_true", help="Enable debug printing")

    args = parser.parse_args()

    papers = fetch_pubmed_papers(args.query, max_results=args.max, debug=args.debug)
    save_to_csv(papers, args.file)

if __name__ == "__main__":
    main()
