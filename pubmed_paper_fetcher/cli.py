# cli.py
# Importing required modules
import argparse    # For building command-line interfaces
import csv        # To save results into CSV files
from .fetcher import fetch_pubmed_papers    # Importing the function to fetch papers

# Function to save the fetched paper data to a CSV file
def save_to_csv(papers, filename):
    if not papers:
        print("No papers found.")
        return

    # Open the specified file for writing CSV content
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        # Define the column headers for the CSV
        fieldnames = [
            "PubmedID",
            "Title",
            "Publication Date",
            "Non-academic Author(s)",
            "Company Affiliation(s)",
            "Corresponding Author Email"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()    # Write the header row

        # Write each paper’s data as a row in the CSV
        for paper in papers:
            writer.writerow(paper)

    print(f"Results saved to {filename}")  # Notify the user

# Main function to handle command-line interaction

def main():
    # Create an argument parser for command-line options
    parser = argparse.ArgumentParser(description="Fetch PubMed papers based on a search query.")
    parser.add_argument("query", help="Search query for PubMed")
    parser.add_argument("--max", type=int, default=10, help="Maximum number of results to fetch")
    parser.add_argument("--file", default="results.csv", help="Output CSV file name")
    parser.add_argument("--debug", action="store_true", help="Enable debug printing")

    args = parser.parse_args()
    
    # Save the fetched results into a CSV file
    papers = fetch_pubmed_papers(args.query, max_results=args.max, debug=args.debug)
    save_to_csv(papers, args.file)
    
# Entry point of the script when run directly
if __name__ == "__main__":
    main()
