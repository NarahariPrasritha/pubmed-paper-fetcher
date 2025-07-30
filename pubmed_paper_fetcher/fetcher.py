# fetcher.py

# Required imports
from typing import List, Dict
from Bio import Entrez # Biopython's module for accessing NCBI databases like PubMed

# Required by NCBI to identify the user accessing their services
Entrez.email = "prasrijyo2428@gmail.com"

# Keywords to identify company affiliations
COMPANY_KEYWORDS = ["pharma", "biotech", "inc", "ltd", "laboratories", "corp", "company", "healthcare", "diagnostics"]

# Keywords to identify academic affiliations
ACADEMIC_KEYWORDS = ["university", "college", "institute", "school", "faculty", "hospital", "department"]

# Function to determine if the affiliation is likely from a company
def is_company_affiliation(affil: str) -> bool:
    affil_lower = affil.lower()
    return any(keyword in affil_lower for keyword in COMPANY_KEYWORDS)
    
# Function to determine if the affiliation is from an academic institute
def is_academic_affiliation(affil: str) -> bool:
    affil_lower = affil.lower()  # Convert to lowercase for keyword matching
    return any(keyword in affil_lower for keyword in ACADEMIC_KEYWORDS)

# Function to search PubMed using a query string and return a list of PubMed IDs
def search_pubmed(query: str, max_results: int = 10) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)  # Read and parse the search results
    handle.close()
    return record["IdList"] # Return a list of PubMed IDs
    
# Function to fetch and filter articles with at least one company-affiliated author
def fetch_pubmed_papers(query: str, max_results: int = 10, debug: bool = False) -> List[Dict[str, str]]:
    ids = search_pubmed(query, max_results)
    if not ids:
        return []  # Return empty list if no results
        
    # Fetch full article details in XML format
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    papers = [] # List to hold filtered paper details
    for article in records["PubmedArticle"]:
        try:
            # Extract the PubMed ID
            pmid = article["MedlineCitation"]["PMID"]
            article_data = article["MedlineCitation"]["Article"]
             # Extract the article title
            title = article_data.get("ArticleTitle", "")
             # Extract and format the publication date
            pub_date = article_data["Journal"]["JournalIssue"]["PubDate"]
            pub_date_str = f"{pub_date.get('Year', '')}-{pub_date.get('Month', '')}-{pub_date.get('Day', '')}".strip("-")
            # Initialize storage for author and affiliation info
            non_academic_authors = []
            company_affiliations = set()
            corresponding_email = ""
             # Loop through all authors and check their affiliations
            for author in article_data.get("AuthorList", []):
                # Combine fore and last names
                full_name = f"{author.get('ForeName', '')} {author.get('LastName', '')}".strip()
                affils = author.get("AffiliationInfo", [])
                
                for affil in affils:
                    affil_text = affil.get("Affiliation", "")
                    if not affil_text:
                        continue

                    # Check if affiliation is from a company but not academic
                    if is_company_affiliation(affil_text) and not is_academic_affiliation(affil_text):
                        non_academic_authors.append(full_name)
                        company_affiliations.add(affil_text)

                    # Try to extract an email address from the affiliation
                    if "@" in affil_text and not corresponding_email:
                        words = affil_text.replace(">", " ").replace("<", " ").split()
                        for word in words:
                            if "@" in word:
                                corresponding_email = word.strip(",.;()<>")
                                break
                                
            # Add the collected data to the results
            papers.append({
                "PubmedID": str(pmid),
                "Title": title,
                "Publication Date": pub_date_str,
                "Non-academic Author(s)": ", ".join(non_academic_authors),
                "Company Affiliation(s)": "; ".join(company_affiliations),
                "Corresponding Author Email": corresponding_email
            })

            # Print each paper if debug mode is on
            if debug:
                print(papers[-1])
                print("-" * 60)

        except Exception as e:
            if debug:
                print(f"Skipping article due to error: {e}")
            continue

    return papers






