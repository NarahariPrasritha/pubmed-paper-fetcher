# fetcher.py

from typing import List, Dict
from Bio import Entrez

Entrez.email = "prasrijyo2428@gmail.com"

COMPANY_KEYWORDS = ["pharma", "biotech", "inc", "ltd", "laboratories", "corp", "company", "healthcare", "diagnostics"]
ACADEMIC_KEYWORDS = ["university", "college", "institute", "school", "faculty", "hospital", "department"]

def is_company_affiliation(affil: str) -> bool:
    affil_lower = affil.lower()
    return any(keyword in affil_lower for keyword in COMPANY_KEYWORDS)

def is_academic_affiliation(affil: str) -> bool:
    affil_lower = affil.lower()
    return any(keyword in affil_lower for keyword in ACADEMIC_KEYWORDS)

def search_pubmed(query: str, max_results: int = 10) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_pubmed_papers(query: str, max_results: int = 10, debug: bool = False) -> List[Dict[str, str]]:
    ids = search_pubmed(query, max_results)
    if not ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    papers = []
    for article in records["PubmedArticle"]:
        try:
            pmid = article["MedlineCitation"]["PMID"]
            article_data = article["MedlineCitation"]["Article"]
            title = article_data.get("ArticleTitle", "")
            pub_date = article_data["Journal"]["JournalIssue"]["PubDate"]
            pub_date_str = f"{pub_date.get('Year', '')}-{pub_date.get('Month', '')}-{pub_date.get('Day', '')}".strip("-")

            non_academic_authors = []
            company_affiliations = set()
            corresponding_email = ""

            for author in article_data.get("AuthorList", []):
                full_name = f"{author.get('ForeName', '')} {author.get('LastName', '')}".strip()
                affils = author.get("AffiliationInfo", [])
                for affil in affils:
                    affil_text = affil.get("Affiliation", "")
                    if not affil_text:
                        continue
                    if is_company_affiliation(affil_text) and not is_academic_affiliation(affil_text):
                        non_academic_authors.append(full_name)
                        company_affiliations.add(affil_text)
                    if "@" in affil_text and not corresponding_email:
                        words = affil_text.replace(">", " ").replace("<", " ").split()
                        for word in words:
                            if "@" in word:
                                corresponding_email = word.strip(",.;()<>")
                                break

            papers.append({
                "PubmedID": str(pmid),
                "Title": title,
                "Publication Date": pub_date_str,
                "Non-academic Author(s)": ", ".join(non_academic_authors),
                "Company Affiliation(s)": "; ".join(company_affiliations),
                "Corresponding Author Email": corresponding_email
            })

            if debug:
                print(papers[-1])
                print("-" * 60)

        except Exception as e:
            if debug:
                print(f"Skipping article due to error: {e}")
            continue

    return papers






