#!/usr/bin/env python3
"""
Extended NCBI GenBank Data Retriever with Filtering, CSV Export, and Visualization
"""

from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import os

class NCBIExtendedRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10Plus'
        self.records_data = []

    def search_taxid(self, taxid):
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        organism = Entrez.read(handle)[0]["ScientificName"]
        print(f"Organism: {organism}")
        term = f"txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y", retmax=1000)
        result = Entrez.read(handle)
        self.webenv = result["WebEnv"]
        self.query_key = result["QueryKey"]
        return int(result["Count"])

    def fetch_filtered_records(self, min_len, max_len, max_records=100):
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retmax=max_records,
            webenv=self.webenv,
            query_key=self.query_key
        )
        records = SeqIO.parse(handle, "gb")
        filtered = []
        for record in records:
            seq_len = len(record.seq)
            if min_len <= seq_len <= max_len:
                filtered.append({
                    "accession": record.id,
                    "length": seq_len,
                    "description": record.description
                })
        self.records_data = filtered
        return filtered

    def export_csv(self, filename="report.csv"):
        df = pd.DataFrame(self.records_data)
        df.to_csv(filename, index=False)
        print(f"Saved CSV report to {filename}")

    def generate_plot(self, filename="length_plot.png"):
        df = pd.DataFrame(self.records_data)
        df_sorted = df.sort_values(by="length", ascending=False)
        plt.figure(figsize=(10, 6))
        plt.plot(df_sorted["accession"], df_sorted["length"], marker="o")
        plt.xticks(rotation=90, fontsize=8)
        plt.ylabel("Sequence Length")
        plt.xlabel("Accession Number")
        plt.title("Sequence Lengths Sorted (Longest to Shortest)")
        plt.tight_layout()
        plt.savefig(filename)
        print(f"Saved plot to {filename}")

def main():
    email = input("Enter your email: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID: ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIExtendedRetriever(email, api_key)
    total = retriever.search_taxid(taxid)
    print(f"Total records found: {total}")
    records = retriever.fetch_filtered_records(min_len, max_len, max_records=200)
    print(f"Filtered records: {len(records)}")

    if records:
        retriever.export_csv("taxid_filtered_report.csv")
        retriever.generate_plot("taxid_sequence_plot.png")
    else:
        print("No records found within specified length range.")

if __name__ == "__main__":
    main()
