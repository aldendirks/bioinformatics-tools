#!/usr/bin/env python3

"""
Download GenBank ITS sequences for a given taxon name.

Author: Alden Dirks
Date: December 8, 2025
Version: 0.2

Given a taxon name, this script downloads all available ITS sequences from
GenBank and saves them in FASTA format. Sequence headers are returned exactly
as provided by GenBank. Optionally, only type material can be fetched.

Features:
- Uses NCBI Entrez via Biopython to search and fetch sequences.
- Reformats FASTA headers and adds geographic information if available.
- Supports type-material-only filtering.
- Prompts the user before downloading large sequence sets.
- Allows user-specified output file path.

Examples:
    fetch_gb_seqs.py Pseudorhizina
    fetch_gb_seqs.py Pseudorhizina --type_only
    fetch_gb_seqs.py Pseudorhizina --output seqs/genbank.fasta
"""

import argparse
import os
import re
import sys
from pathlib import Path
from Bio import Entrez, SeqIO
from urllib.error import HTTPError


# Set email; API key optional (recommended for higher rate limits)
Entrez.email = os.getenv("NCBI_EMAIL", "aldendirks@gmail.com")
api_key = os.getenv("NCBI_API_KEY")
if api_key:
    Entrez.api_key = api_key
BATCH_SIZE = 20


def us_state_abbr(state_name):
    states = {
        "Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR",
        "California": "CA", "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE",
        "Florida": "FL", "Georgia": "GA", "Hawaii": "HI", "Idaho": "ID", "Illinois": "IL",
        "Indiana": "IN", "Iowa": "IA", "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA",
        "Maine": "ME", "Maryland": "MD", "Massachusetts": "MA", "Michigan": "MI",
        "Minnesota": "MN", "Mississippi": "MS", "Missouri": "MO", "Montana": "MT",
        "Nebraska": "NE", "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ",
        "New Mexico": "NM", "New York": "NY", "North Carolina": "NC", "North Dakota": "ND",
        "Ohio": "OH", "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA",
        "Rhode Island": "RI", "South Carolina": "SC", "South Dakota": "SD",
        "Tennessee": "TN", "Texas": "TX", "Utah": "UT", "Vermont": "VT",
        "Virginia": "VA", "Washington": "WA", "West Virginia": "WV",
        "Wisconsin": "WI", "Wyoming": "WY"
    }
    return states.get(state_name.strip())


def can_province_abbr(name):
    provinces = {
        "Alberta": "AB",
        "British Columbia": "BC",
        "Manitoba": "MB",
        "New Brunswick": "NB",
        "Newfoundland and Labrador": "NL",
        "Nova Scotia": "NS",
        "Ontario": "ON",
        "Prince Edward Island": "PE",
        "Quebec": "QC",
        "Saskatchewan": "SK",
        "Northwest Territories": "NT",
        "Nunavut": "NU",
        "Yukon": "YT"
    }
    return provinces.get(name.strip())


def build_query(taxon_name, type_only):
    """Construct an Entrez query string."""
    query = f"{taxon_name}[Organism] AND internal[All Fields]"
    if type_only:
        query += " AND type_material[Properties]"
    return query


def get_sequence_ids(query):
    """Return a list of sequence IDs matching the Entrez query."""
    try:
        with Entrez.esearch(db="nucleotide", term=query, retmax=0) as handle:
            record = Entrez.read(handle)
    except HTTPError as e:
        sys.exit(f"Entrez search failed: {e}")

    count = int(record["Count"])
    print(f"Found {count} sequences.")

    if count > 1000:
        confirm = input(f"Downloading {count} sequences may take a while. Continue? (y/n): ")
        if confirm.strip().lower() != "y":
            sys.exit("Aborting.")

    if count == 0:
        sys.exit("No sequences found for this query.")

    try:
        with Entrez.esearch(db="nucleotide", term=query, retmax=count) as handle:
            full_record = Entrez.read(handle)
    except HTTPError as e:
        sys.exit(f"Entrez ID retrieval failed: {e}")

    return full_record["IdList"]


def fetch_fasta_by_ids(ids):
    """Fetch FASTA records given a list of Entrez IDs."""
    print(f"Fetching {len(ids)} sequences...")

    try:
        with Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text") as handle:
            return handle.read()
    except HTTPError as e:
        sys.exit(f"Entrez efetch failed: {e}")


def parse_accession(record_id):
    """Extract accession without version (AB123456.1 → AB123456)."""
    acc = record_id.split(".")[0]
    return acc


def fetch_geo_info_batch(accessions):
    """Fetch /geo_loc_name or /country metadata."""
    geo_dict = {}

    try:
        with Entrez.efetch(
            db="nucleotide",
            id=",".join(accessions),
            rettype="gb",
            retmode="text"
        ) as handle:
            text = handle.read()

        for acc in accessions:
            # isolate each record
            pattern = rf"LOCUS.*?\n.*?ACCESSION\s+{re.escape(acc)}.*?\n(.*?)(?=LOCUS|\Z)"
            match = re.search(pattern, text, re.DOTALL)
            if not match:
                geo_dict[acc] = "NA"
                continue

            rec = match.group(1)

            # geo_loc_name has priority
            geo = None
            m = re.search(r'/geo_loc_name="([^"]+)"', rec)
            if m:
                geo = m.group(1)
            else:
                m = re.search(r'/country="([^"]+)"', rec)
                geo = m.group(1) if m else "NA"

            if geo.startswith("USA:") or geo.startswith("United States:"):
                state = geo.split(":")[1].strip()
                abbr = us_state_abbr(state)
                geo = f"USA-{state}" if abbr else "USA"
            elif geo.startswith("Canada:"):
                prov = geo.split(":")[1].strip()
                abbr = can_province_abbr(prov)
                geo = f"Canada-{prov}" if abbr else "Canada"
            else:
                geo = geo.split(":", 1)[0].strip().replace(" ", "_")

            geo_dict[acc] = geo

    except Exception:
        for acc in accessions:
            geo_dict[acc] = "NA"

    return geo_dict


def fetch_all_geo_metadata(records):
    """Batch GEO metadata for all FASTA records."""
    accessions = [parse_accession(rec.id) for rec in records]
    geo_info = {}

    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i:i+BATCH_SIZE]
        geo_info.update(fetch_geo_info_batch(batch))

    return geo_info


def reformat_header(original_header, geo, acc):
    """Reformat FASTA header to include species, accession, and geo info."""
    parts = original_header.split()

    if len(parts) < 3:
        # fallback: keep ID only
        species = parts[1] if len(parts) > 1 else "Unknown_sp"
    else:
        species = f"{parts[1]}_{parts[2]}"

    return f"{species}_{acc}_{geo}"


def write_reformatted_fasta(records, geo_info, output_path):
    """Write sequences with reformatted headers."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as out:
        for rec in records:
            acc = parse_accession(rec.id)
            geo = geo_info.get(acc, "NA")

            new_id = reformat_header(rec.description, geo, acc)
            rec.id = new_id
            rec.description = ""

            SeqIO.write(rec, out, "fasta")

    print(f"Sequences written to {output_path.resolve()}")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Download ITS sequences from GenBank (with GEO metadata and cleaned headers)."
    )

    parser.add_argument("taxon_name")
    parser.add_argument("--output", help="Destination FASTA file.", default=None)
    parser.add_argument("--type_only", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    print(f"\nSearching GenBank for '{args.taxon_name}'...\n")

    query = build_query(args.taxon_name, args.type_only)
    ids = get_sequence_ids(query)

    fasta_data = fetch_fasta_by_ids(ids)

    # temporarily write raw fasta
    tmp_path = Path("tmp_download.fasta")
    tmp_path.write_text(fasta_data)

    records = list(SeqIO.parse(tmp_path, "fasta"))
    geo_info = fetch_all_geo_metadata(records)

    # determine output
    if args.output:
        output_path = Path(args.output)
    else:
        safe = args.taxon_name.replace(" ", "_")
        output_path = Path(f"{safe}_ITS_genbank.fasta")

    write_reformatted_fasta(records, geo_info, output_path)

    tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()