#!/usr/bin/env python3

"""
Download iNaturalist sequences for a given taxon ID number. 

Author: Alden Dirks
Date: November 21, 2025
Version: 1.0

Provided a taxon ID number, this script downloads all ITS sequences from
iNaturalist and outputs a formatted FASTA file. Provisional species names are 
used in place of consensus names when available. Sequences on iNaturalist are 
often a bit messier than those on GenBank, so basic cleaning steps are included. 

To reverse complement ITS sequences that are in the wrong orientation, use 
Alan Rockefeller's script `fixfasta.py` at: 
https://github.com/AlanRockefeller/fixfasta.py

Examples
- Get ITS sequences for the genus Pseudorhizina with taxon ID 951406
  fetch_inat_seqs.py 951406

- Specify a page limit to return fewer records. Useful if you just want to see total. 
  fetch_inat_seqs.py 47170 --max-pages 1 # 18,349,888 fungal observations

- Specify output file name and path. 
  fetch_inat_seqs.py 951406 --output seqs/inat.fasta
"""

import argparse
import csv
import re
import requests
import sys
import time


URL = "https://api.inaturalist.org/v1/observations"
PLACE_URL = "https://api.inaturalist.org/v1/places/"
place_cache = {}


def clean_sequence(seq):
    """Automatic cleaning of ITS sequences."""
    if not seq:
        return ""
    # Remove whitespace and newlines and make uppercase
    seq = seq.strip().replace(" ", "").replace("\n", "").upper()
    # Define valid IUPAC DNA characters
    valid = "ACGTRYWSMKHBVDN"
    # Find first stretch of ≥20 valid DNA characters
    pattern = f"[{valid}]{{20,}}"
    match = re.search(pattern, seq)
    if match:
        seq = seq[match.start():]
    else:
        seq = None
    return seq


def get_inat_observations(taxon_id=None, per_page=200, max_pages=None, delay=0.5):
    """Download iNaturalist observations"""
    page = 1
    observations = []
    total = 0

    while True:
        params = {"page": page, "per_page": per_page, "taxon_id": taxon_id,}
        try:
            response = requests.get(URL, params=params)
            response.raise_for_status()
            data = response.json()
            total = data.get("total_results", total)
            results = data.get("results", [])
        except requests.exceptions.RequestException as e:
            print(f"An error occurred on page {page}: {e}")
            break

        # Get only the observations with the field "DNA Barcode ITS"
        sequenced_observations = [
            r for r in results
            if any(
                ofv.get("name") == "DNA Barcode ITS"
                for ofv in r.get("ofvs", [])
            )
        ]

        # Add the sequenced observations to the observations list
        observations.extend(sequenced_observations)
        print(f"Fetched page {page}: {len(results)} observations, "
              f"{len(sequenced_observations)} with ITS")

        # Stop conditions
        if max_pages and page > max_pages:
            print("Max pages limit reached.")
            break
        if len(results) < per_page:
            print("Last page reached.")
            break

        page += 1
        time.sleep(delay)

    return total, observations


def get_place_info(place_id):
    """Look up a place by ID, with caching"""
    if place_id in place_cache:
        return place_cache[place_id]

    try:
        response = requests.get(f"{PLACE_URL}{place_id}")
        response.raise_for_status()
        data = response.json().get("results", [])
        if not data:
            place_cache[place_id] = None
            return None
        place_cache[place_id] = data[0]
        return data[0]
    except requests.RequestException:
        place_cache[place_id] = None
        return None
    

def extract_country_state(obs):
    """Given a single iNat observation JSON object, return (country, state) based on place_ids."""
    place_ids = obs.get("place_ids") or []
    country = None
    state = None

    for pid in place_ids:
        place = get_place_info(pid)
        if not place:
            continue

        level = place.get("admin_level")
        name  = place.get("name")

        if level == 0 and not country:
            country = name
            if country == "United States":
                country = "USA"

        elif level == 10 and not state:
            state = name

        # Stop early if both found
        if country and state:
            break

    return country or "NA", state or "NA"
    

def parse_fasta(observations):
    """Return FASTA headers and sequence data and metadata"""
    headers = []
    sequences = []
    metadata = []

    for obs in observations:
        species = obs.get("taxon", {}).get("name")
        species_rank = obs.get("taxon", {}).get("rank").lower().strip()
        if species_rank == "species":
            species = species.replace(" ", "-")
        elif species_rank == "genus":
            species += "-sp."
        else: 
            species = f"{species_rank}_{species}"
        inat_id = obs.get("id") or "NA"
        country, state = extract_country_state(obs)
        lon, lat = (obs.get("geojson", {}).get("coordinates") or [None, None])[:2]
        observed_on = obs.get("observed_on", "")
        user = obs.get("user", {}).get("login", "")

        # Find the "observation field" entry named "DNA Barcode ITS" and get sequence
        ofv = next(
            (field for field in obs.get("ofvs", []) if field.get("name") == "DNA Barcode ITS"),
            None
        )
        if not ofv:
            continue

        raw_seq = ofv.get("value", "")
        cleaned_seq = clean_sequence(raw_seq)
        if not cleaned_seq:
            print(f"[WARN] Sequence for {header} contains no valid DNA, skipping.")
            continue

        header = f"{species}_iNat{inat_id}_{country}-{state}".replace(" ", "_")
        headers.append(f">{header}")
        sequences.append(cleaned_seq)

        metadata.append({
            "header": header,
            "species": species,
            "species_rank": species_rank,
            "country": country,
            "state": state,
            "inat_id": inat_id,
            "latitude": lat,
            "longitude": lon,
            "observed_on": observed_on,
            "user": user,
            "raw_seq_length": len(raw_seq) if raw_seq else 0,
            "cleaned_seq_length": len(cleaned_seq),
        })
    
    return headers, sequences, metadata


def write_fasta(output_path, headers, sequences):
    """Write FASTA file"""
    count = 0
    with open(output_path, "w") as fasta:
        for header, seq in zip(headers, sequences):
            fasta.write(f"{header}\n")
            fasta.write(f"{seq}\n")
            count += 1
    print(f"\nWrote {count} records to FASTA: {output_path}")


def write_metadata_tsv(tsv_path, metadata):
    """Write metadata TSV file"""
    if not metadata:
        print("[WARN] No metadata to write.")
        return 0
    fieldnames = list(metadata[0].keys())
    with open(tsv_path, "w", newline="") as meta:
        writer = csv.DictWriter(meta, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in metadata:
            writer.writerow(row)
    print(f"Wrote metadata TSV: {tsv_path}")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Download ITS sequences from iNaturalist"
    )
    parser.add_argument("taxon_id", help="iNaturalist taxon ID")
    parser.add_argument("--output", metavar="OUTPUT_PATH", help="Path to write FASTA output", default="inat.fasta")
    parser.add_argument("--per_page", help="Number of observations to request per page (max 200)", type=int, default=200)
    parser.add_argument("--max-pages", help="Maximum number of pages to fetch (for testing)", type=int, default=None)
    parser.add_argument("--delay", help="Delay (seconds) between API page requests", type=float, default=0.5)
    
    # If no arguments are provided, show help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    # Enforce per_page max
    per_page = min(args.per_page, 200)
    if args.per_page > 200:
        print("[WARN] per_page > 200; using 200 (API limit).")

    # Download observations
    print(f"\nDownloading taxon_id={args.taxon_id} (per_page={per_page})...\n")
    total, observations = get_inat_observations(taxon_id=args.taxon_id, per_page=per_page, max_pages=args.max_pages, delay=args.delay)
    print(f"\nTotal number of observations: {total}")
    print(f"Observations with ITS sequence data: {len(observations)}")

    # Parse and clean FASTA
    print(f"\nParsing FASTA and fetching location information.")
    headers, sequences, metadata = parse_fasta(observations)

    # Write FASTA and metadata TSV
    fasta_path = args.output
    metadata_path = fasta_path.replace(".fasta", ".tsv")
    write_fasta(fasta_path, headers, sequences)
    write_metadata_tsv(metadata_path, metadata)

    print("\nDone.\n")


if __name__ == "__main__":
    main()