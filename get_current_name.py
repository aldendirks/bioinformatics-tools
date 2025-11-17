#!/usr/bin/env python3

"""
Get the current MycoBank name for a list of species names (batch querying). 

Author: Alden Dirks
Date: November 14, 2025
Version: 1.0

Script to query MycoBank API for current names of fungal species. Results are 
saved to a tab-delimited text file. Excluded species are saved to a separate file.

Examples: 
- Query MycoBank for a list of species in species.txt.
  get_current_name.py species.txt

- Exclude species listed in excluded_species.txt from being queried. 
  get_current_name.py species.txt --exclude excluded_species.txt

- Specify batch size and output files.
  get_current_name.py species.txt --output results.tsv --excluded-output excluded.tsv --verbose

"""


import argparse
import csv
import os
import requests
import sys
import time
from itertools import islice


ACCESS_TOKEN = os.getenv("MYCOBANK_ACCESS_TOKEN")
if not ACCESS_TOKEN:
    raise EnvironmentError("MYCOBANK_ACCESS_TOKEN environment variable not set.")
ACCESS = {"Authorization": f"Bearer {ACCESS_TOKEN}"}
BASE_URL = "https://webservices.bio-aware.com/cbsdatabase_new/mycobank/taxonnames"
MB_URL_TEMPLATE = "https://www.mycobank.org/page/Name%20details%20page/field/Mycobank%20%23/"


def read_species_list(file_path, excluded_species=None, 
        excluded_output_file="excluded_species.tsv", verbose=False):
    """Read species from a file, excluding provisional/ambiguous or user-specified species."""
    
    if excluded_species is None:
        excluded_species = []
    species_to_check = []
    excluded_characters = ['[', ']', 'sp.', 'cf.', 'aff.']
    excluded_rows = []

    # Summary counters
    excluded_user_count = 0
    excluded_ambig_count = 0

    with open(file_path) as f:
        for line in f:
            sp = line.strip()
            if not sp:
                continue
            if sp in excluded_species:
                reason = "user_excluded"
                excluded_user_count += 1
            elif any(char in sp for char in excluded_characters):
                reason = "ambiguous"
                excluded_ambig_count += 1
            else:
                reason = None

            if reason:
                if verbose:
                    print(f"Skipping '{sp}' ({reason})")
                excluded_rows.append([sp, reason])
            else:
                if verbose:
                    print(f"Will check '{sp}'")
                species_to_check.append(sp)

    # Write excluded species to file if specified
    if excluded_rows:
        with open(excluded_output_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["species", "reason"])
            writer.writerows(excluded_rows)
        if verbose:
            print(f"Excluded species written to {excluded_output_file}\n")

    return species_to_check, excluded_user_count, excluded_ambig_count


def batch(iterable, n=20):
    """Yield successive n-sized batches from iterable."""
    it = iter(iterable)
    while True:
        chunk = list(islice(it, n))
        if not chunk:
            return
        yield chunk


def print_progress(processed_count, total_species, bar_length=20):
    """Prints an in-place progress bar."""
    progress_fraction = processed_count / total_species
    filled_length = int(bar_length * progress_fraction)
    bar = '#' * filled_length + '-' * (bar_length - filled_length)
    print(f"Processed {processed_count}/{total_species} species: [{bar}]", end="\r", flush=True)


def get_current_names_batch(species_list, output_file="mycobank_results.tsv", batch_size=20, verbose=False):
    """Query MycoBank API in batches and write results to TSV."""
    
    results = []
    total_species = len(species_list)
    processed_count = 0

    # Summary counters
    summary = {
        "current": 0,
        "not_current": 0,
        "no_records": 0,
        "no_valid": 0,
        "no_current": 0,
        "multiple_records": 0,
        "error": 0,
    }

    for species_batch in batch(species_list, batch_size):
        filter_str = " or ".join([f"name startWith '{sp}'" for sp in species_batch])
        params = {"filter": filter_str}
        response = requests.get(BASE_URL, headers=ACCESS, params=params)
        
        # Debug: print the parameters being sent
        # print(params)
        
        # Possibility 1: API query failed
        if response.status_code != 200:
            for sp in species_batch:
                results.append([sp, "error", "NA"])
                summary["error"] += 1
                processed_count += 1
            if verbose:
                print(f"❌ API error for batch {species_batch}: {response.status_code}".ljust(120))
            print_progress(processed_count, total_species)
            continue

        data = response.json()
        items = data.get("items", [])

        # Debug: print the items received
        # print(items)

        for species in species_batch:

            # Filter items for exact matches, ignoring illegitimate or invalid names
            species_items = [
                item for item in items
                if item.get("name", "").strip().lower() == species.strip().lower()
                and item.get("nameStatus") not in ["Illegitimate", "Invalid"]
            ]

            # Possibility 2: no records or valid items found
            if not species_items:
                status = "no_valid" if any(
                    sp.lower() in (item.get("name", "").lower() for item in items) for sp in [species]
                    ) else "no_records"
                if verbose:
                    msg = "⚠️ No valid records" if status == "no_valid" else "⚠️ No records found"
                    print(f"{msg} for '{species}'".ljust(120))
                results.append([species, status, "NA"])
                summary[status] += 1
                processed_count += 1
                print_progress(processed_count, total_species)
                continue
            
            # Possibility 3: multiple records found, either identical or differing current names
            if len(species_items) > 1:
                current_ids = [item.get("synonymy", {}).get("currentNameId") for item in species_items]
                first_id = current_ids[0]
                if all(cid == first_id for cid in current_ids):
                    item = species_items[0]
                else:
                    options = [item.get("name") for item in species_items]
                    mb_numbers = [str(item.get("mycobankNr")) for item in species_items]
                    if verbose:
                        msg_lines = [f"⚠️ Multiple valid records with different current names found for '{species}':"]
                        for name, mb in zip(options, mb_numbers):
                            msg_lines.append(f"        {name} ({MB_URL_TEMPLATE}{mb})")
                        msg = "\n".join(msg_lines)
                        print(msg.ljust(120))
                    results.append([species, "multiple_records", "NA"])
                    summary["multiple_records"] += 1
                    processed_count += 1
                    print_progress(processed_count, total_species)
                    continue
            else:
                item = species_items[0]

            _id = item["id"]
            current_id = item.get("synonymy", {}).get("currentNameId")

            # Possibility 4: no current name ID found (never encountered this before)
            if not current_id:
                if verbose:
                    print(f"⚠️ No current name ID found for '{species}'".ljust(120))
                results.append([species, "no_current", "NA"])
                summary["no_current"] += 1
                processed_count += 1
                print_progress(processed_count, total_species)
                continue
            
            # Possibility 5: species name is current name
            if _id == current_id:
                status = "current"
                current_name = item["name"]
                summary["current"] += 1
                if verbose:
                    print(f"✅ '{species}' is the current MycoBank name.".ljust(120))
            
            # Possibility 6: species name is not current
            else:
                status = "not_current"
                response_current = requests.get(f"{BASE_URL}/{current_id}", headers=ACCESS)
                data_current = response_current.json()
                current_name = data_current.get("name", "<unknown>")
                mycobank_number = data_current.get("mycobankNr", "<unknown>")
                summary["not_current"] += 1
                if verbose:
                    print(f"🔄 The current MycoBank name for '{species}' is '{current_name}' "
                               f"with MycoBank number '{mycobank_number}': \n        "
                               f"{MB_URL_TEMPLATE}{mycobank_number}".ljust(120))

            results.append([species, status, current_name])
            processed_count += 1
            print_progress(processed_count, total_species)

    # Write results to TSV
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["species_query", "status", "current_name"])
        writer.writerows(results)
    if verbose:
        print(f"\n\nResults written to {output_file}\n")

    return summary


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        description="Get the current MycoBank name for a list of species."
        )
    parser.add_argument(
        "species_list", 
        help="Path to the species list txt file (one species per line)")
    parser.add_argument(
        "--exclude", 
        metavar="FILE_PATH", 
        help="Path to a txt file containing species names to exclude from query (one species per line)", 
        default=None
        )
    parser.add_argument(
        "--batch-size", 
        help="Number of species per API request", 
        type=int, default=20
        )
    parser.add_argument(
        "--output", 
        metavar="OUTPUT_PATH", 
        help="Path to write MycoBank query output TSV file", 
        default="mycobank_results.tsv"
        )
    parser.add_argument(
        "--excluded-output", 
        metavar="EXCLUDED_OUTPUT_PATH", 
        help="Path to write TSV file for excluded species", 
        default="excluded_species.tsv"
        )
    parser.add_argument(
        "--verbose", 
        help="Print detailed output", 
        action="store_true"
        )
    
    # If no arguments are provided, show help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    excluded_species = []
    if args.exclude:
        with open(args.exclude) as f:
            excluded_species = [line.strip() for line in f if line.strip()]

    species_to_check, excluded_user_count, excluded_ambig_count = read_species_list(
        args.species_list,
        excluded_species=excluded_species,
        excluded_output_file=args.excluded_output, 
        verbose=args.verbose,
    )

    summary = get_current_names_batch(
        species_to_check,
        output_file=args.output,
        batch_size=args.batch_size,
        verbose=args.verbose
    )

    # Final summary
    total = sum(1 for line in open(args.species_list) if line.strip())
    print("\n\n=========== FINAL SUMMARY ===========")
    print(f"Total species in input file:      {total}")
    print(f"Excluded (user list):             {excluded_user_count}")
    print(f"Excluded (ambiguous):             {excluded_ambig_count}")
    print(f"Queried:                          {len(species_to_check)}\n")

    print("Query Results:")
    for key, count in summary.items():
        print(f"{key:>20}: {count}")
    print("=====================================\n")


if __name__ == "__main__":
    main()
