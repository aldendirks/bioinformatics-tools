# Bioinformatics Utilities, Scripts, and Programs

## Scripts

* `fetch_inat_seqs.py` - fetch ITS sequences from iNaturalist for a user-provided taxon ID. 
* `get_current_name.py` - with a user-provided list of species names, queries MycoBank for the current name. 
* `print_seqs.py` - print one or more sequences from a multifasta file based on sequence position. 
* `remove_seqs.py` - remove sequences from a multifasta file that are shorter than a specified length.

## ITS phylogenetics pipeline

These scripts can be strung together for a semi-automatic ITS phylogenetic analysis pipeline. The example belwo shows a pipeline to create an ITS phylogeny for the jelly fungus genus *Phaeotremella*. This genus has an ID of 994028 on iNaturalist. 

Create directory structure. 

```
mkdir -p aln seqs phylo
```

Fetch sequences from iNaturalist and GenBank. 

```
fetch_inat_seqs.py 994028 --output seqs/inat.fasta
# fetch_genbank_seqs.py ...
```

Find ITS sequences that are in the wrong orientation and reverse complement them. See Alan Rockefeller's [`fixfasta.py`](https://github.com/AlanRockefeller/fixfasta.py) script. 

```
fixfasta.py seqs/its_formatted.fasta seqs/genbank_formatted.fasta > seqs/combined.fasta
```

Align and trim sequences. **I highly recommend manually reviewing the alignment and correcting or deleting any misformatted or low-quality sequences**

```
mafft --auto seqs/combined.fasta > aln/combined.aln
trimal -in aln/combined.aln -out aln/trimmed.aln -fasta -automated1
```

Phylogenetic analysis.

```
iqtree3 -s aln/combined.aln --alrt 1000 -B 1000 -pre phylo/combined -redo
iqtree3 -s aln/trimmed.aln --alrt 1000 -B 1000 -pre phylo/trimmed -redo
```

Print formatted ITS trees. 

```
# Rscript ...
```