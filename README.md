A utility to determine ROH based on runs of heterozygous or homozygous calls in arbitrary VCFs.

Requirements:
- Linux Ubunutu 22.04 (tested)
- Python3.6 (tested)
- see requirements.txt for python library versions


# Instructions
- Create a python venv and install required components
- Run the script like so:
  roh-bins-genome.py vcf1 vcf2 vcf3
- The script will process a maximum of three vcfs.
- Run time is approximately 2 minutes per 3x VCFs (for demo data)
- Output is a pdf with ROH indicated in each plot (see example .png).
