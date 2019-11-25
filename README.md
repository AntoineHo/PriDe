# PriDe
PriDe: Primer Design tool

## Installation
```bash
conda create -n pride --file packages.txt
conda activate pride
```

## Usage
- Print usage and help

`python PriDe.py {design, test} -h`

### Design primers and test for cross matches
- Design primers with primer3:

`python PriDe.py design reference.fa regions.bed output_directory --offset 50 --processes 10`

- Blast primers on the reference:

`python PriDe.py test output_directory/primers.fasta reference.fa output_directory --processes 10`

- Grep perfect matches only

`cat design_out/primers.fasta.blast.tsv | awk '{if($6 == $7) print}'`
