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

- Design primers with primer3:

`python PriDe.py design reference.fa regions.bed output_directory --offset 50 --processes 10`

- Blast primers on the reference:

`python PriDe.py test primers.fa reference.fa output_directory --processes 10`
