# gene_info.py - get information about the number of introns and the length of genes from the GFF3 file
## Installation
Before using the script, [PyInterval](https://github.com/taschini/pyinterval/) needs to be installed:
```bash
pip install --upgrade pip
pip install pyinterval
```
Then, download or clone this repository, and you are good to go.
```bash
git clone https://github.com/yangwu91/gene_info.git
```
## Usage
It is very straightforward:
```bash
python gene_info.py <GFF3>
```
The output file will be named as `<GFF3>_gene_structure.csv`, which could be opened with Excel or Notepad etc.