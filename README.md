# Preprocessing of yeast dataset for Kotoura method

## Usage

```bash
python3 yeapreproc.py <CSV file> <name of column condition>
```

Result will be a file named like `processed_<CSV file>`

CSV file is parsed using `pandas` module. The parsed file will have a two-level header, reflecting the raw CSV file header structure. 

- `<name of column condition>` is the name of a first level header you want to select as the condition. 
- The control condition is always looked up in the columns having a first level header named `0perc`, as per the initial CSV file structure.

The result file header is flattened as a one level header only to allow easy loading using C/C++ libraries.

## Example 

```
python3 yeapreproc.py yeast_ko.csv 2perc
```
Computes deltas between 2% ethanol condition and 0% (control)

## Requirements

```
pandas
python-libsbml
```
### Installation
```
pip3 install -r requirements.txt
```
