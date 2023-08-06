# scHi-CNN

This pipeline consists of three main scripts: `knear.py`, `normalize.py`, and `significant-filter.py`. They are designed to work in sequence to produce significant intra-chromatin interactions from single cells.

## Requirements

- Python 3.7+
- Libraries: `numpy`, `traceback`, `scipy`, `mpi4py`, `cooler`, `h5py`
- MPI Runtime Environment for parallel computing (e.g., OpenMPI) is highly recommended for reduced running time
- Input: A directory containing `.cool` files. A script `convert-cool.sh` is provided to convert paired end files to `.cool` files.
- To generate `.cool` files from your Hi-C data, you may use the `cooler` library which provides CLI tools for creating, managing, and manipulating `.cool` files. For more information, please refer to the [cooler documentation](https://cooler.readthedocs.io/en/latest/).

## Input File Format

The pipeline accepts input files in two formats: tab-separated value files(gzipped) or `.cool` files. 

### Tab-separated value files
The file should have the following columns:

| Cell name | chr1   | coord1     | chr2  | coord2     |
|-----------|--------|------------|-------|------------|
| 2_B7_AD004| chr10  | 100011044  | chr10 | 100011145  |
| 2_B7_AD004| chr10  | 100020824  | chr10 | 100020610  |
| 2_B7_AD004| chr10  | 100026356  | chr10 | 100026377  |
| 2_B7_AD004| chr10  | 100037257  | chr4  | 6366322    |
| 2_B7_AD004| chr10  | 100037707  | chr10 | 100037871  |
| 2_B7_AD004| chr10  | 100042086  | chr10 | 100041852  |

Each column corresponds to:

- `Cell name`: The name of the cell.
- `chr1`: The chromosome of the first coordinate.
- `coord1`: The first coordinate.
- `chr2`: The chromosome of the second coordinate.
- `coord2`: The second coordinate.

### COOL Format

The pipeline also accepts `.cool` files as input. A `.cool` file is a format for storing genomic interaction data that includes chromosomal coordinates and the associated interaction data.

Please ensure the input files are formatted correctly to avoid errors during processing. All the cell files should be contained in one directory.

## Script Descriptions

### 1. `knear.py`

This script is designed to run in an MPI environment using multiple processes (`mpirun -np`). It accepts following arguments.

- `resol`: Resolution (required, integer).
- `cutoff`: maximum spatial distance to be imputed (required, integer).
- `pad-local`: Padding for local area (optional, integer, default is 5).
- `topk`: The top k neighbours to consider for imputation (optional, integer, default is 4).
- `type`: Cell type (required).
- `prefix`: Name for the output directory (required).
- `out-suffix`: Suffix for the output sub directory (required).
- `cell-file`: Path to the file containing cell file names (required).
- `genome`: Genome name, currently supported: hg19,mm10 (required).

### 2. `normalize.py`

This script is used for normalizing the data.

- `type`: Cell type (required).
- `prefix`: Name for the output directory (required).
- `res`: Resolution (optional, integer, default is 100000).
- `in-suffix`: Suffix for the input directory (required).
- `out-suffix`: Suffix for the output directory (required).
- `cell-file`: Path to the file containing cell file names (required).
- `genome`: Genome name, currently supported: hg19,mm10 (required).
- `max-dist`: Maximum distance to consider for normalization (optional, float, default is 2e6).

### 3. `significant-filter.py`

This script is used to filter the data based on significance.

- `type`: Cell type (required).
- `genome`: Genome name, currently supported: hg19,mm10 (required).
- `res`: Resolution (optional, integer, default is 100000).
- `prefix`: Name for the output directory (required).
- `in-suffix`: Suffix for the input directory (required).
- `center-dist`: Bin distance to consider for paired t-test (optional, integer, default is 2).
- `max-dist`: Maximum distance to consider for imputation (optional, float, default is 2e6).
- `cell-file`: Path to the file containing cell file names (required).
- `fdr`: False discovery rate threshold (optional, float, default is 0.1).
- `tstat`: T-statistic threshold (optional, float, default is 3.0).

## Running the Pipeline

Before running, make sure you have all necessary inputs and parameters. All three scripts need to be executed in the order: `knear.py` > `normalize.py` > `significant-filter.py`.

- Input: A directory containing `.cool` files. A script `convert-cool.sh` is provided to convert paired end files to `.cool` files.
- Example paired end raw input data is provided in the `sample-data` directory.
- Example .cool files for the corresponding input data are provided in the `sample-cool-out` directory.
- Example pipeline is provided in the `job-scHi-CNN.pbs` file.
- Example output file is provided in the `sample-output` directory. (intermediate directories and files are not included). Please delete or rename the current output directory before running the pipeline.
- Example final output contains in the `sample-significant.txt` file.

## Final Output File Format

The final output of the pipeline is a table with the following columns. Each row represents a significant chromatin interaction.

| chr1   | start1  | end1    | chr2  | start2  | end2    | tstat  | pvalue  | qvalue  |
|--------|---------|---------|-------|---------|---------|--------|---------|---------|
| chr1   | 2400000 | 2500000 | chr1  | 2700000 | 2800000 | 4.11554| 0.00262 | 0.02441 |
| chr1   | 2900000 | 3000000 | chr1  | 3200000 | 3300000 | 3.68953| 0.00500 | 0.03247 |
| chr1   | 8000000 | 8100000 | chr1  | 8300000 | 8400000 | 3.32910| 0.00881 | 0.04141 |
| chr1   | 8100000 | 8200000 | chr1  | 8400000 | 8500000 | 3.93407| 0.00344 | 0.02961 |
| chr1   | 8500000 | 8600000 | chr1  | 8800000 | 8900000 | 3.42525| 0.00756 | 0.03917 |
| chr1   | 8600000 | 8700000 | chr1  | 8900000 | 9000000 | 4.08187| 0.00275 | 0.02464 |

Each column corresponds to:

- `chr1`: Chromosome of the first region.
- `start1`: Start coordinate of the first region.
- `end1`: End coordinate of the first region.
- `chr2`: Chromosome of the second region.
- `start2`: Start coordinate of the second region.
- `end2`: End coordinate of the second region.
- `tstat`: The t-statistic value.
- `pvalue`: The p-value of the test.
- `qvalue`: The q-value (adjusted p-value) of the test.

Please note that the coordinates are 0-based, and the end coordinate is exclusive. This is consistent with the BED file format commonly used in bioinformatics.

## Contact
For any issues or concerns, please open an issue on this GitHub repository or contact the repository owner directly.

Note: This pipeline is for research purposes only. Please cite the relevant work if you use it in your published research.
