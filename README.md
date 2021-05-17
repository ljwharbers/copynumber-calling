# Workflow for (sc)CUTseq preprocessing and copy number calling

## Preprocessing
### Dependencies
Python:
* argparse
* gzip
* re
* os
* multiprocessing
* functools
* collections
* itertools
* operator
* Biopython
* distance
* glob
* pandas

Other software:
* snakemake
* alfred

### Workflow
This workflow expects standard illumina demultiplexed `fastq` files as input (can be either multiple `fastq` files for each lane or merged).
To process the a single multiplex library edit the `config.yaml` file in the `preprocessing` folder with the appropriate values.

A reference genome, the link to the `alfred` binary, a link the the `processBarcode.py` and a list of barcodes are also required. The list of barcodes should be a standard text file with each row containing one barcode. 

The `barcode pattern` depends on the sequencing method used. For standard CUTseq (96 linkers) this is typically `UUUUUUUUBBBBBBBBDDDD` and for scCUTseq (384 linkers) this is typically `UUUUUUUUBBBBBBBBBBBDDDD`.

After successfully editing the `config.yaml` you can run the entire snakemake file from the command line. If you are in the `preprocessing` directory this can be done in the following way: `snakemake -j {number of threads to use}`
