## Needle

### A fast and space-efficient pre-filter for estimating the quantification of very large collections of nucleotide sequences

[![Build Status](https://github.com/seqan/app-template/workflows/App%20CI/badge.svg)](https://github.com/seqan/needle/actions?query=branch%3Amaster+workflow%3A%22App+CI%22) [![codecov](https://codecov.io/gh/seqan/needle/branch/master/graph/badge.svg?token=SJVMYRUKW2)](https://codecov.io/gh/seqan/needle)  [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](#install-with-bioconda-linux)

Needle is a tool for semi-quantitative analysis of very large collections of nucleotide sequences.
Needle stores its data in multiple interleaved Bloom filter, a fast and space efficient probabilistic data structure and uses a windowing scheme (also called minimisers) to reduce the amount of data to store. How many interleaved Bloom filter are used is defined by the user. Each interleaved Bloom filter has a so called expression threshold and stores minimisers with an occurrence greater than or equal to its own expression threshold and smaller than the next biggest expression threshold (if there is no bigger expression threshold, all greater than or equal to the threshold are stored). These expression thresholds are then used during the query (called estimate) to approximate the expression values of given transcripts.

## Citation

Please cite:

Mitra Darvish, Enrico Seiler, Svenja Mehringer, René Rahn, Knut Reinert, Needle: a fast and space-efficient prefilter for estimating the quantification of very large collections of expression experiments, Bioinformatics, Volume 38, Issue 17, 1 September 2022, Pages 4100–4108, https://doi.org/10.1093/bioinformatics/btac492

## Download, Install & Build

<details><summary>Prerequisites (click to expand)</summary>

* CMake >= 3.10
* GCC 10, 11 or 12 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth
information.
</details>

### Install with [bioconda](https://bioconda.github.io/recipes/needle/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge needle
```

### Install via github

Needle can be built by following these commands:

```
git clone --recurse-submodules https://github.com/seqan/needle.git
mkdir build-needle && cd build-needle
cmake ../needle
make
```

Run test to check, if Needle is working as intended. All tests should pass.

```
make test
```

If you are interested in building the documentation, just use the command: `make doc`

## Build an Needle index
In order to build a Needle index a number of sequence files have to be given. All sequence file formats supported by seqan3 are accepted as an input (fasta, fastq, embl,... and their compressed forms). The flag `--paired` in the example below indicates that the given sequence files are paired-end experiments. Furthermore, the false positive rate has to be specified with the parameter `f`.
Use -h/--help for more information and to see further parameters. The flag `-c` can be used to build a compressed Needle index.

The following example creates a compressed Needle index for two paired-end experiments for the expression thresholds 4 and 32.

```
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
```

Although, this works. It is recommended to calculate the minimisers beforehand by using the option `minimisers`. It calculates the minimisers of given experiments and stores their hash values and their occurrences in a binary file named ".minimiser".

The following command calculates the minimisers in the two experiments.
```
./bin/needle minimiser ../needle/test/data/exp_*.fasta --paired
```

A minimiser file is a binary file containing the following data:
- number of minimisers (uint64_t)
- kmer-size (uint8_t)
- window-size (uint32_t)
- seed (uint64_t)
- flag which is true, if shape is ungapped (bool)
- shape (uint64_t), if flag is false
- all minimiser hashes (uint64_t) with their occurrences (uint16_t)

Based on the minimiser files the Needle index can be computed by using the following command:
```
./bin/needle ibfmin exp*.minimiser -e 16 -e 32  -f 0.3 -c -o example
```

## Estimate
To estimate the expression value of one transcript a sequence file has to be given. Use the parameter "-i" to define where the Needle index can be found (should be equal with "-o" in the previous commands).
Use -h/--help for more information and to see further parameters.
The following example searches for one gene, which is expressed in the first experiment with expression 6 and in the second with expression 37. Therefore, it should be found only in the second experiment but not the first when using expression levels of 16 and 32.

```
./bin/needle estimate ../needle/test/data/gene.fasta -i example
```

The created file "expressions.out" (if you prefer a different name, use "-o") should contain the following:
```
GeneA   0      32
```

## Insert into an existing Needle index
It is possible to insert new sequence files into an uncompressed Needle index. Similar to the build step this can be done by either using the sequence files as input directly or the minimiser files outputed by `needle minimiser`. Most options are the same as the ones from the build step, however as the Needle index already exist, neither the false positive rate nor the number of hash functions can be changed. It is necessary to specify `i` to the directory, where the existing Needle index can be found.

The following example inserts into the Needle index build above for two paired-end experiments.

```
./bin/needle ibf ../needle/test/data/exp_0*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example // Create Index
./bin/needle insert ../needle/test/data/exp_1*.fasta --paired -i example                    // Insert into created index
```

Based on minimiser files an insertion to the Needle index can be achieved by using the following command:
```
./bin/needle ibf ../needle/test/data/exp_0*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example // Create Index
./bin/needle insertmin exp*.minimiser -i example                                             // Insert into created index
```

The insert methods based on minimiser or on sequence files is independent of the way the index was created.

## Delete experiments from an existing Needle index
It is possible to delete sequence files from an uncompressed Needle index by specifiying the position of the experiment, which should be deleted.
These deleted experiments won't change the size of the index as the space is kept for later insertions.
```
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example // Create Index
./bin/needle delete  -i example 0                                            // Delete first experiment exp_0 (with position 0) from index
``



## Note

This app was created with the [seqan3 app-template](https://github.com/seqan/app-template).
