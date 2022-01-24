# vgaligner-new-experiments
New experiments for rs-vgaligner](https://github.com/HopedWall/rs-vgaligner). 

The main purpose of this repo is to compare the following aligners:
- rs-vgaligner
- [vg map](https://github.com/vgteam/vg/blob/master/src/subcommand/map_main.cpp)
- [vg giraffe](https://github.com/vgteam/vg/blob/master/src/subcommand/giraffe_main.cpp)
- [graphaligner](https://github.com/maickrau/GraphAligner)
- [graphchainer](https://github.com/algbio/GraphChainer)

## Requirements
First, [snakemake](https://snakemake.readthedocs.io/en/stable/#) must be installed. The aligners listed above must be already installed and available on PATH.

Other requirements include:
- [pbsim](https://github.com/vgteam/vg) for generating the reads
- [odgi](https://github.com/pangenome/odgi) for sorting the graphs
- python3 with pandas installed
- awk

Using conda to handle requirements is a TODO.

## How to run the experiments
In order to execute the pipeline for a specific graph, move into the folder for a specific graph and run:

`snakemake --keep-going --cores 8`

`--keep-going` is required because vgaligner can crash on larger graphs... yet another TODO

## How validation works
For each graph, do the following:
1. Obtain the paths, store each path as a separate FASTA file
2. For each path, simulate reads with PBSIM
3. Map the reads to the graph (repeat for each aligner)
4. Check the quality of the resulting alignments

Step 1-3 are fairly straightforward, but step 4 is where problems start to arise: what is a correct alignment? How can we measure "correctness"?

For each path in the graph, we know the nodes it is comprised of. Since we generate reads for each path separately, we also have a general idea of where (i.e. to which nodes) they should be mapped to i.e. the "ground truth".

The alignments are also returned as paths on the graph. Then we compute the following "overlap": nodes_in_truth_length / total_alignment_length. We consider a read to be correctly mapped if the "overlap" is over a certain threshold (i.e. 0.8).

## Datasets
TODO
