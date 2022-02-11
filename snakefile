#configfile: "config.yaml"
#datasets = config["datasets"].split(",")
#kmersize = config["kmersize"]
#threshold = config["threshold"]

import os
DATASETS = [name for name in os.listdir(".") if os.path.isdir(name) and not name.startswith(".")]
DATASETS.remove("scripts")

'''
PATHS_FOR_DATASET = dict()
for dataset in DATASETS:
	graph_filename = os.path.join(".", dataset, "graph.gfa")
	paths = os.popen("odgi paths -i {} -L".format(graph_filename)).read().split("\n")[:-1]
	PATHS_FOR_DATASET[dataset] = paths
'''

FINAL_DATASETS = []
PATHS = []
for dataset in DATASETS:
	graph_filename = os.path.join(".", dataset, "graph.gfa")
	paths = os.popen("odgi paths -i {} -L".format(graph_filename)).read().split("\n")[:-1]
	FINAL_DATASETS.extend([dataset]*len(paths))
	PATHS.extend(paths)

assert(len(FINAL_DATASETS) == len(PATHS))
#print(FINAL_DATASETS)
#print(PATHS)

kmersize = 11	# VG can't build GCSA2 index with kmer size > 16 (i.e. can't use 21 for comparisons)
threshold = 0.8

rule all:
	input:
		expand("{dataset}/comparisons/{path}_graphaligner.txt",zip,dataset=FINAL_DATASETS,path=PATHS),
		expand("{dataset}/comparisons/{path}_vgaligner.txt",zip,dataset=FINAL_DATASETS, path=PATHS),
		expand("{dataset}/comparisons/{path}_vgmap.txt",zip,dataset=FINAL_DATASETS, path=PATHS),
		expand("{dataset}/comparisons/{path}_giraffe.txt",zip,dataset=FINAL_DATASETS,path=PATHS),
		expand("{dataset}/comparisons/{path}_graphchainer.txt",zip,dataset=FINAL_DATASETS, path=PATHS)

##### Sort the graph #####
rule odgi_sort:
	input:
		"{dataset}/graph.gfa"
	output:
		"{dataset}/sorted_starting_graph.gfa"
	threads: 8
	shell:
		"odgi sort -i {input} -o - -p Ygs -P -t {threads} | odgi view -i - -g -t {threads} > {output}"

#### Graph stats ####
rule graph_stats:
    input:
        "{dataset}/sorted_graph_{graph}.gfa"
    output:
        "{dataset}/stats/stats_{graph}.txt"
    shell:
        "vg stats -zAL {input} > {output}"

##### Split into connected components #####
rule odgi_explode:
	input:
		"{dataset}/sorted_starting_graph.gfa"
	output:
		"{dataset}/component.{i}.og"
	threads: 8
	shell:
		"odgi explode -i {input} -t {threads} -p {dataset}/component"

##### Convert .og into .gfa #####
rule odgi_convert:
	input:
		"{dataset}/component.0.og"
	output:
		"{dataset}/sorted_graph.gfa"
	shell:
		"odgi view -i {input} -g > {output}"

##### Create a fasta file with all the paths #####
rule get_paths:
	input:
		"{dataset}/sorted_graph.gfa"
	output:
		"{dataset}/paths.fasta"
	shell:
		"odgi paths -i {input} -f > {output}"

##### Split fasta into separate files #####
rule split_paths:
	input:
		"{dataset}/paths.fasta"
	output:
		expand("{{dataset}}/paths/{path}.fa", path=PATHS)
	shell:
		"awk -v prefix='{dataset}/paths/' -f scripts/split.awk {input}"

##### Read generation with pbsim #####
rule pbsim_generate:
	input:
		"{dataset}/paths/{path}.fa"
	output:
		"{dataset}/reads/{path}_0001.fastq"
	shell:
		"pbsim '{input}' --model_qc scripts/pbsim-models/model_qc_clr --prefix '{dataset}/reads/{wildcards.path}'"

##### Run vgaligner #####
rule vgaligner_index:
	input:
		"{dataset}/sorted_graph.gfa"
	output:
		index = "{dataset}/sorted_graph.idx",
		mappings = "{dataset}/mappings.json"
	threads: 8
	log:
		log = "{dataset}/logs/vgaligner-index.log",
		time = "{dataset}/logs/vgaligner-index.time"
	shell:
		'''
		/usr/bin/time -v -o '{log.time}' vgaligner index -i {input} -k {kmersize} -t {threads} -o {output.index} --generate-mappings > {log.log}
		mv mappings.json {dataset}/mappings.json
		'''

rule vgaligner_map:
	input:
		index = "{dataset}/sorted_graph.idx",
		reads = "{dataset}/reads/{path}_0001.fastq"
	output:
		"{dataset}/results/{path}_vgaligner.gaf"
	threads: 8
	log:
		log = "{dataset}/logs/{path}_vgaligner-map.log",
		time = "{dataset}/logs/{path}_vgaligner-map.time"
	shell:
		"/usr/bin/time -v -o '{log.time}' vgaligner map -i {input.index} -f '{input.reads}' --also-align -t {threads} -o '{output}' > '{log.log}'"

##### Run VG #####
rule vg_convert_gfa_to_vg:
	input:
		"{dataset}/sorted_graph.gfa"
	output:
		"{dataset}/sorted_graph_temp.vg"
	shell:
		"vg convert -g {input} > {output}"

rule vg_mod:
	input:
		"{dataset}/sorted_graph_temp.vg"
	output:
		"{dataset}/sorted_graph.vg"
	shell:
		"vg mod -X 256 {input} > {output}"

rule vg_index:
	input: 
		"{dataset}/sorted_graph.vg"
	output:
		xg="{dataset}/sorted_graph.xg",
		gcsa="{dataset}/sorted_graph.gcsa"
	threads: 8
	log:
		time="{dataset}/logs/vgindex.time"
	shell:
		"/usr/bin/time -v -o {log.time} vg index {input} -x {output.xg} -g {output.gcsa} -k {kmersize} -t {threads}"

rule vg_map:
	input:
		"{dataset}/sorted_graph.xg",
		"{dataset}/sorted_graph.gcsa",
		reads="{dataset}/reads/{path}_0001.fastq"
	output:
		temporary("{dataset}/results/{path}_vgmap.gam")
	params:
		prefix="{dataset}/sorted_graph"
	threads: 8
	log:
		time="{dataset}/logs/{path}_vgmap.time"
	shell:
		"/usr/bin/time -v -o '{log.time}' vg map -d {params.prefix} -f '{input.reads}' -t {threads} > '{output}'"

rule vg_map_convert:
	input: 
		gam="{dataset}/results/{path}_vgmap.gam",
		graph="{dataset}/sorted_graph.vg"
	output:
		"{dataset}/results/{path}_vgmap.gaf"
	shell:
		"vg convert --gam-to-gaf '{input.gam}' {input.graph} >'{output}'"

##### Run graphaligner #####
rule graphaligner:
	input:
		graph="{dataset}/sorted_graph.gfa",
		reads="{dataset}/reads/{path}_0001.fastq"
	output:
		"{dataset}/results/{path}_graphaligner.gaf"
	threads: 8
	log:
		log="{dataset}/logs/{path}_graphaligner.log",
		time="{dataset}/logs/{path}_graphaligner.time"
	shell:
		"touch '{output}' && /usr/bin/time -v -o '{log.time}' GraphAligner -g {input.graph} -f '{input.reads}' -a '{output}' -x vg -t {threads} > '{log.log}'"

##### Run giraffe #####
# https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe

rule vg_autoindex:
	input:
		"{dataset}/sorted_graph.gfa"
	output:
		"{dataset}/index.dist",
		"{dataset}/index.gg",
		"{dataset}/index.giraffe.gbwt",
		"{dataset}/index.min"
	params:
		prefix="{dataset}/index"
	shell:
		"vg autoindex --workflow giraffe -g {input} -p {params.prefix}"

rule vg_giraffe:
	input:
		gbwt_index="{dataset}/index.giraffe.gbwt",
		gbwt_graph="{dataset}/index.gg",
		minimizer_index="{dataset}/index.min",
		dist_index="{dataset}/index.dist",
		reads="{dataset}/reads/{path}_0001.fastq"
	output:
		"{dataset}/results/{path}_giraffe.gaf"
	threads: 8
	shell:
		"vg giraffe -H {input.gbwt_index} -g {input.gbwt_graph} -m {input.minimizer_index} -d {input.dist_index} -f '{input.reads}' -t {threads} -o gaf > '{output}'"

##### Run graphchainer #####
rule graphchainer:
	input:
		graph="{dataset}/sorted_graph.gfa",
		reads="{dataset}/reads/{path}_0001.fastq"
	output:
		"{dataset}/results/{path}_graphchainer.gaf"
	threads: 8
	log:
		log="{dataset}/logs/{path}_graphchainer.log",
		time="{dataset}/logs/{path}_graphchainer.time"
	shell:
		"touch '{output}' && /usr/bin/time -v -o '{log.time}' GraphChainer -g {input.graph} -f '{input.reads}' -a '{output}' -t {threads} > '{log.log}'"

rule vg_map_convert_graphchainer:
	input: 
		gam="{dataset}/results/{path}_graphchainer.gam",
		graph="{dataset}/sorted_graph.gfa"
	output:
		"{dataset}/results/{path}_graphchainer.gaf"
	shell:
		"vg convert --gam-to-gaf '{input.gam}' {input.graph} >'{output}'"

##### Run comparison scripts #####
rule gamcompare_graphaligner:
	input:
		gaf="{dataset}/results/{path}_graphaligner.gaf",
		mappings="{dataset}/mappings.json"
	output:
		"{dataset}/comparisons/{path}_graphaligner.txt"
	shell:
		"python3 scripts/GAFMAFcomparison2.py '{input.gaf}' {input.mappings} '{wildcards.path}' {threshold} graphaligner > '{output}'"

rule gamcompare_vgaligner:
	input:
		gaf="{dataset}/results/{path}_vgaligner.gaf",
		mappings="{dataset}/mappings.json"
	output:
		"{dataset}/comparisons/{path}_vgaligner.txt"
	shell:
		"python3 scripts/GAFMAFcomparison2.py '{input.gaf}' {input.mappings} '{wildcards.path}' {threshold} vgaligner > '{output}'"

rule gamcompare_vgmap:
	input:
		gaf="{dataset}/results/{path}_vgmap.gaf",
		mappings="{dataset}/mappings.json"
	output:
		"{dataset}/comparisons/{path}_vgmap.txt"
	shell:
		"python3 scripts/GAFMAFcomparison2.py '{input.gaf}' {input.mappings} '{wildcards.path}' {threshold} graphaligner > '{output}'"

rule gamcompare_giraffe:
	input:
		gaf="{dataset}/results/{path}_giraffe.gaf",
		mappings="{dataset}/mappings.json"
	output:
		"{dataset}/comparisons/{path}_giraffe.txt"
	shell:
		"python3 scripts/GAFMAFcomparison2.py '{input.gaf}' {input.mappings} '{wildcards.path}' {threshold} graphaligner > '{output}'"

rule gamcompare_graphchainer:
	input:
		gaf="{dataset}/results/{path}_graphchainer.gaf",
		mappings="{dataset}/mappings.json"
	output:
		"{dataset}/comparisons/{path}_graphchainer.txt"
	shell:
		"python3 scripts/GAFMAFcomparison2.py '{input.gaf}' {input.mappings} '{wildcards.path}' {threshold} graphaligner > '{output}'"