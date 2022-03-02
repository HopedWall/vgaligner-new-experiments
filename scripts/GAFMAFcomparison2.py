import argparse
import pandas as pd
import re
import json
import sys

# Parse MAF to a dataframe
def parse_MAF(path_to_file, genome_name):
        f = open(path_to_file)

        ref_data = []
        read_data = []
        for line in f.readlines():
                if line.startswith("s {}".format(genome_name)):
                        s = line.split()
                        new_row = {"seqname": s[1], "strand": s[4], "start": s[2], "len": s[3]}
                        ref_data.append(new_row)
                elif line.startswith("s S"):
                        s = line.split()
                        new_row = {"seqname": s[1], "strand": s[4], "start": s[2], "len": s[3]}
                        read_data.append(new_row)

        ref_df = pd.DataFrame(ref_data, columns=["seqname", "strand", "start", "len"])
        read_df = pd.DataFrame(read_data, columns=["seqname", "strand", "start", "len"])

        return ref_df, read_df

# Parse GAF to a dataframe
def parse_GAF(path_to_file, tool):
        gaf_fields = {
                'vgaligner' : ["name", "qlen", "qstart", "qend", 
                                "strand", "path", "plen", "pstart", "pend", 
                                "residue", "alblock", "quality", "extra",
                                "extra1", "extra2"],
                'graphaligner' : ["name", "qlen", "qstart", "qend", 
                                "strand", "path", "plen", "pstart", "pend", 
                                "residue", "alblock", "quality", "extra",
                                "extra1", "extra2", "extra3", "extra4"],
                'vgmap' : ["name", "qlen", "qstart", "qend", 
                        "strand", "path", "plen", "pstart", "pend", 
                        "residue", "alblock", "quality", "extra",
                        "extra1", "extra2"]
        }
        my_gaf = pd.read_csv(path_to_file, sep='\t', names=gaf_fields[tool]) 
        return my_gaf



# CLI arguments parsing
parser = argparse.ArgumentParser(description='Convert paths in a GAF to ranges')
parser.add_argument('GAF', help='Path to the GAF file')
parser.add_argument('Mappings', help='Path to the mappings')
#parser.add_argument('MAF', help='Path to the MAF file')
parser.add_argument("Genome-name", help="Name (as a string) of the genome the reads are from")
parser.add_argument("Threshold", help="Threshold that should be used to consider a read mapped correctly")
parser.add_argument('Tool', help='Tool the GAF was obtained from (vgaligner, graphaligner, vgmap)')
args = vars(parser.parse_args())

print("Command used: {}".format(" ".join(sys.argv)))

print("GAF file: {}".format(args["GAF"]))
# Read GAF file
my_gaf = parse_GAF(args["GAF"], args["Tool"])

# Which path (in the original graph) should be used to perform the comparison
#path = "gi|568815569:3979127-3993865"
path = args["Genome-name"]

# Print out some debug information
print("Chosen path is: {}".format(path))
print("Threshold is: {}\n".format(args["Threshold"]))

# Read JSON with mappings
f = open(args["Mappings"])
mappings = json.load(f)

# This should be false only if there's a typo in the chosen path
if path in mappings.keys():

        # Get the mappings for the specific path
        mappings_for_path = mappings[path]

        # Read MAF file
        #ref_df, read_df = parse_MAF(args["MAF"], path)
        
        # This is ok, disable for debugging purposes
        #n_reads = len(read_df.index)
        n_reads = len(my_gaf.index)
        
        reads_mapped_correctly = 0

        # Find the nodes for the current path
        nodes_as_int = []
        for node in mappings_for_path.keys():
                nodes_as_int.append(int(node)) 
        nodes_as_int.sort()
        print("Nodes for current path: {}".format(nodes_as_int))

        # For each GAF record
        for i in range(len(my_gaf)):

                # Get the path it encodes
                my_gaf_row = my_gaf.iloc[i]
                my_gaf_name = my_gaf_row["name"]
                my_gaf_nodes_str = my_gaf_row["path"]

                if my_gaf_nodes_str == '*':
                        n_reads -= 1
                        continue

                my_gaf_tuples = re.findall("(>|<)([0-9]+)", my_gaf_nodes_str)
                
                #my_gaf_int = list(map(lambda x: +int(x[1]) if x[0]=='>' else -int(x[1]), my_gaf_tuples))
                my_gaf_int = list(map(lambda x: int(x[1]), my_gaf_tuples)) # strand should not matter

                my_gaf_path_len = int(my_gaf_row["plen"])

                #print(my_gaf_int)

                # Check if alignment is on reverse strand
                is_rev = any(map(lambda x: x < 0, my_gaf_int))
                # TODO: handle reverse strand... 
                
                # Iterate over each node in the path
                correct_alignment_length = 0
                for i in range(0, len(my_gaf_int)):
                        node_as_str = str(my_gaf_int[i])
                        if node_as_str in mappings_for_path.keys():
                                
                                # Get the mappings for the current node
                                curr_int = mappings_for_path[node_as_str]
                                
                                # Find start and end for the current node
                                interval_start = curr_int["start"]
                                interval_end = curr_int["end"]

                                # Update start (end) if first (last) node
                                '''
                                if i == 0:
                                        # Add pstart to left border
                                        interval_start += int(my_gaf_row["pstart"])
                                elif i==len(my_gaf_int)-1:
                                        # Remove pend from right border
                                        interval_end -= int(my_gaf_row["pend"])
                                # This is not wrong, but should happen in the mappings as well
                                # For the time being, it should not be added to either.
                                '''

                                # Add to overall path length
                                correct_alignment_length += len(range(interval_start, interval_end))


                '''
                for node in my_gaf_int:
                        node_as_str = str(node)
                        if node_as_str in mappings_for_path.keys():
                                curr_int = mappings_for_path[node_as_str]
                                correct_alignment_length += len(range(curr_int["start"], curr_int["end"]))
                        
                        # Add length of nodes not in path (DO NOT ADD BACK)
                        else:
                                for path in mappings.keys():
                                        if node_as_str in mappings[path].keys():
                                                curr_int = mappings[path][node_as_str]
                                                intervals.append(len(range(curr_int["start"], curr_int["end"])))
                '''

                overlap_ratio = correct_alignment_length / my_gaf_path_len
                threshold = float(args["Threshold"])
                if overlap_ratio > threshold:
                        #print("Read correctly mapped")
                        reads_mapped_correctly += 1
                else:
                        #print("Read not mapped correctly")
                        print("Incorrect alignment nodes: {}".format(my_gaf_int))
                        print("Read is: {}".format(my_gaf_name))
                        print("GAF alignment length: {}".format(my_gaf_path_len))
                        print("Correct alignment length: {}\n".format(correct_alignment_length))

        # Compute final stats
        correct_reads_ratio = reads_mapped_correctly / n_reads if n_reads != 0 else 0
        print("\nReads mapped correctly: {}/{} ({:.2f})".format(reads_mapped_correctly, n_reads, correct_reads_ratio))

        reads_not_mapped_correctly = n_reads-reads_mapped_correctly
        wrong_reads_ratio = reads_not_mapped_correctly / n_reads if n_reads != 0 else 0
        print("Reads NOT mapped correctly: {}/{} ({:.2f})".format(reads_not_mapped_correctly, n_reads, wrong_reads_ratio))
else:
       print("Path not found!") 


                        
