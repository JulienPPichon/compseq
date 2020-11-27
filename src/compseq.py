import argparse
from intervaltree import Interval, IntervalTree
from tqdm import tqdm


def get_args(argv = None):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--filename", help="Take fasta or alignment file as input")
	parser.add_argument("-o", "--output", help="Name of the tsv output file with kmer, its start position, and the start position of its reverse complementary sequence", default = "matching_kmer.tsv")
	parser.add_argument("-d", "--minimal_distance", type=int, help="Length of the minimal distance between a kmer and its reverse complementary sequence", default = 100)
	parser.add_argument("-k", "--kmer_length", type=int, help="Length of the matching sequences", default = 10)
	parser.add_argument("-f", "--keep_all_kmer", help="Filter out matching sequences belonging to a greater k length", action="store_false")
	parser.add_argument('-s', "--select_sequence", nargs='*', help="Comparison only between the specified sequences")
	return parser.parse_args(argv)


def read_fasta(input_file):
	seq = ""
	with open(input_file, "r") as filin:
		for line in filin:
			line = line.strip()
			if line.startswith(">"):
				continue
			else:
				seq = seq + line
	return seq


def read_align(input_file, select_sequences):
	dict_seq = {}
	with open(input_file, "r") as filin:
		for n, line in enumerate(filin):
			line = line.strip()
			if n < 7 or line.startswith(";") or line.startswith("END") or len(line) == 0:
				continue
			else:
				line = line.split()
				if select_sequences == None:
					dict_seq[line[0]] = line[1]
				else:
					if line[0] in select_sequences:
						dict_seq[line[0]] = line[1]
	return dict_seq
				

def create_kmer_dict(dict_seq, k):
	dict_seq_kmer = {}
	for seq in dict_seq:
		dict_kmer = {}
		for position in range(len(dict_seq[seq]) - k):
			if dict_seq[seq][position:position + k].count("-") != 0:
				continue
			else:
				if dict_seq[seq][position:position + k] not in dict_kmer:
					dict_kmer[dict_seq[seq][position:position + k]] = [position]
				else:
					dict_kmer[dict_seq[seq][position:position + k]].append(position)
		dict_seq_kmer[seq] = dict_kmer
	return dict_seq_kmer


def find_match(dict_seq_kmer, min_distance, k, file_output, keep_all_kmer):
	dict_interval = {}
	with open(file_output, "w") as filout:
		for seq in dict_seq_kmer:
			match = IntervalTree()
			for kmer in dict_seq_kmer[seq]:
				reverse = kmer[::-1]
				comp_reverse = reverse.replace("T", "z")
				comp_reverse = comp_reverse.replace("A", "T")
				comp_reverse = comp_reverse.replace("z", "A")
				comp_reverse = comp_reverse.replace("G", "z")
				comp_reverse = comp_reverse.replace("C", "G")
				comp_reverse = comp_reverse.replace("z", "C")
				if comp_reverse in dict_seq_kmer[seq]:
					if len(dict_seq_kmer[seq][kmer]) == 1 and len(dict_seq_kmer[seq][comp_reverse]) == 1:
						if int(dict_seq_kmer[seq][kmer][0]) + k <= int(dict_seq_kmer[seq][comp_reverse][0]) and int(dict_seq_kmer[seq][comp_reverse][0]) - int(dict_seq_kmer[seq][kmer][0]) + k <= min_distance:
							if keep_all_kmer is True:
								filout.write("{}\t{}\t{}\t{}\n".format(seq, kmer, dict_seq_kmer[seq][kmer][0] + 1, dict_seq_kmer[seq][comp_reverse][0] + 1))
							else:
								match[dict_seq_kmer[seq][kmer][0]:dict_seq_kmer[seq][comp_reverse][0]] = kmer
					else:
						for position in dict_seq_kmer[seq][kmer]:
							for position_reverse in dict_seq_kmer[seq][comp_reverse]:
								if int(position + k) <= int(position_reverse) and int(position_reverse) - int(position) + k <= min_distance:
									if keep_all_kmer is True:
										filout.write("{}\t{}\t{}\t{}\n".format(seq, kmer, position + 1, position_reverse + 1))
									else:
										match[position:position_reverse] = kmer
			dict_interval[seq] = match
	return dict_interval


def k_filter(dict_interval):
	filtered_tree = IntervalTree()
	count_dict = {}
	for seq in tqdm(dict_interval):
		for matching_kmer in sorted(dict_interval[seq]):
			if Interval(matching_kmer.begin + 1, matching_kmer.end - 1, matching_kmer.data[1:] + "A") not in dict_interval[seq] and Interval(matching_kmer.begin + 1, matching_kmer.end - 1, matching_kmer.data[1:] + "C") not in dict_interval[seq] and Interval(matching_kmer.begin + 1, matching_kmer.end - 1, matching_kmer.data[1:] + "G") not in dict_interval[seq] and Interval(matching_kmer.begin + 1, matching_kmer.end - 1, matching_kmer.data[1:] + "T") not in dict_interval[seq] and Interval(matching_kmer.begin - 1, matching_kmer.end + 1, "A" + matching_kmer.data[:-1]) not in dict_interval[seq] and Interval(matching_kmer.begin - 1, matching_kmer.end + 1, "C" + matching_kmer.data[:-1]) not in dict_interval[seq] and Interval(matching_kmer.begin - 1, matching_kmer.end + 1, "G" +  matching_kmer.data[:-1]) not in dict_interval[seq] and Interval(matching_kmer.begin - 1, matching_kmer.end + 1, "T" + matching_kmer.data[:-1]) not in dict_interval[seq]:
				filtered_tree[matching_kmer.begin + 1:matching_kmer.end + 1] = [seq, matching_kmer.data]			
				if str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1) not in count_dict:
					count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)] = {}
					count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)]["total"] = 1
					count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)][matching_kmer.data] = 1
				else:
					count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)]["total"] += 1
					if matching_kmer.data not in count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)]:
						count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)][matching_kmer.data] = 1
					else:
						count_dict[str(matching_kmer.begin + 1) + "_" + str(matching_kmer.end + 1)][matching_kmer.data] += 1
	return [filtered_tree, count_dict]
					
					
def create_file(filtered_tree, count_dict, file_output, k):
	with open(file_output, "w") as filout:
		for interval_obj in sorted(filtered_tree):
			filout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(interval_obj.data[0], interval_obj.data[1], interval_obj.begin, interval_obj.end, count_dict[str(interval_obj.begin) + "_" + str(interval_obj.end)]["total"], count_dict[str(interval_obj.begin) + "_" + str(interval_obj.end)][interval_obj.data[1]], k))
	

if __name__ == "__main__":

	argvals = None
	args = get_args(argvals)
	dict_seq = read_align(args.filename, args.select_sequence)
	dict_seq_kmer = create_kmer_dict(dict_seq, args.kmer_length)
	dict_interval = find_match(dict_seq_kmer, args.minimal_distance, args.kmer_length, args.output, args.keep_all_kmer)
	if args.keep_all_kmer is False:
		filtered_tree, count_dict = k_filter(dict_interval)
	create_file(filtered_tree, count_dict, args.output, args.kmer_length)

