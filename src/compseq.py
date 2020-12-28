#!/bin/env python3
# -*- coding: utf-8 -*-
# This program is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# A copy of the GNU General Public License is available at http://www.gnu.org/licenses/gpl-3.0.html

"""Analyse the conservation of reversed complementary motifs in an alignment file."""

__author__ = "Julien Pichon"
__contact__ = "julien.pichon@cri-paris.org"
__version__ = "1.0.0"
__license__ = "GPL"

import argparse
from intervaltree import Interval, IntervalTree
from tqdm import tqdm


def get_args(argv = None):
	"""Retrieves the arguments of the program. Returns: An object that contains the arguments."""	
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--filename", help="Take nexus alignement file as input.")
	parser.add_argument("-o", "--output", help="Name of the tsv output file containing the\
	individual's name, the kmer sequence, the start position of the kmer sequence, the start\
	position of its reverse complementary sequence, number of individuals possessing this kmer\
	sequence and its reversed complementary at these positions, number of individuals possessing a\
	kmer sequence and its reverted complementary at these positions, kmer length.",
	default = "matching_kmer.tsv")
	parser.add_argument("-d", "--minimal_distance", type=int, help="Minimal distance between a kmer\
	and its reverse complementary sequence.", default = 100)
	parser.add_argument("-k", "--kmer_length", type=int, help="Complementary sequences length.",
	default = 10)
	parser.add_argument("-f", "--filter_kmer", help="Filter out complementary sequences belonging\
	to a greater kmer length. (e. g: AGGGA TCCCT. AGGG or GGGA will not count because\
	AGGGA have a reversed complementary sequence).", action="store_false")
	parser.add_argument('-s', "--select_sequence", nargs='*', help="Only the specified sequences\
	will count for the conservation.")
	parser.add_argument('-r', "--reject_sequence", nargs='*', help="The specified sequences will\
	not count for the conservation.")
	return parser.parse_args(argv)


def read_align(input_file, select_sequences, reject_sequences):
	"""Read alignment file and returns a dictionnary of sequences.
	
	Parameters
	----------
	arg1 : string
		Name of the input file.
	arg2 : list
		List of selected individuals.
	arg3 : list
		List of rejected individuals.
		
	Returns
	-------
	dict
		Dictionnary with the indivual as key and its sequence as value.
	"""
	dict_seq = {}
	with open(input_file, "r") as filin:
		for n, line in enumerate(filin):
			line = line.strip()
			# Parse nexus file to keep only sequences.
			if n < 7 or line.startswith(";") or line.startswith("END") or len(line) == 0:
				continue
			else:
				line = line.split()
				# Test for optionnal selected or rejected sequences arguments.
				if select_sequences == None and reject_sequences == None:
					dict_seq[line[0]] = line[1]
				elif select_sequences != None:
					if line[0] in select_sequences:
						dict_seq[line[0]] = line[1]
				elif reject_sequences != None:
					if line[0] not in reject_sequences:
						dict_seq[line[0]] = line[1]
	return dict_seq
				

def create_kmer_dict(dict_seq, k):
	"""Create a dictionnary of kmers with positions as values in the dictionnary of sequences.
	
	Parameters
	----------
	arg1 : dict
		Dictionnary with the indivual as key and its sequence as value.
	arg2 : int
		Length of the kmers.
		
	Returns
	-------
	dict
		Dictionnary of dictionnary of list as {individual1:{kmer1:[[start1:end1], [start2:end2]]}}
	"""
	dict_seq_kmer = {}
	for seq in dict_seq:
		dict_kmer = {}
		for position in range(len(dict_seq[seq]) - k):
			# Don't keep kmers containing gaps.
			if dict_seq[seq][position:position + k].count("-") != 0:
				continue
			# Test if the kmer is already present in the dict.
			else:
				if dict_seq[seq][position:position + k] not in dict_kmer:
					dict_kmer[dict_seq[seq][position:position + k]] = [position]
				else:
					dict_kmer[dict_seq[seq][position:position + k]].append(position)
		dict_seq_kmer[seq] = dict_kmer
	return dict_seq_kmer


def find_match(dict_seq_kmer, min_distance, k, file_output):
	"""Create a dictionnary of intervaltree containing the kmers and their reversed complementary.
	
	Parameters
	----------
	arg1 : dict
		Dictionnary of dictionnary of list as {individual1:{kmer1:[[start1:end1], [start2:end2]]}}
	arg2 : int
		Minimal distance between a kmer and its reverse complementary sequence.
	arg3 : string
		Name of the output file.
		
	Returns
	-------
	dict
		Dictionnary with individuals as key and the intervaltree of its kmers as value.
		{Individual1:([kmer_start:complementary_start = kmer], ...), ...}
	"""
	dict_interval = {}
	with open(file_output, "w") as filout:
		for seq in dict_seq_kmer:
			match = IntervalTree()
			for kmer in dict_seq_kmer[seq]:
				# Convert the kmer to its reversed complementary sequence.
				reverse = kmer[::-1]
				comp_reverse = reverse.replace("T", "z")
				comp_reverse = comp_reverse.replace("A", "T")
				comp_reverse = comp_reverse.replace("z", "A")
				comp_reverse = comp_reverse.replace("G", "z")
				comp_reverse = comp_reverse.replace("C", "G")
				comp_reverse = comp_reverse.replace("z", "C")
				if comp_reverse in dict_seq_kmer[seq]:
					# If only one position exists for this kmer and its reversed complementary.
					if (len(dict_seq_kmer[seq][kmer]) == 1 and 
					len(dict_seq_kmer[seq][comp_reverse]) == 1):
						if (int(dict_seq_kmer[seq][kmer][0]) + k <= int(dict_seq_kmer[seq]
						[comp_reverse][0]) and int(dict_seq_kmer[seq][comp_reverse][0]) - 
						int(dict_seq_kmer[seq][kmer][0]) + k <= min_distance):
							match[dict_seq_kmer[seq][kmer][0]:dict_seq_kmer[seq][comp_reverse][0]] = kmer
					# If several positions exist for this kmer and its reversed complementary.
					else:
						for position in dict_seq_kmer[seq][kmer]:
							for position_reverse in dict_seq_kmer[seq][comp_reverse]:
								if (int(position + k) <= int(position_reverse) and
								int(position_reverse) - int(position) + k <= min_distance):
									match[position:position_reverse] = kmer
			dict_interval[seq] = match
	return dict_interval


def k_filter(dict_interval, keep_all_kmer):
	filtered_tree = IntervalTree()
	count_dict = {}
	for seq in tqdm(dict_interval):
		for matching_kmer in sorted(dict_interval[seq]):
			if keep_all_kmer is False:
				if (Interval(matching_kmer.begin + 1, matching_kmer.end - 1, 
				matching_kmer.data[1:] + "A") not in dict_interval[seq] and
				Interval(matching_kmer.begin + 1, matching_kmer.end - 1, 
				matching_kmer.data[1:] + "C") not in dict_interval[seq] and
				Interval(matching_kmer.begin + 1, matching_kmer.end - 1, 
				matching_kmer.data[1:] + "G") not in dict_interval[seq] and
				Interval(matching_kmer.begin + 1, matching_kmer.end - 1, 
				matching_kmer.data[1:] + "T") not in dict_interval[seq] and
				Interval(matching_kmer.begin - 1, matching_kmer.end + 1, 
				"A" + matching_kmer.data[:-1]) not in dict_interval[seq] and
				Interval(matching_kmer.begin - 1, matching_kmer.end + 1, 
				"C" + matching_kmer.data[:-1]) not in dict_interval[seq] and
				Interval(matching_kmer.begin - 1, matching_kmer.end + 1, 
				"G" +  matching_kmer.data[:-1]) not in dict_interval[seq] and
				Interval(matching_kmer.begin - 1, matching_kmer.end + 1, 
				"T" + matching_kmer.data[:-1]) not in dict_interval[seq]):
					filtered_tree[matching_kmer.begin + 1:matching_kmer.end + 1] = [seq, matching_kmer.data]
			else:
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
	dict_seq = read_align(args.filename, args.select_sequence, args.reject_sequence)
	dict_seq_kmer = create_kmer_dict(dict_seq, args.kmer_length)
	dict_interval = find_match(dict_seq_kmer, args.minimal_distance, args.kmer_length, args.output)
	filtered_tree, count_dict = k_filter(dict_interval, args.filter_kmer)
	create_file(filtered_tree, count_dict, args.output, args.kmer_length)

