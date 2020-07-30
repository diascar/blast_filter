#!/usr/bin/env python3

import fileinput
import os
import sys


if __name__ == "__main__":

	identity = int(sys.argv[1])
	coverage = int(sys.argv[2])
	base_dir = os.getcwd()

# creating empty lists to hold lines from input files
	maf_data = []
	fas_data = []
	bla_data = []
	rep_data = []

# reading the input files

	with fileinput.input(files=sys.argv[3:]) as file:
		for line in file:
			if set(['.fas', '.fa', '.fasta']).intersection(set([os.path.splitext(file.filename())[1]])):
				fas_data.append(line)
			elif ".blast" in file.filename():
				bla_data.append(line)
			elif ".db" in file.filename():
				rep_data.append(line)

# building a python class to deal with blast data
	class SeqBlast:
		def __init__(self, name, repName, identity, alnSize, evalue, sequence):
			self.name = name
			self.repName = repName
			self.identity = float(identity)
			self.alnSize = int(alnSize)
			self.evalue = float(evalue)
			self.sequence = sequence

		def get_fasta(self):
			return(">{}\n{}\n".format(self.name, self.sequence))

		def get_tab_data(self):
			return("{}\t{}\t{}\t{}\t{}".format(self.name, self.repName, self.identity, self.alnSize, self.evalue))

# building a dictionary to hold data from fasta files
	def fas_to_dict(fas_list):
		f_dict = {}
		for li in fas_list:
			if li[0] == ">":
				temp_name = li.strip()
			elif li[0] != ">":
				if temp_name not in f_dict:
					f_dict[temp_name] = li
				else:
					f_dict[temp_name] += li
		return(f_dict)

	fasta_dictionary = fas_to_dict(fas_data)

	repdb_dictionary = fas_to_dict(rep_data)

# building a dictionary of repetition names (key = repetition name, value = a list containing the repetition type and origin)
	rep_type_dic = {}

	for key in repdb_dictionary.keys():
		rep_type_dic[key.split()[0].strip().lstrip(">")] = key.split()[1:3]
		print(key.split()[0].strip().lstrip(">"))

# creating a list to hold SeqBlast objects
	work_list = []

	for i in range(len(bla_data)):
		if bla_data[i][0] != "#":
			temp_var = bla_data[i].split()
			work_list.append(SeqBlast(temp_var[0], temp_var[1], temp_var[2], temp_var[3], temp_var[10], fasta_dictionary[">"+temp_var[0]] ))

# creating a dictionary with repetitive DNA as keys and a list of matching sequences as values
	rep_dict = {}

	for j in range(len(work_list)):
		if work_list[j].repName not in rep_dict:
			if work_list[j].identity > identity and work_list[j].alnSize > coverage:
				rep_dict[work_list[j].repName] = {work_list[j].name:work_list[j].sequence}
		else:
			if work_list[j].identity > identity and work_list[j].alnSize > coverage:
				if work_list[j].name not in rep_dict[work_list[j].repName]:
					rep_dict[work_list[j].repName][work_list[j].name] = work_list[j].sequence
				else:
					rep_dict[work_list[j].repName][work_list[j].name+"_"+str(j)] = work_list[j].sequence


	new_dic = {}
	for k in repdb_dictionary:
		ke = k.split()[0].strip().strip(">")
		new_dic[ke] = repdb_dic[k]


	for k in rep_dict:
		if len(rep_dict[k]) > 1:
			fname = os.path.join(base_dir, k+".fa")
			with open(fname, "w") as f:
				for v in rep_dict[k]:
					f.write(">{}\n{}\n".format(v, rep_dict[k][v]))
			ref_name = os.path.join(base_dir, k+"_ref.fa")
			with open(ref_name, "w") as fr:
				fr.write(">{}\n{}".format(k, new_dic[k]))

	table_name = os.path.join(base_dir, "RepDNA_table.txt")
	with open(table_name, "w") as table:
		table.write("query\trep DNA name\tidentity\talignment size\te-value\trep type\trep origin\n")
		for n in range(len(work_list)):
			table.write(work_list[n].get_tab_data()+"\t{}\t{}\n".format(rep_type_dic[work_list[n].repName][0], rep_type_dic[work_list[n].repName][1] ))

