#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import errno
import sys
import shutil
import subprocess
import datetime
import csv
import argparse
from multiprocessing import Pool, TimeoutError
import subprocess

# Need BIOPYTHON SEQ/IO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

# DEBUG
"""
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
"""

# Need primer3 and primer3-py
import primer3

class TFH :
	def __init__(self, ref, primers, outdir) :
		self.outdir = outdir
		self.ref = ref
		self.primers = primers

		if os.path.isfile(ref) :
			self.ref = os.path.abspath(ref)
		else :
			raise Exception("ERROR: Reference genome .fa file does not exist!")

		if os.path.isfile(primers) :
			self.primers = os.path.abspath(primers)
		else :
			raise Exception("ERROR: Primers fasta file does not exist!")

	def make_outdir(self) :
		try :
			os.makedirs(self.outdir)
			self.outdir = os.path.abspath(self.outdir)
		except OSError as e :
			if e.errno != errno.EEXIST :
				raise Exception("ERROR: Cannot create output directory!")
			else :
				self.outdir = os.path.abspath(self.outdir)
				print("WARNING: Directory already exists!")
	def __str__(self) :
		return "Reference: {}\nPrimers: {}\nOut directory: {}".format(self.ref, self.primers, self.outdir)

class FH :
	def __init__(self, ref, bed, outdir) :
		self.outdir = outdir
		self.ref = None
		self.bed = None
		if os.path.isfile(ref) :
			self.ref = os.path.abspath(ref)
		else :
			raise Exception("ERROR: Query .fa file does not exist!")

		if os.path.isfile(bed) :
			self.bed = os.path.abspath(bed)
		else :
			raise Exception("ERROR: Regions .bed file does not exist!")

		self.make_outdir()

	def make_outdir(self) :
		try :
			os.makedirs(self.outdir)
			self.outdir = os.path.abspath(self.outdir)
		except OSError as e :
			if e.errno != errno.EEXIST :
				raise Exception("ERROR: Cannot create output directory!")
			else :
				self.outdir = os.path.abspath(self.outdir)
				print("WARNING: Directory already exists!")

	def __str__(self) :
		return "Reference: {}\nRegions: {}\nOut directory: {}".format(self.ref, self.bed, self.outdir)

class Region :
	def __init__(self, name, start, end, ctg, sequence) :
		self.name = name
		self.start = start
		self.end = end
		self.ctg = ctg
		self.seq = sequence

	def __str__(self) :
		return "Region {} ({} from {} to {}; {}bp)".format(self.name, self.ctg, self.start, self.end, len(self.seq))

def to_chunks(bed, nproc) :
	reglist = []
	chunks = []
	for line in open(bed, "r") :
		reglist.append(line.strip())
	for i in range(0, len(reglist), nproc) :
		chunks.append(reglist[i:i+nproc])
	return chunks

def get_regions(job) :

	ref = job[0]
	chunk = job[1]
	offset = job[2]

	regions = []

	for line in chunk :
		s = line.split("\t")
		ctg = s[0]
		start = int(s[1]) - offset
		end = int(s[2]) + offset
		if len(s) == 3 :
			name = ctg+":"+str(start)+"-"+str(end)
		elif len(s) == 4 :
			name = s[3]
		else :
			name = ctg+":"+str(start)+"-"+str(end)

		for record in SeqIO.parse(ref, "fasta") :
			if record.id == ctg :
				regions.append( Region(name, start, end, ctg, record.seq[start:end]) )

	return regions

def run_get_regions(chunks, nproc, ref, offset) :
	jobs = []
	for chunk in chunks :
		jobs.append([ref, chunk, offset])

	p = Pool(processes=nproc)
	regions = p.map(get_regions, jobs)
	return [r for sublist in regions for r in sublist]

def thermodynamics(job) :
	ref = job[0]
	chunk = job[1]
	tm_offset = job[2]
	tm_size = job[3]

	record_dict = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))

	t_results = []

	for line in chunk :
		if line[0] == "#" :
			continue
		s = line.split("\t")
		if  int(s[6]) < tm_size :
			t_results.append(line + "\t{}\t{}\t{}\t{}\t{}".format("/", "/", "/", "/", "/"))
			continue

		tStart = int(s[2]) - 1 # 0-indexed position
		tEnd = int(s[3])

		seq1 = s[8]
		if tStart > tEnd :
			tStart, tEnd = tEnd, tStart
			tStart -= tm_offset
			tEnd += tm_offset
			seq2 = record_dict[s[1]].seq[tStart:tEnd].reverse_complement()
		else :
			tStart -= tm_offset
			tEnd += tm_offset
			seq2 = record_dict[s[1]].seq[tStart:tEnd]

		# DEBUG
		#alignments = pairwise2.align.globalxx(Seq(seq1), seq2)
		#print(format_alignment(*alignments[0]))

		tm = primer3.bindings.calcTm(seq1)
		tR1 = primer3.bindings.calcHeterodimer(seq1, str(seq2))
		tR2 = primer3.bindings.calcEndStability(seq1, str(seq2))

		t_results.append(line + "\t{}\t{}\t{}\t{}\t{}".format(tm, tR1.tm, tR1.dg, tR2.tm, tR2.dg))

	return t_results

def run_thermodynamics(chunks, nproc, ref, tm_offset, tm_size) :
	jobs = []
	for chunk in chunks :
		jobs.append([ref, chunk, tm_offset, tm_size])

	p = Pool(processes=nproc)
	tm_results = p.map(thermodynamics, jobs)
	return [r for sublist in tm_results for r in sublist]

def get_primers(region, product_size_range, mintm, maxtm, mingc, maxgc) :
	seq = region.seq

	primer_dict = primer3.bindings.designPrimers(
	{
	'SEQUENCE_ID':region.name,
	'SEQUENCE_TEMPLATE':str(seq),
	},
	{
	'PRIMER_OPT_SIZE': 20,
	'PRIMER_PICK_INTERNAL_OLIGO': 1,
	'PRIMER_INTERNAL_MAX_SELF_END': 8,
	'PRIMER_MIN_SIZE': 18,
	'PRIMER_MAX_SIZE': 25,
	'PRIMER_OPT_TM': int((mintm+maxtm)/2),
	'PRIMER_MIN_TM': mintm,
	'PRIMER_MAX_TM': maxtm,
	'PRIMER_MIN_GC': mingc,
	'PRIMER_MAX_GC': maxgc,
	'PRIMER_MAX_POLY_X': 100,
	'PRIMER_INTERNAL_MAX_POLY_X': 5, # The maximum allowable length of a mononucleotide repeat, for example AAAAAA.
	'PRIMER_SALT_MONOVALENT': 50.0,
	'PRIMER_DNA_CONC': 50.0,
	'PRIMER_MAX_NS_ACCEPTED': 0,
	'PRIMER_MAX_SELF_ANY': 8,
	'PRIMER_MAX_SELF_END': 3,
	'PRIMER_MAX_END_GC': 2, # The maximum number of Gs or Cs allowed in the last five 3' bases of a left or right primer
	'PRIMER_PAIR_MAX_COMPL_ANY': 8,
	'PRIMER_PAIR_MAX_COMPL_END': 3,
	'PRIMER_PRODUCT_SIZE_RANGE': [product_size_range]
	})
	primer_dict["CHROM"] = region.ctg
	primer_dict["START"] = region.start
	primer_dict["END"] = region.end
	primer_dict["REGION_ID"] = region.name
	return primer_dict

def return_primers(regions) :
	p = Pool(processes=nproc)
	primers = p.map(get_primers, regions)
	return [r for sublist in primers for r in sublist]


def parse_designed_primers(filename) :
	if os.path.isfile(filename) :
		filename = os.path.abspath(filename)
	else :
		raise Exception("ERROR: {} does not exist!".format(filename))

	path, ext = os.path.splitext(filename)
	clear_output = os.path.join(path + ".clear.tsv")
	fasta_output = os.path.join(path + ".fasta")

	HDcolnum = {}
	parsed = []
	num_pairs = 0
	for n, line in enumerate(open(filename, "r")) :
		if n == 0 :
			headers = line.strip().split("\t")
			for n, col in enumerate(headers) :
				HDcolnum[col] = n
				if "PRIMER_RIGHT" in col and "GC_PERCENT" in col :
					num_pairs += 1
		else :
			s = line.strip().split("\t")
			for n in range(0, num_pairs) :
				if s[HDcolnum["PRIMER_LEFT_"+str(n)+"_SEQUENCE"]] != "" :
					p = {} # "ID":[], "CHROM":[], "START":[], "END":[], "LEFT":[], "RIGHT":[], "LEFT_GC":[], "RIGHT_GC":[], "LEFT_TM":[], "RIGHT_TM":[]
					p["ID"] = s[HDcolnum["REGION_ID"]] + "_" + str(n)
					p["CHROM"] = s[HDcolnum["CHROM"]]
					p["START"] = get_start(s[HDcolnum["START"]], s[HDcolnum["PRIMER_LEFT_"+str(n)]])
					p["END"] = get_end(s[HDcolnum["START"]], s[HDcolnum["PRIMER_RIGHT_"+str(n)]])
					p["PRODUCT_SIZE"] = s[HDcolnum["PRIMER_PAIR_"+str(n)+"_PRODUCT_SIZE"]]
					p["LEFT"] = s[HDcolnum["PRIMER_LEFT_"+str(n)+"_SEQUENCE"]]
					p["RIGHT"] = s[HDcolnum["PRIMER_RIGHT_"+str(n)+"_SEQUENCE"]]
					p["LEFT_GC"] = s[HDcolnum["PRIMER_RIGHT_"+str(n)+"_GC_PERCENT"]]
					p["RIGHT_GC"] = s[HDcolnum["PRIMER_RIGHT_"+str(n)+"_GC_PERCENT"]]
					p["LEFT_TM"] = s[HDcolnum["PRIMER_RIGHT_"+str(n)+"_TM"]]
					p["RIGHT_TM"] = s[HDcolnum["PRIMER_RIGHT_"+str(n)+"_TM"]]
					parsed.append(p)


	with open(clear_output, 'w') as output_file :
		dict_writer = csv.DictWriter(output_file, fieldnames=parsed[0].keys(), delimiter="\t")
		dict_writer.writeheader()
		dict_writer.writerows(parsed)

	sequences = []
	for cd in parsed :
		desc = "left_primer|START:{}-{}|GC:{}|TM:{}".format(cd["CHROM"], cd["START"], cd["LEFT_GC"], cd["LEFT_TM"])
		sequences.append(SeqRecord(Seq(cd["LEFT"], IUPAC.ambiguous_dna), id=cd["ID"] + "_LEFT", description=desc))
		desc = "right_primer|END:{}-{}|GC:{}|TM:{}".format(cd["CHROM"], cd["END"], cd["RIGHT_GC"], cd["RIGHT_TM"])
		sequences.append(SeqRecord(Seq(cd["RIGHT"], IUPAC.ambiguous_dna), id=cd["ID"] + "_RIGHT", description=desc))

	with open(fasta_output, "w") as output_handle :
		SeqIO.write(sequences, output_handle, "fasta")


def get_start(real_start, tuple_region_start_0_indexed_region_length) :
	return int(real_start) + int(tuple_region_start_0_indexed_region_length.split(",")[0][1:]) + 1

def get_end(real_start, tuple_region_end_0_indexed_region_length) :
	return int(real_start) + int(tuple_region_end_0_indexed_region_length.split(",")[0][1:]) + 1

def run(cmd) :
	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	while True :                                            # Waits and prints cout
            line = proc.stdout.readline()                       # Reads line from stdout
            if line.strip() == "" :                             # If line is empty
                pass
            else :                                              # Else prints the line
                print(line.decode().strip())
            if not line :
                break                                           # If there is no piping in anymore
            proc.wait()

def str_to_bool(v) :
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


def main() :
	parser = argparse.ArgumentParser(description='Find primers in a fasta assembly from a bed file and test for alignment in reference assembly.')
	subparsers = parser.add_subparsers()
	des = subparsers.add_parser('design')
	des.add_argument('Reference',nargs=1,type=str,help="<STRING> A fasta file containing the target sequences.")
	des.add_argument('Regions',nargs=1,type=str,help="<STRING> A fasta file containing the query sequences.")
	des.add_argument('Output',nargs=1,type=str,help="<STRING> An output directory path.")
	des.add_argument('-of','--offset',nargs=1,type=int,default=[50],required=False,help="<INT> Offset around the region of interest. Default: 50.")
	des.add_argument('-p','--processes',nargs=1,type=int,default=[4],required=False,help="<INT> Maximum threads to use. Default: 4.")
	des.add_argument('-s1','--min-size',nargs=1,type=int,default=[200],required=False,help="<INT> Minimum product size. Default: 200.")
	des.add_argument('-s2','--max-size',nargs=1,type=int,default=[400],required=False,help="<INT> Maximum product size. Default: 400.")
	des.add_argument('-t1','--min-tm',nargs=1,type=int,default=[57],required=False,help="<INT> Minimum TM of primer. Default: 57.")
	des.add_argument('-t2','--max-tm',nargs=1,type=int,default=[63],required=False,help="<INT> Maximum TM of primer. Default: 63.")
	des.add_argument('-g1','--min-gc',nargs=1,type=int,default=[20],required=False,help="<INT> Minimum %GC of primer. Default: 20.")
	des.add_argument('-g2','--max-gc',nargs=1,type=int,default=[80],required=False,help="<INT> Maximum %GC of primer. Default: 80.")
	#parser.add_argument('--keep-temp', '-kp',type=str_to_bool, nargs='?', const=True, default=False, help="Keep temporary files (single alignments and single fasta). Unset by default.")
	des.set_defaults(func=design_primers)

	test = subparsers.add_parser('test')
	test.add_argument('Primers',nargs=1,type=str,help="<STRING> A fasta file containing the primers.")
	test.add_argument('Reference',nargs=1,type=str,help="<STRING> A fasta file containing the reference genome.")
	test.add_argument('Output',nargs=1,type=str,help="<STRING> An output directory path.")
	test.add_argument('--skip-tm',type=str_to_bool, nargs='?', const=True, default=False, help="<BOOL> Skip thermodynamics analysis of BLAST results. False by default.")
	test.add_argument('-to','--tm-offset',nargs=1,type=int,default=[3], required=False,help="<INT> Numbers of neighbours to consider in TM check of BLAST result. Default: 3.")
	test.add_argument('-ts','--tm-size',nargs=1,type=int,default=[15], required=False,help="<INT> Minimum aligned length to perform thermodynamics analysis. Default: 15.")
	#test.add_argument('-ma','--min-align',nargs=1,type=int,default=[15], required=False,help="<INT> Minimum BLAST aligned length to report. Default: 15.")
	test.add_argument('-p','--processes',nargs=1,type=int,default=[4], required=False,help="<INT> Maximum threads to use. Default: 4.")
	test.set_defaults(func=test_primers)

	args = parser.parse_args()

	args.func(args)

	print("Done")
	sys.exit(0)

def test_primers(args) :
	ref = args.Reference[0]
	primers = args.Primers[0]
	out = args.Output[0]
	nproc = args.processes[0]
	tm_offset = args.tm_offset[0]
	tm_size = args.tm_size[0]
	#min_align = args.min_align[0]
	skip_tm = args.skip_tm

	# File Handler
	iTFH = TFH(ref, primers, out)

	# 2. Run blastmakedb
	db = os.path.join(iTFH.outdir, os.path.basename(iTFH.ref) + ".db")
	cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=iTFH.ref, out=db)
	print("Building BLAST Database...")
	print(cline)
	run(cline.__str__())

	# 3. Run short-blast
	result = os.path.join(iTFH.outdir, os.path.basename(iTFH.primers) + ".blast.tsv")
	result_tmp = os.path.join(iTFH.outdir, os.path.basename(iTFH.primers) + ".tmp")
	cline = NcbiblastnCommandline(query=iTFH.primers, db=db, task="blastn-short", num_threads=nproc, outfmt="6 qseqid sseqid sstart send mismatch qlen length pident qseq sseq", out=result_tmp)
	print("Running short-BLAST...")
	print(cline)
	run(cline.__str__())

	f = open(result, "w")
	f.write("#PrimerName\tTargetName\tTargetStart\tTargetEnd\t#Mismatches\tPrimerLength\tAlignedLength\t%Identity\tPrimerSeq\tContigSeq\n")
	f.writelines(open(result_tmp, "r").readlines())
	f.close()
	os.remove(result_tmp)

	if skip_tm :
		return

	# 4. Thermodynamics of BLAST results
	print("Running thermodynamic check on blast results...")
	tm_result_file = os.path.join(iTFH.outdir, os.path.basename(iTFH.primers) + ".blast.TM.tsv")
	chunks = to_chunks(result, nproc)
	tm_result = run_thermodynamics(chunks, nproc, iTFH.ref, tm_offset, tm_size)

	# PrimerName	TargetName	TargetStart	TargetEnd	#Mismatches	PrimerLength	AlignedLength	%Identity	PrimerSeq	ContigSeq	Struct_found	TM	DG	DH	DS
	f = open(tm_result_file, "w")
	f.write("#PrimerName\tTargetName\tTargetStart\tTargetEnd\t#Mismatches\tPrimerLength\tAlignedLength\t%Identity\tPrimerSeq\tContigSeq\tPrimerTM\tHeteroDimerTM\tHeteroDimerDG\t3EndStabilityTM\t3EndStabilityDG\n")
	for line in tm_result :
		f.write(line + "\n")
	f.close()

def design_primers(args) :
	ref = args.Reference[0]
	bed = args.Regions[0]
	out = args.Output[0]
	off = args.offset[0]
	nproc = args.processes[0]
	minsize = args.min_size[0]
	maxsize = args.max_size[0]
	mintm = args.min_tm[0]
	maxtm = args.max_tm[0]
	mingc = args.min_gc[0]
	maxgc = args.max_gc[0]

	# File Handler
	iFH = FH(ref, bed, out)

	# Get regions
	chunks = to_chunks(iFH.bed, nproc)
	regions = run_get_regions(chunks, nproc, iFH.ref, off)

	# Get primers for each region
	primers = []
	for n, region in enumerate(regions) :
		if region.end - region.start < minsize :
			raise Exception("ERROR: Region is too small for specified PCR product size!")

		primer_dict = get_primers(region, [minsize, maxsize], mintm, maxtm, mingc, maxgc)
		primers.append(primer_dict)

	keys = []
	for d in primers :
		for k in d.keys() :
			if k not in keys :
				keys.append(k)
			continue

	of = os.path.join(iFH.outdir, 'primers.tsv')
	with open(of, 'w') as output_file :
		dict_writer = csv.DictWriter(output_file, fieldnames=keys, delimiter="\t")
		dict_writer.writeheader()
		dict_writer.writerows(primers)

	parse_designed_primers(of)

if __name__ == '__main__':
	main()
