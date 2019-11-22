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

from Bio import SeqIO # Need BIOPYTHON SEQ/IO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import primer3 # Need prrimer3 and primer3-py

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

def read_bed(bed, nproc) :
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
	"""
	for r in regions :
		print(r)
	"""

	return regions

def get_regs(chunks, nproc, ref, offset) :
	jobs = []
	for chunk in chunks :
		jobs.append([ref, chunk, offset])

	p = Pool(processes=nproc)
	regions = p.map(get_regions, jobs)
	return [r for sublist in regions for r in sublist]

def get_primers(region) :
	seq = region.seq
	if len(seq) < 100 :
		PRIMER_PRODUCT_SIZE_RANGE = [75,100]
	elif len(seq) < 200 :
		PRIMER_PRODUCT_SIZE_RANGE = [150,200]
	elif len(seq) < 300 :
		PRIMER_PRODUCT_SIZE_RANGE = [200,300]
	elif len(seq) < 500 :
		PRIMER_PRODUCT_SIZE_RANGE = [300,500]
	elif len(seq) < 750 :
		PRIMER_PRODUCT_SIZE_RANGE = [500,750]
	elif len(seq) < 1000 :
		PRIMER_PRODUCT_SIZE_RANGE = [750,1000]
	elif len(seq) < 2000 :
		PRIMER_PRODUCT_SIZE_RANGE = [1000,2000]
	else :
		PRIMER_PRODUCT_SIZE_RANGE = [2000,5000]

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
	'PRIMER_OPT_TM': 60.0,
	'PRIMER_MIN_TM': 57.0,
	'PRIMER_MAX_TM': 63.0,
	'PRIMER_MIN_GC': 20.0,
	'PRIMER_MAX_GC': 80.0,
	'PRIMER_MAX_POLY_X': 100,
	'PRIMER_INTERNAL_MAX_POLY_X': 100,
	'PRIMER_SALT_MONOVALENT': 50.0,
	'PRIMER_DNA_CONC': 50.0,
	'PRIMER_MAX_NS_ACCEPTED': 0,
	'PRIMER_MAX_SELF_ANY': 12,
	'PRIMER_MAX_SELF_END': 8,
	'PRIMER_PAIR_MAX_COMPL_ANY': 12,
	'PRIMER_PAIR_MAX_COMPL_END': 8,
	'PRIMER_PRODUCT_SIZE_RANGE': [PRIMER_PRODUCT_SIZE_RANGE]
	})
	primer_dict["REGION_ID"] = region.name
	return primer_dict

def return_primers(regions) :
	p = Pool(processes=nproc)
	primers = p.map(get_primers, regions)
	return [r for sublist in primers for r in sublist]

def main() :
	parser = argparse.ArgumentParser(description='Find primers in a fasta assembly from a bed file')
	subparsers = parser.add_subparsers()
	des = subparsers.add_parser('design')
	des.add_argument('Reference',nargs=1,type=str,help="A fasta file containing the target sequences")
	des.add_argument('Regions',nargs=1,type=str,help="A fasta file containing the query sequences")
	des.add_argument('Output',nargs=1,type=str,help="An output directory path")
	des.add_argument('--offset','-o',nargs=1,type=int,default=[50],required=False,help="Offsets around the region of interest. Default: 50")
	des.add_argument('--processes','-p',nargs=1,type=int,default=[4],required=False,help="Maximum threads to use. Default: 4.")
	#parser.add_argument('--keep-temp', '-kp',type=str_to_bool, nargs='?', const=True, default=False, help="Keep temporary files (single alignments and single fasta). Unset by default.")
	des.set_defaults(func=design_primers)

	test = subparsers.add_parser('test')
	test.add_argument('Primers',nargs=1,type=str,help="A fasta file containing the primers")
	test.add_argument('Reference',nargs=1,type=str,help="A fasta file containing the reference genome")
	test.add_argument('Output',nargs=1,type=str,help="An output directory path")
	test.add_argument('--processes','-p',nargs=1,type=int,default=[4],required=False,help="Maximum threads to use. Default: 4.")
	test.set_defaults(func=test_primers)

	args = parser.parse_args()

	args.func(args)

	print("Done")
	sys.exit(0)
	#return parser.parse_args()

def test_primers(args) :
	ref = args.Reference[0]
	primers = args.Primers[0]
	out = args.Output[0]
	nproc = args.processes[0]

	# File Handler
	iTFH = TFH(ref, primers, out)
	"""
	# 1. Get primers sequences
	primers = [line.strip().split("\t") for line in open(iTFH.primers, "r")]
	"""
	# 2. Run blastmakedb
	db = os.path.join(iTFH.outdir, os.path.basename(iTFH.ref) + ".db")
	if not os.path.isfile(db) :
		cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=iTFH.ref, out=db)
		print("Running:")
		print(cline)
		run(cline.__str__())
	else :
		pass
	# 3. Run short-blast
	result = os.path.join(iTFH.outdir, os.path.basename(iTFH.primers) + ".blast")
	result_tmp = os.path.join(iTFH.outdir, os.path.basename(iTFH.primers) + ".tmp")
	cline = NcbiblastnCommandline(query=iTFH.primers, db=db, task="blastn-short", num_threads=nproc, outfmt="10 qseqid sseqid sstart send mismatch length pident qseq sseq", out=result_tmp)
	print("Running:")
	print(cline)
	run(cline.__str__())

	f = open(result, "w")
	f.write("PrimerName,TargetName,TargetStart,TargetEnd,#Mismatches,AlignedLength,%Identity,PrimerSeq,ContigSeq\n")
	f.writelines(open(result_tmp, "r").readlines())
	f.close()
	os.remove(result_tmp)

def design_primers(args) :
	ref = args.Reference[0]
	bed = args.Regions[0]
	out = args.Output[0]
	off = args.offset[0]
	nproc = args.processes[0]

	# File Handler
	iFH = FH(ref, bed, out)
	# Get regions
	chunks = read_bed(iFH.bed, nproc)
	regions = get_regs(chunks, nproc, iFH.ref, off)
	# Get primers for each region
	primers = []
	for n, region in enumerate(regions) :
		primer_dict = get_primers(region)
		primers.append(primer_dict)

	keys = []
	for d in primers :
		for k in d.keys() :
			if k not in keys :
				keys.append(k)
			continue

	with open(os.path.join(iFH.outdir, 'primers.csv'), 'w') as output_file :
		dict_writer = csv.DictWriter(output_file, fieldnames=keys, delimiter="\t")
		dict_writer.writeheader()
		dict_writer.writerows(primers)

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

if __name__ == '__main__':
	main()
