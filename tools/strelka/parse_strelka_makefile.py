import argparse
import sys
import re

if __name__ == '__main__':

	desc = "Convert Strelka Makefile into a sub-Makefile used STRICTLY for parallelization by chromosome in Galaxy"
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument(
		"-m", "--makefile",
		required=True,
		help="Input Must be a Strelka Makefile"
		)
	parser.add_argument(
		"-c", "--chrom",
		type=str,
		required=False,
		help="Chromosome present in the Strelka Makefile"
		)
	parser.add_argument(
		"-cf", "--chrom_file",
		required=False,
		help="Chromosome File separated by newline"
		)
	parser.add_argument(
		'-o', '--output',
		required=True,
		help="Output Makefile"
		)
	args = parser.parse_args()

	if (not args.chrom and not args.chrom_file) or (args.chrom and args.chrom_file):
		sys.exit('Must define a single chromosome entry or a single chromosome file')

	chromosomes = [];

	if args.chrom_file:
		chromosomes = open(args.chrom_file).read().splitlines()

	if args.chrom:
		chromosomes = [args.chrom]
	add_at_end = []
	makefile_in = open(args.makefile, 'r')
	makefile_out = open(args.output, 'w')

	for line in makefile_in:
		if line.startswith('script_dir'):
			makefile_out.write(line)
		elif line.startswith('call_script'):
			makefile_out.write(line)		
		elif line.startswith('filter_script'):
			makefile_out.write(line)
		elif line.startswith('finish_script'):
			makefile_out.write(line)
		elif line.startswith('script_dir'):
			makefile_out.write(line)
		elif line.startswith('config_file'):
			makefile_out.write(line)
		elif line.startswith('analysis_dir'):
			makefile_out.write(line)
		elif line.startswith('results_dir'):
			makefile_out.write(line)		
		elif line.startswith('get_chrom_dir'):
			makefile_out.write(line)
		elif line.startswith('get_chrom_task'):
			makefile_out.write(line)
		elif line.startswith('get_bin_task'):
			makefile_out.write(line)
			makefile_out.write('\n')
		elif line.startswith('all:'):
			makefile_out.write('all:\n')
		elif re.search("--chrom=", line):
			if re.search(''.join(["\s|".join(chromosomes),'\s']), line):
				if re.search('bin', line):
					makefile_out.write(line)
				else:
					add_at_end.append(line)

	makefile_in.close()
	for line in add_at_end:
		makefile_out.write(line)
	makefile_out.close()
