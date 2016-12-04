import argparse
import math
import subprocess

if __name__ == "__main__":
    desc = "Create BAM interval file for easily parallelism in Galaxy"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
    	"--mode", "-m", 
    	required=True,
    	choices=['by_rname', 'by_chunk']
    	)
    parser.add_argument(
    	"--input", "-i",
    	type=str,
    	required=True,
    	help="Input BAM file"
    	)
    parser.add_argument(
    	"--output", "-o", 
    	type=str,
    	required=True,
    	help="Output BED file"
    	)
    parser.add_argument(
      "--order_file", "-of",
      type=str,
      required=False,
      help="Generate a file with the original contig order"
      )
    parser.add_argument(
    	  "--chunk_size", "-s", 
    	  type=int,
    	  required=False,
    	  default=50000000,
    	  help="The number of bases to include in a chunk interval, required when 'by_chunks' specified"
    	  )
    parser.add_argument(
        "--chromosome",
        action="store_true",
        required=False,
        help="Add the Chromosome to the Interval File"
        )
    parser.add_argument(
        "--start_position",
        action="store_true",
        required=False,
        help="Add the Start Position to the Interval File"
        )
    parser.add_argument(
        "--end_position",
        action="store_true",
        required=False,
        help="Add the End Position to the Interval File"
        )
    parser.add_argument(
        "--interval",
        action="store_true",
        required=False,
        help="Add the entire interval to the interval file"
        )
    parser.add_argument(
        "--prefixes_to_ignore",
        nargs="+",
        required=False,
        help="Ignore Intervals that contains the following prefix"
        )
    parser.add_argument(
        "--group_according_to_largest_chromosome",
        action="store_true",
        required=False,
        help="Output a set of files whose length do not sum larger then the largest chromosome"
        )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="outputs",
        help="Output directory for collection output"
        )
    args = parser.parse_args()

 
    order_file = None
    if not args.order_file == None:
        order_file = open(args.order_file, 'w')

    rnames = []
    data = {}
    commands = [
      "samtools", "view", 
      "-H", args.input
      ]
    p = subprocess.Popen(commands, stdout=subprocess.PIPE)
    out, err = p.communicate()
    seqs = []
    lengths = []
    for line in out.split("\n"):
        if not line.find("@SQ") == -1:
            tmp = line.split("\t")
            seqs.append(tmp[1][3:])
            lengths.append(tmp[2][3:])


    chunk_size = args.chunk_size

    for i in range(len(seqs)):
        rname = seqs[i]
        length = lengths[i]

        if not rname in data:
            rnames.append(rname)
            data[rname] = []

        if args.mode == 'by_rname':
            chunk_size = int(length)
        for i in range(int(math.ceil(float(length) / float(chunk_size)))):
            if (i+1) * chunk_size > length:
                data[rname].append([rname,(i * chunk_size) + 1, length])
            else:
                data[rname].append([rname,(i * chunk_size) + 1, (i+1) * chunk_size])

    if args.group_according_to_largest_chromosome :
        ids = []
        id_to_length = []
        bins_to_length = []
        bins_to_id = []
        for rname in rnames:
            ids.append(rname)

        for i in range(len(ids)):
            id_to_length.append(data[ids[i]][-1][2])
        
        max_length = 0
        for i in range(len(id_to_length)):
            if max_length < id_to_length[i]:
                max_length = id_to_length[i]
        
        for i in range(len(ids)):
            rname = ids[i]
            if any([not rname.find(x) == -1 for x in args.prefixes_to_ignore]):
                continue
            if not order_file == None:
                order_file.write(''.join([rname, '\n']))

            if bins_to_id == []:
                bins_to_id.append([rname])
                bins_to_length.append(id_to_length[i])
            else:
                added = False
                for bin_id in range(len(bins_to_length)):
                    if bins_to_length[bin_id] + id_to_length[i] < max_length:
                        bins_to_id[bin_id].append(rname)
                        bins_to_length[bin_id] += id_to_length[i]
                        added = True
                        break
                if added == False:
                    bins_to_id.append([rname])
                    bins_to_length.append(id_to_length[i])
        
        if not order_file == None:
            order_file.close()
         
        for i in range(len(bins_to_id)):
            output = open(''.join(["./outputs/samp", str(i+1),".bed"]), 'w')
            for rname in bins_to_id[i]:
                line = ''.join(['%s\n' % rname])
                output.write(line)
            output.close()
        
        
    else:
        output = open(args.output,'w')
        for i in range(len(seqs)):
            rname = seqs[i]
            for interval in data[rname]:
                if any([not interval[0].find(x) == -1 for x in args.prefixes_to_ignore]):
                    continue
                line = ''
                if args.interval:
                    line = ''.join([line, '%s:%s-%s\t' % (interval[0], interval[1], interval[2])])
                if args.chromosome:
                    line = ''.join([line,'%s\t' % interval[0]])
                if args.start_position:
                    line = ''.join([line,'%s\t' % interval[1]])
                if args.end_position:
                    line = ''.join([line,'%s\t' % interval[2]])
                output.write(''.join([line[:-1], '\n']))
