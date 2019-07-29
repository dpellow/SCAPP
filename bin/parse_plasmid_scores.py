#!/usr/bin/env python

# parse contig scores from PlasFlow output files, write to file
# to be used by Recycler
# usage: python -i <input file> -o <output file>

import re, argparse, os, math

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'parse_plasmid_scores parses the output from PlasFlow and coverts to Recycler input format'
        )
    parser.add_argument('-i','--input',
     help='input file - generated by PlasFlow (v1.1)',
     required=True, type=str
     )
    parser.add_argument('-o','--output',
     help='output file name for contig plasmid scores',
     required=False, type=str
     )

    return parser.parse_args()

def parsePlasFlow(filename):
    """ Parse the file into lists of edge names and length-adjusted PlasFlow scores

    """

    with open(filename) as f:
        names = []
        lengths = []
        scores = []
        # # read header
        # plasmid_inds = []
        # is_header = True
        # for line in f:
        #     if is_header:
        #         is_header = False
        #         for i,v in enumerate(line.split('\t')):
        #             if v == "contig_name":
        #                 name_ind = i
        #             elif v == "contig_length":
        #                 len_ind = i
        #             elif "plasmid" in v:
        #                 plasmid_inds.append(i)
        #         continue
        #     # parse lines and store names, lengths, scores
        #     fields = line.split('\t')
        #     # convert name from fastg header format to single edge name
        #     name = re.split(':|;',fields[name_ind])[0]
        #     names.append(name)
        #     lengths.append(int(fields[len_ind]))
        #     score = sum([float(fields[x]) for x in plasmid_inds])
        name = ''
        for line in f:
            if line[0] == '>':
                edge_name = line.strip()[1:]
                # convert name from fastg header format to single edge name
                name = re.split(':|;',edge_name)[0]
                names.append(name)
                length = int(name.split('_')[3])
                lengths.append(length)
            else:
                line = line.strip()
                if len(line) > 0:
                    score = float(line)
                    scores.append(score)

    return names, lengths, scores


def transformByLength(lengths, scores):
    """ pull scores of short contigs towards 0.5
        For score x, length l:
        (x-0.5)*1/(1+e^(-0.001(l-2000))) + 0.5

    """
    transformed = []
    for i,v in enumerate(scores):
        l = lengths[i]
        weight = 1.0/(1+math.exp(-0.001*float(l-2000)))
        t = (v-0.5)*weight + 0.5
        transformed.append(t)
    return transformed

def writeOutput(outfile,edges,scores):
    with open(outfile,'w') as o:
        for e,s in zip(edges,scores):
            o.write(e+'\t'+str(s)+'\n')


if __name__=='__main__':
    args = parse_user_input()
    infile = args.input
    outfile = args.output
    edges, lengths, scores = parsePlasFlow(infile)
    scores = transformByLength(lengths, scores)
    writeOutput(outfile,edges,scores)
