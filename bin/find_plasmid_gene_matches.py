#!/usr/bin/env python

# get the plasmid gene matches in the assembled contigs, given lists of plasmid-specific genes and contigs
# usage: python -c <contigs file> -o <output dir> -g <plasmid-specific gene files> -p <plasmid-specific protein files> -dg <blast gene database> -dp <blast protein database>

import argparse, os, subprocess, re
from subprocess import CalledProcessError

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'find_plasmid_gene_matches finds the list of contigs with matches to the plasmid specific genes'
        )
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('-f','--fasta',
     help='input file - contigs in fasta format',
     type=str
     )
    g.add_argument('-db', '--blastdb',
      help='path of blast database for the contigs',
      type=str
     )
    parser.add_argument('-o','--output',
     help='path to directory to write output and temporary files',
     required=True, type=str
     )
    parser.add_argument('-g','--genefiles',
     help='csv list of paths to plasmid-specific gene fasta files',
     required=False, type=str
     )
    parser.add_argument('-p','--proteinfiles',
     help='csv list of paths to plasmid-specific protein fasta files',
     required=False, type=str
     )
    parser.add_argument('--ncbi_bin',
     help='path of ncbi executables executable',
     default="", type=str
     )
    parser.add_argument('-n', '--nthreads',
     help='number of threads for blastn and tblastn to use',
     default=1, type=int
     )
    parser.add_argument('-t', '--thresh',
     help='threshold for sequence identity and length of blast hits',
     default=0.75, type=float
     )
    parser.add_argument('-c', '--clean',
     help='Optionally remove all files that were created.',
     action='store_true'
     )

    return parser.parse_args()


def create_db(infile,makeblastdb,outdir):
    ''' create a blast database for the infile
    '''
    outputfile = os.path.basename(infile) + ".blastdb"
    outputpath = os.path.join(outdir,outputfile)

    command = makeblastdb + " -in " + infile + " -dbtype nucl -parse_seqids -out " + outputpath

    print "Running command: " + command
    try:
        subprocess.check_call(command, shell=True)
    except CalledProcessError:
        print "Error creating blast database. Check path to makeblastdb executable."
        raise

    return outputpath


def blast_gene_file(gf, dbpath, blastn, numthreads, outdir):
    ''' run blastn of the genefile against the blast database
    '''

    outputfile = os.path.basename(gf) + "_blastdb.out"
    outputpath = os.path.join(outdir,outputfile)

    command = blastn + " -task megablast -db " + dbpath + " -query " + gf + \
                " -out " + outputpath + " -num_threads " + str(numthreads) + \
                ' -outfmt "6 qseqid sseqid length pident qlen slen evalue"'

    print "Running command: " + command
    try:
        subprocess.check_call(command, shell=True)
    except CalledProcessError:
        print "Error running blastn. Check paths to blastn executable and input files."
        raise

    return outputpath


def blast_protein_file(pf, dbpath, tblastn, numthreads, outdir):
    ''' same, for protein lists
    '''
    outputfile = os.path.basename(pf) + "_blastdb.out"
    outputpath = os.path.join(outdir,outputfile)

    num_tblastn_threads = min(numthreads, 8) # tblastn is buggy and seg faults randomly
                                             # seems to do this less when there are fewer threads

    command = tblastn + " -db " + dbpath + " -db_gencode 11 -query " + pf + \
                " -out " + outputpath + " -num_threads " + str(num_tblastn_threads) + \
                ' -outfmt "6 qseqid sseqid length pident qlen slen evalue"'

    print "Running command: " + command
    try:
        subprocess.check_call(command, shell=True)
    except CalledProcessError:
        print "Error running tblastn. Check paths to tblastn executable and input files."
        raise

    return outputpath


def get_blast_hits(blastfile, hit_contigs_set, thresh=0.75):
    with open(blastfile) as f:
        for line in f:
            splt = line.strip().split('\t')
            contig = splt[1]
            contig = re.split(':|;',contig)[0]
            percentid = float(splt[3])
            genelen = int(splt[4])
            matchlen = int(splt[2])
            if percentid > thresh*100 and float(matchlen)/float(genelen) > thresh:
                hit_contigs_set.add(contig)


def main():
    args = parse_user_input()
    infile = args.fasta
    outdir = args.output
    genefiles = args.genefiles
    proteinfiles = args.proteinfiles
    dbpath = args.blastdb
    makeblastdb = os.path.join(args.ncbi_bin,'makeblastdb')
    blastn = os.path.join(args.ncbi_bin,'blastn')
    tblastn = os.path.join(args.ncbi_bin,'tblastn')
    numthreads = args.nthreads
    thresh = args.thresh
    clean = args.clean

    # create blast databases for the contigs
    if dbpath is None:
        print "No blast db provided, creating one"
        dbpath = create_db(infile, makeblastdb, outdir)

    # run BLAST searches for the genes in the gene lists
    genefiles_blast_results = []
    protfiles_blast_results = []

    if genefiles is not None:
        print "Running blast search for gene (nt) sequences in blast db"
        genefiles_list = genefiles.split(',')
        for gf in genefiles_list:
            gf_blastpath = blast_gene_file(gf, dbpath, blastn, numthreads, outdir)
            genefiles_blast_results.append(gf_blastpath)

    # run BLAST searches for the proteins in the protein lists
    if proteinfiles is not None:
        print "Running blast search for protein (aa) sequences in blast db"
        protfiles_list = proteinfiles.split(',')
        for pf in protfiles_list:
            pf_blastpath = blast_protein_file(pf, dbpath, tblastn, numthreads, outdir)
            protfiles_blast_results.append(pf_blastpath)

    # parse the output files to create set of contig hits
    hit_contigs_set = set()
    for gbr in genefiles_blast_results:
        get_blast_hits(gbr, hit_contigs_set, thresh)
    for pbr in protfiles_blast_results:
        get_blast_hits(pbr, hit_contigs_set, thresh)

    # write out the hit contigs to output file
    hit_contigs_list = list(hit_contigs_set)
    print "{} contigs hit".format(len(hit_contigs_list))
    outputfile = os.path.join(outdir, "hit_contigs.out")
    print "Writing list of hit contigs in: " + outputfile
    with open(outputfile,'w+') as o:
        o.write('\n'.join(hit_contigs_list))

    # remove intermediate files
    if clean:
        print "Removing intermediate files..."
        files_list = ' '.join(genefiles_blast_results) + ' ' + ' '.join(protfiles_blast_results)
        if infile: files_list += ' ' + dbpath + "*" # we created the blast db
        command = "rm " + files_list
        print "Running command: " + command
        subprocess.call(command, shell=True)


if __name__=='__main__':
    main()
