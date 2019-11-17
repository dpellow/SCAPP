# The entire Recycler2 workflow
# Allow flexibility for different parts to be run manually

import argparse
import glob
import json
import logging
import os, subprocess, sys
import time

import recyclelib.utils as utils

import make_fasta_from_fastg # for creating the fasta for read mappings
import find_plasmid_gene_matches # for running BLAST to find matches to the plasmid-specific genes
import parse_plasmid_scores # for parsing and transforming the plasmid classifier scores
import classify_fastg # for running PlasClass to classify the nodes in the graph
import recycle # for running the core algorithms
import create_hits_fasta # for creating fasta of filtered plasmids

def parse_user_input():
    parser = argparse.ArgumentParser(
        description='Recycler2 extracts likely plasmids from de novo assembly graphs'
    )
    parser.add_argument('-g','--graph',
     help='Assembly graph FASTG file to process',
     required=True, type=str
    )
    parser.add_argument('-o','--output_dir',
        help='Output directory',
        required=True, type=str
    )
    parser.add_argument('-k','--max_k',
        help='Integer reflecting maximum k value used by the assembler',
        required=False, type=int, default=55
    )
    parser.add_argument('-l', '--length',
     help='Minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
    )
    parser.add_argument('-m', '--max_CV',
     help='Coefficient of variation used for pre-selection [default: 0.5, higher--> less restrictive]',
      required=False, default=1./2, type=float
    )

    parser.add_argument('-p','--num_processes',
        help='Number of processes to use',
        required=False, type=int, default=1
    )

    parser.add_argument('-sc','--use_scores',
        help='Boolean flag of whether to use sequence classification scores in plasmid assembly',
        required=False, type=bool, default=True
    )

    parser.add_argument('-gh','--use_genes',
        help='Boolean flag of whether to use plasmid-specific gene hits in plasmid assembly',
        required=False, type=bool, default=True
    )
    parser.add_argument('-b','--bam',
        help='BAM file resulting from aligning reads to contigs file, filtering for best matches',
        required=False, type=str
    )
    parser.add_argument('-r1','--reads1',
        help='1st paired-end read file path',
        required=False, type=str
    )
    parser.add_argument('-r2','--reads2',
        help='1st paired-end read file path',
        required=False, type=str
    )
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-pc', '--plasclass', type=str, help="PlasClass score file with scores of the assembly graph nodes" )
    group.add_argument('-pf','--plasflow',type=str, help="PlasFlow score file with scores of the assembly graph nodes")

    return parser.parse_args()

def main():
    # Get command line args
    args = parse_user_input()

    fastg = args.graph
    outdir = args.output_dir
    max_k = args.max_k

    max_CV = args.max_CV
    min_length = args.length

    num_procs = args.num_processes

    # for optional workflow steps
    bamfile = args.bam
    reads1 = args.reads1
    reads2 = args.reads2

    plasclass_file = args.plasclass
    plasflow_file = args.plasflow # these are mutually exclusive

    # flags
    use_scores = args.use_scores
    use_genes = args.use_genes

    parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(parent_path,'data')
    bin_path  = os.path.join(parent_path, 'bin')


    # Get config variables
    try:
        config_path = os.path.join(bin_path,'config.json')
        with open(config_path) as config:
            config = json.load(config)
            bwa_path = config['BWA_PATH']
            ncbi_path = config["NCBI_PATH"]
            samtools_path = config["SAMTOOLS_PATH"]
    except:
        print("Error loading config variables. Please check config.json file")
        raise

    int_dir = os.path.join(outdir, 'intermediate_files')
    if not os.path.exists(int_dir):
        os.makedirs(int_dir)
    logs_dir = os.path.join(outdir, 'logs')
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)

    # Set up logging and write config and options to the log file
    logfile = os.path.join(logs_dir,"recycler2.log")
    logging.basicConfig(filemode='w', filename=logfile, level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%d/%m/%Y %H:%M')
    logger = logging.getLogger("recycle_logger")

    logger.info("Beginning Recycler2 workflow")
#TODO: ADD print of plasclass path
    logger.info("Got parameters:\n\tInput graph: {}\n\tOutput directory: {} \n\tMaximum k value: {}\
    \n\t# processes: {}\n\tMaximum CV: {}\n\tMinimum Length: {}\n\tBamfile: {}\n\tReads file 1: {}\
    \n\tReads file 2: {}\n\tUse scores: {}\n\tUse genes: {}\n\tPath to BWA executables: {}\
    \n\tPath to NCBI executables: {}".format(
        fastg, outdir, max_k, num_procs, max_CV, min_length, bamfile, reads1, reads2,
        use_scores, use_genes, bwa_path, ncbi_path))

    # TODO: Support GFA format: script to convert gfa to fastg

    # Step 1: Map reads and create BAM file
    if bamfile is None:
        # first create fasta of contigs in only one direction
        logger.info('Creating fasta from fastg')
        nodes_fasta = os.path.join(int_dir, 'assembly_graph.nodes.fasta')
        make_fasta_from_fastg.parse_lines(fastg, nodes_fasta)

        print("Creating BAM file of read mappings with BWA")
        time_start = time.time()
        logger.info("Creating BAM file")

        bwa_file = os.path.join(bwa_path,'bwa')
        samtools_file = os.path.join(samtools_path,'samtools')

        # BWA index
        bwa_outfile_name = os.path.join(logs_dir,'bwa_std.log')
        bwa_outfile = open(bwa_outfile_name,'w') # Don't clutter with print files
        cmd = bwa_file + " index " + nodes_fasta
        logger.info("Executing command: '{}'".format(cmd))
        subprocess.check_call(cmd, stderr=subprocess.STDOUT, stdout=bwa_outfile, shell=True)
        bwa_outfile.flush()

        # BWA mem and samtools view
        reads_bam = os.path.join(int_dir,"reads_pe.bam")
        cmd =  bwa_file + " mem -t " + str(num_procs) + " " + nodes_fasta + " " \
                + reads1 + " " + reads2 + " | " + \
                samtools_file + "  view -buS -@ " + str(num_procs-1) + " - > " + reads_bam
        logger.info("Executing command: '{}'".format(cmd))
        subprocess.check_call(cmd,stderr=subprocess.STDOUT, stdout=bwa_outfile, shell=True)
        bwa_outfile.flush()

        # Filter reads with samtools
        primary_bam = os.path.join(int_dir, "reads_pe_primary.bam")
        #TODO: WHAT SHOULD THE FLAG BE SET TO?
        cmd = samtools_file + " view -bF 0x0900 -@ " + str(num_procs-1) + " " + reads_bam + ' > ' + primary_bam
#        cmd = samtools_file + " view -bF 0x0800 -@ " + str(num_procs-1) + " " + reads_bam + ' > ' + primary_bam
        logger.info("Executing command: '{}'".format(cmd))
        subprocess.check_call(cmd,stderr=subprocess.STDOUT, stdout=bwa_outfile, shell=True)
        bwa_outfile.flush()

        # Sort filtered reads with samtools
        sorted_reads = os.path.join(int_dir, "reads_pe_primary.sort")
        cmd = samtools_file + " sort " + primary_bam + ' ' + sorted_reads
        logger.info("Executing command: '{}'".format(cmd))
        subprocess.check_call(cmd,stderr=subprocess.STDOUT, stdout=bwa_outfile, shell=True)
        bwa_outfile.flush()

        # Index filtered reads with samtools
        sorted_bam = os.path.join(int_dir, "reads_pe_primary.sort.bam")
        cmd = samtools_file + " index " + sorted_bam
        logger.info("Executing command: '{}'".format(cmd))
        subprocess.check_call(cmd,stderr=subprocess.STDOUT, stdout=bwa_outfile, shell=True)
        bwa_outfile.flush()

        bamfile = sorted_bam

        # Remove the intermediate files
        os.remove(primary_bam)
        os.remove(reads_bam)
        for f in glob.glob(nodes_fasta+"*"): os.remove(f)
        bwa_outfile.close()

        time_end = time.time()
        logger.info("{} seconds to create indexed sorted bam file".format(
            time_end-time_start))
    else:
        logger.info("Using file {} as the bamfile".format(bamfile))

    # Step 2: Classify contigs using PlasClass (default) and parse scores
    if use_scores and plasclass_file is None and plasflow_file is None:
        print("Getting scores of graph nodes")
        logger.info("Using PlasClass to obtain sequence scores")
        time_start = time.time()
        plasclass_file = os.path.join(int_dir, 'plasclass.out')
        sys_stdout = sys.stdout # don't want to clutter with more print statements
        sys_stderr = sys.stderr
        plasclass_outfile = os.path.join(logs_dir,'plasclass_std.log')
        stdfile = open(plasclass_outfile,'w')
        sys.stdout = stdfile
        sys.stderr = sys.stdout
#TODO: Add try-except blocks around main calls to scripts (AND IN THE SCRIPTS)
        classify_fastg.classify(fastg, plasclass_file, num_procs)
        sys.stdout = sys_stdout
        sys.stderr = sys_stderr
        stdfile.close()
        time_end = time.time()
        logger.info("{} seconds to classify the assembly graph".format(
            time_end-time_start))

    scores_file = None
    if use_scores:
        # Parse and transform the scores
        logger.info("Transforming scores")
        scores_file = os.path.join(int_dir, 'assembly_graph.nodes.scores')
        if plasflow_file:
            parse_plasmid_scores.transformPlasFlow(plasflow_file, scores_file)
        else:
            parse_plasmid_scores.transformPlasClass(plasclass_file, scores_file)

    # Step 3: BLAST for plasmid-specific genes and parse BLAST output
    gene_hits_path = None
    if use_genes:
        print('Finding plasmid-specific genes with BLAST')
        time_start = time.time()
        logger.info('Finding plasmid-specific genes with BLAST')
        # don't want to clutter with more print statements
        # a bit more complicated for BLAST which prints from c
        new_stdout = os.dup(sys.stdout.fileno())
        new_stderr = os.dup(sys.stderr.fileno())
        blast_outfile = os.path.join(logs_dir,'blast_std.log')
        stdfile = open(blast_outfile,'w')
        os.dup2(stdfile.fileno(),sys.stdout.fileno())
        os.dup2(stdfile.fileno(),sys.stderr.fileno())
        genefiles_path = os.path.join(data_path,'nt')
        protfiles_path = os.path.join(data_path,'aa')
        genefiles_list = [os.path.join(genefiles_path,fname) for fname in os.listdir(genefiles_path)]
        protfiles_list = [os.path.join(protfiles_path,fname) for fname in os.listdir(protfiles_path)]
        genefiles = ','.join(genefiles_list)
        protfiles = ','.join(protfiles_list)
        try:
            find_plasmid_gene_matches.find_plasmid_gene_matches(fastg, int_dir, genefiles, protfiles, None, \
                                    ncbi_path, num_procs)
        except:
            os.dup2(new_stdout, sys.stdout.fileno())
            os.dup2(new_stderr, sys.stderr.fileno())
            stdfile.close()
            print("Error finding plasmid genes - check BLAST output file (blast_std.log)")
            raise
        os.dup2(new_stdout, sys.stdout.fileno())
        os.dup2(new_stderr, sys.stderr.fileno())
        stdfile.close()
        gene_hits_path = os.path.join(int_dir,'hit_seqs.out')
        time_end = time.time()
        logger.info("{} seconds to find plasmid-specific gene hits".format(
            time_end-time_start))

    # Step 4: Run Recycler2 annotated-assembly-graph-based plasmid assembly
    print("Starting Recycler2 plasmid finding")
    logger.info("Starting plasmid finding")
    time_start = time.time()
    recycle.run_recycler2(fastg, outdir, bamfile, num_procs, max_k, \
                    gene_hits_path, use_genes, scores_file, use_scores, \
                    max_CV, min_length)
    time_end = time.time()
    logger.info("{} seconds to run Recycler2 plasmid finding".format(
        time_end-time_start))


    basename, _ = os.path.splitext(os.path.basename(fastg))
    fasta_ofile = os.path.join(outdir,basename+'.cycs.fasta')
    self_loops_ofile = os.path.join(outdir,basename+'.self_loops.fasta')

    # Step 5: Post-process filtering: BLAST output plasmids for plasmid-specific genes
    if use_genes:
        print("Filtering plasmids by plasmid-specific genes")
        logger.info("Filtering plasmids by plasmid-specific genes")
        time_start = time.time()


        # don't want to clutter with more print statements
        # a bit more complicated for BLAST which prints from c
        new_stdout = os.dup(sys.stdout.fileno())
        new_stderr = os.dup(sys.stderr.fileno())
        blast_outfile = os.path.join(logs_dir,'blast_std.log')
        stdfile = open(blast_outfile,'a')
        os.dup2(stdfile.fileno(),sys.stdout.fileno())
        os.dup2(stdfile.fileno(),sys.stderr.fileno())

        hit_plasmids_dir = os.path.join(int_dir,"hit_cycs")
        if not os.path.exists(hit_plasmids_dir):
            os.mkdir(hit_plasmids_dir)
        hit_plasmids_fname = os.path.join(hit_plasmids_dir,"hit_seqs.out")
        try:
            find_plasmid_gene_matches.find_plasmid_gene_matches(fasta_ofile, hit_plasmids_dir, genefiles, protfiles, \
                                    None, ncbi_path, num_procs)
        except:
            os.dup2(new_stdout, sys.stdout.fileno())
            os.dup2(new_stderr, sys.stderr.fileno())
            stdfile.close()
            print("Error filtering by plasmid genes. Check BLAST output file (blast_std.log)")
            raise

        os.dup2(new_stdout, sys.stdout.fileno())
        os.dup2(new_stderr, sys.stderr.fileno())
        stdfile.close()

        gene_filtered_ofile = os.path.join(outdir, basename+".gene_filtered_cycs.fasta")
        create_hits_fasta.create_hits(fasta_ofile, hit_plasmids_fname, gene_filtered_ofile)
        time_end = time.time()
        logger.info("{} seconds to filter plasmids using plasmid-specific gene hits".format(
            time_end-time_start))

    # Step 6: Post-process filtering: Classify gene filtered plasmids
    if use_scores:
        print("Getting scores of gene filtered plasmids")
        logger.info("Using PlasClass to obtain scores of cycles")
        time_start = time.time()
        plasclass_filtered_file = os.path.join(int_dir, 'plasclass_filtered.out')
        sys_stdout = sys.stdout # don't want to clutter with more print statements
        sys_stderr = sys.stderr
        plasclass_outfile = os.path.join(logs_dir,'plasclass_std.log')
        stdfile = open(plasclass_outfile,'a')
        sys.stdout = stdfile
        sys.stderr = sys.stdout
        classify_fastg.classify(fasta_ofile, plasclass_filtered_file, num_procs)
        logger.info("Transforming scores")
        plasmid_scores_file = os.path.join(int_dir, 'filtered_plasmids.scores')
        parse_plasmid_scores.transformPlasClass(plasclass_filtered_file, plasmid_scores_file)

        classified_plasmids_fname = os.path.join(hit_plasmids_dir,"classified_cycs.out")
        classification_thresh = 0.5
        with open(plasmid_scores_file) as f, open(classified_plasmids_fname,'w') as o:
            for line in f:
                splt = line.strip().split()
                if float(splt[1]) > classification_thresh: o.write(splt[0] + '\n')
        classification_filtered_ofile = os.path.join(outdir, basename+".classified_cycs.fasta")
        create_hits_fasta.create_hits(fasta_ofile, classified_plasmids_fname, classification_filtered_ofile)
        sys.stdout = sys_stdout
        sys.stderr = sys_stderr
        stdfile.close()
        time_end = time.time()
        logger.info("{} seconds to filter cycles by PlasClass score".format(
            time_end-time_start))
        os.remove(plasmid_scores_file)
        os.remove(plasclass_filtered_file)

    # Step 7: Create set of confident plasmid predictions:
    # Confident plasmid predictions - 2 out of 3 of: gene hits, plasmid classification, self loops
    if use_scores and use_genes:
        time_start = time.time()
        gene_hit_set = set()
        fp = open(gene_filtered_ofile, 'r')
        for name,_,_ in utils.readfq(fp):
            gene_hit_set.add(name)
        fp.close()

        classified_set = set()
        fp = open(classification_filtered_ofile, 'r')
        for name,_,_ in utils.readfq(fp):
            classified_set.add(name)
        fp.close()

        self_loop_set = set()
        fp = open(self_loops_ofile, 'r')
        for name,_,_ in utils.readfq(fp):
            self_loop_set.add(name)
        fp.close()

        classified_loops = classified_set & self_loop_set
        gene_hit_loops = gene_hit_set & self_loop_set
        classified_gene_hit = gene_hit_set & classified_set

        confident_plasmid_set = classified_loops | classified_gene_hit | gene_hit_loops

        confident_plasmids_fname = os.path.join(hit_plasmids_dir,"confident_cycs.out")
        with open(confident_plasmids_fname,'w') as o:
            for cyc in confident_plasmid_set:
                o.write(cyc+'\n')

        confident_plasmid_ofile = os.path.join(outdir, basename+".confident_cycs.fasta")
        create_hits_fasta.create_hits(fasta_ofile, confident_plasmids_fname, confident_plasmid_ofile)

        time_end = time.time()
        logger.info("{} seconds to filter confident plasmids".format(
            time_end-time_start))

    # Step 8: Clean up any remaining intermediate files
    if use_scores:
        os.remove(scores_file)

if __name__ == '__main__':
    main()
