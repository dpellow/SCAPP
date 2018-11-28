#!/usr/bin/env python

import argparse, os
from recyclelib.utils import *
import pysam
import logging
import multiprocessing

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'Recycler extracts likely plasmids (and other circular DNA elements) from de novo assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) assembly graph FASTG file to process; recommended for spades 3.5: before_rr.fastg, for spades 3.6+:assembly_graph.fastg',
     required=True, type=str
     )
    parser.add_argument('-k','--max_k',
        help='integer reflecting maximum k value used by the assembler',
        required=True, type=int, default=55
        )
    parser.add_argument('-b','--bam',
        help='BAM file resulting from aligning reads to contigs file, filtering for best matches',
        required=True, type=str
        )
    parser.add_argument('-l', '--length',
     help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
     )
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection [default: 0.5, higher--> less restrictive]',
      required=False, default=1./2, type=float
      )
    parser.add_argument('-i','--iso',
        help='True or False value reflecting whether data sequenced was an isolated strain',
        required=False, type=bool, default=False
        )
    parser.add_argument('-o','--output_dir',
        help='Output directory',
        required=False, type=str
        )
    parser.add_argument('-p','--num_processes',
        help='Number of processes to use',
        required=False, type=int, default=1
        )
    parser.add_argument('-s','--scores',
            help='Contig plasmid scores file',
            required=False, type=str
        )
    return parser.parse_args()


if __name__ == '__main__':

    ####### entry point  ##############
    # inputs: spades assembly fastg, BAM of reads aligned to unitigs
    # outputs: fasta of cycles found by joining component edges

    ###################################
    # read in fastg, load graph, create output handle
    args = parse_user_input()
    num_procs = args.num_processes
    fastg = args.graph
    max_CV = args.max_CV
    max_k = args.max_k
    min_length = args.length
    fp = open(fastg, 'r')
    files_dir = os.path.dirname(fp.name)
    if args.scores:
        scores_file = args.scores
        use_scores = True
    else: use_scores = False
    # output 1 - fasta of sequences
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        basename = os.path.basename(fp.name)
        out_fp = os.path.join(args.output_dir, basename)
        (root, ext) = os.path.splitext(out_fp)
    else:
        (root,ext) = os.path.splitext(fp.name)

    fasta_ofile = root + ext.replace(".fastg", ".cycs.fasta")
    f_cycs_fasta = open(fasta_ofile, 'w')
    # output 2 - file containing path name (corr. to fasta),
    # path, coverage levels when path is added
    cycs_ofile = root + ext.replace(".fastg", ".cycs.paths_w_cov.txt")
    f_cyc_paths = open(cycs_ofile, 'w')
    bampath = args.bam
    ISO = args.iso

    logfile = root+ext.replace(".fastg", ".log")
    logging.basicConfig(filemode='w', filename=logfile, level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%d/%m/%Y %H:%M')
    logger = logging.getLogger("recycle_logger")

    # graph processing begins

    G = get_fastg_digraph(fastg)
####################
#    logger.info("Removing isolate nodes: %s" % ", ".join(list(nx.isolates(G))))######

    G.remove_nodes_from(list(nx.isolates(G)))
######### NOTE: if don't require circular path, should leave long isolates in

    cov_vals = [get_cov_from_spades_name(n) for n in G.nodes()]
    MED_COV = np.median(cov_vals)
    STD_COV = np.std(cov_vals)
    # set thresholds for max. CV, min
    # path coverage for allowing cross mappings
    if ISO:
        thresh = np.percentile(cov_vals, 95)
    else:
        thresh = np.percentile(cov_vals, 75)

    logger.info("Coverage threshold:%4f" % thresh)
    print(MED_COV, STD_COV, thresh)
    path_count = 0
    SEQS = get_fastg_seqs_dict(fastg,G)

#####################
    # add a score to every node, remove long nodes that are most probably chrom.
    if use_scores:
        get_node_scores(scores_file,G)
        remove_hi_confidence_chromosome(G)
        #TODO: take into account how to discount coverage of these nodes from their neighbours

#####################

    # gets set of long simple loops, removes short
    # simple loops from graph
    long_self_loops = get_long_self_loops(G, min_length, SEQS)

    final_paths_dict = {}

    for nd in long_self_loops:
        name = get_spades_type_name(path_count, nd,
            SEQS, max_k, G, get_cov_from_spades_name(nd[0]))
        final_paths_dict[name] = nd
        path_count += 1

#    comps = (G.subgraph(c) for c in nx.strongly_connected_components(G))
#   #below function is deprecated in nx 2.1....
    comps = nx.strongly_connected_component_subgraphs(G)
    COMP = nx.DiGraph()

    ###################################
    # iterate through SCCs looking for cycles
###########################

    # set up the multiprocessing #########################################
    job_queue = multiprocessing.JoinableQueue()
    result_queue = multiprocessing.Queue()
    workers = [multiprocessing.Process(target=process_component, args=(job_queue,result_queue, G, max_k, min_length, max_CV, SEQS, thresh, bampath, use_scores)) for i in xrange(num_procs)]
    for w in workers:
        w.daemon = True
        w.start()
    njobs = 0
    print("================== path, coverage levels when added ====================")

#    for c in comps:
    VISITED_NODES = set([]) # used to avoid problems due to RC components
    redundant = False
    for c in sorted(comps,key=len,reverse=True): # descending order - schedule large components first ######

	# check if any nodes in comp in visited nodes
        # if so continue
        for node in c.nodes():
             if node in VISITED_NODES: # I assume the below is an error?
        #     if c in VISITED_NODES:
                 redundant = True
                 break
        if redundant:
             redundant = False
             continue # have seen the RC version of component
        COMP = c.copy()

        job_queue.put(COMP)
        rc_nodes = [rc_node(n) for n in COMP.nodes()]
        VISITED_NODES.update(COMP.nodes())
        VISITED_NODES.update(rc_nodes)
        njobs+=1

    for i in xrange(num_procs):
        job_queue.put(None) # signal that the jobs are all done - one for each process
    job_queue.join() # wait till all processes return
####################################

    # done peeling
    # print final paths to screen
    print "%d jobs" % njobs
    print("==================final_paths identities after updates: ================")

    # write out sequences to fasta - long self-loops
    for p in final_paths_dict.keys():
        seq = get_seq_from_path(final_paths_dict[p], SEQS, max_k_val=max_k)
        print(final_paths_dict[p])
        print(" ")
        if len(seq)>=min_length:
            f_cycs_fasta.write(">" + p + "\n" + seq + "\n")
##################
    ## All processes done
    ## Read results queues of each
    for i in xrange(njobs):
        paths_set = result_queue.get()
        for p in paths_set:
            name = get_spades_type_name(path_count, p, SEQS, max_k, G)
            covs = get_path_covs(p,G)
            seq = get_seq_from_path(p, SEQS, max_k_val=max_k)
            print(p)
            if len(seq)>=min_length:
                f_cycs_fasta.write(">" + name + "\n" + seq + "\n")
                f_cyc_paths.write(name + "\n" +str(p)+ "\n" + str(covs)
                    + "\n" + str([get_num_from_spades_name(n) for n in p]) + "\n")
            path_count += 1
