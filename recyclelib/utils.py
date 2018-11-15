from __future__ import division
import numpy as np
#########################
import math ###############################
import networkx as nx
import re, pysam
import logging
import multiprocessing

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
######################
logger = logging.getLogger("recycle_logger")#######################

def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader:
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def get_node_scores(scores_file,G):
    """ Write the PlasFlow scores into each node in the graph
    """
    scores = {}
    with open(scores_file) as f:
        for line in f:
            split = line.strip().split()
            scores[split[0]] = float(split[1])
    for nd in G.nodes():
        G.add_node(nd, score=scores[nd])



def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

def get_num_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[1]
    return int(contig_length)

def get_length_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_spades_name(name):
    name_parts = name.split("_")
    cov = name_parts[5]
    if cov[-1]=="'": cov=cov[:-1]
    return float(cov)

def get_fastg_digraph(fastg_name):
    """ scans through fastg headers as an adjacency list
        builds and returns a nx directed graph using adjacencies
        note: no connections are created between each node and its
        rc node - we need to take care to maintain these
    """
    lines = []
    fp = open(fastg_name, 'r')
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1])
        lines.append(name)
    G = nx.DiGraph()
    return nx.parse_adjlist(lines, create_using=G)

def get_fastg_seqs_dict(fastg_name, G):
    """ returns a dictionary of sequences in graph
        where node names are keys and sequence strings
        are values; useful for saving memory when G
        is a subgraph (e.g., a component)
    """
    fp = open(fastg_name, 'r')
    seqs = {}
    for name,seq,qual in readfq(fp):
        name_parts = re.sub('[:,]'," ", name[:-1]).split()
        node = name_parts[0]
        seqs[node] = seq
    return seqs

def rc_node(node):
    """ gets reverse complement
        spades node label
    """
    if node[-1] == "'": return node[:-1]
    else: return node + "'"


def get_cov_from_spades_name_and_graph(name,G):
    if name not in G:
        return 0
    if 'cov' in G.node[name]:
        return G.node[name]['cov']
    else:
        return get_cov_from_spades_name(name)

def update_node_coverage(G, node, new_cov):
    """ changes coverage value stored in 'cov'
        field on both F and R version of a node
        if new_cov is 0, both node versions are removed
    """
    if node not in G.nodes(): # nothing to be done, perhaps already removed
        return
    if new_cov == 0:
        G.remove_node(node)
        if rc_node(node) in G.nodes():
            G.remove_node(rc_node(node))
    else:
        G.add_node(node, cov=new_cov)
        G.add_node(rc_node(node), cov=new_cov)

def get_spades_base_mass(G, name):
    length = get_length_from_spades_name(name)
    coverage = get_cov_from_spades_name_and_graph(name,G)
    if coverage <= 0.0: coverage = 1.0/float(length) # ensure no division by zero, consider more principled way to do this
    return length * coverage

def get_seq_from_path(path, seqs, max_k_val=55, cycle=True):
    """ retrieves sequence from a path;
        instead of specifying cycles by having the first and
        last node be equal, the user must decide if a cycle is
        the intent to avoid redundant k-mers at the ends
    """
    start = seqs[path[0]]
    if len(path)==1:
        if cycle:
            return start[max_k_val:]
        else:
            return start
    else:
        seq = ''
        for p in path:
            seq += seqs[p][max_k_val:]
        if cycle: return seq
        else: return start[:max_k_val] + seq

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
    if len(path)< 2: return 0
    mean, std = get_path_mean_std(path, G, seqs, max_k_val)
    if mean<=0: return 0
    return std/mean

def get_node_cnts_hist(path):
    d = {}
    for p in path:
        # always count based on positive node
        pos_name = p if (p[-1]!="'") else p[:-1]
        d[pos_name] = d.get(pos_name,0) + 1
    return d

def get_path_covs(path,G):
    covs = [get_cov_from_spades_name_and_graph(n,G) for n in path]
    cnts = get_node_cnts_hist(path)
    for i in range(len(path)):
        p = path[i]
        pos_name = p if (p[-1]!="'") else p[:-1]
        if cnts[pos_name] > 1:
            covs[i] /= cnts[pos_name]
    return covs


def get_path_mean_std(path, G, seqs, max_k_val=55):
    # covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    covs = get_path_covs(path,G)
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, cycle=True))
    if tot_len<=0: return (0,0)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs):
    mean, _ = get_path_mean_std(path, G, seqs)
    # covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    covs = get_path_covs(path,G)
    new_covs = covs - mean
    for i in range(len(path)):
        if new_covs[i] > 0:
            update_node_coverage(G,path[i],new_covs[i])
        else:
            update_node_coverage(G,path[i],0)

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])

def get_long_self_loops(G, min_length, seqs):
    """ returns set of self loop nodes paths that are longer
        than min length; removes those and short self loops
        from G
    """
    self_loops = set([])
    to_remove = []

    for nd in G.nodes_with_selfloops():
#TODO: Consider only doing this for isolated self-loops ?
        nd_path = (nd,)
        if len(get_seq_from_path(nd_path, seqs)) >= min_length \
        and (rc_node(nd),) not in self_loops:
            self_loops.add(nd_path)
#############
            logger.info("Added path: %s  - long self loop" % nd)
#        elif (rc_node(nd,)) not in self_loops:
#            logger.info('\tRemoving short self loop: %s' % nd) #########
        to_remove.append(nd)

    for nd in to_remove:
        update_node_coverage(G, nd, 0)
    return self_loops

def remove_hi_confidence_chromosome(G, len_thresh=10000, score_thresh=0.2):
    """ Remove the long nodes that are predicted to likely be chromosomal
    """
    to_remove = []
    for nd in G.nodes():
        if get_length_from_spades_name(nd) > len_thresh and \
            G.node[nd]['score']<score_thresh:
            to_remove.append(nd)
            to_remove.append(rc_node(nd)) ####################
            ############################# TODO: IS this really necessary?
    print "Removing {} nodes".format(len(to_remove))
    G.remove_nodes_from(to_remove)
    logger.info("Removed %d long likely chromosomal nodes" % len(to_remove))


def get_unoriented_sorted_str(path):
    """ creates unique, orientation-oblivious string representation of path,
        used to make sure node covered whenever rc of node is;
        lets us avoid issue of rc of node having different weight than node
    """
    all_rc_path = []
    for p in path:
        if p[-1] != "'": p = p+"'"
        all_rc_path.append(p)
    return "".join(sorted(all_rc_path))

def enum_high_mass_shortest_paths(G, use_scores=False, seen_paths=None):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes all shortest paths starting at each node n (assigning
        node weights to be 1/(length * coverage)) to each of
        its predecessors, and returning to n
    """
    if seen_paths == None:
        seen_paths = []
    nodes = []
    nodes = list(G.nodes()) # creates a copy

    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []

    logger.info("Getting edge weights")####################################################################

    # use add_edge to assign edge weights to be 1/mass of starting node
    for e in G.edges():
        #######################################################
        if use_scores:
            G.add_edge(e[0], e[1], cost = (1.-(G.node[e[0]]['score']))/get_spades_base_mass(G, e[0]))
        else:
            G.add_edge(e[0], e[1], cost = 1./get_spades_base_mass(G, e[0]))
    #    G.add_edge(e[0], e[1], cost = math.log(get_spades_base_mass(G, e[0])))
        ################################

######################################
    # edge weights on edges (v,x) are -log(\sum_u w(v,u))
#    for n in nodes:
#        outgoing_edges = G.edges(n,data=True)
#        norm_factor = sum([o[2]['cost'] for o in outgoing_edges])
#        for e in outgoing_edges:
#            G.add_edge(e[0],e[1], cost = math.log(norm_factor)-math.log(e[2]['cost'])) ############
############################################

    for node in nodes:

#	logger.info("%s: Paths to node: %s" % (multiprocessing.current_process().name, node))##############################################################

        # if node[-1] == "'": continue
        for pred in G.predecessors(node):
            # needed because some nodes removed on updates
            if pred not in G: continue

            try:
                path = nx.shortest_path(G, source=node,
                    target=pred, weight='cost')
            except nx.exception.NetworkXNoPath:
                continue


            # below: create copy of path with each node as rc version
            # use as unique representation of a path and rc of its whole
            unoriented_sorted_path_str = get_unoriented_sorted_str(path)

            # here we avoid considering cyclic rotations of identical paths
            # by sorting their string representations (all_rc_path above)
            # and comparing against the set already stored
            if unoriented_sorted_path_str not in unq_sorted_paths:
                unq_sorted_paths.add(unoriented_sorted_path_str)
                paths.append(tuple(path))

    return paths

def get_non_repeat_nodes(G, path):
    """ returns a list of all non-repeat (in degree and out-degree
        == 1) nodes in a path; if there are no such nodes,
        returns an empty list
        NB: G input should be whole graph, not specific SCC, to avoid
        disregarding isolated nodes
    """
    sing_nodes = []
    for nd in path:
        if G.out_degree(nd)==1 and G.in_degree(nd)==1:
            sing_nodes.append(nd)
    return sing_nodes


def get_spades_type_name(count, path, seqs, max_k_val, G, cov=None):
    path_len = len(get_seq_from_path(path,seqs,max_k_val))
    if cov==None:
        cov = get_total_path_mass(path,G)/float(path_len)
    info = ["RNODE", str(count+1), "length", str(path_len),
     "cov", '%.5f' % (cov)]
    return "_".join(info)

def get_contigs_of_mates(node, bamfile, G):
    """ retrieves set of nodes mapped to by read pairs
        having one mate on node; discards isolated nodes
        because they tend to reflect irrelevant alignments
    """
    mate_tigs = set([])
    if node[-1] == "'": node=node[:-1]
    try:
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                mate_tigs.add(nref)

    except ValueError:
        pass
    source_name = node #re.sub('NODE_','EDGE_', node)

#    print "before removal", mate_tigs
    to_remove = set([])
    for nd in mate_tigs:
        # flip name from "NODE_" prefix back to "EDGE_"
        # differs between contigs set and graph node names
        nd_name = nd #re.sub('NODE_','EDGE_', nd)
        if (G.in_degree(nd_name)==0 and G.out_degree(nd_name)==0) or \
        (not G.has_node(nd_name)):
            to_remove.add(nd)
        # see if nd reachable by node or vice-versa
        # try both flipping to rc and switching source and target
        elif not any([nx.has_path(G, source_name, nd_name), nx.has_path(G, rc_node(source_name),nd_name),
          nx.has_path(G, nd_name, source_name), nx.has_path(G, nd_name, rc_node(source_name))]):
            to_remove.add(nd)
    mate_tigs -= to_remove
    # print "after removal", mate_tigs

    return mate_tigs

def is_good_cyc(path, G, bamfile):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle, ignoring mappings
        to isolated nodes or non-reachable nodes
    """

    sing_nodes = get_non_repeat_nodes(G,path)
    logger.info("Checking path: %s", path) #######################
    for nd in sing_nodes:
        mate_tigs = get_contigs_of_mates(nd, bamfile, G)  #re.sub('EDGE_', 'NODE_' ,nd), bamfile, G)

        # mate_tigs_fixed_names = [re.sub('NODE_','EDGE_', x) for x in mate_tigs]
        # print mate_tigs_fixed_names
        # need to check against F and R versions of path nodes

        logger.info("\tNode: %s" % nd)
        logger.info("\t\tMates: %s" % ", ".join(mate_tigs))
        in_path = [x in path for x in mate_tigs]
        # print in_path
        path_rc = [rc_node(x) for x in path]
        in_rc_path = [x in path_rc for x in mate_tigs]
        # print in_rc_path
        if any([ (not in_path[i] and not in_rc_path[i]) for i in range(len(mate_tigs))]):
    	    logger.info("Mate not in path")
            return False
    return True

#########################
#TODO: calculate SEQS only for the component
#TODO: get rid of G and use COMP
# debugging - seem to be missing node in COMP??
def process_component(job_queue, result_queue, G, max_k, min_length, max_CV, SEQS, thresh, bampath, use_scores=False):
    """ run recycler for a single component of the graph
        use multiprocessing to process components in parallel
    """
    proc_name = multiprocessing.current_process().name
    bamfile = pysam.AlignmentFile(bampath)
    while True:
        COMP = job_queue.get()
        if COMP is None:
            print proc_name + ' is done'
            job_queue.task_done()
            break # done process
        logger.info("%s: Next comp" % (multiprocessing.current_process().name))
        # initialize shortest path set considered
        paths = enum_high_mass_shortest_paths(COMP,use_scores)
        logger.info("%s: Shortest paths: %s" % (multiprocessing.current_process().name, str(paths)))
        # peeling - iterate until no change in path set from
        # one iteration to next

        last_path_count = 0
        last_node_count = 0


        path_count = 0
        non_self_loops = set([])
        paths_set = set([]) #the set of paths found
        # continue as long as you either removed a low mass path
        # from the component or added a new path to final paths
        while(path_count!=last_path_count or\
            len(COMP.nodes())!=last_node_count):

            last_node_count = len(COMP.nodes())
            last_path_count = path_count

            if(len(paths)==0): break

            # using initial set of paths or set from last iteration
            # sort the paths by CV and test the lowest CV path for removal
            # need to use lambda because get_cov needs 2 parameters

            # make tuples of (CV, path)
            path_tuples = []
            for p in paths:
                path_tuples.append((get_wgtd_path_coverage_CV(p,COMP,SEQS,max_k_val=max_k), p))

            # sort in ascending CV order
            path_tuples.sort(key=lambda path: path[0])

            curr_path = path_tuples[0][1]
            logger.info("%s: Lowest CV path: %s, CV: %f, Total mass: %f" % (multiprocessing.current_process().name, str(curr_path),\
                        path_tuples[0][0],get_total_path_mass(curr_path,G)))

            if get_unoriented_sorted_str(curr_path) not in non_self_loops:
                path_mean, _ = get_path_mean_std(curr_path, COMP, SEQS, max_k_val=max_k)

                ## only report to file if long enough and good
                ## first good case - paired end reads on non-repeat nodes map on cycle
                ## typical or low coverage level
                ## second good case - high coverage (pairs may map outside due to high chimericism),
                ## near constant coverage level
                if (
                    len(get_seq_from_path(curr_path, SEQS, max_k_val=max_k))>=min_length \
                    and is_good_cyc(curr_path,G,bamfile) and \
                    get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k) <= (max_CV/len(curr_path))
                    ) or \
                (
                    len(get_seq_from_path(curr_path, SEQS, max_k_val=max_k))>=min_length and (path_mean > thresh) \
                    and get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k) <= (max_CV/len(curr_path))
                    ):
                    print(curr_path)

                    logger.info("Added path %s" % ", ".join(curr_path))
                    logger.info("\tCV: %4f" % get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k))
                    logger.info("\tCV cutoff: %4f" % (max_CV/len(curr_path)))
                    logger.info("\tPath mean cov: %4f" % path_mean) ####################
                    non_self_loops.add(get_unoriented_sorted_str(curr_path))

                    update_path_coverage_vals(curr_path, COMP, SEQS)
                    path_count += 1
                    paths_set.add(curr_path)

                else:
                    logger.info("%s: Did not add path: %s" % (proc_name, ", ".join(curr_path)))
                    logger.info("\tCV: %4f" % get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k))
                    logger.info("\tCV cutoff: %4f" % (max_CV/len(curr_path)))
                    logger.info("\tPath mean cov: %4f" % path_mean) ####################

                # recalculate paths on the component
                print(proc_name + ': ' + str(len(COMP.nodes())), " nodes remain in component\n")

                paths = enum_high_mass_shortest_paths(COMP,use_scores,non_self_loops)
                logger.info("%s: Shortest paths: %s" % (multiprocessing.current_process().name, str(paths)))

        job_queue.task_done()
        result_queue.put(paths_set)
    #end while
    return
