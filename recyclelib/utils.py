from __future__ import division
import numpy as np
import math
import networkx as nx
import re, pysam
import logging
import multiprocessing

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
logger = logging.getLogger("recycle_logger")

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

def get_gene_nodes(genes_file,G):
    """ Annotate each node in the graph whether it has plasmid gene in it
    """
    gene_nodes = set()
    with open(genes_file) as f:
        for line in f:
            gene_nodes.add(line.strip())
    for nd in G.nodes():
        if nd in gene_nodes:
            G.add_node(nd, gene=True)
        else:
            G.add_node(nd,gene=False)



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
    if 'cov' in G.nodes[name]:
        return G.nodes[name]['cov']
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
        if rc_node(node) in G.nodes(): # HACKy TODO: FIX THIS BUG PROPERLY!
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
    mean, std = get_path_mean_std(path, G, seqs, max_k_val,discount=True)
    if mean<=0: return 0
    return std/mean

def get_node_cnts_hist(path):
    d = {}
    for p in path:
        # always count based on positive node
        pos_name = p if (p[-1]!="'") else p[:-1]
        d[pos_name] = d.get(pos_name,0) + 1
    return d

#######################################
def get_discounted_node_cov(node,path,G):
    """ Return the coverage of the node, discounted by the coverage of neighbouring
        nodes not in the path
    """
    ######################## Another try:
    pred_covs = [(get_cov_from_spades_name_and_graph(p,G),p) for p in G.predecessors(node)]
    succ_covs = [(get_cov_from_spades_name_and_graph(s,G),s) for s in G.successors(node)]

##    pred_in_path_cov = sum(p[0] for p in pred_covs if p[1] in path)
##    pred_non_path_cov = sum(p[0] for p in pred_covs if p[1] not in path)
##    pred_path_weight = pred_in_path_cov/(pred_in_path_cov+pred_non_path_cov)
##    succ_in_path_cov = sum(s[0] for s in succ_covs if s[1] in path)
##    succ_non_path_cov = sum(s[0] for s in succ_covs if s[1] not in path)
##    succ_path_weight = succ_in_path_cov/(succ_in_path_cov+pred_non_path_cov)
    non_path_cov = sum([p[0] for p in pred_covs if p[1] not in path]) + sum([s[0] for s in succ_covs if s[1] not in path])
    in_path_cov = sum([p[0] for p in pred_covs if p[1] in path]) + sum([s[0] for s in succ_covs if s[1] in path])
    node_cov = get_cov_from_spades_name_and_graph(node,G)
    node_cov *= in_path_cov/(non_path_cov + in_path_cov)
    ###################### A possible alternative would be to discount in this way for both the in- and out-neighbours
    ###################### and then average the two discounted - TODO: compare this
    return node_cov

######    pred_cov = sum([get_cov_from_spades_name_and_graph(p,G) for p in G.predecessors(node) if p not in path])
######    succ_cov = sum([get_cov_from_spades_name_and_graph(s,G) for s in G.successors(node) if s not in path])
######    avg_cov_discount = (pred_cov+succ_cov)/2.0
######    return max((get_cov_from_spades_name_and_graph(node,G)-avg_cov_discount),0.0)
#########################################


def get_path_covs(path,G,discount=False): ###############
#############
    if discount:
        covs = [get_discounted_node_cov(n,path,G) for n in path] #################
        # discount weight of nodes that path passes through multiple times
        cnts = get_node_cnts_hist(path)
        for i in range(len(path)):
            p = path[i]
            pos_name = p if (p[-1]!="'") else p[:-1]
            if cnts[pos_name] > 1:
                covs[i] /= cnts[pos_name] ########################################
    else:
        covs = [get_cov_from_spades_name_and_graph(n,G) for n in path]

    return covs


def get_path_mean_std(path, G, seqs, max_k_val=55,discount=True): ######################## TODO:
    ###################################################################################### should this be True?
    # covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    covs = np.array(get_path_covs(path,G,discount))
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, cycle=True))
    if tot_len<=0: return (0,0)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs):
    mean, _ = get_path_mean_std(path, G, seqs) ## NOTE: CAN WE STILL GUARANTEE CONVERGENCE WHEN DISCOUNTING COVERAGE ??!
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
            G.nodes[nd]['score']<score_thresh:
            to_remove.append(nd)
            to_remove.append(rc_node(nd))
    print "Removing {} nodes".format(len(set(to_remove)))
    G.remove_nodes_from(to_remove)
    logger.info("Removed %d long, likely chromosomal nodes" % len(set(to_remove)))

def get_hi_conf_plasmids(G, len_thresh=10000, score_thresh=0.9):
    """ Return a list of nodes that are likely plasmids
    """
#############
    try:
      hi_conf_plasmids = [nd for nd in G.nodes() if (get_length_from_spades_name(nd) > len_thresh and \
                        G.nodes[nd]['score'] > score_thresh)]
    except:
      print G.nodes()
      for nd in G.nodes():
        print nd
        if G.nodes[nd]['score'] > score_thresh: print nd
    logger.info("Found %d long, likely plasmid nodes" % len(hi_conf_plasmids))
    return hi_conf_plasmids

def get_plasmid_gene_nodes(G):
    """ Return list of nodes annotated as having a plasmid gene
    """
    plasmid_gene_nodes = [nd for nd in G.nodes() if G.nodes[nd]['gene']==True]
    logger.info("Found %d nodes with plasmid genes" % len(plasmid_gene_nodes))
    return plasmid_gene_nodes

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

def enum_high_mass_shortest_paths(G, use_scores=False, use_genes=False, seen_paths=None):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes all shortest paths starting at each node n (assigning
        node weights to be 1/(length * coverage)) to each of
        its predecessors, and returning to n
    """
    if seen_paths == None:
        seen_paths = []
    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []

    nodes = []
    nodes = list(G.nodes()) # creates a copy

    logger.info("Getting edge weights")

    # use add_edge to assign edge weights to be 1/mass of starting node
    #TODO: only calculate these if they haven't been/need to be updated
    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0) ##################################################
        elif use_scores==True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

            ########### NOTE: The shortest paths method counts shortest paths on the edges to the incoming neighbours of
            ################# each node. The weight of the source node is the same for all these paths and we need to count
            ################# the weight of the final node in each path. So - it makes more sense to put the weight of the
            ################# end node on each edge.

#############            G.add_edge(e[0], e[1], cost = (1.-(G.node[e[0]]['score']))/get_spades_base_mass(G, e[0]))


    logger.info("Getting shortest paths")
    #TODO: consider first running all pairs-shortest paths - is it more efficient?

    for node in nodes:

        shortest_score = float("inf")
        path = None
        for pred in G.predecessors(node):
            try:
                path_len = nx.shortest_path_length(G, source=node, target=pred, weight='cost')
                if path_len < shortest_score:
                    path = tuple(nx.shortest_path(G, source=node,target=pred, weight='cost'))
                    shortest_score = path_len
            except nx.exception.NetworkXNoPath:
                continue
        if path is None: continue

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



def get_high_mass_shortest_path(node,G,use_scores,use_genes):
    """ Return the shortest circular path back to node
    """
    #TODO: potentially add check for unique paths so that don't check same cycle
    # twice if there are two potential plasmid nodes in it

    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0) ##################################################
        elif use_scores == True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))
##########        G.add_edge(e[0], e[1], cost = (1.-(G.node[e[0]]['score']))/get_spades_base_mass(G, e[0]))

    shortest_score = float("inf")
    path = None
    for pred in G.predecessors(node):
        try:
            path_len = nx.shortest_path_length(G, source=node, target=pred, weight='cost')
            if path_len < shortest_score:
                path = tuple(nx.shortest_path(G, source=node,target=pred, weight='cost'))
                shortest_score = path_len
        except nx.exception.NetworkXNoPath:
            continue

    return path

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

#    print "before removal", mate_tigs
    to_remove = set([])
    for nd in mate_tigs:
        if (G.in_degree(nd)==0 and G.out_degree(nd)==0) or \
        (not G.has_node(nd)):
            to_remove.add(nd)
        # see if nd reachable by node or vice-versa
        # try both flipping to rc and switching source and target
        elif not any([nx.has_path(G, node, nd), nx.has_path(G, rc_node(node),nd),
                nx.has_path(G, nd, node), nx.has_path(G, nd, rc_node(node))]):
            to_remove.add(nd)
    mate_tigs -= to_remove

    return mate_tigs

#TODO: validate choice for ratio threshold
def is_good_cyc(path, G, bamfile,mate_ratio_thresh=0.5):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle
    """

    sing_nodes = set()
    for node in path:
      if node[-1] == "'": node = node[:-1]
      sing_nodes.add(node)
    ############sing_nodes = get_non_repeat_nodes(G,path)
    #TODO: instead of this, just create set of all (forward only) nodes in the path
    if len(sing_nodes)==0: return True ####################
#############    tot_in_path = 0
#############    tot_not_in_path = 0

#######################
    non_path_dominated_nodes = 0
#######################
    for nd in sing_nodes:
        mate_tigs = get_contigs_of_mates(nd, bamfile, G)
        # NOTE: ^ this only gets mates that are reachable from nd in G
        logger.info("\tNode: %s" % nd)
        logger.info("\t\tMates: %s" % ", ".join(mate_tigs))

        # need to check against F and R versions of path nodes
        path_rc = [rc_node(x) for x in path]
        num_mates_in_path = sum([1 for x in mate_tigs if (x in path or x in path_rc)])
        num_mates_not_in_path = len(mate_tigs)-num_mates_in_path
        if len(mate_tigs)>1 and num_mates_in_path < num_mates_not_in_path:
            non_path_dominated_nodes += 1
    if float(non_path_dominated_nodes)/float(len(sing_nodes)) > mate_ratio_thresh:
        logger.info("Too many nodes with majority of mates not on path")
        return False
    else: return True
    ###########################################################################
    #     if num_mates_in_path < num_mates_not_in_path:
    #         logger.info("Fewer mates in path than not")
    #     tot_in_path += num_mates_in_path # Note, this double counts in-path mates
    #     tot_not_in_path += num_mates_not_in_path
    #
    # # avoid divide by zeros
    # if tot_not_in_path == 0: return True
    # elif tot_in_path == 0:
    #     # TODO: or add this threshold in the other case too?
    #     if float(tot_not_in_path)/float(len(sing_nodes)) > 2*mate_ratio_thresh: # want to be stricter in this case
    #         return False
    #     else: return True
    #
    # if float(tot_not_in_path)/(float(tot_in_path)/2.0) > mate_ratio_thresh:
    #     logger.info("Too many mates not in path")
    #     return False
    # else: return True
    ###########################################################################


#########################
#TODO: calculate SEQS only for the component
#TODO: get rid of G and use COMP
def process_component(job_queue, result_queue, G, max_k, min_length, max_CV, SEQS, thresh, bampath, use_scores=False, use_genes=False):
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
        logger.info("%s: Next comp" % (proc_name))
        # initialize shortest path set considered

        path_count = 0
        seen_unoriented_paths = set([])
        paths_set = set([]) #the set of paths found

        # first look for paths starting from the nodes annotated with plasmid genes
        if use_genes:
            plasmid_gene_nodes = get_plasmid_gene_nodes(COMP)
            potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in plasmid_gene_nodes]
            potential_plasmid_mass_tuples.sort(key = lambda n: n[0])
            while potential_plasmid_mass_tuples: # could be removing other nodes from the list
                top_node = potential_plasmid_mass_tuples.pop() # highest mass node
                top_node_name = top_node[1]
##########################################################
##########################################################
                path = get_high_mass_shortest_path(top_node_name,COMP,use_scores,use_genes) #######
                if path is None: continue
                # check coverage variation
                path_CV = get_wgtd_path_coverage_CV(path,COMP,SEQS,max_k_val=max_k)
                logger.info("%s: Plasmid gene path: %s, CV: %4f" % (proc_name, str(path),path_CV))

                if path_CV <= max_CV and is_good_cyc(path,G,bamfile):
                    logger.info("%s: Added plasmid gene path %s" % (proc_name,str(path)))

                    # prevent checking nodes that have been removed
                    # TODO: do this more efficiently
                    # TODO: the right way to do this is remove just the top node and its RC
                    #       and recalculate the masses and sort in each iteration
                    i = 0
                    while i < len(potential_plasmid_mass_tuples):
                        if potential_plasmid_mass_tuples[i][1] in path or \
                            rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                            potential_plasmid_mass_tuples.pop(i)
                        else: i += 1

                    seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                    update_path_coverage_vals(path, COMP, SEQS)
                    path_count += 1
                    paths_set.add(path)

                else:
                    logger.info("%s: Did not add plasmid gene path: %s" % (proc_name, str(path)))



        # then look for circular paths that start from hi confidence plasmid nodes
        # TODO: implement own scoring function
        if use_scores:
            potential_plasmid_nodes = get_hi_conf_plasmids(COMP)
            potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in potential_plasmid_nodes]
            potential_plasmid_mass_tuples.sort(key = lambda n: n[0])
            while potential_plasmid_mass_tuples: # could be removing other nodes from the list
                top_node = potential_plasmid_mass_tuples.pop() # highest mass node
                top_node_name = top_node[1]
##########################################################
##########################################################
                path = get_high_mass_shortest_path(top_node_name,COMP,use_scores,use_genes) #######
                if path is None: continue
                # check coverage variation
                path_CV = get_wgtd_path_coverage_CV(path,COMP,SEQS,max_k_val=max_k)
                logger.info("%s: Hi conf path: %s, CV: %4f" % (proc_name, str(path),path_CV))

                if path_CV <= max_CV and is_good_cyc(path,G,bamfile):
                    logger.info("%s: Added hi conf path %s" % (proc_name,str(path)))

                    # prevent checking nodes that have been removed
                    # TODO: do this more efficiently
                    # TODO: the right way to do this is remove just the top node and its RC
                    #       and recalculate the masses and sort in each iteration
                    i = 0
                    while i < len(potential_plasmid_mass_tuples):
                        if potential_plasmid_mass_tuples[i][1] in path or \
                            rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                            potential_plasmid_mass_tuples.pop(i)
                        else: i += 1

                    seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                    update_path_coverage_vals(path, COMP, SEQS)
                    path_count += 1
                    paths_set.add(path)

                else:
                    logger.info("%s: Did not add hi-conf path: %s" % (proc_name, str(path)))

        # 3rd step. Run Recycler algorithm that looks for circular high mass shortest
        # paths and accept them as plasmid predictions if the coverages and mate pairs
        # match the required thresholds
#######################################################################################
#######################################################################################
        paths = enum_high_mass_shortest_paths(COMP,use_scores,use_genes,seen_unoriented_paths)
        ###################################logger.info("%s: Shortest paths: %s" % (proc_name, str(paths)))

        last_path_count = 0
        last_node_count = 0

        # continue as long as you either removed a low mass path
        # from the component or added a new path to final paths
        while(path_count!=last_path_count or\
            len(COMP.nodes())!=last_node_count):

            last_node_count = len(COMP.nodes())
            last_path_count = path_count

            # make tuples of (CV, path)
            path_tuples = []
            for p in paths:
                if len(get_seq_from_path(p, SEQS, max_k_val=max_k)) < min_length:
                    seen_unoriented_paths.add(get_unoriented_sorted_str(p))
                    logger.info("%s: Num seen paths: %d" % (proc_name, len(seen_unoriented_paths)))
                    continue
                path_tuples.append((get_wgtd_path_coverage_CV(p,COMP,SEQS,max_k_val=max_k), p))

            logger.info("%s: Num path tuples: %d" % (proc_name, len(path_tuples)))
            if(len(path_tuples)==0): break

            # sort in ascending CV order
            path_tuples.sort(key=lambda path: path[0])

            for pt in path_tuples:
                curr_path = pt[1]
                curr_path_CV = pt[0]
                logger.info("%s: Path: %s" % (proc_name, ",".join(curr_path)))
                if get_unoriented_sorted_str(curr_path) not in seen_unoriented_paths:

                ## only report if low CV and matches mate pair info
###################################3
                    if (curr_path_CV <= (max_CV) and \
                        is_good_cyc(curr_path,G,bamfile)):

                        print(curr_path)

                        logger.info("Added path %s" % ", ".join(curr_path))
                        logger.info("\tCV: %4f" % curr_path_CV)
                        seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))
                        logger.info("%s: Num seen paths: %d" % (proc_name, len(seen_unoriented_paths)))

                        update_path_coverage_vals(curr_path, COMP, SEQS)
                        path_count += 1
                        paths_set.add(curr_path)
                        break

                    else:
                        logger.info("%s: Did not add path: %s" % (proc_name, ", ".join(curr_path)))
                        logger.info("\tCV: %4f" % curr_path_CV)
                        if curr_path_CV > max_CV:
                            break # sorted by CV
                        else: # not good mate pairs
                            seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))
                            logger.info("%s: Num seen paths: %d" % (proc_name, len(seen_unoriented_paths)))

            # recalculate paths on the component
            print proc_name + ': ' + str(len(COMP.nodes())) + " nodes remain in component"
            logger.info("%s: Remaining nodes: %d" % (proc_name, len(COMP.nodes())))
            paths = enum_high_mass_shortest_paths(COMP,use_scores,use_genes,seen_unoriented_paths)
            logger.info("%s: Shortest paths: %s" % (proc_name, str(paths)))

        job_queue.task_done()
        result_queue.put(paths_set)
    #end while
    return
