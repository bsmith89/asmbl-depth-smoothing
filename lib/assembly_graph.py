from collections import defaultdict

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


def complement(a):
    return COMPLEMENT[a]


def reverse_complement(kmer):
    return ''.join(complement(a) for a in reversed(kmer))


def canonical(kmer):
    return min(kmer, reverse_complement(kmer))


def build_full_from_seed_graph(downstream_):
    upstream = defaultdict(list)
    downstream = defaultdict(list)
#     downstream.update(downstream_)
    for kmer in downstream_:
        rc_kmer = reverse_complement(kmer)
        for kmer_downstream in downstream_[kmer]:
            downstream[kmer].append(kmer_downstream)  # Kept for consistency; could have been done more efficiently with update above.
            upstream[kmer_downstream].append(kmer)
            rc_kmer_downstream = reverse_complement(kmer_downstream)
            downstream[rc_kmer_downstream].append(rc_kmer)
            upstream[rc_kmer].append(rc_kmer_downstream)
    assert mapping_all_upstream(upstream)
    return downstream, upstream


def add_reverse_complement_depth(depth_):
    depth = defaultdict(lambda: 0)
    for kmer in depth_:
        depth[kmer] = depth_[kmer]
        depth[reverse_complement(kmer)] = depth_[kmer]
    return depth


def is_ordered(upstream, downstream):
    return upstream[1:] == downstream[:-1]


def mapping_all_upstream(graph):
    for k in graph:
        for u in graph[k]:
            if not is_ordered(u, k):
                return False
    return True