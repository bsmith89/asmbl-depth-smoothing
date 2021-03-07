# Assembly Depth Smoothing

We simplify assembly graphs is in order to make longer contigs.
Among other things, one reason we need longer contigs is for accurate
read-mapping (and therefore depth estimation) of these contigs.
Unfortunately, besides the problems with bowtie2 style read-mapping in general,
this process also results in small, fragmented contigs being lost entirely.

Here I'm trying a different approach, which I'm hoping gets rid of both of
these problems.
I'm _not_ simplifying the assembly graph (so it's a _raw_ de Bruijn graph,
instead), and I'm getting depth estimates by counting median k-mer depths
for each segment.
Basically, this means that I'm doing assembly and mapping all in one step!

Using this approach, the problem becomes the estimation error of unitig (or
k-mer) depths, since they are sampled as a Poisson process.
To solve this, I use a algorithm to
estimate the depth of each unitig based on the depth of adjacent unitigs.
This algorithm converges in a way that "smooths" unitig depths based on their
neighbors, so even very short unitigs with few linear k-mers get a smoothed,
per-sample depth estimate.

Here's an illustration:

```
\                          /
 \ 5                      / 1
  \                      /
   \         2          /
    --------------------
   /                    \
  / 7                    \
 /                        \ 3
/                          \

```

Edges are unitigs (alternatively kmers) and number represent the
observed depth, or some estimate there-of[^median of a sample? Bloom-filter? Minimizers?].
The numbers here could reflect the super-position of three underlying genome
sequences:

```
\
 \ 5
  \               /
   \             / 1
    |           /
   /           /
  /           |
 / 5           \
/               \
                 \ 1
                  \

             2
    --------------------
   /                    \
2 /                      \ 2
 /                        \
/                          \

```

Notice that when all three of these are all sequenced together,
linear paths no longer share the same depth, because they're a combination of
multiple genomes.
The exact numbers on each path will depend on the relative coverage of each
underlying genome, and may vary across samples.
It should be possible to find a minimal set of paths and latent coverages
that fully explain the observed unitig depths across all samples.
NMF implies one very simple model for how these coverages might work,
allowing decomposition of the graph/matrix into genome components.
However, you can imagine, that if our sampling was not deterministic (and it
never is), that the coverage numbers will include noise making it harder
for NMF to factorize correctly.
This will be a problem especially in the case of low-coverage, short unitigs
which will have high coefficients of variation.

Instead of trying to build and fit a more sophisticated model,
(perhaps adding a contiguity constraint?)
we can instead adjust coverages heuristically based on their neighbors,
and then use the simple factorization approach.

The key to this heuristic is that each correctly constructed unitig
has two branches on either side.
(Here we're ignoring unitigs at the ends of linear molecules or where neighbors
are not observed in any sample.
--- We might imagine adding "phantom" unitigs at these ends with coverage of 0.
This is also an elegant link between the k-mer and unitig versions of this
explanation.))
For every such unitig, the coverage along it's length should be
(approximately) the same as the "incoming" coverage from each end.
When it's not, we can both update our coverage estimate and
suggest to our neighbors that they update their coverage estimates.
This is the crux of the message-passing approach to coverage smoothing.
(We hope that) The coverages eventually converge over success rounds of updates
and message passing.

The beauty of this approach is that we can build the cDBG just once (e.g. using
BCALM),
and then for each collection of reads (library) estimate the coverage on each
segment based on kmer counts (these counts could be stored efficiently in a
bloom filter or using minimizers).
We would then get coverage estimates for each unitig using the MP algorithm
sketched above.
Counting could be parallelized over libraries, and then
the MP algorithm could be parallelized over both libraries and connected
components in the cdBG (finding connected components is
probably hard[^see e.g. [@Flick2015] and papers citing it]).
The result would be big tables of exact unitig depths for each sample
without any mapping step and where we've incorporated
connectivity information directly.
We could then do NMF directly on these coverages and try to bin
genomes, for instance.
I _think_ one benefit would be that we could treat coverage estimates as
equivilantly precise, because everything would have similar numbers of kmers
informing it?

Not having to do any mapping to get coverage is goal enough on its own,
but is there a way to leverage this idea into longer contigs?
I think I've got to use some of the statistical path-finding approaches
(unzipping) that I had been looking at last year...
But maybe we can "save" the coverage information in the unzipping
process, so we unzip sub-paths _and_ give them coverages at the same time?
