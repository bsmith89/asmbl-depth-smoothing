# Assembly Depth Smoothing

We simplify assembly graphs is in order to make longer contigs.
Among other things, one reason we need longer contigs is for accurate
read-mapping (and therefore depth estimation) of these contigs.
Unfortunately, besides the problems with bowtie2 style read-mapping in general,
this process also results in small, fragmented contigs being lost entirely.

Here I'm trying a different approach, which I'm hoping gets rid of both of
these problems.
I'm not simplifying the assembly graph (so it's a _raw_ de Bruijn graph,
instead), and I'm getting depth estimates by counting median k-mer depths
for each segment.
Basically, this means that I'm doing assembly and mapping all in one step.

Using this approach, the problem becomes the estimation error of k-mer depths,
since they are sampled as a Poisson process.
To solve this, I use a (maybe some-day "message passing") algorithm to
estimate the depth of each kmer based on the depth of adjacent k-mers.
This algorithm converges in a way that "smooths" k-mer depths based on their
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

Edges are unitigs (alternatively kmers) and number represent the (median)
observed depth.
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
that fully explain the observed k-mer counts across all samples.
NMF implies one very simple model for how these coverages might work,
allowing decomposition of the graph/matrix into genome components.
However, you can imagine, that if our sampling was not deterministic (and it
never is), that the coverage numbers will include noise making it harder
for NMF to factorize correctly.
This will be a problem especially in the case of low-coverage, short unitigs
which will have high coefficients of variation.

Instead of trying to build and fit a more sophisticated mode,
(perhaps adding a contiguity constraint?)
we can instead adjust coverages heuristically based on their neighbors,
and then use the simple factorization approach.

The key to this heuristic is that each correctly constructed unitig
has two branches on either side.
(Here we're ignoring unitigs at the ends of linear molecules or where neighbors
are not observed in any sample.
We might imagine adding "phantom" unitigs at these ends with coverage of 0.
This is also an elegant link between the k-mer and unitig versions of this
explanation.)
For every such unitig, the coverage along it's length should be
(approximately) the difference between the coverage at either end.
In the schematic above, notice that for the flanks of the unitig with coverage
have a difference in coverage of 2.
When this is _not_ the case, we should assume that one of the observed
coverages is an imperfect estimate and adjust to account for this.

But what if it's one of the four flanks that has the bad coverage estimate??
The important thing to remember is that for each flank, our original unitig
is itself flanking, therefore the symmetric application of our update
rule will (hopefully) balance these two possibilities.

So what kind of a rule might we use? Here's a proposal:

1. For every unitig consider its two left and two right flanks.
    A. If the difference in coverage between the left flanks and the difference in coverage of the right flanks match our unitigs coverage do nothing.
    B. If the difference in coverage between the left/right flanks does _not_ match our unitigs coverage,
        i. If the two flanking pairs have very similar coverages, update the unitig coverage to be much closer to these two.
        ii. If the two flanking pairs do not have similar coverages, do not update the unitigs coverage as much.
        iii. (Make the adjustments above in a way that respects higher coverage more.)
2. Repeat 1. until convergence.

So I want an update formula based on our unitigs coverage at step $t$, $C_t$, and the difference
in coverage on the left and on the right, $L_t$ and $R_t$.
Will this actually converge? Who knows! Probably depends on our exact formula.

$C_{t+1} \leftarrow C_t + \frac{C_t + L_t + R_t}{3 * }

**NOTE: Just realized a lot of the above is wrong because unitigs that connect to this unitig
(usually) don't connect to each other.**
