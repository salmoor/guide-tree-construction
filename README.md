#  Guide Tree Construction

Guide Tree Construction using [UPGMA](https://en.wikipedia.org/wiki/UPGMA) algorithm.

## Primer

Evolutionary tree is a branching diagram or "tree" showing the evolutionary relationships among various biological species based upon similarities and differences in their genetic characteristics.

UPGMA (unweighted pair group method with arithmetic mean) is a simple agglomerative (bottom-up) hierarchical clustering method.

## Installation

```bash
$ git clone https://github.com/salmoor/guide-tree-construction.git
$ cd guide-tree-construction
$ make
```

## How to use

Example invocation:

```bash
$ ./buildUPGMA --fasta sequences.fa --match 5 --mismatch -3 --gapopen -8 --gapext -1 --out sequences.tree
```

**Input** parameters are as follows:

- **−−fasta:** FASTA-formatted file containing all sequences. Up to 25 sequences
- **−−match:** match score
- **−−mismatch:** mismatch penalty
- **−−gapopen:** gap opening penalty
- **−−gapext:** gap extension penalty



**Output:**

- **−−out: **file that contains final output tree in [Newick](https://en.wikipedia.org/wiki/Newick_format) format. The following is the output for sequences in `test_sequences` folder with the above input options.

```bash
(((((A:1.00, C:1.00):0.75, B:1.75):0.75, F:2.50):2.50, D:5.00):1.00, E:6.00);
```