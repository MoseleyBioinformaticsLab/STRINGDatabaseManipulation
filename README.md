
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STRINGDatabaseManipulation

The goal of STRINGDatabaseManipulation is to provide functions for
reading in and manipulating [STRING](https://string-db.org/)
protein-protein interaction networks.

## Installation

Only the development version of STRINGDatabaseManipulation is available.

``` r
remotes::install_github("moseleyBioinformaticsLab/STRINGDatabaseManipulation")
```

## Issues

If you have a question about how to use the package, a request for
something new to be implemented, or have a bug, please [file an
issue](https://github.com/moseleyBioinformaticsLab/STRINGDatabaseManipulation/issues).

## Getting Data

This package works with preprocessed data from STRING itself.

You can preprocess files from STRING by downloading link detail and
alias files. These files are available from the [STRING download
page](https://string-db.org/cgi/download). I **do not** recommend trying
to work with the **full** interaction network across all species, but
rather download files for your specific species of interest. For
example, you can filter to just the [human
files](https://string-db.org/cgi/download?&species_text=Homo+sapiens).

We can get the files we want at the command line using a tool like
`wget`:

    # get the PPI network itself
    wget https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz
    # get aliases
    wget https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz

And then process them so they are easier to use:

``` r
library(STRINGDatabaseManipulation)
# assuming data is in the file, you would do:
ppi_data = process_string_links("9606.protein.links.detailed.v11.5.txt.gz")

# save it for later
saveRDS(ppi_data, file = "tmp_ppi_data.rda")

# process aliases
protein_aliases = process_string_aliases("9606.protein.aliases.v11.5.txt.gz")

# save for later
saveRDS(protein_aliases, file = "tmp_aliases.rda")
```

This package has smaller versions of the link and data files (10,000
links) saved for examples, both as raw text files, and as package data
files. These files have been filtered to those links with an
**experimental** score \>= 400, as well as links that include EGFR and
TP53.

``` r
library(STRINGDatabaseManipulation)
links_file = system.file("extdata", "STRING11.5_9606_links_raw.txt.gz", package = "STRINGDatabaseManipulation")
ppi_data = process_string_links(links_file)
head(ppi_data)
#>               protein1             protein2 neighborhood fusion cooccurence coexpression
#> 1 9606.ENSP00000001008 9606.ENSP00000354558            0      0           0           69
#> 2 9606.ENSP00000003084 9606.ENSP00000306330            0      0           0            0
#> 3 9606.ENSP00000003084 9606.ENSP00000449404            0      0           0            0
#> 4 9606.ENSP00000005257 9606.ENSP00000353590            0      0           0           51
#> 5 9606.ENSP00000005340 9606.ENSP00000262320            0      0           0           86
#> 6 9606.ENSP00000005340 9606.ENSP00000335677            0      0           0            0
#>   experimental database textmining combined_score
#> 1          835        0        343            890
#> 2          810        0         84            818
#> 3          550        0        875            941
#> 4          653        0         64            664
#> 5          987      900        475            999
#> 6          476        0         87            501

aliases_file = system.file("extdata", "STRING11.5_9606_aliases_raw.txt.gz", package = "STRINGDatabaseManipulation")
ppi_aliases = process_string_aliases(aliases_file)
head(ppi_aliases)
#>                 string   other                                          type
#> 1 9606.ENSP00000001008    2288                   Ensembl_HGNC_Entrez_Gene_ID
#> 2 9606.ENSP00000001008    2288 Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI)
#> 3 9606.ENSP00000001008 5.2.1.8                   BLAST_UniProt_DE_RecName_EC
#> 4 9606.ENSP00000001008   FKBP4                            Ensembl_EntrezGene
#> 5 9606.ENSP00000001008  FKBP51                    Ensembl_EntrezGene_synonym
#> 6 9606.ENSP00000001008  FKBP52                    Ensembl_EntrezGene_synonym
```

## Using Data

Now lets actually do something with the STRING data. The most common
analysis we want to do is find proteins that interact (directly or
indirectly) with one or more query proteins.

We will use the example data from the package.

To filter the data, we can use `dplyr` to choose which evidence or set
of evidences to use a filter. Here we will use `combined_score` \>= 400
(it’s actually already filtered, but this is to show how we can use it).

``` r
ppi_filtered = ppi_data |>
  dplyr::filter(combined_score >= 400)
dim(ppi_data)
#> [1] 21070    10
dim(ppi_filtered)
#> [1] 21070    10
```

In this case it doesn’t change the number of interactions, because the
data was pre-filtered to make it tractable for inclusion in the package.

## TP53

Lets take everyone’s favorite cancer gene, TP53, and look for those
proteins that experimentally are known to interact with it. We will find
the STRING-db ID, and then query everything that is connected to it.

``` r
ppi_graph = string_2_tidygraph(ppi_data)
tp53_alias = ppi_aliases |>
  dplyr::filter(other %in% "TP53")
tp53_alias
#>                 string other               type
#> 1 9606.ENSP00000269305  TP53 Ensembl_EntrezGene
```

Now we can go fetch the neighbors of our query protein (n_hops = 0), and
find everything that it interacts with. Notice, `n_hops = 0`! This is
because the `hops` refers to hops over other proteins. To find just the
interacting pairs, we are doing **0** hops.

``` r
tp53_interactions = find_nodes_n_hops(ppi_graph, n_hops = 0, start_nodes = tp53_alias$string)
tp53_interactions$graph
#> # A tbl_graph: 322 nodes and 479 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 322 × 1 (active)
#>   name                
#>   <chr>               
#> 1 9606.ENSP00000005340
#> 2 9606.ENSP00000011619
#> 3 9606.ENSP00000025008
#> 4 9606.ENSP00000156084
#> 5 9606.ENSP00000212015
#> 6 9606.ENSP00000215754
#> # … with 316 more rows
#> #
#> # Edge Data: 479 × 3
#>    from    to weight
#>   <int> <int>  <dbl>
#> 1     1    51      2
#> 2     1    68      2
#> 3     2    68      2
#> # … with 476 more rows
```

Here we can see that the returned interactions with TP53 includes
another 321 proteins. This is in **stark** contrast to the **maximum of
50** first shell entries returned by the STRING web tool. In addition to
the large number of interacting proteins, notice that many of them are
connected to each other as well, given the number of edges compared to
the number of nodes.

We can also examine, for this limited network, how many genes TP53
interacts with when we allow one protein in between. For this example,
this is likely to be lower than the real number, as it is a limited
dataset. In reality, TP53 interactors with a single hop is an incredibly
large number of interactors that you will want to watch your memory
usage.

``` r
tp53_onehop = find_nodes_n_hops(ppi_graph, n_hops = 1, start_nodes = tp53_alias$string)
tp53_onehop$graph
#> # A tbl_graph: 1275 nodes and 2461 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 1,275 × 1 (active)
#>   name                
#>   <chr>               
#> 1 9606.ENSP00000001008
#> 2 9606.ENSP00000003084
#> 3 9606.ENSP00000005340
#> 4 9606.ENSP00000011619
#> 5 9606.ENSP00000011653
#> 6 9606.ENSP00000025008
#> # … with 1,269 more rows
#> #
#> # Edge Data: 2,461 × 3
#>    from    to weight
#>   <int> <int>  <dbl>
#> 1     1   676      2
#> 2     2   386      2
#> 3     3   149      2
#> # … with 2,458 more rows
```

Even with the limited network, we end up with a whopping 1274 proteins
that interact one hop out from TP53!

## To Do

There are several improvements we have in mind to make to make this
package more useful:

- [functions to easily add protein identifiers to the returned
  network.](https://github.com/MoseleyBioinformaticsLab/STRINGDatabaseManipulation/issues/3)
- [take pathway data from the `graphite` package and work with it, we
  aren’t actually limited to STRING
  data.](https://github.com/MoseleyBioinformaticsLab/STRINGDatabaseManipulation/issues/5)
- [visualization
  examples.](https://github.com/MoseleyBioinformaticsLab/STRINGDatabaseManipulation/issues/2)
- [tests. We really need more
  tests.](https://github.com/MoseleyBioinformaticsLab/STRINGDatabaseManipulation/issues/1)
- [could we use `sqlite` databases to enable analysis of things with
  **really** large numbers of
  interactors?](https://github.com/MoseleyBioinformaticsLab/STRINGDatabaseManipulation/issues/4)

## Code of Conduct

Note that the ‘STRINGDatabaseManipulation’ project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.
