
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STRINGDatabaseManipulation

The goal of STRINGDatabaseManipulation is to provide functions for
working with STRING protein-protein interaction networks.

## Installation

Only the development version of STRINGDatabaseManipulation is available.
You will also need the `graph` package from Bioconductor.

``` r
install.packages("BiocManager")
BiocManager::install("graph")
remotes::install_github("rmflight/STRINGDatabaseManipulation")
```

## Issues

If you have a question about how to use the package, a request for
something new to be implemented, or have a bug, please [file an
issue](https://github.com/rmflight/STRINGDatabaseManipulation/issues).

## Getting Data

This package works with preprocessed data from STRING itself.

You can preprocess files from STRING by downloading detail and alias
files.

For example, for human you can download the detail file by:

    # get the PPI network itself
    wget https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
    # get aliases
    wget https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz

Note: these links become available after filtering by first choosing an
organism.

And then process them so they are easier to use:

``` r
library(STRINGDatabaseManipulation)
# assuming data is in the file, you would do:
ppi_data = process_string_data("9606.protein.links.detailed.v11.0.txt.gz")

# save it for later
saveRDS(ppi_data, file = "tmp_ppi_data.rda")

# process aliases
protein_aliases = process_string_id("9606.protein.aliases.v11.0.txt.gz")

# save for later
saveRDS(protein_aliases, file = "tmp_aliases.rda")
```

I have smaller versions of the data files saved for examples.

``` r
library(STRINGDatabaseManipulation)
links_file = system.file("extdata", "STRING11_9606_links_raw.txt.gz", package = "STRINGDatabaseManipulation")
ppi_data = process_string_links(links_file)
head(ppi_data)
#>               protein1             protein2 neighborhood fusion cooccurence coexpression
#> 1 9606.ENSP00000000233 9606.ENSP00000382239            0      0           0          143
#> 2 9606.ENSP00000000412 9606.ENSP00000463393            0      0           0           44
#> 3 9606.ENSP00000001146 9606.ENSP00000242208            0      0           0            0
#> 4 9606.ENSP00000001146 9606.ENSP00000225235            0      0           0            0
#> 5 9606.ENSP00000002165 9606.ENSP00000485444            0      0           0            0
#> 6 9606.ENSP00000003084 9606.ENSP00000351777            0      0           0           56
#>   experimental database textmining combined_score
#> 1            0        0         58            158
#> 2            0        0        160            162
#> 3            0        0        166            166
#> 4            0        0        223            223
#> 5            0      900          0            900
#> 6          408      900        450            965

aliases_file = system.file("extdata", "STRING11_9606_aliases_raw.txt.gz", package = "STRINGDatabaseManipulation")
head(STRING11_9606_aliases)
#>                 string other
#> 1 9606.ENSP00000361965   100
#> 2 9606.ENSP00000269141  1000
#> 3 9606.ENSP00000263826 10000
#> 4 9606.ENSP00000351484 10004
#> 5 9606.ENSP00000217455 10005
#> 6 9606.ENSP00000365312 10006
#>                                                                                                                                                                                 type
#> 1 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
#> 2 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
#> 3 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
#> 4 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
#> 5 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
#> 6 BLAST_UniProt_DR_GeneID Ensembl_HGNC_Entrez_Gene_ID Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DR_GeneID Ensembl_UniProt_DR_GeneID
```

## Using Data

Now lets actually do something with the STRING data. Lets find all nodes
within so 3 hops of a set of starting nodes.

We will use the example data from the package.

``` r
set.seed(1234)
ppi_graph = string_2_tidygraph(ppi_data)
all_nodes = ppi_graph |>
  tidygraph::activate(nodes) |>
  tibble::as_tibble() |>
  dplyr::pull(name)
start_nodes = all_nodes |>
  sample(size = 10)

n_hops = 3

after_3 = find_nodes_n_hops(ppi_graph, n_hops = n_hops, start_nodes = start_nodes)
after_3$graph
#> # A tbl_graph: 263 nodes and 264 edges
#> #
#> # An undirected simple graph with 3 components
#> #
#> # Node Data: 263 × 1 (active)
#>   name                
#>   <chr>               
#> 1 9606.ENSP00000037502
#> 2 9606.ENSP00000162749
#> 3 9606.ENSP00000204726
#> 4 9606.ENSP00000211998
#> 5 9606.ENSP00000215832
#> 6 9606.ENSP00000217233
#> # … with 257 more rows
#> #
#> # Edge Data: 264 × 3
#>    from    to weight
#>   <int> <int>  <dbl>
#> 1     1    23      2
#> 2     2    32      2
#> 3     2    81      2
#> # … with 261 more rows
```

In this case nothing changed, but depending on the inputs, it might.

## Visualizing

The nice thing about using `tidygraph`, is that we have access to the
`ggraph` library for visualization!

## Code of Conduct

Note that the ‘STRINGDatabaseManipulation’ project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.
