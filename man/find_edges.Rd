% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodes_n_hops.R
\name{find_edges}
\alias{find_edges}
\title{find links between nodes}
\usage{
find_edges(
  in_graph,
  start_nodes,
  n_hop = 1,
  end_nodes = NULL,
  drop_same_after_2 = TRUE
)
}
\arguments{
\item{in_graph}{a graphBAM graph}

\item{start_nodes}{which nodes to start from}

\item{n_hop}{how many hops to go out (default is 1)}

\item{end_nodes}{optional, only keep edges that end at these nodes}

\item{drop_same_after_2}{should edges that only go to a single other node
be dropped? Default is TRUE. See DETAILS for more information.}
}
\value{
list
}
\description{
given a data frame of edges, find nodes within so many edges of initial nodes
with option that final set of edges go to known nodes
}
\details{
One frequently encountered case will be a set of nodes that do:
N1 --> N2 --> N1 (because the network is assumed to be undirected), where
after a single hop from N1 we reach N2, and then a second hop returns to N1.
Even in the case where the end_nodes are the same as the start_nodes, this
is likely \emph{not} useful information. So \code{drop_same_after_1 = TRUE}
will set the results so that these edge paths are not returned and kept
in the final graph.
}
\examples{
library(STRINGDatabaseManipulation)
library(graph)
set.seed(1234)
link_data = STRING10_links
link_data = link_data[sample(nrow(link_data), 10000),]
link_graph = string_2_graphBAM(link_data)
start_nodes = sample(nodes(link_graph), 10)
end_nodes = NULL
n_hop = 3
find_edges(link_graph, start_nodes, n_hop)

find_edges(link_graph, start_nodes, n_hop, end_nodes, drop_same_after_2 = FALSE)
}
