#' process raw STRING data
#'
#' given a STRING data file, process it and save it for future use
#'
#' @param string_file the data file
#' @export
#' @return data.frame
#'
process_string_data = function(string_file){
  stopifnot(file.exists(string_file))

  string_data = read.table(string_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
  return(string_data)
}

#' process raw STRING link data
#'
#' given a STRING links file, process it and save it for future use
#'
#' @param string_links the links file
#' @export
#' @return data.frame
#'
process_string_links = function(string_links){
  stopifnot(file.exists(string_links))

  string_ppi = read.table(string_links, header = TRUE, sep = " ", stringsAsFactors = FALSE)
  return(string_ppi)
}

#' process STRING id files
#'
#' given a STRING id file, generate a data.frame that can be used to map
#' various symbols
#'
#' @param string_file the alias file
#' @export
#' @return data.frame
#'
process_string_id = function(string_file){
  stopifnot(file.exists(string_file))

  string_id = read.table(string_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
  names(string_id) = c("string", "other", "type")
  return(string_id)
}

#' process STRING alias files
#'
#' given a STRING alias file, generate a data.frame that can be used to map
#' various symbols.
#'
#' @details The STRING alias file actually doesn't have a header, so the names we use
#'   here for the columns are:
#'
#'   * **string**: the STRING ID, as ORG1.ENSPXXXXXX, *9606.ENSP00000269305* is human, TP53
#'   * **other**: the other ID, whatever that is, *TP53* for this example
#'   * **type**: the source of the alias, *Ensembl_EntrezGene* in this case
#'
#' @param string_aliases the alias file
#' @export
#' @return data.frame
#'
process_string_aliases = function(string_aliases){
  stopifnot(file.exists(string_aliases))

  string_id = read.table(string_aliases, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
  names(string_id) = c("string", "other", "type")
  return(string_id)
}
