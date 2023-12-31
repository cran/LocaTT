% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trim_sequences.R
\name{trim_sequences}
\alias{trim_sequences}
\title{Trim Target Nucleotide Sequence from DNA Sequences}
\usage{
trim_sequences(
  sequences,
  target,
  anchor = "start",
  fixed = TRUE,
  required = TRUE,
  quality_scores
)
}
\arguments{
\item{sequences}{A character vector of DNA sequences to trim.}

\item{target}{A string specifying the target nucleotide sequence.}

\item{anchor}{A string specifying whether the target nucleotide sequence should be trimmed from the start or end of the DNA sequences. Allowable values are \code{"start"} (the default) and \code{"end"}.}

\item{fixed}{A logical value specifying whether the position of the target nucleotide sequence should be fixed at the ends of the DNA sequences. If \code{TRUE} (the default), then the position of the target nucleotide sequence is fixed at either the start or end of the DNA sequences, depending on the value of the anchor argument. If \code{FALSE}, then the target nucleotide sequence is searched for anywhere in the DNA sequences.}

\item{required}{A logical value specifying whether trimming is required. If \code{TRUE} (the default), then sequences which could not be trimmed are returned as \code{NA}s. If \code{FALSE}, then untrimmed sequences are returned along with DNA sequences for which trimming was successful.}

\item{quality_scores}{An optional character vector of DNA sequence quality scores. If supplied, these will be trimmed to their corresponding trimmed DNA sequences.}
}
\value{
If quality scores are not provided, then a character vector of trimmed DNA sequences is returned. If quality scores are provided, then a list containing two elements is returned. The first element is a character vector of trimmed DNA sequences, and the second element is a character vector of quality scores which have been trimmed to their corresponding trimmed DNA sequences.
}
\description{
Trims a target nucleotide sequence from the front or back of DNA sequences. Ambiguous nucleotides in the target nucleotide sequence are supported.
}
\details{
For each DNA sequence, the target nucleotide sequence is searched for at either the front or back of the DNA sequence, depending on the value of the anchor argument. If the target nucleotide sequence is found, then it is removed from the DNA sequence. If the required argument is set to \code{TRUE}, then DNA sequences in which the target nucleotide sequence was not found will be returned as \code{NA}s. If the required argument is set to \code{FALSE}, then untrimmed DNA sequences will be returned along with DNA sequences for which trimming was successful. Ambiguous nucleotides in the target nucleotide sequence are supported through the internal use of the \code{\link[=substitute_wildcards]{substitute_wildcards}} function on the target nucleotide sequence, and a regular expression with a leading or ending anchor is used to search for the target nucleotide sequence in the DNA sequences. If the fixed argument is set to \code{FALSE}, then any number of characters are allowed between the start or end of the DNA sequences and the target nucleotide sequence. Trimming will fail for DNA sequences which contain ambiguous nucleotides (\emph{e.g.}, Ns) in their target nucleotide sequence region, resulting in \code{NA}s for those sequences if the required argument is set to \code{TRUE}.
}
\examples{
trim_sequences(sequences=c("ATATAGCGCG","TGCATATACG","ATCTATCACCGC"),
               target="ATMTA",
               anchor="start",
               fixed=TRUE,
               required=TRUE,
               quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"))
}
