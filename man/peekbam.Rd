% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergebamsR.R
\name{peekbam}
\alias{peekbam}
\title{Peek into a BAM File with Tag Filtering}
\usage{
peekbam(bam, n = 100, field = c("tag", "name"), TAG = "CB")
}
\arguments{
\item{bam}{A character string specifying the path to a single BAM file.
Only a single file should be specified.}

\item{n}{An integer, defaulting to 100, indicating the number of entries to peek.
\code{n} must be greater than 0.}

\item{TAG}{A character string specifying the tag to filter by within the BAM file.
The default is "CB" (Cell Barcode).}
}
\value{
The function itself does not return a value; it operates through side effects
such as invoking a Rust helper function or printing messages to the console.
}
\description{
This function provides a quick look into the contents of a BAM file, filtered by a specific tag.
It validates the input parameters, checks for the file's existence, and utilizes a helper function to
handle the BAM processing if the file exists. An error is raised if the input conditions are not met.
}
\details{
If the file does not exist, a message will be displayed. If \code{n} is less than 1
or if more than one BAM file is specified, the function will stop with an error.
}
\examples{
# Assuming 'example.bam' is a valid BAM file path:
peekbam("example.bam", n = 10, TAG = "CB")

}
\references{
This documentation was written by ChatGPT v4 - OpenAI, conversation with the author, 5-2-2024.
}
