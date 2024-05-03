#debug.R
#remotes::install_github("extendr/rextendr")

rm(list=ls())


library(rextendr)

library(mergebamsR)



roxygen2::roxygenise()


rextendr::clean()
rextendr::document()

# Run once to configure package to use pkgdown
#usethis::use_pkgdown()
usethis::use_pkgdown_github_pages()
# Run to build the website


mergebams(bams = c("inst/extdata/test/bam1.bam","inst/extdata/test/bam2.bam"), out_path = "inst/extdata/test/out", prefixes = c("test1_","test2_"))
#works
peekbam(bam="inst/extdata/test/bam1.bam", n=3000)
#works

mergebams(bams = c("inst/extdata/test/bam1.bam","inst/extdata/test/bam2.bam"), 
          out_path = "test/out", 
          prefixes = c("test1_","test2_"), 
          names = list(c("VH00738:4:AAAW2TWHV:1:2512:20125:30230", "VH00738:4:AAAW2TWHV:1:1612:47827:5146"), 
                       c("VH00738:4:AAAW2TWHV:2:2114:49683:55957", "VH00738:4:AAAW2TWHV:2:1203:65210:13741")))

mergebams(bams = c("test/bam1.bam","test/bam2.bam"), 
          out_path = "test/out", 
          prefixes = c("test1_","test2_"), 
          names = list(c("VH00738:4:AAAW2TWHV:1:2512:20125:30230", "VH00738:4:AAAW2TWHV:1:1612:47827:5146"), 
                       NULL))



samtools view out_path.bam | tail -n 13012 | head -n 4

