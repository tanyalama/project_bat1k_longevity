# Export the official ICS/CGMW chronostratigraphic interval colours and boundary ages used by the
# geologic time scale under panel A of Figure 1. Everything comes from the deeptime package (2.4.0);
# no colour or boundary is transcribed by hand.
#   Rscript scripts/00_export_ics_intervals.R
library(deeptime)
library(jsonlite)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
if (length(here) == 0) here <- "scripts"
data(periods); data(epochs)
keep <- function(d) d[d$min_age < 100 & d$max_age > 0,
                      c("name", "abbr", "max_age", "min_age", "color", "lab_color")]
out <- list(periods = keep(periods), epochs = keep(epochs),
            deeptime_version = as.character(packageVersion("deeptime")))
write_json(out, file.path(here, "..", "data", "ics_intervals.json"), auto_unbox = TRUE, digits = 6)
cat("wrote data/ics_intervals.json:", nrow(keep(epochs)), "epochs,", nrow(keep(periods)), "periods\n")
