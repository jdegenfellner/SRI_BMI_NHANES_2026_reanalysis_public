# Records the computing environment of a run (hardware, operating system, R and package versions).
# Usage: Rscript system_info.R [output_file]   (default: system_info.txt next to the working directory)
out <- commandArgs(trailingOnly = TRUE)
out <- if (length(out)) out[1] else "system_info.txt"
hw <- tryCatch(system2("system_profiler", "SPHardwareDataType", stdout = TRUE), error = function(e) character())
hw <- grep("Model Name|Model Identifier|Chip|Total Number of Cores|Memory", hw, value = TRUE)
sink(out)
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
cat("Hardware\n"); cat(trimws(hw), sep = "\n")
cat("\nOS:", system2("sw_vers", "-productVersion", stdout = TRUE), "(", Sys.info()[["machine"]], ")\n\n")
print(sessionInfo())
sink()
cat("written:", out, "\n")
