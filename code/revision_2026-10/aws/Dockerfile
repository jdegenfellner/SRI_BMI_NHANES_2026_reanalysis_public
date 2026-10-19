# Rechenumgebung fuer die Ableitung aus den 80-Hz-Rohdaten (GGIR 3.3-6, SWaN_accel 1.18); Versionen wie auf dem lokalen Mac.
FROM rocker/r-ver:4.6.1
ENV DEBIAN_FRONTEND=noninteractive MAKEFLAGS=-j2
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 curl ca-certificates perl python3-pip libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libuv1-dev pkg-config && rm -rf /var/lib/apt/lists/*
RUN pip3 install --break-system-packages uv && uv venv --python 3.8 /opt/swan38
COPY python_pins.txt /tmp/python_pins.txt
RUN uv pip install --python /opt/swan38/bin/python -r /tmp/python_pins.txt
COPY r_packages_pinned.tsv /tmp/r_packages_pinned.tsv
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(2, parallel::detectCores() - 1)); install.packages("remotes"); p <- read.delim("/tmp/r_packages_pinned.tsv"); for (i in seq_len(nrow(p))) tryCatch(remotes::install_version(p$package[i], p$version[i], upgrade = "never", quiet = TRUE), error = function(e) { message("FALLBACK ", p$package[i], ": ", conditionMessage(e)); install.packages(p$package[i]) }); ip <- installed.packages(); w <- ip[match(p$package, rownames(ip)), "Version"]; bad <- is.na(w) | w != p$version; d <- p[bad, ]; write.table(data.frame(package = d$package, wanted = d$version, installed = w[bad]), "/opt/r_version_differences.tsv", sep = "\t", quote = FALSE, row.names = FALSE); cat("Abweichende Versionen:", nrow(d), "\n")'
RUN mkdir -p /opt/pipeline
COPY ggir_one_cloud.R ggir_double_cumsum.R swan_one_cloud.sh swan_prep.R one_participant.sh queue_worker.sh /opt/pipeline/
RUN chmod +x /opt/pipeline/*.sh
ENV NH_PY=/opt/swan38/bin/python NH_ZCAT=zcat NH_SWAN_PREP=/opt/pipeline/swan_prep.R
