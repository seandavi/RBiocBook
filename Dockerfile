# Container image for "The RBioc Book".
#
# One image, three uses:
#   1. CI            — render the book with no per-run package install
#   2. devcontainer  — render/run ANY chapter locally (incl. the single-cell stack)
#   3. Orchestra     — browser-based RStudio workshop (RStudio Server ships in the base)
#
# Pinned to the Bioconductor RELEASE the book is validated against. The base is
# derived from rocker/rstudio, so it provides RStudio Server + R 4.6 + every
# Bioconductor system library + BINARY package installs (no source compiles).
ARG BIOC_VERSION=RELEASE_3_23
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}

# System tools Quarto needs that the base doesn't carry:
#   - librsvg2-bin: rsvg-convert, used to turn SVG figures into PDF
#   - python3-pip:  to install quartobot (the pre-render citation hook)
RUN apt-get update \
  && apt-get install -y --no-install-recommends librsvg2-bin python3-pip \
  && rm -rf /var/lib/apt/lists/*

# Install the book's R package dependencies (CRAN + Bioconductor) as binaries.
# DESCRIPTION Depends: is the single source of truth, kept honest by
# tools/check-deps.R. Version constraints are stripped with regex character
# classes (no backslashes, so this survives every quoting layer).
COPY DESCRIPTION /tmp/DESCRIPTION
RUN Rscript -e 'deps <- read.dcf("/tmp/DESCRIPTION")[1, "Depends"]; \
      deps <- trimws(gsub("[(].*?[)]", "", strsplit(deps, ",")[[1]])); \
      deps <- deps[deps != "" & deps != "R"]; \
      BiocManager::install(deps, update = FALSE, ask = FALSE, Ncpus = parallel::detectCores())' \
  && rm /tmp/DESCRIPTION

# Quarto CLI + TinyTeX (PDF) + quartobot (pre-render hook declared in _quarto.yml).
ARG QUARTO_VERSION=1.9.37
RUN curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
  && apt-get update && apt-get install -y --no-install-recommends /tmp/quarto.deb \
  && rm /tmp/quarto.deb && rm -rf /var/lib/apt/lists/*
# Install TinyTeX to a SYSTEM path so LaTeX (lualatex/tlmgr) is on PATH for any
# user and any $HOME. GitHub Actions sets HOME=/github/home inside container
# jobs, so a per-user `quarto install tinytex` (installed as root at build time)
# isn't found at render time; the RStudio/Orchestra user needs it too.
RUN Rscript -e 'install.packages("tinytex", repos = "https://cloud.r-project.org"); tinytex::install_tinytex(dir = "/opt/TinyTeX")'
ENV PATH="/opt/TinyTeX/bin/x86_64-linux:${PATH}"
RUN pip3 install --no-cache-dir --break-system-packages quartobot

# Node.js + npm (apt 18.x) so the devcontainer can install in-container agent
# tooling (Claude Code) at postCreate time. Appended last as its own layer so it
# doesn't invalidate the heavy dependency layers above. Not used by CI render or
# Orchestra.
RUN apt-get update && apt-get install -y nodejs npm && rm -rf /var/lib/apt/lists/*

# The nmfs-opensci titlepage extension is project-local (_extensions/, gitignored)
# and installed at render time, so it is intentionally not baked in here.
