# Dowser, TIgGER, Alakazam and SHazaM on r-base, plus raxml-ng. Replaces the
# Immcantation suite image, which was 8.73 GB of which ~98% was never executed.
FROM immcantation/suite:4.8.0 AS suite

FROM r-base:4.4.2

ENV RENV_CONFIG_PPM_ENABLED=false
# Without this, every package is stored twice: once in renv's cache and once in the library.
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

RUN apt-get update && apt-get install -y \
    ca-certificates curl dash default-libmysqlclient-dev libbz2-dev libcairo2-dev \
    libcurl4-openssl-dev libfreetype6-dev libfribidi-dev libgit2-dev libglpk-dev \
    libgmp-dev libgsl-dev libharfbuzz-dev libjpeg-dev liblzma-dev libpcre2-dev \
    libpng-dev libpq-dev libsasl2-dev libsqlite3-dev libssh2-1-dev libssl-dev \
    libtiff5-dev libuv1-dev libxml2-dev unixodbc-dev unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# Absolute paths: the working directory changes at run time.
ENV RENV_PATHS_RENV=/app/renv
ENV RENV_PATHS_LOCKFILE=/app/renv.lock
ENV R_LIBS_USER=/app/renv/library
ENV RENV_PATHS_LIBRARY=${R_LIBS_USER}
RUN mkdir -p "${R_LIBS_USER}"

# restore() reads repositories out of the lockfile, so do not override them here.
COPY renv.lock /app/renv.lock
RUN R --no-echo -e "tryCatch(renv::restore(), error = function(e) message('parallel pass finished with errors: ', conditionMessage(e))); options(Ncpus = 1); renv::restore(clean = TRUE)" \
    && find /root/.cache/R -name keys.html -delete \
    && find /usr/local/lib/R -name keys.html -delete \
    && find "${R_LIBS_USER}" -name keys.html -delete \
    && rm -rf /root/.cache/R /usr/share/doc /usr/share/man

ENV RAXML_NG_EXEC=/usr/local/bin/raxml-ng
RUN curl -sSL -o /tmp/raxml.zip \
      https://github.com/amkozlov/raxml-ng/releases/download/2.0.3/raxml-ng_v2.0.3_linux_x86_64.zip \
    && unzip -q /tmp/raxml.zip -d /tmp/raxml \
    && mv /tmp/raxml/raxml-ng "${RAXML_NG_EXEC}" && chmod +x "${RAXML_NG_EXEC}" \
    && rm -rf /tmp/raxml.zip /tmp/raxml

# IgPhyML, built from source. The suite ships a Fedora binary that wants libsatlas
# and libcblas, which Debian does not provide. Recipe and the prototype patch come
# from the suite's Dockerfile; its BLAS variant needs ATLAS's clapack.h, absent on
# Debian, so this uses the OpenMP-only build. Placed after the R layer to keep that cache.
RUN apt-get update && apt-get install -y --no-install-recommends automake autoconf \
    && curl -sSL -o /tmp/igphyml.tar.gz \
       https://github.com/immcantation/igphyml/archive/refs/tags/2.0.0.tar.gz \
    && mkdir -p /usr/local/share/igphyml \
    && tar -zxf /tmp/igphyml.tar.gz -C /usr/local/share/igphyml --strip-components 1 \
    && sed -i 's/t_tree \*Make_Tree();/t_tree *Make_Tree(int n_otu);/' /usr/local/share/igphyml/src/utilities.h \
    && (cd /usr/local/share/igphyml && ./make_phyml_omp) \
    && rm -f /tmp/igphyml.tar.gz \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
ENV IGPHYML_EXEC=/usr/local/share/igphyml/src/igphyml

# IMGT germline reference sets, data only, so no binary compatibility concern.
COPY --from=suite /usr/local/share/germlines /usr/local/share/germlines

# FastTree drafts the topology that RAxML then fits; see raxml-ng-fasttree for
# why. Debian's build is double precision, which is what the short branches of a
# hypermutated lineage need to be distinguishable from zero.
ENV FASTTREE_EXEC=/usr/bin/fasttree
RUN apt-get update && apt-get install -y --no-install-recommends fasttree \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# What dowser is handed as its raxml-ng. It is raxml-ng for every call but the
# topology search, which it serves from FastTree instead.
ENV RAXML_EXEC=/usr/local/bin/raxml-ng-fasttree
COPY raxml-ng-fasttree "${RAXML_EXEC}"
RUN chmod +x "${RAXML_EXEC}"

COPY trees.R /app/trees.R
COPY alleles.R /app/alleles.R
COPY common.R /app/common.R

RUN Rscript -e "stopifnot(packageVersion('dowser') == '2.5.1'); suppressMessages({library(dowser); library(tigger); library(ggtree)})" \
 && "${RAXML_NG_EXEC}" --version \
 && "${FASTTREE_EXEC}" -expert 2>&1 | grep -q "Double precision" \
 && test -x "${RAXML_EXEC}" \
 && test -x "${IGPHYML_EXEC}"

# Grow R's heap aggressively. With millions of live objects every full collection
# walks all of them; measured on a 0.8 GB heap of small objects, a loop of small
# allocations spent 6.6 of its 10.9 s collecting by default and 0.7 of 5.1 s here.
ENV R_GC_MEM_GROW=3

ENTRYPOINT []
CMD ["Rscript", "/app/trees.R"]
