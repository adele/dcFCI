# dcFCI R Package

Implementation of the data-compatible Fast Causal Inference (dcFCI) algorithm by Ribeiro &amp; Heider, 2025, for robust causal discovery under latent confounding, unfaithfulness, and mixed data. 

[Paper on ArXiv](https://www.arxiv.org/pdf/2505.06542)

### Citation

If you use this work, please cite:

Ribeiro, A. H., & Heider, D. (2025). dcFCI: Robust Causal Discovery Under Latent Confounding, Unfaithfulness, and Mixed Data. ArXiv preprint arXiv:2505.06542. https://arxiv.org/abs/2505.06542

```
@misc{ribeiro2025dcfcirobustcausaldiscovery,
      title={dcFCI: Robust Causal Discovery Under Latent Confounding, Unfaithfulness, and Mixed Data}, 
      author={Adèle H. Ribeiro and Dominik Heider},
      year={2025},
      eprint={2505.06542},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2505.06542}, 
}
```

---

## 📦 Step-by-Step Installation Instructions

This guide walks you through all necessary steps to install dependencies and set up the package environment in **R**, including how to install archived versions of required packages.


### 1. Install Bioconductor Dependencies

```r
install.packages("BiocManager")
BiocManager::install(c("RBGL", "graph", "Rgraphviz"))
```

---

### 2. Install CRAN Dependencies

```r
package_list <- c(
  "lmtest", "pscl", "brms", "dagitty", "ggm", "igraph", 
  "pcalg", "SEMgraph", "doFuture", "DOT", "jsonlite", "rsvg"
)

install.packages(package_list, dependencies = TRUE, repos = "http://cran.us.r-project.org")

# Install any missing packages
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if (length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE, repos = "http://cran.us.r-project.org")
}
```

---

### 3. Install Archived Version of `BFF` (v3.0.1)

```bash
wget https://cran.r-project.org/src/contrib/Archive/BFF/BFF_3.0.1.tar.gz
```

```r
install.packages(c("BSDA", "hypergeo", "gsl"), dependencies = TRUE)
install.packages("./BFF_3.0.1.tar.gz", repos = NULL, type = "source")
```

---

### 4. Install Archived Version of `MXM` (v1.5.5)

```bash
wget https://cran.r-project.org/src/contrib/Archive/MXM/MXM_1.5.5.tar.gz
```

```r
mxm_packages <- c(
  "lme4", "doParallel", "relations", "Rfast", "visNetwork", 
  "energy", "geepack", "bigmemory", "coxme", "Rfast2", "Hmisc"
)
mxm_packages <- mxm_packages[!(mxm_packages %in% installed.packages()[,"Package"])]

if (length(mxm_packages)) {
  install.packages(mxm_packages, dependencies = TRUE, repos = "http://cran.us.r-project.org")
}

install.packages("./MXM_1.5.5.tar.gz", repos = NULL, type = "source")
```

---

### 5. Install FCI.Utils

The FCI.Utils R package is available at <https://github.com/adele/FCI.Utils>/

You can install the development version directly from GitHub:

``` r
install.packages("devtools", dependencies=TRUE)
devtools::install_github("adele/FCI.Utils", dependencies=TRUE)
```

---

## 🐧 Optional: System Dependencies for Debian/Ubuntu Linux

The following system libraries may be required for full functionality (e.g., for graphics rendering, parallel processing, or Rcpp-based packages). You can install them using:

```bash
sudo apt-get update
sudo apt-get install -y \
  libmagick++-dev \
  cargo \
  librsvg2-dev \
  libavfilter-dev \
  libharfbuzz-dev \
  libcurl4-openssl-dev \
  libudunits2-dev \
  cmake \
  libsodium-dev \
  libssl-dev \
  libxml2-dev \
  libgdal-dev \
  libfontconfig1-dev \
  libcairo2-dev
```

> ⚠️ These packages are typically needed for rendering graphs, working with web graphics, compiling C++ code, or handling advanced statistical routines.

---

## ✅ Final Notes

* Ensure that your R version is up to date (≥ 4.2 recommended).
* If you're running on Windows or macOS, equivalent system tools may be needed (e.g., Rtools, Xcode).
* For questions or bug reports, please open an issue on the [dcFCI GitHub repository](#).


