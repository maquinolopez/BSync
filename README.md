# BSync

*Bayesian Synchronisation for Proxy Time‑Series*

![BSync logo](logo.png)

BSync is an open‑source R implementation of the **Bayesian Synchronisation** framework described in our forthcoming paper *Bayesian Synchronization of Proxy Paleorecords with Reference Chronologies*.  The method aligns an undated (or imprecisely dated) proxy record to one or more well‑dated reference records **while fully propagating alignment uncertainty**.  BSync supports single‑target, double‑target and target‑with‑age‑uncertainty workflows.

---

## Table of contents

1. [Overview](#overview)
2. [Key features](#key-features)
3. [Installation](#installation)
4. [Quick start](#quick-start)
5. [Function reference](#function-reference)
6. [Reproducing the paper figures](#reproducing-the-paper-figures)
7. [Citing BSync](#citing-bsync)
8. [Contributing](#contributing)
9. [Licence](#licence)
10. [Contact](#contact)

---

## Overview

Synchronising palaeo‑environmental records is essential for multi‑archive climate reconstructions, yet traditional "tie‑point" approaches provide little or no statistical measure of the alignment uncertainty.  BSync tackles this by casting synchronisation as a Bayesian inverse problem:

* A flexible **link function** (inverse‑gamma process) deforms the input age/depth scale.
* Prior information on sedimentation behaviour and memory is encoded through interpretable hyper‑parameters.
* The likelihood couples the rescaled input and target proxy signals using a heavy‑tailed *t* distribution or a kernel density estimator when the target carries its own age uncertainty.
* The full posterior is sampled with the **t‑walk** MCMC algorithm, enabling credible intervals for the alignment, accumulation rates and derived chronologies.

> **Status of the accompanying manuscript**
> A detailed description of the model, benchmarking and application case‑studies is **in preparation and will be submitted for peer review in late‑2025**.  A pre‑print will be posted on ArXiv as soon as it is ready.  Early adopters are encouraged to cite this repository and to contact the author for an advance copy.

---

## Key features

* **Three alignment modes** – single target, double target (optimal mixture of two chronology‑consistent targets) and target with reported age uncertainty.
* **Objective yet configurable** – default priors work out‑of‑the‑box, but every hyper‑parameter is exposed for expert tuning.
* **Heavy‑tailed error model** to reduce the impact of local noise and proxy‑specific variance.
* **Automatic convergence monitoring** through log‑objective traces and ESS diagnostics.
* **Publication‑ready plots** of proxy alignment, age‑depth models and parameter posteriors.

---

## Installation

BSync is written in R (>= 4.2) with small C++ extensions via **Rcpp**.  You will need the following libraries:

```r
install.packages(c("Rcpp", "KernSmooth", "coda"))
```

Clone the repository and source the main script:

```bash
git clone https://github.com/maquinolopez/BSync.git
cd BSync
```

---

## Quick start

```r
# Load the core routine
source("Bsynchv4.R")  # script will load the helpers automatically

# files shall be in a folder name "input target"

# Run Bayesian synchronisation (defaults are usually fine)
Bsynch(input,           # name of the input record
       target,          # name of the target record
       burn      = 5e4,
       iterations = 3e5,
       thinning   = 100)

```

See `examples/` for fully‑worked scripts that reproduce selected figures from the manuscript.

---

## Function reference

The main entry point is `Bsynch()`

| Argument                   | Type         | Default       | Description                                    |
| -------------------------- | ------------ | ------------- | ---------------------------------------------- |
| `input_record`             | `data.frame` | required      | Proxy record to be aligned                     |
| `target_record`            | `data.frame` | required      | Reference proxy record                         |
| `folder`                   | `character`  | `getwd()`     | Output directory for MCMC chains, plots & logs |
| `burn`                     | `integer`    | `5e4`         | Burn‑in iterations                             |
| `iterations`               | `integer`    | `3e5`         | Total iterations                               |
| `thinning`                 | `integer`    | `100`         | Keep every *n*‑th sample                       |
| `sigma.prior`              | `numeric`    | `1`           | Prior SD of proxy residuals                    |
| `n_sections`               | `integer`    | `50`          | Number of inverse‑gamma sections               |
| `shape_acc`, `mean_acc`    | `numeric`    | `1.5`, `50`   | Hyper‑parameters of accumulation‑rate prior    |
| `strength_mem`, `mean_mem` | `numeric`    | `10`, `0.5`   | Memory prior parameters                        |
| `sd_shape`, `sd_scale`     | `numeric`    | `1.5`, `0.01` | Hyper‑parameters of the SD prior               |
| `savefig`                  | `logical`    | `TRUE`        | Save alignment plot as PDF                     |

Full documentation is available in the function header.

---

## Reproducing the paper figures

Scripts in `paper/` replicate every figure and table from the upcoming manuscript, including:

* Iberian Margin case‑study (benthic & planktic δ¹⁸O).
* Synthetic‐core benchmarks for single‑, double‑ and uncertainty‑aware targets.
* Comparison with **BIGMACS**.

Running `make paper` will execute all workflows and save outputs in `results/` (≈ 4 h on a modern laptop).

---

## Citing BSync

If you use BSync in your research, please cite **both** the software and the manuscript (once available):

```bibtex
@software{aquino2025_bsync,
  author  = {Marco A. Aquino‑López},
  title   = {BSync: Bayesian Synchronisation of Proxy Palaeo‑records},
  year    = {2025},
  url     = {https://github.com/maquinolopez/BSync}
}

@unpublished{aquino2025_manuscript,
  author = {Aquino‑López, M. A. and Osman, M. and Muschitiello, F.},
  title  = {Bayesian Synchronization of Proxy Paleorecords with Reference Chronologies},
  note   = {Manuscript in preparation, planned submission 2025}
}
```

---

## Contributing

Pull requests are welcome!  Please open an issue first to discuss major changes.  By contributing you agree to license your work under the terms below.

---

## Licence

BSync is released under the **MIT License** – see `LICENSE` for details.

---

## Contact

Questions, suggestions or requests for an early copy of the manuscript?  Drop me a line:

**Marco A. Aquino‑López**
Centro de Investigación en Matemáticas (CIMAT), Mexico
✉️ [aquino@cimat.mx](mailto:aquino@cimat.mx)

---

*Happy synchronising!*

