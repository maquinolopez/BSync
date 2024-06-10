
# BSynch

![BSynch Logo](logo.png)

BSynch is an innovative tool designed to provide automatic and continuous synchronization between proxy records. With an intuitive setup and robust architecture, BSynch makes the task of maintaining up-to-date proxy records a breeze.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Bsynch Function](#bsynch-function)
- [Support & Feedback](#support--feedback)

## Features

- **Continuous Synchronization**: Continuous and objective synchronization between records.
- **Performance**: Automatic monitoring of convergence and thinning.
- **Bayesian Approach**: The linking function is inferred using Bayesian Inference, which allows for a robust estimate of the uncertainties of the alignment.

## Installation

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/your-username/BSynch.git
    cd BSynch
    ```

2. **Install Dependencies**:

   Ensure you have the necessary dependencies installed. The following R packages are required:

    - **Rcpp**: This package provides a seamless integration of R and C++ and is crucial for the performance of BSynch.
    - **KernSmooth**: This package provides functions for kernel smoothing.
    - **coda**: This package is used for output analysis and diagnostics of MCMC.

    You can install them using the following command in your R console:

    ```r
    # Install the required R packages
    install.packages(c("Rcpp", "KernSmooth", "coda"))
    ```

## Usage

1. Add your input and target records to a folder inside the BSynch folder.

2. Run the synchronization process in R:

    ```r
    source('Bsynchv4.R')

    # Define your input and target records
    input_record <- read.csv('path/to/input_record.csv')
    target_record <- read.csv('path/to/target_record.csv')

    # Perform synchronization
    synchronized_record <- Bsynch(input_record, target_record)

    # View the synchronized record
    print(synchronized_record)
    ```

## Bsynch Function

The `Bsynch` function is the core of the BSynch tool. It synchronizes the input record with the target record using a sophisticated Bayesian inference method. Here’s how you can utilize the `Bsynch` function:

### Function Definition

```r
Bsynch <- function(input_record, target_record, ...) {
    # Function implementation
}
```

### Parameters

- `input_record`: **(Required)** The name of the CSV or TXT file containing the input data that needs to be synchronized.
- `target_record`: **(Required)** The name of the CSV or TXT file containing the input data that needs to be synchronized.
- `folder`: **(required)** The folder where the input_record and target_record are located. By default `~/Documents/BSync/`
- `burn`: **(Optional)** The number of initial iterations to discard (burn-in period). Default is `5e+4`. This helps in allowing the MCMC chain to stabilize before collecting samples.
- `iterations`: **(Optional)** The total number of iterations for the MCMC. Default is `3e+3`. This determines how many samples the MCMC will generate.
- `thinning`: **(Optional)** The thinning interval for the MCMC chain. Default is `100`. This reduces the autocorrelation by only keeping every nth sample.
- `sigma.prior`: **(Optional)** The prior distribution for the standard deviation. Default is `1`. This influences the smoothness of the alignment.
- `savefig`: **(Optional)** A boolean indicating whether to save the alignment plot as a PDF. Default is `TRUE`.
- `folder`: **(Optional)** The directory where the output files and plots will be saved. Default is the current working directory.
- `Input`: **(Optional)** The name of the input record file.
- `Target`: **(Optional)** The name of the target record file.
- `n_sections`: **(Optional)** The number of sections to divide the data into for processing. Default is `50`. This helps in managing large datasets by processing them in smaller chunks.
- `n_sections`: **(Optional)** The number of sections to divide the data into for processing. Default is `50`. 
- - `shape_acc`: **(Optional)** The shape parameter for the accumulation rate of the model. Default is `1.5`.
- `mean_acc`: **(Optional)** The mean parameter for the accumulation rate of the model. Default is `50`.
- `strength_mem`: **(Optional)** The strength parameter for the memory effect in the model. Default is `10`.
- `mean_mem`: **(Optional)** The mean parameter for the memory effect in the model. Default is `0.5`.
- `sd_shape`: **(Optional)** The shape parameter for the prior distribution of the variance in the alignment model. Default is `1.5`.
- `sd_scale`: **(Optional)** The scale parameter for the prior distribution of the variance in the alignment model. Default is `0.01`.

- `savefig`: **(Optional)** A boolean indicating whether to save the alignment plot as a PDF. Default is `TRUE`.

### Usage Example

Here’s a basic example of how to use the `Bsynch` function in your R code:

```r
# Load the BSynch library
source('Bsynchv4.R')

# Define your input and target records
input_record <- read.csv('path/to/input_record.csv')
target_record <- read.csv('path/to/target_record.csv')

# Perform synchronization
synchronized_record <- Bsynch(input_record, target_record, burn=500, iterations=5000, thinning=5, sigma.prior=0.5, savefig=FALSE, folder='output')

# View the synchronized record
print(synchronized_record)
```


## Support & Feedback

For any questions, feedback, or concerns, please raise an issue on the GitHub repository or contact the maintainers directly.

---

Thank you for choosing BSynch for your synchronization needs!
