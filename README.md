# Scalable Variational Bayes for Spatial Generalized Linear Mixed Models

This repository contains R scripts accompanying the paper:

**“A Scalable Variational Bayes Approach to Fit High-Dimensional Spatial Generalized Linear Mixed Models.”**

The code includes simulation studies, implementations of the proposed variational inference methods, and comparisons with alternative inference methods.

## Repository Structure

The `01SimulationStudy` directory contains two main subdirectories:

- **`01FullSGLMM`**: Implementations of the full spatial generalized linear mixed model (SGLMM).
- **`02BasisSGLMM`**: Implementations of SGLMMs using spatial basis functions.

The `01FullSGLMM` directory contains simulation studies for three response types:

- **`01Binary`**: Binary responses.
- **`02Count`**: Count responses.
- **`03Normal_Supplementary`**: Gaussian responses for supplementary analyses.

## Example: Binary Response Simulation

The `01FullSGLMM/01Binary` directory contains the following subdirectories:

| Directory | Description |
| --- | --- |
| `01DataGeneration` | Generates binary responses at 500 spatial locations. |
| `02MCMC` | Implements Metropolis–Hastings MCMC for the full SGLMM. |
| `03INFVBfixphi_ParallelComputing` | Implements INFVB with the spatial range parameter, φ, held fixed within each conditional fit, using parallel computing. |
| `04INFVBfixssqphi_ParallelComputing` | Implements INFVB with both the spatial range parameter, φ, and the partial sill, σ², held fixed within each conditional fit, using parallel computing. |
| `05Summary_Figures` | Produces figures comparing the methods. |
| `06INLA` | Implements integrated nested Laplace approximation (INLA) for comparison. |
| `07Measures` | Summarizes evaluation results, including the continuous ranked probability score (CRPS) and 95% interval assessments. |
| `08INFVBfixphi_SA` | Contains INFVB code for 95% interval calculations, with φ held fixed within each conditional fit. |
| `09INFVBfixssqphi_SA` | Contains INFVB code for 95% interval calculations, with φ and σ² held fixed within each conditional fit. |
| `10Stan` | Implements Hamiltonian Monte Carlo (HMC) using Stan for comparison. |

The basis-function implementations and the simulation studies for count and Gaussian responses follow a similar organizational structure.

## Contact

For questions about the code, please contact Jin Hyung Lee at [lee6024@purdue.edu](mailto:lee6024@purdue.edu).
