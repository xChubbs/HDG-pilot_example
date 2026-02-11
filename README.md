# HDG: pilot example

This project provides a computational framework for exploring and solving a pilot example of a second-order curl problem using the powerful Hybridizable Discontinuous Galerkin (HDG) finite element method. Whether you're a seasoned researcher or a student in numerical analysis, this repository aims to offer a clear, practical example implemented in MATLAB.

## Project Overview

This repository addresses a canonical second-order curl problem, a fundamental challenge in computational electromagnetics, fluid dynamics, and other fields involving vector fields. The problem is formulated as follows: Find $\boldsymbol{u}$ in $\Omega$ such that:

```math
\begin{align*}
    \boldsymbol{u} + \textbf{curl}(\textbf{curl}(\boldsymbol{u})) &= \boldsymbol{f} \quad &\textup{in}~\Omega \\
    \boldsymbol{u} \times \boldsymbol{n} &= \boldsymbol{g} \times \boldsymbol{n} \quad &\textup{on}~\partial\Omega
\end{align*}
```

Our implementation is primarily in MATLAB, designed to be accessible for researchers and practitioners in numerical analysis and scientific computing.

## Key Features

This repository offers a comprehensive look into an HDG implementation, highlighting several crucial components:

*   **HDG Scheme Implementation**: Dive into the core MATLAB source files (`./src`), where you'll find the building blocks of the HDG scheme. This includes the construction of element and face matrices, the systematic imposition of boundary conditions, and robust error estimation techniques.
*   **Second-Order Curl Problem**: We focus on a challenging "second-order curl" problem, which is vital for understanding vector calculus applications. This class of problems is known for its computational complexities, particularly in ensuring divergence-free conditions, and is highly relevant to areas like Maxwell's equations.
*   **Quadrature and Basis Functions**: The `./include/sayas` directory provides specialized tools for numerical integration and polynomial basis functions. Here, you'll discover routines like `TablesQuadForm.m`, `TablesQuadForm3d.m`, `dubiner2d.m`, and `dubiner3d.m`, which are essential for computing high-order quadrature rules and evaluating Dubiner polynomials—common choices for shape functions on triangular or tetrahedral elements.
*   **Matrix Assembly**: Explore files like `./src/matricesFace.m`, which are central to the efficient assembly of the global system matrix within the HDG formulation. These files meticulously handle terms related to face integrals and penalty parameters.
*   **Boundary Conditions and Error Analysis**: The repository includes functionalities for applying 3D boundary conditions (`./src/BC3D.m`) and for performing both element-wise (`./src/errorElem.m`) and face-wise (`./src/errorFace.m`) error analysis. These are critical for verifying the accuracy and convergence of the numerical scheme.

## Getting Started

To explore the HDG scheme in action, we recommend starting with the `TB-pilot_example.mlx` file. This MATLAB Live Script serves as the main entry point, guiding you through the setup, execution, and visualization of the HDG solution for the pilot example.

## Repository Tree

For a quick overview of the project structure, here's the repository tree:

```
HDG/
├── .gitignore
├── LICENSE
├── README.md
├── TB-pilot_example.mlx
├── Tex/
├── include/
│   ├── gmsh/
│   └── sayas/
│       ├── TablesQuadForm.m
│       ├── TablesQuadForm3d.m
│       ├── dubiner2d.m
│       ├── dubiner3d.m
│       ├── testElem.m
│       └── testFaces.m
└── src/
    ├── BC3D.m
    ├── convMatrix.m
    ├── errorElem.m
    ├── errorFace.m
    ├── matricesFace.asv
    └── matricesFace.m
