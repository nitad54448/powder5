# Powder Pattern Decomposition - Help (powder5)

This document serves as a scientific and technical guide for the **powder5** web application, a tool for the analysis of powder X-ray diffraction (PXRD) data.

## Help Topics

* [Introduction](#introduction)
* [Getting Started: Data Input](#getting-started-data-input)
* [Interactive Data Visualization](#interactive-data-visualization)
* [Pattern Decomposition Methodologies](#pattern-decomposition-methodologies)
    * [The Le Bail Method](#the-le-bail-method)
    * [The Pawley Method](#the-pawley-method)
* [Minimization Algorithms](#minimization-algorithms)
    * [Levenberg-Marquardt (LM)](#levenberg-marquardt-lm)
    * [Simulated Annealing (SA)](#simulated-annealing-sa)
    * [Parallel Tempering (Replica Exchange)](#parallel-tempering-replica-exchange)
* [Guide to Refinable Parameters](#guide-to-refinable-parameters)
    * [Crystal System & Space Group](#crystal-system--space-group)
    * [Instrumental Parameters](#instrumental-parameters)
    * [Background Modeling](#background-modeling)
    * [Profile Function #4 (Simple pseudo-Voigt)](#profile-function-4-simple-pseudo-voigt)
    * [Profile Function #3 (Anisotropic TCH)](#profile-function-3-anisotropic-tch)
* [Recommended Refinement Strategy](#recommended-refinement-strategy)
* [Interpretation of Results](#interpretation-of-results)
    * [Williamson-Hall Size-Strain Analysis](#williamson-hall-size-strain-analysis)
* [References & Further Reading](#references--further-reading)
* [About This Tool](#about-this-tool)

---

## Introduction

This document serves as a scientific and technical guide for the **powder5** web application, a tool for the analysis of powder X-ray diffraction (PXRD) data via whole-pattern fitting. This technique, also known as pattern decomposition, is a crucial method in materials science and crystallography for refining structural and microstructural parameters when a complete structural model is either unknown or unnecessary.

The application facilitates the extraction of precise lattice parameters, peak profile information, and integrated intensities of Bragg reflections. It implements two principal decomposition algorithms: the iterative **Le Bail method** for rapid and stable convergence, and the simultaneous **Pawley method** for rigorous, unbiased intensity extraction. Peak profiles are modeled using phenomenological functions adapted from the GSAS software package, including a versatile **Simple pseudo-Voigt (Profile Function #4)** and a more physically rigorous **Anisotropic model (Profile Function #3)** based on the Thompson-Cox-Hastings (TCH) formulation with a Stephens model for anisotropic line broadening.

---

## Getting Started: Data Input

Analysis commences with the loading of a powder diffraction data file. The application is designed to automatically parse numerous common ASCII-based file formats from major instrument manufacturers and standard crystallographic software.

* **Supported Formats:** Built-in parsers are included for Bruker (`.brml`, `.uxd`), PANalytical (`.xrdml`), Rigaku (`.ras`), Philips (`.udf`), and GSAS (`.esd`, `.gsa`, `.xra`).
* **Generic Data:** Standard two-column ASCII files containing $2\theta$ and intensity values are also supported. The parser accommodates space, comma, or semicolon delimiters. Comment lines prefixed with `#` are ignored.
* **Metadata Parsing:** For many instrument-specific formats, instrument parameters such as the X-ray wavelength for Kα1 are read from the file's metadata and used to populate the relevant fields in the user interface. It is incumbent upon the user to verify the correctness of these automatically populated values.

---

## Interactive Data Visualization

The diffraction pattern is rendered in an interactive plot to facilitate detailed inspection of the experimental data and the quality of the model fit. You can hide any of the plots by clicking on its legend at the top of the chart.

* **Pan:** Click and drag on the chart to translate the view along the $2\theta$ and intensity axes.
* **Zoom:** Use the mouse wheel to adjust the zoom level. The zoom behavior is context-dependent:
    * **Over the plot area:** Zooms both axes isotropically.
    * **Over the Intensity (Y) axis:** Zooms the vertical axis exclusively.
    * **Over the $2\theta$ (X) axis:** Zooms the horizontal axis exclusively.
* **Reset View:** A **right-click** on the chart resets the zoom and pan to the range currently defined by the **2θ Min/Max sliders**.
* **Reflection Data:** Hovering the cursor near a Bragg reflection marker (tick mark) displays a tooltip containing the corresponding Miller indices ($hkl$).
* **Background Anchor Points:** To impose constraints on the background model, hold down the **Ctrl** key and **click** on a point in the pattern. This designates the nearest experimental data point as an anchor, a point that is given a significantly higher weight during refinement, thereby compelling the background function to pass through its vicinity.

---

## Pattern Decomposition Methodologies

Pattern decomposition enables the fitting of a powder diffraction pattern based on a unit cell and space group, without requiring a full structural model (atomic coordinates). This is essential for the precise determination of lattice parameters and the extraction of integrated intensities, which are requisite for ab initio structure determination.

> **A Note on Intensity:** In pattern decomposition, the "intensity" parameter ($I_{hkl}$) of a Bragg reflection refers to its total **integrated area**, not its maximum peak height. The program correctly normalizes the calculated peak shape function by its mathematical area, ensuring that the refined $I_{hkl}$ parameter is a physically meaningful quantity representing the total scattering power of a reflection.

### The Le Bail Method

The Le Bail method is an iterative, sequential algorithm known for its computational efficiency and robust convergence. The process is as follows:

1.  **Initialization:** A theoretical pattern is calculated from the user-supplied lattice, profile, and background parameters. The integrated intensities of all reflections ($I_{hkl}$) are initially assumed to be equal.
2.  **Intensity Extraction:** The observed intensity at each data point ($y_{i,obs}$) is partitioned among the calculated Bragg peaks that contribute to that point. The contribution of each peak is proportional to its profile function value at that point. Summing these partitioned intensities for each reflection yields a new set of "observed" integrated intensities.
3.  **Parameter Refinement:** These newly extracted intensities are held constant, and a non-linear least-squares refinement is performed on the lattice, background, and profile parameters to minimize the difference between the calculated and observed patterns.
4.  **Iteration:** The refined parameters from Step 3 are used to restart the process from Step 1. This cycle repeats until the parameters and R-factors converge.

### The Pawley Method

The Pawley method employs a simultaneous, non-iterative approach. It treats the integrated intensity ($I_{hkl}$) of each Bragg reflection as an independent, refinable variable within a single, large-scale least-squares minimization.

This means that all parameters—lattice, profile, background, and all individual integrated intensities—are refined concurrently.

* **Advantages:** The Pawley method is considered more rigorous as it avoids the iterative bias of the Le Bail method, particularly in cases of severe peak overlap. It can yield more accurate and statistically sound integrated intensities and parameter uncertainties.
* **Disadvantages:** The inclusion of hundreds or thousands of intensity variables significantly increases the computational complexity and memory requirements. The refinement is also more susceptible to instability and parameter correlation, especially if the initial model is poor.

---

## Minimization Algorithms

The goal of refinement is to minimize the sum-of-squares objective function, $\chi^2 = \sum w_i (y_{i,obs} - y_{i,calc})^2$, where $w_i$ is the statistical weight of each data point. This application provides several algorithms to navigate the complex parameter space and find the minimum of this function.

### Levenberg-Marquardt (LM)

The LM algorithm is a standard gradient-based method for non-linear least-squares problems. It effectively interpolates between the Gauss-Newton algorithm and the method of gradient descent. By calculating the Jacobian matrix (the matrix of first partial derivatives of the calculated pattern with respect to the parameters), it determines the most efficient path toward the nearest local minimum.

* **Characteristics:** LM is a local minimizer, exhibiting rapid quadratic convergence when the initial parameters are close to the true minimum. It is the preferred method for final, high-precision refinement and is the only algorithm here that can calculate valid estimated standard deviations (ESDs) for the refined parameters from the covariance matrix.
* **Limitations:** It is susceptible to converging to a local minimum if the starting model is far from the global solution.

### Simulated Annealing (SA)

Simulated Annealing is a stochastic, global optimization metaheuristic inspired by the annealing process in metallurgy. The algorithm introduces a "temperature" parameter, $T$, which controls the probability of accepting a parameter step that increases the objective function.

At high temperatures, the algorithm explores the parameter space widely, frequently accepting "uphill" moves that allow it to escape local minima. As the refinement proceeds, $T$ is gradually reduced according to a cooling schedule, reducing the probability of accepting non-improving solutions. The system eventually "freezes" into a low-energy state, which is expected to be the global minimum.

* **Use Case:** Ideal for exploring complex solution spaces or when the initial parameters are poorly known. It can effectively escape local minima where gradient-based methods like LM would become trapped.

### Parallel Tempering (Replica Exchange)

Parallel Tempering, also known as Replica Exchange MCMC, is an advanced stochastic algorithm designed to overcome the slow convergence of traditional Simulated Annealing. Instead of a single system with a decreasing temperature, Parallel Tempering simulates multiple copies (or "replicas") of the system simultaneously, each at a different, fixed temperature in a predefined ladder ($T_1 < T_2 < ... < T_N$).

* **Mechanism:** Each replica evolves independently according to a standard Monte Carlo or SA algorithm at its respective temperature. The high-temperature replicas explore the parameter space broadly (high mobility), while the low-temperature replicas perform a fine-grained search of local minima (high precision).
* **The Swap Move:** Periodically, the algorithm attempts to swap the entire set of parameters between adjacent replicas (e.g., between replica $i$ at temperature $T_i$ and replica $i+1$ at $T_{i+1}$). The swap is accepted with a Metropolis-like probability that depends on the energies and temperatures of the two replicas. This crucial step allows a good solution discovered by a high-temperature replica in a distant valley to "percolate down" to the low-temperature replicas, dramatically improving the efficiency of finding the global minimum.
* **Advantages:** Significantly more efficient at global exploration than standard Simulated Annealing, making it the most robust choice for complex problems or when the starting model is highly uncertain.

---

## Guide to Refinable Parameters

This section provides a detailed breakdown of the parameters you can control and refine.

> #### A Note on Parameter Scaling & GSAS Comparison
>
> Following a long-standing convention in crystallographic software like GSAS, some refinable parameters in this program are internally scaled. This is done for user convenience, allowing you to work with manageable numbers (e.g., 1.0) instead of very small decimals (e.g., 1.0e-4). The documentation below provides the exact formulas used, allowing for direct comparison with physical models and values from other software.

### Crystal System & Space Group

These parameters define the crystallographic symmetry of the material.

* The **System** selection imposes metrical constraints on the lattice parameters (e.g., for Cubic, $a=b=c$, $\alpha=\beta=\gamma=90^\circ$).
* The **Space Group** selection determines the systematic reflection conditions ($hkl$ absences) used to generate the list of Bragg peaks. The underlying logic for these conditions is consistent with established crystallographic libraries and was taken from Computational Crystallography Toolbox (cctbx).

### Instrumental Parameters

Found under the "Sample" tab, these parameters model the diffractometer configuration.

* `Radiation 1/2 (Å) & Ratio`: Defines the X-ray source. For divergent-beam laboratory instruments, a Kα1/Kα2 doublet is typically used. For synchrotron radiation, the Ratio is set to 0.
* `Zero`: A refinable parameter that corrects for instrumental zero-point error in the $2\theta$ axis. It is highly correlated with lattice parameters and must be refined with caution.
* `2θ Min / Max`: These sliders define the refinement range. It is standard practice to exclude regions of low signal-to-noise or known artifacts from the calculation.

### Background Modeling

The background is modeled as a superposition of empirical functions on the "Background" tab.

* **Chebyshev Polynomial:** A 6-term Chebyshev polynomial of the first kind provides a flexible and orthogonal basis for modeling slowly varying background features. It is recommended to refine coefficients sequentially, starting with the lowest-order terms.
* **Amorphous Hump:** A single Lorentzian function can be added to the background model to account for broad scattering features, such as from a glass capillary or an amorphous sample component.
* **Anchor Points:** Manual anchor points can be added by **Ctrl+Clicking** on the chart. These points are assigned a high weight in the least-squares minimization, effectively constraining the background function to pass through specified regions.

### Profile Function #4 (Simple pseudo-Voigt)

This function models the peak shape as a linear combination of a Gaussian and a Lorentzian function: $pV(x) = \eta L(x) + (1-\eta)G(x)$. The angular dependence of the Full Width at Half Maximum (FWHM) for each component is modeled empirically.
$$H_G^2 = GU \tan^2\theta + GV \tan\theta + GW + GP / \cos^2\theta$$
$$H_L = LX / \cos\theta$$

* `GU, GV, GW, GP`: Parameters describing the Gaussian FWHM. The terms are associated with strain ($GU$), instrumental factors ($GV, GW$), and particle size effects ($GP$).
* `LX`: Describes the Lorentzian FWHM, primarily associated with crystallite size broadening.
* `eta`: A simple linear mixing parameter ($\eta=0$ for pure Gaussian, $\eta=1$ for pure Lorentzian).
* `shft & trns`: Corrects for peak position shifts due to sample displacement and transparency, respectively.
    > **Unit and Scaling for the `shft` parameter:**
    >
    > The refined `shft` parameter is a **dimensionless, scaled coefficient**, not a direct physical length. Its relationship to the physical specimen displacement ($s$) and the goniometer radius ($R$) is defined as follows:
    >
    > * The physical peak shift in radians is: $\Delta(2\theta)_{\text{rad}} = -2 \frac{s}{R} \cos(\theta)$
    > * The program calculates this shift (in degrees) using the formula: $\Delta(2\theta)_{\text{deg}} = -(\text{shft} / 1000) \times \cos(\theta) \times (180 / \pi)$
    > * Therefore, the relationship is: $\frac{\text{shft}}{1000} = \frac{2s}{R}$
    > * To find the physical displacement $s$ from the refined parameter, use: $s = R \times (\text{shft} / 2000)$.
    >
    > *Example: For a typical instrument with $R=240$ mm, a refined `shft` value of 1.0 corresponds to a physical displacement $s$ of $240 \times (1/2000) = 0.12$ mm.*

### Profile Function #3 (Anisotropic TCH)

This is a physically rigorous profile function adapted from GSAS, convoluting a **Thompson-Cox-Hastings (TCH) pseudo-Voigt** with a **Stephens model** for anisotropic strain broadening.

#### Isotropic Broadening (TCH Model)

The TCH formulation models the FWHM of the Gaussian ($H_G$) and Lorentzian ($H_L$) components based on physical contributions to line broadening:
$$H_G^2 = U \tan^2\theta + V \tan\theta + W$$
$$H_L = X \tan\theta + Y / \cos\theta$$

The total FWHM ($H$) and pseudo-Voigt mixing parameter ($\eta$) are then derived from these components using polynomial approximations. The final shape is $pV(x) = \eta L(x, H) + (1-\eta)G(x, H)$, where both functions share the same convoluted FWHM.

* `U, V, W`: Gaussian broadening parameters related to strain ($U, V$) and instrumental resolution ($W$), following the Caglioti formulation.
* `X, Y`: Lorentzian broadening parameters related to strain ($X$) and crystallite size ($Y$), following the Scherrer formula.

#### Peak Asymmetry

The `S/L` and `H/L` parameters introduce an angle-dependent asymmetry, primarily correcting for axial divergence effects at low $2\theta$.

#### Anisotropic Broadening (Stephens Model)

Anisotropic microstrain, where broadening varies with crystallographic direction, is modeled by adding terms to the Lorentzian component ($H_L$) that are dependent on the Miller indices ($hkl$). The model is a fourth-order polynomial in the reciprocal lattice vectors.

The refinable parameters (`S400`, `S040`, etc.) are the non-zero, symmetry-unique coefficients of this polynomial. The application automatically applies symmetry constraints based on the Laue class of the selected space group.

> **Unit and Scaling for Stephens `S_hkl` parameters:**
>
> The user-inputted `S_hkl` parameters are scaled for convenience. The dimensionless anisotropic broadening term ($H_{aniso}$) is calculated from these parameters, and its contribution to the total Lorentzian width (in degrees $2\theta$) is scaled by a factor of 1000.
> $$H_L(\text{total}) = H_L(\text{isotropic}) + \frac{|H_{aniso}|}{1000}$$
> This scaling allows the user to refine values in a manageable range (e.g., -10 to +10) rather than requiring input of very small decimals (e.g., 1e-4), a convention common in other refinement software.

---

## Recommended Refinement Strategy

A sequential and hierarchical refinement strategy is crucial for achieving a stable and physically meaningful solution. Attempting to refine all parameters simultaneously from a poor starting model will likely lead to divergence or convergence to a false minimum.

### Phase 1: Initial Model Setup

1.  **Define the Model:** Load data, select the crystal system and space group, and define the refinement range using the **$2\theta$ sliders**.
2.  **Background Refinement:** Model the background by refining the first two or three Chebyshev terms (`B0`, `B1`, `B2`). Use anchor points if necessary to constrain complex background shapes.
3.  **Peak Position Refinement:** Using the **Le Bail** and **LM** algorithms, refine only the **Lattice Parameter(s)** and, if necessary, the instrumental **Zero Shift**. The goal is to align the calculated Bragg positions with the observed peak maxima.

### Phase 2: Peak Profile Refinement

4.  **Isotropic Broadening:** Once positions are correct, refine the primary isotropic peak shape parameters (e.g., **W, Y, U, X** in Profile #3). This will account for the dominant size and strain contributions.
5.  **Asymmetry and Shape:** Introduce asymmetry parameters (**S/L**) if there is a clear misfit at low angles.
6.  **Anisotropic Broadening:** If systematic misfits remain (e.g., some peaks are consistently broader than the model), introduce the anisotropic Stephens parameters (e.g., **S400**). Refine only the symmetry-independent terms.

### Phase 3: Finalization and Intensity Extraction

7.  **Global Optimization:** If the LM algorithm converges to a poor solution, switch to **Parallel Tempering** or **Simulated Annealing** for one or more cycles to perform a global search. Afterwards, switch back to LM for a final, precise local minimization.
8.  **Pawley Refinement:** With a stable and well-refined model from the Le Bail method, a final refinement can be performed using the **Pawley method**. This will provide the most statistically robust set of integrated intensities, suitable for subsequent structure solution.

> #### Technical Note on Calculation Parameters
>
> The theoretical pattern is calculated using a hybrid approach for performance and accuracy. A search window (defined by `CALCULATION_WINDOW_MULTIPLIER`) is used to find all potential contributions to a given data point. However, only points where the calculated peak height is greater than a defined threshold (`PEAK_HEIGHT_CUTOFF`, now set at 0.2%) relative to its maximum are included in the final sum.

---

## Interpretation of Results

Assessing the quality of a refinement requires both statistical analysis and critical visual inspection of the fit.

### Figures of Merit

Standard crystallographic R-factors are provided to quantify the quality of the fit.

* **R-pattern ($R_p$):** The unweighted residual error, sensitive primarily to the fit of high-intensity reflections.
    $$R_p = \frac{\sum |y_{i,obs} - y_{i,calc}|}{\sum y_{i,obs}} \times 100\%$$
* **Weighted R-pattern ($R_{wp}$):** The primary figure of merit, weighted by the inverse of the observed intensities ($w_i = 1/y_{i,obs}$), which properly accounts for the counting statistics across the entire pattern.
    $$R_{wp} = \left[ \frac{\sum w_i (y_{i,obs} - y_{i,calc})^2}{\sum w_i y_{i,obs}^2} \right]^{1/2} \times 100\%$$
* **Reduced Chi-squared ($\chi^2$, Goodness of Fit):** The most statistically rigorous indicator. For a statistically perfect fit where the model correctly describes the data and the weights are accurate, $\chi^2$ should be exactly 1.0.
    $$\chi^2 = \frac{\sum w_i (y_{i,obs} - y_{i,calc})^2}{N - P}$$
    where $N$ is the number of data points and $P$ is the number of refined parameters.

#### Calculating Observed Intensities ($I_{obs}$) for Overlapping Peaks

A simple numerical integration over a fixed angular range is insufficient for accurately determining the observed integrated intensity ($I_{obs}$) of overlapping peaks. This tool employs a more robust **intensity partitioning** method.

At each point in the diffraction pattern, the net observed intensity ($y_{obs} - y_{bkg}$) is distributed among all contributing Bragg reflections. This distribution is proportional to the value of each peak's calculated profile function at that specific point. By integrating these partitioned "slices" of intensity for each reflection across the entire pattern, the method yields a reliable $I_{obs}$ value that correctly deconvolutes contributions from neighboring peaks.

### Visual Inspection

Numerical indicators can be misleading. Visual inspection of the difference plot (observed minus calculated) is the most critical step in evaluating the fit.

* **The Difference Plot:** A successful refinement should yield a difference plot that consists of random, uncorrelated noise centered on zero.
* **Systematic Residuals:** The presence of structured, non-random features in the difference plot (e.g., "M-shaped" residuals around peaks, broad humps, or un-indexed peaks) is a clear indication of systematic errors in the model. These may arise from an incorrect peak shape, unmodeled anisotropy or asymmetry, an inadequate background model, or the presence of an unaccounted-for impurity phase.

> ### Williamson-Hall Size-Strain Analysis
>
> For refinements utilizing **Profile Function #3 (Anisotropic TCH)** and the **Pawley method**, the application can automatically perform a Williamson-Hall analysis to extract microstructural information. This method separates the contributions of crystallite size and microstrain to the total peak broadening by analyzing their different dependencies on the diffraction angle, $\theta$.
>
> The analysis is based on the linear Williamson-Hall equation, where $\beta$ is the total physical peak breadth (FWHM) in radians:
> $$\beta \cos(\theta) = \frac{K\lambda}{L} + 4\epsilon \sin(\theta)$$
>
> This equation describes a straight line when plotting $\beta \cos(\theta)$ vs. $4\sin(\theta)$. The software performs a linear least-squares fit on the refined sample-broadening parameters (`U`, `X`, `Y`) to determine the y-intercept (related to crystallite size, $L$) and the slope (related to microstrain, $\epsilon$).
>
> #### Reported Values
>
> * **Apparent Crystallite Size (nm):** An estimate of the average size of the coherently scattering domains, calculated from the y-intercept of the Williamson-Hall plot.
> * **Apparent Microstrain (%):** An estimate of the root-mean-square strain within the crystallites, calculated from the slope of the plot.
> * **Linear Fit R²:** The coefficient of determination for the linear regression. A value close to 1.0 indicates that the isotropic size/strain model is a good fit for the observed peak broadening. Values significantly less than 1.0 may suggest that broadening is anisotropic or that the model is otherwise inadequate.

### Data Export

* **Save Report:** Generates a comprehensive ASCII text file containing all statistical indicators, refined parameter values with their ESDs (for LM refinements), and a point-by-point list of observed, calculated, and difference intensities.
* **Generate PDF:** Creates a summary PDF document containing a high-resolution plot and tables of final parameters, suitable for archival or reporting.

---

## References & Further Reading

**Pawley Method:**
Pawley, G. S. (1981). "Unit-cell refinement from powder diffraction scans". *Journal of Applied Crystallography*, 14(6), 357-361.

**Le Bail Method:**
Le Bail, A., Duroy, H. & Fourquet, J.L. (1988). "Ab-initio structure determination of LiSbWO6 by X-ray powder diffraction". *Materials Research Bulletin*, 23(3), 447-452.

**GSAS Profile Functions:**
Larson, A. C. & Von Dreele, R. B. (2004). "General Structure Analysis System (GSAS)". *Los Alamos National Laboratory Report LAUR 86-748*.

**TCH Profile Function:**
Thompson, P., Cox, D. E. & Hastings, J. B. (1987). "Rietveld refinement of Debye-Scherrer synchrotron X-ray data from Al2O3". *Journal of Applied Crystallography*, 20(2), 79-83.

**Stephens Anisotropy Model:**
Stephens, P. W. (1999). "Phenomenological model of anisotropic peak broadening in powder diffraction". *Journal of Applied Crystallography*, 32(2), 281-289.

**Parallel Tempering:**
Swendsen, R. H., & Wang, J. S. (1986). "Replica Monte Carlo simulation of spin-glasses". *Physical Review Letters*, 57(21), 2607.

**cctbx - Computational Crystallography Toolbox:**
Grosse-Kunstleve, R. W., Sauter, N. K., Moriarty, N. W., & Adams, P. D. (2002). "The Computational Crystallography Toolbox: crystallographic algorithms in a reusable software framework". *Journal of Applied Crystallography*, 35(1), 126-136.

---

## About This Tool

The **powder5** toolkit was developed by Nita Dragoe from Université Paris-Saclay as a simple browser-based implementation of powder pattern decomposition methods. It is a long-time successor of PowderV2 (Dragoe, N. (2001). *J. Appl. Cryst.*, 34, 535) and has been updated to include the Pawley method and modern global optimization algorithms.

Last Updated: 20 October 2025.
This document was updated with the assistance of an AI.

> **Disclaimer:** This application is provided for educational and research purposes. While it implements standard and robust algorithms, it is not a substitute for fully validated, peer-reviewed software packages (e.g., GSAS-II, FullProf, TOPAS) for analyses intended for publication.
