## Introduction

            
                This document serves as a scientific and technical guide for the **powder5** web application, a tool for the analysis of powder X-ray diffraction (PXRD) data via whole-pattern fitting. This technique, also known as pattern decomposition, is a crucial method in materials science and crystallography for refining structural and microstructural parameters when a complete structural model is either unknown or unnecessary.
            

            
                The application facilitates the extraction of precise lattice parameters, peak profile information, and integrated intensities of Bragg reflections. It implements two principal decomposition algorithms: the iterative **Le Bail method** for rapid and stable convergence, and the simultaneous **Pawley method** for rigorous, unbiased intensity extraction. Peak profiles are modeled using phenomenological functions, including a versatile **Simple pseudo-Voigt**, an **Asymmetric Split pseudo-Voigt**, and a physically rigorous **Anisotropic model (TCH)** based on the Thompson-Cox-Hastings formulation with a Stephens model for anisotropic line broadening. The background is modeled using a monotonic cubic spline interpolation between user-defined points.
            

        

        
---

        
            ## Getting Started: Data Input

            
                Analysis commences with the loading of a powder diffraction data file. The application is designed to automatically parse numerous common ASCII-based file formats from major instrument manufacturers and standard crystallographic software.
            

            
                * **Supported Formats:** Built-in parsers are included for Bruker (`.brml`, `.uxd`), PANalytical (`.xrdml`), Rigaku (`.ras`), Philips (`.udf`, `.rd`, `.sd`), GSAS (`.esd`, `.gsa`, `.std`, `.xra`), and Jade (`.mdi`). Format detection is by content first and file extension second, so a correctly formatted file will usually be recognised even if it has been renamed.

                * **Scientific notation** (`1.5E+03`, `1.2e4`) is accepted in the intensity or 2θ column.

                * **Ordering:** data stored in descending 2θ is detected and sorted ascending on load. Duplicated 2θ values and non-finite points are removed, with a note in the browser console.

                * **Generic Data:** Standard two-column ASCII files (`.xy`, `.csv`, `.txt`, `.dat`, `.asc` etc.) containing $2\theta$ and intensity values are also supported. The parser accommodates space, comma, or semicolon delimiters. Comment lines prefixed with `#`, `!`, `;` or `//` are ignored, as are non-numeric header lines.

                * **Metadata Parsing:** For many instrument-specific formats, instrument parameters such as the X-ray wavelength for Kα1 are read from the file's metadata and used to populate the relevant fields in the user interface. It is incumbent upon the user to verify the correctness of these automatically populated values.

            
            > **Warning: Rigaku `.rasx`:**
> this format is a ZIP archive rather than plain
                text, and is *not* currently readable. Export to `.ras` or to a
                two-column ASCII file instead. Attempting to load one gives a clear parse error
                rather than silently producing wrong data.
            > **Note: A Note on Parameter Scaling & GSAS Comparison**
> 
                    Following a long-standing convention in crystallographic software like GSAS, some refinable parameters in this program are internally scaled. This is done for user convenience, allowing you to work with manageable numbers (e.g., 1.0) instead of very small decimals (e.g., 1.0e-4). The documentation below provides the exact formulas used, allowing for direct comparison with physical models and values from other software.
                

            ### Crystal System & Space Group

            
                These parameters define the crystallographic symmetry of the material.
            

            
                * The **System** selection imposes metrical constraints on the lattice parameters (e.g., for Cubic, $a=b=c$, $\alpha=\beta=\gamma=90^\circ$).

                * The **Space Group** selection determines the systematic reflection conditions ($hkl$ absences) used to generate the list of Bragg peaks. The underlying logic for these conditions is consistent with established crystallographic libraries and was taken from Computational Crystallography Toolbox (cctbx).

                * Space groups are chosen from a searchable modal listing every setting of all 230 groups. Where a group has several settings (for example No. 62: $Pnma$, $Pmnb$, $Pbnm$, $Pcmn$, $Pmcn$, $Pnam$) each is offered separately, since the reflection conditions differ between them. Rhombohedral groups are offered on hexagonal axes only, because the $d$-spacing formula is written for that indexing.

                * **Monoclinic settings.** All three unique-axis choices are supported. The unique
                    axis is read from the symmetry operators of the setting you select, and the lattice
                    panel then offers the one angle that is free: $\beta$ for unique axis $b$
                    ($P12_1/c1$), $\gamma$ for unique axis $c$ ($P112_1/b$), $\alpha$ for unique axis
                    $a$ ($P2_1/b11$). The other two are fixed at $90^\circ$ and are neither shown nor
                    refined. The reciprocal metric, the reflection orbits and the multiplicities all
                    follow that same axis, so the three settings of one lattice give identical peak
                    positions and identical intensities. Switching between settings carries the angle
                    across, so the cell is not lost.

                * **Startup default:** the application opens on $Pnma$ (No. 62, standard setting) with an orthorhombic cell of $a = 8.478$, $b = 5.397$, $c = 6.958$ Å. This is only a starting point — change it to match your material before refining.

            

            ### Reflection List & Multiplicities

            
                The Bragg peak list is built from the *rotation operators* of the selected
                space-group setting. Every candidate $hkl$ inside the resolution sphere is mapped
                onto its orbit under the Laue group; the orbit is kept once, labelled by a canonical
                member, and its size becomes the multiplicity $m$ that scales the intensity,
                $I \propto m \cdot LP \cdot |F|^2$. Systematic absences are then tested on that
                representative, which decides the whole orbit, since an absence is a property of the
                orbit and not of any one index triple.
            

            
                Working from the operators rather than from a fixed index range per crystal system
                matters wherever two Laue classes share a system, because the two classes have
                genuinely different sets of independent reflections:
            

            
| System | Laue classes | What separates them |
| --- | --- | --- |
| Cubic | $m\bar{3}m$ / $m\bar{3}$ | $m\bar{3}$ permutes the axes cyclically only, so $(210)$ and $(201)$ are                             independent reflections. In $m\bar{3}m$ they are equivalent. |
| Tetragonal | $4/mmm$ / $4/m$ | $4/m$ has no mirror exchanging $a$ and $b$, so $(hkl)$ and $(khl)$ are                             independent. |
| Hexagonal | $6/mmm$ / $6/m$ | As above: in $6/m$, $(hkl)$ and $(khl)$ are independent. |
| Trigonal | $\bar{3}m$ / $\bar{3}$ | In $\bar{3}$ the orbit is generated by the 3-fold and inversion alone, and                             is half the size. |
| Trigonal | $\bar{3}m1$ / $\bar{3}1m$ | The two settings differ in which form is special: $\bar{3}m1$ gives                             $h0l$ multiplicity 6 and $hhl$ multiplicity 12, $\bar{3}1m$ the                             reverse. |

            
                In the hexagonal and trigonal systems the special forms are set by the
                Bravais–Miller index $i = -(h+k)$ as much as by $h$ and $k$: a reflection lies on
                a mirror whenever any two of $h$, $k$, $i$ are equal in magnitude, so
                $(1\bar{2}l)$ and $(2\bar{4}l)$ are special even though $|h| \neq |k|$. The orbit
                calculation accounts for this automatically.
            

            
                Operators are taken from the `symops` field of the space-group database
                (or the older `rotations` field, which serves equally well here — only
                the rotation parts affect a reflection orbit). The closed group is checked against the
                `order_p` recorded for the setting before it is used. If the database
                carries no operators, a built-in table of the Laue-class generators is used instead;
                if the Laue class itself cannot be identified, the generator falls back to the lowest
                Laue class of the crystal system. That last case is safe rather than wrong: orbits
                split into smaller ones at identical $d$, so peak positions and total intensities are
                unaffected, and only the number of listed entries changes. The browser console
                records which source was used whenever it is not the first.
            

            > **Note: Completeness**
> 
                    The sum of the multiplicities over the generated list equals the number of
                    reciprocal-lattice points inside the limiting sphere, exactly, for every Laue
                    class. A reflection cannot be missed and cannot be counted twice without
                    breaking that identity, which makes it a usable self-check if you modify the
                    symmetry data.
                

            ### Instrumental Parameters

            
                Found under the "Sample" tab, these parameters model the diffractometer configuration.
            

            
                * `Radiation 1/2 (Å) & Ratio`: Defines the X-ray source. For divergent-beam laboratory instruments, a Kα1/Kα2 doublet is typically used. For synchrotron radiation, the Ratio is set to 0.

                * `Polarisation`: Selects the polarisation model used in the Lorentz–polarisation factor — *Lab* (default), *Synchrotron*, or *None*. A second field appears beside it: the monochromator angle $2\theta_M$ for a laboratory source, or the polarised fraction $f$ for a synchrotron. This setting does **not** affect the fit; it governs how an extracted intensity is converted into $|F|$. See [Lorentz–Polarisation Factor](#parameters-lp).

                * `Zero`: A refinable parameter that corrects for instrumental zero-point error in the $2\theta$ axis. It is highly correlated with lattice parameters and must be refined with caution.

                * `2θ Min / Max`: These sliders define the refinement range. It is standard practice to exclude regions of low signal-to-noise or known artifacts from the calculation.

            

            ### Lorentz–Polarisation Factor (Lp)

            
                The integrated intensity of a powder reflection is not $|F|^2$ alone. It carries two
                geometric weights — the multiplicity $m$ of the reflection and the
                Lorentz–polarisation factor $Lp$ — on top of an arbitrary overall scale $s$:
            

            
$$
I(hkl) \;=\; s \cdot m \cdot Lp(2\theta) \cdot |F(hkl)|^2
$$

            
                $Lp$ is therefore the last thing standing between an extracted Le Bail or Pawley
                intensity and a structure factor, and it is the reason the reflection table in the
                report prints $m$, $Lp$ and $|F_o|$ next to the intensity rather than the intensity
                on its own.
            

            > **Note: Lp does not affect the refinement**
> 
                    In both Le Bail and Pawley the intensity of every peak is a free quantity:
                    it is re-partitioned from the observation each cycle, or refined as a
                    least-squares parameter. Either way it absorbs $Lp$ completely. Changing the
                    polarisation model therefore cannot move a lattice parameter, a profile
                    coefficient, $R_{wp}$ or $\chi^2$ by so much as a digit — and if it appears
                    to, something else changed at the same time. What it changes is the
                    *meaning* of the refined intensities: $|F_o|$, the input to charge
                    flipping, the Wilson prior behind the French–Wilson correction, and the
                    space-group probability test. Because nothing needs re-refining, Powder 5
                    applies a change of model retroactively to the live fit and to every run in the
                    history, so re-exporting an old run uses the model currently on screen.
                

            #### The Lorentz factor

            
                For a conventional $\theta$–$2\theta$ powder scan, Powder 5 uses
            

            
$$
L(2\theta) \;=\; \frac{1}{\sin^2\theta \, \cos\theta}
$$

            
                This lumps together the rate at which a crystallite passes through the reflecting
                condition and the fraction of the Debye–Scherrer ring intercepted by the
                detector. It is often written $1/(4\sin^2\theta\cos\theta)$; the factor of 4 is a
                constant and is swallowed by $s$, so it makes no difference to anything reported
                here. $L$ has no adjustable parameters and is applied for every polarisation
                setting, including *None*.
            

            #### The polarisation factor

            
                All three models are the same expression with one constant $K$, normalised so that
                $P(0) = 1$ in every case:
            

            
$$
P(2\theta) \;=\; \frac{1 + K\cos^2 2\theta}{1 + K}
              \qquad\qquad Lp = L \cdot P
$$

            
| Setting | $K$ | Use it when |
| --- | --- | --- |
| **Lab** (default) | $\cos^2 2\theta_M$, or $1$ if $2\theta_M = 0$ | Any sealed-tube or rotating-anode instrument. Leave $2\theta_M$ at 0 for a                             plain unpolarised source; enter the monochromator angle if there is one. |
| **Synchrotron** | $(1-f)/f$ | A polarised source. $f$ is the fraction of the beam polarised                             *perpendicular* to the diffraction plane. |
| **None** | $0$, so $P \equiv 1$ | Neutron data, intensities that have already been corrected elsewhere, or                             a deliberate test of how much the correction is worth. |

            
                With $2\theta_M = 0$ the Lab setting reduces to $P = \tfrac12(1 + \cos^2 2\theta)$,
                the classical unpolarised expression, and this is what Powder 5 uses by default. Note that a fully perpendicular-polarised synchrotron beam ($f = 1$)
                gives $K = 0$ and therefore exactly the same arithmetic as *None*; the two
                settings are kept distinct because they describe different geometry and are labelled
                differently in the report.
            

            #### Choosing $2\theta_M$ for a laboratory monochromator

            
                $2\theta_M$ is the take-off angle of the monochromator crystal at the wavelength in
                use, not a free parameter. Common values for Cu K$\alpha$:
            

            
| Monochromator | $2\theta_M$ (°) | Resulting $K$ |
| --- | --- | --- |
| Graphite (002), the usual diffracted-beam crystal | 26.6 | 0.800 |
| Ge (111) | 27.3 | 0.790 |
| Si (111) | 28.4 | 0.773 |
| Quartz (10$\bar{1}$1) | 26.6 | 0.800 |
| None — Ni filter, or no filter at all | 0 | 1.000 |

            > **Warning:**
> Mosaic versus perfect crystals. $K = \cos^2 2\theta_M$ is the
> ideally mosaic result, which is the right one for graphite and for
> ordinary pyrolytic monochromators. A genuinely perfect crystal diffracts in the
> dynamical regime and the correct constant is $K = |\cos 2\theta_M|$ instead
> — 0.894 rather than 0.800 for graphite geometry. Powder 5 implements
> the mosaic convention. If you need the perfect-crystal value, enter the
> equivalent angle $2\theta_{M}^{\,\mathrm{eff}} = \arccos\sqrt{K}$ (for
> $K = 0.894$, that is 19.1°). In practice the two differ by well under 2% in
> $P$ at any angle, which is smaller than most other systematic errors in a
> laboratory pattern.

            #### Choosing $f$ for a synchrotron

            
                Bending-magnet and undulator radiation is polarised in the *horizontal* plane,
                with $f$ typically between 0.90 and 0.98 in the plane of the orbit. What matters is
                the orientation of the diffractometer relative to that:
            

            
                * **Vertical scattering plane** (the usual arrangement for
                    high-resolution powder diffraction): the electric field is perpendicular to the
                    scattering plane, so $f \approx 0.95$–$1.0$ and $P$ is nearly flat. The
                    correction is almost pure Lorentz, and this is precisely why the geometry is
                    chosen.

                * **Horizontal scattering plane**: the field lies *in* the
                    scattering plane. Enter $f \approx 0.05$–$0.1$, which gives
                    $P \to \cos^2 2\theta$ — a severe correction that goes to zero at
                    $2\theta = 90^\circ$. Getting this backwards is a much worse error than
                    misjudging $f$ by a few per cent.

            

            #### How large is the effect?

            
                $Lp$ falls steeply with angle whatever the model, and the polarisation term adds a
                further factor of two by $2\theta = 90^\circ$ for an unpolarised source. Values for
                Cu K$\alpha$:
            

            
| $2\theta$ (°) | $L$ | $Lp$, Lab (no mono) | $Lp$, graphite | $Lp$, synchrotron $f=0.95$ | $Lp$, None |
| --- | --- | --- | --- | --- | --- |
| 10 | 132.15 | 130.16 | 130.38 | 131.95 | 132.15 |
| 20 | 33.68 | 31.71 | 31.93 | 33.48 | 33.68 |
| 30 | 15.46 | 13.52 | 13.74 | 15.26 | 15.46 |
| 45 | 7.39 | 5.54 | 5.75 | 7.21 | 7.39 |
| 60 | 4.62 | 2.89 | 3.08 | 4.45 | 4.62 |
| 90 | 2.83 | 1.41 | 1.57 | 2.69 | 2.83 |
| 120 | 2.67 | 1.67 | 1.78 | 2.57 | 2.67 |
| 150 | 4.14 | 3.62 | 3.68 | 4.09 | 4.14 |

            
                $Lp$ spans a factor of about 90 between $2\theta = 10^\circ$ and $90^\circ$. Omitting
                it does not perturb a structure-solution calculation slightly — it weights the
                first few peaks roughly two orders of magnitude too heavily, and a charge-flipping map
                built on uncorrected intensities is dominated by them. The choice *between*
                models is a much smaller effect: at most a factor of two, and concentrated near
                $2\theta = 90^\circ$. Getting the correction on at all matters more than getting the
                model exactly right, but the models are cheap to select correctly.
            

            #### Where it is applied

            
| Consumer | What $Lp$ does there |
| --- | --- |
| Reflection table in the report | Printed as its own column, and used for                             $|F_o| = \sqrt{I_{hkl} / (m \cdot Lp \cdot (1 + I_2/I_1))}$. |
| Charge flipping | Each observed intensity is divided by $Lp$ before it becomes a target                             amplitude. Without it the map is dominated by the low-angle peaks. |
| French–Wilson correction | The Wilson prior is estimated from $I/(m\,Lp)$, which is the quantity                             actually proportional to $|F|^2$, rather than from raw areas. |
| Space-group probability test | Intensities are normalised by $m \cdot Lp$ before the Wilson-like scale                             $\tau_j$ is estimated, so the extinction evidence is not confounded with                             the $Lp$ fall-off. |
| Theoretical (data-free) export | $Lp$ and $m \cdot Lp$ are tabulated per reflection. With no structure                             factors available, $m \cdot Lp$ *is* the whole of the predicted                             relative intensity. |

            
                All of these read one shared calculation, evaluated at the *corrected*
                $2\theta$ — the same angle printed in the reflection table, including the zero
                shift and any displacement or transparency term. The charge-flipping worker takes the
                per-reflection value computed by the refinement rather than recomputing its own, so
                the map, the report and the prior cannot disagree about it. The charge-flipping
                summary panel reports the model the worker actually used, which is the check that
                they have not.
            

            > **Note: A note on $|F_o|$ and the Kα doublet**
> 
                    With a K$\alpha_1$/K$\alpha_2$ source the integrated area of a reflection is the
                    sum over the pair, so it exceeds the single-wavelength intensity by
                    $(1 + I_2/I_1)$. The $|F_o|$ column divides that back out, which makes the printed
                    structure factors comparable with a monochromatic or synchrotron measurement
                    instead of carrying a source-dependent constant. It is a single global factor and
                    so has no effect on *relative* $|F_o|$; the report states the divisor it
                    used. The overall scale remains arbitrary in any case: $|F_o|$ is on the scale of
                    the observed pattern, not on an absolute electron scale.
                

            ### Background Modeling (Spline Interpolation)

            
                The background contribution is modeled using a **monotonic cubic Hermite spline interpolation** between a series of user-defined points (spline points or knots). This approach provides flexibility and ensures a smooth, physically realistic background shape without introducing refinable background parameters into the least-squares minimization. The background shape is therefore considered **fixed** during the refinement process based on the current spline points.
            

            
                Control of the background spline is located under the "Background" tab:
            

             
                * **Auto-generation:** The application **automatically estimates** background points immediately upon loading a data file. You can adjust the density of these points using the **Auto-points slider** (10-40 points). Adjusting the slider automatically recalculates the points based on local intensity minima within intervals distributed across the current **2θ Min/Max slider** range. The points at the exact 2θ Min and Max slider positions are always included and fixed to these $2\theta$ values.

                * **Manual Addition:** Add individual points by holding **Ctrl** and clicking on the chart. The closest experimental point will be added to the list, provided it's within the current slider range and not an edge point.

                * **Editable List:** The generated and manually added points appear in a list below the controls.
                    
                        You can directly edit the **$2\theta$** and **Intensity** values for any point, except for the $2\theta$ values of the first (Min) and last (Max) points, which are fixed by the sliders. Edits trigger recalculation of the spline.

                        * Points can be deleted using the **×** button, except for the first and last points.

                    
                
                 * **Chart Display:** The spline points can be toggled on/off on the chart using the "Show Points on Chart" checkbox. The calculated spline curve is always shown.

             
             > **Note: Isotropic Broadening (TCH Model)**
> The TCH formulation models the FWHM of the Gaussian ($H_G$) and Lorentzian ($H_L$) components based on physical contributions to line broadening:
                
$$
H_G^2 = U \tan^2\theta + V \tan\theta + W
$$

                
$$
H_L = X \tan\theta + Y / \cos\theta
$$

            

            The total FWHM ($H$) and pseudo-Voigt mixing parameter ($\eta$) are then derived from these components using polynomial approximations. The final shape is $pV(x) = \eta L(x, H) + (1-\eta)G(x, H)$, where both functions share the same convoluted FWHM.

            
                * `U, V, W`: Gaussian broadening parameters related to strain ($U, V$) and instrumental resolution ($W$). Physically, $U$ and $W$ should be non-negative.

                * `X, Y`: Lorentzian broadening parameters related to strain ($X$) and crystallite size ($Y$). Physically, $X$ and $Y$ should be non-negative.

            

            #### Peak Asymmetry

             The `S/L` and `H/L` parameters introduce an angle-dependent asymmetry, primarily correcting for axial divergence effects at low $2\theta$.

            #### Anisotropic Broadening (Stephens Model)

            
                Anisotropic microstrain, where broadening varies with crystallographic direction, is modeled by adding terms to the Lorentzian component ($H_L$) that are dependent on the Miller indices ($hkl$). The model is a fourth-order polynomial in the reciprocal lattice vectors.
            

            
                The refinable parameters (`S400`, `S040`, etc.) are the non-zero, symmetry-unique coefficients of this polynomial. The application automatically applies symmetry constraints based on the Laue class of the selected space group (e.g., for cubic, $S400=S040=S004$).
            

            > **Unit and Scaling for Stephens `S_hkl` parameters:**
> 

                The user-inputted `S_hkl` parameters are scaled for convenience. The dimensionless anisotropic broadening term ($H_{aniso}$) is calculated from these parameters, and its contribution to the total Lorentzian width (in degrees $2\theta$) is scaled by a factor of 1000.
                
$$
H_L(\text{total}) = H_L(\text{isotropic}) + \frac{|H_{aniso}|}{1000}
$$

                This scaling allows the user to refine values in a manageable range (e.g., -10 to +10) rather than requiring input of very small decimals (e.g., 1e-4), a convention common in other refinement software.
          ### Split pVoigt (Asymmetric)

            
                This profile function is a modification of the Simple pseudo-Voigt designed to model asymmetric peaks (e.g., from axial divergence or stacking faults). It achieves this by defining independent sets of profile width parameters for the left side (at $2\theta$ values less than the peak center) and the right side of the peak.
            

            
                The shape is still a linear combination $pV(x) = \eta L(x) + (1-\eta)G(x)$, but the $H_G$ and $H_L$ parameters used in the calculation depend on whether $x$ is to the left or right of the peak center.
            

            #### Gaussian Broadening (Left & Right)

            The Gaussian FWHM ($H_G$) for each side is modeled as:

            
$$
H_{G, \text{side}}^2 = GU_{\text{side}} \tan^2\theta + GV_{\text{side}} \tan\theta + GW_{\text{side}}
$$

            (Note: This model does not use the $GP$ term used in the Simple pVoigt profile.)

            
                * `GU-L, GV-L, GW-L`: Parameters describing the Gaussian FWHM for the **left side** of the peak.

                * `GU-R, GV-R, GW-R`: Parameters describing the Gaussian FWHM for the **right side** of the peak.

            

            #### Lorentzian Broadening (Left & Right)

            The Lorentzian FWHM ($H_L$) for each side is modeled as:

            
$$
H_{L, \text{side}} = LX_{\text{side}} / \cos\theta
$$

            
                * `LX-L`: Describes the Lorentzian FWHM for the **left side**, primarily associated with size broadening.

                * `LX-R`: Describes the Lorentzian FWHM for the **right side**, primarily associated with size broadening.

            

            #### Peak Shape & Position

            
                * `eta (Mixing)` (param: `eta_split`): A simple linear mixing parameter ($0 \le \eta \le 1$; $\eta=0$ for pure Gaussian, $\eta=1$ for pure Lorentzian). This single value is used for both sides of the peak.

                * `shft (Displ.)` (param: `shft_split`): Corrects for peak position shifts due to sample displacement.

                * `trns (Transp.)` (param: `trns_split`): Corrects for peak position shifts due to transparency.
                    
                        **Unit and Scaling for the `shft_split` parameter:**

                        This parameter is a **dimensionless, scaled coefficient**, identical in function to the `shft` parameter in the Simple pVoigt profile.
                        
                            The physical peak shift in radians is: $\Delta(2\theta)_{\text{rad}} = -2 \frac{s}{R} \cos(\theta)$

                            * The program calculates this shift (in degrees) using the formula: $\Delta(2\theta)_{\text{deg}} = -(\text{shft\_split} / 1000) \times \cos(\theta) \times (180 / \pi)$

                            * Therefore, the relationship is: $\frac{\text{shft\_split}}{1000} = \frac{2s}{R}$

                            * To find the physical displacement $s$ from the refined parameter, use: $s = R \times (\text{shft\_split} / 2000)$.

                        
                
            
        
        
---

        
            ## Probabilistic Space-Group Determination

            
                Once a Pawley refinement has completed, a small **?** button becomes
                active on the space-group line of the control panel. It scores every space group
                compatible with the current Laue class against the extracted intensities and
                returns a ranked list of posterior probabilities. The method follows
                Markvardsen, David, Johnson & Shankland (2001).
            

            ### The Principle

            
                A candidate space group $S$ makes exactly one testable statement about a powder
                pattern: a particular subset $e$ of the reflections is *systematically absent*,
                so their true integrated intensities are identically zero. Everything else it leaves
                free. The question is therefore not “does this space group fit?” but
                “are the reflections it forbids consistent with zero, given how well the data
                actually determine them?”
            

            
                A Pawley refinement is what makes this answerable. Because every intensity is a free
                least-squares parameter, the refinement returns not only the intensities
                $\hat{\mathbf{I}}$ but the normal matrix behind them, and therefore their full
                *covariance*. That matters enormously: in a powder pattern reflections overlap,
                and overlapping intensities are strongly correlated. A reflection that should be
                absent can carry a large apparent intensity purely by borrowing it from a neighbour
                it cannot be resolved from. Judging each reflection on its own $I/\sigma$ ignores
                this and can be badly misleading in either direction.
            

            > **Note: 1. Orbits: symmetry equivalents share one intensity**
> 
                A Pawley fit reports one intensity per symmetry-unique reflection. Its
                equivalents $hR$ must be placed on the $P1$ grid before a transform
                is possible, and they are generated from the *actual symmetry
                operators* of the selected group. The distinction matters: the Laue
                class is not the holohedry of the crystal system. $P4$ has Laue class
                $4/m$, so a general reflection has 8 equivalents, not the 16 that
                $4/mmm$ would give; the same applies to $\bar 3$ versus $\bar 3 m$,
                $6/m$ versus $6/mmm$ and $m\bar 3$ versus $m\bar 3 m$. Spreading
                one measured intensity over twice as many positions as really exist puts
                half of it onto reflections that are not equivalent at all.
            

            
                Within an orbit the amplitudes are *equal by symmetry*:
                $|F| = \sqrt{I^{obs}/m}$, where $m$ is the orbit size. This is a hard
                constraint, applied every cycle. Only the phases are free.
            

            #### 2. Overlap: distinct reflections may share one peak

            
                The genuinely powder-specific ambiguity is different. Reflections that
                are *not* symmetry equivalent can still fall at the same
                $2\theta$, either exactly (511/333 in a cubic cell) or accidentally.
                The pattern gives only their sum. Powder 5 groups such reflections
                into a *cluster* and re-partitions the cluster's total intensity
                between them each cycle in proportion to the current calculated values
                — the same idea as Le Bail extraction, applied inside the
                flipping loop. The split is free; the total is fixed.
            

            > **These two mechanisms are easy to confuse and they pull in opposite
                directions. **Symmetry equivalents must be equal**
> (constraint);
                overlapping distinct reflections are free to differ**
                (unknown). Letting the members of an orbit float would discard exactly
                the information the space group provides.

            #### 3. Systematic absences are held at zero

            
                Reflections forbidden by lattice centring, screw axes or glide planes are
                forced to $F = 0$ every cycle. They are found directly from the
                operators: $h$ is extinct if some operator $(R,\mathbf{t})$ satisfies
                $hR = h$ with $h \cdot \mathbf{t}$ non-integral. Leaving them free
                lets the density break the lattice centring outright, and for an $F$ or
                $I$ lattice that is most of the reciprocal grid.
            

            #### 4. Symmetrisation of the phases

            
                For a real-space operator $\mathbf{x}' = \mathbf{x}R + \mathbf{t}$, the
                structure factors of a correctly positioned structure satisfy
            

            
$$
F(hR) = F(h)\,\exp(-2\pi i\, h \cdot \mathbf{t})
$$

            
                Each cycle the members of an orbit are reduced to a common
                representative, averaged, and pushed back out. The
                **Symmetry in loop** control sets how strongly:
                $F \leftarrow (1-\lambda)F + \lambda F_{sym}$. At $\lambda = 0$ nothing
                is imposed and the run is classical $P1$ charge flipping. At
                $\lambda = 1$ the projection is exact, which converges faster and has a
                useful side effect: the solution comes out on a standard origin instead
                of a random one, so the automatic origin search in the Structure block
                becomes a verification rather than a necessity.
            

            > **Warning: The trade-off.**
> Strict symmetry cannot recover from a
                wrong space group — it will happily converge to a low $R$ for a
                structure that does not exist. If the group is uncertain, run at
                $\lambda = 0$ as well and compare: a density that is genuinely
                consistent with the group will symmetrise well afterwards (high symmetry
                correlation) without having been forced to.

            #### Intensity scaling

            
                The refined Pawley parameter is a peak *height*, because the
                profile functions Powder 5 uses are normalised to unit height, not
                unit area. What charge flipping needs is the integrated intensity, so
                each height is multiplied by the area of its own profile (summed over
                the Kα1/Kα2 doublet) before it is used.
                Skipping that step underweights the high-angle reflections by the full
                Caglioti broadening, which is roughly a factor of two across a typical
                scan.
            

            
                The integrated intensity then carries
                $I \propto m \cdot LP \cdot |F|^2$. The multiplicity $m$ is absorbed
                by spreading the intensity across the $m$ grid points of the orbit, but
                the Lorentz-polarisation factor
                $LP = (1+\cos^2 2\theta)/(2\sin^2\theta\cos\theta)$ is divided out
                explicitly. Without it the low-angle reflections are weighted several
                times too heavily and the map is dominated by the first few peaks.
            

            > **Note: Choosing the grid**
> 
                The grid must be large enough to hold every reflection in the fit: an
                index $h_{max}$ needs $N \geq 2h_{max}+2$. If it is too small the
                summary reports how many reflections did not fit and what $N$ they
                would need, and those reflections are simply lost. Beyond that
                requirement, aim for a voxel of roughly 0.25–0.35 Å
                (the summary prints the spacing). Finer than ~0.2 Å buys
                nothing when the data stop at 1 Å resolution and costs
                $N^3\log N$: 64³ to 128³ is about an eight-fold increase in
                run time.
            

            > **A 128³ grid holds about 50 MB of GPU memory and pushes several
                dispatches close to the WebGPU ceiling of 65535 workgroups per
                dimension. Powder 5 checks every dispatch against the adapter's
                reported limits before it encodes anything; if the device cannot take
                the grid, the run falls back to the CPU with a message saying so rather
                than producing a map that was never computed. The CPU path at 128³
                is slow — expect minutes per random start — so on a
                constrained device it is usually better to stay at 64³.

            #### Choosing the threshold

            
                $\delta/\sigma$ is the single most sensitive parameter. Too low and
                almost nothing is flipped, so the iteration stalls with $R$ flat and
                high; too high and the perturbation is so violent that the phases never
                settle, giving an $R$ that jumps around without trend. The working
                range is 0.8–1.2 and 1.1 is a good default. If $R$ plateaus above
                0.3 with no oscillation, *lower* it in steps of 0.1; if $R$ is
                erratic from the first cycles, *raise* it. Weak or noisy data
                generally want a slightly lower threshold.
            

            #### Iterations and random starts

            
                Charge flipping does not converge gradually — it stays flat and
                then drops abruptly when the phases lock in, often after several hundred
                cycles. Judge a run by whether that drop happened, not by the final
                $R$. More *starts* is usually a better investment than more
                *iterations*: a run that has not found the solution by ~1000
                cycles is unlikely to find it by 3000, whereas a fresh random start might
                find it immediately. For a difficult problem, 10–20 starts of 1000
                cycles beats 2 starts of 10000.
            

            #### Choosing λ

            
                Use $\lambda = 1$ when the space group is established (a clean
                space-group test, a known compound, an unambiguous extinction pattern):
                it converges in fewer cycles, needs fewer random starts, and delivers the
                map already on a standard origin. Use $\lambda = 0$ when the group is
                in doubt, or as an independent check on a $\lambda = 1$ solution —
                if the $P1$ map symmetrises to a high correlation without having been
                forced to, the group is right. $\lambda = 0.5$ is the compromise and
                the default: it accelerates convergence appreciably while still allowing
                the density to disagree with a slightly wrong group.
            

            #### Choosing the overlap tolerance

            
                Set it to roughly half the FWHM of a mid-angle peak. Too small and truly
                overlapped reflections are treated as independently measured, which
                imposes intensities the data never determined; too large and reflections
                that *were* resolved get to trade intensity freely, discarding
                real information. 0.05° suits sharp laboratory data; broad peaks or a
                high-symmetry cell with many coincidences want 0.10–0.15°. The
                summary reports how many reflections ended up sharing a peak — if
                that number is close to the total, the tolerance is almost certainly too
                large.
            

            > **Note: Peak and site table columns**
> 
| Column | Meaning |
| --- | --- |
| **#** > | Peak rank by height (1 = strongest). |
| x, y, z** | Fractional coordinates. In the raw peak list these are in $P1$ with an arbitrary origin; in the built structure they are placed on the space-group origin. |
| **Height** | Peak density relative to the strongest peak (which is set to 1). A rough proxy for atomic number: heavier atoms scatter more and peak higher. |
| **dmin (Å)** | (Peak list.) Distance to the nearest other peak — more meaningful than the absolute coordinates, since the origin is arbitrary. |
| **Mult** | Site multiplicity: how many symmetry-equivalent copies of this site exist in the unit cell. See the note on the asterisk below. |

            > **Note:**
> The asterisk (*) on a multiplicity marks a
> special position. A site lying on a symmetry element — a
> mirror plane, rotation axis or inversion centre — maps onto itself
> under one or more symmetry operators, so it has fewer
> equivalent copies than a general position. Its multiplicity is therefore
> smaller than the space group's general multiplicity, and the asterisk
> flags exactly that: this atom sits on a special position and its
> coordinates are constrained by symmetry. For example, in $Pnma$ a general
> site has multiplicity 8; an atom on the mirror at $y=\tfrac14$ has
> multiplicity 4, shown as 4*. Special positions are common and
> expected — many atoms in real structures sit on them — and
> recognising them is the first step toward assigning Wyckoff letters.
> Requirements and Limitations
> Data quality. The solution is only as good as the
> extracted intensities. Severe peak overlap, a poor background, or an
> incorrect cell will all degrade or prevent a solution.
> Resolution. Atomic-resolution data (roughly
> $d \lesssim 1.2$ Å) is normally needed to resolve
> individual atoms; lower resolution gives a blurred map.
> The space group is an assumption. With
> $\lambda > 0$ it is imposed, not tested. A wrong group can still
> produce a low $R$, because the constraint makes the calculated
> intensities agree with the observed ones by construction. Cross-check
> with a $\lambda = 0$ run and with the symmetry correlation from
> the Structure block.
> Non-uniqueness. Random starts can land on different
> origins or the inverted structure. This is expected at
> $\lambda = 0$; the structure builder resolves origin and hand
> against the space-group symmetry.
> Light atoms near heavy ones. As with any
> Fourier method, light atoms sitting close to much heavier ones may be
> lost in truncation ripples.
> Overlap is modelled, not solved. Re-partitioning the
> intensity of a cluster in proportion to the calculated values is a
> reasonable guess, not a measurement. Structures whose solution
> depends on resolving a specific pair of overlapped reflections may be
> beyond reach from powder data alone.
> Recommended Refinement Strategy
> A sequential and hierarchical refinement strategy is crucial for achieving a stable and physically meaningful solution. Attempting to refine all parameters simultaneously from a poor starting model will likely lead to divergence or convergence to a false minimum.
> Phase 1: Initial Model Setup
> Define the Model: Load data, select the crystal system and space group, and define the refinement range using the $2\theta$ sliders.
> Set Background Points: The background is automatically estimated upon loading. Adjust the Auto-points slider under the "Background" tab to optimize the density of points so that they reasonably follow the experimental background. You can also edit the points manually or use Ctrl+Click on the chart. This background shape is fixed during refinement.
> Peak Position Refinement: Using the Le Bail and LM algorithms, refine only the Lattice Parameter(s) and, if necessary, the instrumental Zero Shift. The goal is to align the calculated Bragg positions with the observed peak maxima.
> Phase 2: Peak Profile Refinement
> Isotropic Broadening: Once positions are correct, refine the primary isotropic peak shape parameters (e.g., W, Y, U, X in TCH, or GW, LX in Simple/Split pVoigt). This will account for the dominant size and strain contributions.
> Asymmetry and Shape: Introduce asymmetry parameters (S/L in TCH, shft/trns in Simple/Split pVoigt) if there is a clear misfit, particularly at low angles. Refine the mixing parameter (eta) if needed.
> Anisotropic Broadening (TCH only): If systematic misfits remain (e.g., some peaks are consistently broader than the model), introduce the anisotropic Stephens parameters (e.g., S400). Refine only the symmetry-independent terms.
> Phase 3: Finalization and Intensity Extraction
> Global Optimization (Optional): If the LM algorithm converges to a poor solution during the Le Bail steps, switch to Parallel Tempering (PT) for one or more Le Bail cycles to perform a global search. Afterwards, switch back to LM for a final, precise local minimization.
> Pawley Refinement: With a stable and well-refined model from the Le Bail method, perform a final refinement using the Pawley method. It is generally recommended to use the Levenberg-Marquardt (LM) algorithm for stability, though Parallel Tempering (PT) can also be used (potentially requiring more iterations). This will provide the most statistically robust set of integrated intensities, suitable for subsequent structure solution.
> > **Note: Technical Note on Calculation Parameters**
> >
> Each reflection is evaluated only within a finite window around its centre,
> set by CALCULATION_WINDOW_MULTIPLIER (currently 8.0 × the
> peak FWHM, doubled when an asymmetry correction is active). This radius is the
> single authoritative truncation distance, shared by the preview and the
> refinement worker.
> Rather than discarding contributions below a fixed fraction of peak height, the
> profile value at the truncation radius is subtracted from the whole
> peak as a small constant “pedestal”. The profile therefore reaches
> zero continuously at the window edge instead of stepping off a cliff. An
> earlier version used a separate height cutoff, which put a small discontinuity
> in the calculated pattern — and hence in the derivatives — wherever
> a peak was truncated.
> Interpretation of Results
> Assessing the quality of a refinement requires both statistical analysis and critical visual inspection of the fit.
> Figures of Merit
> Standard crystallographic R-factors are provided to quantify the quality of the fit.
> R-pattern ($R_p$): The unweighted residual error based on net intensities, sensitive primarily to the fit of high-intensity reflections.
> 
$$
R_p = \frac{\sum |(y_{i,obs} - y_{i,bkg}) - (y_{i,calc} - y_{i,bkg})|}{\sum |y_{i,obs} - y_{i,bkg}|} \times 100\%
$$

> Weighted R-pattern ($R_{wp}$): The primary figure of merit, weighted by the inverse of the observed gross intensities ($w_i = 1/y_{i,obs}$), which properly accounts for the counting statistics across the entire pattern.
> 
$$
R_{wp} = \left[ \frac{\sum w_i (y_{i,obs} - y_{i,calc})^2}{\sum w_i y_{i,obs}^2} \right]^{1/2} \times 100\%
$$

> Reduced Chi-squared ($\chi^2$, Goodness of Fit): The most statistically rigorous indicator. For a statistically perfect fit where the model correctly describes the data and the weights are accurate, $\chi^2$ should approach 1.0.
> 
$$
\chi^2 = \frac{1}{N - P} \sum w_i (y_{i,obs} - y_{i,calc})^2 = \left(\frac{R_{wp}}{R_{exp}}\right)^2
$$

> where $N$ is the number of data points, $P$ is the number of refined parameters, and $R_{exp}$ is the statistically expected minimum $R_{wp}$.
> Calculating Observed Intensities ($I_{obs}$) for Overlapping Peaks
> A simple numerical integration over a fixed angular range is insufficient for accurately determining the observed integrated intensity ($I_{obs}$) of overlapping peaks. This tool employs a more robust intensity partitioning method.
> At each point in the diffraction pattern, the net observed intensity ($y_{obs} - y_{bkg}$) is distributed among all contributing Bragg reflections. This distribution is proportional to the value of each peak's calculated profile function (including both Kα1 and Kα2 components, scaled by the refined peak height) at that specific point. By integrating these partitioned "slices" of intensity for each reflection across the entire pattern (using the trapezoidal rule), the method yields a reliable $I_{obs}$ value (reported as integrated area) that correctly deconvolutes contributions from neighboring peaks.
> Visual Inspection
> Numerical indicators can be misleading. Visual inspection of the difference plot (observed minus calculated) is the most critical step in evaluating the fit.
> The Difference Plot: A successful refinement should yield a difference plot that consists of random, uncorrelated noise centered on zero. The plot is scaled relative to the main pattern for visibility.
> Systematic Residuals: The presence of structured, non-random features in the difference plot (e.g., "M-shaped" residuals around peaks, broad humps where the spline is inadequate, or un-indexed peaks) is a clear indication of systematic errors in the model. These may arise from an incorrect peak shape, unmodeled anisotropy or asymmetry, an inadequate background model (requiring adjustment of spline points), or the presence of an unaccounted-for impurity phase.
> Williamson-Hall Size-Strain Analysis
> For refinements utilizing the TCH (Size/Strain/Aniso) profile function, the application can automatically perform a Williamson-Hall analysis to extract approximate microstructural information. This method separates the contributions of crystallite size and microstrain to the total peak broadening by analyzing their different dependencies on the diffraction angle, $\theta$.
> The analysis is based on the linear Williamson-Hall equation, where $\beta$ is the total physical peak breadth (FWHM) in radians derived from the refined sample-only broadening parameters (U, X, Y), excluding instrumental contributions (V, W):
> 
$$
\beta \cos(\theta) = \frac{K\lambda}{L} + 4\epsilon \sin(\theta)
$$

> This equation describes a straight line when plotting $\beta \cos(\theta)$ vs. $4\sin(\theta)$. The software performs a linear least-squares fit on the relevant data points (within the fitted 2θ range) to determine the y-intercept (related to crystallite size, $L$) and the slope (related to microstrain, $\epsilon$).
> > **Note: Reported Values**
> >
> Apparent Crystallite Size (nm): An estimate of the average size of the coherently scattering domains, calculated from the y-intercept of the Williamson-Hall plot (using $K=0.9$).
> Apparent Microstrain (%): An estimate of the root-mean-square strain within the crystallites, calculated from the slope of the plot.
> Linear Fit R²: The coefficient of determination for the linear regression. A value close to 1.0 indicates that the isotropic size/strain model is a good fit for the observed peak broadening. Values significantly less than 1.0 may suggest that broadening is anisotropic or that the model is otherwise inadequate.
> Note: This Williamson-Hall analysis provides an approximation based only on the isotropic TCH parameters (U, X, Y). It does not account for anisotropic broadening effects (Stephens parameters) or instrumental contributions (V, W). For rigorous quantitative analysis, dedicated size/strain analysis software should be employed.
> Data Export
> Save Report: Behavior depends on whether experimental data is loaded.
> With data loaded, after a refinement: generates a comprehensive ASCII text file containing all statistical indicators, refined parameter values with their ESDs (for LM refinements), the Williamson-Hall analysis results (if applicable), the declared Lorentz–polarisation model, the list of background spline points used, the reflection table described below, and a point-by-point list of observed, calculated, background, and difference intensities across the fitted range.
> With no data loaded: the button remains active and instead exports a theoretical reflection list for the current space group, lattice parameters, and wavelength(s) — h k l, $d$-spacing, and $2\theta$ for Kα1 (and Kα2 if the intensity ratio is greater than 0), together with the multiplicity, $Lp$ and $m \cdot Lp$ of each reflection. With no structure factors available, $m \cdot Lp$ is the whole of the predicted relative intensity. No refinement is required for this mode; the chart itself also displays this theoretical stick pattern live as you change the space group or lattice, so you can preview expected peak positions before ever loading a scan.
> Generate PDF: Creates a summary PDF document containing a high-resolution plot and tables of final parameters and statistics, suitable for archival or reporting. This requires a completed refinement. The PDF carries the same reflection table as the text report; rows too wide for the page are set in a smaller size rather than being clipped.
> The reflection table
> Both the text and PDF reports list every reflection inside the fitted range with the
> complete chain from peak to structure factor. Le Bail and Pawley share the same
> columns; Pawley adds one.
> ColumnMeaning
> h,k,l
> Canonical member of the reflection orbit. One row per orbit, not per
> index triple.
> 2th_corr
> Peak position after the zero shift and any displacement or transparency
> correction — where the peak actually sits, not the ideal Bragg
> angle.
> dInterplanar spacing in Ångström.
> m
> Powder multiplicity: the size of the orbit under the Laue group. See
> Reflections & Multiplicities.
> Lp
> Lorentz–polarisation factor, evaluated at 2th_corr with the
> declared model. See Lorentz–Polarisation
> Factor.
> I_hkl
> Integrated intensity of the whole reflection — peak height times
> profile area, summed over the Kα1/Kα2 pair, with the scale
> factor applied. This is an area, not a height: a height
> underweights the high-angle reflections by the full profile broadening,
> which typically doubles across a scan.
> sigma(I)
> Standard uncertainty on the area. For Pawley it comes from the
> covariance matrix, converted from a height ESD; for Le Bail it is the
> counting statistics of the observed points propagated through the
> partition fractions.
> |Fo|
> $\sqrt{I_{hkl} / (m \cdot Lp \cdot (1 + I_2/I_1))}$, on an arbitrary
> scale. Blank for a non-positive intensity — see the note
> below.
> I_obs (Pawley only)
> The same reflection re-integrated from the observed pattern by
> partitioning the net counts over the calculated profile. A cross-check on
> the refined parameter, which is what I_hkl holds for a Pawley
> fit. Large disagreement points at an overlap the fit has resolved badly.
> For Le Bail the two are the same integral by construction, so the
> column is omitted.
> A legend above the table restates the polarisation model, the formulae and the value
> of $K$ actually used, so every $|F_o|$ can be reproduced from the printed numbers
> without knowing which convention the program follows.
> > **Note: Negative intensities are printed, not hidden**
> >
> A least-squares or partitioned intensity can come out negative when the
> background runs above the data across a peak's window. That is a real
> measurement and I_hkl shows it; only $|F_o|$ is left blank, because
> there is no square root to take. Clipping such values at zero would bias every
> weak reflection upward, which is exactly what the French–Wilson correction
> used by the charge-flipping path exists to avoid.
> References & Further Reading
> Pawley Method:
> Pawley, G. S. (1981). "Unit-cell refinement from powder diffraction scans". Journal of Applied Crystallography, 14(6), 357-361.
> Le Bail Method:
> Le Bail, A., Duroy, H. & Fourquet, J.L. (1988). "Ab-initio structure determination of LiSbWO6 by X-ray powder diffraction". Materials Research Bulletin, 23(3), 447-452.
> Charge Flipping:
> Oszlányi, G. & Sütő, A. (2004). "Ab initio structure solution by charge flipping". Acta Crystallographica Section A, 60(2), 134-141.
            
                **Charge Flipping — practice and variants:**

                Oszlányi, G. & Sütő, A. (2008). "The charge flipping algorithm". *Acta Crystallographica Section A*, 64(1), 123-134.
            
            
                **Charge Flipping with symmetry and powder data:**

                Palatinus, L. & Chapuis, G. (2007). "SUPERFLIP – a computer program for the solution of crystal structures by charge flipping in arbitrary dimensions". *Journal of Applied Crystallography*, 40(4), 786-790.
            
            
                **Charge Flipping from powder data (intensity repartitioning):**

                Baerlocher, C., McCusker, L. B. & Palatinus, L. (2007). "Charge flipping combined with histogram matching to solve complex crystal structures from powder diffraction data". *Zeitschrift für Kristallographie*, 222(2), 47-53.
            
            
                **Polarisation with a crystal monochromator:**

                Azaroff, L. V. (1955). "Polarization correction for crystal-monochromatized X-radiation". *Acta Crystallographica*, 8(11), 701-704.
            
            
                **Lorentz and polarisation factors for powders:**

                Klug, H. P. & Alexander, L. E. (1974). *X-Ray Diffraction Procedures for Polycrystalline and Amorphous Materials*, 2nd ed. Wiley, New York.
            
            
                **Polarisation conventions in modern refinement software:**

                Toby, B. H. & Von Dreele, R. B. (2013). "GSAS-II: the genesis of a modern open-source all purpose crystallography software package". *Journal of Applied Crystallography*, 46(2), 544-549. (Source of the single-parameter polarisation convention used here.)
            
            
                **Posterior intensities for weak and negative reflections:**

                French, S. & Wilson, K. (1978). "On the treatment of negative intensity observations". *Acta Crystallographica Section A*, 34(4), 517-525.
            
            
                **GSAS Profile Functions:**

                Larson, A. C. & Von Dreele, R. B. (2004). "General Structure Analysis System (GSAS)". *Los Alamos National Laboratory Report LAUR 86-748*.
            
             
                **TCH Profile Function:**

                Thompson, P., Cox, D. E. & Hastings, J. B. (1987). "Rietveld refinement of Debye-Scherrer synchrotron X-ray data from Al2O3". *Journal of Applied Crystallography*, 20(2), 79-83.
            
            
                **Stephens Anisotropy Model:**

                Stephens, P. W. (1999). "Phenomenological model of anisotropic peak broadening in powder diffraction". *Journal of Applied Crystallography*, 32(2), 281-289.
            
            
                **Parallel Tempering:**

                Swendsen, R. H., & Wang, J. S. (1986). "Replica Monte Carlo simulation of spin-glasses". *Physical Review Letters*, 57(21), 2607.
            
            
                **Monotonic Splines:**

                 Fritsch, F. N., & Carlson, R. E. (1980). "Monotone Piecewise Cubic Interpolation". *SIAM Journal on Numerical Analysis*, 17(2), 238–246.
            
            
                **Probabilistic Space-Group Determination:**

                Markvardsen, A. J., David, W. I. F., Johnson, J. C. & Shankland, K. (2001). "A Probabilistic Approach to Space-Group Determination from Powder Diffraction Data". *Acta Crystallographica Section A*, 57(1), 47-54.
            
            
                **cctbx - Computational Crystallography Toolbox:**

                Grosse-Kunstleve, R. W., Sauter, N. K., Moriarty, N. W., & Adams, P. D. (2002). "The Computational Crystallography Toolbox: crystallographic algorithms in a reusable software framework". *Journal of Applied Crystallography*, 35(1), 126-136. (Used for space group systematic absence rules).
            
        

        
---

        
            ## About This Tool

            
                The **powder5** toolkit was developed by Nita Dragoe from Université Paris-Saclay as a simple browser-based implementation of powder pattern decomposition methods. It is a long-time successor of PowderV2 (Dragoe, N. (2001). *J. Appl. Cryst.*, 34, 535) and has been updated to include the Pawley method and modern global optimization algorithms.
            

            
                All numerical work — the linear solver, the Gauss-Jordan/Cholesky routines
                behind the covariance matrix and the ESDs, and the least-squares machinery itself
                — is implemented locally. An earlier version depended on
                [math.js](https://mathjs.org/) for these; it was removed
                so that a refinement can never be blocked by a content-delivery network, and so the
                application works fully offline. Charting uses
                [Chart.js](https://www.chartjs.org/); all pan, zoom and
                pinch behaviour is implemented directly against the chart's scales, with no plugin.
                PDF generation uses [jsPDF](https://github.com/parallax/jsPDF)
                and [html2canvas](https://html2canvas.hertzen.com/).
                Space-group symmetry data is derived from the
                [Computational Crystallography
                Toolbox (cctbx)](https://cctbx.github.io/).
            

            
                This document was updated with the assistance of an AI.
            

            > **Disclaimer:**
> This application is provided for educational and research purposes. While it implements standard and robust algorithms, it is not a substitute for fully validated, peer-reviewed software packages (e.g., GSAS-II, FullProf, TOPAS) for analyses intended for publication.
        
        
            Powder 5, version 147, 08 August 2026