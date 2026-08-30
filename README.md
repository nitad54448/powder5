# powder5

Whole-pattern analysis of powder X-ray diffraction data, in a browser.
[Launch it](https://nitad54448.github.io/powder5/powder5.html) — there is
nothing to install, and your data never leaves your machine.

**powder5** decomposes a powder pattern against a unit cell and a space group,
without needing a structural model. From there it can go on to solve one: the
extracted intensities feed a charge-flipping map or a direct-space Wyckoff
search, both of which run in the browser too.

---

## What it does

**Pattern decomposition** — refine the cell, zero-point, peak profile and
background, and extract integrated intensities, by the **Le Bail** or
**Pawley** method. Three profile functions: simple pseudo-Voigt, the
Thompson–Cox–Hastings model with Stephens anisotropic broadening, and a split
pseudo-Voigt for asymmetric peaks. Refinement runs on a background thread, by
Levenberg–Marquardt or by parallel tempering with an LM polish.

**Space-group determination** — test the settings compatible with the observed
extinctions against the Pawley intensities, and rank them.

**Structure solution** — charge flipping (WebGPU where available, CPU
otherwise), or a Wyckoff-position swarm search that fits atoms directly to the
Pawley intensities. Either produces coordinates, a structure plot, contacts and
bond-valence sums, and a CIF.

**Reporting** — statistics, refined parameters with standard uncertainties,
reflection tables, and export to text, CIF or PDF.

---

## Using it

1. **Load data.** Drop in a file, or use *Open*. Bruker (`.brml`, `.uxd`,
   `.raw`), PANalytical (`.xrdml`), Rigaku (`.ras`, `.rasx`), Philips (`.udf`),
   GSAS (`.gsa`, `.esd`), CIF, and plain two-column ASCII are all read
   directly. Where the file carries the wavelength, it is filled in for you —
   check it.
2. **Set the cell and space group** on the *Sample* tab, and trim the 2θ range
   with the sliders.
3. **Place background anchors** on the *Background* tab, or Ctrl-click the
   plot. The background is a spline through those anchors; it is interpolated,
   not refined.
4. **Choose a profile** on the *Profile* tab. Start with the simple
   pseudo-Voigt.
5. **Refine.** Le Bail first, to get the cell and profile settled — it is
   faster and far more forgiving of a poor starting point. Then Pawley, when
   you want individual intensities and their uncertainties.
6. **Read the difference curve**, not just R<sub>wp</sub>. Structure in the
   residual means the model is wrong somewhere, whatever the number says.
7. **Solve**, if that is where you are going: *Charge Flipping* or *Wyckoff*,
   both of which need a converged Pawley fit.

A single status bar at the foot of the results panel reports whatever is
running, from any tab. One report button, top right between the theme toggle
and the help icon, exports whichever tab is in front as a PDF &mdash; including
the plot before you have loaded any data, which gives you the theoretical
reflection list.

---

## Help

The [help file](powder5_help.html) is the real documentation: the methods, what
every control does, what the numbers mean, and where each one stops being
trustworthy. It opens from the **?** button in the application.

Two things worth knowing before you start:

- **Peak overlap sets the limit, not the algorithm.** Where two reflections are
  exactly coincident the data determines only their sum. powder5 says so
  explicitly rather than printing a split it cannot support — see *Exactly
  Coincident Reflections* in the report and the *Overlap clusters* section of
  the help.
- **A structure solution from powder data is a hypothesis.** Judge it on the
  contacts, the bond-valence sums and the chemistry, not on the correlation
  coefficient.

---

## Requirements

Any current browser. WebGPU accelerates charge flipping and the Wyckoff search
where it is available (Chrome and Edge today); both fall back to the CPU
otherwise. Everything runs locally — no upload, no account, no network needed
after the page has loaded.

---

## Citing and credits

Developed by Nita Dragoe, Université Paris-Saclay. It is the successor to
PowderV2 — Dragoe, N. (2001), *J. Appl. Cryst.* **34**, 535.

Methods implemented here are due to their authors; the help file lists the
references in full. The principal ones are Le Bail *et al.* (1988), Pawley
(1981), Thompson, Cox & Hastings (1987), Stephens (1999), Oszlányi & Sütő
(2004) for charge flipping, and Larson & Von Dreele (2004) for the GSAS profile
parameterisation.

> **Disclaimer.** Provided for research and teaching. The algorithms are
> standard and the implementation is tested, but this is not a substitute for
> GSAS-II, FullProf or TOPAS for work intended for publication. Check anything
> that matters against one of them.

---

## License

[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International](http://creativecommons.org/licenses/by-nc-nd/4.0/)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">
  <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" />
</a>
