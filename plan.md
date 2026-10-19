## Plan: Gaia Open-Cluster Ages

Build a reproducible command-line Python workflow that queries Gaia DR3 `gaia_source_lite` around the Pleiades, Praesepe, and M67, identifies likely members from the joint `(pmra, pmdec, parallax)` distribution, produces Gaia HR diagrams, estimates relative ages from the turn-off, and performs a quantitative isochrone age fit with uncertainties.

**Steps**

1. **Define the scientific contract and configuration.**
   - Create a YAML/TOML configuration containing the three cluster names, ICRS cone centers, initial radii, Gaia table name (`gaiadr3.gaia_source_lite`), ADQL column list, quality thresholds, clustering settings, isochrone source, parameter grids, random seed, and output directory.
   - Use configurable starting cones rather than hard-coded assumptions; sensible initial radii are about 5 degrees for the Pleiades and Praesepe and 2 degrees for M67, with the final values checked against source counts and field contamination.
   - Record units, Gaia reference epoch, model version, query timestamp, configuration hash, and literature priors used for metallicity/extinction.

2. **Set up the project and dependencies.**
   - Create a small Python package/CLI with `pyproject.toml`, a README, configuration example, and separate modules for querying, preprocessing, membership, plotting, isochrones, and reporting.
   - Use `astroquery.gaia`, `astropy`, `numpy`, `pandas` or Astropy Tables, `scipy`, `scikit-learn`, `matplotlib`, and a structured configuration/CLI library.
   - Prefer an isochrone adapter interface. Start with a Gaia-passband PARSEC grid if a maintained grid is available; otherwise use a documented MIST grid and preserve the adapter so the model family can be changed without rewriting the fitter. Do not assume the `isochrones` Python package provides PARSEC, since its documented built-in models are MIST-oriented.
   - Pin compatible versions and provide an environment/installation check.

3. **Implement Gaia table discovery and cone queries.**
   - Validate that the requested `gaiadr3.gaia_source_lite` table exists and that every required column is present before querying. Required data should include `source_id`, `ra`, `dec`, `parallax`, `parallax_error`, `pmra`, `pmra_error`, `pmdec`, `pmdec_error`, `phot_g_mean_mag`, `phot_bp_mean_mag`, `phot_rp_mean_mag`, and available quality/observation fields such as `ruwe`, flux-excess information, and visibility/transit counts.
   - Build an ADQL cone query using `CONTAINS(POINT(...), CIRCLE(...)) = 1` or the `Gaia.cone_search_async` API, and use asynchronous jobs for reproducible larger results. Never use `SELECT *`; select and document the exact fields.
   - Set the row limit explicitly, handle archive timeouts/retries, and save raw results locally as VOTable/FITS plus the ADQL text, job metadata, and a checksum.
   - On later runs, reuse the cached result when the cluster, cone, table, column list, and query version match; expose a `--refresh` option.
   - Do not silently substitute `gaiadr3.gaia_source` if `gaia_source_lite` is unavailable. Either fail with an actionable message or use an explicit configuration fallback.

4. **Normalize and quality-filter each catalogue.**
   - Convert all columns to numeric arrays with units; remove rows missing astrometry or any of the three Gaia magnitudes.
   - Apply documented astrometric/photometric quality filters: finite five-parameter astrometry, reasonable parallax uncertainty/significance, positive photometric fluxes, sensible BP/RP flux-excess limits, and a configurable RUWE/visibility/transit criterion. Keep a less restrictive diagnostic sample so filtering effects can be inspected.
   - Compute `bp_rp = phot_bp_mean_mag - phot_rp_mean_mag`, magnitude uncertainties from flux uncertainties where available, and an astrometric quality summary.
   - Apply the Gaia DR3 parallax zero-point correction using a versioned implementation and record both raw and corrected parallaxes. Treat the correction as a systematic uncertainty rather than implying it removes all distance systematics.
   - Preserve every source and add filter flags instead of deleting information irreversibly.

5. **Make diagnostic plots before membership selection.**
   - For each cluster, plot sky positions, `pmra` versus `pmdec`, parallax histograms, proper-motion histograms, and a preliminary `pmra/pmdec/parallax` scatter view.
   - Mark quality-rejected sources and report the number of queried, valid, and retained sources. These diagnostics are the cheap check that the cone contains a visible overdensity and that the selected table/columns are behaving as expected.

6. **Cluster the astrometry in three dimensions.**
   - Build the feature matrix from corrected parallax, `pmra`, and `pmdec`; use robust scaling based on the field sample so parallax and proper-motion units cannot dominate by accident.
   - Use HDBSCAN as the preferred density-based method if available, with DBSCAN as a deterministic fallback. Configure minimum cluster size, neighborhood scale, and random seed where relevant.
   - Select the cluster component using a combination of spatial concentration, high parallax/proper-motion density, and proximity to the expected cluster locus, rather than assuming label zero is the cluster.
   - Estimate the cluster center and intrinsic scatter robustly, then run an iterative refinement. Where uncertainties and Gaia covariance terms are available, compute a measurement-error-aware Mahalanobis membership distance; otherwise retain the standardized clustering score and label it as approximate.
   - Produce continuous membership probabilities or scores, not only a hard label. Store `member_score`, `member_flag`, cluster-locus residuals, and the model parameters for every source.
   - Add optional validation against known published members or an external membership catalogue if supplied, but keep this outside the primary selection so the method remains Gaia-based.

7. **Build the member catalogue and estimate common distance.**
   - Use high-confidence astrometric members for an initial cluster parallax location, robustly fit a common parallax while accounting for measurement errors and the corrected zero point, and report the intrinsic depth/systematic floor.
   - Use the fitted cluster distance/parallax for the HR diagram rather than independently inverting noisy individual parallaxes. Retain individual parallax-based absolute magnitudes only as a comparison diagnostic for high-S/N sources.
   - Save a member table containing source IDs, coordinates, astrometry, errors, photometry, quality flags, membership score, cluster distance quantities, and provenance.

8. **Create extinction-corrected Gaia HR diagrams.**
   - Calculate absolute `G` magnitude from the cluster distance modulus and Gaia `BP-RP` colour. Apply cluster-specific extinction coefficients consistently to `G`, BP, and RP, with `A_V`/`A_G` and colour excess either from literature priors or a configured fit range.
   - Produce, for each cluster, a two-panel HR diagram: all quality-retained cone sources in the background and members colour-coded by membership score, plus a clean member-only panel with photometric error bars where useful.
   - Overlay the best-fit isochrone, turn-off marker, uncertainty envelope, and labels for likely binaries, outliers, and evolved stars where the data support them. Invert the magnitude axis and document the adopted photometric system/model passbands.

9. **Estimate relative ages from the turn-off.**
   - Define the turn-off operationally: the hottest/brightest termination or bend of the member main sequence, using a reproducible colour/magnitude envelope or the turn-off mass returned by the isochrone grid, not a hand-picked point.
   - Fit a smooth main-sequence ridge/envelope as a diagnostic and identify the turn-off region with bootstrap resampling of members and perturbation by photometric/distance/extinction uncertainties.
   - Compare turn-off luminosity/colour and fitted turn-off mass across the three clusters. Report the ordering with confidence intervals and state the expected qualitative result: Pleiades youngest, Praesepe intermediate, M67 oldest, subject to model/data uncertainty.
   - Keep this as a relative-age diagnostic; do not present turn-off ordering alone as the final age measurement.

10. **Implement quantitative isochrone fitting.**
    - Define the model parameters per cluster: `log10(age/yr)`, distance modulus or cluster parallax, extinction (`A_V` or `A_G`), and `[Fe/H]`; use literature-informed priors/ranges for metallicity and extinction, with an option to fix them. Include an unresolved-binary broadening/mixture term or a documented binary exclusion strategy.
    - Interpolate the chosen isochrone grid into Gaia `G`, `BP`, and `RP` passbands. Apply distance and extinction transformations, and compare observed member colours/magnitudes to the nearest point on the model sequence, including photometric errors, distance uncertainty, intrinsic sequence width, binaries, and membership weights.
    - Start with a deterministic grid likelihood over age, distance, extinction, and metallicity. Normalize likelihoods, marginalize nuisance parameters, and return the maximum-likelihood/maximum-posterior age plus median and credible intervals.
    - Add an optional MCMC refinement only after the grid fitter is validated; use it to sample posterior correlations, especially age-extinction-distance-metallicity degeneracy. Keep this optional to avoid making the baseline workflow unnecessarily fragile.
    - Fit using high-confidence members first, then repeat with probability-weighted members and sensitivity variants. Compare ages from both to expose contamination dependence.
    - Save the complete likelihood/posterior grid, best-fit parameters, uncertainty summaries, residuals, and model metadata.

11. **Quantify systematics and robustness.**
    - Repeat fits under plausible choices of cone radius, quality cuts, cluster-membership threshold, parallax zero-point treatment, extinction prior, metallicity prior, binary model, and isochrone family when available.
    - Bootstrap members and perturb measurements to estimate sampling uncertainty. Separate statistical uncertainty from model/extinction/membership systematic ranges.
    - Check for differential extinction, unresolved binaries, blue stragglers, pre-main-sequence effects in the Pleiades, saturation/bright-star issues, faint incompleteness, and field stars that share similar astrometry.
    - Compare fitted distances, metallicities, and ages against independent literature values and clearly label disagreement as a diagnostic, not an automatic correction.

12. **Produce outputs and a final report.**
    - Write per-cluster raw and processed catalogues, member catalogues, fit grids/results in CSV/JSON/FITS, and publication-quality PNG/PDF figures.
    - Generate one summary table with fitted age, uncertainty, distance, extinction, metallicity, number of members, turn-off estimate, and relative-age rank.
    - Generate a markdown/HTML report containing query provenance, filters, diagnostic plots, membership results, HR diagrams, turn-off comparison, isochrone overlays, fit residuals, sensitivity tests, caveats, and Gaia/model citations.
    - Provide CLI stages such as `query`, `members`, `plot`, `fit`, `report`, and `all`, with a `--offline` mode that uses cached catalogues and model grids.

13. **Test and validate the workflow.**
    - Unit-test ADQL construction, table-column validation, unit conversions, parallax correction plumbing, absolute-magnitude/extinction calculations, cluster-label selection, likelihood normalization, and turn-off extraction using small synthetic fixtures.
    - Add an integration test with a saved miniature Gaia-like catalogue so tests do not depend on the live archive.
    - Run a live smoke test for one small cone/query, verify cached reruns are identical, and verify all three clusters produce non-empty member samples and nonblank figures.
    - Check expected sanity constraints: cluster members form compact astrometric loci, the member HR diagrams show coherent sequences, the fitted distances are plausible, and the age ranking is stable under reasonable membership thresholds.
    - Record failures caused by unavailable archive/model downloads separately from code failures.

**Relevant files**

- `/home/tracer/TEST_CHAMBER2/agent_testing/pyproject.toml` - dependencies, CLI entry point, and test tooling.
- `/home/tracer/TEST_CHAMBER2/agent_testing/config/clusters.yaml` - cluster coordinates, cone radii, quality thresholds, priors, and fitting grids.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/query.py` - Gaia table validation, ADQL/cone queries, retries, and cache metadata.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/preprocess.py` - units, quality flags, flux-excess/RUWE filtering, colours, and parallax-zero-point handling.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/membership.py` - scaled 3D clustering, robust locus refinement, membership scores, and diagnostics.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/photometry.py` - distance modulus, extinction transformations, absolute magnitudes, and HR-diagram data.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/isochrones.py` - model-grid adapter, interpolation, binary treatment, likelihood, grid search, and optional posterior refinement.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/plots.py` - astrometric diagnostics, HR diagrams, turn-off overlays, and fit residual plots.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/turnoff.py` - turn-off definition, ridge/envelope extraction, bootstrap uncertainty, and age ordering.
- `/home/tracer/TEST_CHAMBER2/agent_testing/src/gaia_clusters/cli.py` - staged commands and end-to-end orchestration.
- `/home/tracer/TEST_CHAMBER2/agent_testing/tests/` - unit, synthetic integration, and cache/reproducibility tests.
- `/home/tracer/TEST_CHAMBER2/agent_testing/README.md` - installation, live/offline usage, scientific assumptions, citations, and interpretation limits.
- `/home/tracer/TEST_CHAMBER2/agent_testing/data/` and `/home/tracer/TEST_CHAMBER2/agent_testing/outputs/` - git-ignored caches, downloaded grids, tables, figures, and reports.

**Verification**

1. Validate the exact Gaia DR3 table name and required columns before any production query; fail clearly if `gaia_source_lite` does not expose a required field.
2. Run a one-cluster live cone query with an explicit row limit, inspect row counts and required columns, and save the raw result plus ADQL provenance.
3. Run the same query in offline/cache mode and compare source IDs and derived columns for deterministic equality.
4. Inspect proper-motion/parallax diagnostics to confirm a compact cluster overdensity exists before accepting the clustering configuration.
5. Confirm member loci are compact in all three clustering dimensions and that the member-only HR diagram is more coherent than the full cone diagram.
6. Fit synthetic data generated from a known isochrone and verify recovery of age, distance, extinction, and metallicity within expected uncertainties.
7. Run the three-cluster end-to-end workflow and verify every requested figure, table, fit result, and report section is produced.
8. Perform sensitivity runs and require the qualitative age order to remain stable, while reporting any cluster whose absolute age is strongly model- or prior-dependent.

**Decisions**

- Deliverable: command-line Python workflow with a reusable configuration and optional notebook/report layer, rather than a one-off notebook.
- Primary fit: grid likelihood with nuisance-parameter marginalization; optional MCMC is a later enhancement.
- Isochrones: use a Gaia-passband PARSEC-compatible grid when available, with an adapter and a documented MIST fallback; do not rely on undocumented package support.
- Extinction and metallicity: use literature-informed cluster priors and fit them within bounded ranges by default, with fixed-value options for controlled comparisons.
- Membership: astrometry drives the initial selection, while photometry validates and refines the sample; the HR diagram must not be used to hide a failed astrometric selection.
- Scope excludes radial-velocity membership as a required input, detailed dynamical modeling, spectroscopic abundance reanalysis, and a claim of publication-grade ages without systematic validation.

**Further Considerations**

1. `gaia_source_lite` must be checked at implementation time because archive table names/column coverage can differ by service schema; keep the requested table as the default and make any fallback explicit.
2. Absolute cluster ages are degenerate with extinction, metallicity, binaries, model physics, and membership contamination. The report should present posterior intervals plus a systematic range, not a single over-precise number.
3. Gaia DR3 astrometry carries a global parallax bias and cluster-specific systematics; preserve raw and corrected values and cite the adopted correction/version.
