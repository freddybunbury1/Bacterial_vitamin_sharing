# Compiled figures

Figure S2 uses the full split-rooted tree of 231 QC-retained genomes and its
GTDB taxonomy, without requiring measured B12 phenotypes. Figure 1A retains
the complete-trait subset.

Figure S11 compares initial source-culture OD600 (before B12 addition) with
four-hour OD600 in the 0 mM glucose, 500 pM added-B12 uptake wells. Cultures are
paired by strain, experiment, and biological replicate before averaging;
two-sided paired Wilcoxon tests use strain means within each trait group.
The former Figures S11–S18 are now S12–S19.

Figure S15C includes a pre-culture adjusted OD600 group before its two four-hour
glucose conditions, using the same strains and averaging the initial source
cultures per strain. Lines and paired tests compare only the four-hour 0 and
10 mM C groups. Figure S15 column widths are 32%, 32%, and 36%.

The final figure-assembly notebook and its required helper assets belong here.
Panel letters are bold, including labels inside wrapped plots and the raster
Figure 1 assembly. Figure 3 gives both plotting areas 5% more width within
its existing canvas. Figure 4 transfers 4% of the left-column allocation to
the right column and lowers only panel B's y-axis title by 4% of the plot
height. Figure 5 likewise transfers 4% of panel A's width allocation to B.
Export dimensions and text sizes are unchanged by these layout adjustments.
Figure S1 places the sampling-size legend below the map and reduces panel B's
plotting height by 30%, keeping its width and centring it vertically.

Figure S9 combines raw survey OD600 versus heat-killed DCS1 fluorescence
(panel A, with marginal densities and blanks in red) and DCS1 fluorescence
increase above each strain/experiment/biological-replicate 0%-killed baseline
versus killed proportion times blank-adjusted source-culture OD600 at t = 0
(panel B, coloured by trait group). Panel B excludes blanks and nonpositive
axis values, retaining 642 measurements from 23 strains; the unfiltered
differences are retained in `figure_s9b_all_data`. Both panels use logarithmic axes. The
figure uses 10-point base text and 11-point panel labels and is saved at
4.97 × 6.448 inches; panel B's legend is inside the top-left corner, and its
plotting width is reduced by 8% while retaining its left edge. Its inputs are generated
within the release pipeline, with no dependency on exploratory files.
Figures S12–S14 compare first-order (red) and quota (blue) structures,
with line alpha 0.8, shared facet scales and the blank control. Figure S13
retains the continuous y-axis upper limit of 1.5. Figure S16 shows the GTDB phylogeny of
the 23 selected release/uptake strains alongside broad environment, survey
total B12, uptake, and EC-B12 fraction, trait group, the glucose-zero
first-order release and uptake rates, strain ID, and the established
trait-group ordinal.
Figures 4B, 5B–C and S16 show OD-specific coefficients, normalized per
replicate to mean initial/final blank-adjusted OD600 before geometric
averaging. The additional Figure 5 uptake panels use the same normalization.
Figure 5B statistics are recomputed on normalized rates; Figure 5D's ratio
is unchanged because both coefficients use the same OD denominator.
Figures S17–S19 show candidate-gene associations, bioassay comparison, and
survey fluorescence distributions. The same Table S3 strain ordinals are
used in Figures S12–S14 and S16. Manuscript captions are maintained separately in Overleaf.

The shared genome-linked phenotype cohort is defined in
`R/utils/b12_figure_cohort.R`: QC-passing study-tree membership plus finite
synthesis and uptake, with no minimum sample count or uncertainty requirement.
Figures 1B and S5A use 202 strains; synthesizer-only panels use 79
(`total_b12_gm > 50`). Figures 2 and S6 use the same cohort after canonicalizing
the KC identifier. Missing error bars do not remove valid points.
