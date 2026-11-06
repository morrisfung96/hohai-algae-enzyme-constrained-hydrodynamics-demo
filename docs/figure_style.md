# Project Figure Style Memory

This file is a persistent project-level plotting rule for all new and revised
figures under `<project-root>`.

## Required Axis Frame Style

Use a full rectangular axis frame for every quantitative data panel.

- Keep the top, right, bottom and left spines visible.
- Do not use the half-frame style that hides the top and right spines.
- Apply this rule to line plots, scatter plots, bar plots, heatmaps, maps,
  stacked plots, error-bar plots and validation panels.
- If a panel has a secondary y-axis, the primary axis should still draw the
  full box; the secondary axis may draw only the right-side ticks/spine as
  needed.
- Empty layout placeholders, image-only molecular renderings and schematic
  panels may hide axes when the axes are not data coordinates.

## Matplotlib Pattern

Use this style pattern in new plotting scripts:

```python
def style_full_frame(ax, grid=True):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#9AA7A9")
        spine.set_linewidth(1.1)
    ax.tick_params(labelsize=8, colors="#142A38", direction="out", length=3)
    if grid:
        ax.grid(True, color="#D9E4E2", lw=0.5, alpha=0.75, zorder=0)
```

Avoid patterns like:

```python
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
```

## Visual Reference

The reference style is the user's supplied composite figure: each data subplot
is enclosed by a clean full rectangular border. The desired effect is a
publication-style framed panel, not minimalist half axes.

## Manuscript Asset Naming

Development iterations may use internal counters in working directories, but
submission-facing assets must not expose exploratory names such as `v20`,
`clean`, `highinfo`, `final`, or `new`.

- Use a stable public filename: `FigureN_<scientific-topic>.png/.pdf/.svg`.
- Store public figures under `docs/manuscript_figures/FigureN/`.
- Keep revision identity in `FigureN_audit.json` using restrained labels such
  as `manuscript-revA`; do not put that label on the figure canvas.
- Never call an asset `final`; the stable filename is replaced only after the
  new rendering passes its visual and scientific audit.
- Preserve old development directories as provenance unless the user
  explicitly asks for cleanup.

## Cross-Figure Typography and Panel Grammar

Figure 3 is the editorial reference for Figures 3 and 4.

- Use Arial/Helvetica with DejaVu Sans fallback and editable text in SVG/PDF.
- Keep one hierarchy: panel letters 8 pt bold; panel titles 7.0--7.4 pt bold;
  axis labels 6.5--7.2 pt; tick labels 5.8--6.4 pt; restrained annotations
  5.2--6.0 pt at final two-column size.
- Place panel letters and titles outside the data frame on a shared baseline.
- Use upright lowercase panel letters without circles or colored boxes.
- Use the restrained Figure 3 palette as the default: near-black observations,
  teal primary model/product, vermilion comparison, blue-grey uncertainty,
  pale grey grids, and a very pale teal validation background.
- Use direct labels when there are two series. Avoid a redundant legend plus
  direct labels.
- Keep method details, sample-lineage caveats, and interpretive footnotes in
  the caption or audit file. Only information required to read the marks may
  remain inside the panel.
- Use sentence case and concise scientific titles. Do not put development
  terminology or version identifiers on the canvas.

## Cross-Figure Consistency Feedback (2026-10-12)

- The user found the standalone mist-blue/coral/lavender Figure 4 treatment
  inconsistent with Figures 1–3. Inspect the complete figure family before
  changing a figure's palette; a locally attractive panel is not sufficient.
- Figures 1 and 2 are user-edited PPT artwork. Use them as visual references
  without silently replacing or restyling them.
- The quantitative-panel reference is the retained Figure 3 v18 artwork:
  `docs/est_figure3_component_contributions_v18/figures/Figure3_biomass_component_contributions_v18.svg`.
  This is a style reference only, not authorization to restore historical
  metric definitions or observation mappings.
- Existing reference colors are teal `#007F78`, orange `#D55E00`, text
  `#23343B`, frame `#7E8C90`, and secondary blue-grey `#5F8FAE`.
  At 7.5-inch figure width, the reference full-frame stroke is about 0.72 pt.
  Related shades may distinguish categories, with their meaning stated in
  the legend; continuous and categorical variables need distinct scales.
- The compact paired concentration/state layout remains a candidate. Check
  it alongside the first three figures; do not describe a new color variant
  as user-approved before the user has reviewed it.

## Orange Saturation Feedback (2026-10-12)

- The user requested a slightly less saturated orange family in Figure 4.
  The current map/legend trial reduces the two orange colors' HLS saturation
  by 30%, retaining hue and lightness before RGB rounding: `#D55E00` becomes
  `#B56220`, and `#E8A673` becomes `#D6A885`.
- This preference applies to the revised ecological state maps and their
  legend. It does not authorize changing the other manuscript figures.

## Figure 4a Resolution and Palette Review (2026-10-12)

- The user identified blurred Figure 4a in the composed trials. The confirmed
  source PDF and SVG embed panel a at only 199 x 168 pixels, while the exact
  matching original PNG (4500 x 4200 pixels) preserves sharp artwork.
- Restore the panel-a image from that PNG at its original placement before
  reusing the PDF/SVG: the retained crop is 1194 x 1008 pixels, approximately
  600 dpi. Increasing the export dpi of the old PDF does not fix this defect.
  Restoration assets and the source/crop audit are retained under
  `outputs/Figure4_palette_studies/restored_source/` and
  `outputs/Figure4_palette_studies/panel_a_resolution_audit.json`.
- The user still disliked the muted orange trial and requested multiple
  alternatives. Four review candidates are under `outputs/Figure4_palette_studies/`:
  A teal/apricot, B blue-grey/rose, C sage/slate, and D teal/mauve.
  All four retain the restored panel a and identical quantitative data.
  The user subsequently selected A ("使用 图A吧"); B, C and D remain
  unselected alternatives.

## Selected Figure 4 Palette (2026-10-12)

- Current Figure 4 uses the user-selected A palette: teal + apricot.
  Net gain / other losses: `#A4C7C2`; net gain / turnover: `#398B85`;
  net loss / turnover: `#D59870`; net loss / other losses: `#E9CDB4`.
- Keep the restored 1194 x 1008 pixel original panel-a artwork and the
  existing blue-grey continuous concentration scale.
- The selected PNG/PDF/SVG are published at
  `docs/manuscript_figures/Figure4/Figure4_concentration_ecological_states.*`.
  The audit records the explicit visual selection; it is not a new
  scientific validation or permission to change model states.
- Earlier muted-orange and palette-study outputs remain as provenance.
