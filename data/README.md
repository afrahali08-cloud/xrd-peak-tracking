# Public experimental diffraction data

The included archive contains all 342 JSON patterns under `sutter-fella_heymans_0/` in the opXRD v4 in-situ archive, across 16 source groups. Each JSON member retains its original bytes, labels and metadata. The original records identify these patterns as experimental (`is_simulated=false`) and attribute them to Sutter Fella Lab, Lawrence Berkeley National Laboratory.

## Attribution and license

Hollarek, Daniel, and Friederich, Pascal. **opXRD: Open Experimental Powder X-ray Diffraction Database**, version v4, December 5, 2024. Zenodo. DOI [10.5281/zenodo.14279434](https://doi.org/10.5281/zenodo.14279434).

The selected contribution is associated with Gregoire Heymans and Carolin M. Sutter-Fella in the dataset contributor list. Embedded pattern metadata names Carolin Sutter-Fella as contributor and Sutter Fella Lab, LBNL as institution.

Dataset license is [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). The data's license applies to the redistributed subset. Selection and repackaging were performed for this repository. Individual JSON content is unchanged. Analysis outputs derive from these data after the intensity scaling documented below. No endorsement by the source authors is implied.

## Reproducibility

- Pinned source file is `opxrd_in_situ.zip`, 373,157,285 bytes, from the versioned DOI above.
- Source MD5 is `0ab26d52fedde4dae422ddcce8625df9`, verified against the publisher's record.
- Included subset SHA-256 is recorded in `opxrd-heymans.manifest.json`, along with the original archive paths and each member's SHA-256.
- `python public_data.py --download` downloads the pinned source and deterministically rebuilds the subset, then analyzes it.
- `python public_data.py --archive /path/to/opxrd_in_situ.zip` rebuilds from a local copy.
- `python public_data.py` runs offline using the included subset.

## Selection and interpretation

This is a complete contributor subset with increasing, uniform angle grids, 4,501 samples per pattern from 15° to 60° at 0.01° spacing. It was selected for compatibility with the detector's input contract. Several other contributor examples in this older archive had nonmonotonic angle arrays, so they were excluded rather than repaired speculatively. The chosen source version has a separate, moderately sized in-situ archive and is pinned for repeatability.

Patterns are grouped by the source filename's group number and ordered numerically by its final index. Tracker state resets between groups. These indices are an ordering convention in the demonstration. Physical elapsed time, continuous acquisition and operating conditions are unverified from the available metadata, so the result is described as an in-situ archive demonstration rather than a validated operando time-series benchmark.

The original intensities remain in the included JSON files. For processing, each group's intensities are multiplied by 100 divided by that group's maximum intensity. This makes the detector's fixed prominence threshold interpretable on a common scale while preserving relative amplitudes within the group. Scaling uses the whole group, so this preprocessing is offline. No angle interpolation or invented scan frames are used.

There are no hand-labeled peak centers or track identities in this subset. Phase-label fields are incomplete and are not used as peak truth. Detection counts, associations and runtime are reported. Precision, recall and identity accuracy require additional reference labels.
