"""Site-determining ion analysis for phosphosite localization pairs.
Analyzes site-determining ions for phosphosite localization pairs to distinguish between different modification sites.

This module provides an end-to-end workflow:
1. Read MGF spectra indexed by TITLE.
2. Parse two localization peptide strings with configurable modification masses.
3. Compute b/y fragment ions and select site-determining ions.
4. Match theoretical ions to observed peaks in the paired spectrum.
5. Run batch analysis from CSV with configurable column presets.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from tqdm import tqdm


# Preset schemas for quickly switching between different peptide-column combinations.
SPECTRUM_TITLE_COLUMN = "Title"
METHOD_COLUMN = "Method"
COLUMN_MAPPING_PRESETS: Dict[str, Dict[str, object]] = {
	"2vs4": {
		"method_a": {
			"name": "Method2",
			"peptide_col": "IsoformSequence_autort_pDeep",
		},
		"method_b": {
			"name": "Method4",
			"peptide_col": "IsoformSequence_autort_pDeep_phosSight",
		},
	},
	"1vs4": {
		"method_a": {
			"name": "Method1",
			"peptide_col": "IsoformSequence_PhosphoRS",
		},
		"method_b": {
			"name": "Method4",
			"peptide_col": "IsoformSequence_autort_pDeep_phosSight",
		},
	},
	"6vs8": {
		"method_a": {
			"name": "Method6",
			"peptide_col": "IsoformSequence_autort_pDeep",
		},
		"method_b": {
			"name": "Method8",
			"peptide_col": "IsoformSequence_autort_pDeep_phosSight",
		},
	},
}

# Monoisotopic masses.
PROTON_MASS = 1.00727646677
HYDROGEN_MASS = 1.007825032

# CID/HCD fragment ion terminal mass adjustments.
N_TERM_MASS = HYDROGEN_MASS
C_TERM_MASS = 17.002739654

# Modification delta masses.
PHOSPHO_DELTA_MASS = 79.966331
OXIDATION_DELTA_MASS = 15.994915

# Standard residue masses (neutral residue masses used in peptide mass formulas).
AMINO_ACID_MASSES: Dict[str, float] = {
	"A": 71.037114,
	"C": 103.009185,
	"D": 115.026943,
	"E": 129.042593,
	"F": 147.068414,
	"G": 57.021464,
	"H": 137.058912,
	"I": 113.084064,
	"K": 128.094963,
	"L": 113.084064,
	"M": 131.040485,
	"N": 114.042927,
	"P": 97.052764,
	"Q": 128.058578,
	"R": 156.101111,
	"S": 87.032028,
	"T": 101.047679,
	"V": 99.068414,
	"W": 186.079313,
	"Y": 163.063329
}


# Default modification masses. You can extend/override this dictionary.
MODIFICATION_MASSES: Dict[str, float] = {
	"1": AMINO_ACID_MASSES["M"] + OXIDATION_DELTA_MASS,
	"2": AMINO_ACID_MASSES["S"] + PHOSPHO_DELTA_MASS,
	"3": AMINO_ACID_MASSES["T"] + PHOSPHO_DELTA_MASS,
	"4": AMINO_ACID_MASSES["Y"] + PHOSPHO_DELTA_MASS,
}

@dataclass(frozen=True)
class Peak:
	mz: float
	intensity: float


@dataclass(frozen=True)
class Spectrum:
	title: str
	charge: int
	peaks: Tuple[Peak, ...]


@dataclass(frozen=True)
class FragmentIon:
	ion_type: str
	index: int
	charge: int
	mz: float

	@property
	def label(self) -> str:
		return f"{self.ion_type}{self.index}^{self.charge}+"


@dataclass(frozen=True)
class MatchedIon:
	ion: FragmentIon
	observed_mz: Optional[float]
	observed_intensity: float
	mass_error: Optional[float]


@dataclass(frozen=True)
class RowAnalysisResult:
	row_index: int
	spectrum_title: str
	peptide_a: str
	peptide_b: str
	precursor_charge: int
	site_determining_count_a: int
	site_determining_count_b: int
	site_determining_ions_a: Tuple[FragmentIon, ...]
	site_determining_ions_b: Tuple[FragmentIon, ...]


@dataclass(frozen=True)
class IndexedMethodRow:
	row_index: int
	spectrum_title: str
	method: str
	peptide: str


@dataclass(frozen=True)
class UnmatchedPairRecord:
	spectrum_title: str
	present_method: str
	missing_method: str
	present_row_index: int
	peptide_a: str
	peptide_b: str


def parse_mgf(mgf_path: str | Path) -> Dict[str, Spectrum]:
	"""Parse MGF and return {TITLE: Spectrum}.

	The parser expects standard MGF blocks enclosed by BEGIN IONS / END IONS.
	"""

	spectra: Dict[str, Spectrum] = {}
	in_block = False
	meta: Dict[str, str] = {}
	peaks: List[Peak] = []

	def finalize_block() -> None:
		title = meta.get("TITLE", "").strip()
		if not title:
			raise ValueError("MGF spectrum block is missing TITLE.")
		if title in spectra:
			raise ValueError(f"Duplicate TITLE '{title}' found in MGF file '{mgf_path}'.")
		charge_text = meta.get("CHARGE")
		if not charge_text or not charge_text.strip():
			raise ValueError(f"MGF spectrum '{title}' is missing CHARGE.")
		charge = _parse_charge(charge_text)
		if not peaks:
			raise ValueError(f"MGF spectrum '{title}' has no peak list.")
		spectra[title] = Spectrum(title=title, charge=charge, peaks=tuple(peaks))

	with Path(mgf_path).open("r", encoding="utf-8") as handle:
		for raw_line in tqdm(
			handle,
			desc=f"Parsing spectra from {Path(mgf_path).name}",
			unit="line",
			leave=False,
		):
			line = raw_line.strip()
			if not line:
				continue
			if line == "BEGIN IONS":
				in_block = True
				meta = {}
				peaks = []
				continue
			if line == "END IONS":
				if in_block:
					finalize_block()
				in_block = False
				meta = {}
				peaks = []
				continue
			if not in_block:
				continue

			if "=" in line:
				key, value = line.split("=", 1)
				meta[key.strip().upper()] = value.strip()
				continue

			parts = line.split()
			if len(parts) == 2:
				try:
					peaks.append(Peak(mz=float(parts[0]), intensity=float(parts[1])))
				except ValueError:
					continue

	return spectra


def collect_mgf_files(mgf_dir_path: str | Path) -> Tuple[Path, ...]:
	"""Collect non-recursive .mgf files from the provided directory."""

	mgf_dir = Path(mgf_dir_path)
	if not mgf_dir.exists():
		raise ValueError(f"MGF directory does not exist: '{mgf_dir_path}'.")
	if not mgf_dir.is_dir():
		raise ValueError(f"MGF directory path is not a directory: '{mgf_dir_path}'.")

	mgf_files = sorted(
		(
			path
			for path in mgf_dir.iterdir()
			if path.is_file() and path.suffix.lower() == ".mgf"
		),
		key=lambda p: (p.name.lower(), p.name),
	)
	if not mgf_files:
		raise ValueError(f"No .mgf files found in directory: '{mgf_dir_path}'.")

	return tuple(mgf_files)


def _parse_charge(charge_text: Optional[str]) -> int:
	if charge_text is None:
		raise ValueError("MGF spectrum is missing CHARGE.")
	text = charge_text.strip()
	match = re.fullmatch(r"(\d+)\+", text)
	if not match:
		raise ValueError(f"Invalid CHARGE value: {charge_text!r}. Expected a single integer followed by '+'.")
	return int(match.group(1))


def parse_modified_peptide(
	peptide: str,
) -> List[float]:
	"""Convert modified peptide notation into per-residue neutral masses.

	Supported notation examples:
	- "ACDEFG" (plain sequence)
	- "A1C" (digit token as modified amino acid, e.g. 1=oxidized M)
	- "A2C" / "A3C" / "A4C" (digit token as phospho S/T/Y)
	"""

	if not peptide:
		raise ValueError("Peptide string is empty.")

	masses: List[float] = []
	i = 0
	n = len(peptide)

	while i < n:
		char = peptide[i]

		if char.isdigit():
			if char not in MODIFICATION_MASSES:
				raise ValueError(
					f"Unknown modified-residue key '{char}' in peptide '{peptide}'."
				)
			masses.append(float(MODIFICATION_MASSES[char]))
			i += 1
			continue

		if char.isalpha():
			if char not in AMINO_ACID_MASSES:
				raise ValueError(f"Unknown amino acid '{char}' in peptide '{peptide}'.")

			mass = AMINO_ACID_MASSES[char]
			i += 1

			masses.append(mass)
			continue

		raise ValueError(f"Unsupported character '{char}' in peptide '{peptide}'.")

	if len(masses) < 2:
		raise ValueError(f"Peptide '{peptide}' must have at least two residues.")
	return masses


def determine_fragment_charges(precursor_charge: Optional[int]) -> Tuple[int, ...]:
	"""Return fragment charge states for a required precursor charge.

	The precursor charge must be present. Missing charge is treated as a data error
	instead of being silently replaced with a default value.

	- precursor charge == 2: use fragment z=1 only
	- precursor charge >= 3: use fragment z=1,2
	"""

	if precursor_charge is None:
		raise ValueError("Precursor charge is missing.")
	if precursor_charge >= 3:
		return (1, 2)
	return (1,)


def calculate_b_y_ions(
	residue_masses: Sequence[float],
	fragment_charges: Sequence[int],
) -> Dict[Tuple[str, int, int], FragmentIon]:
	"""Compute charged b/y ions.
	Generates theoretical charged b/y fragment ions based on given residue masses and charge states.

	Neutral formulas:
	- b: [N] + [M] - H
	- y: [C] + [M] + H
	"""

	ions: Dict[Tuple[str, int, int], FragmentIon] = {}
	n = len(residue_masses)

	prefix = 0.0
	for i in range(1, n):
		prefix += residue_masses[i - 1]
		neutral_b = N_TERM_MASS + prefix - HYDROGEN_MASS
		for z in fragment_charges:
			mz = (neutral_b + z * PROTON_MASS) / z
			ion = FragmentIon("b", i, z, mz)
			ions[("b", i, z)] = ion

	suffix = 0.0
	for i in range(1, n):
		suffix += residue_masses[n - i]
		neutral_y = C_TERM_MASS + suffix + HYDROGEN_MASS
		for z in fragment_charges:
			mz = (neutral_y + z * PROTON_MASS) / z
			ion = FragmentIon("y", i, z, mz)
			ions[("y", i, z)] = ion

	return ions


def extract_site_determining_ions(
	ions_a: Mapping[Tuple[str, int, int], FragmentIon],
	ions_b: Mapping[Tuple[str, int, int], FragmentIon],
	ms2_tolerance: float,
) -> Tuple[Dict[Tuple[str, int, int], FragmentIon], Dict[Tuple[str, int, int], FragmentIon]]:
	"""Keep site-determining ions between two localizations.
	Extracts signature fragment ions that can differentiate between peptide A and peptide B.

	A site-determining ion is either:
	- present on one side only (e.g., a peak exists in one peptide but not the other), or
	- present on both sides but with m/z difference greater than tolerance (due to different modification sites).
	"""

	keys = set(ions_a).union(ions_b)
	selected_a: Dict[Tuple[str, int, int], FragmentIon] = {}
	selected_b: Dict[Tuple[str, int, int], FragmentIon] = {}
	for key in keys:
		ion_a = ions_a.get(key)
		ion_b = ions_b.get(key)

		if ion_a is None and ion_b is not None:
			selected_b[key] = ion_b
			continue
		if ion_b is None and ion_a is not None:
			selected_a[key] = ion_a
			continue

		if ion_a is not None and ion_b is not None and abs(ion_a.mz - ion_b.mz) > ms2_tolerance:
			selected_a[key] = ions_a[key]
			selected_b[key] = ions_b[key]
	return selected_a, selected_b


def match_ions_to_spectrum(
	ions: Iterable[FragmentIon],
	spectrum: Spectrum,
	ms2_tolerance: float,
) -> Tuple[MatchedIon, ...]:
	"""Match ions to nearest peak within tolerance and record observed intensity.
	Matches theoretical ions with actual observed MS2 peaks and records intensity within ms2_tolerance.

	Only matched ions are returned; unmatched theoretical ions are omitted.
	Theoretical ions not matched in the actual spectrum are ignored.
	"""

	sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.mz)
	peak_mzs = [p.mz for p in sorted_peaks]
	matched: List[MatchedIon] = []

	for ion in sorted(ions, key=lambda x: (x.ion_type, x.index, x.charge)):
		best = _find_best_peak(ion.mz, sorted_peaks, peak_mzs, ms2_tolerance)
		if best is None:
			continue
		else:
			obs, err = best
			matched.append(
				MatchedIon(
					ion=ion,
					observed_mz=obs.mz,
					observed_intensity=obs.intensity,
					mass_error=err,
				)
			)
	return tuple(matched)


def _find_best_peak(
	target_mz: float,
	sorted_peaks: Sequence[Peak],
	peak_mzs: Sequence[float],
	ms2_tolerance: float,
) -> Optional[Tuple[Peak, float]]:
	import bisect

	left = bisect.bisect_left(peak_mzs, target_mz - ms2_tolerance)
	right = bisect.bisect_right(peak_mzs, target_mz + ms2_tolerance)
	if left >= right:
		return None

	window = sorted_peaks[left:right]
	best_peak = min(window, key=lambda p: abs(p.mz - target_mz))
	return best_peak, (best_peak.mz - target_mz)


def analyze_localization_pair(
	spectrum: Spectrum,
	peptide_a: str,
	peptide_b: str,
	ms2_tolerance: float,
) -> RowAnalysisResult:
	"""Analyze one spectrum + two localization peptides.
	Analyzes a single spectrum and its two predicted localization peptides (including generating
	theoretical ions, comparing differences, and matching against actual observed peaks).

	Only site-determining ions that can be matched to the spectrum within
	ms2_tolerance are retained in the returned result.
	"""

	fragment_charges = determine_fragment_charges(spectrum.charge)

	masses_a = parse_modified_peptide(peptide_a)
	masses_b = parse_modified_peptide(peptide_b)

	if len(masses_a) != len(masses_b):
		raise ValueError(
			"Peptide A and B must have the same length for localization comparison."
		)

	ions_a = calculate_b_y_ions(
		residue_masses=masses_a,
		fragment_charges=fragment_charges,
	)
	ions_b = calculate_b_y_ions(
		residue_masses=masses_b,
		fragment_charges=fragment_charges,
	)

	sd_a, sd_b = extract_site_determining_ions(ions_a, ions_b, ms2_tolerance)
	matched_sd_a = match_ions_to_spectrum(sd_a.values(), spectrum, ms2_tolerance)
	matched_sd_b = match_ions_to_spectrum(sd_b.values(), spectrum, ms2_tolerance)
	matched_ions_a = tuple(match.ion for match in matched_sd_a)
	matched_ions_b = tuple(match.ion for match in matched_sd_b)

	return RowAnalysisResult(
		row_index=-1,
		spectrum_title=spectrum.title,
		peptide_a=peptide_a,
		peptide_b=peptide_b,
		precursor_charge=spectrum.charge,
		site_determining_count_a=len(matched_ions_a),
		site_determining_count_b=len(matched_ions_b),
		site_determining_ions_a=tuple(sorted(matched_ions_a, key=lambda x: (x.ion_type, x.index, x.charge))),
		site_determining_ions_b=tuple(sorted(matched_ions_b, key=lambda x: (x.ion_type, x.index, x.charge))),
	)


def row_result_to_dict(result: RowAnalysisResult) -> Dict[str, object]:
	"""Convert a row result into a flat summary dictionary."""

	return {
		"row_index": result.row_index,
		"spectrum_title": result.spectrum_title,
		"peptide_a": result.peptide_a,
		"peptide_b": result.peptide_b,
		"precursor_charge": result.precursor_charge,
		"site_determining_count_a": result.site_determining_count_a,
		"site_determining_count_b": result.site_determining_count_b,
		"site_determining_ions_a": _serialize_fragment_ions(result.site_determining_ions_a),
		"site_determining_ions_b": _serialize_fragment_ions(result.site_determining_ions_b),
	}


def _serialize_fragment_ions(ions: Sequence[FragmentIon]) -> str:
	payload = [
		{
			"label": ion.label,
			"ion_type": ion.ion_type,
			"index": ion.index,
			"charge": ion.charge,
			"theoretical_mz": ion.mz,
		}
		for ion in ions
	]
	return json.dumps(payload, ensure_ascii=True)


def analyze_pairs_from_csv(
	csv_path: str | Path,
	mgf_dir_path: str | Path,
	method_a_name: str,
	method_b_name: str,
	peptide_a_col: str,
	peptide_b_col: str,
	ms2_tolerance: float,
	output_csv_path: Optional[str | Path] = None,
	output_unmatched_csv_path: Optional[str | Path] = None,
	output_matched_diff_csv_path: Optional[str | Path] = None,
	output_unmatched_recheck_csv_path: Optional[str | Path] = None,
) -> Tuple[List[RowAnalysisResult], List[UnmatchedPairRecord]]:
	"""Batch analyze peptide localization pairs from CSV.

	This function performs cross-method pairing:
	1) index rows for method A/B by spectrum title,
	2) pair by key for analysis,
	3) keep unmatched keys for downstream inspection.

	Fail-fast conditions:
	- duplicate spectrum title inside the same method side,
	- missing required title / peptide values on target method rows,
	- matched pair refers to a spectrum title that does not exist in MGF.
	"""

	mgf_files = collect_mgf_files(mgf_dir_path)
	spectra: Dict[str, Spectrum] = {}
	title_sources: Dict[str, Path] = {}
	mgf_iter: Iterable[Path] = tqdm(mgf_files, total=len(mgf_files), desc="Loading MGF files", unit="file")
	for mgf_file in mgf_iter:
		current = parse_mgf(mgf_file)
		for title, spectrum in current.items():
			if title in spectra:
				first_source = title_sources[title]
				raise ValueError(
					f"Duplicate TITLE '{title}' found across MGF files: "
					f"'{first_source}' and '{mgf_file}'."
				)
			spectra[title] = spectrum
			title_sources[title] = mgf_file

	target_methods = {method_a_name, method_b_name}
	if len(target_methods) != 2:
		raise ValueError("method_a and method_b must be different method names.")

	indexed_a: Dict[str, IndexedMethodRow] = {}
	indexed_b: Dict[str, IndexedMethodRow] = {}
	target_rows_by_index: Dict[int, Mapping[str, str]] = {}

	with Path(csv_path).open("r", encoding="utf-8", newline="") as handle:
		reader = csv.DictReader(handle)
		if not reader.fieldnames:
			raise ValueError("Input CSV has no header row.")
		header_set = set(reader.fieldnames)
		required_columns = {
			METHOD_COLUMN,
			SPECTRUM_TITLE_COLUMN,
			peptide_a_col,
			peptide_b_col,
		}
		missing_columns = sorted(column for column in required_columns if column not in header_set)
		if missing_columns:
			raise ValueError(f"Input CSV missing required columns: {', '.join(missing_columns)}")

		row_iter: Iterable[Mapping[str, str]] = tqdm(
			reader,
			desc="Analyzing CSV rows",
			unit="row",
		)
		for idx, row in enumerate(row_iter):
			if all((value is None) or (str(value).strip() == "") for value in row.values()):
				continue

			row_method = (row.get(METHOD_COLUMN) or "").strip()
			if row_method not in target_methods:
				continue

			title = (row.get(SPECTRUM_TITLE_COLUMN) or "").strip()

			if not title:
				raise ValueError(
					f"Row {idx}: Missing title for target method row. "
					f"method={row_method!r}, title={title!r}."
				)

			if row_method == method_a_name:
				peptide = (row.get(peptide_a_col) or "").strip()
				if not peptide:
					raise ValueError(
						f"Row {idx}: Missing peptide A value for method '{method_a_name}' "
						f"in column '{peptide_a_col}'."
					)
				target_index = indexed_a
			else:
				peptide = (row.get(peptide_b_col) or "").strip()
				if not peptide:
					raise ValueError(
						f"Row {idx}: Missing peptide B value for method '{method_b_name}' "
						f"in column '{peptide_b_col}'."
					)
				target_index = indexed_b

			if title in target_index:
				prev = target_index[title]
				raise ValueError(
					f"Duplicate title '{title}' within method '{row_method}'. "
					f"Rows {prev.row_index} and {idx}."
				)

			target_index[title] = IndexedMethodRow(
				row_index=idx,
				spectrum_title=title,
				method=row_method,
				peptide=peptide,
			)
			target_rows_by_index[idx] = row

	results: List[RowAnalysisResult] = []
	unmatched: List[UnmatchedPairRecord] = []
	all_titles = sorted(set(indexed_a).union(indexed_b))
	for spectrum_title in all_titles:
		entry_a = indexed_a.get(spectrum_title)
		entry_b = indexed_b.get(spectrum_title)

		if entry_a is None and entry_b is not None:
			unmatched.append(
				UnmatchedPairRecord(
					spectrum_title=entry_b.spectrum_title,
					present_method=entry_b.method,
					missing_method=method_a_name,
					present_row_index=entry_b.row_index,
					peptide_a="",
					peptide_b=entry_b.peptide,
				)
			)
			continue
		if entry_b is None and entry_a is not None:
			unmatched.append(
				UnmatchedPairRecord(
					spectrum_title=entry_a.spectrum_title,
					present_method=entry_a.method,
					missing_method=method_b_name,
					present_row_index=entry_a.row_index,
					peptide_a=entry_a.peptide,
					peptide_b="",
				)
			)
			continue

		if entry_a is None or entry_b is None:
			raise RuntimeError(
				f"Internal invariant violated for title '{spectrum_title}': "
				"title exists in union(indexed_a, indexed_b) but both sides are missing."
			)

		spectrum = spectra.get(spectrum_title)
		if spectrum is None:
			raise ValueError(
				f"Spectrum title '{spectrum_title}' not found in MGF."
			)

		result = analyze_localization_pair(
			spectrum=spectrum,
			peptide_a=entry_a.peptide,
			peptide_b=entry_b.peptide,
			ms2_tolerance=ms2_tolerance,
		)
		results.append(
			RowAnalysisResult(
				row_index=min(entry_a.row_index, entry_b.row_index),
				spectrum_title=result.spectrum_title,
				peptide_a=result.peptide_a,
				peptide_b=result.peptide_b,
				precursor_charge=result.precursor_charge,
				site_determining_count_a=result.site_determining_count_a,
				site_determining_count_b=result.site_determining_count_b,
				site_determining_ions_a=result.site_determining_ions_a,
				site_determining_ions_b=result.site_determining_ions_b,
			)
		)

	if output_csv_path:
		_write_results_csv(results, unmatched, output_csv_path)
	if output_unmatched_csv_path:
		_write_unmatched_csv(unmatched, output_unmatched_csv_path)
	if output_matched_diff_csv_path:
		_write_matched_diff_csv(results, output_matched_diff_csv_path)
	if output_unmatched_recheck_csv_path:
		recheck_rows = _build_unmatched_recheck_rows(
			unmatched=unmatched,
			target_rows_by_index=target_rows_by_index,
			spectra=spectra,
			method_a_name=method_a_name,
			method_b_name=method_b_name,
			peptide_a_col=peptide_a_col,
			peptide_b_col=peptide_b_col,
			ms2_tolerance=ms2_tolerance,
		)
		_write_unmatched_recheck_csv(recheck_rows, output_unmatched_recheck_csv_path)

	return results, unmatched


def analyze_pairs_from_csv_with_preset(
	csv_path: str | Path,
	mgf_dir_path: str | Path,
	mapping_name: str,
	ms2_tolerance: float,
	output_csv_path: Optional[str | Path] = None,
	output_unmatched_csv_path: Optional[str | Path] = None,
	output_matched_diff_csv_path: Optional[str | Path] = None,
	output_unmatched_recheck_csv_path: Optional[str | Path] = None,
) -> Tuple[List[RowAnalysisResult], List[UnmatchedPairRecord]]:
	"""Batch analysis entry using predefined column-mapping presets."""

	if mapping_name not in COLUMN_MAPPING_PRESETS:
		available = ", ".join(sorted(COLUMN_MAPPING_PRESETS))
		raise ValueError(
			f"Unknown mapping preset '{mapping_name}'. Available presets: {available}."
		)

	mapping = COLUMN_MAPPING_PRESETS[mapping_name]
	method_a = mapping.get("method_a")
	method_b = mapping.get("method_b")

	if not isinstance(method_a, dict) or not isinstance(method_b, dict):
		raise ValueError(
			f"Invalid mapping preset '{mapping_name}': method_a/method_b must be objects."
		)

	method_a_name = method_a.get("name")
	method_b_name = method_b.get("name")
	peptide_a_col = method_a.get("peptide_col")
	peptide_b_col = method_b.get("peptide_col")
	if not all(isinstance(x, str) for x in [method_a_name, method_b_name, peptide_a_col, peptide_b_col]):
		raise ValueError(
			f"Invalid mapping preset '{mapping_name}': method names and peptide columns must be strings."
		)

	return analyze_pairs_from_csv(
		csv_path=csv_path,
		mgf_dir_path=mgf_dir_path,
		method_a_name=method_a_name,
		method_b_name=method_b_name,
		peptide_a_col=peptide_a_col,
		peptide_b_col=peptide_b_col,
		ms2_tolerance=ms2_tolerance,
		output_csv_path=output_csv_path,
		output_unmatched_csv_path=output_unmatched_csv_path,
		output_matched_diff_csv_path=output_matched_diff_csv_path,
		output_unmatched_recheck_csv_path=output_unmatched_recheck_csv_path,
	)


def _derive_unmatched_output_path(output_csv_path: str | Path) -> Path:
	output_path = Path(output_csv_path)
	return output_path.with_name(f"{output_path.stem}_unmatched{output_path.suffix}")


def _derive_matched_diff_output_path(output_csv_path: str | Path) -> Path:
	output_path = Path(output_csv_path)
	return output_path.with_name(f"{output_path.stem}_matched_diff{output_path.suffix}")


def _derive_unmatched_recheck_output_path(output_csv_path: str | Path) -> Path:
	output_path = Path(output_csv_path)
	return output_path.with_name(f"{output_path.stem}_unmatched_recheck{output_path.suffix}")


def _write_results_csv(
	results: Sequence[RowAnalysisResult],
	unmatched: Sequence[UnmatchedPairRecord],
	output_csv_path: str | Path,
) -> None:
	fields = [
		"pair_status",
		"row_index",
		"spectrum_title",
		"peptide_a",
		"peptide_b",
		"precursor_charge",
		"site_determining_count_a",
		"site_determining_count_b",
		"site_determining_ions_a",
		"site_determining_ions_b",
		"present_method",
		"missing_method",
	]

	with Path(output_csv_path).open("w", encoding="utf-8", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=fields)
		writer.writeheader()
		for result in results:
			row = row_result_to_dict(result)
			row["pair_status"] = "matched"
			row["present_method"] = ""
			row["missing_method"] = ""
			writer.writerow(row)

		for record in unmatched:
			writer.writerow(
				{
					"pair_status": "unmatched",
					"row_index": record.present_row_index,
					"spectrum_title": record.spectrum_title,
					"peptide_a": record.peptide_a,
					"peptide_b": record.peptide_b,
					"precursor_charge": "",
					"site_determining_count_a": "",
					"site_determining_count_b": "",
					"site_determining_ions_a": "",
					"site_determining_ions_b": "",
					"present_method": record.present_method,
					"missing_method": record.missing_method,
				}
			)


def _write_unmatched_csv(unmatched: Sequence[UnmatchedPairRecord], output_csv_path: str | Path) -> None:
	fields = [
		"spectrum_title",
		"present_method",
		"missing_method",
		"present_row_index",
		"peptide_a",
		"peptide_b",
	]

	with Path(output_csv_path).open("w", encoding="utf-8", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=fields)
		writer.writeheader()
		for record in unmatched:
			writer.writerow(
				{
					"spectrum_title": record.spectrum_title,
					"present_method": record.present_method,
					"missing_method": record.missing_method,
					"present_row_index": record.present_row_index,
					"peptide_a": record.peptide_a,
					"peptide_b": record.peptide_b,
				}
			)


def _write_matched_diff_csv(results: Sequence[RowAnalysisResult], output_csv_path: str | Path) -> None:
	fields = [
		"pair_status",
		"row_index",
		"spectrum_title",
		"peptide_a",
		"peptide_b",
		"precursor_charge",
		"site_determining_count_a",
		"site_determining_count_b",
		"site_determining_ions_a",
		"site_determining_ions_b",
	]

	with Path(output_csv_path).open("w", encoding="utf-8", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=fields)
		writer.writeheader()
		for result in results:
			if result.peptide_a == result.peptide_b:
				continue
			row = row_result_to_dict(result)
			row["pair_status"] = "matched"
			writer.writerow({field: row.get(field, "") for field in fields})


def _build_unmatched_recheck_rows(
	unmatched: Sequence[UnmatchedPairRecord],
	target_rows_by_index: Mapping[int, Mapping[str, str]],
	spectra: Mapping[str, Spectrum],
	method_a_name: str,
	method_b_name: str,
	peptide_a_col: str,
	peptide_b_col: str,
	ms2_tolerance: float,
) -> List[Dict[str, object]]:
	rows: List[Dict[str, object]] = []

	for record in unmatched:
		present_peptide = record.peptide_a if record.peptide_a else record.peptide_b

		if record.missing_method == method_a_name:
			missing_col = peptide_a_col
		elif record.missing_method == method_b_name:
			missing_col = peptide_b_col
		else:
			raise ValueError(
				f"Unknown missing_method '{record.missing_method}' for spectrum '{record.spectrum_title}' "
				f"at row {record.present_row_index}."
			)

		row_data = target_rows_by_index.get(record.present_row_index)
		if row_data is None:
			raise ValueError(
				f"Cannot find source CSV row {record.present_row_index} for unmatched spectrum "
				f"'{record.spectrum_title}'."
			)

		candidate_missing_peptide = (row_data.get(missing_col) or "").strip()
		if not candidate_missing_peptide:
			rows.append(
				{
					"pair_status": "unmatched",
					"row_index": record.present_row_index,
					"spectrum_title": record.spectrum_title,
					"peptide_a": record.peptide_a,
					"peptide_b": record.peptide_b,
					"candidate_missing_peptide": "",
					"peptide_equal": "",
					"precursor_charge": "",
					"site_determining_count_a": "",
					"site_determining_count_b": "",
					"has_site_determining_ions": "",
					"site_determining_ions_a": "",
					"site_determining_ions_b": "",
					"present_method": record.present_method,
					"missing_method": record.missing_method,
				}
			)
			continue

		if record.missing_method == method_a_name:
			peptide_a = candidate_missing_peptide
			peptide_b = record.peptide_b
		else:
			peptide_a = record.peptide_a
			peptide_b = candidate_missing_peptide

		if candidate_missing_peptide == present_peptide:
			rows.append(
				{
					"pair_status": "matched",
					"row_index": record.present_row_index,
					"spectrum_title": record.spectrum_title,
					"peptide_a": peptide_a,
					"peptide_b": peptide_b,
					"candidate_missing_peptide": candidate_missing_peptide,
					"peptide_equal": True,
					"precursor_charge": "",
					"site_determining_count_a": "",
					"site_determining_count_b": "",
					"has_site_determining_ions": "",
					"site_determining_ions_a": "",
					"site_determining_ions_b": "",
					"present_method": record.present_method,
					"missing_method": record.missing_method,
				}
			)
			continue

		spectrum = spectra.get(record.spectrum_title)
		if spectrum is None:
			raise ValueError(
				f"Spectrum title '{record.spectrum_title}' not found in MGF while rechecking unmatched "
				f"row {record.present_row_index}."
			)

		result = analyze_localization_pair(
			spectrum=spectrum,
			peptide_a=peptide_a,
			peptide_b=peptide_b,
			ms2_tolerance=ms2_tolerance,
		)

		rows.append(
			{
				"pair_status": "matched",
				"row_index": record.present_row_index,
				"spectrum_title": result.spectrum_title,
				"peptide_a": result.peptide_a,
				"peptide_b": result.peptide_b,
				"candidate_missing_peptide": candidate_missing_peptide,
				"peptide_equal": result.peptide_a == result.peptide_b,
				"precursor_charge": result.precursor_charge,
				"site_determining_count_a": result.site_determining_count_a,
				"site_determining_count_b": result.site_determining_count_b,
				"has_site_determining_ions": (result.site_determining_count_a + result.site_determining_count_b) > 0,
				"site_determining_ions_a": _serialize_fragment_ions(result.site_determining_ions_a),
				"site_determining_ions_b": _serialize_fragment_ions(result.site_determining_ions_b),
				"present_method": record.present_method,
				"missing_method": record.missing_method,
			}
		)

	return rows


def _write_unmatched_recheck_csv(rows: Sequence[Mapping[str, object]], output_csv_path: str | Path) -> None:
	fields = [
		"pair_status",
		"row_index",
		"spectrum_title",
		"peptide_a",
		"peptide_b",
		"candidate_missing_peptide",
		"peptide_equal",
		"precursor_charge",
		"site_determining_count_a",
		"site_determining_count_b",
		"has_site_determining_ions",
		"site_determining_ions_a",
		"site_determining_ions_b",
		"present_method",
		"missing_method",
	]

	with Path(output_csv_path).open("w", encoding="utf-8", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=fields)
		writer.writeheader()
		for row in rows:
			writer.writerow({field: row.get(field, "") for field in fields})


def build_arg_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description="Calculate and match site-determining b/y ions for phosphosite localization pairs.",
	)
	parser.add_argument("--csv", required=True, help="Input CSV with spectrum and peptide columns.")
	parser.add_argument(
		"--mgf-dir",
		required=True,
		help="Input directory containing .mgf files (non-recursive).",
	)
	parser.add_argument(
		"--ms2-tolerance",
		required=True,
		type=float,
		help="Absolute m/z tolerance used for ion filtering and spectrum matching (in Dalton).",
	)
	parser.add_argument("--out", required=True, help="Output CSV path.")

	parser.add_argument(
		"--mapping-name",
		required=True,
		help="Optional preset mapping key, e.g. AB or CB.",
	)
	parser.add_argument(
		"--unmatched-out",
		default=None,
		help="Optional unmatched output CSV path. Defaults to <out_stem>_unmatched.csv.",
	)
	parser.add_argument(
		"--matched-diff-out",
		default=None,
		help="Optional matched-diff output CSV path. Defaults to <out_stem>_matched_diff.csv.",
	)
	parser.add_argument(
		"--unmatched-recheck-out",
		default=None,
		help="Optional unmatched recheck output CSV path. Defaults to <out_stem>_unmatched_recheck.csv.",
	)

	return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
	parser = build_arg_parser()
	args = parser.parse_args(argv)
	output_unmatched = args.unmatched_out
	if output_unmatched is None:
		output_unmatched = str(_derive_unmatched_output_path(args.out))
	output_matched_diff = args.matched_diff_out
	if output_matched_diff is None:
		output_matched_diff = str(_derive_matched_diff_output_path(args.out))
	output_unmatched_recheck = args.unmatched_recheck_out
	if output_unmatched_recheck is None:
		output_unmatched_recheck = str(_derive_unmatched_recheck_output_path(args.out))

	analyze_pairs_from_csv_with_preset(
		csv_path=args.csv,
		mgf_dir_path=args.mgf_dir,
		mapping_name=args.mapping_name,
		ms2_tolerance=args.ms2_tolerance,
		output_csv_path=args.out,
		output_unmatched_csv_path=output_unmatched,
		output_matched_diff_csv_path=output_matched_diff,
		output_unmatched_recheck_csv_path=output_unmatched_recheck,
	)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())



# python '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/cal_site_determining_ions/analyze_site_determ_ions.py' --csv '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/Figure3D_v2_Method1-8_all.csv' --mgf-dir '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/spectra' --ms2-tolerance 0.02 --out '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/output/site_determining_ions/2vs4.csv' --mapping-name 2vs4
# python '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/cal_site_determining_ions/analyze_site_determ_ions.py' --csv '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/Figure3D_v2_Method1-8_all.csv' --mgf-dir '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/spectra' --ms2-tolerance 0.02 --out '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/output/site_determining_ions/1vs4.csv' --mapping-name 1vs4
# python '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/cal_site_determining_ions/analyze_site_determ_ions.py' --csv '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/Figure3D_v2_Method1-8_all.csv' --mgf-dir '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/input/spectra' --ms2-tolerance 0.02 --out '/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/output/site_determining_ions/6vs8.csv' --mapping-name 6vs8