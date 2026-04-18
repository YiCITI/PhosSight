from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "analyze_site_determ_ions.py"
SPEC = importlib.util.spec_from_file_location("analyze_site_determ_ions", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_parse_modified_peptide_success() -> None:
    masses = MODULE.parse_modified_peptide("AS2C")
    assert len(masses) == 4
    assert masses[0] == pytest.approx(MODULE.AMINO_ACID_MASSES["A"])
    assert masses[1] == pytest.approx(MODULE.AMINO_ACID_MASSES["S"])
    assert masses[2] == pytest.approx(MODULE.MODIFICATION_MASSES["2"])
    assert masses[3] == pytest.approx(MODULE.AMINO_ACID_MASSES["C"])


def test_parse_modified_peptide_invalid_symbol() -> None:
    with pytest.raises(ValueError, match="Unsupported character"):
        MODULE.parse_modified_peptide("AC*")


def test_determine_fragment_charges() -> None:
    assert MODULE.determine_fragment_charges(2) == (1,)
    assert MODULE.determine_fragment_charges(3) == (1, 2)
    with pytest.raises(ValueError, match="missing"):
        MODULE.determine_fragment_charges(None)


def test_extract_site_determining_ions_respects_tolerance() -> None:
    ion_a = MODULE.FragmentIon("b", 1, 1, 100.000)
    ion_b = MODULE.FragmentIon("b", 1, 1, 100.030)
    ions_a = {("b", 1, 1): ion_a}
    ions_b = {("b", 1, 1): ion_b}

    site_determining_a, site_determining_b = MODULE.extract_site_determining_ions(
        ions_a, ions_b, ms2_tolerance=0.02
    )
    assert ("b", 1, 1) in site_determining_a
    assert ("b", 1, 1) in site_determining_b

    non_site_determining_a, non_site_determining_b = MODULE.extract_site_determining_ions(
        ions_a, ions_b, ms2_tolerance=0.05
    )
    assert non_site_determining_a == {}
    assert non_site_determining_b == {}


def test_match_ions_to_spectrum_selects_closest_peak() -> None:
    spectrum = MODULE.Spectrum(
        title="spec1",
        charge=2,
        peaks=(
            MODULE.Peak(mz=100.005, intensity=10.0),
            MODULE.Peak(mz=100.012, intensity=50.0),
            MODULE.Peak(mz=150.0, intensity=20.0),
        ),
    )
    ion = MODULE.FragmentIon("b", 1, 1, 100.010)
    matched = MODULE.match_ions_to_spectrum([ion], spectrum, ms2_tolerance=0.01)

    assert len(matched) == 1
    assert matched[0].observed_mz == pytest.approx(100.012)
    assert matched[0].observed_intensity == pytest.approx(50.0)


def test_analyze_localization_pair_keeps_only_matched_site_determining_ions() -> None:
    peptide_a = "AS2C"
    peptide_b = "ASTC"
    masses_a = MODULE.parse_modified_peptide(peptide_a)
    masses_b = MODULE.parse_modified_peptide(peptide_b)
    ions_a = MODULE.calculate_b_y_ions(masses_a, fragment_charges=(1,))
    ions_b = MODULE.calculate_b_y_ions(masses_b, fragment_charges=(1,))
    sd_a, _ = MODULE.extract_site_determining_ions(ions_a, ions_b, ms2_tolerance=0.02)

    # Keep exactly one theoretical site-determining ion matched in the spectrum.
    kept_ion = sorted(sd_a.values(), key=lambda x: (x.ion_type, x.index, x.charge))[0]
    spectrum = MODULE.Spectrum(
        title="spec_matched_only",
        charge=2,
        peaks=(
            MODULE.Peak(mz=kept_ion.mz + 0.005, intensity=100.0),
            MODULE.Peak(mz=999.0, intensity=10.0),
        ),
    )

    result = MODULE.analyze_localization_pair(
        spectrum=spectrum,
        peptide_a=peptide_a,
        peptide_b=peptide_b,
        ms2_tolerance=0.02,
    )

    assert result.site_determining_count_a == 1
    assert result.site_determining_count_b == 0
    assert len(result.site_determining_ions_a) == 1
    assert result.site_determining_ions_a[0].label == kept_ion.label


def test_calculate_all_theoretical_ions_for_as2c_with_precursor_charge_2() -> None:
    masses = MODULE.parse_modified_peptide("AS2C")
    fragment_charges = MODULE.determine_fragment_charges(2)
    ions = MODULE.calculate_b_y_ions(masses, fragment_charges)

    assert fragment_charges == (1,)
    assert set(ions.keys()) == {
        ("b", 1, 1),
        ("b", 2, 1),
        ("b", 3, 1),
        ("y", 1, 1),
        ("y", 2, 1),
        ("y", 3, 1),
    }

    expected_mz = {
        ("b", 1, 1): MODULE.AMINO_ACID_MASSES["A"] + MODULE.PROTON_MASS,
        ("b", 2, 1): MODULE.AMINO_ACID_MASSES["A"] + MODULE.AMINO_ACID_MASSES["S"] + MODULE.PROTON_MASS,
        ("b", 3, 1): MODULE.AMINO_ACID_MASSES["A"] + MODULE.AMINO_ACID_MASSES["S"] + MODULE.MODIFICATION_MASSES["2"] + MODULE.PROTON_MASS,
        ("y", 1, 1): MODULE.AMINO_ACID_MASSES["C"] + MODULE.C_TERM_MASS + MODULE.HYDROGEN_MASS + MODULE.PROTON_MASS,
        ("y", 2, 1): MODULE.MODIFICATION_MASSES["2"] + MODULE.AMINO_ACID_MASSES["C"] + MODULE.C_TERM_MASS + MODULE.HYDROGEN_MASS + MODULE.PROTON_MASS,
        ("y", 3, 1): MODULE.AMINO_ACID_MASSES["S"] + MODULE.MODIFICATION_MASSES["2"] + MODULE.AMINO_ACID_MASSES["C"] + MODULE.C_TERM_MASS + MODULE.HYDROGEN_MASS + MODULE.PROTON_MASS,
    }

    for ion_key, mz in expected_mz.items():
        assert ions[ion_key].mz == pytest.approx(mz)


def test_parse_mgf_and_collect_mgf_files(tmp_path: Path) -> None:
    mgf_dir = tmp_path / "spectra"
    mgf_dir.mkdir()
    mgf_path = mgf_dir / "a.mgf"
    _write_text(
        mgf_path,
        """BEGIN IONS
TITLE=spec1
CHARGE=2+
100.0 10
200.0 20
END IONS
""",
    )

    mgf_files = MODULE.collect_mgf_files(mgf_dir)
    assert mgf_files == (mgf_path,)

    spectra = MODULE.parse_mgf(mgf_path)
    assert "spec1" in spectra
    assert spectra["spec1"].charge == 2
    assert len(spectra["spec1"].peaks) == 2
    assert spectra["spec1"].peaks[0].mz == pytest.approx(100.0)
    assert spectra["spec1"].peaks[0].intensity == pytest.approx(10.0)
    assert spectra["spec1"].peaks[1].mz == pytest.approx(200.0)
    assert spectra["spec1"].peaks[1].intensity == pytest.approx(20.0)


def test_analyze_pairs_from_csv_with_preset(tmp_path: Path) -> None:
    mgf_dir = tmp_path / "spectra"
    mgf_dir.mkdir()
    mgf_path = mgf_dir / "input.mgf"
    _write_text(
        mgf_path,
        """BEGIN IONS
TITLE=spec_match
CHARGE=2+
50.0 5
100.0 20
200.0 10
END IONS
BEGIN IONS
TITLE=spec_unmatched
CHARGE=2+
80.0 7
160.0 14
END IONS
""",
    )

    csv_path = tmp_path / "input.csv"
    rows = [
        {
            "Method": "Method2",
            "Title": "spec_match",
            "IsoformSequence_autort_pDeep": "AS2C",
            "IsoformSequence_autort_pDeep_phosSight": "ASTC",
        },
        {
            "Method": "Method4",
            "Title": "spec_match",
            "IsoformSequence_autort_pDeep": "AS2C",
            "IsoformSequence_autort_pDeep_phosSight": "ASTC",
        },
        {
            "Method": "Method2",
            "Title": "spec_unmatched",
            "IsoformSequence_autort_pDeep": "AS2C",
            "IsoformSequence_autort_pDeep_phosSight": "ASTC",
        },
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Method",
                "Title",
                "IsoformSequence_autort_pDeep",
                "IsoformSequence_autort_pDeep_phosSight",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    out_path = tmp_path / "out.csv"
    unmatched_out = tmp_path / "unmatched.csv"
    matched_diff_out = tmp_path / "matched_diff.csv"

    results, unmatched = MODULE.analyze_pairs_from_csv_with_preset(
        csv_path=csv_path,
        mgf_dir_path=mgf_dir,
        mapping_name="2vs4",
        ms2_tolerance=0.02,
        output_csv_path=out_path,
        output_unmatched_csv_path=unmatched_out,
        output_matched_diff_csv_path=matched_diff_out,
    )

    assert len(results) == 1
    assert results[0].spectrum_title == "spec_match"
    assert len(unmatched) == 1
    assert unmatched[0].spectrum_title == "spec_unmatched"
    assert out_path.exists()
    assert unmatched_out.exists()
    assert matched_diff_out.exists()


def test_main_uses_derived_output_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    csv_path = tmp_path / "in.csv"
    mgf_dir = tmp_path / "spectra"
    out_path = tmp_path / "result.csv"
    mgf_dir.mkdir()
    csv_path.write_text("Method,Title,IsoformSequence_autort_pDeep,IsoformSequence_autort_pDeep_phosSight\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def _fake_analyze(**kwargs: object) -> tuple[list[object], list[object]]:
        captured.update(kwargs)
        return [], []

    monkeypatch.setattr(MODULE, "analyze_pairs_from_csv_with_preset", _fake_analyze)

    exit_code = MODULE.main(
        [
            "--csv",
            str(csv_path),
            "--mgf-dir",
            str(mgf_dir),
            "--ms2-tolerance",
            "0.02",
            "--out",
            str(out_path),
            "--mapping-name",
            "2vs4",
        ]
    )

    assert exit_code == 0
    assert captured["output_unmatched_csv_path"] == str(tmp_path / "result_unmatched.csv")
    assert captured["output_matched_diff_csv_path"] == str(tmp_path / "result_matched_diff.csv")
    assert captured["output_unmatched_recheck_csv_path"] == str(tmp_path / "result_unmatched_recheck.csv")
