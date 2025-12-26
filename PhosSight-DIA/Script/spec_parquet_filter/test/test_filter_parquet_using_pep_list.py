import os
import tempfile
import pyarrow as pa
import pytest
import pyarrow.parquet as pq

from filter_parquet.filter_parquet_using_pep_list import (
    process_peptide_to_modified_sequence,
    _read_and_process_peptides_from_txt,
    filter_parquet_by_peptide_list,
)

def test_process_peptide_to_modified_sequence_basic():
    seq, count = process_peptide_to_modified_sequence("PEPTIDEK")
    assert seq == "PEPTIDEK"
    assert count == 0

def test_process_peptide_to_modified_sequence_phospho():
    seq, count = process_peptide_to_modified_sequence("PEPTsDEK")
    assert seq == "PEPTS(UniMod:21)DEK"
    assert count == 1

def test_process_peptide_to_modified_sequence_multiple_phospho():
    seq, count = process_peptide_to_modified_sequence("MULTyPLEtsITES")
    assert seq == "MULTY(UniMod:21)PLET(UniMod:21)S(UniMod:21)ITES"
    assert count == 3

def test_process_peptide_to_modified_sequence_carbamidomethyl():
    seq, count = process_peptide_to_modified_sequence("ACDC")
    assert seq == "AC(UniMod:4)DC(UniMod:4)"
    assert count == 0

def test_process_peptide_to_modified_sequence_other_lowercase(monkeypatch):
    warnings = []
    monkeypatch.setattr("builtins.print", lambda msg: warnings.append(msg))
    seq, count = process_peptide_to_modified_sequence("ANOhERPEPTIDE")
    assert "Warning: Peptide 'ANOhERPEPTIDE' contains lowercase letter 'h'" in "".join(warnings)
    assert seq == "ANOHERPEPTIDE"
    assert count == 0

def test_read_and_process_peptides_from_txt(tmp_path):
    txt = tmp_path / "peps.txt"
    txt.write_text("PEPTIDEK\nPEPTsDEK\nsty\n")
    result = _read_and_process_peptides_from_txt(str(txt))
    assert "PEPTIDEK" in result
    assert "PEPTS(UniMod:21)DEK" in result
    assert "S(UniMod:21)T(UniMod:21)Y(UniMod:21)" not in result

def test_read_and_process_peptides_from_txt_empty(tmp_path):
    txt = tmp_path / "empty.txt"
    txt.write_text("")
    result = _read_and_process_peptides_from_txt(str(txt))
    assert result == set()

def test_read_and_process_peptides_from_txt_file_not_found():
    result = _read_and_process_peptides_from_txt("nonexistent_file.txt")
    assert result == set()

def make_test_parquet(path, modified_sequences):
    table = pa.table({
        "Modified.Sequence": pa.array(modified_sequences),
        "OtherCol": pa.array([1]*len(modified_sequences)),
    })
    pq.write_table(table, path)

def test_filter_parquet_by_peptide_list(tmp_path):
    # Prepare peptide txt
    pep_txt = tmp_path / "peps.txt"
    pep_txt.write_text("PEPTIDEK\nPEPTsDEK\nsty\n")
    # Prepare parquet
    parquet_path = tmp_path / "test.parquet"
    out_path = tmp_path / "filtered.parquet"
    make_test_parquet(str(parquet_path), [
        "PEPTIDEK",
        "PEPTS(UniMod:21)DEK",
        "S(UniMod:21)T(UniMod:21)Y(UniMod:21)",
        "NOTINLIST"
    ])
    # Run filter
    filter_parquet_by_peptide_list(str(pep_txt), str(parquet_path), str(out_path))
    # Check output
    filtered = pq.read_table(str(out_path))
    filtered_seqs = set(filtered["Modified.Sequence"].to_pylist())
    assert filtered_seqs == {
        "PEPTIDEK",
        "PEPTS(UniMod:21)DEK"
    }
    assert "NOTINLIST" not in filtered_seqs

def test_filter_parquet_by_peptide_list_empty_pep_list(tmp_path, capsys):
    pep_txt = tmp_path / "empty.txt"
    pep_txt.write_text("")
    parquet_path = tmp_path / "test.parquet"
    out_path = tmp_path / "filtered.parquet"
    make_test_parquet(str(parquet_path), ["PEPTIDEK"])
    filter_parquet_by_peptide_list(str(pep_txt), str(parquet_path), str(out_path))
    # Output file should not exist or be empty
    assert not out_path.exists() or pq.read_table(str(out_path)).num_rows == 0

def test_filter_parquet_by_peptide_list_missing_column(tmp_path, capsys):
    pep_txt = tmp_path / "peps.txt"
    pep_txt.write_text("PEPTIDEK\n")
    parquet_path = tmp_path / "test.parquet"
    out_path = tmp_path / "filtered.parquet"
    # Parquet without the required column
    table = pa.table({"OtherCol": pa.array([1, 2])})
    pq.write_table(table, str(parquet_path))
    filter_parquet_by_peptide_list(str(pep_txt), str(parquet_path), str(out_path))
    captured = capsys.readouterr()
    assert "Error: Column 'Modified.Sequence' not found" in captured.out

def test_filter_parquet_by_peptide_list_input_parquet_not_found(tmp_path, capsys):
    pep_txt = tmp_path / "peps.txt"
    pep_txt.write_text("PEPTIDEK\n")
    parquet_path = tmp_path / "doesnotexist.parquet"
    out_path = tmp_path / "filtered.parquet"
    filter_parquet_by_peptide_list(str(pep_txt), str(parquet_path), str(out_path))
    captured = capsys.readouterr()
    assert "Error: Input Parquet file not found" in captured.out