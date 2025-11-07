import os
import tempfile
import pyarrow as pa
import pytest
from filter_parquet.filter_parquet_using_mod import filter_parquet_file

import pyarrow.parquet as pq


@pytest.fixture
def sample_parquet_file(tmp_path):
    # Create a sample table
    data = {
        "Modified.Sequence": [
            "PEPTIDE[UniMod:21]K",
            "PEPTIDEK",
            "PEPTIDE[UniMod:21][UniMod:21]K",
            None,
            "PEPTIDE[UniMod:21]K[UniMod:21]",
            "PEPTIDE[UniMod:35]K"
        ],
        "OtherColumn": [1, 2, 3, 4, 5, 6]
    }
    table = pa.table(data)
    file_path = tmp_path / "input.parquet"
    pq.write_table(table, file_path)
    return file_path

def read_column_from_parquet(file_path, column):
    table = pq.read_table(file_path)
    return table[column].to_pylist()

def test_filter_parquet_file_basic(tmp_path, sample_parquet_file):
    output_file = tmp_path / "output.parquet"
    filter_parquet_file(
        str(sample_parquet_file),
        str(output_file),
        "Modified.Sequence",
        "UniMod:21"
    )
    result = read_column_from_parquet(output_file, "Modified.Sequence")
    # Only rows with at most one "UniMod:21" should remain (including None)
    assert result == [
        "PEPTIDE[UniMod:21]K",
        "PEPTIDEK",
        None,
        "PEPTIDE[UniMod:35]K"
    ]

def test_filter_parquet_file_column_not_found(tmp_path, sample_parquet_file, capsys):
    output_file = tmp_path / "output2.parquet"
    filter_parquet_file(
        str(sample_parquet_file),
        str(output_file),
        "NonExistentColumn",
        "UniMod:21"
    )
    captured = capsys.readouterr()
    assert "Error: Column 'NonExistentColumn' not found" in captured.out
    assert not os.path.exists(output_file)

def test_filter_parquet_file_no_matches(tmp_path):
    # All rows have more than one occurrence
    data = {
        "Modified.Sequence": [
            "A[UniMod:21][UniMod:21]",
            "B[UniMod:21][UniMod:21][UniMod:21]",
        ]
    }
    table = pa.table(data)
    input_file = tmp_path / "input3.parquet"
    output_file = tmp_path / "output3.parquet"
    pq.write_table(table, input_file)
    filter_parquet_file(
        str(input_file),
        str(output_file),
        "Modified.Sequence",
        "UniMod:21"
    )
    result = read_column_from_parquet(output_file, "Modified.Sequence")
    assert result == []

def test_filter_parquet_file_empty_file(tmp_path):
    # No rows
    data = {"Modified.Sequence": pa.array([], pa.string())}
    table = pa.table(data)
    input_file = tmp_path / "input4.parquet"
    output_file = tmp_path / "output4.parquet"
    pq.write_table(table, input_file)
    filter_parquet_file(
        str(input_file),
        str(output_file),
        "Modified.Sequence",
        "UniMod:21"
    )
    result = read_column_from_parquet(output_file, "Modified.Sequence")
    assert result == []