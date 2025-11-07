import pyarrow.parquet as pq
import pyarrow.compute as pc

def filter_parquet_file(input_path, output_path, column_name, filter_string):
    """
    Read a Parquet file, filter rows where the specified column contains the given substring
    at most once, and save the result to a new Parquet file.

    Args:
        input_path (str): Path to the input Parquet file.
        output_path (str): Path to the output Parquet file.
        column_name (str): Name of the column to filter.
        filter_string (str): Substring to search for in the column.
    """
    try:
        # Read the Parquet file
        table = pq.read_table(input_path)

        # Ensure the specified column exists
        if column_name not in table.column_names:
            print(f"Error: Column '{column_name}' not found. Available columns: {table.column_names}")
            return

        # Get the target column data
        column_data = table[column_name]

        # Use PyArrow compute functions for efficient filtering
        # 1. Count occurrences of filter_string in each value
        #    pa.compute.count_substring returns null for null input
        string_counts = pc.count_substring(column_data, pattern=filter_string)

        # 2. Replace nulls with 0 (treat null as zero occurrences)
        counts_no_nulls = string_counts.fill_null(0)

        # 3. Create a boolean mask for rows where count <= 1
        mask = pc.less_equal(counts_no_nulls, 1)

        # Filter the table using the mask
        filtered_table = table.filter(mask)

        # Write the filtered table to a new Parquet file
        pq.write_table(filtered_table, output_path)
        print(f"Filtering complete. Output saved to: {output_path}")
        print(f"Original row count: {len(table)}")
        print(f"Filtered row count: {len(filtered_table)}")

    except Exception as e:
        print(f"Error processing Parquet file: {e}")


