import pandas as pd

def convert_parquet_to_csv_select_columns(input_parquet_path, output_csv_path, columns_to_keep, is_drop_duplicate=False):
    """
    Use Pandas to read a Parquet file, select specified columns, optionally remove duplicates, and save as a CSV file (with header).

    Args:
    input_parquet_path (str): Path to the input Parquet file.
    output_csv_path (str): Path to the output CSV file.
    columns_to_keep (list): List of column names to keep.
    is_drop_duplicate (bool): Whether to drop duplicate rows based on the "Modified.Sequence" column.
    """
    try:
        # Read the Parquet file using Pandas
        df = pd.read_parquet(input_parquet_path)

        # Check if the specified columns exist and only select existing columns
        existing_columns_to_keep = [col for col in columns_to_keep if col in df.columns]
        
        if not existing_columns_to_keep:
            print("Error: None of the specified columns exist in the Parquet file.")
            return

        if len(existing_columns_to_keep) < len(columns_to_keep):
            missing_cols = set(columns_to_keep) - set(existing_columns_to_keep)
            print(f"Warning: The following specified columns do not exist and will be ignored: {missing_cols}")

        # Select the specified columns from the DataFrame
        selected_df = df[existing_columns_to_keep]
        
        # If is_drop_duplicate is True and "Modified.Sequence" column exists, drop duplicates
        if is_drop_duplicate and "Modified.Sequence" in selected_df.columns:
            selected_df = selected_df.drop_duplicates(subset=["Modified.Sequence"])

        # Write the selected DataFrame to a CSV file, without index, but with header
        selected_df.to_csv(output_csv_path, index=False)
        print(f"Successfully saved selected columns from '{input_parquet_path}' to '{output_csv_path}'")

    except Exception as e:
        print(f"An error occurred while processing the file: {e}")

# if __name__ == "__main__":
#     # Please replace with your actual file paths and column names
#     input_file = "test/test.parquet"  # e.g., "test/test.parquet"
#     output_file = "test/output.csv"    # e.g., "test/output.csv"
    
#     # Columns you want to keep
#     columns = ["Protein.Ids", "Protein.Group", "Modified.Sequence"]
    
#     # Call the function
#     convert_parquet_to_csv_select_columns(input_file, output_file, columns, is_drop_duplicate=True)
