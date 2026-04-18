import argparse
import csv


def filter_csv(input_csv: str, output_csv: str, target_present_method: str) -> None:
	with open(input_csv, "r", newline="", encoding="utf-8") as infile:
		reader = csv.DictReader(infile)
		if reader.fieldnames is None:
			raise ValueError("Input CSV has no header.")

		required_columns = {"peptide_equal", "present_method"}
		missing_columns = required_columns - set(reader.fieldnames)
		if missing_columns:
			raise ValueError(
				f"Missing required column(s): {', '.join(sorted(missing_columns))}"
			)

		with open(output_csv, "w", newline="", encoding="utf-8") as outfile:
			writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
			writer.writeheader()

			for row in reader:
				peptide_equal_value = (row.get("peptide_equal") or "").strip().upper()
				present_method_value = (row.get("present_method") or "").strip()
				if peptide_equal_value == "FALSE" and present_method_value == target_present_method:
					writer.writerow(row)


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"Filter CSV rows where peptide_equal is FALSE and present_method equals a target value."
		)
	)
	parser.add_argument("input_csv", help="Path to the input CSV file")
	parser.add_argument("output_csv", help="Path to the output CSV file")
	parser.add_argument("present_method", help="Target value for present_method")
	args = parser.parse_args()

	filter_csv(args.input_csv, args.output_csv, args.present_method)


if __name__ == "__main__":
	main()
