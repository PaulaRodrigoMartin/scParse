import os
import pandas as pd

root_path = "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/starsolo_v1/"

# Samples used in the analysis, modify as necessary:
sp_sample_dict = {
	"CHMP":["Chimpanzee_Stephan","Chimpanzee_marlock"],
	"GUB":["Guinea_Baboon_5","Guinea_Baboon_6"],
	"SGB":["Siamang_Gibbon_Coco"],
	"GLB":["Gelada_Baboon_1"],
	"BWM":["Brown_Wooly_Monkey_1"],
	"PTM":["Pigtailed_macaque_1"],
	"WCG":["Whitecheeked_gibbon"],
	"DIM":["Diana_Monkey_Suacoco"],
	"LRG":["Lar_Gibbon_Tarzan"]
}

sum_rows = [
    "Number of Reads",
    # "Reads Mapped to Genome: Unique+Multiple",
    # "Reads Mapped to Genome: Unique",
    # "Reads Mapped to GeneFull: Unique+Multiple GeneFull",
    # "Reads Mapped to GeneFull: Unique GeneFull",
    "Estimated Number of Cells",
    "Unique Reads in Cells Mapped to GeneFull",
    "UMIs in Cells",
    "Total GeneFull Detected"
]

sublibraries = ["PARSE1_UDI_WT_1", "PARSE2_UDI_WT_2", "PARSE3_UDI_WT_3", "PARSE4_UDI_WT_4", "PARSE6_UDI_WT_4"]

all_samples_data = []

# Loop through each species in the dictionary
for species, samples in sp_sample_dict.items():
    # Loop through each sample for the current species
    for sample_name in samples:
        print(f"Processing sample: {sample_name}")
        
        # Create an empty DataFrame to store combined data for the current sample
        combined_df = pd.DataFrame()
        
        # Loop through each sublibrary
        for sublib in sublibraries:
            # Construct the file path for the summary.csv
            file_path = os.path.join(root_path, species, sample_name, sublib, "Solo.out", "GeneFull", "Summary.csv")
            
            # Check if the summary.csv file exists
            if os.path.exists(file_path):
                # Read the summary.csv file into a DataFrame
                df = pd.read_csv(file_path, header=None, index_col=0)
                
                # Rename the column to match the sublibrary
                df.columns = [sublib]
                
                # Merge into the combined DataFrame
                if combined_df.empty:
                    combined_df = df
                else:
                    combined_df = combined_df.join(df)
            else:
                print(f"File not found: {file_path}")

        if not combined_df.empty:
            # Calculate the sum across all rows and add as a new column named 'all'
            all_column = combined_df.apply(lambda x: x.sum() if x.name in sum_rows else x.mean(), axis=1)

            # depending on the column sum it or make the mean 
            combined_df['all'] = all_column

            #################
            # Append the sample name and the combined data to the list
            sample_row = pd.DataFrame({sample_name: [sample_name]}, index=[""])  # Creates a row with the sample name
            all_samples_data.append(sample_row)            
            all_samples_data.append(combined_df)
        else:
            print(f"No data to save for sample: {sample_name}")
            #################

        # Define output path for the combined summary of the current sample
        output_path = f'{root_path}/Summary_starsolo_{sample_name}_all.csv'
        
        # Save the combined DataFrame to a new CSV file if it's not empty
        if not combined_df.empty:
            combined_df.to_csv(output_path)
            print(f"Combined summary saved to {output_path}")
        else:
            print(f"No data to save for sample: {sample_name}")


# import os
# import pandas as pd

csv_files = []

# Walk through all subdirectories and files in the root directory
for subdir, dirs, files in os.walk(root_path):
    for file in files:
        if 'Summary_starsolo' in file:  # Check if 'Summary_starsolo' is in the filename
            # Get the full path of the file
            file_path = os.path.join(subdir, file)
            # Append the full path to the csv_files list
            csv_files.append(file_path)

combined_df = pd.DataFrame()

# Process each CSV file
for file_path in csv_files:
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path, index_col= 0)
    
    # Extract the columns named 'all'
    all_columns = df.filter(like='all')
    
    # Generate the new column names by removing 'Summary_starsolo' and '_all'
    base_name = os.path.basename(file_path)
    new_column_prefix = base_name.replace('Summary_starsolo_', '').replace('_all.csv', '')
    all_columns.columns = [f"{new_column_prefix}_{col}" for col in all_columns.columns]
    
    # Concatenate the columns to the combined DataFrame
    combined_df = pd.concat([combined_df, all_columns], axis=1)

# Save the combined DataFrame to a new CSV file
combined_df.to_csv(f'{root_path}combined_summary.csv', index=True)

print("Combined CSV file created successfully: combined_summary.csv")

# output_excel_file = 'combined_output.xlsx'


# writer = pd.ExcelWriter(output_excel_file, engine='openpyxl')

# # Loop through each CSV file and write it to a different sheet in the Excel file
# for csv_file in csv_files:
#     # Read the CSV file into a DataFrame
#     df = pd.read_csv(csv_file)
    
#     # Use the CSV filename without extension as the sheet name
#     sheet_name = os.path.splitext(os.path.basename(csv_file))[0]
    
#     # Write the DataFrame to a specific sheet
#     df.to_excel(writer, sheet_name=sheet_name, index=False)

# # Save the Excel file
# writer.save()

