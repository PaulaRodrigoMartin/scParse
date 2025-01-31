
import sys, os, re
from pathlib import Path
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack, vstack, csr_matrix
import shutil

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

sublibraries = ["PARSE1_UDI_WT_1", "PARSE2_UDI_WT_2", "PARSE3_UDI_WT_3", "PARSE4_UDI_WT_4", "PARSE6_UDI_WT_4"]


BASE_DIR = Path(__file__).parent

out_path = BASE_DIR / "starsolo_v1_t2t"
out_path = str(out_path.resolve())

# Get the paths for the matrices for each sublibrary per sample (GeneFull)
path_matrices = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_matrices[f'{sp}/{sample}/combined_fil'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/matrix.mtx'
			path_matrices[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

path_tsvs = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_tsvs[f'{sp}/{sample}/combined_fil'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/barcodes.tsv'
			path_tsvs[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

path_features = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_features[f'{sp}/{sample}/combined_fil'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/features.tsv'
			path_features[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

# # For the raw files
path_tsvs_raw = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_tsvs_raw[f'{sp}/{sample}/combined'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/barcodes.tsv'
			path_tsvs_raw[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)

path_features_raw = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_features_raw[f'{sp}/{sample}/combined'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/features.tsv'
			path_features_raw[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)


# # # For velocyto raw, unspliced
path_unspliced = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_unspliced[f'{sp}/{sample}/Velocyto'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/unspliced.mtx'
			path_unspliced[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# # For velocyto raw, spliced
path_spliced = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_spliced[f'{sp}/{sample}/Velocyto'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/spliced.mtx'
			path_spliced[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# # For velocyto raw, ambiguous
path_ambiguous = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_ambiguous[f'{sp}/{sample}/Velocyto'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/ambiguous.mtx'
			path_ambiguous[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# ## tsvs and barcodes from velocyto
path_tsvs_velocyto = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_tsvs_velocyto[f'{sp}/{sample}/Velocyto'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/barcodes.tsv'
			path_tsvs_velocyto[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

path_features_velocyto = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_features_velocyto[f'{sp}/{sample}/Velocyto'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/features.tsv'
			path_features_velocyto[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

path_unique_and_multi = {}
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:
		path_unique_and_multi[f'{sp}/{sample}/combined'] = []
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		for sublib in sublibraries:
			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx'
			path_unique_and_multi[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)


# #############################################################
# # # Add the sublibrary tag to each cell name
for sublib in sublibraries:
	for key, values in path_tsvs.items():
		# print(f"{key}: {value}")
		for value in values:
			if sublib in value:
				df = pd.read_csv(value, sep='\t', header=None)
	            # Add the sublibrary to the cell name
				df[0] = df[0] + '_' + sublib
	            # Save the modified df back to a TSV
				df.to_csv(value, sep='\t', index=False, header=False)
				print(f"Processed {value} with sublibrary {sublib}")
# print(path_tsvs_raw)
for sublib in sublibraries:
	for key, values in path_tsvs_raw.items():
		# print(f"{key}: {value}")
		for value in values:
			if sublib in value:
				df = pd.read_csv(value, sep='\t', header=None)
	            # Add the sublibrary to the cell name
				df[0] = df[0] + '_' + sublib
	            # Save the modified df back to a TSV
				df.to_csv(value, sep='\t', index=False, header=False)
				print(f"Processed {value} with sublibrary {sublib}")

for sublib in sublibraries:
	for key, values in path_tsvs_velocyto.items():
		# print(f"{key}: {value}")
		for value in values:
			if sublib in value:
				df = pd.read_csv(value, sep='\t', header=None)
	            # Add the sublibrary to the cell name
				df[0] = df[0] + '_' + sublib
	            # Save the modified df back to a TSV
				df.to_csv(value, sep='\t', index=False, header=False)
				print(f"Processed {value} with sublibrary {sublib}")

############################################################
# # # # # # ### Combine matrices
combined_matrices = {}
for key in sorted(path_matrices.keys()):
    matrices = []  # list to hold all matrices for the current key

    for file_path in path_matrices[key]:
        matrix = mmread(file_path).tocsr()  # Read the matrix from the file
        matrices.append(matrix)  # Add the matrix to the list
    # print(matrices)
    # Combine all matrices for the current key into one big matrix
    big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
    print(big_matrix.shape)
	# Create a directory for the current key to store the combined matrix
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined matrix
    output_file_path = os.path.join(key_directory, f'matrix.mtx')
    
    # Save the combined matrix to the file
    mmwrite(output_file_path, big_matrix)
    
    combined_matrices[key] = big_matrix
    print(f"Combined matrix for {key} saved to {output_file_path}")

# # for unspliced
combined_matrices_unspliced = {}
for key in sorted(path_unspliced.keys()):
    matrices = []  # list to hold all matrices for the current key

    for file_path in path_unspliced[key]:
        matrix = mmread(file_path).tocsr()  # Read the matrix from the file
        matrices.append(matrix)  # Add the matrix to the list
    # print(matrices)
    # Combine all matrices for the current key into one big matrix
    big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
    print(big_matrix.shape)
	# Create a directory for the current key to store the combined matrix
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined matrix
    output_file_path = os.path.join(key_directory, f'unspliced.mtx')
    
    # Save the combined matrix to the file
    mmwrite(output_file_path, big_matrix)
    
    combined_matrices_unspliced[key] = big_matrix
    print(f"Combined matrix for {key} saved to {output_file_path}")

# # for spliced
combined_matrices_spliced = {}
for key in sorted(path_spliced.keys()):
    matrices = []  # list to hold all matrices for the current key

    for file_path in path_spliced[key]:
        matrix = mmread(file_path).tocsr()  # Read the matrix from the file
        matrices.append(matrix)  # Add the matrix to the list
    # print(matrices)
    # Combine all matrices for the current key into one big matrix
    big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
    print(big_matrix.shape)
	# Create a directory for the current key to store the combined matrix
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined matrix
    output_file_path = os.path.join(key_directory, f'spliced.mtx')
    
    # Save the combined matrix to the file
    mmwrite(output_file_path, big_matrix)
    
    combined_matrices_spliced[key] = big_matrix
    print(f"Combined matrix for {key} saved to {output_file_path}")


# # for ambiguous
combined_matrices_ambiguous = {}
for key in sorted(path_ambiguous.keys()):
    matrices = []  # list to hold all matrices for the current key

    for file_path in path_ambiguous[key]:
        matrix = mmread(file_path).tocsr()  # Read the matrix from the file
        matrices.append(matrix)  # Add the matrix to the list
    # print(matrices)
    # Combine all matrices for the current key into one big matrix
    big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
    print(big_matrix.shape)
	# Create a directory for the current key to store the combined matrix
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined matrix
    output_file_path = os.path.join(key_directory, f'ambiguous.mtx')
    
    # Save the combined matrix to the file
    mmwrite(output_file_path, big_matrix)
    
    combined_matrices_ambiguous[key] = big_matrix
    print(f"Combined matrix for {key} saved to {output_file_path}")

# # for raw/UniqueandMult-EM.mtx
combined_matrices_mult = {}
for key in sorted(path_unique_and_multi.keys()):
    matrices = []  # list to hold all matrices for the current key

    for file_path in path_unique_and_multi[key]:
        matrix = mmread(file_path).tocsr()  # Read the matrix from the file
        matrices.append(matrix)  # Add the matrix to the list
    # print(matrices)
    # Combine all matrices for the current key into one big matrix
    big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
    print(big_matrix.shape)
	# Create a directory for the current key to store the combined matrix
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined matrix
    output_file_path = os.path.join(key_directory, f'UniqueAndMult-EM.mtx')
    
    # Save the combined matrix to the file
    mmwrite(output_file_path, big_matrix)
    
    combined_matrices_mult[key] = big_matrix
    print(f"Combined matrix for {key} saved to {output_file_path}")

# # # ## Combine tsvs
# ## filtered matrices
for key in sorted(path_tsvs.keys()):
    dfs = []
    # Read each TSV into a df and append to list
    for file_path in path_tsvs[key]:
        df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
        dfs.append(df)  # Add the DataFrame to the list
    
    # Concatenate all dfs for the current key into one big df
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Create a directory for the current key to store the combined TSV
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined TSV
    output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
    # Save the combined DataFrame to a TSV file
    combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
    print(f"Combined TSV for {key} saved to {output_file_path}")

# # for raw matrix
for key in sorted(path_tsvs_raw.keys()):
    dfs = []
    # Read each TSV into a df and append to list
    for file_path in path_tsvs_raw[key]:
        df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
        dfs.append(df)  # Add the DataFrame to the list
    
    # Concatenate all dfs for the current key into one big df
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Create a directory for the current key to store the combined TSV
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined TSV
    output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
    # Save the combined DataFrame to a TSV file
    combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
    print(f"Combined TSV for {key} saved to {output_file_path}")

# # raw/velocyto
for key in sorted(path_tsvs_velocyto.keys()):
    dfs = []
    # Read each TSV into a df and append to list
    for file_path in path_tsvs_velocyto[key]:
        df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
        dfs.append(df)  # Add the DataFrame to the list
    
    # Concatenate all dfs for the current key into one big df
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Create a directory for the current key to store the combined TSV
    key_directory = os.path.join(out_path, key)
    os.makedirs(key_directory, exist_ok=True)
    
    # Define the output file path for the combined TSV
    output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
    # Save the combined DataFrame to a TSV file
    combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
    print(f"Combined TSV for {key} saved to {output_file_path}")

# # # ## Add features (genes) to combined folder
for key, feature_paths in path_features.items():
    # Create the combined directory path
    combined_dir = os.path.join(out_path, key)
    
    # Ensure the combined directory exists
    os.makedirs(combined_dir, exist_ok=True)
    
    # Copy each features.tsv file to the combined directory
    for feature_path in feature_paths:
        if os.path.isfile(feature_path):  # Check if the source file exists
            destination_path = os.path.join(combined_dir, 'features.tsv')
            shutil.copy2(feature_path, destination_path)
            print(f"Copied {feature_path} to {destination_path}")
        else:
            print(f"File not found: {feature_path}")

# # for raw matrix
for key, feature_paths in path_features_raw.items():
    # Create the combined directory path
    combined_dir = os.path.join(out_path, key)
    
    # Ensure the combined directory exists
    os.makedirs(combined_dir, exist_ok=True)
    
    # Copy each features.tsv file to the combined directory
    for feature_path in feature_paths:
        if os.path.isfile(feature_path):  # Check if the source file exists
            destination_path = os.path.join(combined_dir, 'features.tsv')
            shutil.copy2(feature_path, destination_path)
            print(f"Copied {feature_path} to {destination_path}")
        else:
            print(f"File not found: {feature_path}")

for key, feature_paths in path_features_velocyto.items():
    # Create the combined directory path
    combined_dir = os.path.join(out_path, key)
    
    # Ensure the combined directory exists
    os.makedirs(combined_dir, exist_ok=True)
    
    # Copy each features.tsv file to the combined directory
    for feature_path in feature_paths:
        if os.path.isfile(feature_path):  # Check if the source file exists
            destination_path = os.path.join(combined_dir, 'features.tsv')
            shutil.copy2(feature_path, destination_path)
            print(f"Copied {feature_path} to {destination_path}")
        else:
            print(f"File not found: {feature_path}")