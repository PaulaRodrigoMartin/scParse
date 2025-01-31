from gwf import Workflow, AnonymousTarget

gwf=Workflow()
import sys, os, re
import pandas as pd
import fnmatch

import scanpy
import anndata
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack, vstack, csr_matrix 
import shutil
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pathlib import Path


###################
# INPUT DATA #
###################

#string indicating version of parse kit, can be v1, v2, v3
chemistry = "v1" 

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

# wells for each sample
sp_well_dict = {
	"Chimpanzee_Stephan":["A1,A2,A3,A4,D11,D12"],
	"Guinea_Baboon_5" :["A5,A6,A7,A8"],
	"Siamang_Gibbon_Coco":["A9,A10,A11,A12"],
	"Gelada_Baboon_1":["B1,B2,B3,B4"],
	"Brown_Wooly_Monkey_1":["B5,B6,B7,B8"],
	"Pigtailed_macaque_1" :["B9,B10,B11,B12"],
	"Whitecheeked_gibbon" :["C1,C2,C3,C4"],
	"Diana_Monkey_Suacoco" :["C5,C6,C7,C8"],
	"Lar_Gibbon_Tarzan" :["C9,C10,C11,C12"],
	"Guinea_Baboon_6" :["D1,D2,D3,D4"],
	"Chimpanzee_marlock" :["D5,D6,D7,D8,D9,D10"]
}

# names of the parse sublibraries
sublibraries = ["PARSE1_UDI_WT_1", "PARSE2_UDI_WT_2", "PARSE3_UDI_WT_3", "PARSE4_UDI_WT_4", "PARSE6_UDI_WT_4"]


#paths 
BASE_DIR = Path(__file__).parent

in_path = BASE_DIR / "../../backup/PrimaryData/human/ref/t2t"
out_path = BASE_DIR / "starsolo_v1_t2t"
out_path_parse = BASE_DIR / "parse_good_v1_t2t/"
fastq_files = BASE_DIR / "../emma/PARSE/rawdata/split-fq_v1"
path_to_fastq_files = BASE_DIR / "primary_data/"
empty_file = BASE_DIR / "gwf_out_empty/"

# Modify barcode_whitelist accordingly to the chemistry version of parse
barcode_whitelist = [
    BASE_DIR / "data/barcode_whitelists/bc_data_v1.csv",
    BASE_DIR / "data/barcode_whitelists/bc_data_v1.csv",
    BASE_DIR / "data/barcode_whitelists/bc_data_v2.csv"
]
path_to_data = BASE_DIR / "data"
path_to_split = BASE_DIR / "scripts/fastq_sep_groups_v0.5.py"
path_to_plots = BASE_DIR / "plots"
output_path_ref = BASE_DIR / "data/genomes/hg38_t2t" # output path for reference genome in parse pipeline


# Convert paths to absolute paths of your directory and to strings
in_path = str(in_path.resolve())
out_path = str(out_path.resolve())
out_path_parse = str(out_path_parse.resolve())
fastq_files = str(fastq_files.resolve())
path_to_fastq_files = str(path_to_fastq_files.resolve())
empty_file = str(empty_file.resolve())
barcode_whitelist = [str(path.resolve()) for path in barcode_whitelist]
path_to_data = str(path_to_data.resolve())
path_to_split = str(path_to_split.resolve())
path_to_plots = str(path_to_plots.resolve())
output_path_ref = str(output_path_ref.resolve())


###################
# FUNCTIONS #######
###################


# PARSE 
###################

## Function to split the fastqfiles by each well (specify the samples and wells in --group)
def split_fq_by__well(path_to_script,chemistry,path_to_fq,out_path,sublib,sample,group,empty_file):
	inputs = []
	outputs = [f'{empty_file}{sample}.txt',f'{out_path}{sublib}_group_{sample}_R1.fastq.gz',f'{out_path}{sublib}_group_{sample}_R2.fastq.gz']
	options = {"memory": "50g","account":"testis_singlecell","walltime":"23:00:00"}
	spec='''

	python {} \
		--chemistry {} \
        --kit WT \
        --kit_score_skip \
		--fq1 {}{}_R1.fq.gz \
		--fq2 {}{}_R2.fq.gz \
		--opath {} \
		--group {} {}
            
    touch {}

	'''.format(path_to_script,chemistry,path_to_fq,sublib,path_to_fq,sublib,out_path,sample,group,empty_file)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



## Construct reference genome in parse format
def reference_parse(path_to_fasta,path_to_genes,out_dir,empty_file2,empty_file1,sample):
	inputs = [empty_file1]
	outputs = [empty_file2]
	options = {"memory": "50g","account":"testis_singlecell","walltime":"23:00:00"}
	spec='''

    split-pipe \
        --mode mkref \
        --genome_name hg38 \
        --fasta {} \
        --genes {} \
        --output_dir {}

    touch {}

	'''.format(path_to_fasta,path_to_genes,out_dir,empty_file2,empty_file1,sample)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


## Run the parse pipeline
def spipe_parse(chemistry,path_to_genome,path_f1,path_f2,out_dir,empty_file2,sublib,sample,empty_file1):
	inputs = [f'{path_f1}{sublib}_group_{sample}_R1.fastq.gz',f'{path_f1}{sublib}_group_{sample}_R2.fastq.gz',empty_file1] #empty_file1 CHANGE THIS
	outputs = [empty_file2]
	options = {"memory": "300g","account":"testis_singlecell","walltime":"15:00:00","cores":24}
	spec='''

    split-pipe \
        --mode all \
        --chemistry {} \
        --kit WT \
        --kit_score_skip \
        --genome_dir {} \
        --fq1 {}{}_group_{}_R1.fastq.gz \
        --fq2 {}{}_group_{}_R2.fastq.gz \
        --output_dir {} \
        --sample Chimpanzee_Stephan A1-A4,D11-D12 \
        --sample Guinea_Baboon_5 A5-A8 \
        --sample Siamang_Gibbon_Coco A9-A12 \
        --sample Gelada_Baboon_1 B1-B4 \
        --sample Brown_Wooly_Monkey_1 B5-B8 \
        --sample Pigtailed_macaque_1 B9-B12 \
        --sample Whitecheeked_gibbon C1-C4 \
        --sample Diana_Monkey_Suacoco C5-C8 \
        --sample Lar_Gibbon_Tarzan C9-C12 \
        --sample Guinea_Baboon_6 D1-D4 \
        --sample Chimpanzee_marlock D5-D8,D9-D10

    touch {}

	'''.format(chemistry,path_to_genome,path_f1,sublib,sample,path_f2,sublib,sample,out_dir,empty_file2,sublib,sample,empty_file1)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


## Function to get the combined fastqfiles (R1+R2)
def spipe_parse_fastq(chemistry,path_to_genome,path_f1,path_f2,out_dir,empty_file2,sublib,sample,empty_file1):
	inputs = [] 
	outputs = [empty_file2]
	options = {"memory": "100g","account":"testis_singlecell","walltime":"15:00:00","cores":4}
	spec='''

    split-pipe \
        --mode pre \
        --chemistry {} \
        --kit WT \
        --kit_score_skip \
		--keep_temps \
		--genome_dir {} \
        --fq1 {}{}_group_{}_R1.fastq.gz \
        --fq2 {}{}_group_{}_R2.fastq.gz \
        --output_dir {}

    touch {}

	'''.format(chemistry,path_to_genome,path_f1,sublib,sample,path_f2,sublib,sample,out_dir,empty_file2,sublib,sample,empty_file1)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

## Combine alignments of different sublibraries for one sample
def comb_parse(sublib1,sublib2,sublib3,sublib4,sublib5,out_dir,empty_file2,empty_file1):
	inputs = [empty_file1]
	outputs = [empty_file2]
	options = {"memory": "80g","account":"testis_singlecell","walltime":"5:00:00"}
	spec='''

    split-pipe \
        --mode comb \
        --sublibraries {} \
             {} \
             {} \
             {} \
             {} \
        --output_dir {}

    touch {}

	'''.format(sublib1,sublib2,sublib3,sublib4,sublib5,out_dir,empty_file2,empty_file1)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# STARSOLO 
###################

## Construct reference genome to run starsolo pipeline
def reference_starsolo(path_to_ref,genome_fa,gene_annotation,empty_file):
	inputs = [empty_file]
	outputs = [f'{path_to_ref}Log.out']
	options = {"memory": "100g","account":"testis_singlecell","walltime":"05:00:00"}
	spec='''

	cd {} 
	STAR  --runMode genomeGenerate --runThreadN 4 --genomeSAsparseD 3 --genomeDir {} --genomeFastaFiles {}  --sjdbGTFfile {}	

	'''.format(path_to_ref,path_to_ref,genome_fa,gene_annotation,empty_file)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Alignment with starsolo
def mapping_starsolo(path_to_ref, read_2, read_1, CB_whitelist, out_prefix):
	inputs = [f'{path_to_ref}Log.out']
	outputs = [f'{out_prefix}Log.final.out']
	options = {"memory": "200g","walltime":"35:00:00", "account":"testis_singlecell"}
	spec='''

    STAR --runThreadN 16 \
            --genomeDir {} \
            --readFilesIn <(zcat {}) <(zcat {}) \
            --soloType CB_UMI_Complex \
            --outSAMtype BAM SortedByCoordinate \
            --soloCBwhitelist {} \
			--soloCBmatchWLtype EditDist_2 \
            --soloBarcodeReadLength 0 \
            --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
            --soloUMIposition 0_0_0_9 \
            --outFilterScoreMin 30 \
            --outFileNamePrefix {} \
            --soloUMIfiltering MultiGeneUMI \
            --soloUMIdedup 1MM_All \
            --soloMultiMappers EM \
            --soloFeatures Gene GeneFull Velocyto \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
            --soloCellFilter EmptyDrops_CR \
            --soloCellReadStats Standard


	'''.format(path_to_ref, read_1, read_2, CB_whitelist, out_prefix)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# MISC
###################

# function that concatenates fastqfiles from the same sublibrary
def concat_species(in_dir, fastqspecies, out_file):
	inputs = [] 
	outputs = [] #
	options = {"memory": "300g","walltime":"10:00:00", "account":"testis_singlecell"}
	spec='''

	cat {} > {}

	'''.format(fastqspecies, out_file, in_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# function to filter with emptydrops the raw file from starsolo
def filter_emptydrops(path_to_inp, path_out):
	inputs = [] #path_to_inp+'Log.out'
	outputs = [] #f'{path_out}Log.final.out'
	options = {"memory": "50g","walltime":"35:00:00", "account":"testis_singlecell"}
	spec='''

	STAR --runMode soloCellFiltering {} {} --soloCellFilter EmptyDrops_CR

	'''.format(path_to_inp, path_out)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# function to visualize knee plot after empty drops
def plot_combined_knee(matrix_raw, matrix_fil, plot_title = None):
    """
    Function for plotting a combined knee plot of UMI counts vs. rank for two matrices.
    
    Parameters:
    - matrix_raw: AnnData object containing raw single-cell RNA-seq data.
    - matrix_fil: AnnData object containing filtered single-cell RNA-seq data.
    - save_path: Optional path to save the plot. If None, the plot is not saved.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import os

    # Calculate total UMI counts per cell for each matrix
    total_counts_raw = np.array(matrix_raw.X.sum(axis=0)).flatten()  # Convert to 1D array
    total_counts_fil = np.array(matrix_fil.X.sum(axis=0)).flatten()  # Convert to 1D array

    # Sort the total UMI counts in descending order
    sorted_umis_raw = np.sort(total_counts_raw)[::-1]
    sorted_umis_fil = np.sort(total_counts_fil)[::-1]
    
    # Create rank arrays
    ranks_raw = np.arange(1, len(sorted_umis_raw) + 1)
    ranks_fil = np.arange(1, len(sorted_umis_fil) + 1)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(ranks_raw, sorted_umis_raw, label='Raw Data', color='blue')
    plt.plot(ranks_fil, sorted_umis_fil, label='Filtered Data', color='red')
    plt.title(plot_title if plot_title else 'Combined Knee Plot: UMI Counts vs. Rank (Log-Log Scale)')
    plt.xlabel('Rank')
    plt.ylabel('Total UMI Counts')
    plt.xscale('log')  # Log scale for x-axis
    plt.yscale('log')  # Log scale for y-axis
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


def extract_combined_dir_name(path):
    # Split the path into components
    path_parts = path.split(os.sep)
    
    # Get the parts we need (the second last and third last directory)
    species_code = path_parts[-4]  # e.g., 'SGB'
    species_name = path_parts[-3]  # e.g., 'Siamang_Gibbon_Coco'

    # Join with underscore
    combined_dir_name = f"{species_code}_{species_name}"
    
    return combined_dir_name

###################
# PREPROCESSING #
###################

# # Barcode whitelist processing: creates a list with the paths to the 3 rounds of barcoding
barcode_whitelist_use = []
for csv_file in barcode_whitelist:
	df = pd.read_csv(csv_file, delimiter=',')
	# extract the 'sequence' column
	sequence_column = df['sequence']
	path = csv_file.replace(".csv", ".txt")
	barcode_whitelist_use.append(path)
	sequence_column.to_csv(path, index=False, header=None)
barcodes = " ".join(barcode_whitelist_use)
# barcodes: list containing the 3 paths to the 3 sets of barcodes

######################################################################################################################################
# # Split the fastqfiles (R1 and R2) for each sublibrary per sample:

if not os.path.exists(empty_file):
    os.makedirs(empty_file) 


# It calls /scripts/fastq_sep_groups.py
empty_files = None
for sublib in sublibraries:
	for sample,group_list in sp_well_dict.items(): 
		group = group_list[0]
		empty_files = f'{empty_file}empty_gwf{sublib}'
		path_out = f'{fastq_files}/{sublib}/'
		if not os.path.exists(path_out):
			os.makedirs(path_out)
		gwf.target_from_template(f'split_fastq_{sublib}_{sample}',
            split_fq_by__well(path_to_split,chemistry,path_to_fastq_files,path_out,sublib,sample,group,empty_files))


################
# RUN PARSE
################

# ## make reference genome
path_to_fasta = f'{in_path}/genome.fa'
path_to_genes = f'{in_path}/gene_annot_orth_human_mc_ampl.gtf' # For t2t reference
empty_file2 = f'{empty_file}ref_parse.txt'

if not os.path.exists(output_path_ref):
	os.makedirs(output_path_ref) 

for sublib in sublibraries:
	for sample,group_list in sp_well_dict.items(): 
		group = group_list[0]
		empty_files = f'{empty_file}empty_gwf{sublib}{sample}.txt'
		sample = sample

#make reference
gwf.target_from_template(f'parse_ref_annotation',
	reference_parse(path_to_fasta,path_to_genes,output_path_ref,empty_file2,empty_files,sample))

# run parse pipeline
empty_files = None
for sp,name in sp_sample_dict.items():	
	matrix_counts = []
	samples = []
	for sample in sp_sample_dict[sp]:	
		out_path_sample = f'{out_path_parse}{sp}/{sample}/' 
		empty_files = None
		for sublib in sublibraries:
			empty_files = f'{empty_file}empty_gwf_spipe_{sample}_{sublib}.txt'
			matching_files = []
			fastq_files_good = f'{fastq_files}/{sublib}/'
			out_path_sample_sublib = f'{out_path_sample}{sublib}/'
			gwf.target_from_template(f'run_parse_{sublib}_{sample}',
				spipe_parse(chemistry,output_path_ref,fastq_files_good,fastq_files_good,out_path_sample_sublib,empty_files,sublib,sample,empty_file2))

# Combine sublibraries per sample
for sp,name in sp_sample_dict.items():	
	for sample in sp_sample_dict[sp]:	
		sublibrary_files = []
		empty_file1 = f'{empty_file}empty_gwf_combine_{sample}.txt'
		for sublib in sublibraries:
			out_path_sample = f'{out_path_parse}{sp}/{sample}/{sublib}' 
			out_dir = f'{out_path_parse}combined/{sp}/{sample}'
			sublibrary_files.append(out_path_sample)
		gwf.target_from_template(f'comb_parse_{sp}_{sample}',
					comb_parse(sublibrary_files[0],sublibrary_files[1],sublibrary_files[2],sublibrary_files[3],sublibrary_files[4],out_dir,empty_file1,empty_files))


## Don't run, just check for correct combination of R1 and R2 (check that reads are well and version is well used)
# for sp,name in sp_sample_dict.items():	
# 	matrix_counts = []
# 	samples = []
# 	for sample in sp_sample_dict[sp]:	
# 		out_path_sample = f'{out_path_parse}{sp}/{sample}/' 
# 		empty_files = None
# 		for sublib in sublibraries:
# 			empty_files = f'{empty_file}empty_gwf_spipe_fastq_{sample}_{sublib}.txt'
# 			matching_files = []
# 			fastq_files_good = f'{fastq_files}/{sublib}/'
# 			out_path_sample_sublib = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/combined_fastqfiles_v1_t2t/{sublib}/{sample}/'
# 			if not os.path.exists(out_path_sample_sublib):
# 				os.makedirs(out_path_sample_sublib)
# 			gwf.target_from_template(f'parse_fastq_{sublib}_{sample}',
# 				spipe_parse_fastq(chemistry,output_path_ref,fastq_files_good,fastq_files_good,out_path_sample_sublib,empty_files,sublib,sample,empty_file2))


# ################
# # RUN STARSOLO
# ################

# input files for reference
genome_file = f'{in_path}/genome.fa'
gene_annotation_file = f'{in_path}/gene_annot_orth_human_mc_ampl.gtf' # For t2t reference
path_to_ref = f'{out_path}/annotation_starsolo/'

# create the directory to store the reference if it doesn't exist already
if not os.path.exists(path_to_ref):
	os.makedirs(path_to_ref) 

empty_files = None
for sublib in sublibraries: 
    empty_files = f'{empty_file}/empty_gwf{sublib}'

# # construct STARsolo genome reference
gwf.target_from_template(f'run_starsolo_ref_annotation',
	reference_starsolo(path_to_ref,genome_file,gene_annotation_file,empty_files))

# Run starsolo for each species and each sublibrary
for sp,name in sp_sample_dict.items():	
	matrix_counts = []
	samples = []
	for sample in sp_sample_dict[sp]:	
		out_path_sample = f'{out_path}{sp}/{sample}/' 
		
		for sublib in sublibraries:
		# get the specific R1 and R2 fastqfiles for each species
			matching_files = []
			fastq_files_good = f'{fastq_files}/{sublib}'
			out_path_sample_sublib = f'{out_path_sample}{sublib}/'
			for files in os.listdir(fastq_files_good):
				if fnmatch.fnmatch(files, f'*{sample}*'):
					path = os.path.join(fastq_files_good,files)
					matching_files.append(path)
					matching_files.sort()
					# Separate R1 and R2 
					R1s = [file for file in matching_files if '_R1' in file]
					R2s = [file for file in matching_files if '_R2' in file]
					# Converting lists to space-separated strings
					first_three_str = ' '.join(R1s)
					last_three_str = ' '.join(R2s)
		# print(last_three_str)
			# Run starsolo for each sublibrary and sample
			gwf.target_from_template(f'run_starsolo_{sublib}_{sample}',
				mapping_starsolo(path_to_ref, last_three_str, first_three_str, barcodes, out_path_sample_sublib))


# run analysis/combine_files.py to combine all the sublibraries for starsolo
# run analysis/after_mapping_tailored.sh to prepare the dataset
# run scripts/unite_tables.py to create all the summary files for starsolo

#########

# This small code chunk is to run emptydrops 
# don't run emptydrops if we are already taking the filtered matrices (which is this setup)
# for combined_dir in path_matrices.keys():
#     # input path
#     path_to_inp = f'{out_path}{combined_dir}' 

#     # output path
#     path_o = combined_dir.replace('combined', 'combined_filtered')
#     path_out = f'{out_path}{path_o}/'

#     name = combined_dir.replace('/', '_')

#     # Run filter_emptydrops function
#     gwf.target_from_template(f'run_emptydrops_{name}',
# 							 filter_emptydrops(path_to_inp, path_out))