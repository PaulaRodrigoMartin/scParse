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


###################
# INPUT DATA #
###################

#string indicating version of parse kit, can be v1, v2, v3
## parse_version = v2

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

# paths 
# in_path = f'/home/paulilokiestudia/testis_singlecell/backup/PrimaryData/human/ref/GRCh38.p14'
in_path = f'/home/paulilokiestudia/testis_singlecell/backup/PrimaryData/human/ref/t2t' ## modify this path to match yours, maybe copy and move the files to the newvolum
out_path = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/starsolo_v1_t2t/'
out_path_parse = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/parse_good_v1_t2t/' ##modify this to go back to normal
fastq_files = "/home/paulilokiestudia/testis_singlecell/Workspaces/emma/PARSE/rawdata/split-fq_v1" ##modify this to go back to normal
path_to_fastq_files = "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/primary_data/"
chemistry = "v1"
empty_file = "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/gwf_out_empty/"
# /home/paulilokiestudia/testis_singlecell/backup/PrimaryData/human/ref/t2t/gene_annot_orth_human_mc_ampl.gtf
# /home/paulilokiestudia/testis_singlecell/backup/PrimaryData/human/ref/t2t/genome.fa

# Modify barcode_whitelist accordingly to the version
barcode_whitelist = ["/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/barcode_whitelists/bc_data_v1.csv", "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/barcode_whitelists/bc_data_v1.csv", "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/barcode_whitelists/bc_data_v2.csv"]
path_to_data = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/'
path_to_split = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/scripts/fastq_sep_groups_v0.5.py'
path_to_plots = f'/home/paulilokiestudia/testis_singlecell/Workspaces/paula/plots/'
sublibraries = ["PARSE1_UDI_WT_1", "PARSE2_UDI_WT_2", "PARSE3_UDI_WT_3", "PARSE4_UDI_WT_4", "PARSE6_UDI_WT_4"]
# raw_data_folders = ["",]


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

# quitar touch y poner en outputs los nombres
		# --group Chimpanzee_Stephan A1-A4,D11-D12 \
		# --group Guinea_Baboon_5 A5-A8 \
		# --group	Siamang_Gibbon_Coco A9-A12 \
		# --group	Gelada_Baboon_1 B1-B4 \
		# --group	Brown_Wooly_Monkey_1 B5-B8 \
		# --group	Pigtailed_macaque_1 B9-B12 \
		# --group	Whitecheeked_gibbon C1-C4 \
		# --group	Diana_Monkey_Suacoco C5-C8 \
		# --group	Lar_Gibbon_Tarzan C9-C12 \
		# --group	Guinea_Baboon_6 D1-D4 \
		# --group	Chimpanzee_marlock D5-D10

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
	inputs = [] #empty_file1 CHANGE THIS
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
#       --clipAdapterType CellRanger4 \
#		--clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG \ delete it because it's not that one


# MISC
###################
# for postprocessing the files, gets the paths of the matrices, features and barcodes and stores them in a variable
def generate_paths(sp_sample_dict, sublibraries, out_path, file_structure, key_suffix):
    paths = {}
    for sp, name in sp_sample_dict.items():
        for sample in sp_sample_dict[sp]:
            key = f"{sp}/{sample}/{key_suffix}"
            paths[key] = []
            out_path_sample = f"{out_path}{sp}/{sample}/"
            for sublib in sublibraries:
                path = f"{out_path_sample}{sublib}/{file_structure}"
                paths[key].append(path)
    return paths

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

# Plot QC of data
def plotQC(adata, plot_title = None, output_file = None):
    """
    Function for plotting QC values. 
    Needs to run preprocess.calculateQC first.
    """
    from scipy.sparse import issparse
    import scanpy as sc
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import os


    df=pd.DataFrame(index=adata.obs_names)
    df['total_counts'] = adata.obs['total_counts']
    df['n_genes_per_cell'] = adata.obs['n_genes_by_counts']
    if 'percent_mito' in adata.obs_keys():
        df['percent_mito'] = adata.obs['percent_mito']
    else:
        df['percent_mito'] = np.zeros(adata.shape[0])

    plt.rcParams['figure.figsize']=(20,20)
    f, ax = plt.subplots(3,2)

    sns.scatterplot(x='total_counts', y='n_genes_per_cell', hue='percent_mito', data=df, ax=ax[0,0])
    ax[0,0].set_title('UMI vs GENES plot - percent mito genes')

    if 'prop_spl' in adata.obs_keys():
        df['rate_spliced'] = adata.obs['prop_spl']
    else:
        df['rate_spliced'] = np.zeros(adata.shape[0])
    sns.scatterplot(x='total_counts', y='n_genes_per_cell', hue='rate_spliced', data=df, ax=ax[0,1])
    ax[0,1].set_title('UMI vs GENES plot - spliced proportions')

    sns.distplot(df['total_counts'], ax=ax[1,0], kde = True)
    ax[1,0].set_title('UMI counts per cell')

    sns.distplot(df['n_genes_per_cell'], ax=ax[1,1], kde= True)
    ax[1,1].set_title('Genes per cell')

    df['counts_over_genes'] = df['total_counts']/df['n_genes_per_cell']
    sns.scatterplot(x='total_counts', y='counts_over_genes', hue='percent_mito', data=df, ax=ax[2,0])
    ax[2,0].set_title('UMI vs UMI/GENES ratio plot - percent MT genes')
    sns.scatterplot(x='total_counts', y='counts_over_genes', hue='rate_spliced', data=df, ax=ax[2,1])
    ax[2,1].set_title('UMI vs UMI/GENES ratio plot - spliced proportions')

    plt.suptitle(plot_title, fontsize=16)

    plt.rcParams['figure.figsize']=(6,6)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if output_file:
        plt.savefig(output_file, format='png', bbox_inches='tight')
    plt.close() 

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

# if not os.path.exists(empty_file):
#     os.makedirs(empty_file) 


# # Basically calls /scripts/fastq_sep_groups.py
# empty_files = None
# for sublib in sublibraries:
# 	for sample,group_list in sp_well_dict.items(): 
# 		group = group_list[0]
# 		empty_files = f'{empty_file}empty_gwf{sublib}'
# 		path_out = f'{fastq_files}/{sublib}/'
# 		if not os.path.exists(path_out):
# 			os.makedirs(path_out)
# 		gwf.target_from_template(f'split_fastq_{sublib}_{sample}',
#             split_fq_by__well(path_to_split,chemistry,path_to_fastq_files,path_out,sublib,sample,group,empty_files))


################
# RUN PARSE
################

# ## make reference genome

# path_to_fasta = f'{in_path}/genome.fa'
# # path_to_genes = f'{in_path}/genes_amplicons_starsolo.gtf' # For first used reference 
# path_to_genes = f'{in_path}/gene_annot_orth_human_mc_ampl.gtf' # For t2t reference
# output_path_ref = "/home/paulilokiestudia/testis_singlecell/Workspaces/paula/data/genomes/hg38_t2t"
# empty_file2 = f'{empty_file}ref_parse.txt'

# if not os.path.exists(output_path_ref):
# 	os.makedirs(output_path_ref) 

# for sublib in sublibraries:
# 	for sample,group_list in sp_well_dict.items(): 
# 		group = group_list[0]
# 		empty_files = f'{empty_file}empty_gwf{sublib}{sample}.txt'
# 		sample = sample

# gwf.target_from_template(f'parse_ref_annotation',
# 	reference_parse(path_to_fasta,path_to_genes,output_path_ref,empty_file2,empty_files,sample))

		
# empty_files = None
# for sp,name in sp_sample_dict.items():	
# 	matrix_counts = []
# 	samples = []
# 	for sample in sp_sample_dict[sp]:	
# 		out_path_sample = f'{out_path_parse}{sp}/{sample}/' 
# 		empty_files = None
# 		for sublib in sublibraries:
# 			empty_files = f'{empty_file}empty_gwf_spipe_{sample}_{sublib}.txt'
# 			matching_files = []
# 			fastq_files_good = f'{fastq_files}/{sublib}/'
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/'
# 			gwf.target_from_template(f'run_parse_{sublib}_{sample}',
# 				spipe_parse(chemistry,output_path_ref,fastq_files_good,fastq_files_good,out_path_sample_sublib,empty_files,sublib,sample,empty_file2))

# # # Combine sublibraries per sample
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:	
# 		sublibrary_files = []
# 		empty_file1 = f'{empty_file}empty_gwf_combine_{sample}.txt'
# 		for sublib in sublibraries:
# 			out_path_sample = f'{out_path_parse}{sp}/{sample}/{sublib}' 
# 			out_dir = f'{out_path_parse}combined/{sp}/{sample}'
# 			sublibrary_files.append(out_path_sample)
# 		gwf.target_from_template(f'comb_parse_{sp}_{sample}',
# 					comb_parse(sublibrary_files[0],sublibrary_files[1],sublibrary_files[2],sublibrary_files[3],sublibrary_files[4],out_dir,empty_file1,empty_files))


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
# gene_annotation_file = f'{in_path}/genes_amplicons_starsolo.gtf' # For first reference
gene_annotation_file = f'{in_path}/gene_annot_orth_human_mc_ampl.gtf' # For t2t reference
path_to_ref = f'{out_path}annotation_starsolo/'

# create the directory to store the reference if it doesn't exist already
if not os.path.exists(path_to_ref):
	os.makedirs(path_to_ref) 

empty_files = None
for sublib in sublibraries: 
    empty_files = f'{empty_file}empty_gwf{sublib}'

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

# For each species, unify the raw matrices
#############################################################

# # Define different file structures
# file_structures = {
#     "combined_fil_mtx": "Solo.out/GeneFull/filtered/matrix.mtx",
#     "combined_fil_tsv": "Solo.out/GeneFull/filtered/barcodes.tsv",
#     "combined_fil_features": "Solo.out/GeneFull/filtered/features.tsv",
#     "combined_tsv_raw": "Solo.out/GeneFull/raw/barcodes.tsv",
#     "combined_features_raw": "Solo.out/GeneFull/raw/features.tsv",
#     "velocyto_unspliced": "Solo.out/Velocyto/raw/unspliced.mtx",
#     "velocyto_spliced": "Solo.out/Velocyto/raw/spliced.mtx",
#     "velocyto_ambiguous": "Solo.out/Velocyto/raw/ambiguous.mtx",
#     "velocyto_tsvs": "Solo.out/Velocyto/raw/barcodes.tsv",
#     "velocyto_features": "Solo.out/Velocyto/raw/features.tsv",
#     "unique_and_multi": "Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx"
# }

# path_matrices = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["combined_fil_mtx"], "combined_fil")
# path_tsvs = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["combined_fil_tsv"], "combined_fil")
# path_features = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["combined_fil_features"], "combined_fil")
# path_tsvs_raw = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["combined_tsv_raw"], "combined")
# path_features_raw = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["combined_features_raw"], "combined")
# path_unspliced = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["velocyto_unspliced"], "Velocyto")
# path_spliced = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["velocyto_spliced"], "Velocyto")
# path_ambiguous = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["velocyto_ambiguous"], "Velocyto")
# path_tsvs_velocyto = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["velocyto_tsvs"], "Velocyto")
# path_features_velocyto = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["velocyto_features"], "Velocyto")
# path_unique_and_multi = generate_paths(sp_sample_dict, sublibraries, out_path, file_structures["unique_and_multi"], "combined")


# # Get the paths for the matrices for each sublibrary per sample (GeneFull)
# path_matrices = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_matrices[f'{sp}/{sample}/combined_fil'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/matrix.mtx'
# 			path_matrices[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

# path_tsvs = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_tsvs[f'{sp}/{sample}/combined_fil'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/barcodes.tsv'
# 			path_tsvs[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

# path_features = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_features[f'{sp}/{sample}/combined_fil'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/filtered/features.tsv'
# 			path_features[f'{sp}/{sample}/combined_fil'].append(out_path_sample_sublib)

# # # For the raw files
# path_tsvs_raw = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_tsvs_raw[f'{sp}/{sample}/combined'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/barcodes.tsv'
# 			path_tsvs_raw[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)

# path_features_raw = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_features_raw[f'{sp}/{sample}/combined'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/features.tsv'
# 			path_features_raw[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)

# # # ####################################33
# # # ######## do i change the folder names? not combined_fil but GeneFull/filtered/ 
# # # ######## probably that makes more sense
# # # # For velocyto raw, unspliced
# path_unspliced = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_unspliced[f'{sp}/{sample}/Velocyto'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/unspliced.mtx'
# 			path_unspliced[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# # # For velocyto raw, spliced
# path_spliced = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_spliced[f'{sp}/{sample}/Velocyto'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/spliced.mtx'
# 			path_spliced[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# # # For velocyto raw, ambiguous
# path_ambiguous = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_ambiguous[f'{sp}/{sample}/Velocyto'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/ambiguous.mtx'
# 			path_ambiguous[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# # ## tsvs and barcodes from velocyto
# path_tsvs_velocyto = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_tsvs_velocyto[f'{sp}/{sample}/Velocyto'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/barcodes.tsv'
# 			path_tsvs_velocyto[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# path_features_velocyto = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_features_velocyto[f'{sp}/{sample}/Velocyto'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/Velocyto/raw/features.tsv'
# 			path_features_velocyto[f'{sp}/{sample}/Velocyto'].append(out_path_sample_sublib)

# path_unique_and_multi = {}
# for sp,name in sp_sample_dict.items():	
# 	for sample in sp_sample_dict[sp]:
# 		path_unique_and_multi[f'{sp}/{sample}/combined'] = []
# 		out_path_sample = f'{out_path}{sp}/{sample}/' 
# 		for sublib in sublibraries:
# 			out_path_sample_sublib = f'{out_path_sample}{sublib}/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx'
# 			path_unique_and_multi[f'{sp}/{sample}/combined'].append(out_path_sample_sublib)


# # #############################################################
# # # # Add the sublibrary tag to each cell name
# for sublib in sublibraries:
# 	for key, values in path_tsvs.items():
# 		# print(f"{key}: {value}")
# 		for value in values:
# 			if sublib in value:
# 				df = pd.read_csv(value, sep='\t', header=None)
# 	            # Add the sublibrary to the cell name
# 				df[0] = df[0] + '_' + sublib
# 	            # Save the modified df back to a TSV
# 				df.to_csv(value, sep='\t', index=False, header=False)
# 				print(f"Processed {value} with sublibrary {sublib}")
# # print(path_tsvs_raw)
# for sublib in sublibraries:
# 	for key, values in path_tsvs_raw.items():
# 		# print(f"{key}: {value}")
# 		for value in values:
# 			if sublib in value:
# 				df = pd.read_csv(value, sep='\t', header=None)
# 	            # Add the sublibrary to the cell name
# 				df[0] = df[0] + '_' + sublib
# 	            # Save the modified df back to a TSV
# 				df.to_csv(value, sep='\t', index=False, header=False)
# 				print(f"Processed {value} with sublibrary {sublib}")

# for sublib in sublibraries:
# 	for key, values in path_tsvs_velocyto.items():
# 		# print(f"{key}: {value}")
# 		for value in values:
# 			if sublib in value:
# 				df = pd.read_csv(value, sep='\t', header=None)
# 	            # Add the sublibrary to the cell name
# 				df[0] = df[0] + '_' + sublib
# 	            # Save the modified df back to a TSV
# 				df.to_csv(value, sep='\t', index=False, header=False)
# 				print(f"Processed {value} with sublibrary {sublib}")

#############################################################
# # # # # # # ### Combine matrices
# combined_matrices = {}
# for key in sorted(path_matrices.keys()):
#     matrices = []  # list to hold all matrices for the current key

#     for file_path in path_matrices[key]:
#         matrix = mmread(file_path).tocsr()  # Read the matrix from the file
#         matrices.append(matrix)  # Add the matrix to the list
#     # print(matrices)
#     # Combine all matrices for the current key into one big matrix
#     big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
#     print(big_matrix.shape)
# 	# Create a directory for the current key to store the combined matrix
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined matrix
#     output_file_path = os.path.join(key_directory, f'matrix.mtx')
    
#     # Save the combined matrix to the file
#     mmwrite(output_file_path, big_matrix)
    
#     combined_matrices[key] = big_matrix
#     print(f"Combined matrix for {key} saved to {output_file_path}")

# # # for unspliced
# combined_matrices_unspliced = {}
# for key in sorted(path_unspliced.keys()):
#     matrices = []  # list to hold all matrices for the current key

#     for file_path in path_unspliced[key]:
#         matrix = mmread(file_path).tocsr()  # Read the matrix from the file
#         matrices.append(matrix)  # Add the matrix to the list
#     # print(matrices)
#     # Combine all matrices for the current key into one big matrix
#     big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
#     print(big_matrix.shape)
# 	# Create a directory for the current key to store the combined matrix
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined matrix
#     output_file_path = os.path.join(key_directory, f'unspliced.mtx')
    
#     # Save the combined matrix to the file
#     mmwrite(output_file_path, big_matrix)
    
#     combined_matrices_unspliced[key] = big_matrix
#     print(f"Combined matrix for {key} saved to {output_file_path}")

# # # for spliced
# combined_matrices_spliced = {}
# for key in sorted(path_spliced.keys()):
#     matrices = []  # list to hold all matrices for the current key

#     for file_path in path_spliced[key]:
#         matrix = mmread(file_path).tocsr()  # Read the matrix from the file
#         matrices.append(matrix)  # Add the matrix to the list
#     # print(matrices)
#     # Combine all matrices for the current key into one big matrix
#     big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
#     print(big_matrix.shape)
# 	# Create a directory for the current key to store the combined matrix
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined matrix
#     output_file_path = os.path.join(key_directory, f'spliced.mtx')
    
#     # Save the combined matrix to the file
#     mmwrite(output_file_path, big_matrix)
    
#     combined_matrices_spliced[key] = big_matrix
#     print(f"Combined matrix for {key} saved to {output_file_path}")


# # # for ambiguous
# combined_matrices_ambiguous = {}
# for key in sorted(path_ambiguous.keys()):
#     matrices = []  # list to hold all matrices for the current key

#     for file_path in path_ambiguous[key]:
#         matrix = mmread(file_path).tocsr()  # Read the matrix from the file
#         matrices.append(matrix)  # Add the matrix to the list
#     # print(matrices)
#     # Combine all matrices for the current key into one big matrix
#     big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
#     print(big_matrix.shape)
# 	# Create a directory for the current key to store the combined matrix
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined matrix
#     output_file_path = os.path.join(key_directory, f'ambiguous.mtx')
    
#     # Save the combined matrix to the file
#     mmwrite(output_file_path, big_matrix)
    
#     combined_matrices_ambiguous[key] = big_matrix
#     print(f"Combined matrix for {key} saved to {output_file_path}")

# # # for raw/UniqueandMult-EM.mtx
# combined_matrices_mult = {}
# for key in sorted(path_unique_and_multi.keys()):
#     matrices = []  # list to hold all matrices for the current key

#     for file_path in path_unique_and_multi[key]:
#         matrix = mmread(file_path).tocsr()  # Read the matrix from the file
#         matrices.append(matrix)  # Add the matrix to the list
#     # print(matrices)
#     # Combine all matrices for the current key into one big matrix
#     big_matrix = hstack(matrices)  # Change to hstack(matrices) if horizontal stacking is needed
#     print(big_matrix.shape)
# 	# Create a directory for the current key to store the combined matrix
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined matrix
#     output_file_path = os.path.join(key_directory, f'UniqueAndMult-EM.mtx')
    
#     # Save the combined matrix to the file
#     mmwrite(output_file_path, big_matrix)
    
#     combined_matrices_mult[key] = big_matrix
#     print(f"Combined matrix for {key} saved to {output_file_path}")

# # # # ## Combine tsvs
# # ## filtered matrices
# for key in sorted(path_tsvs.keys()):
#     dfs = []
#     # Read each TSV into a df and append to list
#     for file_path in path_tsvs[key]:
#         df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
#         dfs.append(df)  # Add the DataFrame to the list
    
#     # Concatenate all dfs for the current key into one big df
#     combined_df = pd.concat(dfs, ignore_index=True)
    
#     # Create a directory for the current key to store the combined TSV
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined TSV
#     output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
#     # Save the combined DataFrame to a TSV file
#     combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
#     print(f"Combined TSV for {key} saved to {output_file_path}")

# # # for raw matrix
# for key in sorted(path_tsvs_raw.keys()):
#     dfs = []
#     # Read each TSV into a df and append to list
#     for file_path in path_tsvs_raw[key]:
#         df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
#         dfs.append(df)  # Add the DataFrame to the list
    
#     # Concatenate all dfs for the current key into one big df
#     combined_df = pd.concat(dfs, ignore_index=True)
    
#     # Create a directory for the current key to store the combined TSV
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined TSV
#     output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
#     # Save the combined DataFrame to a TSV file
#     combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
#     print(f"Combined TSV for {key} saved to {output_file_path}")

# # # raw/velocyto
# for key in sorted(path_tsvs_velocyto.keys()):
#     dfs = []
#     # Read each TSV into a df and append to list
#     for file_path in path_tsvs_velocyto[key]:
#         df = pd.read_csv(file_path, sep='\t', header=None)  # Read the TSV file
#         dfs.append(df)  # Add the DataFrame to the list
    
#     # Concatenate all dfs for the current key into one big df
#     combined_df = pd.concat(dfs, ignore_index=True)
    
#     # Create a directory for the current key to store the combined TSV
#     key_directory = os.path.join(out_path, key)
#     os.makedirs(key_directory, exist_ok=True)
    
#     # Define the output file path for the combined TSV
#     output_file_path = os.path.join(key_directory, 'barcodes.tsv')
    
#     # Save the combined DataFrame to a TSV file
#     combined_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    
#     print(f"Combined TSV for {key} saved to {output_file_path}")

# # # # ## Add features (genes) to combined folder
# for key, feature_paths in path_features.items():
#     # Create the combined directory path
#     combined_dir = os.path.join(out_path, key)
    
#     # Ensure the combined directory exists
#     os.makedirs(combined_dir, exist_ok=True)
    
#     # Copy each features.tsv file to the combined directory
#     for feature_path in feature_paths:
#         if os.path.isfile(feature_path):  # Check if the source file exists
#             destination_path = os.path.join(combined_dir, 'features.tsv')
#             shutil.copy2(feature_path, destination_path)
#             print(f"Copied {feature_path} to {destination_path}")
#         else:
#             print(f"File not found: {feature_path}")

# # # for raw matrix
# for key, feature_paths in path_features_raw.items():
#     # Create the combined directory path
#     combined_dir = os.path.join(out_path, key)
    
#     # Ensure the combined directory exists
#     os.makedirs(combined_dir, exist_ok=True)
    
#     # Copy each features.tsv file to the combined directory
#     for feature_path in feature_paths:
#         if os.path.isfile(feature_path):  # Check if the source file exists
#             destination_path = os.path.join(combined_dir, 'features.tsv')
#             shutil.copy2(feature_path, destination_path)
#             print(f"Copied {feature_path} to {destination_path}")
#         else:
#             print(f"File not found: {feature_path}")

# for key, feature_paths in path_features_velocyto.items():
#     # Create the combined directory path
#     combined_dir = os.path.join(out_path, key)
    
#     # Ensure the combined directory exists
#     os.makedirs(combined_dir, exist_ok=True)
    
#     # Copy each features.tsv file to the combined directory
#     for feature_path in feature_paths:
#         if os.path.isfile(feature_path):  # Check if the source file exists
#             destination_path = os.path.join(combined_dir, 'features.tsv')
#             shutil.copy2(feature_path, destination_path)
#             print(f"Copied {feature_path} to {destination_path}")
#         else:
#             print(f"File not found: {feature_path}")

##################### 
# run scripts/unite_tables.py to create all the summary files for starsolo

# # #############################################################################################################################################
############################################################################################################################################

# Run emptydrops --> don't run emptydrops if we are already taking the filtered matrices
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