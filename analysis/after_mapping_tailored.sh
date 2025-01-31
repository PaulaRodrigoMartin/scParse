#!/bin/bash

###############################
## MAPPING TO STARSOLO
#  After mapping with starsolo, prepare data to work with scanpy
#  Modify species, samples and paths and
#  Run ./ in terminal
###############################

## species and samples
declare -A get_samples=( ["CHMP"]="Chimpanzee_Stephan Chimpanzee_marlock" ["GUB"]="Guinea_Baboon_5 Guinea_Baboon_6" ["SGB"]="Siamang_Gibbon_Coco" ["GLB"]="Gelada_Baboon_1" ["BWM"]="Brown_Wooly_Monkey_1" ["PTM"]="Pigtailed_macaque_1" ["WCG"]="Whitecheeked_gibbon" ["DIM"]="Diana_Monkey_Suacoco" ["LRG"]="Lar_Gibbon_Tarzan")
path="/home/paulilokiestudia/testis_singlecell/Workspaces/paula/starsolo_v1/"
# CHMP GLB BWM PTM WCG 
for sp in DIM LRG GUB SGB
do

        for sample in ${get_samples[${sp}]}
        do
                # compress files
                gzip ${path}${sp}/${sample}/combined_fil/matrix.mtx #change to combined_filtered if you want to use the matrices combined and then filtered with empty drops instead of 
                gzip ${path}${sp}/${sample}/combined_fil/features.tsv
                gzip ${path}${sp}/${sample}/combined_fil/barcodes.tsv
                cp ${path}${sp}/${sample}/combined/UniqueAndMult-EM.mtx ${path}${sp}/${sample}/combined/matrix.mtx
                # gunzip ${path}${sp}/${sample}/combined/matrix.mtx #change to combined_filtered if you want to use the matrices combined and then filtered with empty drops instead of 
                # gunzip ${path}${sp}/${sample}/combined/features.tsv
                # gunzip ${path}${sp}/${sample}/combined/barcodes.tsv
                # #gunzip ${path}${sp}/${sample}/combined/matrix.mtx
                # # make unspliced folder
                mkdir -p ${path}${sp}/${sample}/Velocyto/raw/unspliced/
                cp ${path}${sp}/${sample}/Velocyto/unspliced.mtx ${path}${sp}/${sample}/Velocyto/raw/unspliced/matrix.mtx
                cp ${path}${sp}/${sample}/Velocyto/features.tsv ${path}${sp}/${sample}/Velocyto/raw/unspliced/
                cp ${path}${sp}/${sample}/Velocyto/barcodes.tsv ${path}${sp}/${sample}/Velocyto/raw/unspliced/
                # # make spliced folder
                mkdir -p ${path}${sp}/${sample}/Velocyto/raw/spliced/
                cp ${path}${sp}/${sample}/Velocyto/spliced.mtx ${path}${sp}/${sample}/Velocyto/raw/spliced/matrix.mtx
                cp ${path}${sp}/${sample}/Velocyto/features.tsv ${path}${sp}/${sample}/Velocyto/raw/spliced/
                cp ${path}${sp}/${sample}/Velocyto/barcodes.tsv ${path}${sp}/${sample}/Velocyto/raw/spliced/
                # # make ambiguous folder
                mkdir -p ${path}${sp}/${sample}/Velocyto/raw/ambiguous/
                cp ${path}${sp}/${sample}/Velocyto/ambiguous.mtx ${path}${sp}/${sample}/Velocyto/raw/ambiguous/matrix.mtx
                cp ${path}${sp}/${sample}/Velocyto/features.tsv ${path}${sp}/${sample}/Velocyto/raw/ambiguous/
                cp ${path}${sp}/${sample}/Velocyto/barcodes.tsv ${path}${sp}/${sample}/Velocyto/raw/ambiguous/
                # # compress all files in the new folders
                gzip ${path}${sp}/${sample}/Velocyto/raw/unspliced/*
                gzip ${path}${sp}/${sample}/Velocyto/raw/spliced/*
                gzip ${path}${sp}/${sample}/Velocyto/raw/ambiguous/*
                # # make folder for files to scanpy
                mkdir -p ${path}${sp}/${sample}/GeneFull/raw/to_scanpy/ 
                cp ${path}${sp}/${sample}/combined/UniqueAndMult-EM.mtx ${path}${sp}/${sample}/GeneFull/raw/to_scanpy/matrix.mtx
                cp ${path}${sp}/${sample}/combined/features.tsv ${path}${sp}/${sample}/GeneFull/raw/to_scanpy/
                cp ${path}${sp}/${sample}/combined/barcodes.tsv ${path}${sp}/${sample}/GeneFull/raw/to_scanpy/
                gzip ${path}${sp}/${sample}/GeneFull/raw/to_scanpy/*

        done
done

