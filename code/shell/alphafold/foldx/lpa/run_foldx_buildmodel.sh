#!/bin/bash
# FoldX BuildModel for LPA variants
cd C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa

# STEP 1: RepairPDB first
foldx --command RepairPDB --pdb LPA_wildtype.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa
foldx --command RepairPDB --pdb LPA_KIV78.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa
foldx --command RepairPDB --pdb LPA_KIV10_Protease.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa

# STEP 2: BuildModel
foldx --command BuildModel --pdb LPA_KIV10_Protease_Repair.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --mutant-file LPA_KIV10_Protease_individual_list.txt --numberOfRuns 5
foldx --command BuildModel --pdb LPA_KIV78_Repair.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --mutant-file LPA_KIV78_individual_list.txt --numberOfRuns 5
foldx --command BuildModel --pdb LPA_wildtype_Repair.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --mutant-file LPA_wildtype_individual_list.txt --numberOfRuns 5
