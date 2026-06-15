#!/bin/bash
# FoldX RepairPDB for LPA structures
cd C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa

foldx --command RepairPDB --pdb LPA_wildtype.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa
foldx --command RepairPDB --pdb LPA_KIV78.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa
foldx --command RepairPDB --pdb LPA_KIV10_Protease.pdb --pdb-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa --output-dir C:/Users/nader/Downloads/calon_ukb_pipeline\alphafold\foldx\lpa
