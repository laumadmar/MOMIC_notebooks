# MOMIC: A Multi-Omics Pipeline for data analysis, integration and interpretation

MOMIC offers a complete analysis environment for analysing and integrating multi-omics data in a single, easy-to-use platform. 

MOMIC currently compiles protocols for whole genome SNP data (GWAS), mRNA expression (both from arrays and from RNAseq experiments) and protein data. Along with enrichment analysis and methods for combining distinct data. The proposed protocols are developed as Jupyter notebooks that guide the user through the tasks of pre-processing and transforming the data and performing the actual analysis, allowing the user to modify any piece of code needed along the process to adequate it to each project.

It is distributed as a docker project and a collection of Jupyter notebooks. Prior to getting the notebooks you need to install the MOMIC server locally. Follow the instrucctions on [this repo](https://github.com/laumadmar/MOMIC_server.git) or check out the [user manual](https://laumadmar.github.io/MOMIC_server).

## MOMIC jupyter Notebooks

This repository contains the collection of MOMIC notebooks. There are two alternatives to get it:

a) clone this repository from within your local MOMIC server, Jupyterhub instance. Once logged in, go to the git tab located on the left menu. Click on the button ‘Clone a Repository’ and provide the url for this repo.

b) clone the repository from the terminal, either using Jupyter terminal or via ssh into the container. CD into your jupyter home directory and type `git clone https://github.com/laumadmar/MOMIC_notebooks.git`.

You should have in your home directory a copy of all the notebooks necessary to carry out the analysis presented in MOMIC.