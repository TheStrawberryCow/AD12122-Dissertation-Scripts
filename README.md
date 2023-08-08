# AD12122 Dissertation Scripts for Metagenomic Analyses

Dissertation Title:  Exploring the Link between the Oral Microbiome and Infective Endocarditis: A Systems Biology Approach

This repository contains a collection of R scripts developed specifically for metagenomic analysis. Each script focuses on a different aspect of the analysis process, from data processing to visualization. 

## Table of Contents

- [Description](#description)
- [Data Requirements](#data-requirements)
- [Usage](#usage)
- [Files in the Repository](#files-in-the-repository)
- [Contribute](#contribute)
- [Acknowledgements](#acknowledgements)

## Description

The AD12122-RScript repository has been crafted to facilitate comprehensive metagenomic analyses using R. Whether you're looking to analyze diversity patterns, differential abundance, or build predictive models, this collection has something to offer.

## Data Requirements

- Data input for these scripts should contain **read counts** for the metagenomic samples.
- It's crucial that the data is correctly formatted and cleaned for optimal results. Ensure there are no missing values or anomalies before processing.

## Usage

1. Ensure you've installed the required packages. This can be done using the `Packages.R` script.
2. Depending on your analysis objective, execute the corresponding R script. For instance, if you wish to analyze alpha diversity patterns, use the `Alpha Diversity.R` script.
3. Always refer to inline comments within each script for specific instructions or notes.

## Files in the Repository

- **Packages.R**: Lists and installs the required R packages for all scripts in this repository.
- **Phyloseq Object.R**: Procedures to create and manipulate a `phyloseq` object, which is often used in metagenomic analyses.
- **Relative Abundance.R**: Provides preliminary visual insights into your data, helping you get a quick overview.
- **Alpha Diversity.R**: Analyzes and visualizes patterns of alpha diversity in your dataset.
- **Beta Diversity.R**: Provides tools for assessing beta diversity across samples or groups.
- **Differential Abundance.R**: Helps in identifying taxa that are differentially abundant between groups.
- **Modelling.R**: Contains tools and functions for building predictive models based on your metagenomic data.

## Contribute

We appreciate any contributions to improve these scripts or extend their functionality. If you wish to contribute, please:

1. Fork the repository.
2. Create a new branch for your feature or fix.
3. Commit your changes.
4. Push to your branch.
5. Open a pull request.

## Acknowledgements

I would like to express my heartfelt thanks to my supervisors, Dr Vanessa Sousa and Joseph Falconer. Their guidance, patience, and expertise have been essential to the completion of this study. I am also grateful to Natasha Oh, whose collaboration and efforts have significantly contributed to this research.
