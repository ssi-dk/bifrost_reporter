#!/usr/bin/env python

import importlib
import importlib.util
import logging
import pandas as pd
from bson import ObjectId
import yaml
import numpy as np
import warnings
import os
import envyaml
import re

# This script includes helper functions to parse and normalize various YAML outputs
# from a bioinformatics pipeline such as Bifrost. These include:
# - MLST (Multi-Locus Sequence Typing)
# - Point mutation detection
# - AMR (Antimicrobial Resistance)
# - Plasmid and virulence detection
# - Species classification
# - Assembly quality metrics
# - QC stamping for pipeline success/failure


# ---------------
# CONFIG LOADING
# ---------------
# Define the core package name used for relative configuration paths
PACKAGE_NAME: str = "bifrost_reporter"

# Dynamically locate the package and its directory for loading config files
try:
    spec = importlib.util.find_spec(PACKAGE_NAME)
    if spec is None:
        raise ModuleNotFoundError(f"Package '{PACKAGE_NAME}' not found.")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    PACKAGE_DIR = os.path.dirname(module.__file__)  # Extract directory location

except ModuleNotFoundError as e:
    print(f"Error: {e}")
    PACKAGE_DIR = None
except AttributeError:
    print(f"Error: Could not determine package directory for '{PACKAGE_NAME}'.")
    PACKAGE_DIR = None
except Exception as e:
    print(f"Unexpected error: {e}")
    PACKAGE_DIR = None



def get_config(config_path: str = None):
    """
    Load specified YAML configuration file. If the path is None falls back to 
    the default config in the package directory 

    Parameters:
    ----------
    config_path : str, optional
        Path to the YAML config file. If None, defaults to 'config.default.yaml' in PACKAGE_DIR.

    Returns:
    -------
    dict
        Configuration parameters as a dictionary.
    """
    if config_path is None:
        config_path = os.path.join(PACKAGE_DIR, "config", "config.default.yaml")

    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except Exception as e:
        raise RuntimeError(f"Failed to load config file: {config_path}. Error: {str(e)}")


def extract_prefix(sample_name):
    """
    Extracts the first three underscore-separated components from a sample name.

    Parameters:
    ----------
    sample_name : str

    Returns:
    -------
    str
        Prefix composed of first 3 elements (e.g., "HER_BTP_WGS")
    """
    return "_".join(sample_name.split("_")[:3])


# ----------------------
# CUSTOM YAML HANDLERS
# ----------------------
# Enable parsing of BSON ObjectIds (common in MongoDB-based YAMLs)
def bson_objectid_constructor(loader, node):
    value = loader.construct_scalar(node)
    return ObjectId(value)

# Register the constructor for custom YAML tags
yaml.add_constructor('!bson.objectid.ObjectId', bson_objectid_constructor)





# ------------------------------------
# PARSERS FOR DIFFERENT BIFROST TOOLS
# ------------------------------------

def parse_mlst(list_files):
    """
    Parse MLST YAML files, extracting the 7 loci alleles and ST (sequence type).
    Returns a DataFrame with one row per sample.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            temp = yaml.load(f, Loader=yaml.Loader)
            if temp["status"] == "Success":
                d[sample_id] = temp["summary"]["mlst_report"]
            else:
                d[sample_id] = "N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A"

    # Split CSV-like MLST reports into list of 8 values (ST + 7 loci)
    for key, value in d.items():
        d[key] = value.split(',')
    
    df = pd.DataFrame.from_dict(d, orient='index', dtype=str)
    df.columns = range(df.shape[1])
    
    return df



def parse_kmapointmutations(list_files):
    """
    Extract point mutation data from KMA-based pointmutations_tsv reports.
    Returns a multi-indexed DataFrame with mutation records per sample.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("results", {}).get("pointmutations_tsv", {}).get("values", [])
                if not summary:
                    warnings.warn(f"Missing pointmutation summary for: {sample_id}")
                    continue
                df = pd.DataFrame.from_dict(summary)
                d[sample_id] = df
    if d:
        df = pd.concat(d, names=['Sample Name']).reset_index(level=0)
        df = df.set_index("Sample Name").drop(columns=["#Sample"])
    else:
        df = pd.DataFrame()
    return df



def check_stampers(list_files):
    """
    Evaluate the quality control (QC) summary stamps from YAML.
    If all QC checks are 'pass', mark overall as 'Pass'.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("results", {})
                df = pd.DataFrame.from_dict(summary)
                d[sample_id] = "Pass" if df["status"].eq("pass").all() else "Fail"
            else:
                d[sample_id] = "Requirement Not Met"
    return pd.DataFrame.from_dict(d, orient='index')



def parse_amrfinder(list_files):
    """
    Parse AMRFinder tool output from YAMLs, extracting AMR gene hits.
    Returns long-form DataFrame of gene hits per sample.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("summary", {}).get("output_tsv", [])
                df = pd.DataFrame(summary)
            else:
                df = pd.DataFrame(np.nan, index=range(1), columns=[
                    '% Coverage of reference sequence', '% Identity to reference sequence',
                    'Accession of closest sequence', 'Alignment length', 'Class', 'Contig id',
                    'Element subtype', 'Element type', 'Gene symbol', 'HMM description', 'HMM id',
                    'Method', 'Name of closest sequence', 'Protein identifier',
                    'Reference sequence length', 'Scope', 'Sequence name', 'Start', 'Stop',
                    'Strand', 'Subclass', 'Target length'])
            d[sample_id] = df
    df = pd.concat(d, names=['Sample Name']).reset_index(level=0).set_index("Sample Name")
    return df



def parse_species(list_files):
    """
    Parse Kraken-style species classification output.
    Adds custom column to compute unclassified + top-species proportion.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("summary", {})
                d[sample_id] = summary
    df = pd.DataFrame.from_dict(d, orient='index')
    df['sum_unclassified_species1'] = df['percent_unclassified'] + df['percent_classified_species_1']
    return df[[
        "name_classified_species_1", "percent_classified_species_1",
        "name_classified_species_2", "percent_classified_species_2",
        'percent_unclassified', 'sum_unclassified_species1']]



def parse_finder_tools(list_files, ariba_type):
    """
    Generic parser for ARIBA-style outputs (e.g., virulencefinder, plasmidfinder).
    Dynamically adapts to the provided ARIBA result type.
    """
    data_df = {}

    def extract_data(info):
        extracted_data = []
        if info:
            for entry in info:
                extracted_data.append({
                    'GENE': entry.get('GENE', np.nan),
                    '%COVERAGE': entry.get('%COVERAGE', np.nan),
                    '%IDENTITY': entry.get('%IDENTITY', np.nan),
                    'SEQUENCE': entry.get('SEQUENCE', np.nan),
                    'START': entry.get('START', np.nan),
                    'END': entry.get('END', np.nan),
                    'DATABASE': entry.get('DATABASE', np.nan),
                    'ACCESSION': entry.get('ACCESSION', np.nan)
                })
        else:
            extracted_data.append({
                'GENE': np.nan, '%COVERAGE': np.nan, '%IDENTITY': np.nan,
                'SEQUENCE': np.nan, 'START': np.nan, 'END': np.nan,
                'DATABASE': np.nan, 'ACCESSION': np.nan})
        return extracted_data

    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            sample_name = sample_id
            if data.get("status") == "Success":
                info = data["summary"].get(ariba_type, [])
                data_df[sample_name] = extract_data(info)
            else:
                data_df[sample_name] = extract_data(None)

    rows = []
    for sample, entries in data_df.items():
        for entry in entries:
            rows.append((sample, entry['GENE'], entry['%COVERAGE'], entry['%IDENTITY'],
                         entry['SEQUENCE'], entry['START'], entry['END'],
                         entry['DATABASE'], entry['ACCESSION']))

    df = pd.DataFrame(rows, columns=[
        'Sample', 'GENE', '%COVERAGE', '%IDENTITY', 'SEQUENCE',
        'START', 'END', 'DATABASE', 'ACCESSION'])

    df[['%COVERAGE', '%IDENTITY', 'START', 'END']] = df[[
        '%COVERAGE', '%IDENTITY', 'START', 'END']].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    df.set_index('Sample', inplace=True)
    return df

def parse_assemblatron(list_files):
    """
    Parse genome assembly metrics from Assemblatron output YAMLs.
    Includes GC content, N50, contig counts, and genome sizes at depth.
    """
    d = {}
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split("__")[0]
        with open(file) as f:
            temp = yaml.load(f, Loader=yaml.Loader)
            if temp["status"] == "Success":
                d[sample_id] = [
                    temp["summary"]["GC"], temp["summary"]["N50"],
                    temp["summary"]["bin_contigs_at_1x"],
                    temp["summary"]["bin_contigs_at_10x"],
                    temp["summary"]["bin_coverage_at_1x"],
                    temp["summary"]["bin_length_at_1x"],
                    temp["summary"]["bin_length_at_10x"],
                    temp["summary"]["bin_length_at_25x"],
                    temp["summary"]["snp_filter_10x_10%"]]

    df = pd.DataFrame.from_dict(d, orient='index', columns=[
        "GC %", "N50", "Number of contigs (1x cov.)", "Number of contigs (10x cov.)",
        "Average coverage (1x)", "Genome size at 1x depth",
        "Genome size at 10x depth", "Genome size at 25x depth", "Ambiguous sites"])
    return df








# ------------------------------------
# PARSERS FOR MINION RESULTS
# ------------------------------------


def load_or_na(list_files):
    """
    Load multiple TSV data files or return a DataFrame filled with "NA"
    if files are empty or unreadable. Adds a SampleID column based on the file name.

    Parameters
    ----------
    list_files : list
        List of file paths to TSV files.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame from all readable files.
    """
    all_dfs = []

    for file in list_files:
        
        try:
            sample_id = os.path.splitext(os.path.basename(file))[0].split(".")[0]
            
            df = pd.read_csv(file, sep="\t")
        
            # If file loads but is empty
            if df.empty:
                df = pd.DataFrame(columns=df.columns)
                df.loc[0] = [np.nan] * len(df.columns)
            
            df["Sample"] = sample_id
            


        except Exception as e:
            continue

        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df = combined_df.set_index("Sample")
    
    return combined_df



def parse_nanostat(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    data = {}
    section = None
    
    
    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith("General summary:"):
            section = "general"
            continue
        elif line.startswith("Number, percentage and megabases of reads above quality cutoffs"):
            section = "cutoffs"
            continue
        elif line.startswith("Top 5 highest mean basecall") or line.startswith("Top 5 longest reads"):
            break

        if section == "general" and ':' in line:
            key, value = map(str.strip, line.split(':', 1))
            value = value.replace(',', '')
            try:
                data[key] = float(value)
            except ValueError:
                data[key] = value

        elif section == "cutoffs" and line.startswith('>Q'):
            parts = re.split(r'\t+', line)
            
            if len(parts) >= 2:
                data[parts[0]] = parts[1]

    if not data:
        raise ValueError(f"No valid data found in file: {filename}")

    df = pd.DataFrame([data])
    
    return df


def parse_fallback_summary(filename):
    # Read TSV ignoring the first line
    d = {"number_of_reads":"Number of reads", 
         "number_of_bases": "Total bases",
         "median_read_length": "Median read length" ,
         "mean_read_length": "Mean read length",
         "read_length_stdev" : "STDEV read length",
         "n50":"Read length N50",
         "mean_qual": "Mean read quality",
         "median_qual": "Median read quality",
         "Reads >Q5:" : ">Q5:",
         "Reads >Q7:" : ">Q7:",
         "Reads >Q10:" : ">Q10:",
         "Reads >Q12:" : ">Q12:",
         "Reads >Q15:" : ">Q15:",
         "active_channels": "Active channels"}
    
    df = pd.read_csv(filename, sep='\t', skiprows=1, header=None,index_col=0)
   
    
    df_transformed = df.T.rename(columns=d)
    return (df_transformed)


def parse_nanoplot_summary(list_files):
    all_dfs = []
    index = []
    for file in list_files:
        sample_id = os.path.splitext(os.path.basename(file))[0].split(".")[0].rstrip("_NanoStats")
        index.append(sample_id)
        try:
            df = parse_nanostat(file)
            
           
        except Exception as e:
            print(f"parse_nanostat failed for '{file}': {e}")
            print(f"Attempting fallback parse as TSV: {file}")
            df = parse_fallback_summary(file)
            
            
           
        all_dfs.append(df)
    
    combined_df = pd.concat(all_dfs)
    combined_df.index = index
    combined_df = combined_df.drop(columns=["Active channels",
                                            "longest_read_(with_Q):1",
                                            "longest_read_(with_Q):2",
                                            "longest_read_(with_Q):3",
                                            "longest_read_(with_Q):4",
                                            "longest_read_(with_Q):5",
                                            "highest_Q_read_(with_length):1",
                                            "highest_Q_read_(with_length):2",
                                            "highest_Q_read_(with_length):3",
                                            "highest_Q_read_(with_length):4",
                                            "highest_Q_read_(with_length):5"])
    
    return combined_df


def parse_mlst_nanopore(list_files):
    all_dfs = []

    for file in list_files: 
        sample_id = os.path.splitext(os.path.basename(file))[0].split(".")[0].rstrip("_mlst")
        
        try:
            df = pd.read_csv(file, sep='\t',header=None)
            
            all_dfs.append(df)
            if df.empty:
                df = pd.DataFrame(columns=df.columns)
                df.loc[0] = [np.nan] * len(df.columns)
        
        except:
            continue
        df["Sample"] = sample_id
    combined_df = pd.concat(all_dfs)
    combined_df = combined_df.drop(combined_df.columns[0], axis=1)
    combined_df = combined_df.set_index("Sample")
   
    return combined_df


