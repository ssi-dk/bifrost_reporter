#!/usr/bin/env python

# Import necessary libraries
from bifrost_reporter import data_processing
import logging
import pandas as pd
import os



def retrieve_samples(sample_sheet_path):
    """
    Reads an Excel sample sheet and builds full paths to each sample directory.

    Parameters:
    ----------
    sample_sheet_path : str
        Path to the Excel file containing sample IDs

    Returns:
    -------
    pd.Series
        Series of full paths to sample directories
    """

    # Load file based on extension
    if sample_sheet_path.endswith(('.xls', '.xlsx')):
        df = pd.read_excel(sample_sheet_path)
    elif sample_sheet_path.endswith('.tsv'):
        df = pd.read_csv(sample_sheet_path, sep='\t')
    else:
        raise ValueError("Unsupported file type. Use .tsv, .xls, or .xlsx.")

    # Create full paths
    base_path = os.path.dirname(os.path.abspath(sample_sheet_path))
    df["SampleID"] = base_path + '/' + df["SampleID"].astype(str)

    return df["SampleID"]



def check_samples(folder_paths, 
                  mode):
    """
    Checks each sample directory for the expected output files from Bifrost.

    Parameters:
    ----------
    folder_paths : list of str
        List of sample directory paths
    mode : str
        Sequencing mode: 'illumina' or 'nanopore'

    Returns:
    -------
    dict
        A dictionary showing whether each sample directory and required files exist
    """

    illumina_files = [
        "__amrfinderplus_fbi.yaml", 
        "__ariba_mlst.yaml", 
        "__ariba_plasmidfinder.yaml",
        "__ariba_resfinder.yaml", 
        "__ariba_virulencefinder.yaml", 
        "__assemblatron.yaml",
        "__kma_pointmutations.yaml", 
        "__min_read_check.yaml", 
        "__reslab_stamper.yaml",
        "__sp_cdiff_fbi.yaml", 
        "__sp_ecoli_fbi.yaml", 
        "__sp_salm_fbi.yaml", 
        "__ssi_stamper.yaml", 
        "__whats_my_species.yaml"
    ]

    nanopore_files = ["_NanoStats.txt",
                      ".unicycler_nponly_resfinder.tsv",
                      ".unicycler_nponly_plasmidfinder.tsv"
    ]

    expected_files = illumina_files if mode == 'illumina' else nanopore_files

    status = {}
    for folder in folder_paths:
        file_status = {}
        sample_name = os.path.basename(folder)

        if os.path.isdir(folder):  # Check if folder exists
            for result in expected_files:
                file_name = sample_name + result # careful with nanostat
                file_path = os.path.join(folder, file_name)
                file_status[file_name] = os.path.isfile(file_path)
        else:
            logging.error(f"The folder : {folder} could not be found. Please check again")

        status[folder] = {
            'exists': os.path.isdir(folder),
            'files': file_status
        }

    return status



def data_collection_from_dict(sample_dict, 
                              mode):
    """
    Gathers and processes data from YAML result files found in each sample directory.

    Parameters:
    ----------
    sample_dict : dict
        Output from check_samples() containing file presence per sample

    Returns:
    -------
    tuple of DataFrames
        Parsed results from different Bifrost analyses
    """
    analysis_files = {}

    for sample_path, sample_info in sample_dict.items():
        if not sample_info["exists"]:
            continue  # Skip if directory doesn't exist

        for file_name, is_present in sample_info["files"].items():
            if is_present:
                try:
                    if mode == "illumina":
                        analysis_name = file_name.strip(".yaml").split("__")[1]
                        
                    elif mode == "nanopore" and file_name.endswith(".txt"):
                        analysis_name = file_name.strip(".txt").split("_")[-1]
                        
                    elif mode == "nanopore" and file_name.endswith(".tsv"):
                        analysis_name = file_name.strip(".tsv").split("_")[-1]

                    analysis_files.setdefault(analysis_name, []).append(os.path.join(sample_path, file_name))
                except:
                    continue  # Ignore parsing errors
    
    
    if mode == "illumina":
        try:
            # Parse MLST analysis
            print("[INFO] Parsing MLST...")
            mlst_df = data_processing.parse_mlst(analysis_files.get("ariba_mlst", []))
        except Exception as e:
            print(f"[ERROR] Failed to parse MLST: {e}")
            mlst_df = pd.DataFrame()  # Return empty dataframe if failure

        try:
            # Parse PlasmidFinder
            print("[INFO] Parsing PlasmidFinder...")
            plasmid_finder_df = data_processing.parse_finder_tools(analysis_files.get("ariba_plasmidfinder", []), "ariba_plasmidfinder")
            plasmid_finder_df = plasmid_finder_df[(plasmid_finder_df["%COVERAGE"] >= 80) & (plasmid_finder_df["%IDENTITY"] >= 80)]
        except Exception as e:
            print(f"[ERROR] Failed to parse PlasmidFinder: {e}")
            plasmid_finder_df = pd.DataFrame()

        try:
            # Parse ResFinder
            print("[INFO] Parsing ResFinder...")
            resfinder_df = data_processing.parse_finder_tools(analysis_files.get("ariba_resfinder", []), "ariba_resfinder")
            resfinder_df = resfinder_df[(resfinder_df["%COVERAGE"] >= 60) & (resfinder_df["%IDENTITY"] >= 90)]
        except Exception as e:
            print(f"[ERROR] Failed to parse ResFinder: {e}")
            resfinder_df = pd.DataFrame()

        try:
            # Parse VirulenceFinder
            print("[INFO] Parsing VirulenceFinder...")
            virulence_df = data_processing.parse_finder_tools(analysis_files.get("ariba_virulencefinder", []), "ariba_virulencefinder")
            virulence_df = virulence_df[(virulence_df["%COVERAGE"] >= 60) & (virulence_df["%IDENTITY"] >= 90)]
        except Exception as e:
            print(f"[ERROR] Failed to parse VirulenceFinder: {e}")
            virulence_df = pd.DataFrame()

        try:
            # Parse Assemblatron
            print("[INFO] Parsing Assemblatron...")
            assemblatron_df = data_processing.parse_assemblatron(analysis_files.get("assemblatron", []))
        except Exception as e:
            print(f"[ERROR] Failed to parse Assemblatron: {e}")
            assemblatron_df = pd.DataFrame()

        try:
            # Parse KMA Point Mutations
            print("[INFO] Parsing KMA Point Mutations...")
            kma_df = data_processing.parse_kmapointmutations(analysis_files.get("kma_pointmutations", []))
        except Exception as e:
            print(f"[ERROR] Failed to parse KMA Point Mutations: {e}")
            kma_df = pd.DataFrame()

        try:
            # Parse AMRFinderPlus
            print("[INFO] Parsing AMRFinderPlus...")
            amr_df = data_processing.parse_amrfinder(analysis_files.get("amrfinderplus_fbi", []))
        except Exception as e:
            print(f"[ERROR] Failed to parse AMRFinderPlus: {e}")
            amr_df = pd.DataFrame()

        try:
            # Check SSI Stamper
            print("[INFO] Checking SSI Stamper...")
            ssi_stamper_df = data_processing.check_stampers(analysis_files.get("ssi_stamper", []))
        except Exception as e:
            print(f"[ERROR] Failed to check SSI Stamper: {e}")
            ssi_stamper_df = pd.DataFrame()

        try:
            # Check ResLab Stamper
            print("[INFO] Checking ResLab Stamper...")
            reslab_stamper_df = data_processing.check_stampers(analysis_files.get("reslab_stamper", []))
        except Exception as e:
            print(f"[ERROR] Failed to check ResLab Stamper: {e}")
            reslab_stamper_df = pd.DataFrame()

        try:
            # Check E. coli Stampers
            print("[INFO] Checking E. coli Stampers...")
            ecoli = data_processing.check_stampers(analysis_files.get("sp_ecoli_fbi", []))
        except Exception as e:
            print(f"[ERROR] Failed to check E. coli Stampers: {e}")
            ecoli = pd.DataFrame()

        return (
            mlst_df,
            plasmid_finder_df,
            resfinder_df,
            virulence_df,
            assemblatron_df,
            kma_df,
            amr_df,
            ssi_stamper_df,
            reslab_stamper_df,
            ecoli
        )

    elif mode == "nanopore":
        plasmid_finder_df = pd.DataFrame()
        resfinder_df = pd.DataFrame()
        nanostat = pd.DataFrame()

        try:
            # Parse PlasmidFinder for Nanopore
            print("[INFO] Parsing PlasmidFinder...")
            plasmid_finder_df = data_processing.load_or_na(analysis_files.get("plasmidfinder", []))
            plasmid_finder_df = plasmid_finder_df[(plasmid_finder_df["%COVERAGE"] >= 80) & (plasmid_finder_df["%IDENTITY"] >= 80)]
        except Exception as e:
            print(f"[WARNING] Failed to parse plasmidfinder: {e}")

        try:
            # Parse ResFinder for Nanopore
            print("[INFO] Parsing ResFinder...")
            resfinder_df = data_processing.load_or_na(analysis_files.get("resfinder", []))
            resfinder_df = resfinder_df[(resfinder_df["%COVERAGE"] >= 60) & (resfinder_df["%IDENTITY"] >= 90)]
        except Exception as e:
            print(f"[WARNING] Failed to parse resfinder: {e}")

        try:
            # Parse NanoPlot Summary for Nanopore
            print("[INFO] Parsing NanoPlot Summary...")
            nanostat = data_processing.parse_nanoplot_summary(analysis_files.get("NanoStats", []))
        except Exception as e:
            print(f"[WARNING] Failed to parse NanoStats: {e}")

        return (
            plasmid_finder_df,
            resfinder_df,
            nanostat
        )
        # missing other stuff


