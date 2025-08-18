# Declare the essential variables

conf_file_suffix: str = "generated_conformers"     # Suffix of input coordinate file (containing coordinates of all conformers)
init_prefix: str = "initial"     # Prefix for the files associated with initially generated conformers
opt_prefix: str = "optimized"     # Prefix for the files associated with geometry optimized conformers
# conf_extension: str = ".xyz"     # Extension of the files with conformer coordinates
conf_extension: str = ".sdf"     # Extension of the files with conformer coordinates

conf_prefix: str = "conformer"     # Prefix given to each conformer name

init_conf_xyz: str = init_prefix + "_" + conf_file_suffix + conf_extension      # XYZ file to store coordinates of generated conformers
opt_conf_xyz: str = opt_prefix + "_" + conf_file_suffix + conf_extension      # XYZ file to store coordinates of geometry optimized conformers

opt_conf_energy_csv: str = "opt_conf_energies.csv"     # CSV file to store total energies of geometry optimized conformers

opt_conf_SMILES_file: str = opt_prefix + "_" + conf_file_suffix + ".smi"     # SMI file to store SMILES of geometry optimized conformers

similarity_output_csv: str = "tanimoto_similarity.csv"        # CSV file with Tanimoto similarity values between the SMILES of geometry optimized conformers and input SMILES
feasible_geometries_csv: str = "feasible_geom_energies.csv"     # CSV file to log conformers with feasible geometries
infeasible_geometries_csv: str = "infeasible_geometries.csv"     # CSV file to log conformers with infeasible geometries

feasible_geometries_xyz: str = "feasible_geometries.xyz"     # XYZ file to store coordinates of conformers with feasible geometries
infeasible_geometries_xyz: str = "infeasible_geometries.xyz"     # XYZ file to store coordinates of conformers with infeasible geometries

pairwise_RMSDs_dat: str = "rmsd_matrix-mx_flags.dat"     # DAT file to store pairwise RMSDs of conformers with feasible geometries
pairwise_RMSDs_csv: str = "pairwise_RMSDs.csv"     # CSV file to store pairwise RMSDs of conformers with feasible geometries

cluster_reps_csv: str = "cluster_rep_conformers.csv"     # CSV file to store details of cluster representative conformers
cluster_reps_xyz: str = "cluster_rep_conformers.xyz"     # XYZ file to store coordinates of cluster representative conformers
cluster_rep_prefix: str = "rep_of_cluster_"     # Prefix to the xyz files of individual cluster representatives
cluster_reps_dir: str = "cluster_rep_conformers"     # Directory in which the xyz files of individual cluster representatives must be stored

clusters_RMSD_stats_csv: str = "cluster_statistics-RMSDs.csv"     # CSV file to store statistics of pairwise RMSDs within each cluster
clusters_energy_stats_csv: str = "cluster_statistics-energies.csv"     # CSV file to store statistics of relative energies within each cluster

opt_cluster_reps_csv: str = "opt_cluster_rep_conformers.csv"     # CSV file to store MMFF and xTB energies of cluster representative conformers


sdf_conf_extension: str = ".sdf"     
cluster_reps_sdf: str = "cluster_rep_conformers.sdf"     # XYZ file to store coordinates of cluster representative conformers
init_conf_sdf: str = init_prefix+"_"+conf_file_suffix+sdf_conf_extension  
opt_conf_sdf: str = opt_prefix+"_"+conf_file_suffix+sdf_conf_extension  

feasible_geometries_sdf: str = "feasible_geometries.sdf"     # CSV file to store coordinates of conformers with feasible geometries
infeasible_geometries_sdf: str = "infeasible_geometries.sdf"     # XYZ file to store coordinates of conformers with infeasible geometries