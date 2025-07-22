import pandas as pd 
import os 
import subprocess
list_of_files = os.listdir('/data/passala/Collaborator_Data/Hagai_zach_circadian/RNA_data')
def decode_file_name(file_name):
    parts = file_name.split('_')
    species = parts[4]
    replicate = parts[5]
    time_point = parts[6]
    if int(time_point) <12:
        time_point = time_point + 'am'
    elif int(time_point) >= 12:
        time_point = time_point + 'pm'
    else: 
        raise Exception("No time point found in file name")

    summary_dictionary = {'species': species, 'replicate': replicate, 'time_point': time_point}
    return summary_dictionary
unique_file_prefixes = pd.Series(list_of_files).str.split('L1').str[0]
list_of_unique_file_prefixes = unique_file_prefixes.unique().tolist()
genome_locations = pd.read_csv('/data/passala/Generated_Tables/Temp_junk/list_of_genome_fastas.csv', sep='\t')
file_location = '/data/passala/Collaborator_Data/Hagai_zach_circadian/RNA_data/'
processed_data_location = '/data/passala/Collaborator_Data/Hagai_zach_circadian/processed_RNA_data/'
for file_prefix in list_of_unique_file_prefixes: 
    r1_file = file_location + file_prefix + 'L1_R1.fastq' 
    r2_file = file_location + file_prefix + 'L1_R2.fastq'
    if file_prefix == 'HT1004_Solanum_SP5G_A07_BLANK_10481877_232FCWLT3_': 
        subprocess.run(f'/data/passala/Module_paper_data/Populus_vs_glycine_drought/STAR --runThreadN 30 --genomeDir /data/passala/Genomes/Solanum_lycopersicum_4_0 --readFilesIn {r1_file} {r2_file} --outSAMtype None --outFileNamePrefix {processed_data_location+file_prefix} --quantMode GeneCounts', shell = True, check = True)
        continue
    decoded_file_characteristics = decode_file_name(file_prefix)
    species = decoded_file_characteristics['species']
    genome_directory = genome_locations.loc[genome_locations['Abbreviation'] == species]['Genome Folder'].item()
    subprocess.run( f'/data/passala/Module_paper_data/Populus_vs_glycine_drought/STAR --runThreadN 30 --genomeDir {genome_directory} --readFilesIn {r1_file} {r2_file} --outSAMtype None --outFileNamePrefix {processed_data_location+file_prefix} --quantMode GeneCounts' , shell = True, check = True)
