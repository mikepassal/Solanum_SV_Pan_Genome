#!/bin/bash
# keep the above line!


# Let's first re-generate the genome using an up-to-date STAR
# ./STAR  --runThreadN 24 --runMode genomeGenerate --genomeDir /home/suresh/salmon/sally --genomeFastaFiles /home/suresh/salmon/GCF_000233375.1_ICSASG_v2_genomic.fna --sjdbGTFfile /home/suresh/salmon/GCF_000233375.1_ICSASG_v2_genomic.gtf --sjdbOverhang 100


dirname="/data/passala/Collaborator_Data/Marchantia_data/Full_marchan_network"
cd $dirname
# comment genDir - /data/shah/coexp_rna/STARindex/yeast/R64
genDir="/data/passala/Genomes/Marchantia_polymorpha"

# list of SRA accessions that need to be realigned i.e. SRP##### 
input="/data/passala/Collaborator_Data/Marchantia_data/Full_marchan_network/marchantia_accs.csv"

nCore=1       # attempted to parallel download fastq files, but ends up running into connection issues on the ftp side if you try to download too many things.
nCoreSTAR=30  

getFastq() {

    local Run=$1  # first position

    link="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    first6="$(cut -b 1-6 <<<"$Run")"
    last="${Run: -1}"
    last2="${Run: -2}"
    L=$(echo -n $Run | wc -m)
    
    if [ $L = 11 ]  ;then # 8 digits
        link="${link}${first6}/0${last2}/${Run}/*"
        wget -nv -nc  "$link"
        
    elif [ $L = 10 ]  ;then # 7 digits
        link="${link}${first6}/00${last}/${Run}/*"
        wget -nv -nc  "$link"

    elif [ $L = 9 ]; then # 6 digits
        link="${link}${first6}/${Run}/*"
        wget -nv -nc  "$link"
    fi 
    
}

# ./STAR-2.6.0c/bin/Linux_x86_64/STAR

while IFS=',' read -r SRA
do
    it=0
    SRAfile="${SRA}_accList_1.txt"  # runs within the SRA 
    echo $SRAfile   
    while IFS= read -r Run
    do
        ((it=it%nCore)) ; ((it++==0)) && wait
        
        (getFastq $Run 
        
        filein=$Run".fastq.gz"
        filesin_paired=$Run"_1.fastq.gz"
        filesin_subread=$Run"_subreads.fastq.gz"
                    
        if test -f $filesin_paired; then
            filein="${filesin_paired} ${Run}_2.fastq.gz"
        elif test -f $filein; then
            filein=$filein        
        fi

        ./STAR --runThreadN $nCoreSTAR --genomeDir $genDir --readFilesIn $filein --readFilesCommand zcat --outSAMtype None --outFileNamePrefix $Run --quantMode GeneCounts
        if test -f "${Run}ReadsPerGene.out.tab" ; then
            rm ${Run}*.fastq.gz
        fi
        # runSTAR $Run $nCoreSTAR $genDir $filein 
        ) &
        

        # foo $Run $SRA &
    done < $SRAfile
    #wait 

    mkdir $SRA
    mv *RR* $SRA

done < $input


