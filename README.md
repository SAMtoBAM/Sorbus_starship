# Sorbus_starship
Walkthrough detailing how the SORBUS cluster in _Penicillium roqueforti_ was confirmed as a _Starship_


[Punt et al. 2022](https://doi.org/10.1371/journal.pgen.1010086) detailed a cluster of genes that conferred Sorbic Acid resistance in a subset of _Pencillium roqueforti_ strains isolated from bread <br/>
The size of this region (~180kb) and presence in several independent lineages made it a good candidate to be a _Starship_ element <br/>
However most analyses from Punt et al. used a strain (DTO006G7) in which only an internal portion of the full region was assembled; additionally the long-read assembly used to confirm deletions was not published <br/>
Therefore this walkthrough details how I extracted the long-reads, re-assembled the genome, found the SORBUS cluster and identified the _Starship_ element I am calling Sorbus <br/>
Additionally, I have found that the this _Starship_ is in two insertion sites therefore indicating this element is actively transposing and horizontally transferring between _P. roqueforti_ strains.


    ##Use the fusemblr conda environment plus a few additions in order to download raw reads and both assemble and evaluate the resulting assembly
    ##fusemblr is an pipeline wrapper (https://github.com/SAMtoBAM/fusemblr) used to assemble fungal genomes and evaluate them automatically using PAQman (https://github.com/SAMtoBAM/PAQman)  
    conda create -n fusemblr samtobam::fusemblr
    conda activate fusemblr
    conda install ncbi-datasets-cli sra-tools

    ##set the number of threads to use
    threads="16"

    ##ok now we are ready to proceed

    ##first step is get the raw reads from sra
    ##the raw long-reads have the accession SRR17178875
    ##this long-read dataset comes from the KU knock-out strain of DTO013F2 before deleting genes within the SORBUS cluster    
    ##and the short reads for the same strain and study
    strain="DTO013F2"
    LR="SRR17178875"
    SR="SRR17211221"

    ##now use fasterq-dump to extract our reads (should be named is in the accession set variables)
    fasterq-dump $LR
    ##this time split into the paired ends for the illumina data
    fasterq-dump --split-3 $SR

    ##zip up the reads to save space
    gzip *.fastq

    ##now we can jump straight to assembling this strain with fusemblr with an estimated genome size of 27Mb
    fusemblr.sh -n ${LR}.fastq.gz -1 ${SR}_1.fastq.gz -2 ${SR}_2.fastq.gz -g 27000000 -p ${strain} -o ${strain}_fusemblr -t ${threads}

    


    
