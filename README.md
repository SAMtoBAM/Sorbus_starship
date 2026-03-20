# Sorbus_starship
Walkthrough detailing how the SORBUS cluster in _Penicillium roqueforti_ was confirmed as a _Starship_


[Punt et al. 2022](https://doi.org/10.1371/journal.pgen.1010086) detailed a cluster of genes that conferred Sorbic Acid resistance in a subset of _Pencillium roqueforti_ strains isolated from bread <br/>
This size of this region (~180kb) and presence in several independent lineages made it a good candidate to be a _Starship_ element <br/>
However most analyses used a strain (DTO006G7) only assembled an internal region of the full region; additionally the long-read assembly used to confirm deletions was not published <br/>
Therefore this walkthrough details how I extracted the long-reads, re-assembled the genome, found the SORBUS cluster and identified the _Starship_ element I am calling Sorbus <br/>
Additionally, I have found that the this _Starship_ is in two insertion sites therefore indicating this element is actively transposing and horizontally transferring between _P. roqueforti_ strains.


    ##Use the fusemblr conda environment in order to download raw reads and both assemble and evaluate the resulting assembly
    ##fusemblr is an pipeline wrapper (https://github.com/SAMtoBAM/fusemblr) used to assemble fungal genomes and evaluate them automatically using PAQman (https://github.com/SAMtoBAM/PAQman)  
    conda create -n fusemblr samtobam::fusemblr
    conda activate fusemblr

    
