# Sorbus_starship
## Walkthrough detailing how the SORBUS cluster in _Penicillium roqueforti_ was confirmed as a _Starship_ I am simply calling sorbus


[Punt et al. 2022](https://doi.org/10.1371/journal.pgen.1010086) detailed a cluster of genes that conferred Sorbic Acid resistance in a subset of _Pencillium roqueforti_ strains isolated from bread <br/>
The size of this region (~180kb) and presence in several independent lineages made it a good candidate to be a _Starship_ element <br/>
However most analyses from Punt et al. used a strain (DTO006G7) in which only an internal portion of the full region was assembled; additionally the long-read assembly used to confirm deletions was not published <br/>
Therefore this walkthrough details how I extracted the long-reads, re-assembled the genome, found the SORBUS cluster and identified the _Starship_ element I am calling Sorbus <br/>
Additionally, I have found that the this _Starship_ is in two insertion sites therefore indicating this element is actively transposing and horizontally transferring between _P. roqueforti_ strains.


### Step 1: Get a long-read assembly of the SORBUS containing strain DTO013F2

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

    ##use the PAQman output to evaluate the best assembly
    ##here we can see the best one is the Flye assembly (no redundant contigs to remove; 4 complete contiguous chromosomes plus 1 contig for MT)
    ##but the Flye assembly is essentially the same with a few additional rDNA contigs 
    ##No telomeres at all detected by PAQman and nothing obvious by visual inspection either; same for the public assembly too...

    ##BIG SIDE NOTE: Does roqueforti have a different telomeric repeat? something like:
    ## '[AG]AGTT[AG]AT[TAC][GAC]TC[TG]' with the degenerecy although generally a little distant from the actual end...need to work on this
    ## using this repeat then the Flye assembly is T2T excluding the rDNA capped chromosome end (and better than the hifiasm version)
    
    assembler="flye"
    
    ##copy into into the current directory for all genomes
    mkdir genomes
    cat ${strain}_fusemblr/${strain}.${assembler}.final.fa | sed "s/>/>${strain}_/g" | sed 's/contig_/contig/g' > ./genomes/${strain}.fa

### Step 2: Use other publicaly available long-read assembly to identify _Starships_ in the new long-read assembly
#### Step 2a: Get public long-read assembly LCP06133
    ##we will use the a new environment for this next step
    ##this is because we want to both download the publicly available long-read assembly for _P. roqueforti_ LCP06133 AND analyse _Starships_
    ##this public long-read assembly has the accession GCA_030518555.1 used below
    reference="GCA_030518555.1"

    ##download the assembly then rename it etc
    datasets download genome accession ${reference}
    unzip ncbi_dataset.zip
    rm ncbi_dataset.zip
    ##rename them as just the ncbi GCA assession and rename the contig headers with the genome name and get rid of softmasking
    ls ncbi_dataset/data/ | grep -v json | while read genome
    do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna | sed "s/>/>${genome2}_/g" | awk -F " " '{print $1}' | awk '{if($0 ~ ">") {print} else {print toupper($0)}}' > genomes/$genome2.fa
    done
    ##clean up
    rm -r ncbi_dataset/ md5sum.txt README.md

    ##now the LCP06133 assembly is in the 'genomes' folder with the new DTO013F2 assembly

#### Step 2b: Identify Starships with just the two assemblies
    ##now we will use a new environment in order to run starfish (we only need starfish as this _Starship_ should be easy to identify)
    ##but we will use the Stargraph wrapper to run all the steps

    conda create -n stargraph samtobam::stargraph
    conda activate stargraph

    ls genomes/*.fa > assemblies_list.txt
    starfish_wrapper.sh -a assemblies_list.txt -t ${threads} -o starfish_wrapper_output

    conda deactivate 

### Step 3: Extract SORBUS cluster from short read assembly and identify region in long-read assembly

    ##we will use the same environment as above 'fusemblr' to download the publicly available assembly for GCA_023141355.1
    ##in this strain (DTO006G7) the cluster was identified in scaffold43 which corresponds to JAKIKL010000044.1
    ##this is a 174kb contig

    sorbusref="GCA_023141355.1"
    sorbuscluster="JAKIKL010000044.1"

    mkdir sorbus_ref
    
    ##download the assembly then rename it etc
    datasets download genome accession ${sorbusref}
    unzip ncbi_dataset.zip
    rm ncbi_dataset.zip
    ##rename them as just the ncbi GCA assession and rename the contig headers with the genome name 
    ls ncbi_dataset/data/ | grep -v json | while read genome
    do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna  > sorbus_ref/$genome2.fa
    done
    ##clean up
    rm -r ncbi_dataset/ md5sum.txt README.md  

    ##extract the sorbus cluster then look for the corresponding region in DTO013F2
    ##use seqkit to grab the contig based on the ID
    seqkit grep -p "${sorbuscluster}" sorbus_ref/*.fa > sorbus_cluster.fa
    ##use blastn to find the region and filter out small matches
    blastn -query sorbus_cluster.fa -subject genomes/DTO013F2.fa -outfmt 6 | awk '{if($4 > 10000) print}'

    JAKIKL010000044.1	DTO013F2_contig1	99.900	77857	13	4	3475	81270	10647787	10569935	0.0	1.433e+05
    JAKIKL010000044.1	DTO013F2_contig1	99.991	64392	5	1	81461	145852	10569933	10505543	0.0	1.189e+05
    JAKIKL010000044.1	DTO013F2_contig1	99.974	34974	7	1	145895	180868	10505552	10470581	0.0	64534

    ##corresponds to a starship found DTO013F2_s00007 (~400kb in total size)
    DTO013F2_s00007	phoenixFam	navis0004-hap0003	DTO013F2_contig1	DTO013F2_tyr4	10297813	10689168	391356	-	insert	GCA030518555_site002	GCA030518555_CP130315.1	46605	46612	ttcttaca

    ##got the captain 
    DTO013F2_contig1	MetaEuk	mRNA	10686206	10688619	1054	-	.	Target_ID=penrub2_KAF3021169.1;Name=DTO013F2_tyr4

    ##with a MYB/SANT gene at the opposite edge
    DTO013F2_contig1	MetaEuk	mRNA	10298717	10299988	782	+	.	Target_ID=IBT25940-HTR6_g15471.t1;Name=DTO013F2_myb6

    
