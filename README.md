# ChIP-seq CUT&Tag and CUT&RUN Processing Pipeline
The analysis pipeline for protein-DNA interaction profiling assays, including ChIP-seq, CUT&Tag, and CUT&RUN, processes raw FASTQ data through sequential steps such as adapter trimming, quality control, genome alignment, peak calling, fingerprint analysis, and principal component analysis (PCA), and it supports multiple samples and various comparison modes.

# Part I Introduction
## i. Workflow
Here stands an throughout workflow of data analysis.
<img width="2048" height="664" alt="ChIPseq" src="https://github.com/user-attachments/assets/04242ba9-70a6-43da-ace5-1334aeea2d88" />

## ii. Features
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow from raw FASTQ input through trimming, quality control, genome alignment, peak calling, fingerprinting, and PCA can be executed reproducibly on any compatible system, supporting multiple samples and various comparison modes.

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **Download Basement Files**:

      * `run_DNAProteinSeq.sh`
      * `DNAProteinSeq.sif` (The Singularity container)
      * `illumina_adapter.fa`

4.  **Reference Data**: A directory containing bowtie2 index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).

      ```bash
      mkdir basement_data
      cd basement_data

      # Download Genome FASTA
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

      # Unzip the files
      gunzip GRCh38.primary_assembly.genome.fa.gz
      gunzip gencode.v46.primary_assembly.annotation.gtf.gz

      # Remove scafford
      awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa

      # Build index
      mkdir hg38
      singularity exec --cleanenv End-seq.sif bowtie2-build --threads 8 -f GRCh38.primary_assembly.genome.chr.fa ./hg38/bowtie2_index

      # get chromatin size (only used when you select SEACR as the peak calling method)
      samtools faidx GRCh38.primary_assembly.genome.chr.fa
      cut -f1,2 GRCh38.primary_assembly.genome.chr.fa.fai > chromatin.size

      # Remove unnecessary files
      rm GRCh38.primary_assembly.genome.chr.fa
      rm GRCh38.primary_assembly.genome.fa
      ```
      
5.   **Required Input Metadata Files**

      * i. Create a tab-seperated file named `SampleInfor.txt`.
        
        `Sample_prefix`: A unique identifier for the sample (e.g., Control). This name will be used for output subdirectories.
        `R1_path`: The absolute path to the Read 1 FASTQ file, which may be gzipped or uncompressed; soft links are not supported.
        `R2_path`: The absolute path to the Read 2 FASTQ file, which may be gzipped or uncompressed; soft links are not supported.
      
        Note: This pipeline supports both paired-end (PE) and single-end (SE) sequencing data. If you have SE data, simply remove the R2_path column. Ensure that there are no extra spaces or empty lines in the files.
      
        PE Example `SampleInfor.txt`:
        ```
        Sample_prefix	R1_path	R2_path
        Input	/path/to/data/Input_1.fastq.gz	/path/to/data/Input_2.fastq.gz
        H3K27ac	/path/to/data/H3K27ac_1.fastq.gz	/path/to/data/H3K27ac_2.fastq.gz
        H3K4me3	/path/to/data/H3K4me3_1.fastq.gz	/path/to/data/H3K4me3_2.fastq.gz
        ```
      
        SE Example `SampleInfor.txt`:
        ```
        Sample_prefix	path
        Input	/path/to/data/Input.fastq.gz
        H3K27ac	/path/to/data/H3K27ac.fastq.gz
        H3K4me3	/path/to/data/H3K4me3.fastq.gz
        ```

      * ii. (Optional) Create a tab-seperated file named `Comparison.txt`.
      
        `Treatment_prefix`: The identifier of the treatment sample as defined in `SampleInfor.txt`.
        
        `Control_prefix`: The identifier of the corresponding control sample as defined in `SampleInfor.txt`, typically an input sample.
      
        Note: Each line represents one treatment-control comparison. Ensure that there are no extra spaces or empty lines in the files.
      
        Example `Comparison.txt`:
        ```
        Treatment_prefix	Control_prefix
        H3K27ac	Input
        H3K4me3	Input
        ```

6.   **Required File Structure**
      ```bash
      basement_data/
      ├── hg38/
            ├── bowtie2_index.1.bt2
            ├── bowtie2_index.2.bt2
            ├── bowtie2_index.3.bt2
            ├── bowtie2_index.4.bt2
            ├── bowtie2_index.rev.1.bt2
            └── bowtie2_index.rev.2.bt2
      ├── chromatin.size
      ├── Comparison.txt
      ├── DNAProteinSeq.sif
      ├── illumina_adapter.fa
      ├── run_DNAProteinSeq.sh
      └── SampleInfor.txt
      ```

# Part III Running

   * **Example code for ChIP-seq (Histomodification)**

      ```bash
      bash ./basement_data/run_DNAProteinSeq.sh \
                           --sampleInfor ./basement_data/SampleInfor.txt \
                           --comparisonInfor ./basement_data/Comparison.txt \
                           --outputdir ./result_chip_histone \
                           --referencedir ./basement_data/hg38/bowtie2_index \
                           --adapterFa ./basement_data/illumina_adapter.fa \
                           --sif ./basement_data/DNAProteinSeq.sif \
                           --threads 8 \
                           --binSize 10 \
                           --g hs \
                           --X 2000 \                               ## fragments are relatively long; not needed for SE
                           --peakcalling MACS3 \
                           --peaktype broad \                       ## broad peaks
                           --peakPerSample yes \                    ## Peaks can be called per sample without a control and can coexist with `--comparisonInfor`
                           --qval 0.05 \                            ## a relaxed threshold
                           --broad_cutoff 0.1 \                     ## Cutoff for broad peaks
                           --llocal 100000 \                        ## larger local background, used for sample without a control
                           --keepdup all \                          ## Keep all duplicates
                           --nomodel on \
                           --nolambda on \
                           --callsummits off \
                           --extsize_val 147                        ## nucleosome size
      ```
      
   * **Example code for ChIP-seq (Transcription Factors)**

      ```bash
      bash ./basement_data/run_DNAProteinSeq.sh \
                           --sampleInfor ./basement_data/SampleInfor.txt \
                           --comparisonInfor ./basement_data/Comparison.txt \
                           --outputdir ./result_chip_tf \
                           --referencedir ./basement_data/hg38/bowtie2_index \
                           --adapterFa ./basement_data/illumina_adapter.fa \
                           --sif ./basement_data/DNAProteinSeq.sif \
                           --threads 8 \
                           --binSize 5 \
                           --g hs \
                           --X 1000 \                               ## fragments are shorter; not needed for SE
                           --peakcalling MACS3 \
                           --peaktype narrow \                      ## narrow peaks
                           --peakPerSample yes \                    ## Peaks can be called per sample without a control and can coexist with `--comparisonInfor`
                           --pval 0.001 \                           ## a more stringent p-value
                           --llocal 50000 \                         ## smaller local background, used for sample without a control
                           --keepdup 1 \                            ## keep uniquely mapped reads for TFs 
                           --nomodel on \
                           --nolambda on \                          ## used for sample without a control
                           --callsummits on \                       ## precise summit calling
                           --extsize_val 200                        ## Typical fragment size for TFs
      ```

   * **Example code for CUT&Tag and CUT&Run**

      ```bash
      bash ./basement_data/run_DNAProteinSeq.sh \
                           --sampleInfor ./basement_data/SampleInfor.txt \
                           --comparisonInfor ./basement_data/Comparison.txt \
                           --outputdir ./result_cut_tag \
                           --referencedir ./basement_data/hg38/bowtie2_index \
                           --adapterFa ./basement_data/illumina_adapter.fa \
                           --sif ./basement_data/DNAProteinSeq.sif \
                           --threads 8 \
                           --binSize 1 \                            ## Higher resolution
                           --g ./basement_data/chromatin.size \
                           --X 700 \                                ## fragments are shorter; not needed for SE
                           --peakcalling SEACR \                    ## use SEACR for peak calling
                           --peaktype stringent \                   ## High signal-to-noise, use stringent mode
                           --peakPerSample yes \
                           --seacr_threshold 0.01                   ## Strict threshold for top signal regions, used for sample without a control
      ```

   * **Command Parameters**

      - `--sampleInfor`:        (required) A tab-separated file containing the sample prefix and the path(s) to the FASTQ file(s)
      - `--comparisonInfor`:    (optinal) A tab-separated file containing the prefixes of treatment and control samples for comparison
      - `--outputdir`:          (required) Path to the directory where the output will be stored
      - `--referencedir`:       (required) Path to the directory where bowtie reference build with prefix
      - `--adapterFa`:          (required) Path to the adapter fasta
      - `--sif`:                (required) Path to the singularity environment file
      - `--threads`:            (optional) Number of threads to use (default: 8)
      - `--binSize`:            (optional) Number of binsize to use (default: 10)
      - `--peakcalling`:        (optional) Number of binsize to use (default: 10)
      - `--g`:                  (optional) If `--peakcalling` is set to `MACS3`, provide the species code accepted by MACS3, for example: `hs` (human), `mm` (mouse), `ce` (C. elegans), `dm` (Drosophila melanogaster), etc; If `--peakcalling` is set to `SEACR`, provide the path to the corresponding `chromatin.size` file for the genome of interest.
      - `--X`:                  (optional) Maximum fragment length for paired-end mapping in Bowtie2. Typically 700 for Cut&Tag and Cut&Run, 1000 for ChIP-seq (TF), 2000 for ChIP-seq (histone marks)
      - `--peakcalling`:        (optional) Peak calling method, either `MACS3` or `SEACR` (default: `MACS3`)
      - `--peaktype`:           (optional) When using MACS3, set `broad` or `narrow`. When using SEACR, set `relaxed` or `stringent` (default: `narrow`)
      - `--peakPerSample`:      (optional) Call peaks for each sample individually, Sst to `yes` if --comparisonInfor is not provided, or optionally yes even when --comparisonInfor is provided, `yes` or `no` (default: `yes`)
      - `--pval`:               (optional) For MACS3, P-value cutoff for peak calling to determine significance of narrow/strong peaks. If specified, MACS3 uses p-value instead of q-value
      - `--qval`:               (optional) For MACS3, Q-value (FDR) cutoff for narrow/strong peaks, controlling false discovery rate (default: 0.01)
      - `--broad_cutoff`:       (optional) For MACS3, P-value or q-value cutoff for broad/weak peaks, effective only when `--peaktype` is `broad` (default: 0.1)
      - `--llocal`:             (optional) For MACS3, `--llocal` value for samples without input control (default: 100000)
      - `--keepdup`:            (optional) For MACS3, `--keep-dup` setting (default: `all`)
      - `--nomodel`:            (optional) For MACS3, Disable model building, set `on` or `off` (default: `on`)
      - `--nolambda`:           (optional) For MACS3, Disable dynamic lambda, set `on` or `off` (default: `on`)
      - `--callsummits`:        (optional) For MACS3, Enable calling of peak summits within enriched regions, set `on` or `off` (default: `off`)
      - `--extsize_val`:        (optional) For MACS3, Set fragment/extension size for MACS3 in base pairs, effective only when `--nomodel` is `on`
      - `--shift_val`:          (optional) For MACS3, Set shift for read 5' ends in base pairs for MACS3, effective only when `--nomodel` is `on`; positive moves 5'->3', negative moves 3'->5'
      - `--seacr_threshold`:    (optional) For SEACR, Percentile-based cutoff for SEACR, e.g., 0.01 = top 1% of signal regions. Recommended 0.01 for sharp marks like H3K4me3 and 0.05 for broad marks like H3K27me3 (default: 0.01)




