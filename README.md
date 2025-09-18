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

		| Sample_prefix | R1_path                                | R2_path                                |
        |--------------|-----------------------------------------|-----------------------------------------|
        | Input        | /path/to/data/Input_1.fastq.gz          | /path/to/data/Input_2.fastq.gz          |
        | H3K27ac      | /path/to/data/H3K27ac_1.fastq.gz        | /path/to/data/H3K27ac_2.fastq.gz        |
        | H3K4me3      | /path/to/data/H3K4me3_1.fastq.gz        | /path/to/data/H3K4me3_2.fastq.gz        |

        SE Example `SampleInfor.txt`:

        | Sample_prefix | path                                   |
        |--------------|----------------------------------------|
        | Input        | /path/to/data/Input.fastq.gz           |
        | H3K27ac      | /path/to/data/H3K27ac.fastq.gz         |
        | H3K4me3      | /path/to/data/H3K4me3.fastq.gz         |

      * ii. (Optional) Create a tab-seperated file named `Comparison.txt`.
      
        `Treatment_prefix`: The identifier of the treatment sample as defined in `SampleInfor.txt`.
        
        `Control_prefix`: The identifier of the corresponding control sample as defined in `SampleInfor.txt`, typically an input sample.
      
        Note: Each line represents one treatment-control comparison. Ensure that there are no extra spaces or empty lines in the files.
      
        Example `Comparison.txt`:

		| Treatment_prefix | Control_prefix |
		|------------------|----------------|
		| H3K27ac          | Input         |
		| H3K4me3          | Input         |

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

# Part IV Output

   * **Output Structure**
      ```bash
      outputdir/
      ├── bam/
            ├── Input.bowtie.stats
            ├── Input.markdup.log
            ├── Input.DeDup.bam
            ├── Input.DeDup.bam.bai
            ├── Input.flagstat.txt
            ├── H3K27ac.bowtie.stats
            ├── H3K27ac.markdup.log
            ├── H3K27ac.DeDup.bam
            ├── H3K27ac.DeDup.bam.bai
            ├── H3K27ac.flagstat.txt
            ├── H3K4me3.bowtie.stats
            ├── H3K4me3.markdup.log
            ├── H3K4me3.DeDup.bam
            ├── H3K4me3.DeDup.bam.bai
            └── H3K4me3.flagstat.txt
      ├── bw/
            ├── Input.DeDup.bw
            ├── H3K27ac.DeDup.bw
            └── H3K4me3.DeDup.bw
      ├── figure/
            ├── Input.peak.pdf
            ├── H3K27ac.peak.pdf
            ├── H3K4me3.peak.pdf
            ├── fingerprints.pdf
            ├── BW_compare_PCA.pdf
            ├── BW_compare_cor.pdf
            ├── H3K27ac.vs.Input.peak.pdf
            └── H3K4me3.vs.Input.peak.pdf
      ├── peak.per.sample/
            ├── Input.macs3.stats
            ├── Input_peaks.xls
            ├── Input_peaks.broadPeak
            ├── Input_peaks.gappedPeak
            ├── H3K27ac.macs3.stats
            ├── H3K27ac_peaks.xls
            ├── H3K27ac_peaks.broadPeak
            ├── H3K27ac_peaks.gappedPeak
            ├── H3K4me3.macs3.stats
            ├── H3K4me3_peaks.xls
            ├── H3K4me3_peaks.broadPeak
            └── H3K4me3_peaks.gappedPeak
      ├── peak.compare/
            ├── H3K27ac.vs.Input.macs3.stats
            ├── H3K27ac.vs.Input_peaks.xls
            ├── H3K27ac.vs.Input_peaks.broadPeak
            ├── H3K27ac.vs.Input_peaks.gappedPeak
            ├── H3K4me3.vs.Input.macs3.stats
            ├── H3K4me3.vs.Input_peaks.xls
            ├── H3K4me3.vs.Input_peaks.broadPeak
            └── H3K4me3.vs.Input_peaks.gappedPeak
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      ```

   * **Output Interpretation**

      - **`*.bowtie.stats`**

        - **Content**: Contains Bowtie alignment summary statistics, including the total number of reads processed, reads aligned, reads discarded, and uniquely mapped reads. It provides an overview of mapping quality and efficiency for each FASTQ file.
        - **Application**: Used to assess alignment quality and sequencing library performance. These statistics help in troubleshooting mapping issues, evaluating experiment success, and can be parsed by downstream tools like MultiQC for visualization and comparison across samples.
          
          <img width="410" height="253" alt="图片" src="https://github.com/user-attachments/assets/958fa19a-a967-4424-8c91-2df3d6a18511" />

      - **`*.DeDup.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.bw`**

        - **Content**: A BigWig file that represents the End-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bigWig.html
        - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.flagstat.txt`**

        - **Content**: Contains alignment statistics generated by samtools flagstat after removing duplicate reads. It reports the total number of reads, mapped reads, properly paired reads, singletons, and the number of duplicate reads removed, providing a summary of the final, deduplicated BAM file.
        - **Application**: Used to evaluate the quality of the deduplicated alignment, check library complexity, and ensure that downstream analyses (e.g., peak calling, coverage calculation) are based on high-quality, non-redundant reads.

          <img width="405" height="329" alt="图片" src="https://github.com/user-attachments/assets/27893e2a-f8ff-43c3-a234-b297f50f0fad" />

      - **`*.markdup.log`**

        - **Content**: Log file generated by `samtools markdup`, summarizing read duplication. It includes READ (total number of input reads), WRITTEN (reads retained after removing duplicates), EXCLUDED, EXAMINED, counts of PAIRED and SINGLE reads, as well as DUPLICATE SINGLE/PAIR and DUPLICATE TOTAL.
        - **Application**: Used to evaluate library complexity and duplication rate. A high WRITTEN/READ ratio indicates low duplication and good library complexity, while a low ratio suggests high PCR duplication or low-complexity sequencing.
        
		  <img width="318" height="104" alt="图片" src="https://github.com/user-attachments/assets/938f588d-b04b-4aa8-aa58-97fcc8aa5a39" />

      - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample:
	  
          <img width="1604" height="600" alt="图片" src="https://github.com/user-attachments/assets/ecde9e5e-17e7-404c-b85d-c3d9836c200b" />

        - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content':
    
          - Sequence counts for each sample. Estimate duplicate read counts:

            <img width="1595" height="557" alt="图片" src="https://github.com/user-attachments/assets/594acb6e-313c-4661-909b-bea2624e133b" />

          - Sequence Quality Histograms: The mean quality value across each base position in the read.
	  
            <img width="1601" height="612" alt="图片" src="https://github.com/user-attachments/assets/67f69cb2-1401-461c-850a-a40f397cbc0e" />

          - Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
	  
            <img width="1593" height="620" alt="图片" src="https://github.com/user-attachments/assets/20120b8e-675a-4f0c-be5c-8bc1dc34feac" />

        - **Samtools**: This module parses the output from samtools flagstat to report the percentage of total, mapped, and properly paired reads, providing a summary of alignment quality. Helps evaluate the effectiveness of deduplication and ensures that downstream analyses (e.g., peak calling, coverage profiling) are based on unique, non-redundant reads.
	  
          <img width="1634" height="779" alt="图片" src="https://github.com/user-attachments/assets/08f72013-30ee-4007-9dee-ea025246b8e9" />

          <img width="1626" height="781" alt="图片" src="https://github.com/user-attachments/assets/4bf9b33c-2ac4-49a7-bf77-2c0e6948e061" />

        - **Bowtie**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:
	  
          <img width="1625" height="507" alt="图片" src="https://github.com/user-attachments/assets/8da0c8a1-7efb-41b5-b47f-641df258307a" />

      - **`*BW_compare_PCA.pdf`**

        - **Content**: PDF file showing the principal component analysis (PCA) of BigWig signal profiles across multiple samples. It visualizes sample-to-sample similarity and variance based on genome-wide coverage or signal intensities.
        - **Application**: Used to assess the overall relationship between samples, detect outliers, and evaluate batch effects or experimental reproducibility.

          <img width="707" height="705" alt="图片" src="https://github.com/user-attachments/assets/e1ed3bc8-d562-4dfd-bf8b-43d64fca4971" />

      - **`*BW_compare_cor.pdf`**

        - **Content**: PDF file showing a heatmap of pairwise correlations between samples based on BigWig signal profiles. It typically includes correlation values (Pearson) and visually represents sample similarity across the genome.
        - **Application**: Used to assess consistency and reproducibility between samples, identify outliers, and evaluate experimental quality in End-seq.

          <img width="714" height="660" alt="图片" src="https://github.com/user-attachments/assets/9beacc3e-35bc-4315-8316-2708de675170" />

      - **`*fingerprints.pdf`**

        - **Content**: PDF file generated by `plotFingerprint` showing the cumulative read coverage across the genome for each BAM file. It visualizes enrichment patterns and sequencing depth consistency among samples.
        - **Application**: In the fingerprint plot, a larger separation between treatment and control curves, indicates stronger enrichment and higher signal-to-noise ratio. Also, the fingerprint plot can help decide whether to call narrow peaks or broad peaks: Narrow peaks are appropriate when the signal is sharp and localized; Broad peaks are used when the signal spans wide genomic regions with diffuse enrichment, such as histone modifications.

          <img width="678" height="511" alt="图片" src="https://github.com/user-attachments/assets/9e338721-680c-452c-8d7f-cd4e4f31de14" />

      - **`*..macs3.stats`**

        - **Content**: Contains summary statistics from MACS3 peak calling, including number of input reads, effective genome size, estimated fragment size, number of peaks called, and other runtime information.
        - **Application**: Used to check if MACS3 ran successfully and to detect any errors or warnings during the peak calling process.

      - **`*_peaks.narrowPeak`**

        - **Content**: BED6+4 format file containing peak locations along with peak summit, p-value, and q-value.  
    Suitable for direct loading into the UCSC Genome Browser when the `--trackline` option is enabled.

        - **Columns**:

        | Column | Description |
        |--------|------------|
        | chrom  | Chromosome name |
        | start  | Start position of the peak (0-based) |
        | end    | End position of the peak (not inclusive) |
        | name   | Peak name or ID |
        | score  | Integer score for display, calculated as `int(-10*log10(pvalue))` or `int(-10*log10(qvalue))` depending on whether `-p` or `-q` was used as the cutoff. |
        | strand | Strand information (‘+’, ‘-’, or ‘.’ if not applicable) |
        | signalValue | Fold enrichment at the peak summit |
        | pValue | -log10 p-value at the peak summit |
        | qValue | -log10 q-value (FDR) at the peak summit |
        | summit | Relative summit position to the peak start |

        - **Application**: Represents high-resolution binding sites and is widely used for downstream visualization and integrative genomic analyses.

          <img width="673" height="305" alt="图片" src="https://github.com/user-attachments/assets/2379edbf-596f-4832-b26d-867317443cc9" />

      - **`*_peaks.broadPeak`**

        - **Content**: BED6+3 format file (similar to narrowPeak, but without the 10th column for peak summits). Only available when `--broad` is enabled. In broad peak mode, the peak summit isn’t called, so the 5th, 7th–9th columns are the mean values across the peak region. Can be loaded directly into UCSC Genome Browser with `--trackline`.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom  | Chromosome name |
        | start  | Start position of the broad peak (0-based) |
        | end    | End position of the broad peak (not inclusive) |
        | name   | Peak name or ID |
        | score  | Mean score across the broad peak (similar to narrowPeak 5th column) |
        | strand | Strand information (‘+’, ‘-’, or ‘.’ if not applicable) |
        | signalValue | Mean enrichment signal across the peak |
        | pValue | Mean -log10 p-value across the peak |
        | qValue | Mean -log10 q-value (FDR) across the peak |

        - **Application**: Used for visualization DNA breadk signals, with bigwig file.

          example of H3K27ac.vs.Input_peaks.broadPeak: 
          <img width="687" height="310" alt="图片" src="https://github.com/user-attachments/assets/4297588f-10b9-4093-981b-451a7fe2bf26" />

          example of H3K27ac.vs.Input_peaks.broadPeak on igv with corresponding bw files:
          <img width="957" height="327" alt="图片" src="https://github.com/user-attachments/assets/a8b5b074-4257-4d6e-ae04-edbc110dbcf2" />

      - **`*_peaks.gappedPeak`**

        - **Content**: BED12+3 format file containing broad regions and narrow peaks within them. Only available when `--broad` is enabled. Can be loaded into UCSC Genome Browser. Columns 5, 7–9 may need adjustment if integrating with narrowPeak conventions.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom       | Chromosome name |
        | start       | Start of the broad region (0-based) |
        | end         | End of the broad region (not inclusive) |
        | name        | Peak name or ID |
        | score       | Score for display in UCSC browser (grey levels, similar to narrowPeak 5th column) |
        | strand      | Strand information (‘+’, ‘-’, or ‘.’) |
        | thickStart  | Start of the first narrow peak within the broad region |
        | thickEnd    | End of the first narrow peak within the broad region |
        | itemRgb     | RGB color for UCSC browser (0 uses default color) |
        | blockCount  | Number of blocks (including 1bp at start and end of broad regions) |
        | blockSizes  | Comma-separated lengths of each block |
        | blockStarts | Comma-separated start positions of each block relative to `start` |
        | foldChange  | Fold-change of enrichment within the peak |
        | -log10(pvalue) | -log10 p-value for the peak |
        | -log10(qvalue) | -log10 q-value (FDR) for the peak |

        - **Application**: Used to analyze subpeak structure, study internal peak features, or visualize complex enrichment patterns in broad regions.

          <img width="1101" height="301" alt="图片" src="https://github.com/user-attachments/assets/989f95f5-942f-4434-ab70-0b44a9333e3e" />

      - **`*_peaks.xls`**

        - **Content**: Tab-delimited summary of all peaks called by MACS3 with detailed metrics.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chr    | Chromosome name |
        | start  | Peak start position (0-based) |
        | end    | Peak end position (not inclusive) |
        | length | Peak length (end - start) |
        | pileup | Maximum pileup (number of overlapping tags) at the peak |
        | -log10(pvalue) | -log10 of p-value for peak significance |
        | fold_enrichment | Fold enrichment of the peak over background |
        | -log10(qvalue) | -log10 of q-value (FDR) for peak significance |
        | name   | Peak name or ID |

        - **Application**: Used for quantitative peak analysis, filtering peaks by significance, fold enrichment, or integrating with downstream functional annotation, motif analysis, and visualization.

      - **`*.peak.pdf`**

        - **Content**: Heatmap visualizing read enrichment over peaks. Generated using `plotHeatmap` from deepTools with a `*_peaks.broadPeak` and `*bw` input. The heatmap shows signal intensity (color-coded, viridis colormap) across all peaks, with missing data represented in white. The height and width of the heatmap are set for clear visualization of peak patterns.
        - **Application**: Used to assess global enrichment patterns across peaks. Peaks with strong enrichment appear as high-intensity bands or curves; if the signal is higher than control samples, it indicates that peak calling was successful and represents true biological enrichment.

          example: peak enrichment without control:
          
          <img width="269" height="769" alt="图片" src="https://github.com/user-attachments/assets/9090167c-4fc4-47a8-91da-12dca50029f1" />

          example: peak enrichment with control:

          <img width="487" height="776" alt="图片" src="https://github.com/user-attachments/assets/48e6e62f-4d5f-484f-9d28-6e327770e9a7" />

# Part V Video Tutorials

   Watch this video tutorial to see a complete walkthrough of running the pipeline:
