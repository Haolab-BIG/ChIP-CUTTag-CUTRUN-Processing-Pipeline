#!/bin/bash

print_help() {
  echo "Usage: $0 --sampleInfor <fastq_path> --comparisonInfor <compare_information> --outputdir <output_directory> --referencedir <reference_index_directory> --adapterFa <illumina_adapter> --sif <singularity_environment> --threads <number_of_threads> --binSize <bin_size> --g <genome_or_chromsize> --X <max_fragment_length> --peakcalling <MACS3_or_SEACR> --peaktype <broad_or_narrow> --peakPerSample <yes_or_no> --pval <p_value> --qval <q_value> --broad_cutoff <cutoff> --llocal <local_region> --keepdup <dup_setting> --nomodel <on_or_off> --nolambda <on_or_off> --callsummits <on_or_off> --extsize_val <extension_size> --shift_val <shift_bp> --seacr_threshold <percentile_cutoff>"

  echo "Options:"
  echo "  --sampleInfor       : (required) Tab-separated file with sample prefix and FASTQ path(s)."
  echo "  --comparisonInfor   : (optional) Tab-separated file listing treatment and control prefixes for comparisons."
  echo "  --outputdir         : (required) Directory to store all output files."
  echo "  --referencedir      : (required) Directory containing Bowtie2 reference index with prefix."
  echo "  --adapterFa         : (required) Path to adapter FASTA file."
  echo "  --sif               : (required) Path to Singularity environment (.sif) file."
  echo "  --threads           : (optional) Number of CPU threads to use (default: 8)."
  echo "  --binSize           : (optional) Bin size for bigWig generation (default: 10)."
  echo "  --g                 : (optional) If --peakcalling=MACS3, provide species code accepted by MACS3 (e.g. hs, mm, ce, dm). If --peakcalling=SEACR, provide path to the genome-specific chromatin.size file."
  echo "  --X                 : (optional) Bowtie2 maximum fragment length for paired-end mapping. Typical values: 700 for Cut&Tag/Cut&Run, 1000 for TF ChIP-seq, 2000 for histone-mark ChIP-seq."
  echo "  --peakcalling       : (optional) Peak calling method: MACS3 or SEACR (default: MACS3)."
  echo "  --peaktype          : (optional) For MACS3 use 'broad' or 'narrow'; for SEACR use 'relaxed' or 'stringent' (default: narrow)."
  echo "  --peakPerSample     : (optional) Call peaks per sample without control. Set to 'yes' if --comparisonInfor is not provided, or optionally 'yes' even when it is (default: yes)."
  echo "  --pval              : (optional) MACS3 p-value cutoff for narrow/strong peaks. If specified, overrides q-value."
  echo "  --qval              : (optional) MACS3 q-value (FDR) cutoff for narrow/strong peaks (default: 0.01)."
  echo "  --broad_cutoff      : (optional) MACS3 p- or q-value cutoff for broad peaks, used only when --peaktype=broad (default: 0.1)."
  echo "  --llocal            : (optional) MACS3 --llocal setting for samples without control (default: 100000)."
  echo "  --keepdup           : (optional) MACS3 --keep-dup setting (default: all)."
  echo "  --nomodel           : (optional) MACS3 disable model building: on/off (default: on)."
  echo "  --nolambda          : (optional) MACS3 disable dynamic lambda: on/off (default: on)."
  echo "  --callsummits       : (optional) MACS3 call summits: on/off (default: off)."
  echo "  --extsize_val       : (optional) MACS3 fragment/extension size (bp); effective only when --nomodel=on."
  echo "  --shift_val         : (optional) MACS3 read 5' end shift (bp); effective only when --nomodel=on. Positive shifts 5'→3'; negative shifts 3'→5'."
  echo "  --seacr_threshold   : (optional) SEACR percentile cutoff (e.g., 0.01 = top 1% of signal regions). Recommended: 0.01 for sharp marks (H3K4me3), 0.05 for broad marks (H3K27me3). Default: 0.01."
}

if [[ "$#" -eq 0 || "$1" == "--help" ]]; then
  print_help
  exit 0
fi

sampleInfor=""
comparisonInfor=""
outputdir=""
referencedir=""
adapterFa=""
sif=""
threads=8
binSize=10
g=hs
X=1000
peakcalling=MACS3
peaktype=narrow
peakPerSample=yes
pval=""
qval=""
broad_cutoff=0.1
llocal=100000
keepdup=all
nomodel=on
nolambda=on
callsummits=off
extsize_val=""
shift_val=""
seacr_threshold=0.01

while [[ $# -gt 0 ]]; do
  case $1 in
    --sampleInfor) 
      if [[ -z "$2" ||! -f "$2" ]]; then
        echo "Error: sample information file $2 not found!"
        exit 1
      fi
      sampleInfor="$2"; shift 2 ;;
    --comparisonInfor)
      if [[ -z "$2" ]]; then
        comparisonInfor=""
        shift 1
      else
        if [[ ! -f "$2" ]]; then
          echo "Error: comparison information file $2 not found!"
          exit 1
        fi
        comparisonInfor="$2"
        shift 2
      fi;;
    --outputdir) 
      if [[ -z "$2" ]]; then
        echo "Error: output directory $2 not found!"
        exit 1
      fi
      outputdir="$2"; shift 2 ;;
    --referencedir) 
      if [[ -z "$2" ]]; then
        echo "Error: reference directory $2 not found!"
        exit 1
      fi
      referencedir="$2"; shift 2 ;;
    --adapterFa) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: adapter fa file $2 not found!"
        exit 1
      fi
      adapterFa="$2"; shift 2 ;;
    --sif) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: singularity sif file $2 not found!"
        exit 1
      fi
      sif="$2"; shift 2 ;;
    --threads) 
      threads="$2"; shift 2 ;;
    --binSize) 
      binSize="$2"; shift 2 ;;
    --g) 
      g="$2"; shift 2 ;;
    --X) 
      X="$2"; shift 2 ;;
    --peakcalling) 
      peakcalling="$2"; shift 2 ;;
    --peaktype) 
      peaktype="$2"; shift 2 ;;
    --peakPerSample) 
      peakPerSample="$2"; shift 2 ;;
    --nolambda) 
      nolambda="$2"; shift 2 ;;
    --pval) 
      pval="$2"; shift 2 ;;
    --qval) 
      qval="$2"; shift 2 ;;
    --broad_cutoff) 
      broad_cutoff="$2"; shift 2 ;;
    --llocal) 
      llocal="$2"; shift 2 ;;
    --keepdup) 
      keepdup="$2"; shift 2 ;;
    --nomodel) 
      nomodel="$2"; shift 2 ;;
    --callsummits) 
      callsummits="$2"; shift 2 ;;
    --extsize_val) 
      extsize_val="$2"; shift 2 ;;
    --shift_val) 
      shift_val="$2"; shift 2 ;;
    --seacr_threshold) 
      seacr_threshold="$2"; shift 2 ;;
    *) 
      echo "Unknown option $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$sampleInfor" || -z "$outputdir" || -z "$referencedir" || -z "$adapterFa" || -z "$sif" ]]; then
  echo "Error: Missing required parameters"
  print_help
  exit 1
fi

########### mkdir
SING_EXEC1="singularity exec --cleanenv $sif"
$SING_EXEC1 mkdir -p ${outputdir}
$SING_EXEC1 mkdir -p ${outputdir}/qc
FastQCdir=${outputdir}/qc
$SING_EXEC1 mkdir -p ${outputdir}/trim
trimedfadir=${outputdir}/trim
$SING_EXEC1 mkdir -p ${outputdir}/bam
bowtieoutdir=${outputdir}/bam
$SING_EXEC1 mkdir -p ${outputdir}/bw
bwoutdir=${outputdir}/bw
$SING_EXEC1 mkdir -p ${outputdir}/figure
FigureDir=${outputdir}/figure
if [[ ${peakPerSample} == "yes" ]]; then
    $SING_EXEC1 mkdir -p ${outputdir}/peak.per.sample
    pcoutdir=${outputdir}/peak.per.sample
fi
$SING_EXEC1 mkdir -p ${outputdir}/multiqc
multiqcdir=${outputdir}/multiqc
if [[ -n ${comparisonInfor} && -f ${comparisonInfor} ]]; then
    $SING_EXEC1 mkdir -p ${outputdir}/peak.compare
    pkdir=${outputdir}/peak.compare
fi

########### singularity command

[[ -n $outputdir ]] && outputdir=$(readlink -f "$outputdir")
[[ -n $referencedir ]] && referencedir=$(readlink -f "$referencedir")
[[ -n $adapterFa ]] && adapterFa=$(readlink -f "$adapterFa")
[[ -n $sif ]] && sif=$(readlink -f "$sif")
if [[ -n $sampleInfor && -f $sampleInfor ]]; then
    while IFS=$'\t' read -r prefix r1 r2; do
        [[ $prefix == "Sample_prefix" || -z "$prefix" ]] && continue
        [[ -n "$r1" ]] && r1=$(readlink -f "$r1") && bind_dirs+=("$(dirname "$r1")")
        [[ -n "$r2" ]] && r2=$(readlink -f "$r2") && bind_dirs+=("$(dirname "$r2")")
    done < $sampleInfor
fi
bind_dirs=()
[[ -n $outputdir ]] && bind_dirs+=("$outputdir")
[[ -n $referencedir ]] && bind_dirs+=("$(dirname "$referencedir")")
[[ -n $adapterFa ]] && bind_dirs+=("$(dirname "$adapterFa")")
bind_dirs_unique=($(printf "%s\n" "${bind_dirs[@]}" | sort -u))
SING_EXEC="singularity exec --cleanenv"
for dir in "${bind_dirs_unique[@]}"; do
    real_dir=$(readlink -f "$dir")
    [[ -n $real_dir ]] && SING_EXEC+=" -B $real_dir:$real_dir"
done
SING_EXEC+=" $sif"

########### processing
while IFS=$'\t\r\n' read -r sample r1 r2; do
    echo "Sample: $sample"
    Sample=$(echo "$sample" | tr -d '\r' | xargs)
    r1=$(echo "$r1" | tr -d '\r' | xargs)
    r2=$(echo "$r2" | tr -d '\r' | xargs)
		## gz
    if [[ -n "$r1" && ! "$r1" =~ \.gz$ ]]; then
        echo "Compressing $r1 ..."
        $SING_EXEC gzip -c "$r1" > "${r1}.gz"
        r1="${r1}.gz"
    fi
    if [[ -n "$r2" && ! "$r2" =~ \.gz$ ]]; then
        echo "Compressing $r2 ..."
        $SING_EXEC gzip -c "$r2" > "${r2}.gz"
        r2="${r2}.gz"
    fi
    ## fastq path check
    #r1=$(readlink -f "$r1")
    #[[ -n "$r2" ]] && r2=$(readlink -f "$r2")
    #echo "  R1: $r1"
    #[[ -n "$r2" ]] && echo "  R2: $r2"
    echo "check: $r1"
    if [[ -f "$r1" ]]; then
        echo "File exists: $r1"
    else
        echo "File does not exist: $r1"
    fi
    #
    if [[ -z "$r2" ]]; then
		    # single-end
		    echo "run fastqc"
		    $SING_EXEC fastqc -o ${FastQCdir} -t ${threads} -q $r1
		    echo "cut adapters and sequence quality filters"
		    $SING_EXEC trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 $r1 -a file:${adapterFa} -o ${trimedfadir} --basename ${Sample} -j ${threads}
        echo "mapped the data to the genome"
        $SING_EXEC bowtie2 -t -k 1 --end-to-end --sensitive -p ${threads} -x ${referencedir} -U ${trimedfadir}/${Sample}_trimmed.fq.gz 2> ${bowtieoutdir}/${Sample}.bowtie.stats | samtools view -q 255 -bS - | samtools sort -n -@ ${threads} - -T $Sample -O BAM -o ${bowtieoutdir}/${Sample}.sort.bam
    else
    echo "check: $r2"
        if [[ -f "$r2" ]]; then
            echo "File exists: $r2"
        else
            echo "File does not exist: $r2"
        fi
        # pair-end
        echo "run fastqc"
		    $SING_EXEC fastqc -o ${FastQCdir} -t ${threads} -q $r1
		    $SING_EXEC fastqc -o ${FastQCdir} -t ${threads} -q $r2
		    echo "cut adapters and sequence quality filters"
		    $SING_EXEC bash -c "trim_galore -q 25 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 --paired '$r1' '$r2' -a file:${adapterFa} -o '${trimedfadir}' --basename '${Sample}' -j ${threads}"
		    echo "mapped the data to the genome"
		    $SING_EXEC bowtie2 -t -k 1 --end-to-end --sensitive -p ${threads} --fr --no-mixed --no-discordant -X $X -x ${referencedir} -1 ${trimedfadir}/${Sample}_val_1.fq.gz -2 ${trimedfadir}/${Sample}_val_2.fq.gz 2> ${bowtieoutdir}/${Sample}.bowtie.stats | samtools view -q 255 -bS - | samtools sort -n -@ ${threads} - -T $Sample -O BAM -o ${bowtieoutdir}/${Sample}.sort.bam
    fi
    echo "rm duplicate reads"
    $SING_EXEC samtools fixmate -r -m ${bowtieoutdir}/${Sample}.sort.bam -| samtools sort -@ ${threads} - | samtools markdup -r -s - ${bowtieoutdir}/${Sample}.DeDup.bam 2> ${bowtieoutdir}/${Sample}.markdup.log
    $SING_EXEC samtools flagstat ${bowtieoutdir}/${Sample}.DeDup.bam > ${bowtieoutdir}/${Sample}.flagstat.txt
    $SING_EXEC samtools index ${bowtieoutdir}/${Sample}.DeDup.bam
    rm ${bowtieoutdir}/${Sample}.sort.bam
    echo "obtain BW files"
    $SING_EXEC samtools index ${bowtieoutdir}/${Sample}.DeDup.bam
    $SING_EXEC bamCoverage -b ${bowtieoutdir}/${Sample}.DeDup.bam -o ${bwoutdir}/${Sample}.DeDup.bw  -p ${threads} --binSize ${binSize} --normalizeUsing CPM
    # peak calling for each sample
    if [[ $peakPerSample == "yes" ]]; then
	    echo "peak calling for each sample"
	    if [[ $peakcalling == "MACS3" ]]; then
		    echo "peak calling using MACS3"
		    extra_opts=""
		    [[ $nomodel == "on" ]] && extra_opts+=" --nomodel"
		    [[ $nolambda == "on" ]] && extra_opts+=" --nolambda"
		    [[ $callsummits == "on" ]] && extra_opts+=" --call-summits"
		    if [[ -z "$r2" ]]; then
			    bam_format=BAM
					echo "running with single-end bam"
		    else
					bam_format=BAMPE
					echo "running with paired-end bam"
		    fi
		    if [[ $nomodel == "on" ]]; then
		      [[ -n "$shift_val" ]] && extra_opts+=" --shift $shift_val"
		      [[ -n "$extsize_val" ]] && extra_opts+=" --extsize $extsize_val"
		    fi
		    if [[ "$peaktype" == "broad" ]]; then
		      extra_opts+=" --broad"
		      [[ -n "$broad_cutoff" ]] && extra_opts+=" --broad-cutoff $broad_cutoff"
		      [[ -n "$pval" ]] && extra_opts+=" -p $pval"
		      [[ -n "$qval" ]] && extra_opts+=" -q $qval"
		    elif [[ "$peaktype" == "narrow" ]]; then
		      [[ -n "$pval" ]] && extra_opts+=" -p $pval"
		      [[ -n "$qval" ]] && extra_opts+=" -q $qval"
		    fi
		    if [[ -z "$pval" && -z "$qval" ]]; then
		      extra_opts+=" -q 0.01"
		    fi
		    $SING_EXEC macs3 callpeak -t ${bowtieoutdir}/${Sample}.DeDup.bam -f $bam_format -g ${g} --outdir ${pcoutdir} -n ${Sample} $extra_opts --llocal ${llocal} --keep-dup ${keepdup} 2> >(tee -a ${pcoutdir}/${Sample}.macs3.stats >&2) 1> ${pcoutdir}/${Sample}.macs3.stdout
		    rm ${pcoutdir}/${Sample}.macs3.stdout
		    peak_file=${pcoutdir}/${Sample}_peaks.${peaktype}Peak
			elif [[ $peakcalling == "SEACR" ]]; then
		    echo "peak calling using SEACR"
		    echo "Generating bedgraph from BAM file..."
		    # normalizing bedgraph
		    if [[ -z "$r2" ]]; then
			    total_reads=$(samtools view -c ${bowtieoutdir}/${Sample}.DeDup.bam)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${Sample}.DeDup.bam -bg -g ${g} > ${pcoutdir}/${Sample}.bedgraph
		    else
			    total_reads=$(echo "$(samtools view ${bowtieoutdir}/${Sample}.DeDup.bam | wc -l) / 2" | bc)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${Sample}.DeDup.bam -bg -pc -g ${g} > ${pcoutdir}/${Sample}.bedgraph
		    fi
		    scaling_factor=$(echo "1000000 / $total_reads" | bc -l)
		    $SING_EXEC awk -v sf="$scaling_factor" '{print $1 "\t" $2 "\t" $3 "\t" $4 * sf}' ${pcoutdir}/${Sample}.bedgraph > ${pcoutdir}/${Sample}.normalized.bedgraph
        threshold=$seacr_threshold
        $SING_EXEC SEACR_1.3.sh ${pcoutdir}/${Sample}.normalized.bedgraph $threshold norm $peaktype ${pcoutdir}/${Sample}_seacr_top${threshold}
        peak_file=${pcoutdir}/${Sample}_seacr_top${threshold}.${peaktype}.bed
        rm ${pcoutdir}/${Sample}.normalized.bedgraph
        rm ${pcoutdir}/${Sample}.bedgraph
			else
		    echo "Error: peakcalling must be MACS3 or SEACR."
			fi
			echo "peak enrichment"
			if [[ -s "$peak_file" ]] && [[ $(wc -l < "$peak_file") -gt 0 ]]; then
		    if [[ $peaktype == "broad" ]]; then
	        $SING_EXEC computeMatrix scale-regions -p ${threads} -S ${bwoutdir}/${Sample}.DeDup.bw \
                              -R $peak_file -a 5000 -b 5000 --regionBodyLength 1000 \
                              --missingDataAsZero -o ${FigureDir}/${Sample}.peak.gz
			    $SING_EXEC plotHeatmap -m ${FigureDir}/${Sample}.peak.gz -out ${FigureDir}/${Sample}.peak.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4  --startLabel "Start" --endLabel "End"
		    elif [[ "$peaktype" == "narrow" || "$peaktype" == "stringent" || "$peaktype" == "relaxed" ]]; then
	        $SING_EXEC computeMatrix scale-regions -p ${threads} -S ${bwoutdir}/${Sample}.DeDup.bw \
                              -R $peak_file -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${FigureDir}/${Sample}.peak.gz
				  $SING_EXEC plotHeatmap -m ${FigureDir}/${Sample}.peak.gz -out ${FigureDir}/${Sample}.peak.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4 --refPointLabel "Peak"
		    fi
		  else
		    echo "Skip: $peak_file is missing or empty."
		  fi
		else
	    echo ""
		fi
done < <(tail -n +2 "$sampleInfor")

########### correlation
echo "Check the correlation"
bwfiles=("${bwoutdir}"/*.bw)
labels=($(basename -s .bw "${bwfiles[@]}"))
$SING_EXEC multiBigwigSummary bins -b "${bwfiles[@]}" --labels "${labels[@]}" -out ${FigureDir}/BW_compare_PCA.npz -p ${threads}
$SING_EXEC plotPCA -in ${FigureDir}/BW_compare_PCA.npz -o ${FigureDir}/BW_compare_PCA.pdf
$SING_EXEC plotCorrelation -in ${FigureDir}/BW_compare_PCA.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers
$SING_EXEC rm ${FigureDir}/BW_compare_PCA.npz

echo "check the seq quality"
bamfiles=("${bowtieoutdir}"/*.DeDup.bam)
$SING_EXEC plotFingerprint -b "${bamfiles[@]}" --labels "${labels[@]}" --skipZeros --plotFile ${FigureDir}/fingerprints.pdf

########### peak compare
echo "peak calling with comparison"

while IFS=$'\t' read -r treatment control; do
    treatment=$(echo "$treatment" | tr -d '\r' | xargs)
    control=$(echo "$control" | tr -d '\r' | xargs)
    echo "Treatment: ${treatment}; Control: ${control}"
    treat_bam="${bowtieoutdir}/${treatment}.DeDup.bam"
    ctrl_bam="${bowtieoutdir}/${control}.DeDup.bam"
    echo ${treat_bam}
    echo ${ctrl_bam}
    r2_path=$($SING_EXEC bash -c "awk -F '\t' -v sample=\"$treatment\" '\$1==sample {print \$3}' \"$sampleInfor\" | tr -d '\r' | xargs")
    if [[ -f "$treat_bam" && -f "$ctrl_bam" ]]; then
        echo "Calling peaks for $treatment vs $control ..."
        if [[ $peakcalling == "MACS3" ]]; then
		    echo "peak calling using MACS3"
		    extra_opts=""
		    [[ $nomodel == "on" ]] && extra_opts+=" --nomodel"
		    [[ $callsummits == "on" ]] && extra_opts+=" --call-summits"
		    if [[ -z "$r2_path" ]]; then
					bam_format=BAM
					echo "running with single-end bam"
		    else
					bam_format=BAMPE
					echo "running with paired-end bam"
		    fi
		    if [[ $nomodel == "on" ]]; then
		      [[ -n "$shift_val" ]] && extra_opts+=" --shift $shift_val"
		      [[ -n "$extsize_val" ]] && extra_opts+=" --extsize $extsize_val"
		    fi
		    if [[ "$peaktype" == "broad" ]]; then
		    	extra_opts+=" --broad"
		    	[[ -n "$broad_cutoff" ]] && extra_opts+=" --broad-cutoff $broad_cutoff"
		    	[[ -n "$pval" ]] && extra_opts+=" -p $pval"
		     	[[ -n "$qval" ]] && extra_opts+=" -q $qval"
		    elif [[ "$peaktype" == "narrow" ]]; then
		    	[[ -n "$pval" ]] && extra_opts+=" -p $pval"
		    	[[ -n "$qval" ]] && extra_opts+=" -q $qval"
		    fi
		    if [[ -z "$pval" && -z "$qval" ]]; then
		      extra_opts+=" -q 0.01"
		    fi
		    echo $extra_opts
		    $SING_EXEC macs3 callpeak -t ${treat_bam} -c ${ctrl_bam} -f ${bam_format} -g ${g} --outdir ${pkdir} -n ${treatment}.vs.${control} $extra_opts --keep-dup ${keepdup} 2> >(tee -a ${pkdir}/${treatment}.vs.${control}.macs3.stats >&2) 1> ${pkdir}/${treatment}.vs.${control}.macs3.stdout
		    rm ${pkdir}/${treatment}.vs.${control}.macs3.stdout
		    peak_fileC=${pkdir}/${treatment}.vs.${control}_peaks.${peaktype}Peak
		elif [[ $peakcalling == "SEACR" ]]; then
		    echo "peak calling using SEACR"
		    echo "Generating bedgraph from BAM file..."
		    # normalizing bedgraph
		    if [[ -z "$r2_path" ]]; then
		    	echo "running with single-end bam"
			    total_readsT=$(samtools view -c ${bowtieoutdir}/${treatment}.DeDup.bam)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${treatment}.DeDup.bam -bg -g ${g} > ${pkdir}/${treatment}.bedgraph
			    total_readsC=$(samtools view -c ${bowtieoutdir}/${control}.DeDup.bam)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${control}.DeDup.bam -bg -g ${g} > ${pkdir}/${control}.bedgraph
		    else
					echo "running with paired-end bam"
			    total_readsT=$(echo "$(samtools view ${bowtieoutdir}/${treatment}.DeDup.bam | wc -l) / 2" | bc)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${treatment}.DeDup.bam -bg -pc -g ${g} > ${pkdir}/${treatment}.bedgraph
			    total_readsC=$(echo "$(samtools view ${bowtieoutdir}/${control}.DeDup.bam | wc -l) / 2" | bc)
			    $SING_EXEC bedtools genomecov -ibam ${bowtieoutdir}/${control}.DeDup.bam -bg -pc -g ${g} > ${pkdir}/${control}.bedgraph
		    fi
		    scaling_factorT=$(echo "1000000 / $total_readsT" | bc -l)
		    $SING_EXEC awk -v sf="$scaling_factorT" '{print $1 "\t" $2 "\t" $3 "\t" $4 * sf}' ${pkdir}/${treatment}.bedgraph > ${pkdir}/${treatment}.normalized.bedgraph
		    scaling_factorC=$(echo "1000000 / $total_readsC" | bc -l)
		    $SING_EXEC awk -v sf="$scaling_factorC" '{print $1 "\t" $2 "\t" $3 "\t" $4 * sf}' ${pkdir}/${control}.bedgraph > ${pkdir}/${control}.normalized.bedgraph
			threshold=$seacr_threshold
			$SING_EXEC SEACR_1.3.sh ${pkdir}/${treatment}.normalized.bedgraph ${pkdir}/${control}.normalized.bedgraph norm $peaktype ${pkdir}/${treatment}.vs.${control}_seacr
			peak_fileC=${pkdir}/${treatment}.vs.${control}_seacr.${peaktype}.bed
			rm ${pkdir}/${control}.normalized.bedgraph
			rm ${pkdir}/${treatment}.normalized.bedgraph
			rm ${pkdir}/${treatment}.bedgraph
			rm ${pkdir}/${control}.bedgraph
		else
			echo "Error: peakcalling must be MACS3 or SEACR."
		fi
		echo "peak enrichment"
		if [[ -s "$peak_fileC" ]] && [[ $(wc -l < "$peak_fileC") -gt 0 ]]; then
		  if [[ $peaktype == "broad" ]]; then
			  $SING_EXEC computeMatrix scale-regions -p ${threads} -S ${bwoutdir}/${treatment}.DeDup.bw ${bwoutdir}/${control}.DeDup.bw \
                              -R $peak_fileC -a 5000 -b 5000 --regionBodyLength 1000 \
                              --missingDataAsZero -o ${FigureDir}/${treatment}.vs.${control}.peak.gz
			  $SING_EXEC plotHeatmap -m ${FigureDir}/${treatment}.vs.${control}.peak.gz -out ${FigureDir}/${treatment}.vs.${control}.peak.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4  --startLabel "Start" --endLabel "End"
		  elif [[ "$peaktype" == "narrow" || "$peaktype" == "stringent" || "$peaktype" == "relaxed" ]]; then
			  $SING_EXEC computeMatrix scale-regions -p ${threads} -S ${bwoutdir}/${treatment}.DeDup.bw ${bwoutdir}/${control}.DeDup.bw \
                              -R $peak_fileC -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${FigureDir}/${treatment}.vs.${control}.peak.gz
			  $SING_EXEC plotHeatmap -m ${FigureDir}/${treatment}.vs.${control}.peak.gz -out ${FigureDir}/${treatment}.vs.${control}.peak.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4 --refPointLabel "Peak"
		  fi
		else
		  echo "Skip: $peak_file is missing or empty."
		fi
	else
        echo "Warning: Missing BAM file for $treatment or $control"
    fi
done < <(tail -n +2 "$comparisonInfor")

########### multiqc
$SING_EXEC multiqc ${FastQCdir}/* ${trimedfadir}/* ${bowtieoutdir}/*.bowtie.stats ${bowtieoutdir}/*flagstat.txt -o ${multiqcdir} --force

########### remove
$SING_EXEC rm -r ${FastQCdir}
$SING_EXEC rm -r ${trimedfadir}
$SING_EXEC rm -f ${FigureDir}/*.peak.gz


echo "Finished."


