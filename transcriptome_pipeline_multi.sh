#!/bin/bash
set -e
set -o pipefail

VERSION="1.1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THREADS=32
MODE="illumina"
PLOIDY=1
MIN_QUAL=20
MIN_DEPTH=10
SKIP_FILTER=false
MULTI_INPUT=""

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Transcriptome Assembly Pipeline v${VERSION} (Multi-sample support)

Usage: $0 -i <input1[,input2,input3,...]> -r <reference> -o <output_prefix> [options]

Required:
    -i  Input fastq file(s) - comma-separated for multiple inputs
        Single: -i sample1.fastq
        Multiple: -i sample1.fastq,sample2.fastq,sample3.fastq
    -r  Reference genome fasta
    -o  Output prefix/sample name

Optional:
    -m  Mode: illumina (default), nanopore, pacbio
    -t  Threads (default: 32)
    -p  Ploidy (default: 1)
    -q  Min quality (default: 20)
    -d  Min depth (default: 10)
    -s  Skip quality filtering (use raw data directly)
    -h  Help

Multi-sample workflow:
    1. Each input file is processed independently (filter + map)
    2. All BAM files are merged using samtools merge
    3. Merged BAM is used for transcript assembly and downstream analysis

Example:
    # Single sample
    $0 -i 01_raw_data/sample.fastq -r 03_ref/genome.fna -o sample1

    # Multiple samples (will be merged)
    $0 -i 01_raw_data/s1.fastq,01_raw_data/s2.fastq,01_raw_data/s3.fastq \\
       -r 03_ref/genome.fna -o merged_sample -m nanopore -s

EOF
    exit 1
}

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) MULTI_INPUT="$2"; shift 2 ;;
        -r) REFERENCE="$2"; shift 2 ;;
        -o) OUTPUT_PREFIX="$2"; shift 2 ;;
        -m) MODE="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -p) PLOIDY="$2"; shift 2 ;;
        -q) MIN_QUAL="$2"; shift 2 ;;
        -d) MIN_DEPTH="$2"; shift 2 ;;
        -s) SKIP_FILTER=true; shift ;;
        -h) usage ;;
        *) log_error "Unknown option: $1" ;;
    esac
done

[ -z "$MULTI_INPUT" ] || [ -z "$REFERENCE" ] || [ -z "$OUTPUT_PREFIX" ] && usage
[ ! -f "$REFERENCE" ] && log_error "Reference not found: $REFERENCE"

# Split input by comma
IFS=',' read -ra INPUT_ARRAY <<< "$MULTI_INPUT"
NUM_INPUTS=${#INPUT_ARRAY[@]}

log_info "Starting pipeline for ${OUTPUT_PREFIX}"
log_info "Mode: $MODE, Threads: $THREADS, Ploidy: $PLOIDY, Skip filter: $SKIP_FILTER"
log_info "Number of input files: $NUM_INPUTS"

# Validate all input files exist
for INPUT_FILE in "${INPUT_ARRAY[@]}"; do
    [ ! -f "$INPUT_FILE" ] && log_error "Input not found: $INPUT_FILE"
    log_info "  - $INPUT_FILE"
done

mkdir -p 02_clean 04_mapping

# Process each input file
BAM_FILES=()
for i in "${!INPUT_ARRAY[@]}"; do
    INPUT="${INPUT_ARRAY[$i]}"
    SAMPLE_NAME="${OUTPUT_PREFIX}_sample$((i+1))"

    log_info "========================================"
    log_info "Processing sample $((i+1))/$NUM_INPUTS: $(basename $INPUT)"
    log_info "========================================"

    # Step 1: Filter reads (or skip)
    if [ "$SKIP_FILTER" = true ]; then
        log_info "Step 1.$((i+1)): Skipping quality filtering (using raw data)"
        if [[ "$INPUT" == *.gz ]]; then
            cp "$INPUT" "02_clean/${SAMPLE_NAME}.fq.gz"
        else
            gzip -c "$INPUT" > "02_clean/${SAMPLE_NAME}.fq.gz"
        fi
        CLEAN_R1="02_clean/${SAMPLE_NAME}.fq.gz"
        CLEAN_R2=""
    else
        log_info "Step 1.$((i+1)): Quality filtering"

        case $MODE in
            illumina)
                R1="$INPUT"
                R2="${INPUT/_R1/_R2}"
                R2="${R2/_1./_2.}"

                if [ -f "$R2" ] && [ "$R1" != "$R2" ]; then
                    fastp -i "$R1" -I "$R2" \
                          -o "02_clean/${SAMPLE_NAME}_R1.fq.gz" \
                          -O "02_clean/${SAMPLE_NAME}_R2.fq.gz" \
                          -q $MIN_QUAL -l 50 -w $THREADS \
                          -j "02_clean/${SAMPLE_NAME}.json" \
                          -h "02_clean/${SAMPLE_NAME}.html"
                    CLEAN_R1="02_clean/${SAMPLE_NAME}_R1.fq.gz"
                    CLEAN_R2="02_clean/${SAMPLE_NAME}_R2.fq.gz"
                else
                    fastp -i "$R1" -o "02_clean/${SAMPLE_NAME}.fq.gz" \
                          -q $MIN_QUAL -l 50 -w $THREADS \
                          -j "02_clean/${SAMPLE_NAME}.json" \
                          -h "02_clean/${SAMPLE_NAME}.html"
                    CLEAN_R1="02_clean/${SAMPLE_NAME}.fq.gz"
                    CLEAN_R2=""
                fi
                ;;
            nanopore|pacbio)
                if command -v NanoFilt &> /dev/null; then
                    cat "$INPUT" | NanoFilt -q $MIN_QUAL -l 200 | gzip > "02_clean/${SAMPLE_NAME}.fq.gz"
                else
                    if [[ "$INPUT" == *.gz ]]; then
                        cp "$INPUT" "02_clean/${SAMPLE_NAME}.fq.gz"
                    else
                        gzip -c "$INPUT" > "02_clean/${SAMPLE_NAME}.fq.gz"
                    fi
                fi
                CLEAN_R1="02_clean/${SAMPLE_NAME}.fq.gz"
                CLEAN_R2=""
                ;;
        esac
    fi

    # Step 2: Mapping
    log_info "Step 2.$((i+1)): Mapping to reference"

    case $MODE in
        illumina)
            [ ! -f "${REFERENCE}.1.ht2" ] && hisat2-build "$REFERENCE" "$REFERENCE" -p $THREADS
            if [ -n "$CLEAN_R2" ]; then
                hisat2 -x "$REFERENCE" -1 "$CLEAN_R1" -2 "$CLEAN_R2" -p $THREADS --dta -S "04_mapping/${SAMPLE_NAME}.sam"
            else
                hisat2 -x "$REFERENCE" -U "$CLEAN_R1" -p $THREADS --dta -S "04_mapping/${SAMPLE_NAME}.sam"
            fi
            ;;
        nanopore|pacbio)
            minimap2 -ax splice -uf -k14 -t $THREADS "$REFERENCE" "$CLEAN_R1" > "04_mapping/${SAMPLE_NAME}.sam"
            ;;
    esac

    samtools view -@ $THREADS -bS "04_mapping/${SAMPLE_NAME}.sam" | samtools sort -@ $THREADS -o "04_mapping/${SAMPLE_NAME}.bam"
    samtools index -@ $THREADS "04_mapping/${SAMPLE_NAME}.bam"
    rm "04_mapping/${SAMPLE_NAME}.sam"
    samtools flagstat "04_mapping/${SAMPLE_NAME}.bam" > "04_mapping/${SAMPLE_NAME}_flagstat.txt"

    BAM_FILES+=("04_mapping/${SAMPLE_NAME}.bam")
    log_info "Completed processing sample $((i+1)): ${SAMPLE_NAME}.bam"
done

# Step 3: Merge BAM files (if multiple inputs)
log_info "========================================"
if [ $NUM_INPUTS -gt 1 ]; then
    log_info "Step 3: Merging $NUM_INPUTS BAM files"
    samtools merge -@ $THREADS "04_mapping/${OUTPUT_PREFIX}_merged.bam" "${BAM_FILES[@]}"
    samtools index -@ $THREADS "04_mapping/${OUTPUT_PREFIX}_merged.bam"
    samtools flagstat "04_mapping/${OUTPUT_PREFIX}_merged.bam" > "04_mapping/${OUTPUT_PREFIX}_merged_flagstat.txt"
    FINAL_BAM="04_mapping/${OUTPUT_PREFIX}_merged.bam"
    log_info "Merged BAM created: $FINAL_BAM"
else
    log_info "Step 3: Single input, skipping merge"
    FINAL_BAM="${BAM_FILES[0]}"
    ln -sf "$(basename $FINAL_BAM)" "04_mapping/${OUTPUT_PREFIX}.bam"
    FINAL_BAM="04_mapping/${OUTPUT_PREFIX}.bam"
fi

# Step 4: StringTie assembly
log_info "Step 4: Transcript assembly"
mkdir -p 05_stringtie
stringtie "$FINAL_BAM" -o "05_stringtie/${OUTPUT_PREFIX}.gtf" -p $THREADS -v

# Step 5: Variant calling
log_info "Step 5: Variant calling"
mkdir -p 06_variant
[ ! -f "${REFERENCE}.fai" ] && samtools faidx "$REFERENCE"

bcftools mpileup -f "$REFERENCE" "$FINAL_BAM" -Ou | \
bcftools call -mv --ploidy $PLOIDY -Oz -o "06_variant/${OUTPUT_PREFIX}.vcf.gz"

bcftools index "06_variant/${OUTPUT_PREFIX}.vcf.gz"
bcftools view -v snps -i "QUAL>=${MIN_QUAL} && DP>=${MIN_DEPTH}" "06_variant/${OUTPUT_PREFIX}.vcf.gz" -Oz -o "06_variant/${OUTPUT_PREFIX}_filtered.vcf.gz"
bcftools index "06_variant/${OUTPUT_PREFIX}_filtered.vcf.gz"
bcftools stats "06_variant/${OUTPUT_PREFIX}_filtered.vcf.gz" > "06_variant/${OUTPUT_PREFIX}_stats.txt"

# Step 6: Consensus genome
log_info "Step 6: Generating consensus"
mkdir -p 07_consensus
bcftools consensus -f "$REFERENCE" "06_variant/${OUTPUT_PREFIX}_filtered.vcf.gz" > "07_consensus/${OUTPUT_PREFIX}_consensus.fasta"

# Step 7: Extract transcripts
log_info "Step 7: Extracting transcripts"
mkdir -p 08_output_transcript

if command -v gffread &> /dev/null; then
    gffread -w "08_output_transcript/${OUTPUT_PREFIX}_transcripts.fasta" -g "07_consensus/${OUTPUT_PREFIX}_consensus.fasta" "05_stringtie/${OUTPUT_PREFIX}.gtf"
else
    awk '$3=="exon"' "05_stringtie/${OUTPUT_PREFIX}.gtf" | \
        awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t0\t"$7}' | tr -d '";' | sort -k1,1 -k2,2n > "08_output_transcript/${OUTPUT_PREFIX}.bed"
    bedtools getfasta -fi "07_consensus/${OUTPUT_PREFIX}_consensus.fasta" -bed "08_output_transcript/${OUTPUT_PREFIX}.bed" -s -name > "08_output_transcript/${OUTPUT_PREFIX}_transcripts.fasta"
fi

awk '/^>/ {if (seq && length(seq) >= 200) print header"\n"seq; header=$0; seq=""; next} {seq=seq$0} END {if (length(seq) >= 200) print header"\n"seq}' \
    "08_output_transcript/${OUTPUT_PREFIX}_transcripts.fasta" > "08_output_transcript/${OUTPUT_PREFIX}_transcripts_filtered.fasta"

# Step 8: Protein prediction
log_info "Step 8: Protein prediction"
mkdir -p 09_proteins
cd 09_proteins
TransDecoder.LongOrfs -t "../08_output_transcript/${OUTPUT_PREFIX}_transcripts_filtered.fasta" -m 100
TransDecoder.Predict -t "../08_output_transcript/${OUTPUT_PREFIX}_transcripts_filtered.fasta" --cpu $THREADS
mv "${OUTPUT_PREFIX}_transcripts_filtered.fasta.transdecoder.pep" "${OUTPUT_PREFIX}_proteins.fasta" 2>/dev/null || true
mv "${OUTPUT_PREFIX}_transcripts_filtered.fasta.transdecoder.cds" "${OUTPUT_PREFIX}_cds.fasta" 2>/dev/null || true
cd ..

log_info "Pipeline completed for ${OUTPUT_PREFIX}!"

# Report
NTRANS=$(grep -c '^>' "08_output_transcript/${OUTPUT_PREFIX}_transcripts_filtered.fasta" 2>/dev/null || echo 0)
NPROT=$(grep -c '^>' "09_proteins/${OUTPUT_PREFIX}_proteins.fasta" 2>/dev/null || echo 0)
echo ""
echo "========================================"
echo "           Final Summary"
echo "========================================"
echo "Sample: ${OUTPUT_PREFIX}"
echo "Input files processed: $NUM_INPUTS"
if [ $NUM_INPUTS -gt 1 ]; then
    echo "Merged BAM: $FINAL_BAM"
fi
echo "Transcripts: $NTRANS"
echo "Proteins: $NPROT"
echo "========================================"
