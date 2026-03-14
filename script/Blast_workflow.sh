#!/usr/bin/env bash
set -e

# Setting out Variables

ACC_NUM="GCF_024496165.1 GCF_011045595.1 GCF_040047695.1 GCF_046619595.1 GCF_052906805.1"
K_ZIP="data/kleb_genome/scripted_kleb_dataset.zip"
KLEB_DATASET_EXTRACT="data/kleb_genome_dataset"
KAG="data/scripted_kleb_assembled_genome.fasta"
KG="data/kleb_genome"
CD="data/card_db"
CARD_TAR="data/card_db/broadstreet-v3.3.0.tar.bz2"
CARD_FASTA="data/card_db/scripted_nucleotide_fasta_protein_homolog_model.fasta"
ARO_INDEX="data/card_db/scripted_aro_index.tsv"
DB="data/card_db/databases/scripted_card_nucl_db"

BO="analysis/scripted_blast_output"
FINAL="analysis/scripted_blast_result.tsv"
CLASSIFIED="analysis/scripted_classified_AMR_result.tsv"
FINAL_SORTED="analysis/scripted_finally_classified_sorted_AMR_result.tsv"

# Creating Directories
cd "$HOME/Bioinfo_UI/Blast_assessment"
mkdir -p data analysis script
mkdir -p "KG" "$CD" "$CD/databases" "$BO"

# Downloading the % assembled genomes
if [ -f data/datasets.exe ]; then
  data/datasets.exe download genome accession $ACC_NUM --filename "$K_ZIP" && unzip -o -q "$K_ZIP" -d "$KLEB_DATASET_EXTRACT"
else
  echo "ERROR: ../data/datasets.exe not found"
fi

# Unzipping any .gz genomes if present, then copy assembled genomes into $KAG
find "$KLEB_DATASET_EXTRACT/ncbi_dataset/data" -name "*genomic.fna.gz" -exec gunzip -f {} \; 
cat "$KLEB_DATASET_EXTRACT/ncbi_dataset/data"/*/*_genomic.fna > "$KAG"
echo "done unzipping"

# 2) Downloading Card Database

curl -o "$CARD_TAR" https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2

# Extract card database and copy files into fixed paths.
tar -xvf "$CARD_TAR"
cp "$(find "$CD" -name nucleotide_fasta_protein_homolog_model.fasta | head -n 1)" "$CARD_FASTA"
cp "$(find "$CD" -name aro_index.tsv | head -n 1)" "$ARO_INDEX"
echo "done copying"

# 3) Creating Blast Database
export PATH="$PATH:$HOME/Bioinfo_UI/Blast_assessment/data/ncbi-blast-2.17.0+/bin"
makeblastdb -in "$CARD_FASTA" -title scripted_card_nucl_db -dbtype nucl -out "$DB"
echo "done creating database"

# 4) Extracting relevant columns inside aro_index tsv file and adding a new header
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "ARO_Index","Variant","Family","Resistance_Class"; next} {print $1,$6,$9,$10}' "$ARO_INDEX" > analysis/scripted_aro_reference.tsv
echo "done extraction"

# 5) Running Blast for the 5 assembled genomes, Filtering %identity, Splitting the ARO_Index and Gene column, Then Mapping the Ref tsv + Classifica

# Blasting & Filtering
blastn -query "$KAG" -db "$DB" -evalue 1e-30 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
| awk -F'\t' 'BEGIN{OFS="\t"} $3>=90 {n=split($2,a,"|"); print $1,a[n-1],a[n],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > "$FINAL"
echo "done blasting"

# Mapping and AMR Classification using the aro index tabl
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR{m[$1]=$2"\t"$3"\t"$4; next} {aro=$2; add=(aro in m)?m[aro]:"NA\tNA\tNA"; print $0, add}
' analysis/scripted_aro_reference.tsv "$FINAL" > "$CLASSIFIED"
echo "done classifying"

# Appending and Sorting relevent columns
echo -e "ID\tGene_Variants\tGene_Families\tResistance_Class" > "$FINAL_SORTED"
awk 'NR>1' "$CLASSIFIED" | cut -f1,14,15,16 | sort -t$'\t' -k2,2 >> "$FINAL_SORTED"

echo "Done sorting: $FINAL_SORTED"
