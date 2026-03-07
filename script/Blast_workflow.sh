#!/usr/bin/env bash
set -e

# Setting out Variables

ACC_NUM="GCF_024496165.1 GCF_011045595.1 GCF_040047695.1 GCF_046619595.1 GCF_052906805.1"

K_ZIP="data/kleb_genome/scripted_kleb.zip"
KLEB_EXTRACT="data/kleb_genome/scripted_kleb"
KAG="data/scripted_kleb_assembled_genome"

CD="data/card_db"
CARD_TAR="data/card_db/broadstreet-v3.3.0.tar.bz2"
CARD_FASTA="data/card_db/scripted_nucleotide_fasta_protein_homolog_model.fasta"
ARO_INDEX="data/card_db/scripted_aro_index.tsv"
DB="data/card_db/databases/scripted_card_nucl_db"

BO="analysis/scripted_blast_output"
FINAL="analysis/scripted_final_result"
CLASSIFIED="analysis/scripted_AMR_result_classified"
FINAL_SORTED="analysis/scripted_AMR_finally_classified"

# Creating Directories
cd "$HOME/Bioinfo_UI/Blast_assessment"
mkdir -p data analysis script
mkdir -p data/kleb_genome "$KAG" "$CD" "$CD/databases" "$BO" "$FINAL" "$CLASSIFIED" "$FINAL_SORTED"

# Downloading ncbi dataset tools
[ -f data/datasets.exe ] || curl -o data/datasets.exe https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe
chmod +x data/datasets.exe 2>/dev/null || true

# 1) Dowloading the 5 assembled genomes
./data/datasets.exe download genome accession $ACC_NUM --filename "$K_ZIP" 
unzip -o -q "$K_ZIP" -d "$KLEB_EXTRACT"

# Unzipping any .gz genomes if present, then copy assembled genomes into $KAG
find "$KLEB_EXTRACT/ncbi_dataset/data" -name "*genomic.fna.gz" -exec gunzip -f {} \; 2>/dev/null || true
cp "$KLEB_EXTRACT/ncbi_dataset/data"/*/*_genomic.fna "$KAG"/ 2>/dev/null || true

# 2) Downloading Card Database
[ -f "$CARD_TAR" ] || curl -o "$CARD_TAR" https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2

# Extract card and copy files into fixed paths.
tar -xvf "$CARD_TAR" -C "$CD" >/dev/null
cp "$(find "$CD" -name nucleotide_fasta_protein_homolog_model.fasta | head -n 1)" "$CARD_FASTA"
cp "$(find "$CD" -name aro_index.tsv | head -n 1)" "$ARO_INDEX"
echo "done"
# 3) Creating Blast Database
export PATH="$PATH:$HOME/Bioinfo_UI/Blast_assessment/data/ncbi-blast-2.17.0+/bin"
makeblastdb -in "$CARD_FASTA" -title scripted_card_nucl_db -dbtype nucl -out "$DB"
echo "done creating database"

# 4) Extracting relevant columns inside aro_index tsv file and adding a new header
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "ARO_Index","Variant","Family","Resistance_Class"; next} {print $1,$6,$9,$10}' "$ARO_INDEX" > analysis/scripted_aro_reference.tsv
echo "done extraction"

# 5) Running Blast for the 5 assembled genomes, Filtering %identity, Splitting the ARO_Index and Gene column, Then Mapping the Ref tsv + Classification

for f in "$KAG"/*_genomic.fna; do
  base=$(basename "$f" _genomic.fna)

blastn -query "$f" -db "$DB" -evalue 1e-30 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
| awk -F'\t' 'BEGIN{OFS="\t"} $3>=90 {n=split($2,a,"|"); print $1,a[n-1],a[n],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > "$FINAL/${base}_final.tsv"
echo "done running blast"

# Mapping and AMR Classification
awk -F'\t' 'BEGIN{OFS="\t"} 
NR==FNR{m[$1]=$2"\t"$3"\t"$4; next} {aro=$2; add=(aro in m)?m[aro]:"NA\tNA\tNA"; print $0, add} 
' analysis/scripted_aro_reference.tsv "$FINAL/${base}_final.tsv" > "$CLASSIFIED/${base}_AMR_classified.tsv"

# Filtering the Relevant columns and sorting
OUTF="$FINAL_SORTED/${base}_AMR_sorted.tsv"
echo -e "ID\tGene_Variants\tGene_Families\tResistance_Class" > "$OUTF"
awk 'NR>1' "$CLASSIFIED/${base}_AMR_classified.tsv" | cut -f1,14,15,16 | sort -t$'\t' -k2,2 >> "$OUTF"
 echo "done running all steps"
done
 echo "Done. Final outputs are in: $FINAL_SORTED"
