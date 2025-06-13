#!/bin/bash
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -J sra_download_qc
#SBATCH --time=190:00:00
#SBATCH --output=sra_process_%j.log

# ===== Configurações =====
THREADS=20
MIN_LEN=200
MIN_COMP=50
MAX_CONT=10
MEMORY=800
ENV_NAME="metawrap-env"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Configurações específicas para SRA Toolkit
export TMPDIR="/temporario2/12689050/tmp_12"
mkdir -p "$TMPDIR"
export VDB_CONFIG="$HOME/.ncbi/user-settings.mkfg"

# ===== Lista de SRA IDs =====
SRA_IDS=(
    "SRR31832930"
    "SRR31832945"
)

# Diretorios base
BASE_DIR="/temporario2/12689050/PRJNA1191028/sra_12"
RESULTS_DIR="/temporario2/12689050/PRJNA1191028/resultados"

# ===== Diretórios =====
OUT_DIR="$BASE_DIR/raw_reads"
QC_DIR="$BASE_DIR/read_qc"
CLEAN="$BASE_DIR/clean_reads"
ASSEMBLY="$BASE_DIR/assembly"
BINNING="$BASE_DIR/initial_binning"
REFINEMENT="$BASE_DIR/bin_refinement"
REASSEMBLY="$BASE_DIR/bin_reassembly"
BLOBOLOGY="$BASE_DIR/blobology"
QUANT_BINS="$BASE_DIR/quant_bins"
CLASSIFICATION="$BASE_DIR/bin_classification"
ANNOTATION="$BASE_DIR/funct_annot"

# Criar todos os diretórios necessários
mkdir -p "$OUT_DIR" "$QC_DIR" "$CLEAN" "$ASSEMBLY" "$BINNING" "$REFINEMENT" \
         "$REASSEMBLY" "$BLOBOLOGY" "$QUANT_BINS" "$CLASSIFICATION" "$ANNOTATION" "$RESULTS_DIR"

# ===== Configurar PATH =====
if [[ ! "$PATH" =~ "/temporario2/12689050/metaWRAP/bin" ]]; then
    export PATH="$PATH:/temporario2/12689050/metaWRAP/bin"
fi

if [[ ! "$PATH" =~ "$HOME/micromamba/envs/$ENV_NAME/bin" ]]; then
    export PATH="$PATH:$HOME/micromamba/envs/$ENV_NAME/bin"
fi

# Função para verificar e mover arquivos/diretórios
safe_move() {
    local src=$1
    local dest_dir=$2
    
    if [ -e "$src" ]; then
        mkdir -p "$dest_dir"
        mv -v "$src" "$dest_dir/"
    else
        echo "AVISO: Item não encontrado para mover: $src" >&2
    fi
}

# Função para processar cada SRA ID
process_sra() {
    local SRA_ID=$1
    
    echo "========================================"
    echo "Processando SRA ID: $SRA_ID"
    echo "========================================"

    # ===== 1. Download SRA =====
    echo "Baixando $SRA_ID..."
    MAX_RETRIES=3
    RETRY_DELAY=60
    attempt=0
    
    while [ $attempt -lt $MAX_RETRIES ]; do
        if prefetch "$SRA_ID" --output-directory "$OUT_DIR"; then
            break
        else
            attempt=$((attempt+1))
            echo "Tentativa $attempt falhou. Aguardando $RETRY_DELAY segundos antes de tentar novamente..."
            sleep $RETRY_DELAY
        fi
    done
    
    if [ $attempt -eq $MAX_RETRIES ]; then
        echo "Falha no download do SRA $SRA_ID após $MAX_RETRIES tentativas"
        return 1
    fi

    # ===== 2. Conversão para FASTQ =====
    echo "Convertendo SRA para FASTQ..."
    
    # Verifica se o download foi bem sucedido (formato VDB moderno)
    if [ -d "$OUT_DIR/$SRA_ID" ]; then
        # Usando fasterq-dump com o formato VDB
        if ! fasterq-dump "$OUT_DIR/$SRA_ID" \
            --bufsize 2 \
            --curcache 20 \
            --mem 200 \
	    --outdir "$OUT_DIR/$SRA_ID" \
            --threads 15 \
            --temp "$TMPDIR" \
            --skip-technical \
            --progress; then
            
            # Tentativa alternativa com parallel-fastq-dump
            echo "Tentando parallel-fastq-dump como alternativa..."
            if ! parallel-fastq-dump --split-files --threads "$((THREADS/2))" \
                --outdir "$OUT_DIR/$SRA_ID" \
                --sra-id "$OUT_DIR/$SRA_ID"; then
                echo "Falha na conversão para FASTQ do $SRA_ID"
                return 1
            fi
        fi
    else
        echo "Erro: Não foi possível encontrar o arquivo/diretório baixado para $SRA_ID"
        return 1
    fi
    
     
    # ===== 3. Controle de Qualidade =====
    echo "Executando metaWRAP read_qc..."
    if ! micromamba run -n "$ENV_NAME" metawrap read_qc \
        -1 "$OUT_DIR/$SRA_ID/${SRA_ID}_1.fastq" \
        -2 "$OUT_DIR/$SRA_ID/${SRA_ID}_2.fastq" \
        -t "$THREADS" \
        -o "$QC_DIR/$SRA_ID"; then
        echo "Falha no read_qc do $SRA_ID"
        return 1
    fi

    # ===== 4. Mover reads limpas =====
    echo "Transferindo reads limpas..."
    mkdir -p "$CLEAN/$SRA_ID"
    mv "$QC_DIR/$SRA_ID/final_pure_reads_1.fastq" "$CLEAN/$SRA_ID/${SRA_ID}_1.fastq" || return 1
    mv "$QC_DIR/$SRA_ID/final_pure_reads_2.fastq" "$CLEAN/$SRA_ID/${SRA_ID}_2.fastq" || return 1

    # 4.1 Salvar relatorios de QC
    mkdir -p "$RESULTS_DIR/$SRA_ID/read_qc"
    safe_move "$QC_DIR/$SRA_ID/post-QC_report" "$RESULTS_DIR/$SRA_ID/read_qc"
    safe_move "$QC_DIR/$SRA_ID/pre-QC_report" "$RESULTS_DIR/$SRA_ID/read_qc"

    # ===== 5. Assembly =====
    echo "Iniciando assembly..."
    mkdir -p "$ASSEMBLY/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap assembly \
        -1 "$CLEAN/$SRA_ID/${SRA_ID}_1.fastq" \
        -2 "$CLEAN/$SRA_ID/${SRA_ID}_2.fastq" \
        -m "$MIN_LEN" \
        -t "$THREADS" \
        -o "$ASSEMBLY/$SRA_ID"; then
        echo "Falha no assembly do $SRA_ID"
        return 1
    fi

    # 5.1 Salvando resultados do assembly
    mkdir -p "$RESULTS_DIR/$SRA_ID/assembly"
    safe_move "$ASSEMBLY/$SRA_ID/QUAST_out/basic_stats" "$RESULTS_DIR/$SRA_ID/assembly"
    safe_move "$ASSEMBLY/$SRA_ID/QUAST_out/icarus_viewers" "$RESULTS_DIR/$SRA_ID/assembly"
    safe_move "$ASSEMBLY/$SRA_ID/QUAST_out/report.html" "$RESULTS_DIR/$SRA_ID/assembly"
    safe_move "$ASSEMBLY/$SRA_ID/assembly_report.html" "$RESULTS_DIR/$SRA_ID/assembly"

    # ===== 6. Binning =====
    echo "Iniciando binning..."
    mkdir -p "$BINNING/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap binning \
        -o "$BINNING/$SRA_ID" \
        -t "$THREADS" \
        -a "$ASSEMBLY/$SRA_ID/final_assembly.fasta" \
        --metabat2 --maxbin2 --concoct \
        "$CLEAN/$SRA_ID/${SRA_ID}"*.fastq; then
        echo "Falha no binning do $SRA_ID"
        return 1
    fi

    # ===== 7. Refinamento =====
    echo "Refinando bins..."
    mkdir -p "$REFINEMENT/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap bin_refinement \
        -o "$REFINEMENT/$SRA_ID" \
        -t "$THREADS" \
        -A "$BINNING/$SRA_ID/metabat2_bins/" \
        -B "$BINNING/$SRA_ID/maxbin2_bins/" \
        -C "$BINNING/$SRA_ID/concoct_bins/" \
        -c "$MIN_COMP" \
        -x "$MAX_CONT"; then
        echo "Falha no refinamento de bins do $SRA_ID"
        return 1
    fi

    # 7.1 Salvar os resultados do refinamento
    mkdir -p "$RESULTS_DIR/$SRA_ID/bin_refinement"
    safe_move "$REFINEMENT/$SRA_ID/figures" "$RESULTS_DIR/$SRA_ID/bin_refinement"
    safe_move "$REFINEMENT/$SRA_ID"/*.stats "$RESULTS_DIR/$SRA_ID/bin_refinement"
    cp -r "$REFINEMENT/$SRA_ID/concoct_bins" "$RESULTS_DIR/$SRA_ID/bin_refinement"
    cp -r "$REFINEMENT/$SRA_ID/maxbin2_bins" "$RESULTS_DIR/$SRA_ID/bin_refinement"
    cp -r "$REFINEMENT/$SRA_ID/metabat_bins" "$RESULTS_DIR/$SRA_ID/bin_refinement"
    cp -r "$REFINEMENT/$SRA_ID/work_files" "$RESULTS_DIR/$SRA_ID/bin_refinement"
    cp -r "$REFINEMENT/$SRA_ID/metawrap_${MIN_COMP}_${MAX_CONT}_bins" "$RESULTS_DIR/$SRA_ID/bin_refinement"

    # ===== 8. Blobology =====
    echo "Executando blobology..."
    mkdir -p "$BLOBOLOGY/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap blobology \
        -a "$ASSEMBLY/$SRA_ID/final_assembly.fasta" \
        -t "$THREADS" \
        -o "$BLOBOLOGY/$SRA_ID" \
        --bins "$REFINEMENT/$SRA_ID/metawrap_${MIN_COMP}_${MAX_CONT}_bins" \
        "$CLEAN/$SRA_ID/${SRA_ID}"*.fastq; then
        echo "Falha no blobology do $SRA_ID"
        return 1
    fi

    # 8.1 Salvar os resultados do blobology
    mkdir -p "$RESULTS_DIR/$SRA_ID/blobology"
    safe_move "$BLOBOLOGY/$SRA_ID/blobplot_figures" "$RESULTS_DIR/$SRA_ID/blobology"
    safe_move "$BLOBOLOGY/$SRA_ID/blobplot_figures_only_binned_contigs" "$RESULTS_DIR/$SRA_ID/blobology"

    # ===== 9. Quantificação =====
    echo "Quantificando bins..."
    mkdir -p "$QUANT_BINS/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap quant_bins \
        -b "$REFINEMENT/$SRA_ID/metawrap_${MIN_COMP}_${MAX_CONT}_bins" \
        -o "$QUANT_BINS/$SRA_ID" \
        -a "$ASSEMBLY/$SRA_ID/final_assembly.fasta" \
        "$CLEAN/$SRA_ID/${SRA_ID}"*.fastq; then
        echo "Falha na quantificação de bins do $SRA_ID"
        return 1
    fi

    # 9.1 Salvar dados do quant_bins
    mkdir -p "$RESULTS_DIR/$SRA_ID/quant_bins"
    cp -r "$QUANT_BINS/$SRA_ID" "$RESULTS_DIR/$SRA_ID/quant_bins"

    # ===== 10. Reassembly =====
    echo "Reassemblando bins..."
    mkdir -p "$REASSEMBLY/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap reassemble_bins \
        -o "$REASSEMBLY/$SRA_ID" \
        -1 "$CLEAN/$SRA_ID/${SRA_ID}_1.fastq" \
        -2 "$CLEAN/$SRA_ID/${SRA_ID}_2.fastq" \
        -t "$THREADS" \
        -m "$MEMORY" \
        -c "$MIN_COMP" \
        -x "$MAX_CONT" \
        -b "$REFINEMENT/$SRA_ID/metawrap_${MIN_COMP}_${MAX_CONT}_bins"; then
        echo "Falha no reassembly de bins do $SRA_ID"
        return 1
    fi

    # 10.1 Obtendo os dados do reassembly_bins
    mkdir -p "$RESULTS_DIR/$SRA_ID/reassembled_bins"
    for item in "$REASSEMBLY/$SRA_ID"/*.stats "$REASSEMBLY/$SRA_ID"/*.png "$REASSEMBLY/$SRA_ID"/original_bins; do
    	cp -r "$item" "$RESULTS_DIR/$SRA_ID/reassembled_bins"
    done


    # ===== 11. Classificação =====
    echo "Classificando bins..."
    mkdir -p "$CLASSIFICATION/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap classify_bins \
        -b "$REASSEMBLY/$SRA_ID/reassembled_bins" \
        -o "$CLASSIFICATION/$SRA_ID" \
        -t "$THREADS"; then
        echo "Falha na classificação de bins do $SRA_ID"
        return 1
    fi

    # 11.1 Obtendo os dados do classify_bins
    mkdir -p "$RESULTS_DIR/$SRA_ID/classify_bins"
    cp -r "$CLASSIFICATION/$SRA_ID" "$RESULTS_DIR/$SRA_ID/classify_bins"

    # ===== 12. Renomeação de contigs =====
    echo "Renomeando contigs..."
    mkdir -p "$REASSEMBLY/$SRA_ID/manual_renamed_bins"
    for bin in "$REASSEMBLY/$SRA_ID/reassembled_bins"/*.fa; do
        bin_name=$(basename "$bin" .fa)
        echo "Renomeando $bin_name"
        awk '/^>/ {print ">contig_" ++i; next} {print}' "$bin" > "$REASSEMBLY/$SRA_ID/manual_renamed_bins/${bin_name}.manual_renamed.fa"
    done
    echo "Renomeação completa"

    # ===== 13. Anotação =====
    echo "Anotando bins..."
    mkdir -p "$ANNOTATION/$SRA_ID"
    if ! micromamba run -n "$ENV_NAME" metawrap annotate_bins \
        -o "$ANNOTATION/$SRA_ID" \
        -t "$THREADS" \
        -b "$REASSEMBLY/$SRA_ID/manual_renamed_bins"; then
        echo "Falha na anotação de bins do $SRA_ID"
        return 1
    fi

    # 13.1 Obtendo os dados do funct_annot
    mkdir -p "$RESULTS_DIR/$SRA_ID/funct_annot"
    cp -r "$ANNOTATION/$SRA_ID" "$RESULTS_DIR/$SRA_ID/funct_annot"

    echo "Pipeline concluído com sucesso para $SRA_ID!"
    return 0
}

# Processar cada SRA ID
for SRA_ID in "${SRA_IDS[@]}"; do
    # Processa o SRA_ID atual
    if ! process_sra "$SRA_ID"; then
        echo "Erro no processamento de $SRA_ID, continuando para o próximo..."
        continue  # Pula a limpeza se falhar
    fi

    # Limpeza APENAS do SRA_ID atual (após sucesso)
    echo "Removendo dados temporários de $SRA_ID..."
    rm -rf "$BASE_DIR/assembly/$SRA_ID"
    rm -rf "$BASE_DIR/raw_reads/$SRA_ID"
    rm -rf "$BASE_DIR/read_qc/$SRA_ID"
    rm -rf "$BASE_DIR/clean_reads/$SRA_ID"
    rm -rf "$BASE_DIR/quant_bins/$SRA_ID"
    rm -rf "$BASE_DIR/initial_binning/$SRA_ID"
    rm -rf "$BASE_DIR/bin_classification/$SRA_ID"
    rm -rf "$BASE_DIR/bin_reassembly/$SRA_ID"
    rm -rf "$BASE_DIR/bin_refinement/$SRA_ID"
    rm -rf "$BLOBOLOGY/$SRA_ID"
    rm -rf "$ANNOTATION/$SRA_ID"
    rm -rf "$TMPDIR"
done

echo "Processamento de todos os SRA IDs concluído!"
