
#
# Run examples from https://github.com/michauhl/RBPBench/README.md
#
# rbpbench should be installed / the env activated.
# test folder needs to contain:
# eclip_clipper_idr
#
# bash run_examples.sh 2>&1 /path/to/hg38.fa /path/to/Homo_sapiens.GRCh38.xxx.gtf.gz /path/to/hg38.phastCons100way.bw /path/to/hg38.phyloP100way.bw | tee run_examples.log
#



genome_fa_path="$1"
gtf_path="$2"
pc_bbw_path="$3"
pp_bbw_path="$3"


# Check if files exist.
if [ ! -f "$genome_fa_path" ]; then
    echo "Error: File '$genome_fa_path' does not exist."
    exit 1
fi
if [ ! -f "$gtf_path" ]; then
    echo "Error: File '$gtf_path' does not exist."
    exit 1
fi
if [ ! -f "$pc_bbw_path" ]; then
    echo "Error: File '$pc_bbw_path' does not exist."
    exit 1
fi
if [ ! -f "$pp_bbw_path" ]; then
    echo "Error: File '$pp_bbw_path' does not exist."
    exit 1
fi

# Search with selected RBPs.
rbpbench search --in eclip_clipper_idr/PUM2_K562_IDR_peaks.bed --genome $genome_fa_path --gtf $gtf_path --out test_search_pum2_ex1_out --rbps PUM2 PUM1 RBFOX2 --ext 10 --regex AATAAA

# Search with all RBPs.
rbpbench search --in eclip_clipper_idr/SLBP_K562_IDR_peaks.bed --genome $genome_fa_path --gtf $gtf_path --out test_search_slbp_tr_rsd_tep_out --functions TR RSD TEP --rbps ALL --ext 20 --goa

# With greatest hits.
rbpbench search --in eclip_clipper_idr/SLBP_K562_IDR_peaks.bed --genome $genome_fa_path --gtf $gtf_path --out test_search_slbp_tr_rsd_tep_gh_out --functions TR RSD TEP  --rbps ALL --ext 20 --set-rbp-id SLBP --greatest-hits

# Select individual motifs.
rbpbench search --in eclip_clipper_idr/PUM2_K562_IDR_peaks.bed --genome $genome_fa_path --gtf $gtf_path --out test_search_pum2_1_2_out --rbps ALL --ext 10 --motifs PUM2_1 PUM2_2

# Motif search with multiple input datasets.
rbpbench batch --bed eclip_clipper_idr.k562.batch_in.txt --genome $genome_fa_path --gtf $gtf_path --ext 10 --out batch_K562_eclip_clipper_idr_out --regex AATAAA

# Comparisons between search results
rbpbench batch --bed batch_compare_test.batch_in.txt --genome $genome_fa_path --ext 10 --out test_batch_out
rbpbench compare --in test_batch_out --out test_compare_out

# Single motif enrichment and co-occurrences: ENMO mode
rbpbench enmo --in eclip_clipper_idr/PUM2_K562_IDR_peaks.bed --genome $genome_fa_path --gtf $gtf_path --out test_enmo_pum2_out --rbps ALL --ext 10 --min-motif-dist 10 --motif-sim-thr 2

# Single motif enrichment and co-occurrences: NEMO mode
gtf_extract_mpt_region_bed.py --gtf $gtf_path --out mrna_regions_hg38_out --mrna-only
bed_print_last_n_pos.py --in mrna_regions_hg38_out/mpt_regions.bed --ext 1 > mrna_region_end_pos.bed
rbpbench nemo --in mrna_region_end_pos.bed --genome $genome_fa_path --gtf $gtf_path --out test_nemo_mrna_ends_out --rbps ALL --ext 40 --min-motif-dist 10 --motif-sim-thr 2 --allow-overlaps --functions TEP RSD TR

# Motif preferences in spliced full transcripts: CSTF2_1 (AATAAA)
rbpbench searchlongrna --genome $genome_fa_path --gtf $gtf_path --out test_searchlongrna_mrna_pas_out --rbps ALL --motifs CSTF2_1 --mrna-only
# Motif preferences in spliced full transcripts: DDX3X
rbpbench searchlongrna --genome $genome_fa_path --gtf $gtf_path --out test_searchlongrna_mrna_ddx3x_out --rbps DDX3X --mrna-only

# Motif preferences in long genomic regions.
rbpbench searchlong --in tr_list.txt --genome $genome_fa_path --gtf $gtf_path --rbps DDX3X --out searchlong_test_ddx3x_out

# Plot nucleotide distribution at genomic positions.
gtf_extract_tr_feat_bed.py --feat stop_codon --gtf $gtf_path --out stop_codons.Homo_sapiens.GRCh38.112.bed --uniq-reg
rbpbench dist --in stop_codons.Homo_sapiens.GRCh38.112.bed --genome $genome_fa_path --out test_dist_out --ext 5

# Comparing conservation scores between two sets of genomic sites.
bed_shift_regions.py --in mrna_region_end_pos.bed --num 50 > mrna_region_end_pos.50ds_shift.bed
rbpbench con --in mrna_region_end_pos.bed --ctrl-in mrna_region_end_pos.50ds_shift.bed --phylop $pp_bbw_path --phastcons $pc_bbw_path --out mrna_region_end_pos_con_sc_out

# Searching for sponge transcripts.
rbpbench sponge --regex 'TGTA[ACGT]ATA' --out test_sponge_search_gtf_out --genome $genome_fa_path --gtf $gtf_path --min-seq-len 1000

# Comparing motif hits between transcript isoforms.
rbpbench isocomp --regex "TGTA[ACGT]ATA" --out test_isocomp_gtf_out --genome $genome_fa_path --gtf $gtf_path --select-mode 1


