This Repository contains the complete software and documentation to execute the Long-Read-Proteogenomics Workflow.

## Data Flow Summary

![1.drawio](D:\AAA研究生课题项目\1.drawio.png)



## Digital Object Identifiers

| Sequence Read Archive (SRA) Project Reference                | Description                                          |
| ------------------------------------------------------------ | ---------------------------------------------------- |
| [PRJNA783347](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA783347) | Long-Read RNA Sequencing Project for Jurkat Samples  |
| [PRJNA193719](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA193719) | Short-Read RNA Sequencing Project for Jurkat Samples |



## Test Module and Results

- Module

  # Transcriptomic Processes

  ## PacBio Iso‑Seq data analysis

  ### SMARTLink CCS 

  - 原始reads处理为CCS reads

  ```shell
  ccs 0data/*.subreads.bam 1SMARTLink_ccs/movie.ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75 --num-threads 12
  ```

  - find and remove adapters/barcodes

  ```shell
  # barcode序列文件
  # 针对RNA建库: 
  # 5'primer - 后缀为_5p
  # 5'primer - 后缀为_3p
  lima --isoseq --dump-clips --peek-guess -j 24 1SMARTLink_ccs/movie.ccs.bam 0data/barcode.fasta 1SMARTLink_ccs/movie.fl_adpt.bam
  ```

  ### IsoSeq3

  ```shell
  # ensure that only qv10 reads from ccs are input
  bamtools filter -tag 'rq':'>=0.90' -in 1SMARTLink_ccs/movie.ccs.bam -out 2IsoSeq3/filtered.movie.ccs.bam 
  
  # create an index for the ccs bam
  pbindex 2IsoSeq3/filtered.movie.ccs.bam 
  
  # find and remove adapters/barcodes
  lima --isoseq --dump-clips --peek-guess -j 24 2IsoSeq3/filtered.movie.ccs.bam 0data/barcode.fasta 2IsoSeq3/movie.demult.bam
  
  # filter for non-concatamer, polya containing reads
  isoseq3 refine --require-polya 2IsoSeq3/movie.demult.primer_5p--primer_3p.bam 0data/barcode.fasta 2IsoSeq3/movie.flnc.bam
  
  # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
  isoseq3 cluster 2IsoSeq3/movie.flnc.bam 2IsoSeq3/movie.clustered.bam --verbose --use-qvs
  
  # align reads to the genome, takes few minutes (40 core machine)
  pbmm2 align 0data/reference.fasta 2IsoSeq3/movie.clustered.hq.bam 2IsoSeq3/movie.aligned.bam --preset ISOSEQ --sort -j 24 --log-level INFO
  
  # collapse redundant reads
  isoseq3 collapse 2IsoSeq3/movie.aligned.bam  2IsoSeq3/movie.collapsed.gff
  ```

  

  ## Transcript isoform classification and filtering

  ### SQANTI3

  ```
  input:
  GENCODE version 35 annotations(GTF file)
  ref genome(GRCh38)
  ```

  ```shell
  sqanti3_qc.py \
  2IsoSeq3/movie.collapsed.gff \
  0data/gencode.gtf \
  0data/reference.fasta \
  --short_reads 0data/short_reads.fofn \
  --skipORF \
  -o sqanti3 -d 3SQANTI3 \
  --fl_count 2IsoSeq3/movie.collapsed.abundance.txt \
  --cpus 24 --report skip
  
  # filter
  # generate_reference_tables
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/generate_reference_tables/src/prepare_reference_tables.py \
    --gtf 0data/gencode.gtf \
    --fa 0data/gencode_transcript.fasta \
    --ensg_gene 3SQANTI3/ensg_gene.tsv \
    --enst_isoname 3SQANTI3/enst_isoname.tsv \
    --gene_ensp 3SQANTI3/gene_ensp.tsv \
    --gene_isoname 3SQANTI3/gene_isoname.tsv \
    --isoname_lens 3SQANTI3/isoname_lens.tsv \
    --gene_lens 3SQANTI3/gene_lens.tsv \
    --protein_coding_genes 3SQANTI3/protein_coding_genes.txt
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/filter_sqanti/src/filter_sqanti.py \
  --sqanti_classification 3SQANTI3/sqanti3_classification.txt \
  --sqanti_corrected_fasta 3SQANTI3/sqanti3_corrected.fasta \
  --sqanti_corrected_gtf 3SQANTI3/sqanti3_corrected.gtf \
  --protein_coding_genes 3SQANTI3/protein_coding_genes.txt \
  --ensg_gene 3SQANTI3/ensg_gene.tsv \
  --filter_protein_coding yes \
  --filter_intra_polyA yes \
  --filter_template_switching yes \
  --percent_A_downstream_threshold 95 \
  --structural_categories_level strict \
  --minimum_illumina_coverage 3
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/filter_sqanti/src/collapse_isoforms.py \
  --name movie \
  --sqanti_gtf 3SQANTI3/sqanti3_corrected.gtf \
  --sqanti_fasta 3SQANTI3/sqanti3_corrected.fasta
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/filter_sqanti/src/collapse_classification.py \
  --name movie \
  --collapsed_fasta movie_corrected.5degfilter.fasta \
  --classification filtered_sqanti3_classification.tsv 
  
  # remove 5′ transcript degradation products
  filter_away_subset.py
  ```

  

  ## Generation of a full‑length protein isoform database from long‑read RNA‑seq ORF prediction

  ### CPAT

  ```shell
    cpat.py \
    -x 0data/Human_Hexamer.tsv \
    -d 0data/Human_logitModel.RData\
    -g SQANTI3/movie_corrected.5degfilter.fasta \
    --min-orf=50 \
    --top-orf=50 \
    -o 4CPAT/4CPAT
  ```

  

  # Proteogenomic Processes

  ### ORF calling

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/transcriptome_summary/src/transcriptome_summary.py \
  --sq_out ../3SQANTI3/movie_classification.5degfilter.tsv \
  --tpm /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/data/jurkat_gene_kallisto.tsv \
  --ribo /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/data/kallist_table_rdeplete_jurkat.tsv \
  --ensg_to_gene ../3SQANTI3/ensg_gene.tsv \
  --enst_to_isoname ../3SQANTI3/enst_isoname.tsv \
  --len_stats ../3SQANTI3/gene_lens.tsv
  
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/orf_calling/src/orf_calling.py \
  --orf_coord ../4CPAT/movie.ORF_prob.tsv \
  --orf_fasta ../4CPAT/movie.ORF_seqs.fa \
  --gencode ../0data/gencode.gtf \
  --sample_gtf ../3SQANTI3/movie_corrected.5degfilter.gff \
  --pb_gene ./pb_gene.tsv \
  --classification ../3SQANTI3/movie_classification.5degfilter.tsv \
  --sample_fasta ../3SQANTI3/movie_corrected.5degfilter.fasta \
  --num_cores 24 \
  --output orfcalling_best_orf.tsv
  ```

  ### DB Generation

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/refine_orf_database/src/refine_orf_database.py \
  --name movie \
  --orfs ../5orf_calling/orfcalling_best_orf.tsv \
  --pb_fasta ../3SQANTI3/movie_corrected.5degfilter.fasta \
  --coding_score_cutoff 0.0
  ```

  ### DB Annotation

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/make_pacbo_cds_gtf/src/make_pacbio_cds_gtf.py \
  --sample_gtf ../3SQANTI3/movie_corrected.5degfilter.gff \
  --agg_orfs ../6DB_Generation/movie_orf_refined.tsv \
  --refined_orfs ../5orf_calling/orfcalling_best_orf.tsv \
  --pb_gene ../5orf_calling/pb_gene.tsv \
  --output_cds pacbio_database.gtf
  
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/rename_cds_to_exon/src/rename_cds_to_exon.py \
  --sample_gtf pacbio_database.gtf \
  --sample_name movie \
  --reference_gtf ../0data/gencode.gtf  \
  --reference_name gencode \
  --num_cores 24
  ```

  ### SQANTI Protein

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/sqanti_protein/src/sqanti3_protein.py \
  ../7DB_Annotation/movie.transcript_exons_only.gtf \
  ../7DB_Annotation/movie.cds_renamed_exon.gtf  \
  ../5orf_calling/orfcalling_best_orf.tsv \
  ../7DB_Annotation/gencode.transcript_exons_only.gtf \
  ../7DB_Annotation/gencode.cds_renamed_exon.gtf \
  -d ./ \
  -p movie
  
  ```

  ### Protein Classification

  ```shell
  # 5' UTR Status
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/5p_utr_status/src/1_get_gc_exon_and_5utr_info.py \
  --gencode_gtf ../0data/gencode.gtf \
  --odir ./
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/5p_utr_status/src/2_classify_5utr_status.py \
  --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
  --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
  --sample_cds_gtf ../7DB_Annotation/pacbio_database.gtf \
  --odir ./
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/5p_utr_status/src/3_merge_5utr_info_to_pclass_table.py \
  --name movie \
  --utr_info pb_5utr_categories.tsv \
  --sqanti_protein_classification ../8SQANTI_Protein/movie.sqanti_protein_classification.tsv  \
  --odir ./
  
  
  # Protein Classification
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_classification/src/protein_classification_add_meta.py \
  --protein_classification movie.sqanti_protein_classification_w_5utr_info.tsv \
  --best_orf ../5orf_calling/orfcalling_best_orf.tsv  \
  --refined_meta ../6DB_Generation/movie_orf_refined.tsv \
  --ensg_gene ../3SQANTI3/ensg_gene.tsv \
  --name movie \
  --dest_dir ./
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_classification/src/protein_classification.py \
  --sqanti_protein movie.protein_classification_w_meta.tsv \
  --name movie_unfiltered \
  --dest_dir ./
  
  # Protein Gene Rename
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_gene_rename/src/protein_gene_rename.py \
  --sample_gtf ../7DB_Annotation/pacbio_database.gtf \
  --sample_protein_fasta ../6DB_Generation/movie_orf_refined.fasta \
  --sample_refined_info ../6DB_Generation/movie_orf_refined.tsv \
  --pb_protein_genes movie_genes.tsv \
  --name movie
  
  # Protein Filtering
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_filter/src/protein_filter.py \
  --protein_classification movie_unfiltered.protein_classification.tsv \
  --gencode_gtf ../0data/gencode.gtf \
  --protein_fasta  movie.protein_refined.fasta \
  --sample_cds_gtf movie_with_cds_refined.gtf \
  --min_junctions_after_stop_codon 2 \
  --name movie
  ```

  ### Protein Hybrid Database

  ```shell
  # GENCODE Database
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/make_gencode_database/src/make_gencode_database.py \
  --gencode_fasta ../0data/gencode_translations.fasta \
  --protein_coding_genes ../3SQANTI3/protein_coding_genes.txt \
  --output_fasta gencode_protein.fasta \
  --output_cluster gencode_isoname_clusters.tsv
  
  # Protein Hybrid Database
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/make_hybrid_database/src/make_hybrid_database.py \
  --protein_classification ../9Protein_Classification/movie.classification_filtered.tsv \
  --gene_lens ../3SQANTI3/gene_lens.tsv \
  --pb_fasta ../9Protein_Classification/movie.filtered_protein.fasta \
  --gc_fasta gencode_protein.fasta \
  --refined_info ../9Protein_Classification/movie_orf_refined_gene_update.tsv \
  --pb_cds_gtf ../9Protein_Classification/movie_with_cds_filtered.gtf \
  --name movie \
  --lower_kb 1 \
  --upper_kb 4 \
  --lower_cpm 3
  ```

  

  # Proteomics

  ## Metamorpheus

  ### Mass Spec File Conversion

  ```shell
  sudo docker run -it --rm -v /mnt/usb1/chenhu/work/0data/PeptideAtlas:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ./*.raw --filter "peakPicking true 1-"
  
  sudo docker run -it --rm -v /mnt/usb1/chenhu/work/0data/PeptideAtlas:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert *.mzML --merge -o combined
  ```

  ### Metamorpheus GENCODE

  ```shell
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -d ../10Protein_Hybrid_Database/gencode_protein.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv
  ```

  ### Metamorpheus UniProt

  ```shell
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -d ../0data/uniprot_reviewed_canonical_and_isoform.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv
  ```

  ### MetaMorpheus with Sample Specific Database - Refined

  ```shell
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -d ../9Protein_Classification/movie.protein_refined.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.movie.refined.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.refined.tsv
  ```

  ### MetaMorpheus with Sample Specific Database - Filtered

  ```shell
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -d ../9Protein_Classification/movie.filtered_protein.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.movie.filtered.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.filtered.tsv
  
  ```

  ### MetaMorpheus with Sample Specific Database - Hybrid

  ```shell
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/metamorpheus/CMD.dll -d ../10Protein_Hybrid_Database/movie_hybrid.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.movie.hybrid.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.hybrid.tsv
  ```

  ### MetaMorpheus using "Rescue and Resolve" Algorithm

  ```shell
  dotnet /mnt/usb1/chenhu/program/MetaMorpheus-LongReadProteogenomics/CMD/bin/Debug/netcoreapp3.1/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet /mnt/usb1/chenhu/program/MetaMorpheus-LongReadProteogenomics/CMD/bin/Debug/netcoreapp3.1/CMD.dll -d ../10Protein_Hybrid_Database/movie_hybrid.fasta settings/Contaminants/MetaMorpheusContaminants.xml -s ../0data/PeptideAtlas/combined/120426_Jurkat_highLC_Frac.mzML -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results --orf ../10Protein_Hybrid_Database/movie_refined_high_confidence.tsv --cpm 25
  
  mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.movie.rescue_resolve.psmtsv
  mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.rescue_resolve.tsv
  ```

  ## Analysis

  ### Peptide Analysis

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/peptide_analysis/src/peptide_analysis.py \
  -gmap ../3SQANTI3/gene_isoname.tsv \
  --gencode_peptides ../11Metamorpheus_GENCODE/search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv \
  --pb_refined_fasta ../9Protein_Classification/movie.protein_refined.fasta \
  --pb_filtered_fasta ../9Protein_Classification/movie.filtered_protein.fasta \
  --pb_hybrid_fasta ../10Protein_Hybrid_Database/movie_hybrid.fasta \
  --pb_gene ../5orf_calling/pb_gene.tsv \
  -odir ./
  ```

  ### Reference Track Visualization

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/gencode_filter_protein_coding.py \
  --reference_gtf ../../0data/gencode.gtf
  
  gtfToGenePred gencode.filtered.gtf gencode.filtered.genePred
  genePredToBed gencode.filtered.genePred gencode.filtered.bed12
  
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/gencode_add_rgb_to_bed.py \
  --gencode_bed gencode.filtered.bed12 \
  --rgb 0,0,140 \
  --version V35
  ```

  ### Peptide Track Visualization

  ```shell
  # Convert GTF to Bed
  # Refined
  gtfToGenePred ../../9Protein_Classification/movie_with_cds_refined.gtf movie_refined_cds.genePred
  genePredToBed movie_refined_cds.genePred movie_refined_cds.bed12
  # Filtered
  gtfToGenePred ../../9Protein_Classification/movie_with_cds_filtered.gtf movie_filtered_cds.genePred
  genePredToBed movie_filtered_cds.genePred movie_filtered_cds.bed12
  # Hybrid
  gtfToGenePred ../../10Protein_Hybrid_Database/movie_cds_high_confidence.gtf movie_hybrid_cds.genePred
  genePredToBed movie_hybrid_cds.genePred movie_hybrid_cds.bed12
  
  # Add RGB colors
  # Refined
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/track_add_rgb_colors_to_bed.py \
  --name movie_refined \
  --bed_file movie_refined_cds.bed12
  # Filtered
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/track_add_rgb_colors_to_bed.py \
  --name movie_filtered \
  --bed_file movie_filtered_cds.bed12
  # Hybrid
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/track_add_rgb_colors_to_bed.py \
  --name movie_hybrid \
  --bed_file movie_hybrid_cds.bed12
  ```

  ### Multiregion BED generation

  ```shell
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_region_bed_for_ucsc.py \
  --name movie_refined \
  --sample_gtf ../9Protein_Classification/movie_with_cds_refined.gtf \
  --reference_gtf ../18Track_Visualization/reference/gencode.filtered.gtf
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_region_bed_for_ucsc.py \
  --name movie_filtered \
  --sample_gtf ../9Protein_Classification/movie_with_cds_filtered.gtf \
  --reference_gtf ../18Track_Visualization/reference/gencode.filtered.gtf
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_region_bed_for_ucsc.py \
  --name movie_high_confidence \
  --sample_gtf ../10Protein_Hybrid_Database/movie_cds_high_confidence.gtf \
  --reference_gtf ../18Track_Visualization/reference/gencode.filtered.gtf
  ```

  ### Peptide Track Visualization

  ```shell
  # Make peptide gtf files
  # Refined
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_peptide_gtf_file.py \
  --name movie_refined \
  --sample_gtf ../../9Protein_Classification/movie_with_cds_refined.gtf \
  --reference_gtf ../../18Track_Visualization/reference/gencode.filtered.gtf \
  --peptides ../../13Metamorpheus_Sample/search_results/Task1SearchTask/AllPeptides.movie.refined.psmtsv \
  --pb_gene ../../5orf_calling/pb_gene.tsv \
  --gene_isoname ../../3SQANTI3/gene_isoname.tsv \
  --refined_fasta ../../9Protein_Classification/movie.protein_refined.fasta
  # Filtered
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_peptide_gtf_file.py \
  --name movie_filtered \
  --sample_gtf ../../9Protein_Classification/movie_with_cds_filtered.gtf \
  --reference_gtf ../../18Track_Visualization/reference/gencode.filtered.gtf \
  --peptides ../../14Metamorpheus_Sample_filtered/search_results/Task1SearchTask/AllPeptides.movie.filtered.psmtsv \
  --pb_gene ../../5orf_calling/pb_gene.tsv \
  --gene_isoname ../../3SQANTI3/gene_isoname.tsv \
  --refined_fasta ../../9Protein_Classification/movie.filtered_protein.fasta
  # Hybrid
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/make_peptide_gtf_file.py \
  --name movie_hybrid \
  --sample_gtf ../../10Protein_Hybrid_Database/movie_cds_high_confidence.gtf \
  --reference_gtf ../../18Track_Visualization/reference/gencode.filtered.gtf \
  --peptides ../../15Metamorpheus_Sample_Hybrid/search_results/Task1SearchTask/AllPeptides.movie.hybrid.psmtsv \
  --pb_gene ../../5orf_calling/pb_gene.tsv \
  --gene_isoname ../../3SQANTI3/gene_isoname.tsv \
  --refined_fasta ../../10Protein_Hybrid_Database/movie_hybrid.fasta
  
  # Convert GTF to bed and add RGB
  # Refined
  gtfToGenePred movie_refined_peptides.gtf movie_refined_peptides.genePred
  genePredToBed movie_refined_peptides.genePred movie_refined_peptides.bed12
  # add rgb to colorize specific peptides 
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/finalize_peptide_bed.py \
  --bed movie_refined_peptides.bed12 \
  --name movie_refined
  # Filtered
  gtfToGenePred movie_filtered_peptides.gtf movie_filtered_peptides.genePred
  genePredToBed movie_filtered_peptides.genePred movie_filtered_peptides.bed12
  # add rgb to colorize specific peptides 
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/finalize_peptide_bed.py \
  --bed movie_filtered_peptides.bed12 \
  --name movie_filtered
  # High Confidence
  gtfToGenePred movie_hybrid_peptides.gtf movie_hybrid_peptides.genePred
  genePredToBed movie_hybrid_peptides.genePred movie_hybrid_peptides.bed12
  # add rgb to colorize specific peptides 
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/visualization_track/src/finalize_peptide_bed.py \
  --bed movie_hybrid_peptides.bed12 \
  --name movie_hybrid
  ```

  ### Accession Mapping 

  ```shell
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/accession_mapping/src/accession_mapping.py \
  --gencode_fasta ../10Protein_Hybrid_Database/gencode_protein.fasta \
  --pacbio_fasta ../9Protein_Classification/movie.protein_refined.fasta  \
  --uniprot_fasta ../9Protein_Classification/movie.protein_refined.fasta 
  ```

  ### Protein Group Comparison

  ```shell
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_groups_compare/src/protein_groups_compare.py \
  --pg_fileOne ../11Metamorpheus_GENCODE/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv \
  --pg_fileTwo ../15Metamorpheus_Sample_Hybrid/search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.hybrid.tsv  \
  --mapping ../20Accession_Mapping/accession_map_gencode_uniprot_pacbio.tsv \
  --output ./
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_groups_compare/src/protein_groups_compare.py \
  --pg_fileOne ../11Metamorpheus_GENCODE/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv \
  --pg_fileTwo ../12Metamorpheus_UniProt/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv \
  --mapping ../20Accession_Mapping/accession_map_gencode_uniprot_pacbio.tsv \
  --output ./
  #
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/protein_groups_compare/src/protein_groups_compare.py \
  --pg_fileOne ../12Metamorpheus_UniProt/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv \
  --pg_fileTwo ../15Metamorpheus_Sample_Hybrid/search_results/Task1SearchTask/AllQuantifiedProteinGroups.movie.hybrid.tsv \
  --mapping ../20Accession_Mapping/accession_map_gencode_uniprot_pacbio.tsv \
  --output ./
  ```

  ### Novel Peptides

  ```shell
  # Refined
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
  --pacbio_peptides ../13Metamorpheus_Sample/search_results/Task1SearchTask/AllPeptides.movie.refined.psmtsv \
  --gencode_fasta ../10Protein_Hybrid_Database/gencode_protein.fasta \
  --uniprot_fasta ../0data/uniprot_reviewed_canonical_and_isoform.fasta \
  --name movie_refined
  # Filtered
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
  --pacbio_peptides ../14Metamorpheus_Sample_filtered/search_results/Task1SearchTask/AllPeptides.movie.filtered.psmtsv \
  --gencode_fasta ../10Protein_Hybrid_Database/gencode_protein.fasta \
  --uniprot_fasta ../0data/uniprot_reviewed_canonical_and_isoform.fasta \
  --name movie_filtered
  # High Confidence
  python /mnt/usb1/chenhu/work/Long-Read-Proteogenomics-main/modules/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
  --pacbio_peptides ../15Metamorpheus_Sample_Hybrid/search_results/Task1SearchTask/AllPeptides.movie.hybrid.psmtsv \
  --gencode_fasta ../10Protein_Hybrid_Database/gencode_protein.fasta \
  --uniprot_fasta ../0data/uniprot_reviewed_canonical_and_isoform.fasta \
  --name movie_hybrid
  ```

- Results

```tex
.
├── 0data
│   ├── barcode.fasta
│   ├── gencode.gtf
│   ├── gencode_transcript.fasta
│   ├── gencode_translations.fasta
│   ├── gencode.v35.annotation.gff3
│   ├── Human_Hexamer.tsv
│   ├── Human_logitModel.RData
│   ├── PeptideAtlas
│   ├── reference.fasta
│   ├── reference.fasta.fai
│   ├── short_reads2.fofn
│   ├── short_reads.fofn
│   └── uniprot_reviewed_canonical_and_isoform.fasta
├── 10Protein_Hybrid_Database
│   ├── gencode_isoname_clusters.tsv
│   ├── gencode_protein.fasta
│   ├── movie_cds_high_confidence.gtf
│   ├── movie_high_confidence_genes.tsv
│   ├── movie_hybrid.fasta
│   └── movie_refined_high_confidence.tsv
├── 11Metamorpheus_GENCODE
│   ├── search_results
│   ├── settings
│   └── toml
├── 12Metamorpheus_UniProt
│   ├── search_results
│   ├── settings
│   └── toml
├── 13Metamorpheus_Sample
│   ├── search_results
│   ├── settings
│   └── toml
├── 14Metamorpheus_Sample_filtered
│   ├── search_results
│   ├── settings
│   └── toml
├── 15Metamorpheus_Sample_Hybrid
│   ├── search_results
│   ├── settings
│   └── toml
├── 16Metamorpheus_Algorithm
│   ├── search_results
│   ├── settings
│   └── toml
├── 17Peptide_Analysis
├── 18Track_Visualization
│   ├── Peptide
│   ├── Protein
│   └── reference
├── 19Multiregion_BED_generation
│   ├── movie_filtered_ucsc_multiregion.bed
│   ├── movie_high_confidence_ucsc_multiregion.bed
│   └── movie_refined_ucsc_multiregion.bed
├── 1SMARTLink_ccs
│   └── jurkat_merged.ccs.bam
├── 20Accession_Mapping
│   ├── accession_map_gencode_uniprot_pacbio.tsv
│   └── accession_map_stats.tsv
├── 21Protein_Group_Comparison
│   ├── ProteinInference_GENCODE_PacBio_comparisons.xlsx
│   ├── ProteinInference_GENCODE_UniProt_comparisons.xlsx
│   └── ProteinInference_UniProt_PacBio_comparisons.xlsx
├── 22Noveli_Peptides
│   ├── movie_filtered.pacbio_novel_peptides_to_gencode.tsv
│   ├── movie_filtered.pacbio_novel_peptides_to_uniprot.tsv
│   ├── movie_filtered.pacbio_novel_peptides.tsv
│   ├── movie_hybrid.pacbio_novel_peptides_to_gencode.tsv
│   ├── movie_hybrid.pacbio_novel_peptides_to_uniprot.tsv
│   ├── movie_hybrid.pacbio_novel_peptides.tsv
│   ├── movie_refined.pacbio_novel_peptides_to_gencode.tsv
│   ├── movie_refined.pacbio_novel_peptides_to_uniprot.tsv
│   └── movie_refined.pacbio_novel_peptides.tsv
├── 2IsoSeq3
│   ├── jurkat.collapsed.abundance.txt
│   ├── jurkat.collapsed.fasta
│   ├── jurkat.collapsed.gff
│   ├── jurkat.collapsed.renamed.fasta
│   └── movie.collapsed.gff
├── 3SQANTI3
│   ├── ensg_gene.tsv
│   ├── enst_isoname.tsv
│   ├── filtered_sqanti3_classification.tsv
│   ├── filtered_sqanti3_corrected.fasta
│   ├── filtered_sqanti3_corrected.gtf
│   ├── gene_ensp.tsv
│   ├── gene_isoname.tsv
│   ├── gene_lens.tsv
│   ├── GMST
│   ├── isoname_lens.tsv
│   ├── kallisto_output
│   ├── movie_classification.5degfilter.tsv
│   ├── movie_corrected.5degfilter.fasta
│   ├── movie_corrected.5degfilter.gff
│   ├── protein_coding_genes.txt
│   ├── refAnnotation_sqanti3.genePred
│   ├── RTS
│   ├── sqanti3_classification.txt
│   ├── sqanti3_corrected.fasta
│   ├── sqanti3_corrected.genePred
│   ├── sqanti3_corrected.gtf
│   ├── sqanti3_corrected.gtf.cds.gff
│   ├── sqanti3_junctions.txt
│   ├── sqanti3.params.txt
│   ├── STAR_index
│   └── STAR_mapping
├── 4CPAT
│   ├── 4CPAT.no_ORF.txt
│   ├── 4CPAT.ORF_prob.best.tsv
│   ├── 4CPAT.ORF_prob.tsv
│   ├── 4CPAT.ORF_seqs.fa
│   └── 4CPAT.r
├── 5orf_calling
│   ├── all_orfs_mapped.tsv
│   ├── gene_level_tab.tsv
│   ├── orfcalling_best_orf.tsv
│   ├── pb_gene.tsv
│   └── sqanti_isoform_info.tsv
├── 6DB_Generation
│   ├── movie_orf_refined.fasta
│   └── movie_orf_refined.tsv
├── 7DB_Annotation
│   ├── gencode.cds_renamed_exon.gtf
│   ├── gencode.transcript_exons_only.gtf
│   ├── movie.cds_renamed_exon.genePred
│   ├── movie.cds_renamed_exon.gtf
│   ├── movie.sqanti_protein_classification.tsv
│   ├── movie.transcript_exons_only.genePred
│   ├── movie.transcript_exons_only.gtf
│   ├── pacbio_database.gtf
│   └── refAnnotation_movie.genePred
├── 8SQANTI_Protein
│   ├── movie.sqanti_protein_classification.tsv
│   └── refAnnotation_movie.genePred
└── 9Protein_Classification
    ├── gc_exon_chain_strings_for_cds_containing_transcripts.tsv
    ├── gencode_exons_for_cds_containing_ensts.bed
    ├── gencode_exons_for_cds_containing_ensts.bed_merged.bed
    ├── movie.classification_filtered.tsv
    ├── movie.filtered_protein.fasta
    ├── movie_genes.tsv
    ├── movie_orf_refined_gene_update.tsv
    ├── movie.protein_classification_w_meta.tsv
    ├── movie.protein_refined.fasta
    ├── movie.sqanti_protein_classification_w_5utr_info.tsv
    ├── movie_unfiltered.protein_classification.tsv
    ├── movie_w_class_info.tsv
    ├── movie_with_cds_filtered.gtf
    ├── movie_with_cds_refined.gtf
    └── pb_5utr_categories.tsv
```

## Usage

```
module load nextflow/21.10.5 singularity
nextflow run main.nf --config conf/Test.config
```
