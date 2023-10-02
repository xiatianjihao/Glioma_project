cd /data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma_vs_PBMC_diff

for id in Monocyte DC MDM
do
/software/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -gmx Msigdb.7.5.1.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk $id.rnk -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out $id
done

cd /data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma_vs_PBMC_diff

for id in CD4_Naive CD4_Stress Treg CD8_Cytotoxicity
do
/software/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -gmx Msigdb.7.5.1.symbols_v2.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk $id.rnk -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out $id
done

