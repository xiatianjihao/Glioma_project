
fdir="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/files"

#V1
cellphonedb method statistical_analysis $fdir/GBM_meta.txt $fdir/GBM_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db

cellphonedb method statistical_analysis $fdir/PBMC_meta.txt $fdir/PBMC_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db

#V2
cellphonedb method statistical_analysis $fdir/GBM_meta_v2.txt $fdir/GBM_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db

cellphonedb method statistical_analysis $fdir/PBMC_meta_v2.txt $fdir/PBMC_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db


#V3
cellphonedb method statistical_analysis $fdir/GBM_meta_v3.txt $fdir/GBM_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db

cellphonedb method statistical_analysis $fdir/PBMC_meta_v3.txt $fdir/PBMC_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db


#V4
cellphonedb method statistical_analysis $fdir/GBM_meta_v4.txt $fdir/GBM_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db

cellphonedb method statistical_analysis $fdir/PBMC_meta_v4.txt $fdir/PBMC_cells_count.txt --iterations=10 --threads=18 --counts-data gene_name --database /data/active_data/ranxiaojuan/cellphonedb_database/CellTalkDB/cellphonedb_user_2022-07-12-14_12.db
