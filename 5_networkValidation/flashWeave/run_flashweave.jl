using FlashWeave

data_path = "muo_datampa4.tsv"
meta_data_path = "muo_meta_data_mpa4.tsv"

netw_results = learn_network(data_path, meta_data_path, sensitive=true, heterogeneous=false, alpha=0.01)

save_network("muo_mv.edgelist", netw_results)
save_network("muo_mv.gml", netw_results)