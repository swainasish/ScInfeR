

process_marker_file <- function(ct_marker_df,subtype_stat,subtype_info){
  if (subtype_stat==F) {
    return(list(ct=ct_marker_df,st=NULL))
  }else if (subtype_stat==T){
    uniq_cts <- base::unique(ct_marker_df$celltype)
    subtype_pre_in <- names(subtype_info)
    subtype_dic = list()
    all_sts = c()
    for (i in uniq_cts) {
      if (i %in% subtype_pre_in) {
        subtype_dic[[i]]$st_status=T
        subtype_dic[[i]]$st_df = subset(ct_marker_df,celltype %in% subtype_info[[i]])
        for (j in subtype_info[[i]]){
          all_sts = append(all_sts,j)
        }
      }else{
        subtype_dic[[i]]$st_status=F
        subtype_dic[[i]]$st_df=NULL
      }
    }
  }
  ct_marker_df_no_st = subset(ct_marker_df, celltype %notin% all_sts)
  return(list(ct=ct_marker_df_no_st,st=subtype_dic))
}

subtype_classification <- function(mat_st,umap_proj,st_info,own_weight=0.5,n_neighbor=30){
  mat_st <- t(data.frame(mat_st))
  ncells <- dim(mat_st)[1]
  ngenes <- dim(mat_st)[2]
  #calculate hybrid matrix of own and neighbor expression
  #Message passing layer from graph neural net using M_adj
  n_near_coor <- knnx.index(umap_proj, umap_proj, k=n_neighbor, algo="kd_tree")
  st_mat_neighbor <- data.frame(matrix(nrow = ncells,ncol = ngenes))
  for (ce in 1:ncells) {
    index_nums <-  n_near_coor[ce,]
    st_mat_neighbor[ce,] <- colMeans(mat_st[index_nums,])
  }
  p=3
  w=0.5
  for (itr in 1:p){
    for (ce in 1:ncells) {
      index_nums <-  n_near_coor[ce,]
      st_mat_neighbor[ce,] <- st_mat_neighbor[ce,] + w * colMeans(mat_st[index_nums,])
      w = w/itr
    }
  }
  nei_weight = 1-own_weight
  # cat("nei_weight:",nei_weight,"n_nei:",n_neighbor,"\n")
  hybrid_mat <- own_weight*mat_st + nei_weight*st_mat_neighbor
  colnames(hybrid_mat) <- st_info[,"celltype"]
  #minmax transformation
  hybrid_mat_norm <- (t(t(hybrid_mat)/colSums(hybrid_mat)))*10000
  #calculate colmean
  colmean_ct <- t(apply(hybrid_mat_norm,1, function(x) tapply(x,colnames(hybrid_mat_norm),mean)))
  st_call <- colnames(colmean_ct)[apply(colmean_ct,1,which.max)]
  return(st_call)
}

core_func<- function(expression_matrix ,
                     group_annt,
                     ct_marker_df,
                     umap_embd,
                     subtype_present = F,
                     subtype_info=NULL,
                     own_weightage = 0.5,
                     n_neighbor=10){
  #create the log file
  tmp <- file.path(getwd(),"scinfer.log")
  log_file <- log_open(tmp)
  log_print("Log file working1",console = F)
  # log_close()

  #create_expression_matrix
  n_cell = dim(expression_matrix)[2]
  n_gene = dim(expression_matrix)[1]
  cat("Expression matrix extracted from Seurat object having Gene_num: ",n_gene,", Sample_num: ",n_cell,"\n")
  #differential expression analysis using presto
  clus_marker_df <- wilcoxauc(expression_matrix, group_annt)
  uniq_clus = base::unique(group_annt)
  all_gene_names <- rownames(expression_matrix)
  ct_marker_df <- data.frame(ct_marker_df)
  common_genes <- intersect(all_gene_names,ct_marker_df[,"marker"])
  input_m_length <- length(ct_marker_df[,"marker"])
  detected_m_length <- length(common_genes)
  cat("Detected",detected_m_length,"/",length(unique(ct_marker_df[,"marker"])),"number of genes in the input marker file","\n")
  ct_marker_df <- subset(ct_marker_df,marker %in% common_genes)
  processed_marker_file <- process_marker_file(ct_marker_df,subtype_present,subtype_info)
  result_df = data.frame(matrix(nrow=n_cell,ncol=2))
  colnames(result_df) = c("celltype","subtype")
  #cluster wise annotation
  for (clus in uniq_clus) {
    clus_bool = clus == group_annt
    clus_marker = clus_marker_df[clus_marker_df$"group"==clus,]
    top_n_auc_genes = clus_marker[order(clus_marker$"auc",decreasing = T),]$"feature"[1:10]
    top_genes = clus_marker[order(clus_marker$"auc",decreasing = T),][1:10,]
    valid_top_genes = sum(top_genes$"auc">0.8)
    log_print(paste("cluster number: ",clus),console = F)
    log_print(paste("Genes having > auc of 0.80",valid_top_genes),console = F)
    log_print(top_genes,console = F)
    ct_genes = processed_marker_file$ct$"marker"
    concate_both = unique(c(top_n_auc_genes,ct_genes))
    master_df = processed_marker_file$ct
    colnames(master_df) = c("celltype","feature","weight")
    #retreieve expression mat
    clus_exp_mat = expression_matrix[,clus_bool]
    clus_exp_mat = as.data.frame(clus_exp_mat[concate_both,])
    #step1: cosine similarity
    cosine_simi = cosine(as.matrix(t(clus_exp_mat)))
    cosine_simi[is.na(cosine_simi)] <- 0 #fill NA's with 0
    rownames(cosine_simi) <- concate_both
    colnames(cosine_simi) <- concate_both
    cosine_simi = rowMeans(cosine_simi[ct_genes,top_n_auc_genes])
    cosine_simi = data.frame(names(cosine_simi),cosine_simi)
    colnames(cosine_simi) = c("feature","cosine")
    #step2: find AUC value
    auc_ct_gene = subset(clus_marker, feature %in% ct_genes )[,c("feature","auc")]
    rownames(auc_ct_gene) <- auc_ct_gene$feature
    auc_ct_gene <- auc_ct_gene[master_df$feature,]
    #merge_all_information
    master_df <- cbind(master_df,cosine_simi[,"cosine"],auc_ct_gene[,"auc"])
    colnames(master_df) <- c("celltype","feature","weight","cosine","auc")
    log_print("printing masterdf",console = F)
    log_print(master_df,console = F)
    if(valid_top_genes<9){
      master_df$cosine<- 1
    }
    log_print("printing masterdf:corrected",console = F)
    log_print(master_df,console = F)
    #   master_df = merge(master_df,cosine_simi,by="feature")
    #   master_df = merge(master_df,auc_ct_gene,by="feature")
    #   master_df = master_df[!duplicated(master_df),]

    #   #step3 multiply all
    master_df_summarise = master_df %>% dplyr::group_by(celltype) %>% dplyr::summarize(weight_func=sum(cosine*auc*weight))
    #print(master_df_summarise)
    celltype_sums = master_df %>% dplyr::group_by(celltype) %>% dplyr::summarize(weight_sum=sum(weight))
    aggrgate_df = merge(master_df_summarise,celltype_sums,by="celltype")
    aggrgate_df$"norm_weight" = aggrgate_df[,"weight_func"] / aggrgate_df[,"weight_sum"]
    #checking for subtype present or not in the clusters
    max_val <- max(aggrgate_df[,"norm_weight"])
    threshold <- max_val - 0.02
    st_bool <- aggrgate_df[,"norm_weight"] > threshold
    st_avl <-  sum(st_bool)
    log_print(aggrgate_df,console = F)
    if (st_avl >1) {
      log_print(cat("Subtype possibility in",clus,"\n"),console = F)
      st_clus_cts <- aggrgate_df[,"celltype"][st_bool]
      st_m_df <- subset(ct_marker_df, celltype %in% st_clus_cts)
      st_mat <- clus_exp_mat
      st_mat <- st_mat[st_m_df[,"marker"],]
      umap_embd_st <- umap_embd[clus_bool,]
      log_print(unique(st_m_df$celltype),console = F)

      st_predictions_clus <- subtype_classification(st_mat,umap_embd_st,
                                                    st_m_df,own_weigh=own_weightage,
                                                    n_neighbor=n_neighbor)
      result_df[which(clus_bool),1] <- st_predictions_clus
    }else if(st_avl==1){
      final_call = aggrgate_df$"celltype"[which.max(aggrgate_df$norm_weight)]
      result_df[which(clus_bool),1] = final_call}}
  #subtype analysis-------------------------------------------------
  if(subtype_present == F){
    log_print("No subtype information was provided",console = F)
  }else if(subtype_present == T){
    ct_res = result_df[,1]
    for (ct in names(processed_marker_file$'st')) {
      st_status <- processed_marker_file[["st"]][[ct]][['st_status']]
      if (st_status==T) {
        st_m_df <- processed_marker_file[["st"]][[ct]][['st_df']]
        st_bool <- which(ct_res==ct)
        st_mat <- expression_matrix[st_m_df[,"marker"],st_bool]
        umap_embd_st <- umap_embd[st_bool,]
        st_predictions <- subtype_classification(st_mat,umap_embd_st,st_m_df,own_weigh=own_weightage,n_neighbor=n_neighbor)
        result_df[st_bool,2] <- st_predictions
      }else if (st_status==F){
        st_bool <- which(ct_res==ct)
        result_df[st_bool,2] <- result_df[st_bool,1]
      }
    }}
  log_close()
  return(a=result_df)
}
get_marker_from_ref_matrix <- function(exp_mat,
                                       annotations,
                                       umap_cor,
                                       num_marker_per_ct=15,
                                       Local_weightage = 0.5,
                                       n_local = 2,
                                       auc_threshold=0.75){
  # expr = GetAssayData(object = seurat_obj, assay = assay_name, slot = slot_name)
  expr = exp_mat
  Global_weightage = 1-Local_weightage
  #get global markers
  global_markers = wilcoxauc(expr, annotations)
  global_markers = global_markers[,c("feature","group","auc")]
  #get local markers
  uniq_cts = unique(annotations)
  #find neighbors in umap space
  #adjancency matrix construction
  umap_cor=umap_cor
  mean_centroid_umap = data.frame(matrix(nrow=length(uniq_cts),ncol=2))
  rownames(mean_centroid_umap) = uniq_cts
  for (ct in uniq_cts) {
    umap_mean = colMeans(umap_cor[annotations %in% ct, ])
    mean_centroid_umap[ct,]=umap_mean
  }
  n_near_cell <- knnx.index(mean_centroid_umap, mean_centroid_umap, k=n_local+1, algo="kd_tree")
  rownames(n_near_cell) <- rownames(mean_centroid_umap)
  for (ct in uniq_cts) {
    n_near_cell[ct,] = rownames(n_near_cell)[as.integer(n_near_cell[ct,])]
  }
  #get local marker
  # marker_df_scinfer = data.frame(matrix(nrow=N,ncol=length(uniq_cts)))
  # colnames(marker_df_scinfer)=uniq_cts
  ct_marker <- data.frame(matrix(nrow =0,ncol=3 ))
  colnames(ct_marker) <- c("celltype","marker","weight")
  for (ct in uniq_cts) {
    N = num_marker_per_ct
    local_markers = wilcoxauc(expr, annotations,n_near_cell[ct,])
    local_markers_ct = local_markers[local_markers$"group" %in% ct,c("feature","group","auc")]
    global_markers_ct = global_markers[global_markers$"group" %in% ct,]
    #remove duplicate gene names if present
    local_markers_ct = local_markers_ct[!duplicated(local_markers_ct$"feature"),]
    colnames(local_markers_ct)=c("feature","group","L-auc")
    global_markers_ct = global_markers_ct[!duplicated(global_markers_ct$"feature"),]
    colnames(global_markers_ct)=c("feature","group","G-auc")
    concate_markers = merge(local_markers_ct,global_markers_ct,by="feature")
    weightage_auc = concate_markers$"L-auc" * Local_weightage + concate_markers$"G-auc"*Global_weightage
    concate_markers[,"weightage_auc"]=weightage_auc
    num_pass = sum(concate_markers$weightage_auc>auc_threshold)
    concate_markers = concate_markers[concate_markers$weightage_auc>auc_threshold,]
    #print(ct)
    #print(dim(concate_markers))
    n_pass_genes = dim(concate_markers)[1]
    if(n_pass_genes<N){
      N=n_pass_genes
    }
    top_genes = sort(concate_markers$weightage_auc,decreasing = T,
                     index.return=T)$ix[0:N]
    # print(ct)
    # print(concate_markers[sort(concate_markers$weightage_auc,decreasing = T,
    #                            index.return=T)$ix[0:5],])
    #print( concate_markers[top_genes,])
    top_genes = concate_markers[top_genes,"feature"]
    top_gene_len = length(top_genes)
    if(top_gene_len>0){
      temp_mat = data.frame(matrix(nrow = top_gene_len,ncol = 3))
      colnames(temp_mat) = c("celltype","marker","weight")
      temp_mat$celltype = ct
      temp_mat$marker = top_genes
      temp_mat$weight=1
      ct_marker = rbind(ct_marker,temp_mat)
    }
    #marker_df_scinfer[,ct]=top_genes
  }
  #make the format for scInfeR
  # ct_marker <- data.frame(matrix(nrow =N*length(uniq_cts),ncol=3 ))
  # colnames(ct_marker) <- c("celltype","marker","weight")
  # index_num=1
  # for (ct in uniq_cts) {
  #   index_end = index_num+(N-1)
  #   ct_marker[index_num:index_end,1]=ct
  #   gene_names = marker_df_scinfer[,ct]
  #   ct_marker[index_num:index_end,2]=gene_names
  #   ct_marker[index_num:index_end,3]=1
  #   index_num=index_end+1
  # }
  return(ct_marker)
}
predict_celltype_scRNA_seurat <- function(s_object,
                                          group_annt,
                                          ct_marker_df,
                                          subtype_present = F,
                                          subtype_info = F,
                                          assay_name="RNA",
                                          slot_name="counts",
                                          reduction="umap",
                                          own_weightage = 0.5,
                                          n_neighbor=10){
  t1 = Sys.time()
  DefaultAssay(s_object) <- assay_name
  expression_matrix <- GetAssayData(object = s_object, assay = assay_name, slot = slot_name)
  umap_embd <- Embeddings(s_object,reduction=reduction)
  ct_output <- core_func(expression_matrix ,
                         group_annt,
                         ct_marker_df,
                         umap_embd,
                         subtype_present = subtype_present,
                         subtype_info=subtype_info,
                         own_weightage = own_weightage,
                         n_neighbor=n_neighbor)
  t2= Sys.time()
  cat("scInfeR took:",t2-t1,"seconds","\n")
  return(ct_output)
}

predict_celltype_scATAC <- function(gene_act_mat,
                                    group_annt,
                                    ct_marker_df,
                                    umap_cord,
                                    subtype_present = F,
                                    subtype_info = F,
                                    own_weightage = 0.5,
                                    n_neighbor=10){
  t1 = Sys.time()
  expression_matrix <- gene_act_mat
  umap_embd <- umap_cord
  ct_output <- core_func(expression_matrix ,
                         group_annt,
                         ct_marker_df,
                         umap_embd,
                         subtype_present = subtype_present,
                         subtype_info=subtype_info,
                         own_weightage = own_weightage,
                         n_neighbor=n_neighbor)
  t2= Sys.time()
  cat("scInfeR took:",t2-t1,"seconds","\n")
  return(ct_output)
}


predict_celltype_spatial <- function(expression_matrix,
                                    group_annt,
                                    ct_marker_df,
                                    spatial_cord,
                                    subtype_present = F,
                                    subtype_info = F,
                                    own_weightage = 0.5,
                                    n_neighbor=10){
  t1 = Sys.time()
  expression_matrix <- expression_matrix
  umap_embd <- spatial_cord
  ct_output <- core_func(expression_matrix ,
                         group_annt,
                         ct_marker_df,
                         umap_embd,
                         subtype_present = subtype_present,
                         subtype_info=subtype_info,
                         own_weightage = own_weightage,
                         n_neighbor=n_neighbor)
  t2= Sys.time()
  cat("scInfeR took:",t2-t1,"seconds","\n")
  return(ct_output)
}

