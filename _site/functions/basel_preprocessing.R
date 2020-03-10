
# This script is to read and preprocess the Nature Jackson et al. data

library(data.table)
sc_data <- fread("../../sc-targeted-proteomics/data/Nature_cytof/Data_publication/BaselTMA/SC_dat.csv")

table(sc_data$channel)


table(sc_data$core)

length(table(sc_data$core))

length(table(sc_data$id))


sc_mat <- dcast(sc_data, id~channel)


sc_mat_id <- sc_mat$id
sc_mat <- sc_mat[, -1]
sc_mat <- as.matrix(sc_mat)

rownames(sc_mat) <- sc_mat_id
dim(sc_mat)

sc_mat_core <- unlist(lapply(strsplit(rownames(sc_mat), "_"), function(x) {
  paste(x[-length(x)], collapse = "_")
}))

length(sc_mat_core)


sc_mat_core <- data.frame(core = sc_mat_core,
                          core_id = unlist(lapply(strsplit(sc_mat_core, "_"), 
                                                  function(x) {
                                                    paste(x[length(x)], collapse = "_")
                                                  })))
dim(sc_mat_core)
dim(sc_mat)







# Nearest neighbour


all_files <- list.files("../../sc-targeted-proteomics//data/Nature_cytof/results/Basel_NN")

NN <- list()
for (i in 1:length(all_files)) {
  print(all_files[i])
  NN[[i]] <- read.csv(file.path("../../sc-targeted-proteomics/data/Nature_cytof/results/Basel_NN/", all_files[i]))
  idx <- grep(strsplit(all_files[i], "_")[[1]][9], 
              rownames(sc_mat), value = TRUE)
  
  # scell_id <- unlist(lapply(strsplit(idx, "_"), function(x) x[length(x)]))
  # idx <- idx[order(as.numeric(scell_id))]
  
  if (nrow(NN[[i]]) != length(idx)) {
    idx1 <- grepl(paste(strsplit(all_files[i], "_")[[1]][9], collapse = "_"), 
                  rownames(sc_mat)) 
    idx2 <- grepl(paste(strsplit(all_files[i], "_")[[1]][2], collapse = "_"), 
                  rownames(sc_mat))
    idx <- rownames(sc_mat)[idx1 & idx2]
    
    if (nrow(NN[[i]]) != length(idx)) {
      idx3 <- grepl(paste0(gsub("000", "", strsplit(all_files[i], "_")[[1]][8]), "_"), 
                    rownames(sc_mat))
      
      idx <- rownames(sc_mat)[idx1 & idx2 & idx3]
      
      if (nrow(NN[[i]]) != length(idx)) {
        
        idx3 <- grepl(paste0(gsub("000", "", strsplit(all_files[i], "_")[[1]][10]), "_"), 
                      rownames(sc_mat))
        idx <- rownames(sc_mat)[idx1 & idx2 & idx3]
        if (nrow(NN[[i]]) != length(idx)) {
          
          stop("wrong number of cells")
        }
      }
    }
    
    
  } 
  
  scell_id <- unlist(lapply(strsplit(idx, "_"), function(x) x[length(x)]))
  idx <- idx[order(as.numeric(scell_id))]
  rownames(NN[[i]]) <- idx
  
  
}


all(unlist(lapply(NN, rownames)) %in% rownames(sc_mat))



NN_dict <- reshape2::melt(lapply(NN, rownames))




# only keep the cells with NN information
sc_mat <- sc_mat[as.character(NN_dict$value), ]
dim(sc_mat)

all(unlist(lapply(NN, rownames)) %in% rownames(sc_mat))




sc_mat_core <- unlist(lapply(strsplit(rownames(sc_mat), "_"), function(x) {
  paste(x[-length(x)], collapse = "_")
}))

length(sc_mat_core)




## meta data


meta <- fread("../../sc-targeted-proteomics/data/Nature_cytof/Data_publication/BaselTMA/Basel_PatientMetadata.csv")


dim(meta)

all(rownames(sc_mat) == sc_data$id)


all(meta$core %in% sc_mat_core)

all(sc_mat_core %in% meta$core)

# the samples with meta data
common <- intersect(unique(sc_mat_core), unique(meta$core))
length(common)



sc_mat <- sc_mat[sc_mat_core %in% common, ]

sc_mat_core <- sc_mat_core[sc_mat_core %in% common]


meta <- meta[meta$core %in% common, ]

table(meta$clinical_type)



## cell type meta



pg <- fread("../../sc-targeted-proteomics/data/Nature_cytof/Data_publication/BaselTMA/PG_final_k20.csv")

pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(25)] = 1 #B cells
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(19)] = 2 #B and T cells
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(2)] = 3 #T cell
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(6)] = 4 #Macrophages
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(38)] = 5 #T cell
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(70)] = 6 #Macrophage
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(10)] = 7 #Endothelial
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(3)] = 8 #Vimentin hi
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(4)] = 9 #small circular
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(1)] = 10 #small elongated
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(15)] = 11 #Fibronectin hi
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(71)] = 12 #larger elongated
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(36)] = 13 #SMA hi Vimentin

#Tumor cell metaclusters
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(11)] = 14 #hypoxic
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(66,43,32,49,56)] = 15 #apoptotic
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(35,68,8,60)] = 16 #proliferative
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(23,33,26)] = 17 #p53 EGFR
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(47,57,17,50)] = 18 #Basal CK
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(16,69,67)] = 27 #Myoepithalial
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(52,28,64,45)] = 19 #CK7 CK hi Cadherin
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(55,13,40,51,42)] = 20 #CK7 CK
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(59,61,14,48,62,20,24)] = 23 #HR hi CK
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(39,9,12,29,34,22)] = 25 #HR low CK
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(21,65,7)] = 21 #Epithelial low
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(18,46,53,37,31)] = 24 #CK HR
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(58,41,5,27,30)] = 22 #CK low HR low
pg$PhenoGraphBasel_new[pg$PhenoGraphBasel %in% c(44,54,63)] = 26 #CK low HR hi p53



table(pg$PhenoGraphBasel)


table(pg$PhenoGraphBasel_new)

cluster_match <- readRDS("../../data/Nature_cytof/results/cluster_match.rds")
cluster_match <- cluster_match[!duplicated(cluster_match$Basel), ]


# cluster_match <- cluster_match[order(cluster_match$zuri),]
# 
# 
# cluster_name_zuri <- cluster_match$cluster_name
# names(cluster_name_zuri) <- as.character(cluster_match$zuri)
# pg$cluster_name <- cluster_name_zuri[as.character(pg$PhenoGraph)]
# 
# cluster_type_zuri <- cluster_match$cluster_type
# names(cluster_type_zuri) <- as.character(cluster_match$zuri)
# pg$cluster_type <- cluster_type_zuri[as.character(pg$PhenoGraph)]

cluster_name <-  c("B cell",
                   "T & B cells",
                   "T cell",
                   "Marcrophage", 
                   "T cell",
                   "Marcrophage",
                   "Endothelial",
                   "Vim+",
                   "Small circular",
                   "Small elongated",
                   "Fibronectin",
                   "Large elongated",
                   "SMA+Vim-",
                   "Hypoxic",
                   "Apoptotic",
                   "Proliferative",
                   "p53+EGFR+",
                   "Basal CK",
                   "CK7+CK+cadherin+",
                   "CK7+CK+",
                   "Epithelial-",
                   "CK-HR-",
                   "CK+HR++",
                   "CK+HR+",
                   "CK+HR-",
                   "CK-HR+p53+",
                   "Myoepithelial")

cluster_type <- c(rep("Immune", 6),
                  "Endothelial",
                  rep("Stromal", 6),
                  rep("Epithelial", 14))


names(cluster_type) <- as.character(1:27)
pg$cluster_type <- cluster_type[as.character(pg$PhenoGraphBasel_new)]


names(cluster_name) <- as.character(1:27)
pg$cluster_name <- cluster_name[as.character(pg$PhenoGraphBasel_new)]

pg <- as.data.frame(pg)
rownames(pg) <- pg$id

pg <- pg[rownames(sc_mat),]
dim(pg)

mycols_basel_meta <- colors()[c(258,257,86,85,259,81, #green
                                652, #yellow
                                636,520,430,109,68,43,#light blue
                                132,#other dark blue
                                26,#dark blue
                                8,12, #turquouis
                                33,#red
                                551,#dark purple
                                95,#light purple
                                419,#light pink
                                117,#pink
                                52,#burnt orange
                                500,#orange
                                568,#pink orange
                                624, #brown
                                133)]


## To rename scmat's chanel names



rename_chanel <- colnames(sc_mat)

channel_exclude = c("ImageId" ,"CellId" ,"In115 115InIn115Di","Xe134 134XeXe134Di",
                    "Hg202 202HgHg202Di","Pb204 204PbPb204Di","Pb206 206PbPb206Di",
                    "ArAr80 80ArArArAr80Di","phospho Erk12", "10311239Ru96Di Rutheni",
                    "10311240Ru98Di Rutheni","10311241Ru99Di Rutheni", 
                    "10311242Ru100Di Rutheni","10311243Ru101Di Rutheni", 
                    "10311244Ru102Di Rutheni","10311245Ru104Di Rutheni",
                    "Xe126 126XeXe126Di","I127 127II127Di","Xe131 131XeXe131Di",
                    "Pb207 207PbPb207Di","Pb208 208PbPb208Di","EulerNumber",
                    "MajorAxisLength","MinorAxisLength", "Orientation",
                    "10331253Ir191Di Iridium","2971330Dy161Di EpCAM",
                    "Perimeter","1971527Ho165Di bCaten","483739Yb171Di Sox9","Solidity",
                    "Number_Neighbors", "Percent_Touching", "Area", "Eccentricity", "Extent")

# Adhesion markers (done)
rename_chanel[rename_chanel == "6967Gd160Di CD44"] <- "CD44"
rename_chanel[rename_chanel == "1031747Er167Di ECadhe"] <- "ECadherin"


Adhesion_markers <- c("CD44", "ECadherin")


# RTK markers (done)
rename_chanel[rename_chanel == "1021522Tm169Di EGFR"] <- "EGFR"
rename_chanel[rename_chanel == "201487Eu151Di cerbB"] <- "HER2"


RTK_markers <- c("HER2", "EGFR")

# endothelial markers (lack CD31??)
rename_chanel[rename_chanel == "378871Yb172Di vWF"] <- "vWF"


endothelial_markers <- c("CD31", "vWF")

# mesenchymal markers (done)
rename_chanel[rename_chanel == "174864Nd148Di SMA"] <- "SMA"
rename_chanel[rename_chanel == "1921755Sm149Di Vimenti"] <- "Vimentin"
rename_chanel[rename_chanel == "3281668Nd142Di Fibrone"] <- "Fibronectin"

mesenchymal_markers <- c("SMA", "Vimentin", "Fibronectin")

# immune contex (done)
rename_chanel[rename_chanel == "71790Dy162Di CD45"] <- "CD45"
rename_chanel[rename_chanel == "8001752Sm152Di CD3epsi"] <- "CD3"
rename_chanel[rename_chanel == "77877Nd146Di CD68"] <- "CD68"
rename_chanel[rename_chanel == "361077Dy164Di CD20"] <- "CD20"

immune_markers <- c("CD45", "CD3", "CD68", "CD20")


# TF (done)
rename_chanel[rename_chanel == "207736Tb159Di p53"] <- "p53"
rename_chanel[rename_chanel == "117792Dy163Di GATA3"] <- "GATA3"
rename_chanel[rename_chanel == "322787Nd150Di cMyc"] <- "cMyc"
rename_chanel[rename_chanel == "3521227Gd155Di Slug"] <- "Slug"
rename_chanel[rename_chanel == "Nd145Di Twist"] <- "Twist"

cellgrowth_markers <- c("p53", "GATA3", "cMyc", "Slug", "Twist")


# Cell growth (done)
rename_chanel[rename_chanel == "1441101Er168Di Ki67"] <- "Ki67"
rename_chanel[rename_chanel == "phospho S6"] <- "pS6"
rename_chanel[rename_chanel == "phospho mTOR"] <- "pmTOR"
rename_chanel[rename_chanel == "phospho Histone"] <- "pHH3"


cellgrowth_markers <- c("Ki67", "pS6", "pmTOR", "pHH3")


# Cell growth (done)
rename_chanel[rename_chanel == "1441101Er168Di Ki67"] <- "Ki67"
rename_chanel[rename_chanel == "phospho S6"] <- "pS6"
rename_chanel[rename_chanel == "phospho mTOR"] <- "pmTOR"
rename_chanel[rename_chanel == "phospho Histone"] <- "pHH3"


cellgrowth_markers <- c("Ki67", "pS6", "pmTOR", "pHH3")


# hypoxia (done)
rename_chanel[rename_chanel == "92964Er166Di Carboni"] <- "CAIX"
hypoxia_markers <- c("CAIX")

# Epi mark (done)
rename_chanel[rename_chanel == "473968La139Di Histone"] <- "H3K27me3"
epi_markers <- c("H3K27me3")

# Nuclei mark (done)
rename_chanel[rename_chanel == "1261726In113Di Histone"] <- "HistoneH3"
epi_markers <- c("HistoneH3")


# Hormone mark (done)
rename_chanel[rename_chanel == "312878Gd158Di Progest"] <- "PR"
rename_chanel[rename_chanel == "112475Gd156Di Estroge"] <- "ER"
hormone_markers <- c("PR", "ER")


# celldeath
rename_chanel[rename_chanel == "198883Yb176Di cleaved"] <- "cleaved"
celldeath_markers <- c("cleaved")

# CK (done)
rename_chanel[rename_chanel == "651779Pr141Di Cytoker"] <- "CK5"
rename_chanel[rename_chanel == "98922Yb174Di Cytoker"] <- "K7"
rename_chanel[rename_chanel == "3111576Nd143Di Cytoker"] <- "CK19"
rename_chanel[rename_chanel == "971099Nd144Di Cytoker"] <- "CK8/18"
rename_chanel[rename_chanel == "346876Sm147Di Keratin"] <- "K14"
rename_chanel[rename_chanel == "234832Lu175Di panCyto"] <- "panCK"

ck_markers <- c("CK5", "K7", "CK19", "CK8/18", "K14", "panCK")

rename_chanel[rename_chanel == "10331254Ir193Di Iridium"] <- "DNA2"

selected_chanel <- rename_chanel[!rename_chanel %in% channel_exclude]

colnames(sc_mat) <- rename_chanel

saveRDS(sc_mat, file = "../../sc-targeted-proteomics/output/basel_sc_mat.rds")
saveRDS(pg, file = "../../sc-targeted-proteomics/output/basel_pg.rds")
saveRDS(meta, file = "../../sc-targeted-proteomics/output/basel_meta.rds")
saveRDS(selected_chanel, file = "../../sc-targeted-proteomics/output/basel_selected_chanel.rds")


tiff_name_list <- list.files("../../data/Nature_cytof/OMEnMasks/Basel_Zuri_masks/", pattern = "Basel")


map_names <- NULL
skip = FALSE

library(raster)
for (s in 1:length(tiff_name_list)) {
  
  
  str_name <- paste("data/Nature_cytof/OMEnMasks/Basel_Zuri_masks/", tiff_name_list[s], sep = "")
  
  r <- raster(str_name)
  cell_label <- raster::values(r)
  
  ncells <- length(unique(cell_label)) - 1
  
  
  idx <- grep(strsplit(tiff_name_list[s], "_")[[1]][9], 
              rownames(sc_mat), value = TRUE)
  
  
  if (ncells != length(idx)) {
    idx1 <- grepl(paste(strsplit(tiff_name_list[s], "_")[[1]][9], collapse = "_"), 
                  rownames(sc_mat)) 
    idx2 <- grepl(paste(strsplit(tiff_name_list[s], "_")[[1]][2], collapse = "_"), 
                  rownames(sc_mat))
    idx <- rownames(sc_mat)[idx1 & idx2]
    
    if (ncells != length(idx)) {
      idx3 <- grepl(paste0(gsub("000", "", strsplit(tiff_name_list[s], "_")[[1]][8]), "_"), 
                    rownames(sc_mat))
      
      idx <- rownames(sc_mat)[idx1 & idx2 & idx3]
      
      if (ncells != length(idx)) {
        
        idx3 <- grepl(paste0(gsub("000", "", strsplit(tiff_name_list[s], "_")[[1]][10]), "_"), 
                      rownames(sc_mat))
        idx <- rownames(sc_mat)[idx1 & idx2 & idx3]
        if (ncells != length(idx)) {
          
          warning("wrong number of cells")
          skip = TRUE
        }
      }
    }
    
    
  } 
  
  if (!skip) {
    pg_sub <- pg[idx, ]
    sc_mat_subset <- sc_mat[idx, ]
    sample_id <- unique(pg_sub$core)
    
    
    
    print(sample_id)
    
    map_names <- rbind(map_names, c(sample_id, tiff_name_list[s]))
    skip = FALSE
  } else {
    map_names <- rbind(map_names, c(NA, tiff_name_list[s]))
    skip = FALSE
  }
  
}

rownames(map_names) <- map_names[, 1]
saveRDS(map_names, file = "output/basel_map_file_core_names.rds")
