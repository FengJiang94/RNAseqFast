#' Run GO and KEGG analysis
#'
#' This function takes a lits of gene symbols (gene names) and perform GO and KEGG analysis.
#'
#' This function create separate lists of BP, MF, CC, and KEGG enriched terms. This function also generate plots
#' showing the top enriched terms in each ontology (BP, MF, CC, KEGG) and a combined plot showing top enriched
#' GO terms (BP, MF, CC).
#'
#' @param sign_de_gene A list of gene name (gene symbols).
#' @param species "human" or "mouse".
#' @param TermsToPlot Number of top enriched terms to plot. If number of enriched term < *TermsToPlot*, all enriched
#' terms will be plotted.
#' @param qcutoff qvalue cutoff on enrichment tests to report as significant.
#' @param pcutoff Adjusted pvalue cutoff on enrichment tests to report.
#' @param minGSSize Minimal size of genes annotated by Ontology term for testing. Only applied to GO analysis
#' @param maxGSSize Minimal size of genes annotated by Ontology term for testing. Only applied to GO analysis
#' @param Method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", method to adjust
#' pvalue. Only applied to GO analysis.
#' @param simcutoff Similarity cutoff to remove redundant GO terms. GO term with smallest padj value will be kept.
#' @param filename File prefix of the outputs.
#' @import RColorBrewer
#' @import pheatmap
#' @import ggplot2
#' @import RColorBrewer
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @export

############## GO analyis supporting function #############
# looping for each col
GOnKEGG <- function(sign_de_gene,
                    species = "human",
                    TermsToPlot = 10,
                    qcutoff = 0.01,
                    pcutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 500,
                    Method = "BH",
                    simcutoff = 0.7,
                    filename = "filename")
{
  # create a data.frame containing the significance genes
  print("start analysis")
  print("outputs stored in GOnKEGG")
  dir.create(path = "GOnKEGG")
  sign_de_gene <- as.data.frame(sign_de_gene)
  names(sign_de_gene) <- "ID"
  print(paste("perform analysis for ", nrow(sign_de_gene), " genes", sep = ""))
  #convert symbol into entrezid and ensembl id
  if (species == "human"){
    print("use org.Hs.eg.db as species set to human")
    gene_id <- bitr(sign_de_gene$ID, fromType = "SYMBOL", toType = c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
    #GO anakysis for all parameters
    print("perform Go analysis")
    go_ALL <- enrichGO(sign_de_gene$ID, OrgDb = org.Hs.eg.db, ont = 'ALL', minGSSize = minGSSize, maxGSSize = maxGSSize, pAdjustMethod = Method, keyType = 'SYMBOL', pvalueCutoff = pcutoff, qvalueCutoff = qcutoff)
    #Go analysis for biological process
    go_BP <- enrichGO(sign_de_gene$ID, OrgDb = org.Hs.eg.db, ont = "BP", keyType = 'SYMBOL', minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #Go analysis for molecule function
    go_MF <- enrichGO(sign_de_gene$ID, OrgDb = org.Hs.eg.db, ont = "MF", keyType = 'SYMBOL', minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #GO analysis for cellular compartments
    go_CC <- enrichGO(sign_de_gene$ID, OrgDb = org.Hs.eg.db, ont = "CC", keyType = 'SYMBOL', minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #kegg pathway
    print("perform kegg analysis")
    kegg <- enrichKEGG(gene_id$ENTREZID, organism = "hsa", keyType = 'kegg',
                       pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff,
                       use_internal_data = FALSE)
  }
  if (species == "mouse"){
    print("use org.Mm.eg.db as species set to mouse")
    gene_id <- bitr(sign_de_gene$ID, fromType = "SYMBOL", toType = c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
    #GO anakysis for all parameters
    print("perform Go analysis")
    go_ALL <- enrichGO(sign_de_gene$ID, OrgDb = org.Mm.eg.db, ont = 'ALL', minGSSize = minGSSize, maxGSSize = maxGSSize, pAdjustMethod = 'BH', keyType = 'SYMBOL', pvalueCutoff = pcutoff, qvalueCutoff = qcutoff)
    #Go analysis for biological process
    go_BP <- enrichGO(sign_de_gene$ID, OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = minGSSize, maxGSSize = maxGSSize, keyType = 'SYMBOL', pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #Go analysis for molecule function
    go_MF <- enrichGO(sign_de_gene$ID, OrgDb = org.Mm.eg.db, ont = "MF", minGSSize = minGSSize, maxGSSize = maxGSSize, keyType = 'SYMBOL', pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #GO analysis for cellular compartments
    go_CC <- enrichGO(sign_de_gene$ID, OrgDb = org.Mm.eg.db, ont = "CC", minGSSize = minGSSize, maxGSSize = maxGSSize, keyType = 'SYMBOL', pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff)
    #kegg pathway
    print("perform kegg analysis")
    kegg <- enrichKEGG(gene_id$ENTREZID, organism = "mmu", keyType = 'kegg',
                       pvalueCutoff = pcutoff, pAdjustMethod = Method, qvalueCutoff = qcutoff,
                       use_internal_data = FALSE)
  }
  #save the GO results
  write.csv(as.data.frame(go_ALL@result), file = paste("GOnKEGG/",filename,"_GO_ALL_result.csv", sep = ""), quote = F)
  #simplify the GO results,select_fun =min: for similar go terms, only keep the one with minium p.adjust
  CC_simp <- simplify(go_CC, cutoff = simcutoff, by = "p.adjust", select_fun = min)
  BP_simp <- simplify(go_BP, cutoff = simcutoff, by = "p.adjust", select_fun = min)
  MF_simp <- simplify(go_MF, cutoff = simcutoff, by = "p.adjust", select_fun = min)
  #save the simplified go results
  write.csv(as.data.frame(CC_simp), file = paste("GOnKEGG/", filename,"_GO_CC_simp.csv", sep = ''))
  write.csv(as.data.frame(BP_simp), file = paste("GOnKEGG/", filename,"_GO_BP_simp.csv", sep = ''))
  write.csv(as.data.frame(MF_simp), file = paste("GOnKEGG/", filename,"_GO_MF_simp.csv", sep = ''))

  #make dotplots for the top enriched go terms, Defaut: top 10
  print(paste("plot top ", TermsToPlot, " terms for each ontology", sep = ""))
  if (nrow(CC_simp@result) > 0){
    print(paste("plot top ", TermsToPlot, " CC terms", sep = ""))
    pdf(file = paste("GOnKEGG/", filename,"_GO_CC_enrichplot.pdf", sep = ''), width = 12, height = 8)
    p1 <- dotplot(CC_simp, x = "count", showCategory = min(TermsToPlot, nrow(CC_simp)))
    print(p1)
    dev.off()
  }else{
    print("no enriched CC term found")
  }
  #make dotplots for the top 10 enriched go terms,barplots can be made in the same way
  if (nrow(BP_simp@result) > 0){
    print(paste("plot top ", TermsToPlot, " BP terms", sep = ""))
    pdf(file = paste("GOnKEGG/", filename,"_GO_BP_enrichplot.pdf", sep = ''), width = 12, height = 8)
    p1 <- dotplot(BP_simp, x = "count", showCategory = min(TermsToPlot, nrow(BP_simp)))
    print(p1)
    dev.off()
  }else{
    print("no enriched BP term found")
  }
  #make dotplots for the top 10 enriched go terms,barplots can be made in the same way
  if (nrow(MF_simp@result) > 0){
    print(paste("plot top ", TermsToPlot, " MF terms", sep = ""))
    pdf(file = paste("GOnKEGG/", filename,"_GO_MF_enrichplot.pdf", sep = ''), width = 12, height = 8)
    p1 <- dotplot(MF_simp, x = "count", showCategory = min(TermsToPlot, nrow(MF_simp)))
    print(p1)
    dev.off()
  }else{
    print("no enriched MF term found")
  }
  #kegg has entrezIDs for each pathway storted in kegg@result$geneID. now we convert the entrezID into symbol
  #loop each enriched pathway in kegg result
  if (!is.null(kegg)){
    for (m in 1:nrow(kegg@result)) {
      #get the entrezIDs for genes enriched in this pathway
      kegg_gene_id <- as.data.frame(strsplit(kegg@result$geneID[m], "/"))
      names(kegg_gene_id) <- c("ENTREZID")
      #map the entrezIDs to symbols
      kegg_symbol_id <- merge(kegg_gene_id, gene_id, by = "ENTREZID")
      #convert the entrezIDs in kegg into symbols
      kegg@result$geneID[m] <- paste(kegg_symbol_id$SYMBOL, collapse = "/")
      write.csv(as.data.frame(kegg@result), file = paste("GOnKEGG/", filename,"_KEGG_result.csv", sep = ""), quote = F)
    }
    #make dotplots for the top enriched kegg pathways
    if (nrow(kegg@result[kegg@result$p.adjust < 0.05,]) > 0){
      print(paste("plot top ", TermsToPlot, " KEGG terms", sep = ""))
      pdf(file = paste("GOnKEGG/", filename,"_KEGG_enrichplot.pdf", sep = ''), width = 12, height = 8, title = "KEGG Enrichment")
      p2 <- dotplot(kegg, showCategory = min(TermsToPlot, nrow(kegg@result[kegg@result$p.adjust < 0.05,])))
      print(p2)
      dev.off()
    }else{
      print("no enriched KEGG term found")
    }
  }else{
    print("no enriched KEGG term found")
  }
  ##function to make barplot containing top 10 od BP,CC and MF
  print(paste("generate a combined plot for top ", TermsToPlot, " BP, CC, and MF terms", sep = ""))
  go_CC <- as.data.frame(CC_simp@result)
  go_BP <- as.data.frame(BP_simp@result)
  go_MF <- as.data.frame(MF_simp@result)
  if (nrow(go_CC)>0){
    go_CC$ontology <- "CC"
    go_CC_10 <- go_CC[1:min(TermsToPlot,nrow(go_CC)),]}else{go_CC_10 <- go_CC}
  if (nrow(go_BP)>0){
    go_BP$ontology <- "BP"
    go_BP_10 <- go_BP[1:min(TermsToPlot,nrow(go_BP)),]}else{go_BP_10 <- go_BP}
  if (nrow(go_MF)>0){
    go_MF$ontology <- "MF"
    go_MF_10 <- go_MF[1:min(TermsToPlot,nrow(go_MF)),]}else{go_MF_10 <- go_MF}
  #create a data.frame to store the top 10 go terms
  go_enrich_df <- data.frame(ID = c(go_BP_10$ID, go_CC_10$ID, go_MF_10$ID),
                             Description = c(go_BP_10$Description, go_CC_10$Description, go_MF_10$Description),
                             adjusted.P = c(go_BP_10$p.adjust, go_CC_10$p.adjust, go_MF_10$p.adjust),
                             type=factor(c(rep("biological processes",min(TermsToPlot,nrow(go_BP))),
                                           rep("cellular components", min(TermsToPlot,nrow(go_CC))),
                                           rep("molecular function", min(TermsToPlot,nrow(go_MF)))),
                                         levels=c("biological processes",
                                                  "cellular components",
                                                  "molecular function")))
  go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$ID),]
  #number the go terms for barplot as x-axis
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

  #make shortened labels for the barplot x-axis
  go_enrich_df <- go_enrich_df[!is.na(go_enrich_df$Description),]
  go_enrich_df$Description<-as.factor(go_enrich_df$Description)
  labels=(sapply( levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)], shorten_names))
  names(labels) = rev(1:nrow(go_enrich_df))
  head(labels)

  #barplot
  CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
  pdf(file = paste("GOnKEGG/", filename,"GO_enrichplot_bar.pdf", sep = ''), width = 12, height = 8)
  p <- ggplot(data=go_enrich_df, aes(x=number, y=-log10(adjusted.P), fill=type)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() + scale_fill_manual(values = CPCOLS) + theme_classic() + scale_x_discrete(labels = labels) +
    xlab("GO term") + guides(fill=guide_legend(reverse = T)) + labs(title = "The Most Enriched GO Terms")
  print(p)
  dev.off()
}

########## function to shorten the GO term description ####################
shorten_names <- function(x, n_word=4, n_char=40)
{
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40){
      x <- substr(x, 1, 40)}
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]),n_word)],collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}






