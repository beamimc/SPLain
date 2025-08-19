

### extract the gene symbol from the DB

get_description <-function(gene_symbol, ref_assembly){

  if (ref_assembly == "hg38") {
    db <- org.Hs.eg.db
  } else if (ref_assembly == "mm39") {
    db <- org.Mm.eg.db
  }
  
  # Get gene symbol
  gene_ensembl <- mapIds(db,
                        keys = gene_symbol,
                        column = "ENSEMBL",
                        keytype = "SYMBOL",
                        multiVals = "first")
  gene_ensembl <- gene_ensembl[[1]]
  
  # Get full gene name/description
  gene_description <- mapIds(db,
                             keys = gene_symbol,
                             column = "GENENAME",
                             keytype = "SYMBOL",
                             multiVals = "first")
  gene_description <-  gene_description[[1]]
  
  description <-paste(gene_ensembl, gene_description)
  return(description)
}
  

get_GO <- function(gene_symbol, ref_assembly){

  if (ref_assembly == "hg38") {
    db <- org.Hs.eg.db
  } else if (ref_assembly == "mm39") {
    db <- org.Mm.eg.db
  }
  
  
  # Map SYMBOL -> GO terms (returns a list if multiVals="list")
  go_terms <- mapIds(db,
                     keys = gene_symbol,
                     column = "GO",
                     keytype = "SYMBOL",
                     multiVals = "list")
  
  # Convert list to data.frame
  go_df <- data.frame(
    SYMBOL = rep(names(go_terms), sapply(go_terms, length)),
    GO = unlist(go_terms)
  )
  go_df$GO_name <- Term(go_df$GO)
  
  # Add ontology if desired
  go_df$ONTOLOGY <- mapIds(db,
                           keys = go_df$GO,
                           column = "ONTOLOGY",
                           keytype = "GO",
                           multiVals = "first")  # Ontology is 1-to-1 with GO terms
  
  return(go_df)
}
