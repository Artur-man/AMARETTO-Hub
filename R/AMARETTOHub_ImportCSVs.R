#' AMARETTOHub_ImportCSVs
#' 
#' This function import CSV files to Neo4j server  
#'
#' @param AMARETTODirectory The AMARETTO Directory
#' @param con_info The list of necessary Neo4j server information: url, username and password 
#' @param Community_Info A list of information on Community AMARETTO: (1) a vector of cohort names and (2) type of nodes in Community AMARETTO
#' 
AMARETTOHub_ImportCSVs <- function(AMARETTODirectory, con_info, Community_Info){
  
  # establish communication with neo4j server
  con <- neo4r::neo4j_api$new(
    url = con_info$url,
    user = con_info$user,
    password = con_info$password)
  
  # Initiate Entity Data 
  Genes <- Modules <- GeneSet <- Phenotype <- Drug <- NULL 
  NBCT <- 10
  
  # Detect existing files in the Neo4j folder
  file_type <- c('Communities','ModuleOverview','Genes','hgtest','phenotype')
  Neo4j_dir <- paste0(AMARETTODirectory,'/Neo4j_import/')
  CSVfile_list <- list.files(Neo4j_dir)
  file_type_exists <- apply(sapply(file_type, function(x) return(stringr::str_detect(CSVfile_list,x))),2,any)
  
  # Import query_header
  query_header <- paste0('USING PERIODIC COMMIT LOAD CSV WITH HEADERS FROM "file:',getwd(),'/',Neo4j_dir)
  
  # Set uniqueness constraints of Entities
  query <- paste0('CREATE CONSTRAINT ON (module:Module) ASSERT module.ModuleName IS UNIQUE ',
                  'CREATE CONSTRAINT ON (gene:Gene) ASSERT gene.GeneName IS UNIQUE ',
                  'CREATE CONSTRAINT ON (geneset:GeneSet) ASSERT geneset.GeneSetName IS UNIQUE ',
                  'CREATE CONSTRAINT ON (phenotype:Phenotype) ASSERT phenotype.PhenotypeName IS UNIQUE',
                  'CREATE CONSTRAINT ON (community:Community) ASSERT community.CommunityName IS UNIQUE')
  suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))
  
  # Importing begins 
  cat('----------------- \n')
  cat('Importing to Neo4j Server at ',con_info$url,' begins ... \n')
  cat('----------------- \n')
  
  # Import Entities for each Cohort
  for(cohort in Community_Info$Cohorts){

    # Import Cohort
    cat('Importing Entities of Cohort:', cohort, '\n')

    # Import Modules
    cat('Importing Modules \n')

    query <- paste0(query_header, cohort,'_ModuleOverview.csv" AS row MERGE ',
                    '(m:Module {ModuleName: "', cohort, ' " + row.ModuleNr, ',
                    'cohort: "', cohort, '"',
                    ', Genes: row.GeneList, DriverList: row.DriverList, TargetGenes: row.TargetGenes, DriverGenes: row.DriverGenes, link: row.Link, namelink: row.ModuleNameLink})')
    suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))

    # Save Modules
    'MATCH (n:Module) RETURN n.ModuleName' %>% neo4r::call_neo4j(con) -> Modules
    Modules <- Modules$n.ModuleName$value

    # Import Genes
    cat('Importing Driver and Target Genes \n')

    query <- paste0(query_header, cohort,'_GenesData.csv" AS row ',
                    'MERGE (g:Gene {GeneName: row.Genes, link: row.Link, namelink: row.GeneNameLink}) ',
                    'MERGE (m:Module {ModuleName: "', cohort, ' " + row.Module}) ',
                    'MERGE (m)<-[:PART_OF {Type: row.Type, Coefficients: toInteger(row.Coef), CNV: row.CNV, MET: row.MET}]-(g);')
    suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))

    # Save Genes
    'MATCH (n:Gene) RETURN n' %>% neo4r::call_neo4j(con) -> Genes
    Genes <- c(Genes$n$GeneName)

    if(file_type_exists['hgtest']){

      # Import Genesets
      cat('Importing Functional Categories \n')

      query <- paste0(query_header, cohort,'_hgtest_tbl.csv" AS row ',
                      'MERGE (module:Module {ModuleName: "', cohort, ' " + row.Testset}) ',
                      'MERGE (geneset:GeneSet {GeneSetName: row.Geneset, NumberOfGenes: toInteger(row.Geneset_length), link: row.Link, namelink: row.GenesetNameLink}) ',
                      'MERGE (module)-[:MODULE_TO_GENESET {GeneOverlap: row.Overlapping_genes, NumberGeneOverlap: toInteger(row.n_Overlapping), Pvalue: toFloat(row.p_value), FDRQvalue: toFloat(row.padj)}]->(geneset);')
      suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))

      # Save Genesets
      'MATCH (n:GeneSet) RETURN n' %>% neo4r::call_neo4j(con) -> GeneSet
      GeneSet <- c(GeneSet$n$GeneSetName)

    }

    if(file_type_exists['phenotype']){

      # Import Phenotypes
      cat('Importing Clinical Characterizations \n')

      query <- paste0(query_header, cohort,'_phenotype_tests.csv" AS row ',
                      'MERGE (module:Module {ModuleName: "', cohort, ' " + row.ModuleNr}) ',
                      'MERGE (phenotype:Phenotype {PhenotypeName: row.Phenotypes}) ',
                      'MERGE (module)-[:MODULE_TO_PHENOTYPE {Test: row.Statistical_Test, Descriptive: row.Descriptive_Statistics, Pvalue: toFloat(row.pvalue), FDRQvalue: toFloat(row.qvalue), Type: row.Type}]->(phenotype);')
      suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))

      # Save Phenotypes
      'MATCH (n:Phenotype) RETURN n' %>% neo4r::call_neo4j(con) -> Phenotype
      Phenotype <- c(Phenotype$n$PhenotypeName)

    }
  }
  
  # Import Communities
  if(file_type_exists['Communities']){
    
    cat('Importing Communities ...\n')
    
    # Import Community Summaries
    query <- paste0(query_header,'RegulatoryModulesOfCommunities.csv" AS row MERGE (m:Community {CommunityName: row.Community, ',
                    paste(paste0(Community_Info$Community_NodeGroups,': row.', Community_Info$Community_NodeGroups), collapse = ', '),
                    ', link: row.Link, namelink: row.CommunityNameLink})')
    suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))
    
    # Import Module vs Community relationships
    # COMMUNITIES AND MODULES
    query <- paste0(query_header,'CommunitiesToModule.csv" AS row
                    MERGE (c:Community {CommunityName: row.Community})
                    MERGE (m:Module {ModuleName: row.ModuleNr})
                    MERGE (c)-[:CONTAINING]->(m);')
    suppressMessages(query_holder <- query %>% neo4r::call_neo4j(con))
    
    # Save Communities 
    'MATCH (n:Community) RETURN n.CommunityName' %>% neo4r::call_neo4j(con) -> Communities
    Communities <- Communities$n.CommunityName$value
  }
  
  # save Entire Entities 
  save(Genes,Modules, GeneSet, Communities, Community_Info, Phenotype, Drug, NBCT, file='AllData.RData')
}
