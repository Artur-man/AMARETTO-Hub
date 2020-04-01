#' AMARETTOHub_InitiateEntities 
#' 
#' This function creates necessary CSV files to import Entities to Neo4j server
#'
#' @param Neo4j_Dir the path to directory that Neo4j files to be stored.
#' @param AMARETTOlink The link to the HTML reports 
#' @param cAMARETTO_Results Community AMARETTO result
#' @param AMARETTO_Results The list of AMARETTO results 
#' @param AMARETTOinits The list of AMARETTO init objects 
#' @param hgtest_tbl_all The list of all Functional Enrichments of Modules tables
#' @param phenotype_tests_all The list of all Phenotype association tables 
#' 
#' @importFrom rlang .data
#'
AMARETTOHub_InitiateEntities <- function(Neo4j_Dir, AMARETTOlink, cAMARETTO_Results, AMARETTO_Results, AMARETTOinits,
                                         hgtest_tbl_all, phenotype_tests_all) 
{
  # If mandatory objects are not provided, stop!
  if(is.null(cAMARETTO_Results) | is.null(AMARETTOinits) | is.null(AMARETTO_Results)){
    stop('cAMARETTO_Results, AMARETTOinits and cAMARETTO_Results objects should be provided!')
  }
  
  # If any of the cohorts of AMARETTO are not provided, stop
  
  # Cohort Names
  Cohorts <- cAMARETTO_Results$cAMARETTOresults$runnames
  if(!all(names(AMARETTOinits)==Cohorts) | !all(names(AMARETTO_Results)==Cohorts)){
    stop('Some cohorts are not provided, please provide AMARETTO results and init objects of all cohorts')
  }
  
  # Websites for Genes and Genesets
  Genesetwebsite <- 'https://www.gsea-msigdb.org/gsea/msigdb/cards/'
  Genewebsite <- 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
  
  # Neo4j directory
  if(!dir.exists(Neo4j_Dir))
    dir.create(Neo4j_Dir)
  
  # Parse Modules, Drivers and Targets
  cat('Preparing CSVs for Modules, Drivers and Targets ... \n')
  for(cohort in Cohorts){

    # Collecting init and results file of AMARETTO of a single cohort
    AMARETTOinit <- AMARETTOinits[[cohort]]
    AMARETTO_Result <- AMARETTO_Results[[cohort]]

    # Extract Modules and Target Genes summaries
    ModuleOverview <- data.frame(Genes = rownames(AMARETTO_Result$ModuleMembership),
                                 Module = AMARETTO_Result$ModuleMembership)
    ModuleOverview <- tibble::as_tibble(ModuleOverview) %>% group_by(ModuleNr) %>%
      summarize(GeneList = paste(sort(.data$Genes),collapse = ', '),
                NumberOfGens = n())
    ModuleOverview <- ModuleOverview %>% dplyr::mutate(ModuleNr = paste0('Module ',ModuleNr))

    # Extract Regulators of Each Module
    RegulatorOverview <- AMARETTO_Result$RegulatoryPrograms
    rownames(RegulatorOverview) <- gsub('_',' ',rownames(RegulatorOverview))
    RegulatorSummary <- apply(RegulatorOverview, 1, nnzero)
    RegulatorList <- apply(RegulatorOverview, 1, function(x){
      temp <- sort(colnames(RegulatorOverview)[which(x!=0)])
      return(paste(temp,collapse = ', '))
    })

    # Collect all and create Module Overview csv files
    ModuleOverview <- ModuleOverview %>% dplyr::mutate(DriverGenes = RegulatorSummary)
    ModuleOverview <- ModuleOverview %>% dplyr::mutate(DriverList = RegulatorList)
    ModuleOverview <- ModuleOverview %>% dplyr::mutate(TargetGenes = .data$NumberOfGens-.data$DriverGenes)
    modules <- gsub('Module ','module', ModuleOverview$ModuleNr)
    Modules_link <- paste0(AMARETTOlink,cohort,'/AMARETTOhtmls/modules/',modules,'.html')
    ModuleOverview <- data.frame(ModuleOverview, Link = Modules_link,
                               ModuleNameLink = paste0('<a href="', Modules_link, '" target="_blank">', ModuleOverview$ModuleNr, '</a>'))
    cohort_file <- paste0(Neo4j_Dir, cohort, '_ModuleOverview.csv')
    write.csv(ModuleOverview, file = cohort_file)

    # Collect Drivers
    DriverList <- apply(RegulatorOverview,1,function(x) x[which(x!=0)])
    names(DriverList) <- NULL
    listnumbers <- lapply(DriverList,length)
    ModuleRow <- rep(ModuleOverview$ModuleNr,unlist(listnumbers))
    DriverGeneData <- data.frame(Genes = names(unlist(DriverList)), Coef = unlist(DriverList),
                                 Module = ModuleRow)

    # Indicate activators (Activator) and suppressors (Repressor)
    DriverGeneData <- tibble::as_tibble(DriverGeneData) %>% dplyr::mutate(Type = ifelse(.data$Coef > 0,'Activator','Repressor'))
    RegulatorAlterations <- AMARETTOinit$RegulatorAlterations$Summary
    RegulatorAlterations <- RegulatorAlterations[as.character(DriverGeneData$Genes),] + 1

    # Indicate CVN or MET Drivers
    DriverGeneData <- tibble::as_tibble(DriverGeneData) %>% dplyr::mutate(CNV = c('No','Yes')[RegulatorAlterations[,'CNV']],
                                                           MET = c('No','Yes')[RegulatorAlterations[,'MET']])

    # Extract Target genes
    TargetGeneData <- data.frame(Genes = rownames(AMARETTO_Result$ModuleMembership),
                                 Coef = 0,
                                 Module = AMARETTO_Result$ModuleMembership,
                                 Type = 'Target')
    rownames(TargetGeneData) <- NULL
    TargetGeneData <- tibble::as_tibble(TargetGeneData) %>% dplyr::mutate(ModuleNr = paste0('Module ',ModuleNr),
                                                           CNV = NA,
                                                           MET = NA)
    colnames(TargetGeneData)[3] <- 'Module'

    # Collect all and create csv for Genes
    GenesData <- rbind(DriverGeneData,TargetGeneData)
    Genes_link <- paste0(Genewebsite,GenesData$Genes)
    GenesData <- data.frame(GenesData, Link = Genes_link,
                            GeneNameLink = paste0('<a href="', Genes_link, '" target="_blank">', GenesData$Genes, '</a>'))
    cohort_file <- paste0(Neo4j_Dir, cohort, '_GenesData.csv')
    write.csv(GenesData, file = cohort_file)
  }

  # Parse Gene sets
  if(!is.null(hgtest_tbl_all)){
    
    # if cohort names do not match with Gene sets of cohorts, stop!
    if(!all(names(hgtest_tbl_all) %in% Cohorts)){
      stop('all cohorts of Gene Sets associations should be from cAMARETTO cohorts')
    }
    
    # Prepare Gene
    cat('Preparing CSVs for Functional Categories ... \n')
    for(cohort in names(hgtest_tbl_all)){
      
      # Collecting hypergeometric test results of AMARETTO modules vs Functional Categories
      # with respect to a single cohort
      hgtest_tbl <- hgtest_tbl_all[[cohort]]
      
      # Collect all and create csv for Functional Categories
      hgtest_tbl$Testset <- gsub('_',' ',hgtest_tbl$Testset)
      hgtest_tbl$Geneset <- gsub('<.*?>','',hgtest_tbl$Geneset)
      Genesets_link <- paste0(Genesetwebsite, gsub(' ','_',hgtest_tbl$Geneset))
      hgtest_tbl <- data.frame(hgtest_tbl, Link = Genesets_link,
                               GenesetNameLink = paste0('<a href="', Genesets_link, '" target="_blank">', hgtest_tbl$Geneset, '</a>'))
      cohort_file <- paste0(Neo4j_Dir, cohort, '_hgtest_tbl.csv')
      write.csv(hgtest_tbl, file = cohort_file)
    }
  }

  # Parse Phenotypes
  if(!is.null(phenotype_tests_all)){
  
    # if cohort names do not match with Gene sets of cohorts, stop!
    if(!all(names(phenotype_tests_all) %in% Cohorts)){
      stop('all cohorts of Phenotype associations should be from cAMARETTO cohorts')
    }
    
    # Prepare
    cat('Preparing CSVs for Clinical Characterizations ... \n')
    for(cohort in names(phenotype_tests_all)){
  
      phenotype_tests <- phenotype_tests_all[[cohort]]
  
      # Collecting hypergeometric test results of AMARETTO modules vs Clinical Characterizations (Phenotypes)
      # with respect to a single cohort
      colnames(phenotype_tests) <- gsub('\\.','',colnames(phenotype_tests))
      levels_phenotype <- levels(factor(phenotype_tests$Phenotypes))
      levels_phenotype <- gsub("\\s*\\([^\\)]+\\)","",as.character(levels_phenotype))
      levels_phenotype <- gsub(",","",as.character(levels_phenotype))
      temp <- factor(phenotype_tests$Phenotypes)
      levels(temp) <- levels_phenotype
      phenotype_tests$Phenotypes <- temp
  
      # Tag survival phenotypes as Worse and Better survival based on Beta statistics
      phenotype_tests$Type <- 0
      beta <- strsplit(phenotype_tests$Descriptive_Statistics[grepl('Beta',phenotype_tests$Descriptive_Statistics)],split = ' ')
      beta <- as.numeric(gsub(',','',sapply(beta,function(x) return(x[2]),simplify = TRUE)))
      beta <- ifelse(beta > 0,'Worse','Better')
      phenotype_tests$Type[grepl('Beta',phenotype_tests$Descriptive_Statistics)] <- beta
  
      # Collect and create csv for Clinical Characterizations
      cohort_file <- paste0(Neo4j_Dir, cohort, '_phenotype_tests.csv')
      write.csv(phenotype_tests, file = cohort_file)
    }
  }
  
  # Parse Communities 
  cat('Preparing CSVs for Communities ... \n')
  Community_key <- Community <- Community_type <- AMARETTOres <- Run_Names <- ModuleNr <- NULL
  ModuleLink <- numTotalEdgesInCommunity <- fractEdgesInVsOut <- CommsizeFrac <- NULL
  cm_gene_df <- suppressWarnings(CommunityAMARETTO::ComRunModGenInfo(cAMARETTO_Results$cAMARETTOresults, 
                                                  cAMARETTO_Results$cAMARETTOnetworkM, cAMARETTO_Results$cAMARETTOnetworkC)) %>% dplyr::select(Community_key, 
                                                                                                                                               Community, Community_type, AMARETTOres, Run_Names, ModuleNr) %>% distinct()
  cm_gene_df$ModuleNr <- gsub('Module_','Module ',cm_gene_df$ModuleNr)
  ComModule <- cm_gene_df
  ComModule <- ComModule %>% group_by(Community_key, Run_Names) %>% 
    summarise(ModuleLinks = paste(ModuleNr, collapse = ", "))
  ComModule <- suppressMessages(reshape2::dcast(ComModule, 
                                                Community_key ~ Run_Names, fill = 0))
  ComModule <- suppressMessages(ComModule %>% left_join(cm_gene_df %>% 
                                                          dplyr::select(Community_key, Community, Community_type) %>% 
                                                          distinct()))
  ComModule <- ComModule %>% dplyr::select(Community, sort(cAMARETTO_Results$cAMARETTOresults$runnames), 
                                           everything()) %>% dplyr::select(-Community_key, -Community_type) 
  ComModule$Community[!grepl('Network',ComModule$Community)] <- paste0('Community ',ComModule$Community[!grepl('Network',ComModule$Community)])
  
  # Collect all and create csv for Community summaries
  communities <- paste0('Community_',1:length(ComModule$Community))
  Communities_link <- paste0(AMARETTOlink,'communities/',communities,'.html')
  ComModule <- data.frame(ComModule, Link = Communities_link, 
                          CommunityNameLink = paste0('<a href="', Communities_link, '" target="_blank">', ComModule$Community, '</a>'))
  cohort_file <- paste0('Neo4j_import/RegulatoryModulesOfCommunities.csv')
  write.csv(ComModule, file = cohort_file)
  
  # Return headers of community tables for importing tables later
  Community_NodeGroups <- colnames(ComModule)[-1]
  Community_NodeGroups <- Community_NodeGroups[!grepl('Link', Community_NodeGroups)]
  
  # Collect all and create csv for Community vs Modules table
  cm_gene_df <- cm_gene_df %>% dplyr::mutate(Community = ifelse(Community_type, 
                                                                Community, paste0('Community ',Community)), 
                                             ModuleNr = paste(Run_Names,ModuleNr, sep = ' ')) %>% 
    dplyr::select(Community, ModuleNr)
  cohort_file <- paste0('Neo4j_import/CommunitiesToModule.csv')
  write.csv(cm_gene_df, file = cohort_file)
  
  # return
  return(Community_Info <- list(Cohorts = Cohorts, Community_NodeGroups = Community_NodeGroups))
}