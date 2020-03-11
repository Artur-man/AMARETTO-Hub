getQueryRToCypher <- function(Gene, Module, Phenotype, GeneSet, Community, TransFactor, DriverPert, Drug, Community_Info)
{
  # Index variables for queries
  RsId <- GeneId <- PhenoId <- GeneSetId <- ModuleId <- CommId <- 1
  TransFactorId <- DriverPertId <- DrugId <- 1
  GTMId <- MTPId <- MTGId <- MTDId <- MTTId <- MTDrId <- CTMId <- 1
  inc <- function(x) eval.parent(substitute(x <- x + 1))
  RelationShipVar <- function(RsId) return(paste0('rs',RsId,'=')) 
  RelationShipInd <- function(RsId) return(paste0('rs',RsId)) 
  ModuleInd <- function(ModuleId) return(paste0('module',ModuleId)) 
  CommInd <- function(CommId) return(paste0('community',CommId)) 
  GeneInd <- function(GeneId) return(paste0('gene',GeneId)) 
  PhenoInd <- function(PhenoId) return(paste0('pheno',PhenoId)) 
  GeneSetInd <- function(GeneSetId) return(paste0('gs',GeneSetId)) 
  DriverPertInd <- function(DriverPertId) return(paste0('dp',DriverPertId)) 
  TransFactorInd <- function(TransFactorId) return(paste0('tf',TransFactorId)) 
  DrugInd <- function(DrugInd) return(paste0('drug',DrugInd)) 
  GTMInd <- function(GTMId) return(paste0('gtm',GTMId)) 
  CTMInd <- function(CTMId) return(paste0('ctm',CTMId)) 
  MTPInd <- function(MTPId) return(paste0('mtp',MTPId)) 
  MTGInd <- function(MTGId) return(paste0('mtg',MTGId)) 
  MTDInd <- function(MTDId) return(paste0('mtd',MTDId)) 
  MTTInd <- function(MTTId) return(paste0('mtt',MTTId)) 
  MTDrInd <- function(MTDrId) return(paste0('mtdr',MTDrId))
  
  ################################################################################
  # initiate Clauses, a MATCH clause should at least have 1 relationship query
  # Setting WHERE clauses
  ################################################################################
  
  clause_match <- NULL
  clause_where <- ' WHERE '
  clause_with <- ' WITH '
  clause_table <- list()
  clause_with_list <- NULL

  ###########
  # Communities 
  ###########
  
  clause_match <- paste(clause_match, 'MATCH', RelationShipVar(RsId),
                        getClauseCommunityToModule(CTMInd = CTMInd(CTMId),CaseStudy = Community$Casestudy))
  if(!is.null(Community$Names))
    clause_where <- paste0(clause_where, 'community', ".CommunityName IN ['", paste(Community$Names, sep='', collapse = "','"),"']",' AND ')
  
  clause_table[['communities']] <- paste0('RETURN DISTINCT ',
                                          'community.namelink, ',
                                          paste(paste0('community.',Community_Info$Community_NodeGroups), collapse = ', '))
  
  clause_with_list <- c(clause_with_list, RelationShipInd(RsId))
  clause_with_list <- c(clause_with_list, 'module', 'community')
  inc(RsId); inc(CommId); inc(CTMId)
  
  ###########
  # Modules
  ###########
  
  if(!is.null(Module$Names))
    clause_where <- paste0(clause_where, 'module', ".ModuleName IN ['", paste(Module$Names, sep='', collapse = "','"),"']",' AND ')
  
  clause_where <- paste0(clause_where, 'module', ".cohort IN ['", paste(Module$Cohort, sep='', collapse = "','"),"']",' AND ')
  clause_table[['modules']] <- paste0('RETURN DISTINCT ',
                                      'community.namelink, ',
                                      'module.namelink, ',
                                      'module.Genes AS V12, ',
                                      'module.DriverList AS V13, ',
                                      ' toInt(module.TargetGenes) AS V14, ',
                                      ' toInt(module.DriverGenes) AS V15')
  
  ###########
  # Genes
  ###########
  
  # Set Gene to Module Clause, If Gene(s) are provided
  if(!is.null(Gene$Names)){
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),getClauseGeneToModule(GeneInd = GeneInd(GeneId), 
                                                                                   GTMInd = GTMInd(GTMId)))
    clause_where <- paste0(clause_where, GeneInd(GeneId), 
                           ".GeneName IN ['", paste(Gene$Names, sep='', collapse = "','"),"']",' AND ')
    clause_table[['genes']] <- paste0('RETURN DISTINCT ', 
                                      GeneInd(GeneId), '.namelink, ',
                                      'module.namelink, ',
                                 'community.namelink, ',
                                 GTMInd(GTMId), '.Type AS V3, ',
                                 GTMInd(GTMId), '.MET AS V4, ',
                                 GTMInd(GTMId), '.CNV AS V5')
    clause_with_list <- c(clause_with_list, RelationShipInd(RsId), GeneInd(GeneId), GTMInd(GTMId))
    inc(RsId); inc(GeneId); inc(GTMId)
  }
  
  # If the user requests that all genes should be given, visualize all genes
  if(Gene$ShowAll) {
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),getClauseGeneToModule(GeneInd = GeneInd(GeneId),
                                                                 GTMInd = GTMInd(GTMId)))
    clause_where <- paste0(clause_where, GTMInd(GTMId), 
                     ".Type IN ['", paste(c('Driver(+)','Driver(-)'), sep='', collapse = "','"),"']",' AND ')
    clause_table[['allgenes']] <- paste0('RETURN DISTINCT ',GeneInd(GeneId), '.namelink, ',
                                         'module.namelink, ',
                                    'community.namelink, ',
                                    GTMInd(GTMId), '.Type AS V8, ',
                                    GTMInd(GTMId), '.MET AS V9, ',
                                    GTMInd(GTMId), '.CNV AS V10')
    clause_with_list <- c(clause_with_list, RelationShipInd(RsId), GeneInd(GeneId), GTMInd(GTMId))
    inc(RsId); inc(GeneId); inc(GTMId)
  }
  
  # add a WITH clause and pass all variables to next clauses
  # create a list of variables passed to other clauses
  clause_with <- paste(clause_with, paste(clause_with_list, collapse=', '), sep = ' ')
  
  # delete last AND and add WHERE Clause
  clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
  clause_match <- paste0(clause_match, clause_where, clause_with)
  
  ###########
  # Gene Sets
  ###########
  
  # Set Module to GeneSet Clause, If GeneSet(s) are provided
  if(!is.null(GeneSet$Names)){
    
    # initiate WITH and WHERE clauses for new clauses
    clause_where <- ' WHERE '
    clause_with <- ' WITH '
    
    if(!GeneSet$ShowAllModules){
      
      # match and where clauses
      clause_match <- paste(clause_match, 'OPTIONAL MATCH', RelationShipVar(RsId), getClauseModuleToGeneSet(GeneSetInd = GeneSetInd(GeneSetId),
                                                                                                            MTGInd = MTGInd(MTGId)))
      clause_where <- paste0(clause_where, GeneSetInd(GeneSetId), 
                             ".GeneSetName IN ['", paste(GeneSet$Names, sep='', collapse = "','"),"']",' AND ')
      clause_where <- paste0(clause_where, MTGInd(MTGId), ".FDRQvalue < ", GeneSet$Qvalue,' AND ')
      clause_where <- paste0(clause_where, MTGInd(MTGId), ".Pvalue < ", GeneSet$Pvalue,' AND ')
      
      clause_table[['genesets']] <- paste0('RETURN DISTINCT ','community.namelink, ',
                                           'module.namelink, ',
                                           GeneSetInd(GeneSetId), '.namelink AS V17, ',
                                           MTGInd(MTGId), '.NumberGeneOverlap AS V18, ',
                                           MTGInd(MTGId), '.GeneOverlap AS V19, ',
                                           MTGInd(MTGId), '.Pvalue AS V20, ',
                                           MTGInd(MTGId), '.FDRQvalue AS V21')
      clause_with_list <- c(clause_with_list, RelationShipInd(RsId), MTGInd(MTGId), GeneSetInd(GeneSetId))
      inc(RsId); inc(MTGId); inc(GeneSetId);
      
    } else {
      
      # match and where clauses
      clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),getClauseModuleToGeneSet(ModuleInd = ModuleInd(ModuleId),
                                                                                                 GeneSetInd = GeneSetInd(GeneSetId),
                                                                                                 MTGInd = MTGInd(MTGId)))
      clause_where <- paste0(clause_where, GeneSetInd(GeneSetId), 
                             ".GeneSetName IN ['", paste(GeneSet$Names, sep='', collapse = "','"),"']",' AND ')
      clause_where <- paste0(clause_where, MTGInd(MTGId), ".FDRQvalue < ", GeneSet$Qvalue,' AND ')
      clause_where <- paste0(clause_where, MTGInd(MTGId), ".Pvalue < ", GeneSet$Pvalue,' AND ')
      
      clause_table[['allgenesets']] <- paste0('RETURN DISTINCT ','community.namelink, ',
                                              ModuleInd(ModuleId),'.namelink, ',
                                              GeneSetInd(GeneSetId), '.GeneSetName AS V23, ',
                                              MTGInd(MTGId), '.NumberGeneOverlap AS V24, ',
                                              MTGInd(MTGId), '.GeneOverlap AS V25, ',
                                              MTGInd(MTGId), '.Pvalue AS V26, ',
                                              MTGInd(MTGId), '.FDRQvalue AS V27')
      
      clause_table[['allgenesets_module']] <-  paste0('RETURN DISTINCT ','community.namelink, ',
                                                      ModuleInd(ModuleId),'.namelink, ',
                                                      ModuleInd(ModuleId),'.Genes AS V29,',
                                                      ModuleInd(ModuleId),'.DriverList AS V30,',
                                                      ' toInt(',ModuleInd(ModuleId),'.TargetGenes) AS V31,',
                                                      ' toInt(',ModuleInd(ModuleId),'.DriverGenes) AS V32')
      clause_with_list <- c(clause_with_list, RelationShipInd(RsId), MTGInd(MTGId), GeneSetInd(GeneSetId), ModuleInd(ModuleId))
      inc(RsId); inc(MTGId); inc(ModuleId)
      
    }
    
    # add a WITH clause and pass all variables to next clauses
    # create a list of variables passed to other clauses
    clause_with <- paste(clause_with, paste(clause_with_list, collapse=', '), sep = ' ')
    
    # delete last AND and add WHERE Clause
    clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
    clause_match <- paste0(clause_match, clause_where, clause_with)
    
  }

  #############
  # Phenotypes
  #############
  
  # Set Module to Phenotype Clause, If Phenotype(s) are provided
  if(!is.null(Phenotype$Names) & !Phenotype$ShowAll){
    
    # initiate WITH and WHERE clauses for new clauses
    clause_where <- ' WHERE '
    clause_with <- ' WITH '
    
    if(!Phenotype$ShowAllModules){
      
      # match and where clauses 
      clause_match <- paste(clause_match, 'OPTIONAL MATCH',RelationShipVar(RsId), getClauseModuleToPhenotype(PhenoInd = PhenoInd(PhenoId), 
                                                                                                             MTPInd = MTPInd(MTPId)))
      clause_where <- paste0(clause_where, PhenoInd(PhenoId), 
                             ".PhenotypeName IN ['", paste(Phenotype$Names, sep='', collapse = "','"),"']",' AND ')
      clause_where <- paste0(clause_where, MTPInd(MTPId), ".FDRQvalue < ", Phenotype$Qvalue,' AND ')
      clause_where <- paste0(clause_where, MTPInd(MTPId), ".Pvalue < ", Phenotype$Pvalue,' AND ')
      
      clause_table[['phenotypes']] <- paste0('RETURN DISTINCT ',
                                             'community.namelink, ',
                                             'module.namelink, ',
                                             PhenoInd(PhenoId), '.PhenotypeName AS V34, ',
                                             MTPInd(MTPId), '.Test AS V35, ',
                                             MTPInd(MTPId), '.Pvalue AS V36, ',
                                             MTPInd(MTPId), '.FDRQvalue AS V37, ',
                                             MTPInd(MTPId), '.Descriptive AS V38')
      clause_with_list <- c(clause_with_list, RelationShipInd(RsId), MTPInd(MTPId), PhenoInd(PhenoId))
      inc(RsId); inc(MTPId); inc(PhenoId)
      
    } else {
      
      clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),getClauseModuleToPhenotype(ModuleInd = ModuleInd(ModuleId),
                                                                                                   PhenoInd = PhenoInd(PhenoId),
                                                                                                   MTPInd = MTPInd(MTPId)))
      clause_where <- paste0(clause_where, PhenoInd(PhenoId), 
                             ".PhenotypeName IN ['", paste(Phenotype$Names, sep='', collapse = "','"),"']",' AND ')
      clause_where <- paste0(clause_where, MTPInd(MTPId), ".FDRQvalue < ", Phenotype$Qvalue,' AND ')
      clause_where <- paste0(clause_where, MTPInd(MTPId), ".Pvalue < ", Phenotype$Pvalue,' AND ')
      clause_table[['allphenotypes']] <- paste0('RETURN DISTINCT ','community.namelink, ',
                                                ModuleInd(ModuleId),'.namelink, ',
                                                PhenoInd(PhenoId), '.PhenotypeName AS V46, ',
                                                MTPInd(MTPId), '.Test AS V47, ',
                                                MTPInd(MTPId), '.Pvalue AS V48, ',
                                                MTPInd(MTPId), '.FDRQvalue AS V49, ',
                                                MTPInd(MTPId), '.Descriptive AS V50')
      clause_table[['allphenotypes_module']] <-  paste0('RETURN DISTINCT ','community.namelink, ',
                                                        ModuleInd(ModuleId),'.namelink, ',
                                                        ModuleInd(ModuleId),'.Genes AS V52,',
                                                        ModuleInd(ModuleId),'.DriverList AS V53,',
                                                        ' toInt(',ModuleInd(ModuleId),'.TargetGenes) AS V54,',
                                                        ' toInt(',ModuleInd(ModuleId),'.DriverGenes) AS V55')
      clause_with_list <- c(clause_with_list, RelationShipInd(RsId), MTPInd(MTPId), PhenoInd(PhenoId), ModuleInd(ModuleId))
      inc(RsId); inc(MTPId); inc(ModuleId)
    }
    
    # add a WITH clause and pass all variables to next clauses
    # create a list of variables passed to other clauses
    clause_with <- paste(clause_with, paste(clause_with_list, collapse=', '), sep = ' ')
    
    # delete last AND and add WHERE Clause
    clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
    clause_match <- paste0(clause_match, clause_where, clause_with)
    
  } 
  
  if(Phenotype$ShowAll){
    
    # initiate WITH and WHERE clauses for new clauses
    clause_where <- ' WHERE '
    clause_with <- ' WITH '
    
    # match and where clauses
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId), getClauseModuleToPhenotype(PhenoInd = PhenoInd(PhenoId),
                                                                                                  MTPInd = MTPInd(MTPId)))
    clause_where <- paste0(clause_where, MTPInd(MTPId), ".FDRQvalue < ", Phenotype$Qvalue,' AND ')
    clause_where <- paste0(clause_where, MTPInd(MTPId), ".Pvalue < ", Phenotype$Pvalue,' AND ')
    clause_table[['allphenotypes']] <- paste0('RETURN DISTINCT ','community.namelink, ',
                                              'module.namelink, ',
                                              PhenoInd(PhenoId), '.PhenotypeName AS V40, ',
                                              MTPInd(MTPId), '.Test AS V41, ',
                                              MTPInd(MTPId), '.Pvalue AS V42, ',
                                              MTPInd(MTPId), '.FDRQvalue AS V43, ',
                                              MTPInd(MTPId), '.Descriptive AS V44')
    clause_with_list <- c(clause_with_list, RelationShipInd(RsId), MTPInd(MTPId), PhenoInd(PhenoId))
    inc(RsId); inc(MTPId); inc(PhenoId)
    
    # add a WITH clause and pass all variables to next clauses
    # create a list of variables passed to other clauses
    clause_with <- paste(clause_with, paste(clause_with_list, collapse=', '), sep = ' ')
    
    # delete last AND and add WHERE Clause
    clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
    clause_match <- paste0(clause_match, clause_where, clause_with)
    
  }
    
  ##########################
  # Transcription Factors
  ##########################
    
  clause_where <- ' WHERE '
  
  if(TransFactor$ShowAll | TransFactor$ShowAllNonValid){
    
    # # initiate WITH and WHERE clauses for new clauses
    # clause_where <- ' WHERE '
    # clause_with <- ' WITH '
    
    # match and where clauses
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),
                          getClauseModuleToTransFactor(GeneInd = GeneInd(GeneId), GTMInd = GTMInd(GTMId),
                                                       MTTInd = MTTInd(MTTId),
                                                       TransFactorInd = TransFactorInd(TransFactorId)))
    
    clause_where <- paste0(clause_where, MTTInd(MTTId), ".FDRQvalue < ", TransFactor$Qvalue,' AND ')
    clause_where <- paste0(clause_where, MTTInd(MTTId), ".Pvalue < ", TransFactor$Pvalue,' AND ')
    if(TransFactor$ShowAll){
      clause_where <- paste0(clause_where, GTMInd = GTMInd(GTMId), 
                             ".Type IN ['", paste(c('Driver(+)','Driver(-)'), sep='', collapse = "','"),"']",' AND ')
    }
    
    clause_table[['transfactors']] <- paste0('RETURN DISTINCT ', TransFactorInd(TransFactorId), '.TransFactorName AS V56, ',
                                             'community.namelink, ',
                                             'module.namelink, ',
                                        MTTInd(MTTId), '.NumberGeneOverlap AS V58, ',
                                        MTTInd(MTTId), '.GeneOverlap AS V59, ',
                                        MTTInd(MTTId), '.Pvalue AS V60, ',
                                        MTTInd(MTTId), '.FDRQvalue AS V61')
    
    clause_table[['transfactors_genes']] <- paste0('RETURN DISTINCT ',GeneInd(GeneId), '.namelink, ',
                                                   'module.namelink, ',
                                      'community.namelink, ',
                                      GTMInd(GTMId), '.Type AS V64, ',
                                      GTMInd(GTMId), '.MET AS V65, ',
                                      GTMInd(GTMId), '.CNV AS V66')
    # clause_with_list <- c(clause_with_list, RelationShipInd(RsId), GeneInd(GeneId), GTMInd(GTMId), MTTInd(MTTId), TransFactorInd(TransFactorId))
    inc(RsId); inc(GeneId); inc(GTMId); inc(MTTId); inc(TransFactorId)
    
    # # add a WITH clause and pass all variables to next clauses
    # # create a list of variables passed to other clauses
    # clause_with <- paste(clause_with, paste(clause_with_list, collapse=', '), sep = ' ')
    # 
    # # delete last AND and add WHERE Clause
    # clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
    # clause_match <- paste0(clause_match, clause_where, clause_with)
  }
    
  ##########################
  # Driver Perturbations
  ##########################
    
  if(DriverPert$ShowAll | DriverPert$ShowAllNonValid){
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),
                          getClauseModuleToDriverPert(GeneInd = GeneInd(GeneId), GTMInd = GTMInd(GTMId),
                                                      MTDInd = MTDInd(MTDId), DriverPertInd = DriverPertInd(DriverPertId)))

    if(DriverPert$ValStat){
      ValidStat <- 'escore-pval-padj-zscore'
    } else {
      ValidStat <- c('escore-pval-padj-zscore','escore-pval-padj')
    }
    clause_where <- paste0(clause_where, MTDInd(MTDId), ".ValidationStatus IN ['", 
                           paste(ValidStat, sep='', collapse = "','"),"'] AND ")
    
    if(DriverPert$ShowAll){
      clause_where <- paste0(clause_where, GTMInd = GTMInd(GTMId),
                             ".Type IN ['", paste(c('Driver(+)','Driver(-)'), sep='', collapse = "','"),"']",' AND ')
    }

    clause_table[['driverperts']] <- paste0('RETURN DISTINCT ',
                                            'community.namelink, ',
                                            'module.namelink, ',
                                       DriverPertInd(DriverPertId), '.PerturbationName AS V68, ',
                                       DriverPertInd(DriverPertId), '.PerturbationEntrezID AS V69, ',
                                       DriverPertInd(DriverPertId), '.PerturbationGene AS V70, ',
                                       DriverPertInd(DriverPertId), '.PerturbationConfidence AS V71, ',
                                       DriverPertInd(DriverPertId), '.PerturbationType AS V72, ',
                                       MTDInd(MTDId), '.Pvalue AS V73, ',
                                       MTDInd(MTDId), '.FDRQvalue AS V74, ',
                                       MTDInd(MTDId), '.ES AS V75, ',
                                       MTDInd(MTDId), '.NES AS V76, ',
                                       MTDInd(MTDId), '.NMoreExtreme AS V77, ',
                                       MTDInd(MTDId), '.Size AS V78, ',
                                       MTDInd(MTDId), '.LeadingEdge AS V79, ',
                                       MTDInd(MTDId), '.DriverWeight AS V80, ',
                                       MTDInd(MTDId), '.Zscore AS V81, ',
                                       MTDInd(MTDId), '.GeneType AS V82, ',
                                       MTDInd(MTDId), '.GeneDetailedType AS V83, ',
                                       MTDInd(MTDId), '.PvalueCheck AS V84, ',
                                       MTDInd(MTDId), '.PadjustedCheck AS V85, ',
                                       MTDInd(MTDId), '.ZscoreCheck AS V86, ',
                                       MTDInd(MTDId), '.ESCheck AS V87, ',
                                       MTDInd(MTDId), '.ValidationStatus AS V88'
                                       )

    clause_table[['driverperts_genes']] <- paste0('RETURN DISTINCT ', GeneInd(GeneId), '.namelink, ',
                                             'module.namelink, ',
                                             'community.namelink, ',
                                             GTMInd(GTMId), '.Type AS V91, ',
                                             GTMInd(GTMId), '.MET AS V92, ',
                                             GTMInd(GTMId), '.CNV AS V93')
    inc(RsId); inc(GeneId); inc(GTMId); inc(MTDId); inc(DriverPertId)
  }
    
  #################################
  # Chemical Perturbations (Drugs)
  #################################
  
  if(Drug$ShowAll){
    clause_match <- paste(clause_match, 'MATCH',RelationShipVar(RsId),
                          getClauseModuleToDrug(MTDrInd = MTDrInd(MTDrId), DrugInd = DrugInd(DrugId)))
    if(!is.null(Drug$Name))
      clause_where <- paste0(clause_where, DrugInd = DrugInd(DrugId),
                             ".DrugName IN ['", paste(Drug$Name, sep='', collapse = "','"),"']",' AND ')
    
    clause_where <- paste0(clause_where, DrugInd = DrugInd(DrugId),
                           ".NbClinicalTrials >= ", Drug$ShowNbclinicaltrial[1],' AND ')
    clause_where <- paste0(clause_where, DrugInd = DrugInd(DrugId),
                           ".NbClinicalTrials <= ", Drug$ShowNbclinicaltrial[2],' AND ')
    
    clause_where <- paste0(clause_where, MTDrInd = MTDrInd(MTDrId),
                           ".NbDGIdbDrugGeneInteractions >= ", Drug$ShowDgidb[1],' AND ')
    clause_where <- paste0(clause_where, MTDrInd = MTDrInd(MTDrId),
                           ".NbDGIdbDrugGeneInteractions <= ", Drug$ShowDgidb[2],' AND ')
    
    clause_table[['drugs']] <- paste0('RETURN DISTINCT ','community.namelink, ',
                                      'module.namelink, ',
                                      DrugInd(DrugId), '.DrugName AS V96, ',
                                      MTDrInd(MTDrId), '.DGIDB_Genes AS V105, ', 
                                      MTDrInd(MTDrId), '.NbDGIdbDrugGeneInteractions AS V106, ',
                                      DrugInd(DrugId), '.NbClinicalTrials AS V107, ',
                                      DrugInd(DrugId), '.NbClinicalTrialsPhase1 AS V108, ', 
                                      DrugInd(DrugId), '.NbClinicalTrialsPhase2 AS V109, ',
                                      DrugInd(DrugId), '.NbClinicalTrialsPhase3 AS V110, ',
                                      DrugInd(DrugId), '.NbClinicalTrialsPhase4 AS V111'
    )
    inc(RsId); inc(DrugId); inc(MTDrId);
  }     
  
  # delete last AND and add WHERE Clause
  clause_where <- substr(clause_where, 1, nchar(clause_where)-5)
  if(nchar(clause_where) > 5)
    clause_match <- paste0(clause_match,clause_where)
  
  # Add RETURN clauses of tables before returning
  clause_table <- lapply(clause_table,function(x) return(paste(clause_match,x)))

  # Add RETURN clauses for the graph and return
  clause_match <- paste0(clause_match,' RETURN ',paste('rs',1:(RsId-1), sep='', collapse = ','))
  
  print(clause_match)
  print(clause_table)
  
  return(list(graph = clause_match, table = clause_table))
}

QuerytoTable <- function(G_table, query, con, Community_Info){
  
  # Module Table
  ModuleTable <- data.frame(Community = character(0), Module = character(0), Genes = character(0), `Driver Genes`= character(0), 
                            `# of Target Genes` = numeric(0), `# of Driver Genes` = numeric(0))
  colnames(ModuleTable) <- c('Community', 'Module', 'Target Genes', 'Driver Genes','# of Target Genes', '# of Driver Genes')
  
  # Community Table
  # CommTable <- data.frame(Community = character(0), TCGA_BLCA = character(0),
  #                         TCGA_CESC = character(0), TCGA_ESCA = character(0),
  #                         TCGA_HNSC = character(0), TCGA_LUSC = character(0),
  #                         ImmuneSignatures = character(0), StemSignatures = character(0))
  CommTable <- data.frame(matrix(character(0), ncol = length(Community_Info$Community_NodeGroups) + 1))
  colnames(CommTable) <- c('Community', Community_Info$Community_NodeGroups)
  if(!is.null(query$table[['communities']])) {
    G_table$CommTable <- data.frame(query$table[['communities']]%>%neo4r::call_neo4j(con,type="row"))
    G_table$CommTable <- G_table$CommTable[order(G_table$CommTable[,1]),]
    colnames(G_table$CommTable) <- colnames(CommTable)
  } else {
    G_table$CommTable <- CommTable
  }
  
  # Gene Table
  GeneTable <- data.frame(Gene = character(0), Module = character(0), 
                          Community = character(0), 
                          Type = character(0), MET = character(0),
                          CNV = character(0))
  if(!is.null(query$table[['genes']])) {
    G_table$GeneTable <- data.frame(query$table[['genes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$GeneTable) <- colnames(GeneTable)
  } else {
    G_table$GeneTable <- GeneTable
  }
  # All Gene Table
  if(!is.null(query$table[['allgenes']])) {
    G_table$AllGeneTable <- data.frame(query$table[['allgenes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$AllGeneTable) <- colnames(GeneTable)
  } else {
    G_table$AllGeneTable <- GeneTable
  }
  # TF Gene Table
  if(!is.null(query$table[['transfactors_genes']])) {
    G_table$TransFactorGeneTable <- data.frame(query$table[['transfactors_genes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$TransFactorGeneTable) <- colnames(GeneTable)
  } else {
    G_table$TransFactorGeneTable <- GeneTable
  }
  # DriverPert Gene Table
  if(!is.null(query$table[['driverperts_genes']])) {
    G_table$DriverPertGeneTable <- data.frame(query$table[['driverperts_genes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$DriverPertGeneTable) <- colnames(GeneTable)
  } else {
    G_table$DriverPertGeneTable <- GeneTable
  }
  G_table$GeneTable <- rbind(G_table$GeneTable,G_table$AllGeneTable,G_table$TransFactorGeneTable,G_table$DriverPertGeneTable)
  G_table$GeneTable <- tibble::as_tibble(G_table$GeneTable) %>% dplyr::distinct()
  colnames(G_table$GeneTable) <- colnames(GeneTable)
  
  # Phenotype Table
  PhenotypeTable <- data.frame(Community = character(0), Module = character(0), Phenotype = character(0),
                               `Statistics Test` = character(0), `P-value` = numeric(0),
                               `FDR Q-value` = numeric(0), `Descriptive Statistics` = character(0))
  colnames(PhenotypeTable) <- c('Community', 'Module', 'Phenotype', 'Statistics Test', 'P-value','FDR Q-value',
                                'Descriptive Statistics')
  if(!is.null(query$table[['phenotypes']])) {
    G_table$PhenotypeTable <- data.frame(query$table[['phenotypes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$PhenotypeTable) <- colnames(PhenotypeTable)
    G_table$PhenotypeTable <- na.omit(G_table$PhenotypeTable)
  } else {
    G_table$PhenotypeTable <- PhenotypeTable
  }
  # All Phenotype Table
  if(!is.null(query$table[['allphenotypes']])){
    G_table$AllPhenotypeTable <- data.frame(query$table[['allphenotypes']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$AllPhenotypeTable) <- colnames(PhenotypeTable)
  } else {
    G_table$AllPhenotypeTable <- PhenotypeTable
  }
  G_table$PhenotypeTable <- rbind(G_table$PhenotypeTable,G_table$AllPhenotypeTable)
  G_table$PhenotypeTable <- tibble::as_tibble(G_table$PhenotypeTable) %>% dplyr::distinct()
  colnames(G_table$PhenotypeTable) <- colnames(PhenotypeTable)
  
  # All modules to those phenotypes
  if(!is.null(query$table[['allphenotypes_module']])){
    G_table$AllPhenotypeModuleTable <- data.frame(query$table[['allphenotypes_module']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$AllPhenotypeModuleTable) <- colnames(ModuleTable)
  } else {
    G_table$AllPhenotypeModuleTable <- ModuleTable
  }
  
  # Gene Set Table
  GeneSetsTable <- data.frame(Community = character(0), Module = character(0), `Gene Set` = character(0),
                              `# Genes in Overlap` = numeric(0), `Genes in Overlap` = character(0), 
                              `P-value` = numeric(0),`FDR Q-value` = numeric(0))
  colnames(GeneSetsTable) <- c('Community', 'Module', 'Gene Set', '# Genes in Overlap', 'Genes in Overlap',
                               'P-value','FDR Q-value')
  if(!is.null(query$table[['genesets']])) {
    G_table$GeneSetsTable <- data.frame(query$table[['genesets']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$GeneSetsTable) <- colnames(GeneSetsTable)
    G_table$GeneSetsTable <- na.omit(G_table$GeneSetsTable)
  } else {
    G_table$GeneSetsTable <- GeneSetsTable
  }
  # All Gene Set Table
  if(!is.null(query$table[['allgenesets']])){
    G_table$AllGeneSetsTable <- data.frame(query$table[['allgenesets']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$AllGeneSetsTable) <- colnames(GeneSetsTable)
  } else {
    G_table$AllGeneSetsTable <- GeneSetsTable
  }
  G_table$GeneSetsTable <- rbind(G_table$GeneSetsTable,G_table$AllGeneSetsTable)
  G_table$GeneSetsTable <- tibble::as_tibble(G_table$GeneSetsTable) %>% dplyr::distinct()
  colnames(G_table$GeneSetsTable) <- colnames(GeneSetsTable)
  
  # All modules to those genesets
  if(!is.null(query$table[['allgenesets_module']])){
    G_table$AllGeneSetsModuleTable <- data.frame(query$table[['allgenesets_module']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$AllGeneSetsModuleTable) <- colnames(ModuleTable)
  } else {
    G_table$AllGeneSetsModuleTable <- ModuleTable
  }
  
  # Transcription Factor Table
  TransFactorTable <- data.frame(`Trans. Factor.` = character(0), Community = character(0), Module = character(0), 
                              `# Genes in Overlap` = numeric(0), `Genes in Overlap` = character(0), 
                              `P-value` = numeric(0),`FDR Q-value` = numeric(0))
  colnames(TransFactorTable) <- c('Trans. Factor.', 'Community',  'Module', '# Genes in Overlap', 'Genes in Overlap',
                               'P-value','FDR Q-value')
  if(!is.null(query$table[['transfactors']])){
    G_table$TransFactorTable <- data.frame(query$table[['transfactors']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$TransFactorTable) <- colnames(TransFactorTable)
  } else {
    G_table$TransFactorTable  <- TransFactorTable
  }
  G_table$TransFactorTable <- tibble::as_tibble(G_table$TransFactorTable) %>% dplyr::distinct()
  colnames(G_table$TransFactorTable) <- colnames(TransFactorTable)

  # Driver Discovery Table
  DriverPertTable <- data.frame(Community = character(0), Module = character(0), PerturbationID = character(0), PerturbationEntrezID = character(0),
                                 PerturbationName = character(0), PerturbationType = character(0),
                                 PerturbationConfidence = character(0), Pvalue = character(0), Padjusted = character(0),
                                 ES = character(0), NES = character(0), NMoreExtreme = character(0), Size = character(0),
                                 LeadingEdge = character(0), DriverWeight = character(0), Zscore = character(0),
                                 GeneType = character(0), GeneDetailedType = character(0), PvalueCheck = character(0),
                                 PadjustedCheck = character(0), ZscoreCheck = character(0), ESCheck = character(0),
                                 ValidationStatus = character(0))
  if(!is.null(query$table[['driverperts']])){
    G_table$DriverPertTable <- data.frame(query$table[['driverperts']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$DriverPertTable) <- colnames(DriverPertTable)
  } else {
    G_table$DriverPertTable  <- DriverPertTable
  }
  G_table$DriverPertTable <- tibble::as_tibble(G_table$DriverPertTable) %>% dplyr::distinct()
  colnames(G_table$DriverPertTable) <- colnames(DriverPertTable)
  
  # Drug Table
  DrugTable <- data.frame(Community = character(0), Module = character(0), Compound = character(0),
                                DGIdbDrugGeneInteractions = character(0),
                                NbDGIdbDrugGeneInteractions = numeric(0), NbClinicalTrials = numeric(0), 
                                NbClinicalTrialsPhase1 = numeric(0),NbClinicalTrialsPhase2 = numeric(0),
                                NbClinicalTrialsPhase3 = numeric(0),NbClinicalTrialsPhase4 = numeric(0))
  
  if(!is.null(query$table[['drugs']])){
    G_table$DrugTable <- data.frame(query$table[['drugs']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$DrugTable) <- colnames(DrugTable)
  } else {
    G_table$DrugTable  <- DrugTable
  }
  G_table$DrugTable <- tibble::as_tibble(G_table$DrugTable) %>% dplyr::distinct()
  colnames(G_table$DrugTable) <- colnames(DrugTable)
  
  # All other modules
  if(!is.null(query$table[['modules']])) {
    G_table$ModuleTable <- data.frame(query$table[['modules']]%>%neo4r::call_neo4j(con,type="row"))
    colnames(G_table$ModuleTable) <- colnames(ModuleTable)
  } else {
    G_table$ModuleTable <- ModuleTable
  }
  G_table$ModuleTable <- rbind(G_table$ModuleTable,G_table$AllPhenotypeModuleTable,
                               G_table$AllGeneSetsModuleTable)
  G_table$ModuleTable <- tibble::as_tibble(G_table$ModuleTable) %>% dplyr::distinct()
  colnames(G_table$ModuleTable) <- colnames(ModuleTable)
  return(G_table)
}

getClauseModuleToDrug <- function(ModuleInd = 'module', DrugInd = '', MTDrInd = ''){
  return(paste(
    getNodewithProperty(ModuleInd,'Module'),
    getRelationshipWithProperty(MTDrInd,'MODULE_TO_DRUG'),
    getNodewithProperty(DrugInd,'Drug')
  ))
}

getClauseModuleToDriverPert <- function(GeneInd, GTMInd, ModuleInd ='module',
                                         MTDInd = '', DriverPertInd = '', DTGInd = ''){
  return(paste(
    getNodewithProperty(GeneInd,'Gene'),
    getRelationshipWithProperty(GTMInd,'PART_OF'),
    getNodewithProperty(ModuleInd,'Module'),
    getRelationshipWithProperty(MTDInd,'MODULE_TO_DRIVERPERT'),
    getNodewithProperty(DriverPertInd,'DriverPert'),
    getRelationshipWithProperty(DTGInd,'PERTURBATION'),
    getNodewithProperty(GeneInd,'Gene')
  ))
}

getClauseModuleToTransFactor <- function(GeneInd, GTMInd, Type = NULL, ModuleInd ='module',
                                         MTTInd = '', TransFactorInd = ''){
  return(paste(
    getNodewithProperty(GeneInd,'Gene'),
    getRelationshipWithProperty(GTMInd,'PART_OF','Type',Type),
    getNodewithProperty(ModuleInd,'Module'),
    getRelationshipWithProperty(MTTInd,'MODULE_TO_TRANSFACTOR'),
    getNodewithProperty(TransFactorInd,'TransFactor'),
    '-[]->',
    getNodewithProperty(GeneInd,'Gene')
  ))
}

getClauseModuleToGeneSet <- function(ModuleInd = 'module', GeneSetInd = '', MTGInd = ''){
  return(paste(
    getNodewithProperty(ModuleInd,'Module'),
    getRelationshipWithProperty(MTGInd,'MODULE_TO_GENESET'),
    getNodewithProperty(GeneSetInd,'GeneSet')
  ))
}

getClauseModuleToPhenotype <- function(ModuleInd = 'module', PhenoInd = '', MTPInd = ''){
  return(paste(
    getNodewithProperty(ModuleInd,'Module'),
    getRelationshipWithProperty(MTPInd,'MODULE_TO_PHENOTYPE'),
    getNodewithProperty(PhenoInd,'Phenotype')
  ))
}

getClauseCommunityToModule <- function(CommInd = 'community', CTMInd = '', CaseStudy = NULL, Type = NULL, ModuleInd = 'module'){
  
  return(paste(
    # getNodewithProperty(GeneInd,'Gene','GeneName',GeneName),
    getNodewithProperty(CommInd,'Community','Case',CaseStudy),
    getRelationshipWithProperty(CTMInd,'CONTAINING'),
    getNodewithProperty(ModuleInd,'Module')
  ))
}

getClauseGeneToModule <- function(GeneInd = '', GTMInd = '', Type = NULL, ModuleInd = 'module'){
  
  return(paste(
    # getNodewithProperty(GeneInd,'Gene','GeneName',GeneName),
    getNodewithProperty(GeneInd,'Gene'),
    ifelse(is.null(Type),
           getRelationshipWithProperty(GTMInd,'PART_OF'),
           getRelationshipWithProperty(GTMInd,'PART_OF','Type',Type)
    ),
    getNodewithProperty(ModuleInd,'Module')
  ))
}

getRelationshipWithProperty <- function(edge, label, property = NULL, property_value = NULL){
  if(is.null(label)) 
    stop('Please enter edge bel')
  clause <- paste('-[',edge,':',label, sep="")
  if(is.null(property) | is.null(property_value)){
    clause <- paste(clause,']->',sep='')
  } else {
    clause <- paste(clause,' {',property,": '",property_value,"'}]->",sep="")
  }
  return(clause)
}

getNodewithProperty <- function(node, label, property = NULL, property_value = NULL){
  if(is.null(label)) 
    stop('Please enter node label')
  clause <- paste('(',node,':',label, sep="")
  if(is.null(property) | is.null(property_value)){
    clause <- paste(clause,')',sep='')
  } else {
    clause <- paste(clause,' {',property,": '",property_value,"'})",sep="")
  }
  return(clause)
}