#' RemoveEntities
#' 
#' This function deletes entities of the Neo4j server
#'
#' @param label Label of the Entity to be removed 
#' @param bulk The number of nodes to be deleted at a time
#' @param con_info The list of necessary Neo4j server information: url, username and password 
#'
RemoveEntities <- function(label, bulk, con_info) {
  
  # establish communication with neo4j server
  con <- neo4r::neo4j_api$new(
    url = con_info$url,
    user = con_info$user,
    password = con_info$password)
  
  # checking existing entities
  listofEntities <- c('Gene','Module','Community','GeneSet','DriverPert','Phenotype','TransFactor')
  if(!label %in% listofEntities)
    stop('You have to choose a valid Entity label')
  
  # remove bulks
  paste0('MATCH (n:', label,') RETURN count(n) as count') %>% neo4r::call_neo4j(con) -> cn
  cn <- ceiling(cn$count$value/bulk)
  
  # start removing entities
  for(i in 1:cn){
    paste0('MATCH (r:', label, ') with r LIMIT ', bulk,' detach DELETE r') %>% neo4r::call_neo4j(con)
  }
  
}