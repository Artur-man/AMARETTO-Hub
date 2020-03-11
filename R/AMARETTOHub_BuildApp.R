#' AMARETTOHub_BuildApp
#' 
#' Building the AMARETTOHub shiny+rmarkdown app 
#'
#' @param AppTitle The title of the App
#' @param AMARETTODirectory The AMARETTO Directory
#' @param con_info The list of necessary Neo4j server information: url, username and password 
#' @param Community_Info A list of information on Community AMARETTO: (1) a vector of cohort names and (2) type of nodes in Community AMARETTO
#'
AMARETTOHub_BuildApp <- function(AppTitle, AMARETTODirectory, con_info, Community_Info){
  
  # AMARETTO-Hub run directory
  if(!dir.exists('run/'))
    dir.create('run/')
  
  # copy templates
  file.copy(system.file('templates/index.Rmd', package = 'AMARETTOHub'),'index.Rmd')
  file.copy(system.file('templates/navbar.html', package = 'AMARETTOHub'),'run/navbar.html')
  file.copy(system.file('templates/tag.html', package = 'AMARETTOHub'), 'run/tag.html')

  # prepare index.Rmd
  index_file  <- readLines("index.Rmd")
  index_file <- gsub('con_url', paste0('"',con_info$url,'"'),index_file)
  index_file <- gsub('con_user', paste0('"',con_info$user,'"'), index_file)
  index_file <- gsub('con_password', paste0('"',con_info$password,'"'), index_file)
  
  # enter title
  index_file <- gsub('Enter Title Here', AppTitle, index_file)
  writeLines(index_file, con="index.Rmd")
}