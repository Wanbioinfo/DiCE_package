# ----------- Functions to calculate different network centralities ------------


#' Run calc_betweenness: Betweenness calculation function on undirected graph based on the distance
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their betweeness values 
#' @noRd
calc_betweenness <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = (1-(interactions$weight))+0.001)
  E(graph)$weight[is.na(E(graph)$weight)] <- 1
  
  btw_vals <- betweenness(graph,
                          normalized = TRUE)
  
  btw_df <- data.frame(value = btw_vals)
  colnames(btw_df)[1] <- "Betweenness"
  
  return(btw_df)
}

#' Run calc_eigenVC: Eigen vector centrality calculation function on undirected graph based on the |correlation coeff|
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their eigen-vector centrality values 
#' @noRd
calc_eigenVC <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = interactions$weight)
  
  # Replace missing edge weights with 0 (min weight -weakest link)
  E(graph)$weight[is.na(E(graph)$weight)] <- 0
  
  eig_vals <- eigen_centrality(graph)$vector
  
  eig_df <- data.frame(value = eig_vals)
  colnames(eig_df)[1] <- "EigenVector"
  
  return(eig_df)
}

#' Run calc_authority: Authority calculation function on undirected graph based on the |correlation coeff|
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their authority values 
#' @noRd
calc_authority <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = interactions$weight)
  
  # Replace missing edge weights with 0 (min weight -weakest link)
  E(graph)$weight[is.na(E(graph)$weight)] <- 0
  
  auth_vals <- authority_score(graph)$vector
  
  auth_df <- data.frame(value = auth_vals)
  colnames(auth_df)[1] <- "Authority"

  return(auth_df)
}

#' Run calc_strength: Strength calculation function on undirected graph based on the |correlation coeff|
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their strength values 
#' @noRd
calc_strength <- function(interactions){

  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = interactions$weight)
  
  # Replace missing edge weights with 0 (min weight -weakest link)
  E(graph)$weight[is.na(E(graph)$weight)] <- 0
  
  strength_vals <- strength(graph)
  
  strength_df <- data.frame(value = strength_vals)
  colnames(strength_df)[1] <- "Strength"
  
  return(strength_df)
}

#' Run calc_closeness: Closeness calculation function on undirected graph based on the distance
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their closeness values 
#' @noRd
calc_closeness <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = (1-(interactions$weight))+0.001)
  E(graph)$weight[is.na(E(graph)$weight)] <- 1
  
  close_vals <- closeness(graph,
                          normalized =TRUE)
  
  close_df <- data.frame(value = close_vals)
  colnames(close_df)[1] <- "Closeness"
  
  return(close_df)
}

#' Run calc_pagerank: pagerank calculation function on undirected graph based on the |correlation coeff|
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their pagerank values 
#' @noRd
calc_pagerank <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = interactions$weight)
  
  # Replace missing edge weights with 0 (min weight -weakest link)
  E(graph)$weight[is.na(E(graph)$weight)] <- 0
  
  pgr_vals <- page_rank(graph)$vector
  
  pgr_df <- data.frame(value = pgr_vals)
  colnames(pgr_df)[1] <- "PageRank"
  
  return(pgr_df)
}

#' Run calc_harmonic: Harmonic calculation function on undirected graph based on the distance
#' Helper function (not for users)
#'
#' @param interactions Dataframe of the protein-protein interaction pairs with corresponding edge weight
#'
#' @return Dataframe of the genes and their harmonic values 
#' @noRd
calc_harmonic <- function(interactions){
  
  graph <- graph_from_data_frame(interactions,directed = FALSE)
  graph <- set_edge_attr(graph,
                         "weight",
                         value = (1-(interactions$weight))+0.001)
  E(graph)$weight[is.na(E(graph)$weight)] <- 1
  
  harmonic_vals <- harmonic_centrality(graph,
                                       normalized = TRUE)
  
  harmonic_df <- data.frame(value = harmonic_vals)
  colnames(harmonic_df)[1] <- "Harmonic"
  
  return(harmonic_df)
  
}

