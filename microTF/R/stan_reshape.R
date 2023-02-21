
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#' @importFrom stringr str_extract
#' @importFrom dplyr mutate select
#' @importFrom posterior summarise_draws
#' @export
posterior_means <- function(draws, index_names = c("V", "K")) {
  summarise_draws(draws, "mean") %>%
    mutate(index = str_extract(variable, "[0-9,]+")) %>%
    separate(index, index_names, convert = TRUE) %>%
    select(-variable)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select any_of
#' @importFrom purrr map
#' @importFrom tidyr pivot_wider
pivot_list <- function(x, var1, var2) {
  x %>%
    split(.$P, ~ pivot_wider(names_from = any_of(var1), values_from = mean)) %>%
    map(~ select(., -P)) %>%
    map(~ pivot_wider(., names_from = any_of(var2), values_from = mean)) %>%
    map(~ select(., -any_of(var1))) %>%
    map(~ as.matrix(.))
}

#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer separate pivot_wider
#' @importFrom dplyr mutate select
#' @importFrom purrr map
reshape_forecast <- function(draws) {
  draws %>%
    as_tibble() %>%
    rownames_to_column("draw") %>%
    pivot_longer(-draw, names_to = "variable") %>% 
    mutate(index = str_extract(variable, "[0-9,]+")) %>%
    separate(index, c("K", "n"), convert = TRUE) %>%
    select(-variable) %>%
    split(.$draw) %>%
    map(~ pivot_wider(., names_from = "n", values_from = "value")) %>%
    map(~ as.matrix(select(., -K:-draw)))
}

#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom posterior summarise_draws
summarize_posterior <- function(params) {
  L <- posterior_means(params$draws("L"), c("K", "V")) %>%
    pivot_wider(names_from = V, values_from = mean) %>%
    arrange(K) %>%
    column_to_rownames("K")
  A <- posterior_means(params$draws("A"), c("P", "K1", "K2")) %>%
    pivot_list("K1", "K2")
  B <- posterior_means(params$draws("B"), c("P", "K", "D")) %>%
    pivot_list("K", "D")
  sigma_e <- summarise_draws(params$draws("sigma_e"), "mean") %>%
    pull(mean)
  sigma_z <- summarise_draws(params$draws("sigma_z"), "mean") %>%
    pull(mean)
  
  list(L = L, A = A, B = B, sigma_e = sigma_e, sigma_z = sigma_z, K = nrow(A[[1]]))
}