##### NET RECLASSIFICATIONS IMROVEMENT INDEX Function ####
# Returns vector with continuous NRI, event-NRI and non-event-NRI from sample data

NRIcont <- function (data, indices) {
  data <- data[indices,]
  n.event <- data %>%
    filter(complications.dicho.fct.Severe == 1) %>%
    nrow()
  
  n.nonevent <- data %>%
    filter(complications.dicho.fct.Severe == 0) %>%
    nrow()
  
  n_up.event <- data %>%
    filter(complications.dicho.fct.Severe == 1) %>%
    filter(up) %>%
    nrow()
  
  n_down.event <- data %>%
    filter(complications.dicho.fct.Severe == 1) %>%
    filter(down) %>%
    nrow()
  
  n_up.nonevent <- data %>%
    filter(complications.dicho.fct.Severe == 0) %>%
    filter(up) %>%
    nrow()
  
  n_down.nonevent <- data %>%
    filter(complications.dicho.fct.Severe == 0) %>%
    filter(down) %>%
    nrow()
  
  eventNRI = (n_up.event - n_down.event) / n.event
  noneventNRI = (n_down.nonevent - n_up.nonevent) / n.nonevent
  NRI = eventNRI + noneventNRI
  
  return(c(NRI, eventNRI, noneventNRI))
}