library(deSolve)
library(ggplot2)
library(reshape2)

# Define the total population
N <- 5070000     # Total population
N_u <- 0.07 * N  # Total population of unvaccinated
N_v <- 0.93 * N  # Total population of vaccinated


# Initial conditions
S_u <- N_u - 11107  # Susceptible unvaccinated population
E_u <- 4031         # Exposed (infected) unvaccinated population
I_u <- 7076         # Infectious unvaccinated population
R_u <- 0            # Recovered (removed) unvaccinated population
S_v <- N_v - 2469   # Susceptible vaccinated population
E_v <- 700          # Exposed (infected) vaccinated population
I_v <- 1769         # Infectious vaccinated population
R_v <- 0            # Recovered (removed) vaccinated population


# Contact rates and physical distancing
c_u <- 2.11    # Average contacts per day for unvaccinated
c_v <- 2.2     # Average contacts per day for vaccinated
l_u <- 0.45   # Level of physical distancing for unvaccinated
l_v <- 0.134    # Level of physical distancing for vaccinated
V_e <- 0.8     # Vaccine efficacy against infection


# Contact matrix
adjust_contact_proportions <- function(p_vv, p_uu) {
  p_vu <- 1 - p_vv
  p_uv <- 1 - p_uu
  return(list(p_vu = p_vu, p_vv = p_vv, p_uv = p_uv, p_uu = p_uu))
}

contact_proportions <- adjust_contact_proportions(0.72, 0.72)  # p_vv $ p_uu
p_vv <- contact_proportions$p_vv  
p_vu <- contact_proportions$p_vu  
p_uv <- contact_proportions$p_uv  
p_uu <- contact_proportions$p_uu 

p_vv
p_vu
p_uv
p_uu


# Function to plot the contact matrix
plot_contact_matrix <- function(p_vu, p_uu, p_vv, p_uv) {
  # Create the contact matrix
  contact_matrix <- matrix(c(p_vu, p_uu, p_vv, p_uv), nrow = 2, ncol = 2, byrow = TRUE)
  colnames(contact_matrix) <- c("Vaccinated", "Unvaccinated")
  rownames(contact_matrix) <- c("Unvaccinated", "Vaccinated")
  
  # Melt the matrix to a data frame for ggplot
  contact_df <- melt(contact_matrix)
  
  # Create the plot
  ggplot(contact_df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "#FFC0CB", high = "red") +
    labs(title = "Low Homophily in Vaccinated and Low Homophily in Unvaccinated",
         x = "Recipient",
         y = "Sender",
         fill = "Proportion") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
}

# Call the function to plot the contact matrix
plot_contact_matrix(p_vu, p_uu, p_vv, p_uv)


# Model parameters
alpha <- 0           # Birth rate
beta <- 0.6          # Transmission rate from vaccinated to unvaccinated
sigma <- 1/2         # Incubation rate
gamma <- 1/7         # Recovery rate
omega <- 0           # Vaccination rate
mu <- 0              # Mortality rate
delta <- 0           # Disease-induced mortality
kappa_u <- 1/60      # Waning immunity for unvaccinated group
kappa_v <- 1/180     # Waning immunity for vaccinated group

# Define model parameters
parms <- list(
  alpha = alpha,
  beta = beta,
  sigma = sigma,
  gamma = gamma,
  omega = omega,
  mu = mu,
  delta = delta,
  kappa_u = kappa_u,
  kappa_v = kappa_v,
  c_u = c_u,
  c_v = c_v,
  p_uv = p_uv,
  p_vu = p_vu,
  p_uu = p_uu,
  p_vv = p_vv,
  l_u = l_u,
  l_v = l_v,
  V_e = V_e,
  N_u = N_u,
  N_v = N_v
)

# Define the initial conditions vector
y <- c(S_u = S_u, E_u = E_u, I_u = I_u, R_u = R_u, S_v = S_v, E_v = E_v, I_v = I_v, R_v = R_v)

# SEIR model function
SEIR_model <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    # Recalculate lambda values at each time step
    lambda_uv <- p_vu * c_v * beta * (1 - V_e) * (1 - l_u) * (I_u / N_u)
    lambda_vu <- p_uv * c_u * beta * (1 - l_v) * (I_v / N_v)
    lambda_uu <- p_uu * c_u * beta * (1 - l_u) * (I_u / N_u)
    lambda_vv <- p_vv * c_v * beta * (1 - V_e) * (1 - l_v) * (I_v / N_v)
    
    # Force of infection
    lambda_u <- lambda_uv + lambda_uu
    lambda_v <- lambda_vu + lambda_vv
    
    N_v <- S_v + E_v + I_v + R_v
    N_u <- S_u + E_u + I_u + R_u
    
    # Differential equations
    dS_u <- alpha * N - lambda_u * S_u - mu * S_u - omega * S_u + kappa_u * R_u
    dE_u <- lambda_u * S_u - sigma * E_u - mu * E_u
    dI_u <- sigma * E_u - mu * I_u - delta * I_u - gamma * I_u
    dR_u <- gamma * I_u - mu * R_u - kappa_u * R_u
    
    dS_v <- omega * S_u - lambda_v * S_v - mu * S_v + kappa_v * R_v
    dE_v <- lambda_v * S_v - sigma * E_v - mu * E_v
    dI_v <- sigma * E_v - mu * I_v - delta * I_v - gamma * I_v
    dR_v <- gamma * I_v - mu * R_v - kappa_v * R_v
    
    return(list(c(dS_u, dE_u, dI_u, dR_u, dS_v, dE_v, dI_v, dR_v)))
  })
}

# Time points for the output
time.out <- seq(0, 2 * 90, 1)

# Solving the ODE
ts_SEIR_model <- data.frame(lsoda(
  y = y,
  times = time.out,
  func = SEIR_model,
  parms = parms
))

# Add the Sum column
ts_SEIR_model$Sum <- rowSums(ts_SEIR_model[, c("S_u", "E_u", "I_u", "R_u", "S_v", "E_v", "I_v", "R_v")])

# Display the updated data frame
print(head(ts_SEIR_model))
print(tail(ts_SEIR_model))


# Plots -------------------------------------------------------------------


# # Combined plot of all the state for Unvaccinated and Vaccinated groups
# ggplot(ts_SEIR_model, aes(x = time.out)) +
#   geom_line(aes(y = S_u, color = "S_u")) +
#   geom_line(aes(y = E_u, color = "E_u")) +
#   geom_line(aes(y = I_u, color = "I_u")) +
#   geom_line(aes(y = R_u, color = "R_u")) +
#   geom_line(aes(y = S_v, color = "S_v")) +
#   geom_line(aes(y = E_v, color = "E_v")) +
#   geom_line(aes(y = I_v, color = "I_v")) +
#   geom_line(aes(y = R_v, color = "R_v")) +
#   xlab("Time (days)") +
#   ylab("Number of Cases") +
#   ggtitle("COVID-19 in British Columbia") +
#   scale_color_manual(values = c("S_u" = "darkgreen", "E_u" = "orange", "I_u" = "black", "R_u" = "pink",
#                             	"S_v" = "darkblue", "E_v" = "grey", "I_v" = "red", "R_v" = "brown")) +
#   theme_classic()+theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# 

# Create the plot for Vaccinated populationHigh Homophily in Vaccinated vs. Low Homophily in Unvaccinated with Varying Compliance Levels

# ggplot(ts_SEIR_model, aes(x = time.out))+
#   geom_line(aes(y = S_v, color = "S_v")) +
#   geom_line(aes(y = E_v, color = "E_v")) +
#   geom_line(aes(y = I_v/N_v, color = "I_v")) +
#   geom_line(aes(y = R_v, color = "R_v")) +
#   xlab("Time (days)") +
#   ylab("Population") +
#   ggtitle("COVID-19 in British Columbia") +
#   scale_color_manual(values = c("S_v" = "darkblue", "E_v" = "grey", "I_v" = "red", "R_v" = "brown")) +
#   theme_classic()+ theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# Create the plot For infections
ggplot(ts_SEIR_model, aes(x = time.out)) +
  geom_line(aes(y = I_u/N_u, color = "I_u")) +
  geom_line(aes(y = I_v/N_v, color = "I_v")) +
  geom_line(aes(y = I_u/N_u + I_v/N_v, color = "Total"))+
  xlab("Time (days)") +
  ylab("Proportion of Cases") +
  ggtitle("High Homophily in Vaccinated vs. High Homophily in Unvaccinated Varying Compliance Levels") +
  scale_color_manual(values = c("I_v" = "orange", "I_u" = "red","Total" = "blue")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.01, face = "bold"))


