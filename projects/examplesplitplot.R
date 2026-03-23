# Simulate a split plot design, with two factors: A (whole plot factor) and B (subplot factor)
# Four replicates of each
# Factor A (tillage): has two levels: CT and NT
# Factor B (fertilizer): has three levels: 0, 140, and 280 kg/ha of N

# Set seed for reproducibility
set.seed(123)
# Define factors and levels
tillage <- factor(c("CT", "NT"))
fertilizer <- factor(c("0", "140", "280"))
# Number of replicates
replicates <- 4
# Create a data frame to hold the simulated data
data <- expand.grid(
  Tillage = tillage,
  Fertilizer = fertilizer,
  Replicate = 1:replicates
)
# Simulate response variable (e.g., crop yield) with some random noise. Assume some effects of tillage and fertilizer.
data$Yield <- with(
  data,
  50 +
    ifelse(Tillage == "CT", 10, 0) + # CT has a positive effect on yield
    as.numeric(as.character(Fertilizer)) * 0.2 + # Fertilizer increases yield
    rnorm(nrow(data), mean = 0, sd = 5) # Random noise
)

# View the first few rows of the simulated data
head(data)
# plot the data to visualize the effects of tillage and fertilizer on yield
# Add a line showing mu (main mean)
# x axis should be main plot
# Within each group (main plot) colors represent different subplots (fertilizer levels)
# Add horizontal lines representing the group means at both levels. And lines from the global mean to the group means. and fro mthe group means to the subgroup means:

# Packages
library(ggplot2)
library(dplyr)

# prepare positions
tillage_centers <- tibble::tibble(
  Tillage = factor(c("CT", "NT"), levels = levels(tillage)),
  center = c(1, 3)
)

fert_offsets <- tibble::tibble(
  Fertilizer = factor(c("0", "140", "280"), levels = levels(fertilizer)),
  offset = c(-0.32, 0, 0.32) # adjust spacing here
)

plot_data <- data %>%
  left_join(tillage_centers, by = "Tillage") %>%
  left_join(fert_offsets, by = "Fertilizer") %>%
  mutate(x = center + offset)

# means
global_mean <- mean(plot_data$Yield)
group_means <- plot_data %>%
  group_by(Tillage, center) %>%
  summarize(
    till_mean = mean(Yield),
    xmin = min(x) - 0.08,
    xmax = max(x) + 0.08,
    .groups = "drop"
  )

subplot_means <- plot_data %>%
  group_by(Tillage, Fertilizer, x) %>%
  summarize(sub_mean = mean(Yield), .groups = "drop") %>%
  left_join(group_means %>% select(Tillage, till_mean), by = "Tillage")


# half-width for subgroup horizontal ticks
tick_half_width <- 0.12

p <- ggplot(plot_data, aes(x = x, y = Yield, color = Fertilizer)) +
  geom_jitter(
    aes(group = Fertilizer),
    width = 0.06,
    height = 0,
    size = 2,
    alpha = 0.9
  ) +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "black") +
  # tillage mean segments (span only subgroup region)
  geom_segment(
    data = group_means,
    aes(x = xmin, xend = xmax, y = till_mean, yend = till_mean),
    color = "blue",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  # global -> tillage mean
  geom_segment(
    data = group_means,
    aes(x = center, xend = center, y = global_mean, yend = till_mean),
    color = "blue",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  # tillage mean -> subgroup mean (vertical)
  geom_segment(
    data = subplot_means,
    aes(x = x, xend = x, y = till_mean, yend = sub_mean),
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  # subgroup horizontal ticks (centered at x, length = 2 * tick_half_width)
  geom_segment(
    data = subplot_means,
    aes(
      x = x - tick_half_width,
      xend = x + tick_half_width,
      y = sub_mean,
      yend = sub_mean
    ),
    color = "red",
    linewidth = 0.9
  ) +
  # optional point marker at subgroup mean
  geom_point(
    data = subplot_means,
    aes(x = x, y = sub_mean),
    color = "red",
    size = 2,
    shape = 18
  ) +
  scale_x_continuous(
    breaks = tillage_centers$center,
    labels = tillage_centers$Tillage,
    limits = c(0, 4)
  ) +
  scale_color_manual(
    values = c("0" = "green", "140" = "orange", "280" = "red"),
    name = "Fertilizer (kg/ha)"
  ) +
  labs(
    title = "Simulated Split-Plot Design",
    x = "Tillage",
    y = "Yield (kg/ha)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank())

p


# run anova using the following:
# fit.lm <- fxn(y ~ Tillage * Fertilizer + (1|Block/WUerror))
# Need to add WUerror and Block to the data frame first

data$WU <- interaction(data$Tillage, data$Replicate) # whole plot error term
data$Block <- data$Replicate # block is the same as replicate in this case

# Fit the mixed model using lme4
library(lme4)
model <- lmer(Yield ~ Tillage * Fertilizer + (1 | Block / WU), data = data)
# Get ANOVA table
library(car)
anova(model)
