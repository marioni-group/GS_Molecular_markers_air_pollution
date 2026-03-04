#Figure 1: Multi-panel plot 
#A = Geographical distribution for 2005-2011 for each pollutant
#B = Exceedance days for each of 4 pollutants for 2005 - 2011
#C = Exceedance over time -> use median not mean, try weekly instead of daily max/min?

library(data.table)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(dplyr)
library(terra)
library(ggplot2)


#_______________________________
#Figure 1: Panel A = Geographical distribution of pollutants (average between 2005 - 2011)

#file list of all pollutants
file_list <- list.files(path = "/Pollution/EMEP4uk", pattern = "*.nc", full.names = TRUE)
file_list <- file_list[6:12] #subset to 2005 - 2012

#list of pollutants
pollutants <- c("SURF_ug_NO2", "SURF_ug_NO", "SURF_ug_PM25_rh50", "SURF_ppb_O3", "SURF_ug_PM10_rh50", "SURF_ug_SO2", "SURF_ug_NO3_C", "SURF_ug_NO3_F")

#Load results
results_df <- fread("/Pollution/EMEP4uk/EMEP4UK_yearlysummary.csv")


results_df_tiles <- results_df %>%
  group_by(Pollutant, x, y) %>%
  summarise(
    overall_mean = mean(mean, na.rm = TRUE),
    n_years = n_distinct(Year[!is.na(Year)]),
    .groups = "drop"
  ) %>%
  group_by(Pollutant) %>%
  mutate(
    dx = median(diff(sort(unique(x))), na.rm = TRUE),
    dy = median(diff(sort(unique(y))), na.rm = TRUE)
  ) %>%
  ungroup()

#Create outline map for Scotland
# 1) Scotland admin-1 units, then dissolve to one polygon
scotland_sf <- ne_states(geounit = "Scotland", returnclass = "sf")
scot_outline_sf <- scotland_sf |>
  st_make_valid() |>
  summarise()


# 2) Make sure in matching crs
file <- file_list[1]
stackcheck <- terra::rast(file, subds = pollutants)


# 2) Match the raster CRS (stackcheck is your SpatRaster)
target_crs <- st_crs(crs(stackcheck, proj = TRUE))
scot_outline_sf <- st_transform(scot_outline_sf, target_crs) #the outline of the map should now be in the same coordinates as the model data
crs(scot_outline_sf)
target_crs

#make sure x is numeric
results_df_tiles$x <- as.numeric(results_df_tiles$x)

#create labels for the pollutants
lab_pollutant <- c(
  SURF_ug_NO2        = "NO2 (µg/m³)",
  SURF_ug_PM25_rh50  = "PM2.5 (µg/m³)",
  SURF_ug_SO2        = "SO2 (µg/m³)",
  SURF_ug_PM10_rh50  = "PM10 (µg/m³)",
  SURF_ug_NO         = "NO (µg/m³)",
  SURF_ppb_O3        = "O3 (ppb)",
  SURF_ug_NO3_C      = "NO3_C (µg/m³)",
  SURF_ug_NO3_F      = "NO3_F (µg/m³)"
)

library(ggnewscale)
library(grid)
library(patchwork)

df <- fread("/Pollution/EMEP4uk/daily_exceedances_bypoll_acrossscot.csv")
results2007 <- df %>% filter(Year == "2007" & Pollutant %in% c("SURF_ug_NO2", "SURF_ug_PM25_rh50", "SURF_ug_SO2", "SURF_ug_PM10_rh50"))
xlim_scot <- range(results2007$x, na.rm = TRUE)
ylim_scot <- range(results2007$y, na.rm = TRUE)
  x_cut <- quantile(results_df_tiles$x, 0.005, na.rm = TRUE)


make_plot_bottom <- function(dfi) {
  pol <- unique(dfi$Pollutant)
    pol_lab <- if (pol %in% names(lab_pollutant)) unname(lab_pollutant[pol]) else pol


dfi <- dfi |>
  dplyr::mutate(
    overall_mean = ifelse(
      x < x_cut, NA, overall_mean
    )
  )

  ggplot() +
    geom_tile(data = dfi, aes(x = x, y = y, fill = overall_mean)) +
    geom_sf(data = scot_outline_sf, fill = NA, color = "black", linewidth = 0.1) +
     coord_sf(crs = target_crs, default_crs = target_crs, 
     xlim = c(-650000, 45400),
     ylim = c(-3787791, -3053791),
     expand = FALSE, datum = NA) +
    scale_fill_viridis_c(
      name   = pol_lab,
      limits = range(dfi$overall_mean, na.rm = TRUE),
      option = "turbo",
      na.value = "transparent",
      guide = guide_colorbar(
        barheight = unit(8, "pt"),
        barwidth  = unit(60, "pt"),
        title.position = "top"
      )
    ) +
    theme(
      plot.margin = margin(10,2,8,2),
      legend.key.size = unit(1, "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.position      = "bottom",
      legend.box.margin    = margin(t = 2),
      legend.background    = element_rect(fill = "white", colour = "grey85"),
      legend.margin        = margin(4, 8, 6, 2)
    )
}

plots <- lapply(split(results_df_tiles, results_df_tiles$Pollutant), make_plot_bottom)

clean_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title       = element_blank(),
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.line        = element_blank()
)

panelafixed <- lapply(plots, function(p) {
 p <- p +
    theme(
      # remove panel background & gridlines
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),

      # remove plot background
      plot.background  = element_blank(),

      # remove axes
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),

      # legend at bottom, horizontal, no box
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(2, "mm"),
      legend.text       = element_text(size = 6),
      legend.title      = element_text(size = 6),

      # small margins
      plot.margin = margin(2, 2, 2, 2)
    )
  # Fix panel size
  egg::set_panel_size(p, width = unit(3.5, "cm"), height = unit(3.5, "cm"))

})


panela <- do.call(
  gridExtra::arrangeGrob,
  c(panelafixed, list(ncol = 4, nrow = 2))
)


png(
  filename = "/Pollution/EMEP4uk/panela_3.png",
  width = 18,
  height = 12,
  units = "cm",
  res = 300
)

grid::grid.draw(panela)
dev.off()

#_______________________________
#Figure 1: Panel B = Geographical distribution of exceedances of WHO air quality limits

#Load exceedances data 
results_df <- fread("/Pollution/EMEP4uk/daily_exceedances_bypoll_acrossscot.csv")

#Create summary exceedances plot. Given exceedances are defined based on a one-year period -> select one year to represent (2007)
#set 0s to transparent
results_df <- results_df %>% mutate(Exceedances = ifelse(Exceedances == 0, NA, Exceedances))

#check projection of raster
stackcheck <- terra::rast("/Pollution/EMEP4uk/EMEP4UK_2005_day.nc", subds = "SURF_ug_PM10_rh50")
terra::crs(stackcheck)


# 1) Scotland admin-1 units, then dissolve to one polygon
scotland_sf <- ne_states(geounit = "Scotland", returnclass = "sf")
scot_outline_sf <- scotland_sf %>%
  st_make_valid() %>%
  summarise()  # dissolve internal borders

# 2) Match the raster CRS (stackcheck is your SpatRaster)
target_crs <- st_crs(crs(stackcheck, proj = TRUE))
scot_outline_sf <- st_transform(scot_outline_sf, target_crs)
max_ex <- max(results_df$Exceedances, na.rm = TRUE)

library(scales) #trim scale so more easily interpretable
vals  <- results_df$Exceedances
q95   <- quantile(vals, 0.95, na.rm = TRUE) 

fill_scale <- scale_fill_viridis_c(
  option = "rocket", direction = -1,
  trans = "sqrt",
  limits = c(0, q95),
  na.value = "transparent",
  oob = scales::squish,
  name = "Exceedances"
)


# df should have x,y in the same CRS as stackcheck (e.g., from as.data.frame(stackcheck, xy = TRUE))
results2007 <- results_df %>% filter(Year == "2007" & Pollutant %in% c("SURF_ug_NO2", "SURF_ug_PM25_rh50", "SURF_ug_SO2", "SURF_ug_PM10_rh50"))

exceedanceplot <- ggplot() +
  geom_raster(data = results2007, aes(x = x, y = y, fill = Exceedances)) +
  geom_sf(data = scot_outline_sf, fill = NA, color = "black", linewidth = 0.1) +
  # scale_fill_viridis_c(option = "turbo", na.value = "transparent") +
  fill_scale + 
  coord_sf(crs = target_crs, default_crs = target_crs, expand = FALSE, datum = NA) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.1, "lines")
  ) +
  facet_grid(~Pollutant,
  labeller = labeller(
      Pollutant = as_labeller(c(
      SURF_ug_NO2  = "NO2",
        SURF_ug_PM25_rh50 = "PM2.5",
       SURF_ug_SO2 = "SO2",
        SURF_ug_PM10_rh50 = "PM10"
      ))), switch = "y") +
  labs(title = "",
       fill = "Exceedances") +
   theme(
      legend.key.size = unit(2, "mm"),
      legend.key.height = unit(4, "mm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.position      = "right",
      legend.background    = element_rect(fill = "white", colour = "white"),
      legend.margin        = margin(4, 4, 4, 4)
    )



#_______________________________
#Figure 1: #C = Exceedance over time -> use median 

df <- fread("/Pollution/EMEP4uk/daily_poll_across_scotland.csv")

df <- df %>%
  group_by(pollutant) %>%
  mutate(day2 = row_number()) %>%
  ungroup()


#all together
df <- df %>% mutate(year2 = as.Date(paste0(year, "-01-01")))
df <- df %>% mutate(year2 = lubridate::year(year2))
year_positions <- df %>% group_by(year2) %>% summarise(pos = min(day2)) %>% ungroup()
df$pollutant <- as.factor(df$pollutant)
df <- df %>% rename(poll = pollutant)
df <- df %>%
  mutate(limit = ifelse(poll == "SURF_ug_PM25_rh50", 15,
                  ifelse(poll == "SURF_ug_PM10_rh50", 45,
                    ifelse(poll == "SURF_ug_NO2", 25, 40))))


#match colours to other plots
clist <- paletteer_d("MetBrewer::Tam")
clst <- c("#FFB242FF", "#9F2D55FF", "#62205FFF", "#341648FF")


# Fit smooths and extract predicted values

df_smooth <- df %>%
  group_by(poll) %>%                     # smooth separately per facet
  mutate(
    max_s = predict(loess(max ~ day2, span = 0.01)),
    min_s = predict(loess(min ~ day2, span = 0.01))
  )

combplot <- df_smooth %>%
ggplot(aes(x = day2, y = median)) +
geom_smooth(aes(color = poll, fill = poll), method = "loess", span = 0.05, se = FALSE) +
geom_line(aes(y = max_s, color = poll), size = 0.3) +
geom_line(aes(y = min_s, color = poll), size = 0.3) +
  geom_ribbon(aes(ymin = min_s, ymax = max_s, fill = poll),
              alpha = 0.25) +
geom_hline(aes(yintercept = limit), linetype = "dashed", colour = "black") +
scale_fill_manual(values = clst) +
scale_colour_manual(values = clst) + 
scale_x_continuous(breaks = year_positions$pos, labels = year_positions$year2) + 
labs(x = "Year", y = "Median & range of pollutant") +
theme_minimal() +
theme(legend.position = "none") +
facet_wrap(~poll, ncol = 1, scales = "free_y", strip.position = "left",
      labeller = labeller(
      poll = as_labeller(c(
      SURF_ug_NO2  = "NO2",
        SURF_ug_PM25_rh50 = "PM2.5",
       SURF_ug_SO2 = "SO2",
        SURF_ug_PM10_rh50 = "PM10"
      )))) +
theme(
  strip.placement = "outside", # move the strip outside the axis
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 4)
)


all <- gridExtra::arrangeGrob(
  panela,
  exceedanceplot,
  combplot,
  ncol = 1,
  heights = c(10,4,8) # relative panel heights
) 


png("/Pollution/EMEP4uk/multi_panel_fig1.png", height = 25, width = 15, units = "cm", res = 300)
grid.newpage()
grid.draw(all)

grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold")) 
grid.text("B", x = unit(0.02, "npc"), y = unit(0.52, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold")) 
grid.text("C", x = unit(0.02, "npc"), y = unit(0.36, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()




