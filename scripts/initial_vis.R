##
##### Initial visualisation of data #####
##
ggplot(data = data, aes(x = diff_sunrise)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Activity onset - sunrise time", y = "Count")


ggplot(data = data, aes(x = diff_sunrise, y = factor(area))) +
  geom_density_ridges(scale = 0.8, 
                      alpha = 0.8, 
                      aes(fill = area)) +
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB", "#5AAE61")) +
  labs(x = "Activity onset - sunrise time", y = "Count")


ggplot(data = data, aes(x = day_before_hatch, y = date_aprildays, group = box)) +
  geom_point() +
  geom_line() +
  facet_grid(year~area) +
  theme_bw()

ggplot(data = data, aes(x = day_before_hatch, y = mean_days_april, group = box)) +
  geom_point() +
  geom_line() +
  facet_grid(year~area) +
  theme_bw()
