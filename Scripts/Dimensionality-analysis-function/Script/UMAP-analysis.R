#####  UMAP (Uniform Manifold Approximation and Projection) analysis #####
# Uniform Manifold Approximation and Projection. 
# 假设可用数据样本均匀分布在拓扑空间中，可以从这些有限数据样本中近似并映射到低维空间。 
library(tidyverse)
library(Rtsne)
library(openxlsx)
library(umap)
library(palmerpenguins)



penguins <- penguins%>%
  drop_na() %>%
  select(-year) %>%
  mutate(ID=row_number())

penguins_meta <- penguins %>%
  select(ID,species,island,sex)

umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>%
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta,by="ID")

## 可视化 ggplot
umap_df %>% ggplot(aes(x = UMAP1,
                       y = UMAP2,
                       color = species,
                       shape =sex))+
  geom_point()+
  labs(x="UMAP1",
       y="UMAP2",
       subtitle = "UMAP plot")





