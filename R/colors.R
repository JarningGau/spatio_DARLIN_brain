array.colors = c(
  "CA" = "#3273A0",
  "TA" = "#E0802C",
  "RA" = "#3A913A"
)

domain.colors = c(
  "Crypt" = "#FCBA7E",
  "Villus" = "#BFAED4",
  "Lamina" = "#3B6DB0",
  "Muscle" = "#BE5812",
  "Meso" = "#7AC57A",
  "Adipose" = "#BCBA2A",
  "ILF" = "#DA292E"
)

# celltype.colors = c(
#   "Paneth" = "#1f77b4",         # 深蓝
#   "ISC" = "#0c3b99",            # 靛蓝
#   "TA" = "#1ca9c9",             # 湖蓝
#   "Ent (bottom)" = "#2ca02c",   # 草绿色
#   "Ent (mid)" = "#6ab344",      # 浅绿
#   "Ent (tip)" = "#b2df8a",      # 青绿色
#   "Ent (ISG+)" = "#154734",     # 深绿
#   "EC" = "#d4b9da",             # 粉紫
#   "EEC" = "#984ea3",            # 紫罗兰
#   "Goblet" = "#ffdd57",         # 柠檬黄
#   "Tuft" = "#00ced1",           # 暗青色
#   "Plasma" = "#e377c2",         # 粉红
#   "T cell" = "#ff7f0e",         # 橘色
#   "Macrophage" = "#bcbd22",     # 暗黄绿
#   "Endo" = "#8c564b",           # 褐色
#   "Fibro" = "#c49c94",          # 浅褐
#   "RBC" = "#d62728",            # 红色
#   "Myocyte" = "#f781bf",        # 粉紫
#   "Neuron" = "#7f7f7f",         # 中性灰
#   "Adipocyte" = "#c7c7c7",      # 浅灰
#   "ILF" = "#17becf",            # 青蓝
#   "Meso" = "#8c6d31"            # 土褐
# )

celltype.colors = c(
  # Crypt (蓝色系)
  "Paneth" = "#a6cee3",   # 浅蓝
  "ISC" = "#0056A4",      # 深蓝
  # Villus Ent (绿色系)
  "TA" = "#2ECC40",       # 浅绿
  "Ent (bottom)" = "#28A745", # 绿
  "Ent (mid)" = "#218838",    # 中绿
  "Ent (tip)" = "#1C7430",    # 深绿
  "Ent (ISG+)" = "#155724",   # 更深绿
  # Villus Others (独特颜色)
  "EC" = "#D8BFD8",       # 浅紫色
  "EEC" = "#9400D3",      # 深紫
  "Goblet" = "#FFD700",   # 金色
  "Tuft" = "#00FFFF",     # 青色
  # Lamina (独特颜色)
  "Plasma" = "#FF1493",   # 深粉红
  "T cell" = "#FF8C00",   # 深橙
  "Macrophage" = "#8B008B",# 深紫
  "Endo" = "#00BFFF",     # 深天蓝
  "Fibro" = "#32CD32",    # 石灰绿
  "RBC" = "#FF4500",      # 红橙
  # Muscle (紫色系)
  "Myocyte" = "#ffffb3",  # 紫色
  # Neuron (青绿色系，独特颜色)
  "Neuron" = "blue",   # 青绿色，明亮且独特
  # Adipose (灰色系)
  "Adipocyte" = "#A9A9A9",# 灰色
  # ILF (青色系)
  "ILF" = "#e31a1c",      # 青色
  # Meso (棕色系)
  "Meso" = "#8B4513"      # 棕色
)

get_colors = function(n) {
  rep(ggsci::pal_d3("category20")(20), 20)[1:n]
}

# color_list = [
#   "#1f78b4", "#ff7f00", "#33a02c", "#e31a1c", "#6a3d9a",
#   "#8c564b", "#e377c2", "#969696", "#17becf", "#253494",
#   "#ffcfcf", "#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6",
#   "#99d8c9", "#8c96c6", "#66c2a5", "#636363", "#ffd92f",
#   "#b3b3b3", "#db5f57"
# ]
# 
# color_list_2 = [
#   "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
#   "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
#   "#ffed6f", "#ffffb3", "#99d8c9", "#8c96c6", "#66c2a5",
#   "#636363", "#ffd92f", "#b3b3b3", "#db5f57"
# ]