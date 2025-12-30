library(readxl)
library(dplyr)
library(writexl)
options(scipen = 999)
# 讀取原始檔
io35_raw <- read_excel("hw35.xlsx")

#取出 rownames 與部門欄位
rownames_35 <- io35_raw[[1]]
dep_cols <- names(io35_raw)[2:36]

#整理 35×35 矩陣
io35_mat <- io35_raw %>%
  select(all_of(dep_cols)) %>%
  mutate(across(everything(), ~as.numeric(gsub(",", "", .)))) %>%
  as.matrix()
rownames(io35_mat) <- rownames_35

# 35→12 mapping
map12 <- c(
  5, 1, 2, 6, 7, 7, 3, 6, 6, 6, 6, 8, 8, 8, 7, 7, 4, 8, 8, 8, 7, 9, 9, 10, 11,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12
)

# 12 部門名稱
new_names12 <- c(
  "天然氣", "煤", "油品", "電力及蒸汽",
  "農林漁牧業", "礦業與基礎材料製造業", "一般製造業",
  "電子與精密製造業", "建築與公共工程部門",
  "批發與零售業", "運輸部門", "其他服務業"
)

#建立 12×(12 + 10 + 其他列) 矩陣
#前 12 = 投入產出，後面每列對應原始列的指標
extra_cols <- c("民間消費","政府消費","固定資本形成","存貨變動",
                "出口","進口","C+I+G+X-M","產值AK+AT","產值","平衡狀態")
other_rows <- rownames_35[36:length(rownames_35)] # 中間投入、勞動報酬…產值
io12_all <- matrix(0, nrow=12 + length(other_rows), ncol=12 + length(extra_cols))
rownames(io12_all) <- c(new_names12, other_rows)
colnames(io12_all) <- c(new_names12, extra_cols)

#整併 12×12 投入產出
for(i in 1:35){
  for(j in 1:35){
    io12_all[map12[i], map12[j]] <- io12_all[map12[i], map12[j]] + io35_mat[i,j]
  }
}

# 整併後面欄位
io35_extra <- io35_raw %>%
  select(all_of(extra_cols)) %>%
  mutate(across(everything(), ~as.numeric(gsub(",", "", .)))) %>%
  as.matrix()

for(i in 1:35){
  for(j in 1:length(extra_cols)){
    io12_all[map12[i], 12 + j] <- io12_all[map12[i], 12 + j] + io35_extra[i,j]
  }
}

#將原始列的經濟指標列依 12 部門加總
for(r in 36:length(rownames_35)){
  for(c in 1:35){
    io12_all[12 + (r-35), map12[c]] <- io12_all[12 + (r-35), map12[c]] + as.numeric(io35_mat[r,c])
  }
}
io12_all<-io12_all[, -ncol(io12_all)]
io12_all <- cbind(投入 = rownames(io12_all), io12_all)
io12_all <- as.data.frame(io12_all)

# 只轉換數值欄位 (排除第一欄)
io12_all <- io12_all %>%
  mutate(across(-投入, as.numeric))

# 輸出 Excel
#write_xlsx(io12_all, "io12_complete.xlsx")

#####################

energy35 <- read_excel("energy_data.xlsx")
energy_cols <- names(energy35)[2:7]

# 數值轉 numeric
energy35_mat <- energy35 %>%
  select(all_of(energy_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

#----------------------------------------
# 建立 12 × 6 能源表
#----------------------------------------
energy12 <- matrix(0, nrow=12, ncol=length(energy_cols))
rownames(energy12) <- new_names12
colnames(energy12) <- energy_cols

#----------------------------------------
# 依 mapping 加總
#----------------------------------------
for(i in 1:35){
  for(j in 1:length(energy_cols)){
    energy12[map12[i], j] <- energy12[map12[i], j] + energy35_mat[i, j]
  }
}

# 轉 data frame
energy12_df <- as.data.frame(energy12) %>%
  mutate(部門 = rownames(energy12)) %>%
  select(部門, everything())

#----------------------------------------
# 輸出
#----------------------------------------
# write_xlsx(energy12_df, "energy12.xlsx")

energy12_df
###########################
energy_co2_35 <- read_excel("energy_co2_data.xlsx")
co2_cols <- names(energy_co2_35)[2:8]

# 數值全部轉 numeric
co2_35_mat <- energy_co2_35 %>%
  select(all_of(co2_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# 建立 12 × 7 (含總排放) 的矩陣
co2_12 <- matrix(0, nrow = 12, ncol = length(co2_cols))
rownames(co2_12) <- new_names12
colnames(co2_12) <- co2_cols

# 用 map12 加總 35→12
for(i in 1:35){
  for(j in 1:length(co2_cols)){
    co2_12[ map12[i], j ] <- co2_12[ map12[i], j ] + co2_35_mat[i, j]
  }
}

# 轉成 data frame
co2_12_df <- as.data.frame(co2_12) %>%
  mutate(部門 = rownames(co2_12)) %>%
  select(部門, everything())

# 檢視結果
co2_12_df

# 輸出
# write_xlsx(co2_12_df, "energy_co2_12.xlsx")

