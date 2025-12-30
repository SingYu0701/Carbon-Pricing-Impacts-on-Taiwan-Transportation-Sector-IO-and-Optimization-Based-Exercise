library(readxl)
library(dplyr)
library(Matrix)
library(tidyr)
library(writexl)
options(scipen = 999)

io12 <- read_excel("io12_complete_n46144056.xlsx") 
co2_12 <- read_excel("energy_co2_12_n46144056.xlsx") 


# -------- 2. 取得 12x12 中間投入矩陣 Z 與部門總產出 x、最終需求 f --------
# io12 的前欄為 "投入" 且前 12 列、前 12 欄就是中間投入

io12_df <- as.data.frame(io12)
names(io12_df)[1] <- "sector_name"
sectors <- io12_df$sector_name[1:12]

# Extract Z: assume columns 2..13 correspond to the 12 sectors' flows
Z <- as.matrix(io12_df[1:12, 2:13])
rownames(Z) <- sectors
colnames(Z) <- sectors

# Try to find "產值" row (could be at position >12). We'll search for it.
pv_row_idx <- which(grepl("產值", io12_df$sector_name))
if(length(pv_row_idx) == 0){
  stop("Couldn't find '產值' row in io12_complete.xlsx. Please check row names and adjust.")
}
# take the first occurrence
pv_row <- io12_df[pv_row_idx[1], 2:(1+12)]
x <- as.numeric(pv_row)
names(x) <- sectors

# Try to find final demand vector f: if io12 has columns after 12 that are final demand aggregates
# get column "C+I+G+X-M" 
fd_col <- which(names(io12_df) == "C+I+G+X-M")
if(length(fd_col) == 1){
  f <- as.numeric(io12_df[1:12, fd_col])
} else {
  # fallback: sum final demand columns if present (民間消費, 政府消費, 固定資本形成, 存貨變動, 出口, 進口)
  fd_names <- c("民間消費","政府消費","固定資本形成","存貨變動","出口","進口")
  present <- fd_names[fd_names %in% names(io12_df)]
  if(length(present) > 0){
    f <- rowSums(as.matrix(io12_df[1:12, present]), na.rm = TRUE)
  } else {
    # fallback to zeros
    f <- rep(0, 12)
  }
}
names(f) <- sectors

# Expect co2_12 has a column named "總排放" or "Total_CO2" or similar; try to find numeric column likely total
co2_df <- as.data.frame(co2_12)

possible_tot_cols <- c("總排放","總排放量","Total","Total_CO2","CO2_total","CO2")
tot_col_idx <- which(names(co2_df) %in% possible_tot_cols)
if(length(tot_col_idx) == 0){
  # fallback: take last numeric column
  numeric_cols <- which(sapply(co2_df, is.numeric))
  if(length(numeric_cols) == 0) stop("No numeric column found in CO2 file.")
  tot_col_idx <- tail(numeric_cols, 1)
}
E <- as.numeric(co2_df[[tot_col_idx]])
names(E) <- co2_df[[1]]  # assume first col is sector name
# Align E to sectors order
E12 <- numeric(12)
for(i in 1:12){
  nm <- sectors[i]
  idx <- which(co2_df[[1]] == nm)
  if(length(idx)==1) E12[i] <- as.numeric(co2_df[idx, tot_col_idx]) else E12[i] <- NA
}
if(any(is.na(E12))) warning("Some sectors missing CO2 values; check names.")
E12[is.na(E12)] <- 0
E <- E12
names(E) <- sectors

# -------- technic coefficient A and Leontief L -------
Xmat <- diag(as.numeric(x))
A <- Z %*% solve(Xmat) # or divide columns by x (column j)
I12 <- diag(12)
L <- solve(I12 - A)

# -------- Function to compute impacts under scenarios ----------
run_scenario <- function(carbon_price, pass_through = 1.0, price_elasticity = NULL){
  # carbon_price: NTD per 10^6 ton CO2
  # pass_through: fraction [0,1] of cost that is passed to product price
  # price_elasticity: if NULL, no demand response; if numeric, used to estimate final demand change (approx)
  
  c_direct <- carbon_price * E            # NTD per sector (annual)
  c_passed <- pass_through * c_direct     # amount passed to price
  # percent price increase approx = (c_passed / x)  (assumes c_passed distributed over gross output)
  price_increase_rate <- ifelse(x>0, c_passed / x, 0) 
  
  # If price elasticity provided, approximate final demand change: Delta f = f * (1 + g)^{epsilon} - f ~ f * epsilon * g
  if(!is.null(price_elasticity)){
    eps <- price_elasticity
    delta_f <- f * ( (1 + price_increase_rate)^eps - 1 )
  } else {
    delta_f <- rep(0, 12)
  }
  
  # If interpret c_passed as required additional revenue -> add to final demand shock that compensates producers:
  # convert c_passed into equivalent final demand increase: delta_f_compensate = solve(L) %*% c_passed
  delta_f_compensate <- solve(L) %*% c_passed
  
  # Net final demand change: if consumers reduce demand, apply that; otherwise set to compensate (or a combination)
  # For reporting, give both: (A) full compensation scenario (no demand reduction): use delta_f_compensate
  # (B) with demand reduction: delta_f_compensate + delta_f
  delta_f_no_demand <- as.numeric(delta_f_compensate)
  delta_f_with_demand <- as.numeric(delta_f_compensate + delta_f)
  
  # Outputs
  delta_x_no_demand <- L %*% delta_f_no_demand
  delta_x_with_demand <- L %*% delta_f_with_demand
  
  # value-added effect: if have value-added shares (VA/X), compute reduction in VA from non-pass-through portion
  # Suppose (1-pass_through) * c_direct is absorbed by firms (reduce VA). Distribute it by VA structure if available.
  # For now just report total VA reduction:
  va_absorbed <- sum((1-pass_through) * c_direct)
  
  list(
    carbon_price = carbon_price,
    c_direct = c_direct,
    c_passed = c_passed,
    price_increase_rate = price_increase_rate,
    delta_f_no_demand = delta_f_no_demand,
    delta_x_no_demand = as.numeric(delta_x_no_demand),
    delta_f_with_demand = delta_f_with_demand,
    delta_x_with_demand = as.numeric(delta_x_with_demand),
    va_absorbed = va_absorbed
  )
}
# --- Extract Value-Added components (Labour, Operating Surplus, Capital Consumption) ---

labour_idx <- which(grepl("勞動報酬", io12_df$sector_name))
surplus_idx <- which(grepl("營業盈餘", io12_df$sector_name))
capital_idx <- which(grepl("固定資本消耗", io12_df$sector_name))

labour_va <- if(length(labour_idx) > 0) as.numeric(io12_df[labour_idx[1], 2:13]) else rep(NA, 12)
surplus_va <- if(length(surplus_idx) > 0) as.numeric(io12_df[surplus_idx[1], 2:13]) else rep(NA, 12)
capital_va <- if(length(capital_idx) > 0) as.numeric(io12_df[capital_idx[1], 2:13]) else rep(NA, 12)

names(labour_va) <- names(surplus_va) <- names(capital_va) <- sectors

VA_sum <- labour_va + surplus_va + capital_va
VA_sum[VA_sum == 0] <- NA

# 各部門 VA 占比
labour_share <- labour_va / VA_sum
surplus_share <- surplus_va / VA_sum
capital_share <- capital_va / VA_sum

# -------- Run scenarios (examples) ----------
scenarios <- list(
  list(p=300, tau=1, eps=NULL),      # 現行碳價，全額轉嫁，不考慮需求彈性
  list(p=300, tau=1, eps=-0.5),      # 現行碳價，全額轉嫁，考慮需求彈性
  list(p=300, tau=0.5, eps=NULL),    # 現行碳價，部分轉嫁，不考慮需求彈性
  list(p=300, tau=0.5, eps=-0.5),    # 現行碳價，部分轉嫁，考慮需求彈性
  list(p=3000, tau=1, eps=NULL),     # 高碳價，全額轉嫁，不考慮需求彈性
  list(p=3000, tau=1, eps=-0.5),     # 高碳價，全額轉嫁，考慮需求彈性
  list(p=3000, tau=0.5, eps=NULL),   # 高碳價，部分轉嫁，不考慮需求彈性
  list(p=3000, tau=0.5, eps=-0.5)    # 高碳價，部分轉嫁，考慮需求彈性
)

run_scenario <- function(carbon_price, pass_through = 1.0, price_elasticity = NULL){
  
  c_direct <- carbon_price * E
  c_passed <- pass_through * c_direct
  c_absorb <- (1 - pass_through) * c_direct
  
  price_increase_rate <- ifelse(x > 0, c_passed / x, 0)
  
  # Demand elasticity effect
  if(!is.null(price_elasticity)){
    delta_f_demand <- f * ((1 + price_increase_rate)^price_elasticity - 1)
  } else {
    delta_f_demand <- rep(0, 12)
  }
  
  # Compensation demand
  delta_f_comp <- solve(L) %*% c_passed
  
  # Final demand shocks
  delta_f_no <- delta_f_comp
  delta_f_yes <- delta_f_comp + delta_f_demand
  
  # Output changes
  delta_x_no <- L %*% delta_f_no
  delta_x_yes <- L %*% delta_f_yes
  
  # === NEW: Value-added 分配 ===
  # 若產業需吸收成本，依 VA 結構拆分
  va_labor_loss   <- labour_share  * c_absorb
  va_surplus_loss <- surplus_share * c_absorb
  va_capital_loss <- capital_share * c_absorb
  
  return(list(
    carbon_price = carbon_price,
    c_direct = c_direct,
    c_passed = c_passed,
    c_absorb = c_absorb,
    price_increase_rate = price_increase_rate,
    delta_x_no_demand = as.numeric(delta_x_no),
    delta_x_with_demand = as.numeric(delta_x_yes),
    va_loss_labor   = va_labor_loss,
    va_loss_surplus = va_surplus_loss,
    va_loss_capital = va_capital_loss
  ))
}
results <- list()
for(i in seq_along(scenarios)){
  s <- scenarios[[i]]
  results[[i]] <- run_scenario(carbon_price = s$p, pass_through = s$tau, price_elasticity = s$eps)
}
# -------- Summary --------

scenario_prices <- sapply(results, function(r) r$carbon_price)

# === Wide Summary ===
summ_wide <- data.frame(
  sector = sectors,
  emissions_ton = E,
  gross_output = x,
  stringsAsFactors = FALSE
)

for (i in seq_along(results)) {
  r <- results[[i]]
  p <- r$carbon_price
  
  summ_wide[[paste0("c_direct_p", p)]]    <- r$c_direct
  summ_wide[[paste0("c_passed_p", p)]]    <- r$c_passed
  summ_wide[[paste0("c_absorb_p", p)]]    <- r$c_absorb
  
  summ_wide[[paste0("dx_no_p", p)]]       <- r$delta_x_no_demand
  summ_wide[[paste0("dx_yes_p", p)]]      <- r$delta_x_with_demand
  
  summ_wide[[paste0("va_labor_p", p)]]    <- r$va_loss_labor
  summ_wide[[paste0("va_surplus_p", p)]]  <- r$va_loss_surplus
  summ_wide[[paste0("va_capital_p", p)]]  <- r$va_loss_capital
}


summ_long_list <- list()

for (i in seq_along(results)) {
  r <- results[[i]]
  p <- r$carbon_price
  
  df_tmp <- data.frame(
    sector = sectors,
    price = p,
    c_direct = r$c_direct,
    c_passed = r$c_passed,
    c_absorb = r$c_absorb,
    dx_no = r$delta_x_no_demand,
    dx_yes = r$delta_x_with_demand,
    va_labor = r$va_loss_labor,
    va_surplus = r$va_loss_surplus,
    va_capital = r$va_loss_capital
  )
  
  summ_long_list[[i]] <- df_tmp
}

summ_long <- dplyr::bind_rows(summ_long_list)

head(summ_wide)
head(summ_long)

##############
library(ggplot2)
library(dplyr)
library(scales) 
library(forcats)

ggplot(summ_long, aes(x = fct_reorder(sector, c_direct, .fun = sum, .desc =FALSE),
                      y = c_direct, fill = factor(price))) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_y_continuous(labels = comma) + 
  labs(title = "各部門直接碳成本（不同碳價情境）",
       x = "", y = "直接碳成本花費(單位：百萬元)",
       fill = "碳價") +
  theme_bw()

##########

summ_long_ext <- summ_long %>%
  mutate(
    tau = case_when(
      c_absorb == 0 ~ "全額轉嫁給產品價格",
      TRUE ~ "部分成本企業自己吸收"
    ),
    elasticity = case_when(
      dx_no == dx_yes ~ "elasticity=0",
      TRUE ~ "elasticity!=0"
    )
  )

ggplot(summ_long_ext, aes(x=factor(price), y=dx_no, fill=sector)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid( ~tau ) +
  scale_y_continuous(labels = comma) + 
  labs(title="產出變動 Δx (單位：百萬元)", x="碳價", y="Δx", fill="部門") +
  theme_bw()

###########
all_sectors_elastic <- summ_long_ext %>%
  filter(elasticity != "elasticity=0") %>%
  select(sector, price, tau, elasticity, dx_no, dx_yes) %>%
  pivot_longer(cols = c(dx_no, dx_yes),
               names_to = "type",
               values_to = "dx") %>%
  mutate(type = recode(type,
                       dx_no = "需求不變",
                       dx_yes = "考慮需求彈性"))

# 2️⃣ 作圖
ggplot(all_sectors_elastic, aes(x = factor(price), y = dx, fill = type)) +
  geom_col(position = "dodge") +
  facet_grid(tau ~ sector) +          # 行：轉嫁比例 τ，列：部門
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("需求不變" = "#1f77b4", "考慮需求彈性" = "#ff7f0e")) +
  labs(title = "各部門 Δx（考慮需求彈性情境）",
       x = "碳價",
       y = "產出變動 Δx (百萬元)",
       fill = "需求狀態") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))
###########

va_long <- summ_long_ext %>%
  select(sector, price, tau, elasticity, va_labor, va_surplus, va_capital) %>%
  pivot_longer(cols = c(va_labor, va_surplus, va_capital),
               names_to = "VA_type",
               values_to = "VA_loss") %>%
  mutate(VA_type = recode(VA_type,
                          va_labor = "勞動報酬",
                          va_surplus = "營業盈餘",
                          va_capital = "固定資本消耗"))

va_long_price <- va_long %>%
  group_by(sector, price, VA_type) %>%
  summarise(VA_loss = mean(VA_loss), .groups="drop")

# 繪圖
ggplot(va_long_price, aes(x=factor(price), y=abs(VA_loss), fill=VA_type)) +
  geom_col(position="stack") +
  scale_y_continuous(labels = comma) +
  facet_wrap(~sector, ncol=4, scales="free_y") +
  labs(title="各部門原始投入吸收的碳成本（僅比較碳價差異）",
       x="碳價", y="吸收碳成本 (單位：百萬元)", fill="原始投入分項") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###

delta_f_transport <- rep(0, 12)   # 12個部門
names(delta_f_transport) <- sectors
delta_f_transport["運輸部門"] <- 1  # 假設衝擊值=1

# Leontief 矩陣看 Δx
delta_x_transport <- L %*% delta_f_transport
delta_x_transport <- as.numeric(delta_x_transport)
names(delta_x_transport) <- sectors

df_transport <- data.frame(
  sector = sectors,
  delta_x = delta_x_transport
)

library(ggplot2)
ggplot(df_transport, aes(x=fct_reorder(sector, delta_x), y=delta_x, fill=sector)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = comma) + 
  labs(title="運輸部門衝擊帶動其他部門產出變動(單位：百萬元)",
       x="", y="Δx") +
  theme_bw() +
  theme(legend.position = "none")
scenarios <- summ_long_ext %>%
  distinct(price, tau, elasticity)

#####################
library(lpSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

n <- 12
sectors <- names(x)

# 碳價與轉嫁比例
carbon_price <- 300
tau <- 0.5
c_direct <- carbon_price * E
c_passed <- tau * c_direct
c_absorb <- (1 - tau) * c_direct
va_labor_loss <- labour_share * c_absorb
va_surplus_loss <- surplus_share * c_absorb
va_capital_loss <- capital_share * c_absorb

# 產出上下限
x_min <- f * 0.5   # 至少滿足一半的最終需求
x_max <- x * 2               # 假設上限為 2 倍原始產出

# =====  約束矩陣 =====
# x >= x_min, x <= x_max
const.mat <- rbind(diag(n), -diag(n))
const.dir <- c(rep(">=", n), rep(">=", n))
const.rhs <- c(as.numeric(x_min), -as.numeric(x_max))

# =====目標函數 =====
# 1：最小化企業吸收碳成本
# 2：最小總成本 (企業+價格傳導)
# 3：最大總產出
# 4：最小勞動損失

f_obj_list <- list(
  cost_min   = c_absorb ,
  total_cost = (c_absorb + c_passed) ,
  output_max = -x ,
  labor_min  = va_labor_loss
)

# =====  LP =====
results_lp <- list()
for(name in names(f_obj_list)){
  lp_sol <- lp(direction="min",
               objective.in = f_obj_list[[name]],
               const.mat = const.mat,
               const.dir = const.dir,
               const.rhs = const.rhs)
  results_lp[[name]] <- lp_sol$solution
  names(results_lp[[name]]) <- sectors
}


results_lp <- list()
for(name in names(f_obj_list)){
  lp_sol <- lp(direction="min",
               objective.in = f_obj_list[[name]],
               const.mat = const.mat,
               const.dir = const.dir,
               const.rhs = const.rhs)
  results_lp[[name]] <- lp_sol$solution
  names(results_lp[[name]]) <- sectors
}

results_df <- data.frame(
  sector = sectors,
  x_base = x,
  x_cost_min   = results_lp$cost_min,
  x_total_cost = results_lp$total_cost,
  x_output_max = results_lp$output_max,
  x_labor_min  = results_lp$labor_min
)

results_va <- lapply(names(results_lp), function(name){
  x_scenario <- results_lp[[name]]
  
  # 假設用 f_obj_list[[name]] 對應的成本作為吸收量
  c_absorb_s <- f_obj_list[[name]]
  
  data.frame(
    sector = sectors,
    scenario = name,
    va_labor   = labour_share   * c_absorb_s,
    va_surplus = surplus_share  * c_absorb_s,
    va_capital = capital_share  * c_absorb_s
  )
})

va_df <- do.call(rbind, results_va)

va_long <- va_df %>%
  pivot_longer(cols = c(va_labor, va_surplus, va_capital),
               names_to = "VA_type", values_to = "VA_loss") %>%
  mutate(VA_type = recode(VA_type,
                          va_labor="勞動報酬",
                          va_surplus="營業盈餘",
                          va_capital="固定資本消耗"))


va_long <- va_long %>%
  mutate(scenario_label = recode(scenario,
                                 cost_min   = "企業吸收碳成本最小",
                                 total_cost = "總成本最小",
                                 output_max = "總產出最大",
                                 labor_min  = "勞動損失最小"))

ggplot(va_long, aes(x=fct_reorder(sector, VA_loss), y=abs(VA_loss), fill=VA_type)) +
  geom_col(position="stack") +
  scale_y_continuous(labels = scales::comma) +
  coord_flip() +
  facet_wrap(~scenario_label, ncol=2, scales="free_x") +
  labs(title="各部門原始投入承擔之碳成本（線性規劃不同目標）",
       x="", y="承擔碳成本(百萬元)", fill="原始投入分類") +
  theme_bw()
