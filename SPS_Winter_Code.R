rm(list = ls())

# -- Install packages (run once, then comment out) ----------
# install.packages("tidyverse")
# install.packages("fixest")
# install.packages("haven")
# install.packages("car")
# install.packages("kableExtra")
# install.packages("modelsummary")
# install.packages("gt")
# install.packages("broom")
# install.packages("qrcode")
# install.packages("base64enc")
# install.packages("multcomp")
# install.packages("webshot2")
# webshot2::install_chromium()   # only needed once for webshot2
# -----------------------------------------------------------

library(tidyverse)
library(fixest)
library(haven)
library(car)
library(kableExtra)
library(modelsummary)
library(gt)
library(broom)
library(qrcode)
library(base64enc)
library(multcomp)
library(webshot2)
select <- dplyr::select

# -- USER SETTINGS: adjust these two paths ------------------
# data_dir : folder containing all .rds SOEP files
# plot_dir : folder where all output plots/tables are saved
data_dir <- "~/Desktop/UNI/Schwerpunkt Seminar/Winter/Datensatz/soepdata"
plot_dir  <- "~/Desktop/UNI/Schwerpunkt Seminar/Winter/Plots"
# -----------------------------------------------------------

setwd(data_dir)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# 1. Data Loading --------------------------------------------

pl <- readRDS("pl_selected_tom.rds") |>
  zap_labels() |>
  select(pid, syear, plb0022_h)

pgen <- readRDS("pgen.rds") |>
  zap_labels() |>
  select(pid, syear, pgexpft, pgexppt, pgbilzeit, pgfamstd, pglabgro, pgemplst)

ppfad <- readRDS("ppathl.rds") |>
  zap_labels() |>
  select(pid, hid, syear, gebjahr, sex) |>
  distinct(pid, syear, .keep_all = TRUE)

region <- readRDS("regionl.rds") |>
  zap_labels() |>
  select(hid, syear, bula) |>
  mutate(east = as.integer(bula %in% c(12, 13, 14, 15, 16)))

equiv <- readRDS("pequiv.rds") |>
  zap_labels() |>
  select(pid, syear, i11102)

health <- readRDS("health.rds") |>
  zap_labels() |>
  select(pid, syear, pcs, mcs)

df_raw <- pl |>
  left_join(pgen,   by = c("pid", "syear")) |>
  left_join(ppfad,  by = c("pid", "syear")) |>
  left_join(equiv,  by = c("pid", "syear")) |>
  left_join(region, by = c("hid", "syear")) |>
  left_join(health, by = c("pid", "syear"))


# 2. Sample Construction – Robustness Control Groups ---------

create_sample <- function(data, lower_bound_control) {
  data |>
    mutate(
      age = syear - gebjahr,
      lfp = case_when(
        plb0022_h %in% c(1, 2)  ~ 1L,
        plb0022_h %in% c(3:11)  ~ 0L,
        TRUE                     ~ NA_integer_
      ),
      contrib_years = pgexpft + 0.5 * pgexppt,
      treat = case_when(
        contrib_years >= 45                                        ~ 1L,
        contrib_years >= lower_bound_control & contrib_years < 45 ~ 0L,
        TRUE                                                       ~ NA_integer_
      ),
      post = case_when(
        syear %in% 2015:2019 ~ 1L,
        syear %in% 2009:2013 ~ 0L,
        TRUE                 ~ NA_integer_
      ),
      did = treat * post
    ) |>
    filter(syear != 2014) |>
    filter(age >= 62, age <= 65) |>
    filter(!is.na(treat), !is.na(post), !is.na(lfp))
}

df_35_44 <- create_sample(df_raw, 35)
df_40_44 <- create_sample(df_raw, 40)
df_41_44 <- create_sample(df_raw, 41)
df_42_44 <- create_sample(df_raw, 42)

m_35_44 <- feols(lfp ~ treat + post + did, data = df_35_44, cluster = ~pid)
m_40_44 <- feols(lfp ~ treat + post + did, data = df_40_44, cluster = ~pid)
m_41_44 <- feols(lfp ~ treat + post + did, data = df_41_44, cluster = ~pid)
m_42_44 <- feols(lfp ~ treat + post + did, data = df_42_44, cluster = ~pid)

summary(m_35_44)
summary(m_40_44)
summary(m_41_44)
summary(m_42_44)


# 3. Parallel Trends Test ------------------------------------

create_sample_es <- function(data, lower_bound_control) {
  data |>
    mutate(
      age = syear - gebjahr,
      lfp = case_when(
        plb0022_h %in% c(1, 2) ~ 1L,
        plb0022_h %in% c(3:11) ~ 0L,
        TRUE                    ~ NA_integer_
      ),
      contrib_years = pgexpft + 0.5 * pgexppt,
      treat = case_when(
        contrib_years >= 45                                        ~ 1L,
        contrib_years >= lower_bound_control & contrib_years < 45 ~ 0L,
        TRUE                                                       ~ NA_integer_
      )
    ) |>
    filter(age >= 62, age <= 65) |>
    filter(!is.na(treat), !is.na(lfp)) |>
    filter(syear %in% 2009:2013)
}

es_35_44 <- feols(lfp ~ i(syear, treat, ref = 2013) | pid + syear,
                  data = create_sample_es(df_raw, 35), cluster = ~pid)
es_40_44 <- feols(lfp ~ i(syear, treat, ref = 2013) | pid + syear,
                  data = create_sample_es(df_raw, 40), cluster = ~pid)
es_41_44 <- feols(lfp ~ i(syear, treat, ref = 2013) | pid + syear,
                  data = create_sample_es(df_raw, 41), cluster = ~pid)
es_42_44 <- feols(lfp ~ i(syear, treat, ref = 2013) | pid + syear,
                  data = create_sample_es(df_raw, 42), cluster = ~pid)

pre_years <- c("syear::2009:treat = 0",
               "syear::2010:treat = 0",
               "syear::2011:treat = 0",
               "syear::2012:treat = 0")

pt_results <- list(
  "35-44" = linearHypothesis(es_35_44, pre_years),
  "40-44" = linearHypothesis(es_40_44, pre_years),
  "41-44" = linearHypothesis(es_41_44, pre_years),
  "42-44" = linearHypothesis(es_42_44, pre_years)
)

pt_table <- data.frame(
  Modell   = names(pt_results),
  Chisq    = sapply(pt_results, function(x) round(x[2, "Chisq"], 3)),
  df       = sapply(pt_results, function(x) x[2, "Df"]),
  p_Wert   = sapply(pt_results, function(x) round(x[2, "Pr(>Chisq)"], 4)),
  Ergebnis = sapply(pt_results, function(x)
    ifelse(x[2, "Pr(>Chisq)"] > 0.05, "✓ Parallel", "✗ Verletzt"))
)

print(pt_table)

pval_35 <- round(pt_results[["35-44"]][2, "Pr(>Chisq)"], 3)
pval_40 <- round(pt_results[["40-44"]][2, "Pr(>Chisq)"], 3)
pval_41 <- round(pt_results[["41-44"]][2, "Pr(>Chisq)"], 3)
pval_42 <- round(pt_results[["42-44"]][2, "Pr(>Chisq)"], 3)


# 3a. Parallel Trends Table – HTML --------------------------

pt_table |>
  kbl(
    format    = "html",
    row.names = FALSE,
    col.names = c("Model", "Chi²", "df", "p-Value", "Result"),
    caption   = "Parallel Trends Test – Pre-Treatment Periode (2009–2012) | Age 62–65",
    digits    = 4
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width        = FALSE,
    position          = "center"
  ) |>
  column_spec(5,
              color = ifelse(pt_table$Ergebnis == "✓ Parallel", "green", "red"),
              bold  = TRUE
  ) |>
  add_footnote(
    "H0: Alle Pre-Treatment Koeffizienten = 0 | Clustered SE auf pid-Ebene | Sample: Age 62–65",
    notation = "none"
  ) |>
  save_kable(file.path(plot_dir, "parallel_trends_table.html"))
cat("Gespeichert: parallel_trends_table.html\n")


# 3b. Parallel Trends Plots ----------------------------------

plot_trends <- function(data, lower_label, filename, pt_pval) {
  means <- data |>
    group_by(syear, treat) |>
    summarise(lfp_mean = mean(lfp, na.rm = TRUE), .groups = "drop") |>
    mutate(gruppe = ifelse(treat == 1,
                           "Treatment (≥ 45 years of contribution)",
                           paste0("Control (", lower_label, " years of contribution)")))

  ggplot(means, aes(x = syear, y = lfp_mean, color = gruppe, group = gruppe)) +
    geom_vline(xintercept = 2014, linetype = "dashed",
               color = "#CC7C2E", linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    annotate("text", x = 2014.15, y = max(means$lfp_mean),
             label = "reform (July 2014)", hjust = 0, color = "#CC7C2E", size = 3.5) +
    annotate("rect", xmin = 2009, xmax = 2013, ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = "#1E2B5C") +
    annotate("rect", xmin = 2015, xmax = 2019, ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = "green") +
    annotate("text", x = 2011, y = min(means$lfp_mean) - 0.01,
             label = "pre-period", color = "#1E2B5C", size = 3) +
    annotate("text", x = 2017, y = min(means$lfp_mean) - 0.01,
             label = "post-period", color = "darkgreen", size = 3) +
    scale_color_manual(
      name   = "",
      values = setNames(
        c("#35746E", "#1E2B5C"),
        c("Treatment (≥ 45 years of contribution)",
          paste0("Control (", lower_label, " years of contribution)"))
      )
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(2009:2013, 2015:2019)) +
    labs(
      title    = "Employed: treatment- vs. controlgroup (2009–2019)",
      subtitle = paste0("Control: ", lower_label, " contribution years | 2014 excluded | Parallel Trends: p = ", pt_pval),
      x        = "Year", y = "Share employed",
      caption  = "SOEP-Data | Age 62–65 | Böhm & Wilhelmi 2026"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  ggsave(file.path(plot_dir, filename), width = 10, height = 5.5, dpi = 300)
  cat("Gespeichert:", filename, "\n")
}

plot_trends(df_35_44, "35–44", "parallel_trends_35_44.png", pval_35)
plot_trends(df_40_44, "40–44", "parallel_trends_40_44.png", pval_40)
plot_trends(df_41_44, "41–44", "parallel_trends_41_44.png", pval_41)
plot_trends(df_42_44, "42–44", "parallel_trends_42_44.png", pval_42)


# 3c. Parallel Trends – Combined HTML with Plots -------------

png_to_base64 <- function(filepath) {
  paste0("data:image/png;base64,", base64encode(filepath))
}

img_35 <- png_to_base64(file.path(plot_dir, "parallel_trends_35_44.png"))
img_40 <- png_to_base64(file.path(plot_dir, "parallel_trends_40_44.png"))
img_41 <- png_to_base64(file.path(plot_dir, "parallel_trends_41_44.png"))
img_42 <- png_to_base64(file.path(plot_dir, "parallel_trends_42_44.png"))

html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
  body { font-family: sans-serif; padding: 2rem; max-width: 1200px; margin: auto; }
  h1 { color: #1E2B5C; }
  h2 { color: #35746E; margin-top: 2rem; }
  table { width: 100%; border-collapse: collapse; margin-bottom: 2rem; }
  th { background: #1E2B5C; color: white; padding: 10px; text-align: left; }
  td { padding: 8px 10px; border-bottom: 1px solid #eee; }
  tr:nth-child(4) { background: #f0f7f0; font-weight: bold; }
  .sig-yes { color: green; font-weight: bold; }
  .sig-no  { color: red;   font-weight: bold; }
  .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
  img { width: 100%; border: 1px solid #eee; border-radius: 4px; }
  .footnote { font-size: 11px; color: #888; margin-top: 0.5rem; }
</style>
</head>
<body>

<h1>Parallel Trends – Full Overview</h1>
<p style="color:#888; font-size:13px;">
  Sample: Age 62–65 | Treatment: ≥45 contribution years |
  Pre-period: 2009–2013 | Böhm & Wilhelmi 2026
</p>

<h2>Summary Table</h2>
<table>
<thead>
  <tr>
    <th>Model</th><th>Control Group</th>
    <th>N Treat</th><th>N Control</th><th>Ratio (C:T)</th>
    <th>DiD Coeff.</th><th>p (DiD)</th>
    <th>Chi²</th><th>df</th><th>p (PT)</th><th>PT Result</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>35–44</td><td>35–44 contrib. years</td>
    <td>572</td><td>5,981</td><td style="color:red;font-weight:bold">10.5:1</td>
    <td>', round(coef(m_35_44)["did"], 3), '</td>
    <td>', round(tidy(m_35_44)$p.value[tidy(m_35_44)$term=="did"], 3), '</td>
    <td>', round(pt_results[["35-44"]][2,"Chisq"], 3), '</td>
    <td>4</td>
    <td>', pval_35, '</td>
    <td class="sig-yes">✓ Parallel</td>
  </tr>
  <tr>
    <td>40–44</td><td>40–44 contrib. years</td>
    <td>572</td><td>2,866</td><td>5.0:1</td>
    <td>', round(coef(m_40_44)["did"], 3), '</td>
    <td>', round(tidy(m_40_44)$p.value[tidy(m_40_44)$term=="did"], 3), '</td>
    <td>', round(pt_results[["40-44"]][2,"Chisq"], 3), '</td>
    <td>4</td>
    <td>', pval_40, '</td>
    <td class="sig-yes">✓ Parallel</td>
  </tr>
  <tr>
    <td>41–44</td><td>41–44 contrib. years</td>
    <td>572</td><td>2,159</td><td>3.8:1</td>
    <td>', round(coef(m_41_44)["did"], 3), '</td>
    <td>', round(tidy(m_41_44)$p.value[tidy(m_41_44)$term=="did"], 3), '</td>
    <td>', round(pt_results[["41-44"]][2,"Chisq"], 3), '</td>
    <td>4</td>
    <td>', pval_41, '</td>
    <td class="sig-yes">✓ Parallel</td>
  </tr>
  <tr style="background:#f0f7f0; font-weight:bold">
    <td>42–44 ★</td><td>42–44 contrib. years (Main)</td>
    <td>572</td><td>1,496</td><td>2.6:1</td>
    <td>', round(coef(m_42_44)["did"], 3), '</td>
    <td>', round(tidy(m_42_44)$p.value[tidy(m_42_44)$term=="did"], 3), '</td>
    <td>', round(pt_results[["42-44"]][2,"Chisq"], 3), '</td>
    <td>4</td>
    <td>', pval_42, '</td>
    <td class="sig-yes">✓ Parallel</td>
  </tr>
</tbody>
</table>
<p class="footnote">
  H0: Alle Pre-Treatment Koeffizienten = 0 | Clustered SE auf pid-Ebene |
  ★ = Main specification | Chi² joint test auf alle Pre-Treatment Koeffizienten
</p>

<h2>Parallel Trends Plots</h2>
<div class="grid">
  <div><img src="', img_35, '" alt="PT 35-44"></div>
  <div><img src="', img_40, '" alt="PT 40-44"></div>
  <div><img src="', img_41, '" alt="PT 41-44"></div>
  <div><img src="', img_42, '" alt="PT 42-44"></div>
</div>

</body>
</html>
')

writeLines(html_content, file.path(plot_dir, "parallel_trends_full.html"))
cat("Gespeichert: parallel_trends_full.html\n")


# 4. Main Sample (df_main: 42–44, Age 62–65) ----------------
# Control group 42–44: best common support + PT p = 0.097

df_main <- df_raw |>
  mutate(
    age           = syear - gebjahr,
    female        = as.integer(sex == 2),
    married       = as.integer(pgfamstd == 1),
    bildung       = pgbilzeit,
    income        = i11102,
    health_phys   = ifelse(pcs < 0, NA, pcs),
    health_ment   = ifelse(mcs < 0, NA, mcs),
    lfp = case_when(
      plb0022_h %in% c(1, 2)  ~ 1L,
      plb0022_h %in% c(3:11)  ~ 0L,
      TRUE                     ~ NA_integer_
    ),
    fulltime = case_when(
      pgemplst == 1             ~ 1L,
      pgemplst %in% c(2, 4, 7) ~ 0L,
      TRUE                      ~ NA_integer_
    ),
    parttime = case_when(
      pgemplst %in% c(2, 4, 7) ~ 1L,
      pgemplst == 1             ~ 0L,
      TRUE                      ~ NA_integer_
    ),
    contrib_years = pgexpft + 0.5 * pgexppt,
    treat = case_when(
      contrib_years >= 45                      ~ 1L,
      contrib_years >= 42 & contrib_years < 45 ~ 0L,
      TRUE                                     ~ NA_integer_
    ),
    post = case_when(
      syear %in% 2015:2019 ~ 1L,
      syear %in% 2009:2013 ~ 0L,
      TRUE                 ~ NA_integer_
    ),
    did         = treat * post,
    east        = as.integer(bula %in% c(12, 13, 14, 15, 16)),
    high_income = as.integer(income  > median(income,  na.rm = TRUE)),
    high_educ   = as.integer(bildung > median(bildung, na.rm = TRUE))
  ) |>
  filter(syear != 2014) |>
  filter(age >= 62, age <= 65) |>
  filter(!is.na(treat), !is.na(post), !is.na(lfp))


# 5. Main Models ---------------------------------------------

m1 <- feols(lfp ~ treat + post,        data = df_main, cluster = ~pid)
m2 <- feols(lfp ~ treat + post + did,  data = df_main, cluster = ~pid)
m3 <- feols(lfp ~ treat + did | syear, data = df_main, cluster = ~pid)
m4 <- feols(lfp ~ post  + did | pid,   data = df_main, cluster = ~pid)  # main model

summary(m1)
summary(m2)
summary(m3)
summary(m4)


# 5a. Main Models – HTML Table -------------------------------

modelsummary(
  list("Model (1)" = m1, "Model (2)" = m2,
       "Model (3)" = m3, "Model (4)" = m4),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c("treat" = "Treat", "post" = "Post", "did" = "Treat × Post"),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~`Model (1)`, ~`Model (2)`, ~`Model (3)`, ~`Model (4)`,
    "Year FE",       "Nein",       "Nein",        "Ja",         "Nein",
    "Individual FE", "Nein",       "Nein",        "Nein",       "Ja"
  )
) |>
  tab_header(
    title    = "Main Results – Effect of 'Rente mit 63'",
    subtitle = "Dep. Var.: LFP | Control: 42–44 contribution years | Age 62–65"
  ) |>
  tab_spanner(label = "Without FE",         id = "sp_ofe",
              columns = c("Model (1)", "Model (2)")) |>
  tab_spanner(label = "With Fixed Effects", id = "sp_mfe",
              columns = c("Model (3)", "Model (4)")) |>
  tab_source_note("Clustered SE auf pid-Ebene | Sample: Age 62–65 | Main control group: 42–44 contribution years") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(75)) |>
  gtsave(file.path(plot_dir, "final_hauptmodelle.html"))
cat("Gespeichert: final_hauptmodelle.html\n")

modelsummary(
  list("Model (1)" = m1, "Model (2)" = m2,
       "Model (3)" = m3, "Model (4)" = m4),
  output = "table_main.html",
  stars = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj|RMSE|R2",
  coef_rename = c("treat" = "Treat", "post" = "Post", "did" = "Treat × Post"),
  gof_map = list(
    list(raw = "nobs", clean = "Observations", fmt = 0)
  ),
  add_rows = tribble(
    ~term,           ~`Model (1)`, ~`Model (2)`, ~`Model (3)`, ~`Model (4)`,
    "Year FE",       "No",         "No",          "Yes",        "No",
    "Individual FE", "No",         "No",          "No",         "Yes"
  ),
  notes = "Clustered SE at individual level. Control: 42-44 contribution years. Sample: men aged 62-65, 2009-2019 (excl. 2014)."
)

webshot2::webshot(
  "table_main.html",
  file.path(plot_dir, "table_main.png"),
  zoom = 2,
  vwidth = 700
)


# 5b. Robustness – Alternative Control Groups – HTML ---------

modelsummary(
  list("35–44" = m_35_44, "40–44" = m_40_44,
       "41–44" = m_41_44, "42–44" = m_42_44),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c("treat" = "Treat", "post" = "Post", "did" = "Treat × Post"),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  )
) |>
  tab_header(
    title    = "Robustness – Alternative Control Groups",
    subtitle = "Dep. Var.: LFP | Varying lower bound of control group | Age 62–65"
  ) |>
  tab_spanner(label = "Control group (contribution years)", id = "sp_rob",
              columns = c("35–44", "40–44", "41–44", "42–44")) |>
  tab_source_note("Clustered SE auf pid-Ebene | 42–44 is main specification") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(75)) |>
  gtsave(file.path(plot_dir, "final_robustness.html"))
cat("Gespeichert: final_robustness.html\n")


# 6. Heterogeneity Analysis ----------------------------------

# RQ1 – East vs. West
m_ow_interact <- feols(lfp ~ treat + post + did + east + did:east,
                       data = df_main, cluster = ~pid)
m_ow_fe       <- feols(lfp ~ post + did + did:east | pid,
                       data = df_main, cluster = ~pid)

summary(m_ow_interact)
summary(m_ow_fe)

modelsummary(
  list("Interaction" = m_ow_interact, "Interaction FE" = m_ow_fe),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c(
    "treat"    = "Treat", "post" = "Post", "did" = "Treat × Post",
    "east"     = "East",  "did:east" = "Treat × Post × East"
  ),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~Interaction, ~`Interaction FE`,
    "Individual FE", "Nein",       "Ja"
  )
) |>
  tab_header(
    title    = "RQ1 – Regional Heterogeneity: East vs. West",
    subtitle = "Dep. Var.: LFP | Control: 42–44 contribution years | Age 62–65"
  ) |>
  tab_source_note("Clustered SE auf pid-Ebene | east absorbiert durch Individual FE im FE-Modell") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(60)) |>
  gtsave(file.path(plot_dir, "final_rq1_ostwest.html"))
cat("Gespeichert: final_rq1_ostwest.html\n")

# RQ2 – Income
m_inc_interact <- feols(lfp ~ treat + post + did + high_income + did:high_income,
                        data = df_main, cluster = ~pid)
m_inc_fe       <- feols(lfp ~ post + did + did:high_income | pid,
                        data = df_main, cluster = ~pid)

summary(m_inc_interact)
summary(m_inc_fe)

modelsummary(
  list("Interaction" = m_inc_interact, "Interaction FE" = m_inc_fe),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c(
    "treat"           = "Treat", "post" = "Post", "did" = "Treat × Post",
    "high_income"     = "High Income",
    "did:high_income" = "Treat × Post × High Income"
  ),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~Interaction, ~`Interaction FE`,
    "Individual FE", "Nein",       "Ja"
  )
) |>
  tab_header(
    title    = "RQ2 – Socioeconomic Heterogeneity: Income",
    subtitle = "Dep. Var.: LFP | Split at median household income | Age 62–65"
  ) |>
  tab_source_note("Clustered SE auf pid-Ebene") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(60)) |>
  gtsave(file.path(plot_dir, "final_rq2_einkommen.html"))
cat("Gespeichert: final_rq2_einkommen.html\n")

# RQ2 – Education
m_edu_interact <- feols(lfp ~ treat + post + did + high_educ + did:high_educ,
                        data = df_main, cluster = ~pid)
m_edu_fe       <- feols(lfp ~ post + did + did:high_educ | pid,
                        data = df_main, cluster = ~pid)

summary(m_edu_interact)
summary(m_edu_fe)

modelsummary(
  list("Interaction" = m_edu_interact, "Interaction FE" = m_edu_fe),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c(
    "treat"         = "Treat", "post" = "Post", "did" = "Treat × Post",
    "high_educ"     = "High Education",
    "did:high_educ" = "Treat × Post × High Education"
  ),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~Interaction, ~`Interaction FE`,
    "Individual FE", "Nein",       "Ja"
  )
) |>
  tab_header(
    title    = "RQ2 – Socioeconomic Heterogeneity: Education",
    subtitle = "Dep. Var.: LFP | Split at median education years | Age 62–65"
  ) |>
  tab_source_note("Clustered SE auf pid-Ebene") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(60)) |>
  gtsave(file.path(plot_dir, "final_rq2_bildung.html"))
cat("Gespeichert: final_rq2_bildung.html\n")

# Combined Heterogeneity Table
modelsummary(
  list(
    "Ost/West"     = m_ow_interact,
    "Ost/West FE"  = m_ow_fe,
    "Einkommen"    = m_inc_interact,
    "Einkommen FE" = m_inc_fe,
    "Bildung"      = m_edu_interact,
    "Bildung FE"   = m_edu_fe
  ),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c(
    "treat"           = "Treat",
    "post"            = "Post",
    "did"             = "Treat × Post",
    "east"            = "East",
    "did:east"        = "Treat × Post × East",
    "high_income"     = "High Income",
    "did:high_income" = "Treat × Post × High Income",
    "high_educ"       = "High Education",
    "did:high_educ"   = "Treat × Post × High Education"
  ),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~`Ost/West`, ~`Ost/West FE`, ~`Einkommen`, ~`Einkommen FE`, ~`Bildung`, ~`Bildung FE`,
    "Individual FE", "Nein",      "Ja",            "Nein",       "Ja",            "Nein",     "Ja"
  )
) |>
  tab_header(
    title    = "Heterogeneity Analysis – With and Without Fixed Effects",
    subtitle = "Dep. Var.: LFP | Control: 42–44 contribution years | Age 62–65"
  ) |>
  tab_spanner(label = "RQ1: East/West", id = "sp_rq1",
              columns = c("Ost/West", "Ost/West FE")) |>
  tab_spanner(label = "RQ2: Income",    id = "sp_rq2inc",
              columns = c("Einkommen", "Einkommen FE")) |>
  tab_spanner(label = "RQ2: Education", id = "sp_rq2edu",
              columns = c("Bildung", "Bildung FE")) |>
  tab_source_note("Clustered SE auf pid-Ebene") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_source_note("Individual FE: pid | Zeitkonstante Variablen werden vom FE absorbiert") |>
  tab_options(table.font.size = 14, table.width = pct(100)) |>
  gtsave(file.path(plot_dir, "interaktionsmodelle_fe.html"))
cat("Gespeichert: interaktionsmodelle_fe.html\n")


# 7. Appendix – Dynamic DiD ----------------------------------

m_dynamic <- feols(lfp ~ treat + i(syear, treat, ref = 2013) | syear,
                   data = df_main, cluster = ~pid)

coef_dynamic <- tidy(m_dynamic, conf.int = TRUE) |>
  filter(grepl("syear::", term)) |>
  mutate(
    year   = as.numeric(sub("syear::(\\d+):treat", "\\1", term)),
    sig    = ifelse(p.value < 0.05, "Significant", "Not significant"),
    period = ifelse(year < 2014, "Pre-Treatment", "Post-Treatment")
  ) |>
  bind_rows(tibble(year = 2013, estimate = 0, conf.low = 0, conf.high = 0,
                   sig = "Reference", period = "Pre-Treatment")) |>
  arrange(year)

ggplot(coef_dynamic, aes(x = year, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 2014, linetype = "dashed",
             color = "#CC7C2E", linewidth = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = sig),
              alpha = 0.15, color = NA) +
  geom_line(aes(group = 1), color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  annotate("text", x = 2014.1, y = 0.18, label = "Reform (July 2014)",
           hjust = 0, color = "#CC7C2E", size = 3.5) +
  annotate("rect", xmin = 2009, xmax = 2013, ymin = -Inf, ymax = Inf,
           alpha = 0.04, fill = "#1E2B5C") +
  annotate("rect", xmin = 2015, xmax = 2019, ymin = -Inf, ymax = Inf,
           alpha = 0.04, fill = "green") +
  annotate("text", x = 2011, y = -0.37, label = "pre-period",
           color = "#1E2B5C", size = 3) +
  annotate("text", x = 2017, y = -0.37, label = "post-period",
           color = "darkgreen", size = 3) +
  scale_color_manual(
    name   = "",
    values = c("Significant"     = "#c0392b",
               "Not significant" = "#3266ad",
               "Reference"       = "gray40")
  ) +
  scale_fill_manual(
    values = c("Significant"     = "#c0392b",
               "Not significant" = "#3266ad",
               "Reference"       = "gray40"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = c(2009:2013, 2015:2019)) +
  scale_y_continuous(limits = c(-0.40, 0.22)) +
  labs(
    title    = "Dynamic DiD - Year-by-Year Effects",
    subtitle = "Reference year: 2013 | Control: 42-44 contribution years",
    x        = "Year",
    y        = "Coefficient (Treat x Year)",
    caption  = "SOEP Data | Age 62-65 | Clustered SE at individual level | Boehm & Wilhelmi 2026"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "dynamic_did.png"), width = 10, height = 5.5, dpi = 300)
cat("Saved: dynamic_did.png\n")


# 7a. Appendix – Parallel Trends Table ----------------------

pt_table |>
  kbl(
    format    = "html",
    row.names = FALSE,
    col.names = c("Model", "Chi²", "df", "p-Value", "Result"),
    caption   = "Appendix – Parallel Trends Test (2009–2012) | Age 62–65",
    digits    = 4
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width        = FALSE,
    position          = "center",
    font_size         = 14
  ) |>
  column_spec(5,
              color = ifelse(pt_table$Ergebnis == "✓ Parallel", "green", "red"),
              bold  = TRUE
  ) |>
  add_footnote(
    "H0: Alle Pre-Treatment Koeffizienten = 0 | Clustered SE auf pid-Ebene | Sample: Age 62–65",
    notation = "none"
  ) |>
  save_kable(file.path(plot_dir, "appendix_parallel_trends.html"))
cat("Gespeichert: appendix_parallel_trends.html\n")


# 7b. Model Overview HTML ------------------------------------

modelle_html <- '<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
  body { font-family: sans-serif; font-size: 12px; padding: 2rem; }
  table { width: 100%; border-collapse: collapse; }
  th { padding: 8px 10px; text-align: left; font-weight: 500;
       background: #f5f5f5; border-bottom: 1px solid #ddd; }
  td { padding: 7px 10px; border-bottom: 0.5px solid #eee;
       font-family: monospace; font-size: 11px; }
  .section { background: #f0f0f0; font-weight: 500; font-size: 11px; color: #555; }
  .sig-yes { color: #3B6D11; font-weight: 500; }
  .sig-no  { color: #A32D2D; font-weight: 500; }
  .sig-na  { color: #888; }
  .badge-haupt { background: #E6F1FB; color: #0C447C;
                 padding: 2px 6px; border-radius: 4px; font-size: 10px; }
  .badge-aux   { background: #F1EFE8; color: #5F5E5A;
                 padding: 2px 6px; border-radius: 4px; font-size: 10px; }
</style>
</head>
<body>
<h2>Übersicht aller Regressionsmodelle – Sample: Age 62–65 | Hauptkontrollgruppe: 42–44</h2>
<table>
<thead>
  <tr><th>#</th><th>Modell</th><th>Typ</th><th>Regressionsgleichung</th><th>Sig.</th></tr>
</thead>
<tbody>
  <tr class="section"><td colspan="5">Kontrollgruppen-Auswahl (Robustness)</td></tr>
  <tr><td>1</td><td>m_35_44</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD | KG: 35–44</td><td class="sig-yes">✓</td></tr>
  <tr><td>2</td><td>m_40_44</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD | KG: 40–44</td><td class="sig-yes">✓</td></tr>
  <tr><td>3</td><td>m_41_44</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD | KG: 41–44</td><td class="sig-yes">✓</td></tr>
  <tr><td>4</td><td>m_42_44</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD | KG: 42–44</td><td class="sig-yes">✓</td></tr>
  <tr class="section"><td colspan="5">Hauptmodelle (42–44, Age 62–65)</td></tr>
  <tr><td>5</td><td>m1</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post</td><td class="sig-na">—</td></tr>
  <tr><td>6</td><td>m2</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD</td><td class="sig-yes">✓</td></tr>
  <tr><td>7</td><td>m3</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·DiD | Year FE</td><td class="sig-yes">✓</td></tr>
  <tr><td>8</td><td>m4</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Post + β₂·DiD | Individual FE</td><td class="sig-yes">✓</td></tr>
  <tr class="section"><td colspan="5">RQ1 – Regionale Heterogenität</td></tr>
  <tr><td>9</td><td>m_ow_interact</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD + β₄·East + β₅·DiD×East</td><td class="sig-no">✗</td></tr>
  <tr><td>10</td><td>m_ow_fe</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Post + β₂·DiD + β₃·DiD×East | Individual FE</td><td class="sig-no">✗</td></tr>
  <tr class="section"><td colspan="5">RQ2 – Einkommen</td></tr>
  <tr><td>11</td><td>m_inc_interact</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD + β₄·HighInc + β₅·DiD×HighInc</td><td class="sig-no">✗</td></tr>
  <tr><td>12</td><td>m_inc_fe</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Post + β₂·DiD + β₃·DiD×HighInc | Individual FE</td><td class="sig-no">✗</td></tr>
  <tr class="section"><td colspan="5">RQ2 – Bildung</td></tr>
  <tr><td>13</td><td>m_edu_interact</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + β₁·Treat + β₂·Post + β₃·DiD + β₄·HighEduc + β₅·DiD×HighEduc</td><td class="sig-no">✗</td></tr>
  <tr><td>14</td><td>m_edu_fe</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Post + β₂·DiD + β₃·DiD×HighEduc | Individual FE</td><td class="sig-no">✗</td></tr>
  <tr class="section"><td colspan="5">Langfristige Effekte</td></tr>
  <tr><td>15</td><td>m_dynamic</td><td><span class="badge-haupt">Haupt</span></td><td>LFP = β₀ + β₁·Treat + Σβₜ·(Year×Treat) | Year FE</td><td class="sig-yes">✓</td></tr>
  <tr class="section"><td colspan="5">Parallel Trends – Event Study</td></tr>
  <tr><td>16</td><td>es_35_44 – es_42_44</td><td><span class="badge-aux">Aux</span></td><td>LFP = β₀ + Σβₜ·(Year×Treat) | pid + Year FE | ref=2013</td><td class="sig-na">p&gt;0.05 ✓</td></tr>
</tbody>
</table>
<p style="font-size:11px; color:#888; margin-top:1rem;">
  Hauptkontrollgruppe: 42–44 Beitragsjahre | Sample: Age 62–65 | Begründung: bester Common Support + PT p=0.097
</p>
</body>
</html>'

writeLines(modelle_html, file.path(plot_dir, "modell_uebersicht.html"))
cat("Gespeichert: modell_uebersicht.html\n")


# 8. Robustness Checks ---------------------------------------

# 8a. Two-way FE model
m5 <- feols(lfp ~ did | pid + syear, data = df_main, cluster = ~pid)
summary(m5)

# 8b. Covariate Controls
m4_cov       <- feols(lfp ~ post + did + age | pid,
                      data = df_main, cluster = ~pid)
m_ow_fe_cov  <- feols(lfp ~ post + did + did:east + age | pid,
                      data = df_main, cluster = ~pid)
m_inc_fe_cov <- feols(lfp ~ post + did + did:high_income + age | pid,
                      data = df_main, cluster = ~pid)
m_edu_fe_cov <- feols(lfp ~ post + did + did:high_educ + age | pid,
                      data = df_main, cluster = ~pid)

summary(m4_cov)
summary(m_ow_fe_cov)
summary(m_inc_fe_cov)
summary(m_edu_fe_cov)


# 8c. Income Quintiles ---------------------------------------

df_main <- df_main |>
  mutate(
    income_quintile = ntile(income, 5),
    inc_q2 = as.integer(income_quintile == 2),
    inc_q3 = as.integer(income_quintile == 3),
    inc_q4 = as.integer(income_quintile == 4),
    inc_q5 = as.integer(income_quintile == 5)
  )

m_inc_quintile <- feols(lfp ~ post + did + did:inc_q2 + did:inc_q3 +
                          did:inc_q4 + did:inc_q5 | pid,
                        data = df_main, cluster = ~pid)
summary(m_inc_quintile)

# Compute absolute effects per quintile with delta-method SEs
V     <- vcov(m_inc_quintile)
coefs <- coef(m_inc_quintile)

se_q1 <- sqrt(V["did", "did"])
se_q2 <- sqrt(V["did","did"] + V["did:inc_q2","did:inc_q2"] + 2*V["did","did:inc_q2"])
se_q3 <- sqrt(V["did","did"] + V["did:inc_q3","did:inc_q3"] + 2*V["did","did:inc_q3"])
se_q4 <- sqrt(V["did","did"] + V["did:inc_q4","did:inc_q4"] + 2*V["did","did:inc_q4"])
se_q5 <- sqrt(V["did","did"] + V["did:inc_q5","did:inc_q5"] + 2*V["did","did:inc_q5"])

est_q1 <- coefs["did"]
est_q2 <- coefs["did"] + coefs["did:inc_q2"]
est_q3 <- coefs["did"] + coefs["did:inc_q3"]
est_q4 <- coefs["did"] + coefs["did:inc_q4"]
est_q5 <- coefs["did"] + coefs["did:inc_q5"]

p_q1 <- 2 * pnorm(abs(est_q1 / se_q1), lower.tail = FALSE)
p_q2 <- 2 * pnorm(abs(est_q2 / se_q2), lower.tail = FALSE)
p_q3 <- 2 * pnorm(abs(est_q3 / se_q3), lower.tail = FALSE)
p_q4 <- 2 * pnorm(abs(est_q4 / se_q4), lower.tail = FALSE)
p_q5 <- 2 * pnorm(abs(est_q5 / se_q5), lower.tail = FALSE)

coef_plot <- tibble(
  quintile  = 1:5,
  label     = c("Q1\n(lowest)", "Q2", "Q3", "Q4", "Q5\n(highest)"),
  estimate  = c(est_q1, est_q2, est_q3, est_q4, est_q5),
  se        = c(se_q1, se_q2, se_q3, se_q4, se_q5),
  p_value   = c(p_q1, p_q2, p_q3, p_q4, p_q5)
) |>
  mutate(
    conf.low  = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    sig_label = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
    est_label = paste0(round(estimate * 100, 1), "%")
  )

ggplot(coef_plot, aes(x = factor(quintile), y = estimate, color = sig_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.12, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_line(aes(group = 1), linewidth = 0.6, linetype = "dotted", color = "gray60") +
  geom_text(aes(label = est_label), vjust = -1.5, size = 3.2, color = "gray30") +
  scale_x_discrete(labels = coef_plot$label) +
  scale_color_manual(
    name   = "",
    values = c("p < 0.05" = "#c0392b", "p >= 0.05" = "#3266ad")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(-0.65, 0.20)
  ) +
  labs(
    title    = "Treatment Effect by Income Quintile",
    subtitle = "Absolute DiD estimates per quintile with 95% CI | Individual FE",
    x        = "Income Quintile",
    y        = "DiD Coefficient (Treat x Post)",
    caption  = "SOEP Data | Age 62-65 | Control: 42-44 contribution years | Boehm & Wilhelmi 2026"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "rq2_income_quintiles.png"),
       width = 8, height = 5.5, dpi = 300)
cat("Saved: rq2_income_quintiles.png\n")

img_quintile <- paste0("data:image/png;base64,",
                       base64encode(file.path(plot_dir, "rq2_income_quintiles.png")))

quintil_table <- coef_plot |>
  select(quintile, label, estimate, conf.low, conf.high, sig_label) |>
  mutate(
    label     = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"),
    estimate  = round(estimate * 100, 1),
    conf.low  = round(conf.low  * 100, 1),
    conf.high = round(conf.high * 100, 1)
  )

html_quintile <- paste0('
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
  body  { font-family: sans-serif; padding: 2rem; max-width: 1000px; margin: auto; }
  h1    { color: #1E2B5C; }
  h2    { color: #35746E; margin-top: 2rem; }
  table { width: 100%; border-collapse: collapse; margin-bottom: 1.5rem; }
  th    { background: #1E2B5C; color: white; padding: 10px; text-align: left; }
  td    { padding: 8px 10px; border-bottom: 1px solid #eee; }
  tr:nth-child(1) { background: #fdf0f0; font-weight: bold; }
  tr:nth-child(5) { background: #f5f5f5; }
  img      { width: 100%; border: 1px solid #eee; border-radius: 4px; margin-top: 1rem; }
  .footnote { font-size: 11px; color: #888; margin-top: 0.5rem; }
  .badge-sig { background: #fde8e8; color: #c0392b; padding: 2px 8px;
               border-radius: 10px; font-size: 11px; font-weight: bold; }
  .badge-ns  { background: #e8eef8; color: #3266ad; padding: 2px 8px;
               border-radius: 10px; font-size: 11px; }
</style>
</head>
<body>
<h1>RQ2 – Treatment Effect by Income Quintile</h1>
<p style="color:#888; font-size:13px;">
  Dep. Var.: LFP | Individual FE | Control: 42–44 contribution years |
  Age 62–65 | Böhm & Wilhelmi 2026
</p>
<h2>Regression: LFP = β₀ + β₁·Post + β₂·DiD + β₃·DiD×Q2 + β₄·DiD×Q3 + β₅·DiD×Q4 + β₆·DiD×Q5 | Individual FE</h2>
<h2>Coefficient Table</h2>
<table>
<thead>
  <tr>
    <th>Quintile</th><th>DiD Effect</th>
    <th>95% CI (low)</th><th>95% CI (high)</th><th>Significance</th>
  </tr>
</thead>
<tbody>',
                        paste0(sapply(1:5, function(i) {
                          sig_badge <- if(quintil_table$sig_label[i] == "p < 0.05") {
                            '<span class="badge-sig">p &lt; 0.05</span>'
                          } else {
                            '<span class="badge-ns">p &ge; 0.05</span>'
                          }
                          paste0('
  <tr>
    <td>', quintil_table$label[i], '</td>
    <td>', quintil_table$estimate[i], '%</td>
    <td>', quintil_table$conf.low[i], '%</td>
    <td>', quintil_table$conf.high[i], '%</td>
    <td>', sig_badge, '</td>
  </tr>')
                        }), collapse = ""),
                        '
</tbody>
</table>
<p class="footnote">
  Reference category: Q1 (lowest income quintile) | Clustered SE at pid level |
  Effects calculated as linear combinations of did + did:inc_qX
</p>
<h2>Coefficient Plot</h2>
<img src="', img_quintile, '" alt="Income Quintile Plot">
<p class="footnote">
  N = 1,674 | Individual FE (pid) | Error bars = 95% CI |
  Red = p &lt; 0.05, Blue = p &ge; 0.05
</p>
</body>
</html>
')

writeLines(html_quintile, file.path(plot_dir, "rq2_income_quintiles.html"))
cat("Gespeichert: rq2_income_quintiles.html\n")

# Pairwise tests on quintile differences (Wald + Bonferroni)
wald_test <- linearHypothesis(
  m_inc_quintile,
  c("did:inc_q2 = 0", "did:inc_q3 = 0", "did:inc_q4 = 0", "did:inc_q5 = 0")
)
print(wald_test)

pairwise_tests <- tryCatch(
  linearHypothesis(
    m_inc_quintile,
    c("did:inc_q2 = did:inc_q3", "did:inc_q2 = did:inc_q4",
      "did:inc_q2 = did:inc_q5", "did:inc_q3 = did:inc_q4",
      "did:inc_q3 = did:inc_q5", "did:inc_q4 = did:inc_q5")
  ),
  error = function(e) {
    cat("Hinweis: Paarweiser Wald-Test nicht moeglich (linear abhaengige Constraints).\n")
    NULL
  }
)
if (!is.null(pairwise_tests)) print(pairwise_tests)

V     <- vcov(m_inc_quintile)
coefs <- coef(m_inc_quintile)

quintile_pairs <- list(
  c("Q2", "Q3", "did:inc_q2", "did:inc_q3"),
  c("Q2", "Q4", "did:inc_q2", "did:inc_q4"),
  c("Q2", "Q5", "did:inc_q2", "did:inc_q5"),
  c("Q3", "Q4", "did:inc_q3", "did:inc_q4"),
  c("Q3", "Q5", "did:inc_q3", "did:inc_q5"),
  c("Q4", "Q5", "did:inc_q4", "did:inc_q5")
)

pairwise_results <- lapply(quintile_pairs, function(p) {
  diff_coef <- coefs[p[3]] - coefs[p[4]]
  diff_se   <- sqrt(V[p[3], p[3]] + V[p[4], p[4]] - 2 * V[p[3], p[4]])
  t_stat    <- diff_coef / diff_se
  p_val     <- 2 * pnorm(abs(t_stat), lower.tail = FALSE)
  tibble(
    Comparison = paste0(p[1], " vs ", p[2]),
    Difference = round(diff_coef * 100, 2),
    SE         = round(diff_se * 100, 2),
    t_stat     = round(t_stat, 3),
    p_raw      = round(p_val, 4)
  )
}) |> bind_rows() |>
  mutate(
    p_bonferroni = round(pmin(p_raw * n(), 1), 4),
    sig_raw  = case_when(p_raw < 0.001 ~ "***", p_raw < 0.01 ~ "**",
                         p_raw < 0.05  ~ "*",   TRUE          ~ "n.s."),
    sig_bonf = case_when(p_bonferroni < 0.001 ~ "***", p_bonferroni < 0.01 ~ "**",
                         p_bonferroni < 0.05  ~ "*",   TRUE                 ~ "n.s.")
  )

print(pairwise_results)

cat("\n=== Pairwise Tests: Quintile Differences (DiD Effects) ===\n")
cat("Comparison | Diff (pp) | SE  | t     | p (raw) | p (Bonf.) | sig\n")
cat(rep("-", 70), "\n", sep = "")
for (i in 1:nrow(pairwise_results)) {
  r <- pairwise_results[i, ]
  cat(sprintf("%-12s | %6.2f%%   | %4.2f | %5.3f | %7.4f | %9.4f  | %s\n",
              r$Comparison, r$Difference, r$SE, r$t_stat,
              r$p_raw, r$p_bonferroni, r$sig_bonf))
}


# 8d. Placebo Tests ------------------------------------------

run_placebo <- function(data, placebo_year) {
  df_tmp <- data |>
    mutate(
      post_placebo = case_when(
        syear >= placebo_year & syear <= 2013 ~ 1L,
        syear < placebo_year                  ~ 0L,
        TRUE                                  ~ NA_integer_
      ),
      did_placebo = treat * post_placebo
    ) |>
    filter(syear <= 2013) |>
    filter(!is.na(post_placebo))

  feols(lfp ~ post_placebo + did_placebo | pid,
        data = df_tmp, cluster = ~pid)
}

m_pl_2010 <- run_placebo(df_main, 2010)
m_pl_2011 <- run_placebo(df_main, 2011)
m_pl_2012 <- run_placebo(df_main, 2012)

summary(m_pl_2010)
summary(m_pl_2011)
summary(m_pl_2012)

placebo_models <- list(
  "Placebo 2010" = m_pl_2010,
  "Placebo 2011" = m_pl_2011,
  "Placebo 2012" = m_pl_2012
)

modelsummary(
  placebo_models,
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  coef_rename = c("post_placebo" = "Post (Placebo)",
                  "did_placebo"  = "Treat × Post (Placebo)"),
  gof_omit = "IC|Log|Adj",
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3)
  )
) |>
  tab_header(
    title    = "Placebo Tests – Fake Reform Years",
    subtitle = "Dependent variable: LFP | Sample: Age 62–65"
  ) |>
  tab_source_note("Clustered SE auf pid-Ebene | Alle Placebos nutzen nur Pre-Period (2009–2013)") |>
  gtsave(file.path(plot_dir, "placebo_tests.html"))

placebo_plot_data <- bind_rows(
  tidy(m_pl_2010) |> mutate(year = 2010),
  tidy(m_pl_2011) |> mutate(year = 2011),
  tidy(m_pl_2012) |> mutate(year = 2012)
) |>
  filter(term == "did_placebo") |>
  mutate(
    conf_low  = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error,
    signif    = ifelse(p.value < 0.05, "Significant", "Not significant")
  )

ggplot(placebo_plot_data, aes(x = factor(year), y = estimate, color = signif)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                width = 0.15, linewidth = 0.8) +
  scale_color_manual(
    values = c("Significant" = "#c0392b", "Not significant" = "#35746E"),
    name = ""
  ) +
  labs(
    title    = "Placebo Test – Estimated DiD Effects",
    subtitle = "Fake reform years (pre-2014 period only)",
    x        = "Placebo Reform Year",
    y        = "Estimated Treatment Effect (DiD)",
    caption  = "SOEP data | Age 62–65 | Clustered SE (pid)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "placebo_test_plot.png"),
       width = 7, height = 5, dpi = 300)


# 8e. Donut Design -------------------------------------------
# Exclude persons with contrib_years in [44, 45)

df_donut <- df_main |>
  filter(!(contrib_years >= 44 & contrib_years < 45))

cat("Ursprüngliches Sample:", nrow(df_main), "\n")
cat("Donut Sample:         ", nrow(df_donut), "\n")
cat("Ausgeschlossen:       ", nrow(df_main) - nrow(df_donut), "\n")

m4_donut    <- feols(lfp ~ post + did | pid,
                     data = df_donut, cluster = ~pid)
m_ow_donut  <- feols(lfp ~ post + did + did:east | pid,
                     data = df_donut, cluster = ~pid)
m_inc_donut <- feols(lfp ~ post + did + did:high_income | pid,
                     data = df_donut, cluster = ~pid)
m_edu_donut <- feols(lfp ~ post + did + did:high_educ | pid,
                     data = df_donut, cluster = ~pid)

summary(m4_donut)
summary(m_ow_donut)
summary(m_inc_donut)
summary(m_edu_donut)

coef_comparison <- bind_rows(
  tibble(
    model    = c("Main", "RQ1: East/West", "RQ2: Income", "RQ2: Education"),
    estimate = c(coef(m4)["did"], coef(m_ow_fe)["did"],
                 coef(m_inc_fe)["did"], coef(m_edu_fe)["did"]),
    se       = c(se(m4)["did"], se(m_ow_fe)["did"],
                 se(m_inc_fe)["did"], se(m_edu_fe)["did"]),
    spec     = "Baseline"
  ),
  tibble(
    model    = c("Main", "RQ1: East/West", "RQ2: Income", "RQ2: Education"),
    estimate = c(coef(m4_donut)["did"], coef(m_ow_donut)["did"],
                 coef(m_inc_donut)["did"], coef(m_edu_donut)["did"]),
    se       = c(se(m4_donut)["did"], se(m_ow_donut)["did"],
                 se(m_inc_donut)["did"], se(m_edu_donut)["did"]),
    spec     = "Donut"
  )
) |>
  mutate(
    conf.low  = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    model     = factor(model, levels = c("Main", "RQ1: East/West",
                                         "RQ2: Income", "RQ2: Education"))
  )

ggplot(coef_comparison,
       aes(x = model, y = estimate, color = spec, shape = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  scale_color_manual(name = "Specification",
                     values = c("Baseline" = "#1E2B5C", "Donut" = "#35746E")) +
  scale_shape_manual(name = "Specification",
                     values = c("Baseline" = 16, "Donut" = 17)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.55, 0.15)) +
  labs(
    title    = "Robustness Check – Donut Design",
    subtitle = "DiD coefficient (Treat × Post) | Baseline vs. excluding contrib. years 44–45",
    x        = "",
    y        = "DiD Coefficient",
    caption  = "SOEP-Data | Age 62–65 | Control: 42–44 contribution years | Böhm & Wilhelmi 2026"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "robustness_donut.png"),
       width = 9, height = 5.5, dpi = 300)
cat("Gespeichert: robustness_donut.png\n")

modelsummary(
  list(
    "Main (Baseline)"    = m4,       "Main (Donut)"       = m4_donut,
    "RQ1 (Baseline)"     = m_ow_fe,  "RQ1 (Donut)"        = m_ow_donut,
    "RQ2 Inc (Baseline)" = m_inc_fe, "RQ2 Inc (Donut)"    = m_inc_donut,
    "RQ2 Edu (Baseline)" = m_edu_fe, "RQ2 Edu (Donut)"    = m_edu_donut
  ),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c(
    "post"            = "Post",
    "did"             = "Treat × Post",
    "did:east"        = "Treat × Post × East",
    "did:high_income" = "Treat × Post × High Income",
    "did:high_educ"   = "Treat × Post × High Education"
  ),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~`Main (Baseline)`, ~`Main (Donut)`,
    ~`RQ1 (Baseline)`,  ~`RQ1 (Donut)`,
    ~`RQ2 Inc (Baseline)`, ~`RQ2 Inc (Donut)`,
    ~`RQ2 Edu (Baseline)`, ~`RQ2 Edu (Donut)`,
    "Individual FE", "Ja", "Ja", "Ja", "Ja", "Ja", "Ja", "Ja", "Ja",
    "Donut",         "Nein", "Ja", "Nein", "Ja", "Nein", "Ja", "Nein", "Ja"
  )
) |>
  tab_header(
    title    = "Robustness Check – Donut Design",
    subtitle = "Excluding contribution years 44–45 | Dep. Var.: LFP | Age 62–65"
  ) |>
  tab_spanner(label = "Main",          id = "sp1",
              columns = c("Main (Baseline)", "Main (Donut)")) |>
  tab_spanner(label = "RQ1",           id = "sp2",
              columns = c("RQ1 (Baseline)", "RQ1 (Donut)")) |>
  tab_spanner(label = "RQ2 Income",    id = "sp3",
              columns = c("RQ2 Inc (Baseline)", "RQ2 Inc (Donut)")) |>
  tab_spanner(label = "RQ2 Education", id = "sp4",
              columns = c("RQ2 Edu (Baseline)", "RQ2 Edu (Donut)")) |>
  tab_source_note("Clustered SE auf pid-Ebene | Donut: contrib. years 44–45 excluded") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_options(table.font.size = 14, table.width = pct(100)) |>
  gtsave(file.path(plot_dir, "robustness_donut.html"))
cat("Gespeichert: robustness_donut.html\n")


# 8f. Flexirente Robustness (January 2017) -------------------

df_main <- df_main |>
  mutate(flexi = as.integer(syear >= 2017))

m4_flexi <- feols(lfp ~ post + did + flexi | pid,
                  data = df_main, cluster = ~pid)
summary(m4_flexi)

df_main |>
  filter(syear %in% 2015:2019) |>
  group_by(syear, treat) |>
  summarise(lfp_mean = mean(lfp, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat, values_from = lfp_mean,
              names_prefix = "treat_") |>
  mutate(diff = treat_1 - treat_0) |>
  print()

df_preflexi <- df_main |>
  filter(syear %in% c(2009:2013, 2015:2016)) |>
  mutate(
    post = case_when(
      syear %in% 2015:2016 ~ 1L,
      syear %in% 2009:2013 ~ 0L,
      TRUE                 ~ NA_integer_
    ),
    did = treat * post
  ) |>
  filter(!is.na(post))

m4_preflexi    <- feols(lfp ~ post + did | pid,
                        data = df_preflexi, cluster = ~pid)
m_ow_preflexi  <- feols(lfp ~ post + did + did:east | pid,
                        data = df_preflexi, cluster = ~pid)
m_inc_preflexi <- feols(lfp ~ post + did + did:high_income | pid,
                        data = df_preflexi, cluster = ~pid)
m_edu_preflexi <- feols(lfp ~ post + did + did:high_educ | pid,
                        data = df_preflexi, cluster = ~pid)

summary(m4_preflexi)
summary(m_ow_preflexi)
summary(m_inc_preflexi)
summary(m_edu_preflexi)

diff_plot <- df_main |>
  filter(syear %in% 2015:2019) |>
  group_by(syear, treat) |>
  summarise(lfp_mean = mean(lfp, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat, values_from = lfp_mean,
              names_prefix = "treat_") |>
  mutate(diff = treat_1 - treat_0)

ggplot(diff_plot, aes(x = syear, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 2017, linetype = "dashed",
             color = "#c0392b", linewidth = 0.8) +
  geom_line(linewidth = 1.2, color = "#1E2B5C") +
  geom_point(size = 4, color = "#1E2B5C") +
  annotate("text", x = 2017.1, y = max(diff_plot$diff),
           label = "Flexirente (Jan. 2017)", hjust = 0, color = "#c0392b", size = 3.5) +
  geom_text(aes(label = paste0(round(diff * 100, 1), "%")),
            vjust = -1.2, size = 3.2, color = "gray30") +
  scale_x_continuous(breaks = 2015:2019) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.15, 0.10)) +
  labs(
    title    = "Treatment–Control Differential by Year (Post-Period)",
    subtitle = "No systematic break at Flexirente introduction (2017)",
    x        = "Year",
    y        = "Difference in LFP (Treatment − Control)",
    caption  = "SOEP-Data | Age 62–65 | Control: 42–44 contribution years | Böhm & Wilhelmi 2026"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "flexirente_diff_plot.png"),
       width = 8, height = 5, dpi = 300)
cat("Gespeichert: flexirente_diff_plot.png\n")

modelsummary(
  list("Baseline" = m4, "Flexi-Dummy" = m4_flexi, "2015–16 only" = m4_preflexi),
  output  = "gt",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c("post" = "Post", "did" = "Treat × Post",
                  "flexi" = "Post-2017 (Flexirente)"),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~Baseline, ~`Flexi-Dummy`, ~`2015–16 only`,
    "Individual FE", "Ja",      "Ja",            "Ja",
    "Post-Period",   "2015–19", "2015–19",       "2015–16"
  )
) |>
  tab_header(
    title    = "Robustness – Flexirente Reform (January 2017)",
    subtitle = "Dep. Var.: LFP | Control: 42–44 contribution years | Age 62–65"
  ) |>
  tab_source_note("Clustered SE auf pid-Ebene | Individual FE: pid") |>
  tab_source_note("* p < 0.05, ** p < 0.01, *** p < 0.001") |>
  tab_source_note("2015–16 only: reduced power due to shorter post-period") |>
  tab_options(table.font.size = 14, table.width = pct(70)) |>
  gtsave(file.path(plot_dir, "flexirente_robustness.html"))
cat("Gespeichert: flexirente_robustness.html\n")

img_flexi <- paste0("data:image/png;base64,",
                    base64encode(file.path(plot_dir, "flexirente_diff_plot.png")))

flexi_table_html <- modelsummary(
  list("Baseline" = m4, "Flexi-Dummy" = m4_flexi, "2015–16 only" = m4_preflexi),
  output  = "html",
  stars   = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  gof_omit = "IC|Log|Adj",
  coef_rename = c("post" = "Post", "did" = "Treat × Post",
                  "flexi" = "Post-2017 (Flexirente)"),
  gof_map = list(
    list(raw = "nobs",      clean = "Num. Obs.", fmt = 0),
    list(raw = "r.squared", clean = "R²",        fmt = 3),
    list(raw = "rmse",      clean = "RMSE",       fmt = 3)
  ),
  add_rows = tribble(
    ~term,           ~Baseline, ~`Flexi-Dummy`, ~`2015–16 only`,
    "Individual FE", "Ja",      "Ja",            "Ja",
    "Post-Period",   "2015–19", "2015–19",       "2015–16"
  )
)

html_flexi <- paste0('
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
  body      { font-family: sans-serif; padding: 2rem; max-width: 1000px; margin: auto; }
  h1        { color: #1E2B5C; }
  h2        { color: #35746E; margin-top: 2rem; }
  img       { width: 100%; border: 1px solid #eee; border-radius: 4px; margin-top: 1rem; }
  .footnote { font-size: 11px; color: #888; margin-top: 0.5rem; }
  .verdict  { background: #f0f7f0; border-left: 4px solid #35746E;
               padding: 1rem; margin: 1.5rem 0; border-radius: 4px; }
  table     { width: 100%; border-collapse: collapse; }
  th        { background: #1E2B5C; color: white; padding: 8px 10px; }
  td        { padding: 7px 10px; border-bottom: 1px solid #eee; }
</style>
</head>
<body>
<h1>Robustness – Flexirente Reform (January 2017)</h1>
<p style="color:#888; font-size:13px;">
  Testing for confounding effects of the Flexirente reform |
  Age 62–65 | Control: 42–44 contribution years | Böhm & Wilhelmi 2026
</p>
<div class="verdict">
  <strong>Verdict:</strong> No systematic confounding from the Flexirente reform detected.
  The treatment-control differential shows no structural break in 2017.
  Insignificance in the restricted sample (2015–16) is attributable to reduced statistical power.
</div>
<h2>Treatment–Control Differential by Year</h2>
<img src="', img_flexi, '" alt="Flexirente Differential Plot">
<p class="footnote">
  No systematic break visible at 2017 Flexirente introduction.
  Effect strongest in 2016, fading thereafter – consistent with dynamic DiD results.
</p>
<h2>Regression Table</h2>
', flexi_table_html, '
<p class="footnote">
  Clustered SE at pid level | Individual FE: pid |
  * p &lt; 0.05, ** p &lt; 0.01, *** p &lt; 0.001 |
  2015–16 only: reduced power due to shorter post-period
</p>
</body>
</html>
')

writeLines(html_flexi, file.path(plot_dir, "flexirente_full.html"))
cat("Gespeichert: flexirente_full.html\n")


# 9. Descriptive – Occupation Analysis -----------------------

pgen_full <- readRDS("pgen.rds") |>
  zap_labels() |>
  select(pid, syear, pgexpft, pgexppt, pgbilzeit, pgfamstd,
         pglabgro, pgemplst, pgisco08, pgstib)

df_main_occ <- df_main |>
  left_join(
    pgen_full |> select(pid, syear, pgisco08, pgstib),
    by = c("pid", "syear")
  ) |>
  filter(post == 0, treat == 1) |>
  mutate(
    income_quintile = ntile(income, 5),
    manual_job = as.integer(pgisco08 >= 6000 & pgisco08 < 10000)
  ) |>
  filter(!is.na(pgisco08), !is.na(income_quintile))

df_main_occ |>
  group_by(income_quintile) |>
  summarise(
    n            = n(),
    share_manual = round(mean(manual_job, na.rm = TRUE) * 100, 1),
    mean_isco    = round(mean(pgisco08, na.rm = TRUE), 0)
  ) |>
  print()

pgen_full2 <- readRDS("pgen.rds") |>
  zap_labels() |>
  select(pid, syear, pgexpft, pgexppt, pgbilzeit,
         pgfamstd, pglabgro, pgemplst, pgisco08,
         pgstib, pgnace, pgegp08, pgkldb2010)

df_main |>
  left_join(
    pgen_full2 |> select(pid, syear, pgstib, pgnace, pgegp08),
    by = c("pid", "syear")
  ) |>
  filter(post == 0, treat == 1, pgstib > 0) |>
  count(pgstib) |>
  arrange(pgstib) |>
  print()

df_main_stib <- df_main |>
  left_join(
    pgen_full2 |> select(pid, syear, pgstib),
    by = c("pid", "syear")
  ) |>
  filter(post == 0, treat == 1, pgstib > 0) |>
  mutate(
    income_quintile = ntile(income, 5),
    manual_worker   = as.integer(pgstib >= 10 & pgstib < 200)
  ) |>
  filter(!is.na(income_quintile))

df_main_stib |>
  group_by(income_quintile) |>
  summarise(
    n            = n(),
    n_manual     = sum(manual_worker, na.rm = TRUE),
    share_manual = round(mean(manual_worker, na.rm = TRUE) * 100, 1)
  ) |>
  print()

chisq_result <- chisq.test(
  table(df_main_stib$income_quintile, df_main_stib$manual_worker)
)
print(chisq_result)

df_main |>
  left_join(
    pgen_full2 |> select(pid, syear, pgstib),
    by = c("pid", "syear")
  ) |>
  filter(post == 0, treat == 1) |>
  summarise(
    n_total        = n(),
    n_valid_stib   = sum(pgstib > 0, na.rm = TRUE),
    n_missing_stib = sum(pgstib <= 0 | is.na(pgstib), na.rm = TRUE),
    share_valid    = round(mean(pgstib > 0, na.rm = TRUE) * 100, 1)
  ) |>
  print()


# 10. Altenquotient Germany 1950–2070 ------------------------
# Source: Destatis, 16. koordinierte Bevoelkerungsvorausberechnung (2025)

historical <- data.frame(
  year  = c(1950, 1960, 1970, 1979, 1985, 1991,
             1995, 2000, 2005, 2010, 2015, 2020, 2024),
  ratio = c(16, 18, 22, 27, 24, 24,
             25, 26, 29, 31, 34, 37, 39),
  type  = "historical"
)

projection <- data.frame(
  year       = c(2024, 2030, 2035, 2038, 2045, 2055, 2070),
  ratio_mid  = c(39,   43,   46,   47,   46,   45,   50),
  ratio_low  = c(39,   41,   43,   44,   42,   41,   43),
  ratio_high = c(39,   45,   49,   51,   53,   56,   61)
)

p <- ggplot() +
  geom_ribbon(data = projection,
              aes(x = year, ymin = ratio_low, ymax = ratio_high),
              fill = "#1E2B5C", alpha = 0.12) +
  geom_line(data = projection,
            aes(x = year, y = ratio_mid),
            linetype = "dashed", linewidth = 0.9, color = "#1E2B5C") +
  geom_line(data = historical,
            aes(x = year, y = ratio),
            linewidth = 1.2, color = "#1E2B5C") +
  geom_vline(xintercept = 2014, color = "#CC7C2E",
             linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = 2014.5, y = 57,
           label = "Rente mit 63\n(July 2014)",
           hjust = 0, color = "#CC7C2E", size = 3.2) +
  annotate("text", x = 1950,  y = 14.5, label = "16",
           size = 3, color = "#1E2B5C", hjust = 0.5) +
  annotate("text", x = 2024,  y = 37.5, label = "39",
           size = 3, color = "#1E2B5C", hjust = 1.2) +
  annotate("text", x = 2070,  y = 63,   label = "up to 61",
           size = 3, color = "#1E2B5C", hjust = 0.5) +
  annotate("text", x = 2070,  y = 41,   label = "min. 43",
           size = 3, color = "#1E2B5C", hjust = 0.5) +
  geom_vline(xintercept = 2025, color = "gray60",
             linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = 2025.5, y = 20, label = "Projection",
           hjust = 0, color = "gray50", size = 3) +
  scale_x_continuous(
    breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2014, 2024, 2038, 2070),
    labels = c("1950","1960","1970","1980","1990","2000","2010","2014","2024","2038","2070")
  ) +
  scale_y_continuous(limits = c(10, 68), breaks = seq(10, 65, by = 10),
                     labels = function(x) paste0(x)) +
  labs(
    title    = "Old-Age Dependency Ratio in Germany, 1950 to 2070",
    subtitle = paste0("Persons aged 65+ per 100 persons aged 20-64 | ",
                      "Shaded band: projection range (favourable to unfavourable scenario)"),
    x        = "Year",
    y        = "Altenquotient",
    caption  = paste0("Source: Destatis, 16. koordinierte Bevoelkerungsvorausberechnung (2025). ",
                      "Historical values: Destatis (2024). Boehm & Wilhelmi 2026")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
    plot.caption       = element_text(size = 7, color = "gray50"),
    plot.subtitle      = element_text(size = 9, color = "gray40"),
    legend.position    = "none"
  )

ggsave(file.path(plot_dir, "altenquotient_germany.pdf"),
       plot = p, width = 10, height = 5.5, device = cairo_pdf)
ggsave(file.path(plot_dir, "altenquotient_germany.png"),
       plot = p, width = 10, height = 5.5, dpi = 300)
cat("Gespeichert: altenquotient_germany.pdf / .png\n")


# 11. DiDiD – Flexirente Robustness --------------------------

df_ddd <- df_raw |>
  mutate(
    age           = syear - gebjahr,
    contrib_years = pgexpft + 0.5 * pgexppt,
    lfp = case_when(
      plb0022_h %in% c(1, 2)  ~ 1L,
      plb0022_h %in% c(3:11)  ~ 0L,
      TRUE                     ~ NA_integer_
    )
  ) |>
  filter(
    sex == 1,
    age >= 62 & age <= 65,
    syear >= 2009 & syear <= 2019,
    syear != 2014,
    !is.na(lfp)
  ) |>
  filter(contrib_years >= 45 | contrib_years < 35) |>
  mutate(
    treat = as.integer(contrib_years >= 45),
    post  = as.integer(syear >= 2015),
    flexi = as.integer(syear >= 2017),
    did   = treat * post,
    ddd   = treat * post * flexi
  )

df_ddd |> count(treat)
df_ddd |> count(treat, syear)

df_ddd <- df_ddd |>
  mutate(
    period    = case_when(
      syear <= 2013              ~ "pre",
      syear %in% c(2015, 2016)  ~ "rente63",
      syear >= 2017              ~ "flexi"
    ),
    p_rente63 = as.integer(period == "rente63"),
    p_flexi   = as.integer(period == "flexi")
  )

m_ddd <- feols(
  lfp ~ p_rente63 + p_flexi + p_rente63:treat + p_flexi:treat | pid,
  data    = df_ddd,
  cluster = ~pid
)

summary(m_ddd)


# 12. URLs & QR Codes ----------------------------------------

url_rq1_ostwest         <- "https://fan05.github.io/SPS_Winter_Regression_tables/final_rq1_ostwest.html"
url_rq2_bildung         <- "https://fan05.github.io/SPS_Winter_Regression_tables/final_rq2_bildung.html"
url_rq2_einkommen       <- "https://fan05.github.io/SPS_Winter_Regression_tables/final_rq2_einkommen.html"
url_hauptmodelle        <- "https://fan05.github.io/SPS_Winter_Regression_tables/final_hauptmodelle.html"
url_interaktionsmodelle <- "https://fan05.github.io/SPS_Winter_Regression_tables/interaktionsmodelle_fe.html"
url_modell_uebersicht   <- "https://fan05.github.io/SPS_Winter_Regression_tables/modell_uebersicht.html"
url_placebo             <- "https://fan05.github.io/SPS_Winter_Regression_tables/placebo_tests.html"
url_donut               <- "https://fan05.github.io/SPS_Winter_Regression_tables/robustness_donut.html"
url_parallel_trends     <- "https://fan05.github.io/SPS_Winter_Regression_tables/parallel_trends_full.html"
url_rq2_quintiles       <- "https://fan05.github.io/SPS_Winter_Regression_tables/rq2_income_quintiles.html"
url_flexirente          <- "https://fan05.github.io/SPS_Winter_Regression_tables/flexirente_full.html"
url_exe                 <- "https://raw.githubusercontent.com/FAN05/SPS_Winter_Regression_tables/main/SPS_Winter_Executive_Summary.pdf"

make_qr <- function(url, filename) {
  qr <- qr_code(url)
  png(file.path(plot_dir, filename), width = 300, height = 300)
  plot(qr)
  dev.off()
  cat("Gespeichert:", filename, "\n")
}

make_qr(url_donut,               "qr_donut.png")
make_qr(url_modell_uebersicht,   "qr_modell_uebersicht.png")
make_qr(url_placebo,             "qr_placebo.png")
make_qr(url_interaktionsmodelle, "qr_interaktionsmodelle.png")
make_qr(url_rq1_ostwest,         "qr_rq1_ostwest.png")
make_qr(url_rq2_bildung,         "qr_rq2_bildung.png")
make_qr(url_rq2_einkommen,       "qr_rq2_einkommen.png")
make_qr(url_hauptmodelle,        "qr_hauptmodelle.png")
make_qr(url_parallel_trends,     "qr_parallel_trends.png")
make_qr(url_rq2_quintiles,       "qr_rq2_quintiles.png")
make_qr(url_flexirente,          "qr_flexirente.png")
make_qr(url_exe,                 "qr_exe.png")



