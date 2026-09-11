# ================================================================
#  README  (English)
# ----------------------------------------------------------------
#  PURPOSE
#    Phenotypic analysis of Alternaria dauci disease scores on carrot,
#    for the phenotyping section of a paper that also reports RNA-seq
#    and metabolomics results.
#
#  DESIGN
#    Treatment (Control = Water / Treated = PRI)
#      x Temperature (NH, HS-7, HS-2, HS)
#      x Genotype (6: DEEP, NANT, NEVA, OXHELLA, PRESTO, ROBILA)
#      x Rep (4 blocks).  Disease scored at 7 dates (11-54 DAI).
#
#  RESPONSE VARIABLES
#    - AUDPC (trapezoidal)  -> PRIMARY variable
#    - D54 (final score)    -> SECONDARY variable
#    - Progression slope    -> DESCRIPTIVE ONLY (no treatment signal in
#                              the data; excluded from inference)
#
#  STATISTICAL APPROACH
#    - One linear mixed model per variable, full 3-way interaction
#      Treatment*Temperature*Genotype + (1|Rep). The SAME fit provides
#      both the type-III ANOVA and the emmeans contrasts (consistent
#      variance / df).
#    - contr.sum contrasts (required for meaningful type-III tests).
#    - Control vs Treated contrasts within each Temperature x Genotype
#      cell; raw p-values kept, plus Benjamini-Hochberg (FDR) across
#      the 24 cells. Figure stars use the FDR-corrected p-values.
#    - Compact-letter displays (Tukey) for the Temperature and Genotype
#      main effects.
#    - Repeated-measures model (AR1) is used ONLY to describe the time
#      course (progression figure), not for per-date significance.
#    - Residual normality checked (Shapiro / Anderson-Darling) -> the
#      Gaussian LMM is justified for AUDPC and D54.
#
#  BLOCK CAVEAT
#    Rep is modelled as a random block (RCBD-style). With only 4 levels
#    its variance component is imprecise; if Rep is NOT a shared
#    physical block, treat it as a fixed effect instead (see section 5).
#
#  OUTPUTS
#    figures/{pdf,png,svg}/  ->  FigA, FigB (paper candidates),
#                                FigS1-S3 (supplementary)
#    tables/                 ->  T1-T11 (summaries, ANOVAs, contrasts, CLD)
#    diagnostics/            ->  DHARMa residual diagnostics
#
#  INPUT
#    The raw phenotype workbook is configured through `PHENOTYPE_FILE` in
#    `config/config.R` (or the `CARROT_PHENOTYPE_FILE` environment variable).
#
#  DEPENDENCIES
#    See section 0; packages auto-install from CRAN if missing.
#    Note: MASS (pulled in by multcomp) masks dplyr::select(); the
#    script forces the dplyr verbs right after loading (section 0).
# ================================================================

# ================================================================
#  ANALYSE PHENOTYPIQUE — SCORES DE MALADIE (Alternaria dauci / carotte)
#  Version 2 — réoptimisée statistiquement
#
#  Design : Treatment (Control / Treated) x Temperature (NH / HS-7 / HS-2 / HS)
#           x Genotype (6) x Rep (4 blocs)   |   Dates : 11..54 DAI
#
#  Changements clés vs v1 (justifiés par les résultats obtenus) :
#   1. UN SEUL modèle par variable : interaction complète ->
#      anova(type = 3) ET emmeans issus du MEME fit (cohérence variance).
#   2. Correction de multiplicité : FDR (Benjamini-Hochberg) par famille
#      de contrastes. p brutes conservées en colonne à côté.
#   3. La PENTE (Score ~ DAI) est retirée de l'inférence : dans les
#      données, aucun effet traitement (T9/T13 : tout ns). Gardée en
#      descriptif uniquement, clairement étiquetée.
#   4. AUDPC = variable primaire ; D54 = variable secondaire.
#      Normalité vérifiée et OK (Shapiro p>0.05) -> LMM gaussien légitime.
#   5. Repeated-measures (AR1) conservé UNIQUEMENT pour décrire la
#      trajectoire temporelle (figure de progression), pas pour les
#      étoiles par date (168 tests non fiables).
# ================================================================

# ----------------------------------------------------------------
# PUBLIC-REPOSITORY INTEGRATION
# ----------------------------------------------------------------
# This is the supplied phenotype-analysis script. Statistical logic is
# preserved. Only filesystem handling was made portable for public Git use.
# Run from the repository root with:
#   Rscript scripts/phenotype/01_phenotype_analysis.R
#
# config/config.R contains generic local paths only; no HPC/cluster path is
# stored in this repository.
source(file.path("config", "config.R"))

dir.create(PHENOTYPE_OUTPUT, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(PHENOTYPE_OUTPUT)

# ----------------------------------------------------------------
# 0. PACKAGES
# ----------------------------------------------------------------
pkgs <- c("readxl","tidyverse","ggplot2","lme4","lmerTest","emmeans",
          "nlme","car","multcomp","multcompView","DHARMa","nortest",
          "patchwork","RColorBrewer","svglite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
options(stringsAsFactors = FALSE)
# contr.sum indispensable pour des type-3 ANOVA interprétables
options(contrasts = c("contr.sum", "contr.poly"))
emm_options(lmerTest.limit = 5000, pbkrtest.limit = 5000)

# --- Résolution de conflits de namespace -------------------------------
# multcomp attache MASS (via TH.data), et MASS::select() masque
# dplyr::select(). On force les verbes dplyr pour tout le script.
select    <- dplyr::select
filter    <- dplyr::filter
recode    <- dplyr::recode
summarise <- dplyr::summarise
mutate    <- dplyr::mutate

# ----------------------------------------------------------------
# 1. IMPORT + RECODAGE
# ----------------------------------------------------------------
raw <- read_excel(PHENOTYPE_FILE)
colnames(raw) <- c("Treatment","Temperature","Genotype","Rep",
                   "D11","D18","D26","D32","D40","D47","D54",
                   "RootWeight","FreshWeight","DryWeight")

score_cols <- c("D11","D18","D26","D32","D40","D47","D54")
dates_num  <- c(11, 18, 26, 32, 40, 47, 54)

raw <- raw %>%
  mutate(
    Treatment = recode(Treatment, "Water" = "Control", "PRI" = "Treated"),
    Temperature = recode(Temperature,
                         "T0" = "NH", "T1" = "HS-7", "T2" = "HS-2", "T3" = "HS"),
    Treatment   = factor(Treatment,   levels = c("Control","Treated")),
    Temperature = factor(Temperature, levels = c("NH","HS-7","HS-2","HS")),
    Genotype    = factor(Genotype),
    Rep         = factor(Rep)
  ) %>%
  select(Treatment, Temperature, Genotype, Rep, all_of(score_cols)) %>%
  mutate(Unit = interaction(Treatment, Temperature, Genotype, Rep, drop = TRUE))

long <- raw %>%
  pivot_longer(all_of(score_cols), names_to = "DAI_chr", values_to = "Score") %>%
  mutate(DAI = as.numeric(str_remove(DAI_chr, "D")),
         Score = as.numeric(Score))

# ----------------------------------------------------------------
# 2. DOSSIERS + HELPERS
# ----------------------------------------------------------------
for (d in c("figures/pdf","figures/png","figures/svg","tables","diagnostics"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

save_plot_all <- function(plot_obj, filename, width, height, dpi = 300) {
  ggsave(file.path("figures/pdf", paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path("figures/png", paste0(filename, ".png")), plot_obj,
         width = width, height = height, dpi = dpi)
  ggsave(file.path("figures/svg", paste0(filename, ".svg")), plot_obj,
         width = width, height = height, device = svglite::svglite)
}

p_to_stars <- function(p) dplyr::case_when(
  is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
  p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "ns")

# ----------------------------------------------------------------
# 3. VARIABLES DÉRIVÉES
#    AUDPC (trapèzes) = variable primaire ; D54 = variable finale.
# ----------------------------------------------------------------
audpc_trap <- function(scores, times)
  sum(diff(times) * (head(scores, -1) + tail(scores, -1)) / 2, na.rm = TRUE)

score_mat <- as.matrix(raw[, score_cols])
audpc_df <- raw %>%
  mutate(AUDPC     = apply(score_mat, 1, function(x) audpc_trap(as.numeric(x), dates_num)),
         AUDPC_rel = AUDPC / (9 * (54 - 11)) * 100)

# Pente : DESCRIPTIVE uniquement (non testée — voir en-tête, point 3)
progression_rate <- long %>%
  group_by(Treatment, Temperature, Genotype, Rep) %>%
  summarise(slope = coef(lm(Score ~ DAI))[2],
            r2    = summary(lm(Score ~ DAI))$r.squared, .groups = "drop")

# ----------------------------------------------------------------
# 4. RÉSUMÉS DESCRIPTIFS
# ----------------------------------------------------------------
score_summary <- long %>%
  group_by(Treatment, Temperature, Genotype, DAI) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE),
            se_score   = sd(Score, na.rm = TRUE)/sqrt(sum(!is.na(Score))),
            n = sum(!is.na(Score)), .groups = "drop")

audpc_summary <- audpc_df %>%
  group_by(Treatment, Temperature, Genotype) %>%
  summarise(mean_AUDPC = mean(AUDPC, na.rm = TRUE),
            se_AUDPC   = sd(AUDPC, na.rm = TRUE)/sqrt(sum(!is.na(AUDPC))),
            n = sum(!is.na(AUDPC)), .groups = "drop")

# Réduction (%) — descriptif, calculé sur moyennes de cellules
reduction_df <- audpc_summary %>%
  select(Treatment, Temperature, Genotype, mean_AUDPC) %>%
  pivot_wider(names_from = Treatment, values_from = mean_AUDPC) %>%
  mutate(Reduction_pct = (Control - Treated) / Control * 100,
         Efficacy = case_when(Reduction_pct >= 50 ~ "Effective",
                              Reduction_pct >= 25 ~ "Moderate",
                              Reduction_pct >=  0 ~ "Low",
                              TRUE ~ "Negative"))

# ================================================================
# 5. MODÈLES  —  UN SEUL FIT PAR VARIABLE (interaction complète)
#    Bloc = Rep en effet aléatoire (analyse type RCBD).
#    NB : 4 blocs -> composante de variance peu précise ; si Rep n'est
#    pas un vrai bloc physique partagé, le passer en fixe (voir doc).
# ================================================================
mod_audpc <- lmer(AUDPC ~ Treatment * Temperature * Genotype + (1 | Rep),
                  data = audpc_df, REML = TRUE)
mod_d54   <- lmer(D54   ~ Treatment * Temperature * Genotype + (1 | Rep),
                  data = raw,      REML = TRUE)

# Repeated-measures (AR1) — DESCRIPTION de la trajectoire temporelle
mod_rm <- tryCatch(
  nlme::lme(Score ~ Treatment * Temperature * Genotype + DAI,
            random = ~1 | Unit, correlation = corAR1(form = ~ DAI | Unit),
            data = long, method = "REML", na.action = na.omit),
  error = function(e) {
    message("AR1 non convergent -> modèle sans AR1.")
    nlme::lme(Score ~ Treatment * Temperature * Genotype + DAI,
              random = ~1 | Unit, data = long, method = "REML",
              na.action = na.omit)
  })

# ----------------------------------------------------------------
# 6. ANOVA TYPE 3  (issue du MEME fit que les emmeans)
# ----------------------------------------------------------------
anova_audpc <- anova(mod_audpc, type = 3)
anova_d54   <- anova(mod_d54,   type = 3)
anova_rm    <- anova(mod_rm)

# ----------------------------------------------------------------
# 7. CONTRASTES  Control vs Treated, FDR par famille
#    -> emmeans depuis le modèle d'interaction primaire (AUDPC/D54).
#    FDR appliqué sur l'ensemble des 24 cellules (correction honnête).
# ----------------------------------------------------------------
fdr_cellwise <- function(model, label) {
  emm <- emmeans(model, ~ Treatment | Temperature * Genotype)
  as.data.frame(pairs(emm, adjust = "none")) %>%
    mutate(variable = label,
           p_raw = p.value,
           p_fdr = p.adjust(p.value, method = "BH"),   # BH sur les 24 cellules
           stars_raw = p_to_stars(p_raw),
           stars_fdr = p_to_stars(p_fdr))
}
contr_audpc <- fdr_cellwise(mod_audpc, "AUDPC")
contr_d54   <- fdr_cellwise(mod_d54,   "D54")

# Effet global du traitement par température (génotypes regroupés)
emm_audpc_temp <- emmeans(mod_audpc, ~ Treatment | Temperature)
contr_audpc_temp <- as.data.frame(pairs(emm_audpc_temp, adjust = "none")) %>%
  mutate(p_fdr = p.adjust(p.value, "BH"),
         stars_fdr = p_to_stars(p_fdr))

# Effets principaux : lettres CLD (Tukey) pour Température et Génotype
cld_temp <- multcomp::cld(emmeans(mod_audpc, ~ Temperature),
                          Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))
cld_geno <- multcomp::cld(emmeans(mod_audpc, ~ Genotype),
                          Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))

# ----------------------------------------------------------------
# 8. DIAGNOSTICS  (normalité déjà OK -> on documente)
# ----------------------------------------------------------------
norm_tab <- bind_rows(
  AUDPC = { r <- residuals(mod_audpc); tibble(SW_p = shapiro.test(r)$p.value,
                                              AD_p = nortest::ad.test(r)$p.value) },
  D54   = { r <- residuals(mod_d54);   tibble(SW_p = shapiro.test(r)$p.value,
                                              AD_p = nortest::ad.test(r)$p.value) },
  .id = "Variable") %>% mutate(across(where(is.numeric), ~round(.x, 4)))

sim <- DHARMa::simulateResiduals(mod_audpc, n = 1000, plot = FALSE)
pdf("diagnostics/DHARMa_AUDPC.pdf", width = 10, height = 6)
plot(sim); dev.off()

# ----------------------------------------------------------------
# 9. THÈME + PALETTES
# ----------------------------------------------------------------
pal_trt  <- c("Control" = "#4575B4", "Treated" = "#D73027")
pal_temp <- c("NH" = "#313695", "HS-7" = "#74ADD1", "HS-2" = "#F46D43", "HS" = "#A50026")
theme_pub <- function(base = 11) theme_classic(base_size = base) +
  theme(strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# ================================================================
# 10. FIGURES
#     [PAPER]  = candidates pour l'article
#     [SUPPL]  = matériel supplémentaire
# ================================================================

# --- [PAPER] FIG A : AUDPC par génotype x température, étoiles FDR --------
stars_A <- contr_audpc %>%
  left_join(audpc_summary %>% group_by(Temperature, Genotype) %>%
              summarise(star_y = max(mean_AUDPC + se_AUDPC) + 6, .groups = "drop"),
            by = c("Temperature","Genotype")) %>%
  filter(stars_fdr != "ns")

figA <- ggplot(audpc_summary,
               aes(Genotype, mean_AUDPC, fill = Treatment)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_AUDPC - se_AUDPC, ymax = mean_AUDPC + se_AUDPC),
                position = position_dodge(0.75), width = 0.2) +
  geom_text(data = stars_A, aes(Genotype, star_y, label = stars_fdr),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ Temperature, ncol = 4) +
  scale_fill_manual(values = pal_trt) +
  labs(title = "AUDPC by genotype and temperature",
       subtitle = "Mean \u00b1 SE; stars = Control vs Treated (FDR-corrected)",
       x = "Genotype", y = "AUDPC") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_all(figA, "FigA_AUDPC_by_cell", 14, 6)

# --- [PAPER] FIG B : AUDPC global par température (effet chaleur + trt) ----
audpc_global <- audpc_df %>% group_by(Treatment, Temperature) %>%
  summarise(mean_AUDPC = mean(AUDPC), se_AUDPC = sd(AUDPC)/sqrt(n()), .groups = "drop")
stars_B <- contr_audpc_temp %>%
  left_join(audpc_global %>% group_by(Temperature) %>%
              summarise(star_y = max(mean_AUDPC + se_AUDPC) + 6, .groups = "drop"),
            by = "Temperature") %>% filter(stars_fdr != "ns")

figB <- ggplot(audpc_global, aes(Temperature, mean_AUDPC, colour = Treatment,
                                 group = Treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_AUDPC - se_AUDPC, ymax = mean_AUDPC + se_AUDPC),
                width = 0.15, linewidth = 0.8) +
  geom_text(data = stars_B, aes(Temperature, star_y, label = stars_fdr),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_colour_manual(values = pal_trt) +
  labs(title = "Global AUDPC response",
       subtitle = "Mean \u00b1 SE across genotypes; stars = Control vs Treated (FDR)",
       x = "Temperature", y = "AUDPC") +
  theme_pub()
save_plot_all(figB, "FigB_global_AUDPC", 7, 5)

# --- [SUPPL] FIG S1 : progression temporelle -----------------------------
figS1 <- ggplot(score_summary, aes(DAI, mean_score, colour = Treatment, fill = Treatment)) +
  geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_grid(Genotype ~ Temperature) +
  scale_colour_manual(values = pal_trt) + scale_fill_manual(values = pal_trt) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 9, 3)) +
  labs(title = "Disease progression over time", subtitle = "Mean \u00b1 SE",
       x = "Days after inoculation (DAI)", y = "Disease score (0\u20139)") +
  theme_pub(10)
save_plot_all(figS1, "FigS1_progression", 12, 9)

# --- [SUPPL] FIG S2 : score final D54 ------------------------------------
stars_S2 <- contr_d54 %>%
  left_join(raw %>% group_by(Temperature, Genotype) %>%
              summarise(star_y = max(D54, na.rm = TRUE) + 0.6, .groups = "drop"),
            by = c("Temperature","Genotype")) %>% filter(stars_fdr != "ns")
figS2 <- ggplot(raw, aes(Genotype, D54, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.6, alpha = 0.6) +
  geom_text(data = stars_S2, aes(Genotype, star_y, label = stars_fdr),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ Temperature, ncol = 4) + scale_fill_manual(values = pal_trt) +
  scale_y_continuous(limits = c(0, 10.5), breaks = seq(0, 9, 3)) +
  labs(title = "Final disease score at DAI 54",
       subtitle = "Stars = Control vs Treated (FDR)", x = "Genotype",
       y = "Final score (0\u20139)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_all(figS2, "FigS2_final_score_D54", 14, 6)

# --- [SUPPL] FIG S3 : réduction (%) --------------------------------------
pal_eff <- c("Effective"="#1B7837","Moderate"="#78C679","Low"="#FEC44F","Negative"="#D73027")
figS3 <- ggplot(reduction_df, aes(Genotype, Reduction_pct, fill = Efficacy)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0,  linetype = "solid",  colour = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
  geom_text(aes(y = Reduction_pct + ifelse(Reduction_pct >= 0, 2.5, -2.5),
                label = paste0(round(Reduction_pct, 1), "%")),
            size = 2.8, fontface = "bold") +
  facet_wrap(~ Temperature, ncol = 4) + scale_fill_manual(values = pal_eff) +
  labs(title = "Disease reduction due to treatment",
       subtitle = "(AUDPC Control \u2212 AUDPC Treated) / AUDPC Control \u00d7 100",
       x = "Genotype", y = "Disease reduction (%)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_all(figS3, "FigS3_reduction", 14, 6)

# ----------------------------------------------------------------
# 11. EXPORT TABLES
# ----------------------------------------------------------------
write.csv(audpc_summary,  "tables/T1_AUDPC_summary.csv", row.names = FALSE)
write.csv(reduction_df,   "tables/T2_reduction.csv",     row.names = FALSE)
write.csv(as.data.frame(anova_audpc), "tables/T3_ANOVA_type3_AUDPC.csv", row.names = TRUE)
write.csv(as.data.frame(anova_d54),   "tables/T4_ANOVA_type3_D54.csv",   row.names = TRUE)
write.csv(as.data.frame(anova_rm),    "tables/T5_repeated_scores.csv",   row.names = TRUE)
write.csv(contr_audpc,    "tables/T6_AUDPC_ControlvsTreated_FDR.csv", row.names = FALSE)
write.csv(contr_d54,      "tables/T7_D54_ControlvsTreated_FDR.csv",   row.names = FALSE)
write.csv(contr_audpc_temp,"tables/T8_AUDPC_global_by_temp_FDR.csv",  row.names = FALSE)
write.csv(cld_temp,       "tables/T9_CLD_Temperature.csv",  row.names = FALSE)
write.csv(cld_geno,       "tables/T10_CLD_Genotype.csv",    row.names = FALSE)
write.csv(norm_tab,       "tables/T11_normality.csv",       row.names = FALSE)

# ----------------------------------------------------------------
# 12. CONSOLE
# ----------------------------------------------------------------
cat("\n", strrep("=", 60), "\nANALYSE v2 — résumé\n", strrep("=", 60), "\n")
cat("\nType-3 ANOVA AUDPC :\n");        print(anova_audpc)
cat("\nContrastes AUDPC significatifs après FDR :\n")
print(contr_audpc %>% filter(p_fdr < 0.05) %>%
        select(Temperature, Genotype, estimate, p_raw, p_fdr, stars_fdr))
cat("\nCLD Température (AUDPC) :\n"); print(cld_temp[, c("Temperature",".group")])
cat("\nCLD Génotype (AUDPC) :\n");    print(cld_geno[, c("Genotype",".group")])
cat("\nNormalité (résidus) :\n");      print(norm_tab)



# ================================================================
#  SUPPLEMENTARY BLOCK — Two-phase temperature effect on disease
# ----------------------------------------------------------------
#  MESSAGE
#    Beyond the treatment (PRI) effect, TEMPERATURE alone (visible in
#    UNTREATED / Control plants) shapes disease development in two
#    opposite phases:
#      * EARLY (DAI 11-26): heat LOWERS the disease level (delayed /
#        reduced onset).
#      * LATE  (DAI 32-54): heat STEEPENS the progression rate (faster
#        late climb from a lower base).
#    Net AUDPC is the resultant of these two opposing forces.
#
#  DESIGN OF THIS BLOCK
#    - FORMAL TESTS  : Control plants only  -> isolate the pure
#                      temperature effect, independent of PRI.
#    - FIGURE        : treatments POOLED    -> full-data visual overview.
#      (CLD letters on the figure come from the Control-only models;
#       stated in the caption. Direction is identical pooled vs Control.)
#    - Late-phase slope is fitted over DAI 32-54 only, so it reflects
#      genuine late progression and is far less confounded by the
#      early baseline than a full 11-54 slope.
#
#  REQUIRES (from the main script analyse_phenotype_v2.R):
#    long, raw, audpc_df, theme_pub(), pal_temp, save_plot_all(),
#    p_to_stars(), and packages lme4/lmerTest/emmeans/multcomp/patchwork.
#    Append this block AFTER the main script, or source it once the
#    main objects exist.
# ================================================================

library(patchwork)

# ----------------------------------------------------------------
# 0. Phase definition + Control subset
# ----------------------------------------------------------------
early_dai <- c(11, 18, 26)
late_dai  <- c(32, 40, 47, 54)

long_ctrl <- long %>% filter(Treatment == "Control")

# ----------------------------------------------------------------
# 1. EARLY ONSET LEVEL  (Control) — per-plant mean score over DAI 11-26
#    LMM: Temperature * Genotype + (1|Rep)  -> type-3 ANOVA + Tukey CLD
# ----------------------------------------------------------------
early_ctrl <- long_ctrl %>%
  filter(DAI %in% early_dai) %>%
  group_by(Temperature, Genotype, Rep) %>%
  summarise(early_level = mean(Score, na.rm = TRUE), .groups = "drop")

mod_early <- lmer(early_level ~ Temperature * Genotype + (1 | Rep),
                  data = early_ctrl, REML = TRUE)
anova_early <- anova(mod_early, type = 3)
cld_early <- multcomp::cld(emmeans(mod_early, ~ Temperature),
                           Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))

# ----------------------------------------------------------------
# 2. LATE PROGRESSION RATE (Control) — per-plant slope over DAI 32-54
#    LMM: Temperature * Genotype + (1|Rep)  -> type-3 ANOVA + Tukey CLD
# ----------------------------------------------------------------
late_ctrl <- long_ctrl %>%
  filter(DAI %in% late_dai) %>%
  group_by(Temperature, Genotype, Rep) %>%
  summarise(late_rate = coef(lm(Score ~ DAI))[2], .groups = "drop")

mod_late <- lmer(late_rate ~ Temperature * Genotype + (1 | Rep),
                 data = late_ctrl, REML = TRUE)
anova_late <- anova(mod_late, type = 3)
cld_late <- multcomp::cld(emmeans(mod_late, ~ Temperature),
                          Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))

# ----------------------------------------------------------------
# 3. OMNIBUS TWO-PHASE TEST (Control) — Temperature x Phase interaction
#    A significant Temperature:Phase term = the temperature effect on
#    disease level DIFFERS between phases (the gap present early closes
#    late as heat catches up). emmeans Temperature | Phase shows it.
# ----------------------------------------------------------------
long_phase <- long_ctrl %>%
  mutate(Phase = factor(ifelse(DAI %in% early_dai, "Early", "Late"),
                        levels = c("Early", "Late")))

mod_phase <- tryCatch(
  nlme::lme(Score ~ Temperature * Phase * Genotype,
            random = ~1 | Unit,
            correlation = corAR1(form = ~ DAI | Unit),
            data = long_phase, method = "REML", na.action = na.omit),
  error = function(e) {
    message("AR1 non convergent (phase) -> sans AR1.")
    nlme::lme(Score ~ Temperature * Phase * Genotype,
              random = ~1 | Unit, data = long_phase,
              method = "REML", na.action = na.omit)
  })
anova_phase <- anova(mod_phase)
emm_phase <- emmeans(mod_phase, ~ Temperature | Phase)
cld_phase <- multcomp::cld(emm_phase, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))

# ----------------------------------------------------------------
# 4. NET AUDPC (Control) by temperature — the resultant
# ----------------------------------------------------------------
audpc_ctrl <- audpc_df %>% filter(Treatment == "Control")
mod_audpc_ctrl <- lmer(AUDPC ~ Temperature * Genotype + (1 | Rep),
                       data = audpc_ctrl, REML = TRUE)
cld_audpc_ctrl <- multcomp::cld(emmeans(mod_audpc_ctrl, ~ Temperature),
                                Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>% mutate(.group = trimws(.group))

# ================================================================
# 5. FIGURE (treatments POOLED) — 4 panels
#    A : progression curves by temperature, early/late windows shaded
#    B : early onset level by temperature   (CLD from Control model)
#    C : late progression rate by temperature (CLD from Control model)
#    D : net AUDPC by temperature            (CLD from Control model)
# ================================================================

# ---- Panel A : pooled progression curves --------------------------------
curve_pool <- long %>%
  group_by(Temperature, DAI) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE),
            se_score   = sd(Score, na.rm = TRUE)/sqrt(sum(!is.na(Score))),
            .groups = "drop")

split_dai <- 29   # frontière visuelle Early (<=26) / Late (>=32)
pA <- ggplot(curve_pool, aes(DAI, mean_score, colour = Temperature, fill = Temperature)) +
  annotate("rect", xmin = 9, xmax = split_dai, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.35) +
  annotate("text", x = mean(c(11, 26)), y = 8.6, label = "Early (11\u201326)",
           fontface = "italic", size = 3.3, colour = "grey30") +
  annotate("text", x = mean(c(32, 54)), y = 8.6, label = "Late (32\u201354)",
           fontface = "italic", size = 3.3, colour = "grey30") +
  geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.2) +
  scale_colour_manual(values = pal_temp) + scale_fill_manual(values = pal_temp) +
  scale_x_continuous(breaks = c(11, 18, 26, 32, 40, 47, 54)) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 3)) +
  labs(title = "A  Disease progression by temperature",
       subtitle = "Mean \u00b1 SE; genotypes and treatments pooled",
       x = "Days after inoculation (DAI)", y = "Disease score (0\u20139)") +
  theme_pub(11) + theme(legend.position = "top")

# ---- helper : bar panel with CLD letters (pooled bars, Control CLD) ------
bar_panel <- function(df, yvar, sevar, cld, ylab, title, digits = 0) {
  ymax <- max(df[[yvar]] + df[[sevar]], na.rm = TRUE)
  cld2 <- cld %>% transmute(Temperature, .group,
                            y = ymax * 1.08)
  ggplot(df, aes(Temperature, .data[[yvar]], fill = Temperature)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(ymin = .data[[yvar]] - .data[[sevar]],
                      ymax = .data[[yvar]] + .data[[sevar]]),
                  width = 0.2) +
    geom_text(data = cld2, aes(Temperature, y, label = .group),
              inherit.aes = FALSE, fontface = "bold", size = 4) +
    scale_fill_manual(values = pal_temp, guide = "none") +
    labs(title = title, x = "Temperature", y = ylab) +
    coord_cartesian(ylim = c(0, ymax * 1.18)) +
    theme_pub(11)
}

# Pooled descriptives for bars
early_pool <- long %>% filter(DAI %in% early_dai) %>%
  group_by(Temperature, Genotype, Rep, Treatment) %>%
  summarise(v = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(Temperature) %>%
  summarise(m = mean(v), se = sd(v)/sqrt(n()), .groups = "drop")

late_pool <- long %>% filter(DAI %in% late_dai) %>%
  group_by(Temperature, Genotype, Rep, Treatment) %>%
  summarise(v = coef(lm(Score ~ DAI))[2], .groups = "drop") %>%
  group_by(Temperature) %>%
  summarise(m = mean(v), se = sd(v)/sqrt(n()), .groups = "drop")

audpc_pool <- audpc_df %>%
  group_by(Temperature) %>%
  summarise(m = mean(AUDPC), se = sd(AUDPC)/sqrt(n()), .groups = "drop")

pB <- bar_panel(early_pool, "m", "se", cld_early,
                "Early level (0\u20139)", "B  Early onset level")
pC <- bar_panel(late_pool, "m", "se", cld_late,
                "Late rate (score/day)", "C  Late progression rate")
pD <- bar_panel(audpc_pool, "m", "se", cld_audpc_ctrl,
                "AUDPC", "D  Net AUDPC")

figS5 <- pA / (pB | pC | pD) +
  plot_annotation(
    caption = paste0(
      "Bars pooled over treatments and genotypes. Letters = Tukey groups ",
      "from Control-only models (Table S-tp). Heat lowers early onset (B) ",
      "but steepens late progression (C); net AUDPC (D) is the resultant."),
    theme = theme(plot.caption = element_text(hjust = 0, colour = "grey35"))) +
  plot_layout(heights = c(1.25, 1))

save_plot_all(figS5, "FigS5_two_phase_temperature", 12, 9)

# ----------------------------------------------------------------
# 6. EXPORT TABLES
# ----------------------------------------------------------------
write.csv(as.data.frame(anova_early), "tables/Ttp1_ANOVA_early_level_Control.csv", row.names = TRUE)
write.csv(as.data.frame(anova_late),  "tables/Ttp2_ANOVA_late_rate_Control.csv",  row.names = TRUE)
write.csv(as.data.frame(anova_phase), "tables/Ttp3_ANOVA_TemperatureByPhase_Control.csv", row.names = TRUE)
write.csv(cld_early,      "tables/Ttp4_CLD_early_level.csv",  row.names = FALSE)
write.csv(cld_late,       "tables/Ttp5_CLD_late_rate.csv",    row.names = FALSE)
write.csv(cld_phase,      "tables/Ttp6_CLD_TemperatureByPhase.csv", row.names = FALSE)
write.csv(cld_audpc_ctrl, "tables/Ttp7_CLD_AUDPC_Control.csv", row.names = FALSE)

# ----------------------------------------------------------------
# 7. CONSOLE
# ----------------------------------------------------------------
cat("\n", strrep("=", 60), "\nTWO-PHASE TEMPERATURE EFFECT (Control only)\n",
    strrep("=", 60), "\n")
cat("\n[Early onset level] Temperature effect (type-3):\n"); print(anova_early)
cat("Tukey groups:\n"); print(cld_early[, c("Temperature", "emmean", ".group")])
cat("\n[Late progression rate] Temperature effect (type-3):\n"); print(anova_late)
cat("Tukey groups:\n"); print(cld_late[, c("Temperature", "emmean", ".group")])
cat("\n[Temperature x Phase interaction]:\n"); print(anova_phase)
cat("emmeans by phase (gap early -> closing late):\n")
print(cld_phase[, c("Phase", "Temperature", "emmean", ".group")])