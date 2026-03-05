# Parameters

Registry from `fortran/Databases/Parameters.txt`.

Columns:
- `name`: FABM parameter key used in code
- `default`: default value in the registry

| name             | default   | unit                    | description                                                             | group   |
|:-----------------|:----------|:------------------------|:------------------------------------------------------------------------|:--------|
| q10_gel          | 3.5       | -                       | Temperature dependency for GEL (Kremer, 1977)                           | scale   |
| q10_phy          | 2.0       | -                       | Temperature factor (Soetaert et al., 2001)                              | scale   |
| q10_dia          | 1.8       | -                       | Temperature factor for DI                                               | scale   |
| q10_che          | 2.0       | -                       | Temperature factor for chemical processes                               | scale   |
| q10_bac          | 2.0       | -                       | Temperature factor for BAC                                              | scale   |
| q10_si_diss      | 3.3       | -                       | Temperature factor for chemical processes                               | scale   |
| q10_zoo          | 2.0       | -                       | Temperature factor Soetart et al., 2001                                 | scale   |
| eff_mic_fla      | 0.0       | -                       | Capture efficiency of MIC on FL                                         | effic   |
| eff_mic_emi      | 1.0       | -                       | Capture efficiency of MIC on EM                                         | effic   |
| eff_mic_dia      | 0.0       | -                       | Capture efficiency of MIC on DI                                         | effic   |
| eff_mic_mic      | 0.0       | -                       | Capture efficiency of MIC on MIC                                        | effic   |
| eff_mic_mes      | 0.0       | -                       | Capture efficiency of MIC on MES                                        | effic   |
| eff_mic_pom      | 0.0       | -                       | Capture efficiency of MIC on POM                                        | effic   |
| eff_mic_bac      | 0.7       | -                       | Capture efficiency of MIC on BAC                                        | effic   |
| eff_mes_fla      | 0.4       | -                       | Capture efficiency of MES on FL                                         | effic   |
| eff_mes_emi      | 0.4       | -                       | Capture efficiency of MES on EM                                         | effic   |
| eff_mes_dia      | 1.0       | -                       | Capture efficiency of MES on DI                                         | effic   |
| eff_mes_mic      | 1.0       | -                       | Capture efficiency of MES on MIC                                        | effic   |
| eff_mes_mes      | 0.0       | -                       | Capture efficiency of MES on MES                                        | effic   |
| eff_mes_pom      | 0.8       | -                       | Capture efficiency of MES on POM                                        | effic   |
| eff_mes_bac      | 0.0       | -                       | Capture efficiency of MES on BAC                                        | effic   |
| eff_noc_fla      | 0.5       | -                       | Capture efficiency of NOC on FL                                         | effic   |
| eff_noc_emi      | 1.0       | -                       | Capture efficiency of NOC on EM                                         | effic   |
| eff_noc_fla      | 1.0       | -                       | Capture efficiency of NOC on DI                                         | effic   |
| eff_noc_mic      | 1.0       | -                       | Capture efficiency of NOC on MIC                                        | effic   |
| eff_noc_mes      | 0.0       | -                       | Capture efficiency of NOC on MES                                        | effic   |
| eff_noc_pom      | 1.0       | -                       | Capture efficiency of NOC on POM                                        | effic   |
| eff_gel_fla      | 0.0       | -                       | Capture efficiency of GEL on FL                                         | effic   |
| eff_gel_emi      | 0.0       | -                       | Capture efficiency of GEL on EM                                         | effic   |
| eff_gel_dia      | 0.0       | -                       | Capture efficiency of GEL on DI                                         | effic   |
| eff_gel_mic      | 0.0       | -                       | Capture efficiency of GEL on MIC                                        | effic   |
| eff_gel_mes      | 1.0       | -                       | Capture efficiency of GEL on MES                                        | effic   |
| eff_gel_pom      | 0.0       | -                       | Capture efficiency of GEL on POM                                        | effic   |
| eff_ass_zoo_n    | 0.77      | -                       | ZOO assimilation efficiencies on N (Anderson and Pondhaven, 2003)       | effic   |
| eff_ass_zoo_c    | 0.64      | -                       | ZOO assimilation efficiency on C (Anderson and Pondhaven, 2003)         | effic   |
| eff_ass_noc_prey | 0.75      | -                       | NOC assimilation efficiency on prey (Lancelot et al., 2002)             | effic   |
| eff_ass_gel_prey | 0.75      | -                       | GEL assimilation efficiency on prey (Lancelot et al., 2002)             | effic   |
| eff_gr_mic_c     | 0.8       | -                       | MIC net growth efficiency on C (Anderson and Pondhaven, 2003)           | effic   |
| eff_gr_mes_c     | 0.8       | -                       | MES net growth efficiency on C                                          | effic   |
| eff_gr_bac_c     | 0.17      | -                       | BAC gross growth efficiency on C (Anderson and Pondhaven, 2003)         | effic   |
| eff_gr_gel_c     | 0.2       | -                       | Part of the assimil. food used for GEL growth (Lancelot et al., 2002)   | effic   |
| eff_gr_noc_c     | 0.2       | -                       | Part of the assimil. food used for GEL growth (Lancelot et al., 2002)   | effic   |
| rmin_n_c_fla     | 0.05      | molN molC-1             | Minimum N:C ratio in FL (Soetaert et al., 2001)                         | ratio   |
| rmin_n_c_emi     | 0.05      | molN molC-1             | Minimum N:C ratio in EM (Soetaert et al., 2001)                         | ratio   |
| rmin_n_c_dia     | 0.05      | molN molC-1             | Minimum N:C ratio in DI (Soetaert et al., 2001)                         | ratio   |
| rmax_n_c_fla     | 0.2       | mol N molC-1            | Maximum N:C ratio in FL (Soetaert et al., 2001)                         | ratio   |
| rmax_n_c_emi     | 0.2       | molN molC-1             | Maximum N:C ratio in EM (Soetaert et al., 2001)                         | ratio   |
| rmax_n_c_dia     | 0.2       | molN molC-1             | Maximum N:C ratio in DI (Soetaert et al., 2001)                         | ratio   |
| rmin_chl_n_fla   | 1.0       | g Chla molN-1           | Minimum Chl:N ratio in FL (Soetaert et al., 2001)                       | ratio   |
| rmin_chl_n_emi   | 1.0       | g Chla molN-1           | Minimum Chl:N ratio in EM (Soetaert et al., 2001)                       | ratio   |
| rmin_chl_n_dia   | 1.0       | g Chla molN-1           | Minimum Chl:N ratio in DI (Soetaert et al., 2001)                       | ratio   |
| rmax_chl_n_fla   | 2.0       | g Chla molN-1           | Maximum Chl:N ratio in FL (Soetaert et al., 2001)                       | ratio   |
| rmax_chl_n_emi   | 2.0       | g Chla molN-1           | Maximum Chl:N ratio in EM (Soetaert et al., 2001)                       | ratio   |
| rmax_chl_n_dia   | 2.0       | g Chla molN-1           | Maximum Chl:N ratio in DI (Soetaert et al., 2001)                       | ratio   |
| r_n_c_noc        | 0.21      | molN molC-1             | N:C molar ratio in NOC (Nakamura, 1998, JPR)                            | ratio   |
| r_n_c_gel        | 0.25      | molN molC-1             | N:C molar ratio in GEL                                                  | ratio   |
| r_n_c_mic        | 0.18      | molN molC-1             | N:C molar ratio in MIC (Anderson and Pondhaven, 2003)                   | ratio   |
| r_n_c_mes        | 0.21      | molN molC-1             | N:C molar ratio in MES                                                  | ratio   |
| r_n_c_bac        | 0.196     | molN molC-1             | N:C (Goldman) ratio in BAC (Anderson and Pondhaven, 2003)               | ratio   |
| r_si_n_dia       | 0.83      | molSi molN-1            | Si:N ratio in DI (Aksnes et al. 1994)                                   | ratio   |
| r_p_n_redfield   | 0.0625    | molP molN-1             | N:P Redfield ratio in PHY                                               | ratio   |
| umax_nos_fla     | 0.50      | molN molC-1 d-1         | Maximal NOS uptake rate by FL                                           | urate   |
| umax_nos_emi     | 1.5       | molN molC-1 d-1         | Maximal NOS uptake rate by EM                                           | urate   |
| umax_nos_dia     | 1.0       | molN molC-1 d-1         | Maximal NOS uptake rate by DI                                           | urate   |
| umax_nhs_fla     | 0.5       | molN molC-1 d-1         | Maximal NHS uptake rate by FL                                           | urate   |
| umax_nhs_emi     | 1.5       | molN molC-1 d-1         | Maximal NHS uptake rate by EM                                           | urate   |
| umax_nhs_dia     | 1.0       | molN molC-1 d-1         | Maximal NHS uptake rate by DI                                           | urate   |
| umax_si_dia      | 0.5       | molSi molC-1 d-1        | Maximal SiOs uptake rate by DI                                          | urate   |
| umax_po4_emi     | 0.09375   | molP molC-1 d-1         | Maximal PO4 uptake rate by EM                                           | urate   |
| umax_po4_dia     | 0.0625    | molP molC-1 d-1         | Maximal PO4 uptake rate by DI                                           | urate   |
| umax_po4_fla     | 0.03125   | molP molC-1 d-1         | Maximal PO4 uptake rate by FL                                           | urate   |
| ki_nhs_phy       | 0.5       | mmolN m-3               | Inhib. constant of NHS for NOS uptake by PHY (Soetaert et al., 2001)    | ihfcs   |
| ks_nos_fla       | 3.0       | mmolN m-3               | Half-saturation constant for NOS uptake by FL                           | uhfcs   |
| ks_nos_emi       | 0.05      | mmolN m-3               | Half-saturation constant for NOS uptake by EM                           | uhfcs   |
| ks_nos_dia       | 1.0       | mmolN m-3               | Half-saturation constant for NOS uptake by DI                           | uhfcs   |
| ks_nhs_fla       | 3.0       | mmolN m-3               | Half-saturation constant for NHS uptake by FL                           | uhfcs   |
| ks_nhs_emi       | 0.05      | mmolN m-3               | Half-saturation constant for NHS uptake by EM                           | uhfcs   |
| ks_nhs_dia       | 1.0       | mmolN m-3               | Half-saturation constant for NHS uptake by DI                           | uhfcs   |
| ks_po4_fla       | 0.2       | mmolP m-3               | Half-saturation constant for PO4 uptake by FL                           | uhfcs   |
| ks_po4_dia       | 0.1       | mmolP m-3               | Half-saturation constant for PO4 uptake by DI                           | uhfcs   |
| ks_po4_emi       | 0.02      | mmolP m-3               | Half-saturation constant for PO4 uptake by EM                           | uhfcs   |
| ks_sio_dia       | 3.5       | mmolSi m-3              | Half-saturation constant for SiOs uptake by DI (Paasche, 1980)          | uhfcs   |
| ks_dsc_bac       | 417.0     | mmolC m-3               | Half-sat. constant for DSC uptake by BAC (Anderson and Pondhaven, 2003) | uhfcs   |
| ks_dls_bac       | 25.0      | mmolC m-3               | Half-sat. constant for DLC uptake by BAC (Anderson and Pondhaven, 2003) | uhfcs   |
| ks_nhs_bac       | 0.5       | mmolN m-3               | Half-sat. constant for NHS uptake by BAC (Anderson and Pondhaven, 2003) | uhfcs   |
| ks_po4_bac       | 0.031     | mmolP m-3               | Half-saturation constant for PO4 uptake by BAC                          | uhfcs   |
| ks_prey_mic      | 5.0       | mmolC m-3               | Half-saturation constant for MIC grazing (Soetart et al., 2001)         | ghfcs   |
| ks_prey_mec      | 5.0       | mmolC m-3               | Half-saturation constant for MEC grazing (Soetart et al., 2001)         | ghfcs   |
| ks_mort_mic      | 1.0       | mmolC m-3               | (?) Mortality half-saturation rate of MIC                               | mhfcs   |
| ks_mort_mes      | 1.0       | mmolC m-3               | (?) Mortality half-saturation rate of MES                               | mhfcs   |
| ks_mort_noc      | 0.0       | mmolC m-3               | (?) Mortality half-saturation rate of NOC                               | mhfcs   |
| ks_mort_gel      | 0.0       | mmolC m-3               | (?) Mortality half-saturation rate of GEL                               | mhfcs   |
| ks_hydr_o2       | 2.7       | mmolO2 m-3              | Half-saturation constant for oxic hydrolysis rate                       | ohfcs   |
| ks_oxic_o2       | 3.0       | mmolO2 m-3              | Half-sat. constant for O2 lim. in oxic min. (Soetaert et al., 1996)     | chfcs   |
| ks_denitr_nos    | 0.3       | mmolN m-3               | Half-sat. constant for NOS lim. in denitrif. (Soetaert et al., 1996)    | chfcs   |
| ks_odu_iron      | 100.0     | mmolFe m-3              | Half-sat. constant for iron lim. in solid ODU formation                 | chfcs   |
| ks_odu_nos       | 2.0       | mmolN m-3               | Half-sat. constant for NOS lim. in ODU oxidation by NOS                 | chfcs   |
| ks_nhs_o2        | 3.0       | mmolO2 m-3              | Half-sat. constant for O2 lim. in NHS oxidation by O2                   | chfcs   |
| ks_odu_o2        | 1.0       | mmolO2 m-3              | Half-sat. constant for O2 lim. in ODU oxidation (Soetaert et al., 1996) | chfcs   |
| ki_anox_o2       | 0.0005    | mmolO2 m-3              | Half-sat. constant for O2 inhibition in anoxic remineralization         | ihfcs   |
| ki_anox_nos      | 0.0005    | mmolN m-3               | Half-sat. constant for NOS inhibition in anoxic remineralization        | ihfcs   |
| ki_nhs_o2        | 8.0       | mmolO2 m-3              | Half-sat. constant for O2 inhibition in NHS oxidation by NOS            | ihfcs   |
| ki_denit_o2      | 0.5       | mmolO2 m-3              | Half-sat. constant for O2 inhibition in denitrification                 | ihfcs   |
| ki_odu_o2        | 5.0       | mmolO2 m-3              | Half-sat. constant for O2 inhibition in ODU oxidation by NOS            | ihfcs   |
| ki_nhs_odu       | 0.5       | mmolO2 m-3              | Half-sat. constant for O2 inhibition in NHS oxidation by ODU            | ihfcs   |
| pi_fla           | 0.2153    | m2 W-1 d-1              | Initial slope of photosynthesis-light curve for FL                      | light   |
| pi_emi           | 0.3       | m2 W-1 d-1              | Initial slope of photosynthesis-light curve for EM                      | light   |
| pi_dia           | 0.3312    | m2 W-1 d-1              | Initial slope of photosynthesis-light curve for DI                      | light   |
| mo_fla           | 0.03      | d-1                     | Mortality rate of FL (Asknes et al., 1994)                              | mrate   |
| mo_emi           | 0.03      | d-1                     | Mortality rate of EM (Asknes et al., 1994)                              | mrate   |
| mo_dia           | 0.03      | d-1                     | Mortality rate of DI (Asknes et al., 1994)                              | mrate   |
| momax_mic        | 0.3       | d-1                     | Maximum mortality rate of MIC (Anderson and William)                    | mrate   |
| momax_mes        | 0.3       | d-1                     | Maximum mortality rate of MES (Anderson and William)                    | mrate   |
| momax_noc        | 0.06      | d-1                     | Maximum mortality rate of NOC (Lancelot et al., 2002)                   | mrate   |
| momax_gel        | 0.009     | d-1                     | Maximum mortality rate of GEL (Lancelot et al., 2002)                   | mrate   |
| moexp_mic        | 2.0       | -                       | Order of the non-linearity of mortality rate for MIC                    | order   |
| moexp_mes        | 2.0       | -                       | Order of the non-linearity of mortality rate for MES                    | order   |
| mo_bac           | 0.05      | d-1                     | Bacteria natural mortality (Anderson and Pondhaven, 2003)               | mrate   |
| mo_anox_pred     | 0.25      | d-1                     | Mortality rate in anoxia                                                | mrate   |
| mumax_fla        | 1.0       | d-1                     | Maximum specific growth rate of FL                                      | grate   |
| mumax_emi        | 2.5       | d-1                     | Maximum specific growth rate of EM                                      | grate   |
| mumax_dia        | 3.5       | d-1                     | Maximum specific growth rate of DI                                      | grate   |
| mumax_bac        | 0.000154  | d-1                     | Maximum labile DOC or NHS uptake by BAC (A&P, 2003)                     | urate   |
| mess_prey_mic    | 0.23      | -                       | Messy feeding fraction of MIC grazing (Anderson and Pondhaven, 2003)    | effic   |
| mess_prey_mes    | 0.23      | -                       | Messy feeding fraction of MES grazing                                   | effic   |
| hmax_poc         | 0.04      | d-1                     | POC hydrolysis rate (Anderson and Pondhaven, 2003)                      | hrate   |
| hmax_pon         | 0.055     | d-1                     | PON hydrolysis rate (Anderson and Pondhaven, 2003)                      | hrate   |
| hmax_dsl         | 4.0       | d-1                     | Maximum DSL hydrolysis (Anderson and Pondhaven, 2003)                   | hrate   |
| hmax_sid         | 0.08      | d-1                     | Rate of dissolution of silicious detritus (Tusseau, 1996)               | hrate   |
| gmax_mic         | 3.6       | d-1                     | Maximum grazing rate of MIC (Strom and Morello, 1998)                   | grate   |
| gmax_mes         | 1.2       | d-1                     | Maximum grazing rate of MES                                             | grate   |
| gmax_gel         | 0.3       | d-1                     | Maximum grazing rate of GEL                                             | grate   |
| gmax_noc         | 0.06      | d-1                     | Maximum grazing rate of NOC                                             | grate   |
| respb_fla        | 0.009     | d-1                     | Basal respiration rate of FL (Lancelot et al., 2002)                    | rrate   |
| respb_emi        | 0.009     | d-1                     | Basal respiration rate of EM (Lancelot et al., 2002)                    | rrate   |
| respb_dia        | 0.009     | d-1                     | Basal respiration rate of DI (Lancelot et al., 2002)                    | rrate   |
| respb_noc        | 0.0001    | d-1                     | Basal respiration rate of NOC                                           | rrate   |
| respb_gel        | 0.0001    | d-1                     | Basal respiration rate of GEL                                           | rrate   |
| f_pp_resp_emi    | 0.1       | -                       | Part of primary production used for respiration by EM                   | effic   |
| f_pp_resp_dia    | 0.1       | -                       | Part of primary production used for respiration by DI                   | effic   |
| f_pp_resp_fla    | 0.1       | -                       | Part of primary production used for respiration by FL                   | effic   |
| iron             | 10.0      | mmolFe m-3              | Concentration of iron in surface water                                  | conc    |
| i1_curve         | 25000.0   | -                       | Parameter of the curve simulating the iron concentration                | scale   |
| i2_curve         | 50.0      | -                       | Parameter of the curve simulating the iron c                            | scale   |
| f_leak_phy       | 0.02      | -                       | Phytoplankton leakage fraction (Van der Molen et al., 2004)             | effic   |
| f_dl_dom         | 0.7       | -                       | Labile fraction of PHY- and nonPHY-produced DOM (A&P, 2003)             | effic   |
| f_dl_phy_mo      | 0.34      | -                       | DOM fraction of phytoplankton mortality                                 | effic   |
| f_dl_phy_ex      | 0.65      | -                       | Labile fraction phytoxcreted DOC (Anderson and Pondhaven, 2003)         | effic   |
| exc_extra_doc    | 0.05      | -                       | Extra-photosynthetic DOC excretion (Van der Molen et al, 2004)          | effic   |
| doxsatmort       | 7.8125    | mmolO2 m-3 (?)          | Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.  | ?       |
| t_g_gel          | 0         | mmolC m-3               | Feeding threshold for GEL grazing (Lancelot et al., 2002)               | ?       |
| t_g_noc          | ?         | mmolC m-3               | Feeding threshold for NOC grazing                                       | ?       |
| rox_nhs_o2       | 0.03      | d-1                     | Maximum NHS oxidation rate of NHS by O2                                 | orate   |
| rox_nhs_nos      | 0.05      | d-1                     | Maximum NHS oxidation rate by NOS                                       | orate   |
| rox_odu_o2       | 0.1       | d-1                     | Maximum ODU oxidation rate by O2 (Oguz et al. 2000)                     | orate   |
| rox_odu_nos      | 0.05      | d-1                     | Maximum ODU oxidation rate by NOS                                       | orate   |
| f_solid_odu      | 0.2       | -                       | Percentage of solid ODU formation                                       | effic   |
| r_n_c_denit      | 0.8       | molN molC-1             | N:C ratio of denitrification                                            | ratio   |
| r_o2_c_resp      | 1.0       | molO2 molC-1            | O2:C ratio of respiration process                                       | ratio   |
| r_odu_c_anox     | 1.0       | molODU molC-1           | ODU:C ratio in anoxic remineralization                                  | ratio   |
| r_o2_nhs_nitr    | 2.0       | molO2 molNS-1           | O2:NHS ratio in NHS oxidation in nitrification                          | ratio   |
| r_o2_odu_oxid    | 1.0       | molO2 molODU-1          | O2:ODU ratio in ODU oxidation                                           | ratio   |
| r_nos_odu_oxid   | 0.8       | molNOS molODU-1         | NOS:ODU ratio in ODU oxidation                                          | ratio   |
| r_nos_nhs_oxid   | 0.6       | molNOS molNHS-1         | NOS:NHS ratio in NHS oxidation                                          | ratio   |
| w_dia            | -1.0      | m d-1                   | Sinking velocity of DI                                                  | speed   |
| w_emi            | -1.5      | m d-1                   | Sinking velocity of EM                                                  | speed   |
| w_fla            | -2.0      | m d-1                   | Sinking velocity of FL                                                  | speed   |
| w_dom            | -1.0      | m d-1                   | Sinking velocity of dissolved detritus                                  | speed   |
| w_pom            | -2.0      | m d-1                   | Sinking velocity of particulate detritus                                | speed   |
| w_sid            | -2.0      | m d-1                   | Sinking velocity of silicious detritus                                  | speed   |
| k_d              | 0.03      | m-1                     | Background light attanuation coefficient                                | light   |
| dr_sid           | 5.5e-06   | mol d-1                 | Deposition rate of silicious detritus                                   | drate   |
| dr_dom           | 5.5e-06   | mol d-1                 | Deposition rate of dissolved detritus                                   | drate   |
| dr_pom           | 5.6e-06   | mol d-1                 | Deposition rate of particulate detritus                                 | drate   |
| dr_emi           | 5.7e-06   | mol d-1                 | Deposition rate of EM                                                   | drate   |
| dr_fla           | 5.8e-06   | mol d-1                 | Deposition rate of FL                                                   | drate   |
| dr_fla           | 5.9e-06   | mol d-1                 | Deposition rate of DI                                                   | drate   |
| ks_nhs_nos       | 0.3       | mmolN m-3               | Half-sat. constant for NOS lim. in NHS oxidation by NOS                 | obsol   |
| q10_fla          | 1.8       | -                       | Temperature factor for FL                                               | obsol   |
| q10_emi          | 1.8       | -                       | Temperature factor for EM                                               | obsol   |
|                  | 0.28      | -                       | Percentage of the sediment C flux which is burried                      | effic   |
|                  | 0.28      | -                       | Percentage of the sediment N flux which is burried                      | effic   |
|                  | 0.1       | -                       | Respirated fraction (linked to activity) Soetaert et al., 2001          | obsol   |
|                  | 0.1       | -                       | Max respirated fraction                                                 | obsol   |
|                  | 0.8       | mmolC (mg Chl dW m-2)-1 | Maximum quantum yield of DI (Soetaert et al., 2001)                     | obsol   |
|                  | 0.6       | mmolC (mg Chl dW m-2)-1 | Maximum quantum yield of FL (Soetaert et al., 2001)                     | obsol   |
|                  | 0.6       | mmolC (mg Chl dW m-2)-1 | Maximum quantum yield of EM (Soetaert et al., 2001)                     | obsol   |
|                  | 2.0       | -                       | Quadratic mortality of NOC                                              | obsol   |
|                  | 2.0       | -                       | Quadratic mortality of GEL                                              | obsol   |
|                  | 0.8       | ?                       | average downward cosine                                                 | other   |
|                  | 0.84      | ?                       | Shear rate (Kriest et al., 2002)                                        | other   |
|                  | 0.62      | -                       | Sinking exponent (Kriest et al., 2002)                                  | other   |
|                  | 2294      | m^(etabio-1) d-1        | Sinking factor (Kriest et al., 2002)                                    | other   |
|                  | 0.55      | -                       | Stickness (Kriest et al., 2002)                                         | other   |
|                  | 2e-05     | m                       | Minimal cell size (Kriest et al., 2002)                                 | other   |
|                  | 0.015     | m                       | Maximum cell size (Kriest et al., 2002)                                 | other   |
|                  | 1.62      | ?                       | N content exponent (Kriest et al., 2002)                                | other   |
|                  | 0.4744    | mmolN (m^dzetabio)-1    | N content coefficient (Kriest et al., 2002)                             | other   |
|                  | 3.5       | (?)                     | Epsilon initial value (Kriest et al., 2002)                             | other   |
|                  | 0.29      | -                       | Fraction of fast C flux                                                 | other   |
|                  | 0.5       | -                       | Fraction of slow Si flux                                                | other   |
|                  | 0.0753    | d-1                     | Degradation rate of the fast degrated C                                 | other   |
|                  | 0.003     | d-1                     | Degradation rate of the slow degrated C                                 | other   |
|                  | 0.014     | d-1                     | Dissolution rate of the fast degrated Si                                | other   |
|                  | 0.0014    | d-1                     | Dissolution rate of the slow degrated Si                                | other   |
|                  | 0.11      | -                       | Fraction of precipiting ODU                                             | bent    |
|                  | 1.0       | molO molC-1             | Oxygen consumption for POM degradation                                  | other   |
|                  | 2.0       | molO molN-1             | Oxygen consumption for nitrification                                    | other   |
|                  | 0.15      | molO molN-1             | Oxygen consumption for nitrification                                    | other   |
|                  | 2.0       | -                       | Temperature factor on carbon remineralization                           | bent    |
|                  | 1e-05     | m d-1                   | (?) Diatom sunking velocity                                             | other   |
|                  | 0.05      | N m-2                   | Bottom stress threshold for deposition of Si                            | bent    |
|                  | 0.05      | N m-2                   | Bottom stress threshold for erosion of Si                               | bent    |
|                  | 0.00037   | g (m2s)-1               | Erosion constant for mineral particles                                  | other   |
|                  | 1e-13     | s-1                     | ?                                                                       | other   |
|                  | 1.0       | ?                       | ?                                                                       | other   |
|                  | 1.0       | ?                       | ?                                                                       | other   |
|                  | 287.355   | ?                       | ?                                                                       | other   |
|                  | 0.54      | ?                       | ?                                                                       | other   |
|                  | 0.63      | ?                       | ? Part of PAR with a long wavelegth                                     | other   |
|                  | 4.0       | ?                       | ?                                                                       | other   |
|                  | 0.23      | ? m-1                   | ? Background attenuation coefficient for longwave radiation             | other   |
|                  | 0.02      | ?                       | ?                                                                       | other   |
|                  | 0.03      | (mg Chl a m2)-1         | Self-shading extinction coefficient by chlorophyll                      | other   |
|                  | 0.03      | (mmolC m2)-1            | Self-shading extinction coefficient by detritus                         | other   |
|                  | 0.0196    | m-1                     | Absorption by pure sea water (Smith and Baker, 1981)                    | other   |
|                  | 0.0015    | m-1                     | Backscattering by pure sea water (Smith and Baker, 1981)                | other   |
|                  | 0.029     | m-1                     | A*CHL^B (Bricaud et al., 1995, Dmitriev et al., 2007)                   | other   |
|                  | 0.6       | m-1                     | A*CHL^B (Bricaud et al., 1995, Dmitriev et al., 2007)                   | other   |
|                  | 5.8e-06   | m2 (mg POC_DIA)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 8.17e-06  | m2 (mg POC_FLA)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 1.008e-05 | m2 (mg POC_EMI)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 0.01      | m2 g-1                  | Absorption by suspended minerals (Neuckermans et al., 2012)             | other   |
|                  | 0.0155    | m-1                     | Backscattering by suspended minerals (Neuckermans et al., 2012)         | other   |
|                  | 0.05      | m2 g-1                  | Absorption by organic matter (Neuckermans et al., 2012)                 | other   |
|                  | 0.0055    | m-1                     | Backscattering by organic matter (Neuckermans et al., 2012)             | other   |
|                  | 0.2522    | ?                       | ?                                                                       | other   |
|                  | -0.0122   | ?                       | ?                                                                       | other   |
|                  | 0.24      | m-1                     | Absorption by pure sea water (Smith and Baker, 1981)                    | other   |
|                  | 0.0007    | m-1                     | Backscattering by pure sea water (Smith and Baker, 1981)                | other   |
|                  | 0.0066    | m-1                     | A*CHL^B (Dmitriev et al., 2007)                                         | other   |
|                  | 0.8       | m-1                     | A*CHL^B (Bricaud et al., 1995, Dmitriev et al., 2007)                   | other   |
|                  | 4.41e-06  | m2 (mg POC_DIA)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 6.12e-06  | m2 (mg POC_FLA)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 7.56e-06  | m2 (mg POC_EMI)-1       | Backscattering by DI (Vaillancourt et al., 2004)                        | other   |
|                  | 0.001     | m2 g-1                  | Absorption by suspended minerals (Neuckermans et al., 2012)             | other   |
|                  | 0.014     | m-1                     | Backscattering by suspended minerals (Neuckermans et al., 2012)         | other   |
|                  | 0.0195    | m2 g-1                  | Absorption by organic matter (Neuckermans et al., 2012)                 | other   |
|                  | 0.005     | m-1                     | Backscattering by OM (Neuckermans et al., 2012)                         | other   |
|                  | 0.2155    | ?                       | ?                                                                       | other   |
|                  | -0.0113   | ?                       | ?                                                                       | other   |
|                  | 1e-05     | mmolSi (m2s)-1          | Erosion constant for slow Si sediment                                   | other   |
|                  | 1e-05     | mmolSi (m2s)-1          | Erosion constant for fast Si sediment                                   | other   |
|                  | 1e-05     | mmolC (m2s)-1           | Erosion constant for slow C sediment                                    | other   |
|                  | 1e-05     | mmolC (m2s)-1           | Erosion constant for fast C sediment                                    | other   |