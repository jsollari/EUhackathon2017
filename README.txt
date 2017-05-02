#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     17.03.2017
#modificado: 02.05.2017
+bin
  |data_analysis_v1.r                      #R script to run the full analysis
  |misc_v2.1.r                             #Miscellaneous functions for Data Mining
+data
  +gz_files
  |ESCO_v3.csv                             #ad hoc mapping between occupations (NACE_Rev2 and ISCO_08) and skills (ISCED-F_13)
  |ESCO_v3.xls                             #ad hoc mapping between occupations (NACE_Rev2 and ISCO_08) and skills (ISCED-F_13)
  |Metadata.xls                            #Description of Metada
  |README.xls                              #Description of Eurostat data sets
  |demo_r_d2jan.tsv                        #Eurostat data set "demo_r_d2jan"
  |earn_ses_hourly.tsv                     #Eurostat data set "earn_ses_hourly"
  |educ_uoe_enra03.tsv                     #Eurostat data set "educ_uoe_enra03"
  |educ_uoe_fine06.tsv                     #Eurostat data set "educ_uoe_fine06"
  |educ_uoe_grad02.tsv                     #Eurostat data set "educ_uoe_grad02"
  |ilc_li41.tsv                            #Eurostat data set "ilc_li41"
  |ilc_lvhl21.tsv                          #Eurostat data set "ilc_lvhl21" 
  |ilc_lvho04n.tsv                         #Eurostat data set "ilc_lvho04n"
  |ilc_mddd21.tsv                          #Eurostat data set "ilc_mddd21" 
  |ilc_peps11.tsv                          #Eurostat data set "ilc_peps11" 
  |jvs_a_nace2.tsv                         #Eurostat data set "jvs_a_nace2"
  |lfso_14leeow.tsv                        #Eurostat data set "lfso_14leeow"
  |lfst_r_lfe2eedu.tsv                     #Eurostat data set "lfst_r_lfe2eedu"
  |lfst_r_lfe2eftpt.tsv                    #Eurostat data set "lfst_r_lfe2eftpt"
  |lfst_r_lfe2emp.tsv                      #Eurostat data set "lfst_r_lfe2emp"
  |lfst_r_lfe2en2.tsv                      #Eurostat data set "lfst_r_lfe2en2"
  |lfst_r_lfp2acedu.tsv                    #Eurostat data set "lfst_r_lfp2acedu"
  |lfst_r_lfu3pers.tsv                     #Eurostat data set "lfst_r_lfu3pers"
  |migr_emi2.tsv                           #Eurostat data set "migr_emi2"
  |nama_10r_2gdp.tsv                       #Eurostat data set "nama_10r_2gdp"
  |nama_10r_2gvagr.tsv                     #Eurostat data set "nama_10r_2gvagr"
  |nama_10r_2hhinc.tsv                     #Eurostat data set "nama_10r_2hhinc"
  |trng_lfse_04.tsv                        #Eurostat data set "trng_lfse_04"
+docs
  |ESS_BigData_Hackathon_data_catalogue.pdf#ESS: Document from "https://ec.europa.eu/eurostat/cros/content/hackathon-data-catalogue_en"
  |ESS_BigData_Hackathon_description.pdf   #ESS: Document from "https://ec.europa.eu/eurostat/cros/content/big-data_en"
  |Methodology_handouts_20170427.pdf       #Methodology presentation
  |Methodology_summary.pdf                 #Methodology description
+results
  +analyse_classification
  +csv
  +data_overview
  +json
  +regression
  +wcna
  |class_nnet_EU_groups_nuts0.xls          #Summary of classification analysis for NUTS0
  |class_nnet_EU_groups_nuts1.xls          #Summary of classification analysis for NUTS1
  |class_nnet_EU_groups_nuts2.xls          #Summary of classification analysis for NUTS2
  |descr_x_nuts0.xls                       #Summary description of x variables for NUTS0
  |descr_x_nuts1.xls                       #Summary description of x variables for NUTS1
  |descr_x_nuts2.xls                       #Summary description of x variables for NUTS2
  |descr_y_nuts0.xls                       #Summary description of y variables for NUTS0
  |descr_y_nuts1.xls                       #Summary description of y variables for NUTS1
  |descr_y_nuts2.xls                       #Summary description of y variables for NUTS2
  |ora_clst_nuts0.xls                      #Summary of Over-representation analysis for NUTS0
  |ora_clst_nuts1.xls                      #Summary of Over-representation analysis for NUTS1
  |ora_clst_nuts2.xls                      #Summary of Over-representation analysis for NUTS2
  |regr_lm_lmktm_ED0-2_nuts0.xls           #Summary of linear regression analysis of "lmktm_ED0-2" for NUTS0
  |regr_lm_lmktm_ED3_4_nuts0.xls           #Summary of linear regression analysis of "lmktm_ED3_4" for NUTS0
  |regr_lm_lmktm_ED5-8_nuts0.xls           #Summary of linear regression analysis of "lmktm_ED5-8" for NUTS0
  |regr_lm_lmktm_OC0_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC0" for NUTS0
  |regr_lm_lmktm_OC1_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC1" for NUTS0
  |regr_lm_lmktm_OC2_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC2" for NUTS0
  |regr_lm_lmktm_OC3_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC3" for NUTS0
  |regr_lm_lmktm_OC4_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC4" for NUTS0
  |regr_lm_lmktm_OC5_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC5" for NUTS0
  |regr_lm_lmktm_OC6_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC6" for NUTS0
  |regr_lm_lmktm_OC7_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC7" for NUTS0
  |regr_lm_lmktm_OC8_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC8" for NUTS0
  |regr_lm_lmktm_OC9_nuts0.xls             #Summary of linear regression analysis of "lmktm_OC9" for NUTS0
  |regr_lm_lmktm_Total_nuts0.xls           #Summary of linear regression analysis of "lmktm_Total" for NUTS0
  |regr_lm_migrt_Total_nuts0.xls           #Summary of linear regression analysis of "migrt_Total" for NUTS0
  |regr_lm_migrt_Y15-24_nuts0.xls          #Summary of linear regression analysis of "migrt_Y15-24" for NUTS0
  |regr_lm_migrt_YGE25-64_nuts0.xls        #Summary of linear regression analysis of "migrt_YGE25-64" for NUTS0
  |regr_lm_mismatch_nuts0.xls              #Summary of linear regression analysis of "mismatch" for NUTS0
  |regr_lm_mismatch_nuts1.xls              #Summary of linear regression analysis of "mismatch" for NUTS1
  |regr_lm_mismatch_nuts2.xls              #Summary of linear regression analysis of "mismatch" for NUTS2
  |regr_nnet_lmktm_ED0-2_nuts0.xls         #Summary of non-linear regression analysis of "lmktm_ED0-2" for NUTS0
  |regr_nnet_lmktm_ED3_4_nuts0.xls         #Summary of non-linear regression analysis of "lmktm_ED3_4" for NUTS0
  |regr_nnet_lmktm_ED5-8_nuts0.xls         #Summary of non-linear regression analysis of "lmktm_ED5-8" for NUTS0
  |regr_nnet_lmktm_OC0_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC0" for NUTS0
  |regr_nnet_lmktm_OC1_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC1" for NUTS0
  |regr_nnet_lmktm_OC2_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC2" for NUTS0
  |regr_nnet_lmktm_OC3_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC3" for NUTS0
  |regr_nnet_lmktm_OC4_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC4" for NUTS0
  |regr_nnet_lmktm_OC5_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC5" for NUTS0
  |regr_nnet_lmktm_OC6_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC6" for NUTS0
  |regr_nnet_lmktm_OC7_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC7" for NUTS0
  |regr_nnet_lmktm_OC8_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC8" for NUTS0
  |regr_nnet_lmktm_OC9_nuts0.xls           #Summary of non-linear regression analysis of "lmktm_OC9" for NUTS0
  |regr_nnet_lmktm_Total_nuts0.xls         #Summary of non-linear regression analysis of "lmktm_Total" for NUTS0
  |regr_nnet_migrt_Total_nuts0.xls         #Summary of non-linear regression analysis of "migrt_Total" for NUTS0
  |regr_nnet_migrt_Y15-24_nuts0.xls        #Summary of non-linear regression analysis of "migrt_Y15-24" for NUTS0
  |regr_nnet_migrt_YGE25-64_nuts0.xls      #Summary of non-linear regression analysis of "migrt_YGE25-64" for NUTS0
  |regr_nnet_mismatch_nuts0.xls            #Summary of non-linear regression analysis of "mismatch" for NUTS0
  |regr_nnet_mismatch_nuts1.xls            #Summary of non-linear regression analysis of "mismatch" for NUTS1
  |regr_nnet_mismatch_nuts2.xls            #Summary of non-linear regression analysis of "mismatch" for NUTS2
  |summ_clst_nuts0.xls                     #Summary description of clustering analysis for NUTS0
  |summ_clst_nuts1.xls                     #Summary description of clustering analysis for NUTS1
  |summ_clst_nuts2.xls                     #Summary description of clustering analysis for NUTS2
  |wcna_nuts0.xls                          #Summary of Weighted Correlation Network Analysis for NUTS0
  |wcna_nuts1.xls                          #Summary of Weighted Correlation Network Analysis for NUTS1
  |wcna_nuts2.xls                          #Summary of Weighted Correlation Network Analysis for NUTS2
  |dat1.RDS                                #Eurostat data sets broken by NUTS1
  |dat2.RDS                                #Eurostat data sets broken by NUTS2
README.txt                                 #This file