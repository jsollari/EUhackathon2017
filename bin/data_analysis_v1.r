#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     13.03.2017
#modificado: 09.06.2017

#setwd("2017.03.17_BigData/bin/")
options("width"=130)
options(contrasts=c("contr.sum","contr.poly"))

## 0. INDEX

# 1. FUNCTIONS
# 2. READ DATA
## 2.1. DATA OVERVIEW
## 2.2. READ DATA
# 3. EDIT DATA
## 3.1. RESHAPE TO WIDE FORMAT
## 3.2. PREPARE VARIABLES
### 3.2.1 Predictors at NUTS2-level
#### 3.2.1.1. Merge datasets
#### 3.2.1.2. Normalize for population size accross regions
### 3.2.2. Predictors at NUTS0-level
### 3.2.3. All predictors
### 3.2.4. Response variable (Skills mismatch)
#### 3.2.4.1. Skills supply
#### 3.2.4.2. Skills demand
#### 3.2.4.3. Skills mismatch
### 3.2.5. Response variable (Mobility)
### 3.2.6. Response variable (Emigration rate)
# 4. ANALYSE DATA
## 4.1. NUTS0-LEVEL (COUNTRY-LEVEL)
### 4.1.1. Summary statistics
### 4.1.2. Scatterplots
### 4.1.3. Barplots
### 4.1.4. Distances between regions
### 4.1.5. Correlation between variables
### 4.1.6. Classification analysis
### 4.1.7. Study classification analysis
#### 4.1.7.1. Descriptive statistics
##### 4.1.7.1.1. Summary statistics on clusters
##### 4.1.7.1.2. Summary statistics on variables
##### 4.1.7.1.3. T-tests on variables
##### 4.1.7.1.4. Boxplots
#### 4.1.7.2. Modeling
##### 4.1.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
### 4.1.8. Enrichment analysis
### 4.1.9. Weighted correlation network analysis (WCNA)
### 4.1.10. Multivariate linear [General LM] and non-linear [ANN] regression
## 4.2. NUTS1-LEVEL
### 4.2.1. Summary statistics
### 4.2.2. Scatterplots
### 4.2.3. Barplots
### 4.2.4. Distances between regions
### 4.2.5. Correlation between variables
### 4.2.6. Classification analysis
### 4.2.7. Study classification analysis
#### 4.2.7.1. Descriptive statistics
##### 4.2.7.1.1. Summary statistics on clusters
##### 4.2.7.1.2. Summary statistics on variables
##### 4.2.7.1.3. T-tests on variables
##### 4.2.7.1.4. Boxplots
#### 4.2.7.2. Modeling
##### 4.2.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
### 4.2.8. Enrichment analysis
### 4.2.9. Weighted correlation network analysis (WCNA)
### 4.2.10. Multivariate linear [General LM] and non-linear [ANN] regression
## 4.3. NUTS1-LEVEL
### 4.3.1. Summary statistics
### 4.3.2. Scatterplots
### 4.3.3. Barplots
### 4.3.4. Distances between regions
### 4.3.5. Correlation between variables
### 4.3.6. Classification analysis
### 4.3.7. Study classification analysis
#### 4.3.7.1. Descriptive statistics
##### 4.3.7.1.1. Summary statistics on clusters
##### 4.3.7.1.2. Summary statistics on variables
##### 4.3.7.1.3. T-tests on variables
##### 4.3.7.1.4. Boxplots
#### 4.3.7.2. Modeling
##### 4.3.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
### 4.3.8. Enrichment analysis
### 4.3.9. Weighted correlation network analysis (WCNA)
### 4.3.10. Multivariate linear [General LM] and non-linear [ANN] regression

## 1. FUNCTIONS
{
# Key:
# a - auxiliar
# d - dataframe
# f - file name
# i - iterator
# l - list
# n - vector size
# m - matrix
# s - selected samples
# t - table
# v - vector
# w - weights
source("misc_v2.1.r")

# Check variables present in datasets
# files - files of datasets to check
# geo_cat - geographical tags to use for particular datasets [e.g. single-year datasets]
# f1 - path for files
overview_data1 = function(files,geo_cat=NULL,f1){
    nfiles = length(files)
    cnam1 = vector("list",length=nfiles)
    names(cnam1) = files
    for(i1 in files){
        f2 = paste(f1,i1,sep="")
        cnam1[[i1]] = read.table(f2,header=TRUE,sep="\t",dec=".",na.strings=": ",
          nrows=2,stringsAsFactors=FALSE,encoding="UTF-8")
        nams1 = strsplit(colnames(cnam1[[i1]])[1],".",fixed=TRUE)[[1]] #get colnames from 1st row
        if("time" %in% nams1){
            nams1 = nams1[-match("time",nams1)] #remove tag "time" 
        }
        n1 = length(nams1)
        col.list = strsplit(cnam1[[i1]][,1],",") #turn the 1st row into a list, split by ,
        cnam1[[i1]] = cbind(data.frame(t(sapply(col.list,function(x)x)),
          stringsAsFactors=FALSE),cnam1[[i1]][,-1])
        n2 = ncol(cnam1[[i1]])
        colnames(cnam1[[i1]]) = c(nams1,colnames(cnam1[[i1]][(n1+1):n2]))
        if(i1 %in% c("lfso_14leeow.tsv")){
            if(is.null(geo_cat)){
                geo_cat = c("AT","BE","BG","CY","CZ","DE","DK","EE","EL","ES",
                  "FI","FR","HR","HU","IE","IT","LT","LU","LV","MT","NL","PL",
                  "PT","RO","SE","SI","SK","UK")
            }
            s1 = colnames(cnam1[[i1]]) %in% c("unit","worktime","isco08","isced11",
              "mgstatus",geo_cat)
            cnam1[[i1]] = cnam1[[i1]][,s1]
            a1 = geo_cat[geo_cat %in% colnames(cnam1[[i1]])]
            d1 = reshape(cnam1[[i1]],varying=a1,v.names="X2014",direction="long",
              times=a1,timevar="geo")
            d1 = d1[,-ncol(d1)]
            cnam1[[i1]] = d1
        }      
    }
    return(cnam1)
}

# Check which dataset contains certain variable
# nams1 - names of variables to check
# cnam1 - list of the variables of data frames
overview_data2 = function(nams1,cnam1){
    for(i1 in nams1){
        tmp = c()
        for(i2 in names(cnam1)){
            if(i1 %in% colnames(cnam1[[i2]])){
                tmp = c(tmp,i2)
            }
        }
        print(i1)
        print(tmp)
        cat("\n")
    }
}

# Read data from files to a lsit of data frames [Use demographic data from 2014]
# files - files to read data from
# geo_cat - geographical tags to keep
# sex_cat - sex categories to keep {"F","M","T"}
# rem_col - columns to remove from data sets
read_data_v2 = function(files,geo_cat,sex_cat,rem_col){
    age_cat0 = c("TOTAL")
    age_cat1 = c("Y_LT1","Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9","Y10",
      "Y11","Y12","Y13","Y14")
    age_cat2 = c("Y15","Y16","Y17","Y18","Y19","Y20","Y21","Y22","Y23","Y24")
    age_cat3 = c("Y25","Y26","Y27","Y28","Y29","Y30","Y31","Y32","Y33","Y34",
      "Y35","Y36","Y37","Y38","Y39","Y40","Y41","Y42","Y43","Y44","Y45","Y46",
      "Y47","Y48","Y49","Y50","Y51","Y52","Y53","Y54","Y55","Y56","Y57","Y58",
      "Y59","Y60","Y61","Y62","Y63","Y64")
    age_cat4 = c("Y65","Y66","Y67","Y68","Y69","Y70","Y71","Y72","Y73","Y74")
    age_cat5 = c("Y75","Y76","Y77","Y78","Y79","Y80","Y81","Y82","Y83","Y84",
      "Y85","Y86","Y87","Y88","Y89","Y90","Y91","Y92","Y93","Y94","Y95","Y96",
      "Y97","Y98","Y99","Y_OPEN","Y_GE100")
    age_cat6 = c("Y15-24","Y25-64")
    age_cat7 = c("Y15-24","Y25-64","Y_GE65")
    age_cat8 = c("Y15-24","Y_GE25")
    age_cat9 = c("Y25-64")
    agedef_cat = c("REACH","COMPLET")
    geo_cat1 = geo_cat[sapply(geo_cat,nchar)==2]
    isco_cat1 = c("OC0","OC1","OC2","OC3","OC4","OC5","OC6","OC7","OC8","OC9")
    isco_cat2 = c("OC0","OC1","OC2","OC3","OC4","OC5","OC6","OC7","OC8","OC9",
      "TOTAL")
    isced11_cat1 = c("ED0-2","ED3_4","ED5-8")
    isced11_cat2 = c("ED5-8")
    isced11_cat3 = c("ED0-2","ED3_4","ED5-8","TOTAL")
    iscedf13_cat1 = c("TOTAL","F00","F01","F02","F03","F04","F05","F06","F07",
      "F08","F09","F10")
    indic_em_cat = c("JOBVAC","JOBOCC")
    indic_se_cat = c("ERN_MN_PPS")
    mgstatus_cat = c("IMG","TOTAL")
    nace_cat1 = c("A-S","B-S","A","B","C","D","E","F","G","H","I","J","K","L",
      "M","N","O","P","Q","R","S")
    nace_cat2 = c("A","B-E","F","G-I","J","K","L","M_N","O-Q","R-U")
    nace_cat3 = c("B-F","G-N","P-S")
    na_item_cat = c("B6N")
    sizeclas_cat = c("TOTAL")
    time_cat1 = c(paste("X",1990:2012,sep=""),paste("X",2014:2017,sep=""))
    time_cat2 = c(paste("X",1990:2013,sep=""),paste("X",2015:2017,sep=""))
    unit_cat1 = c("NR")                     #Number
    unit_cat2 = c("PC_POP")                 #Percentage of total population
    unit_cat3 = c("PC_Y_LT60")              #Percentage of total population aged less than 60
    unit_cat4 = c("AVG")                    #Average
    unit_cat5 = c("PPS_HAB")                #Purchasing power standard per inhabitant
    unit_cat6 = c("PCH_PRE")                #Percentage change on previous period 
    unit_cat7 = c("PPCS_HAB")               #Purchasing power standard on finals consumption per inhabitant
    unit_cat8 = c("THS")                    #Thousand
    unit_cat9 = c("PC")                     #Percentage
    unit_cat10 = c("PC_GDP")                #Percentage per GDP
    unit_cat11 = c("THS_PER")               #Thousand people
    worktime_cat1 = c("PT","FT")
    worktime_cat2 = c("TOTAL")
    nfiles = length(files)
    dat1 = vector("list",length=nfiles)
    names(dat1) = files
    for(i1 in files){
        f1 = paste("../data/",i1,sep="")
        print(paste("Loading file:",f1))
        dat1[[i1]] = read.table(f1,header=TRUE,sep="\t",dec=".",na=c(": ",": b",
          ": c",": d",": u",": z",": bc",": bd",": bu"),stringsAsFactors=FALSE,
          encoding="UTF-8")
        nams1 = strsplit(colnames(dat1[[i1]])[1],".",fixed=TRUE)[[1]] #get colnames from 1st row
        if("time" %in% nams1){
            nams1 = nams1[-match("time",nams1)]                       #remove tag "time" 
        }
        n1 = length(nams1)
        col.list = strsplit(dat1[[i1]][,1],",")                       #turn the 1st row into a list, split by ","
        dat1[[i1]] = cbind(data.frame(t(sapply(col.list,function(x)x)),
          stringsAsFactors=FALSE),dat1[[i1]][,-1])
        n2 = ncol(dat1[[i1]])
        colnames(dat1[[i1]]) = c(nams1,colnames(dat1[[i1]][(n1+1):n2]))
        if(i1 %in% c("lfso_14leeow.tsv")){
            s1 = colnames(dat1[[i1]]) %in% c("unit","worktime","isco08","isced11",
              "mgstatus",geo_cat1)
            dat1[[i1]] = dat1[[i1]][,s1]
            a1 = geo_cat1[geo_cat1 %in% colnames(dat1[[i1]])]
            d1 = reshape(dat1[[i1]],varying=a1,v.names="X2014",direction="long",
              times=a1,timevar="geo")
            d1 = d1[,-ncol(d1)]
            dat1[[i1]] = d1
        }
        #Age {...}
        if(!"age" %in% colnames(dat1[[i1]])){
            a1 = rep("TOTAL",nrow(dat1[[i1]]))
            dat1[[i1]] = cbind(a1,dat1[[i1]])
            colnames(dat1[[i1]]) = c("age",colnames(dat1[[i1]])[-1])
            s1 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("earn_ses_hourly.tsv")){
            s1 = dat1[[i1]][,"age"] %in% c(age_cat0)
        }else if(i1 %in% c("demo_r_d2jan.tsv","migr_emi2.tsv")){
            s1 = dat1[[i1]][,"age"] %in% c(age_cat0,age_cat1,age_cat2,age_cat3,
              age_cat4,age_cat5)
        }else if(i1 %in% c("lfst_r_lfe2eedu.tsv","lfst_r_lfe2en2.tsv")){
            s1 = dat1[[i1]][,"age"] %in% age_cat6
        }else if(i1 %in% c("lfst_r_lfe2emp.tsv")){
            s1 = dat1[[i1]][,"age"] %in% age_cat7
        }else if(i1 %in% c("lfst_r_lfu3pers.tsv")){
            s1 = dat1[[i1]][,"age"] %in% age_cat8
        }else if(i1 %in% c("trng_lfse_04.tsv")){
            s1 = dat1[[i1]][,"age"] %in% age_cat9
        }
        #Geo {...}
        if(!"geo" %in% colnames(dat1[[i1]])){
            a1 = rep("EU28",nrow(dat1[[i1]]))
            dat1[[i1]] = cbind(a1,dat1[[i1]])
            colnames(dat1[[i1]]) = c("geo",colnames(dat1[[i1]])[-1])
            s2 = rep(TRUE,nrow(dat1[[i1]]))
        }else{
            s2 = dat1[[i1]][,"geo"] %in% c(geo_cat)
        }
        #Sex {F; M; T}
        if(!"sex" %in% colnames(dat1[[i1]])){
            a1 = rep("T",nrow(dat1[[i1]]))
            dat1[[i1]] = cbind(a1,dat1[[i1]])
            colnames(dat1[[i1]]) = c("sex",colnames(dat1[[i1]])[-1])
            s3 = rep(TRUE,nrow(dat1[[i1]]))
        }else{
            s3 = dat1[[i1]][,"sex"] %in% sex_cat
        }
        #Unit {...}
        if(!"unit" %in% colnames(dat1[[i1]])){
            a1 = rep("INDIC",nrow(dat1[[i1]]))
            dat1[[i1]] = cbind(a1,dat1[[i1]])
            colnames(dat1[[i1]]) = c("unit",colnames(dat1[[i1]])[-1])
            s4 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("demo_r_d2jan.tsv","educ_uoe_grad02.tsv",
          "migr_emi2.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat1
        }else if(i1 %in% c("ilc_li41.tsv","ilc_mddd21.tsv","ilc_peps11.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat2
        }else if(i1 %in% c("ilc_lvhl21.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat3
        }else if(i1 %in% c("ilc_lvho04n.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat4
        }else if(i1 %in% c("nama_10r_2gdp.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat5
        }else if(i1 %in% c("nama_10r_2gvagr.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat6
        }else if(i1 %in% c("nama_10r_2hhinc.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat7
        }else if(i1 %in% c("lfst_r_lfe2eedu.tsv","lfst_r_lfe2eftpt.tsv",
          "lfst_r_lfe2emp.tsv","lfst_r_lfe2en2.tsv","lfst_r_lfp2acedu.tsv",
          "lfst_r_lfu3pers.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat8
        }else if(i1 %in% c("trng_lfse_04.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat9
        }else if(i1 %in% c("educ_uoe_fine06.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat10
        }else if(i1 %in% c("lfso_14leeow.tsv")){
            s4 = dat1[[i1]][,"unit"] %in% unit_cat11
        }
        #Age definition {COMPLET=Completed years; REACH=Reached during the year}
        if("agedef" %in% colnames(dat1[[i1]])){
            s5 = dat1[[i1]][,"agedef"] %in% agedef_cat
        }else{
            s5 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #Employment indicator {...}
        if("indic_em" %in% colnames(dat1[[i1]])){
            s6 = dat1[[i1]][,"indic_em"] %in% indic_em_cat
        }else{
            s6 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #Structure of earnings indicator {...}
        if("indic_se" %in% colnames(dat1[[i1]])){
            s7 = dat1[[i1]][,"indic_se"] %in% indic_se_cat
        }else{
            s7 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #ISCED 11 {...}
        if(!"isced11" %in% colnames(dat1[[i1]])){
            s8 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("lfst_r_lfe2eedu.tsv","lfst_r_lfp2acedu.tsv")){
            s8 = dat1[[i1]][,"isced11"] %in% isced11_cat1
        }else if(i1 %in% c("educ_uoe_grad02.tsv","educ_uoe_fine06.tsv")){
            s8 = dat1[[i1]][,"isced11"] %in% isced11_cat2
        }else if(i1 %in% c("lfso_14leeow.tsv")){
            s8 = dat1[[i1]][,"isced11"] %in% isced11_cat3
        }
        #ISCED-F 13 {...}
        if("iscedf13" %in% colnames(dat1[[i1]])){
            s9 = dat1[[i1]][,"iscedf13"] %in% iscedf13_cat1
        }else{
            s9 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #ISCO 08 {...}
        if(!"isco08" %in% colnames(dat1[[i1]])){
            s10 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("earn_ses_hourly.tsv")){
            s10 = dat1[[i1]][,"isco08"] %in% isco_cat1
        }else if(i1 %in% c("lfso_14leeow.tsv","jvs_a_nace2.tsv")){
            s10 = dat1[[i1]][,"isco08"] %in% isco_cat2
        }
        #Migration status
        if("mgstatus" %in% colnames(dat1[[i1]])){
            s11 = dat1[[i1]][,"mgstatus"] %in% mgstatus_cat
        }else{
            s11 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #National account indicator {B5N=Balance of priamry icomes; B6N=Disposable income}
        if("na_item" %in% colnames(dat1[[i1]])){
            s12 = dat1[[i1]][,"na_item"] %in% na_item_cat
        }else{
            s12 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #NACE Rev.2 {...}
        if(!"nace_r2" %in% colnames(dat1[[i1]])){
            s13 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("jvs_a_nace2.tsv")){
            s13 = dat1[[i1]][,"nace_r2"] %in% nace_cat1
        }else if(i1 %in% c("lfst_r_lfe2en2.tsv")){
            s13 = dat1[[i1]][,"nace_r2"] %in% nace_cat2
        }else if(i1 %in% c("earn_ses_hourly.tsv")){
            s13 = dat1[[i1]][,"nace_r2"] %in% nace_cat3
        }
        #Size classes in number of employees {TOTAL=Total; GE10=10 employees or more}
        if("sizeclas" %in% colnames(dat1[[i1]])){
            s14 = dat1[[i1]][,"sizeclas"] %in% sizeclas_cat
        }else{
            s14 = rep(TRUE,nrow(dat1[[i1]]))
        }
        #Working time {FT=Full-time; PT=Part-time}
        if(!"worktime" %in% colnames(dat1[[i1]])){
            s15 = rep(TRUE,nrow(dat1[[i1]]))
        }else if(i1 %in% c("lfst_r_lfe2eftpt.tsv")){
            s15 = dat1[[i1]][,"worktime"] %in% worktime_cat1
        }else if(i1 %in% c("lfso_14leeow.tsv","earn_ses_hourly.tsv")){
            s15 = dat1[[i1]][,"worktime"] %in% worktime_cat2
        }
        #Time {X1990,...,X2017}
        if(i1 %in% c("nama_10r_2hhinc.tsv","educ_uoe_fine06.tsv")){
            s16 = colnames(dat1[[i1]]) %in% time_cat1
            dat1[[i1]] = dat1[[i1]][s1 & s2 & s3 & s4 & s5 & s6 & s7 & s8 & s9 &
              s10 & s11 & s12 & s13 & s14 & s15,!s16]
            dat1[[i1]][,"X2013"] = as.numeric(gsub("[b,c,d,e,f,i,n,p,r,s,u,z]",
              "",dat1[[i1]][,"X2013"]))
        }else{
            s16 = colnames(dat1[[i1]]) %in% time_cat2
            dat1[[i1]] = dat1[[i1]][s1 & s2 & s3 & s4 & s5 & s6 & s7 & s8 & s9 &
              s10 & s11 & s12 & s13 & s14 & s15,!s16]
            dat1[[i1]][,"X2014"] = as.numeric(gsub("[b,c,d,e,f,i,n,p,r,s,u,z]",
              "",dat1[[i1]][,"X2014"]))
        }
        #Remove unneeded columns
        cnames1 = colnames(dat1[[i1]])
        if(i1 %in% c("lfst_r_lfe2eftpt.tsv")){
            s17 = cnames1 %in% rem_col & cnames1 != "worktime"
            dat1[[i1]] = dat1[[i1]][,!s17]
        }else if(i1 %in% c("lfst_r_lfe2eedu.tsv","lfst_r_lfp2acedu.tsv",
          "lfso_14leeow.tsv")){
            s17 = cnames1 %in% rem_col & cnames1 != "isced11"
            dat1[[i1]] = dat1[[i1]][,!s17]
        }else{
            s17 = cnames1 %in% rem_col
            dat1[[i1]] = dat1[[i1]][,!s17]
        }
        #Edit Employment indicator
        if(i1 == c("jvs_a_nace2.tsv")){
            s18_cat0 = dat1[[i1]][,"indic_em"] %in% indic_em_cat
            geo1 = names(table(dat1[[i1]][,"geo"]))
            isco1 = names(table(dat1[[i1]][,"isco08"]))
            nace1 = names(table(dat1[[i1]][,"nace_r2"]))
            d1 = data.frame(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(),V6=c(),V7=c(),
              V8=c())
            for(i2 in geo1){
                s19 = dat1[[i1]][,"geo"] == i2
                for(i3 in isco1){
                    s20 = dat1[[i1]][,"isco08"] == i3
                    for(i4 in nace1){
                        s21 = dat1[[i1]][,"nace_r2"] == i4
                        a1_cat0 = ifelse(sum(s18_cat0 & s19 & s20 & s21) > 0,
                          sum(dat1[[i1]][s18_cat0 & s19 & s20 & s21,"X2014"],
                          na.rm=FALSE),NA)
                        a1 = data.frame("NR","T","TOTAL",i4,i3,"JOBTOT",i2,
                          a1_cat0,stringsAsFactors=FALSE)
                        d1 = rbind(d1,a1)
                    }
                }
            }
            colnames(d1) = colnames(dat1[[i1]])
            d1[,ncol(d1)] = as.numeric(d1[,ncol(d1)])
            dat1[[i1]] = d1
        }
        #Edit "age"
        if(i1 %in% c("migr_emi2.tsv")){
            s18_cat0 = dat1[[i1]][,"age"] %in% age_cat0
            s18_cat1 = dat1[[i1]][,"age"] %in% age_cat1
            s18_cat2 = dat1[[i1]][,"age"] %in% age_cat2
            s18_cat3 = dat1[[i1]][,"age"] %in% age_cat3
            s18_cat4 = dat1[[i1]][,"age"] %in% age_cat4
            s18_cat5 = dat1[[i1]][,"age"] %in% age_cat5
            geo1 = names(table(dat1[[i1]][,"geo"]))
            sex1 = names(table(dat1[[i1]][,"sex"]))
            agedef1 = names(table(dat1[[i1]][,"agedef"]))
            d1 = data.frame(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(),V6=c())
            for(i2 in geo1){
                s19 = dat1[[i1]][,"geo"] == i2
                for(i3 in sex1){
                    s20 = dat1[[i1]][,"sex"] == i3
                    for(i4 in agedef1){
                        s21 = dat1[[i1]][,"agedef"] == i4
                        a1_cat0 = sum(dat1[[i1]][s18_cat0 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1_cat1 = sum(dat1[[i1]][s18_cat1 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1_cat2 = sum(dat1[[i1]][s18_cat2 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1_cat3 = sum(dat1[[i1]][s18_cat3 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1_cat4 = sum(dat1[[i1]][s18_cat4 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1_cat5 = sum(dat1[[i1]][s18_cat5 & s19 & s20 & s21,
                          "X2014"],na.rm=FALSE)
                        a1 = cbind(c("TOTAL","Y0-14","Y15-24","Y25-64","Y65-74",
                          "YGE75"),rep(i4,6),rep("NR",6),rep(i3,6),rep(i2,6),
                          c(a1_cat0,a1_cat1,a1_cat2,a1_cat3,a1_cat4,a1_cat5))
                        d1 = rbind(d1,a1)
                    }
                }
            }
            colnames(d1) = colnames(dat1[[i1]])
            d1[,ncol(d1)] = as.numeric(as.character(d1[,ncol(d1)]))
            dat1[[i1]] = d1
        }else if(i1 %in% c("demo_r_d2jan.tsv")){
            s18_cat0 = dat1[[i1]][,"age"] %in% age_cat0
            s18_cat1 = dat1[[i1]][,"age"] %in% age_cat1
            s18_cat2 = dat1[[i1]][,"age"] %in% age_cat2
            s18_cat3 = dat1[[i1]][,"age"] %in% age_cat3
            s18_cat4 = dat1[[i1]][,"age"] %in% age_cat4
            s18_cat5 = dat1[[i1]][,"age"] %in% age_cat5
            geo1 = names(table(dat1[[i1]][,"geo"]))
            sex1 = names(table(dat1[[i1]][,"sex"]))
            d1 = data.frame(V1=c(),V2=c(),V3=c(),V4=c(),V5=c())
            for(i2 in geo1){
                s19 = dat1[[i1]][,"geo"] == i2
                for(i3 in sex1){
                    s20 = dat1[[i1]][,"sex"] == i3
                    a1_cat0 = sum(dat1[[i1]][s18_cat0 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1_cat1 = sum(dat1[[i1]][s18_cat1 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1_cat2 = sum(dat1[[i1]][s18_cat2 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1_cat3 = sum(dat1[[i1]][s18_cat3 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1_cat4 = sum(dat1[[i1]][s18_cat4 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1_cat5 = sum(dat1[[i1]][s18_cat5 & s19 & s20,"X2014"],
                      na.rm=FALSE)
                    a1 = cbind(rep("NR",6),rep(i3,6),c("TOTAL","Y0-14",
                      "Y15-24","Y25-64","Y65-74","YGE75"),rep(i2,6),
                      c(a1_cat0,a1_cat1,a1_cat2,a1_cat3,a1_cat4,a1_cat5))
                    d1 = rbind(d1,a1)
                }
            }
            colnames(d1) = colnames(dat1[[i1]])
            d1[,ncol(d1)] = as.numeric(as.character(d1[,ncol(d1)]))
            dat1[[i1]] = d1
        }
        rownames(dat1[[i1]]) = 1:nrow(dat1[[i1]])
    }
    return(dat1)
}

# Reshape long format to wide format [Use demographic data from 2014]
# dat1 - data frame of variables to be analysed
# sexbreak - break data by "sex"
reshape2wide_v2 = function(dat1,sexbreak=FALSE){
    dat1w = dat1
    files = names(dat1w)
    for(i1 in files){
        #Age
        nams1 = colnames(dat1w[[i1]])
        if("age" %in% nams1){
            if(i1 %in% c("demo_r_d2jan.tsv","lfst_r_lfu3pers.tsv",
              "lfst_r_lfe2eedu.tsv","lfst_r_lfe2emp.tsv","lfst_r_lfe2en2.tsv",
              "migr_emi2.tsv")){
                vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
                idvar1 = nams1[-match(c("age",vnams1),nams1)]
                dat1w[[i1]] = reshape(dat1w[[i1]],timevar="age",v.names=vnams1,
                  idvar=idvar1,direction="wide",sep="_age_")
            }
        }
        #Sex
        if(sexbreak){
            nams1 = colnames(dat1w[[i1]])
            if("sex" %in% nams1){
                if(i1 %in% c("nama_10r_2hhinc.tsv","educ_uoe_fine06.tsv")){       
                    vnams1 = grep("X2013",colnames(dat1w[[i1]]),value=TRUE)
                    idvar1 = nams1[-match(c("sex",vnams1),nams1)]
                    dat1w[[i1]] = reshape(dat1w[[i1]],timevar="sex",v.names=vnams1,
                      idvar=idvar1,direction="wide",sep="_sex_")
                }else{
                    vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
                    idvar1 = nams1[-match(c("sex",vnams1),nams1)]
                    dat1w[[i1]] = reshape(dat1w[[i1]],timevar="sex",v.names=vnams1,
                      idvar=idvar1,direction="wide",sep="_sex_")
                }
            }
        }
        #Age definition
        nams1 = colnames(dat1w[[i1]])
        if("agedef" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("agedef",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="agedef",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_agedef_")
        }
        #Employment indicator
        nams1 = colnames(dat1w[[i1]])
        if("indic_em" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("indic_em",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="indic_em",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_indicem_")
        }
        #ISCED 11
        nams1 = colnames(dat1w[[i1]])
        if("isced11" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("isced11",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="isced11",v.names=vnams1,
            idvar=idvar1,direction="wide",sep="_isced11_")
        }
        #ISCED-F 13
        nams1 = colnames(dat1w[[i1]])
        if("iscedf13" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("iscedf13",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="iscedf13",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_iscedf13_")
        }
        #ISCO 08
        nams1 = colnames(dat1w[[i1]])
        if("isco08" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("isco08",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="isco08",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_isco08_")
        }
        #Migration Status
        nams1 = colnames(dat1w[[i1]])
        if("mgstatus" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("mgstatus",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="mgstatus",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_mgstatus_")
        }
        #NACE Rev.2
        nams1 = colnames(dat1w[[i1]])
        if("nace_r2" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("nace_r2",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="nace_r2",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_nace_")
        }
        #Working time
        nams1 = colnames(dat1w[[i1]])
        if("worktime" %in% nams1){
            vnams1 = grep("X2014",colnames(dat1w[[i1]]),value=TRUE)
            idvar1 = nams1[-match(c("worktime",vnams1),nams1)]
            dat1w[[i1]] = reshape(dat1w[[i1]],timevar="worktime",v.names=vnams1,
              idvar=idvar1,direction="wide",sep="_worktime_")
        }   
    }
    return(dat1w) 
}

# Map info from "jvs_a_nace2" (i.e. NACE Rev2 and ISCO 08) to ISCED-F 13
# mat1 - Job skills broken by {A,...,S} per {OC0,...,OC9}
# mat2 - Job skills broken by {"A-S"} per {"OC0",...,"OC9"}
# mat3 - Job skills broken by {"B-S"} per {"OC0",...,"OC9"}
# mat4 - Job skills {"Total"}
# map1 - mapping between nace_r2 and isco08 ("skills1") and iscedf13 ("skills2")
# geo_cat - geographical tags [if geo_cat==NULL then NUTS0]
# rmF00 - remove "F00" (i.e. education fields without information)
# prop - get proportion of total
# usetot - Use total from data [Use of sum(na.rm=TRUE) may give incorrect proportions] 
mapSCQO = function(mat1,mat2,mat3,mat4,map1,geo_cat=NULL,rmF00=TRUE,prop=TRUE,
  usetot=TRUE){
    if(is.null(geo_cat)){
        geo_cat = c("AT","BE","BG","CY","CZ","DE","DK","EE","EL","ES","FI","FR",
          "HR","HU","IE","IT","LT","LU","LV","MT","NL","PL","PT","RO","SE","SI",
          "SK","UK")
    }
    y1 = matrix(NA,nrow=length(geo_cat),ncol=11)
    rownames(y1) = geo_cat
    colnames(y1) = c("F00","F01","F02","F03","F04","F05","F06","F07","F08",
      "F09","F10")
    for(i1 in geo_cat){
        if(i1 %in% rownames(mat1)){
            v1 = data.frame(cbind(names(mat1),t(mat1[i1,])))
            v2 = data.frame(cbind(names(mat2),t(mat2[i1,])))
            v3 = data.frame(cbind(names(mat3),t(mat3[i1,])))
            colnames(v1) = c("skills1","value")
            colnames(v2) = c("skills1","value")
            colnames(v3) = c("skills1","value")
            rownames(v1) = 1:nrow(v1)
            rownames(v2) = 1:nrow(v2)
            rownames(v3) = 1:nrow(v3)
            v4 = tmp4[i1,1]
            if(any(!is.na(v1[,"value"]))){
                v1[,"skills2"] = map1[match(v1[,"skills1"],map1[,"skills1"]),
                  "skills2"]
                v5 = tapply(as.numeric(v1[,"value"]),v1[,"skills2"],FUN=sum)
                if(!usetot | is.na(v4)){ #Be careful with na.rm=TRUE (incorrect proportions)
                    v4 = sum(as.numeric(v1[,"value"]),na.rm=TRUE)
                }
            }else if(any(!is.na(v2[,"value"]))){
                v2[,"skills2"] = map1[match(v2[,"skills1"],map1[,"skills1"]),
                  "skills2"]
                v5 = tapply(as.numeric(v2[,"value"]),v2[,"skills2"],FUN=sum)
                if(!usetot | is.na(v4)){ #Be careful with na.rm=TRUE (incorrect proportions)
                    v4 = sum(as.numeric(v2[,"value"]),na.rm=TRUE)
                }
            }else if(any(!is.na(v3[,"value"]))){
                v3[,"skills2"] = map1[match(v2[,"skills1"],map1[,"skills1"]),
                  "skills2"]
                v5 = tapply(as.numeric(v3[,"value"]),v3[,"skills2"],FUN=sum)
                if(!usetot | is.na(v4)){ #Be careful with na.rm=TRUE (incorrect proportions)
                    v4 = sum(as.numeric(v3[,"value"]),na.rm=TRUE)
                }
            }else{
                v4 = v5 = NA
            }
            if(prop){
                if(rmF00 & !is.na(v5["F00"])){
                    v5 = 100*v5/(v4 - v5["F00"])
                    v5["F00"] = NA
                }else{
                    v5 = 100*v5/v4
                }
            }else{
                if(rmF00 & !is.na(v5["F00"])){
                    v5["F00"] = NA
                }
            }
            s1 = match(names(v5),colnames(y1))
            y1[i1,s1] = v5
        }
    }
    return(y1)
}

# Plot value~time*area (wide format)
# dat1 - data frame of variables to be analysed
# col1 - colors associated to subj codes
# units1 - 
# legend1 - 
# k1 - number of plots per row to print
# k2 - number of plots per columns to print
# rm_xlab - remove x label information
# rm_ylab - remove y label information
# rm_main - remove main label information
# horiz - plot as landscape
# f1 - output filename
plot_data = function(dat1,col1=NULL,units1,legend1=NULL,k1=4,k2=2,
  rm_xlab=TRUE,rm_ylab=TRUE,rm_main=TRUE,horiz=FALSE,f1){
    nplot = k1*k2
    if(horiz){
        iwidth = 10
        iheight = 7
    }else{
        iwidth = 7
        iheight = 10
    }
    for(i1 in c("eps","tiff")){
        f2 = paste(f1,"_1",sep="")
        if(i1 == "eps"){
            postscript(file=paste(f2,".eps",sep=""),width=iwidth,height=iheight,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f2,".tif",sep=""),units="in",width=iwidth,
              height=iheight,res=300,compression="lzw",family="sans")
        }
        if(horiz){
            oldpar = par(mfrow=c(k2,k1),mar=c(1.7,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        }else{
            oldpar = par(mfrow=c(k1,k2),mar=c(1.7,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
              tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        }
        iplot = 0
        for(i2 in 1:ncol(dat1)){
            iplot = iplot + 1
            if(rm_xlab){
                xlab1 = ""
            }else{
                xlab1 = colnames(dat1)[i2]
            }
            if(rm_ylab){
                ylab1 = ""
                yaxt1 = "n"
            }else{
                ylab1 = units1[i2]
                yaxt1 = "s"
            }
            if(rm_main){
                main1 = ""
            }else{
                main1 = colnames(dat1)[i2]
            }
            plot(x=jitter(rep(1,nrow(dat1)),factor=1),dat1[,i2],
              bg=col1,col=1,pch=22,xlim=c(0.95,1.05),xaxt="n",xlab=xlab1,
              yaxt=yaxt1,ylab=ylab1,main=main1,cex=1.1)
            if(is.null(legend1)){
                legend("topright",rownames(dat1),col=col1,fill=col1,ncol=2,
                  bty="n",cex=0.7,pt.cex=0.8)
            }else{
                legend("topright",names(legend1),col=legend1,fill=legend1,
                  ncol=2,bty="n",cex=0.7,pt.cex=0.8)
            }
            if(iplot%%nplot==0 & iplot < ncol(dat1)){
                invisible(dev.off())
                f2 = paste(f1,"_",1+iplot/nplot,sep="")
                if(i1 == "eps"){
                    postscript(file=paste(f2,".eps",sep=""),width=iwidth,
                      height=iheight,colormodel="rgb",horizontal=FALSE,
                      onefile=FALSE,paper="special",family="Helvetica")
                }else{
                    tiff(file=paste(f2,".tif",sep=""),units="in",width=iwidth,
                      height=iheight,res=300,compression = "lzw",family="sans")
                }
                if(horiz){
                    oldpar = par(mfrow=c(k2,k1),mar=c(1.7,2.3,1.0,0.2),
                      mgp=c(1.0,0.2,0),tcl=-0.2,cex=0.8,cex.lab=0.8,
                      cex.axis=0.8,cex.main=0.9)
                }else{
                    oldpar = par(mfrow=c(k1,k2),mar=c(1.7,2.3,1.0,0.2),
                      mgp=c(1.0,0.2,0),tcl=-0.2,cex=0.8,cex.lab=0.8,
                      cex.axis=0.8,cex.main=0.9)
                }
            }
        }
        invisible(dev.off())
    }
}

# Compare two data sets using barplots
# dat1 - data set 1
# dat2 - data set 2
# col1 - colors associated to subj codes
# border - add border to columns
# ylab1 - y label for plot 1
# ylab2 - y label for plot 2
# main - title for plots
# cex_names - size of x labels
# f1 - output filename
plot_barplot = function(dat1,dat2,col1=NULL,border=1,ylab1,ylab2,main=NULL,
  cex_names=0.9,f1){
    n1 = nrow(dat1)
    if(is.null(col1)){
        col1 = rainbow(n1)
    }
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=10.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=10.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfrow=c(2,1),mar=c(1.7,2.3,1.0,0.2),mgp=c(1.0,0.2,0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        barplot(dat1,beside=TRUE,cex.names=cex_names,las=2,col=col1,
          border=border,ylim=c(0,50),ylab=ylab1)
        legend("topright",rownames(dat1),fill=col1[-1],ncol=6,bty="n",cex=0.7,
          pt.cex=0.8)
        barplot(dat2,beside=TRUE,cex.names=cex_names,las=2,col=col1,
          border=border,ylim=c(0,50),ylab=ylab2)
        legend("topright",rownames(dat2),fill=col1[-1],ncol=6,bty="n",cex=0.7,
          pt.cex=0.8)
        par(oldpar)
        invisible(dev.off())
    }
}

}
## 2. READ DATA
{
#FILES BROKEN BY NUTS 1:
#  earn_ses_hourly     - Structure of earnings survey: hourly earnings
#  educ_uoe_fine06     - Total public expenditure on education by education level and programme orientation - as % of GDP
#  educ_uoe_grad02     - Graduates by education level, programme orientation, sex and field of education
#  lfso_14leeow        - Employees by migration status, educational attainment level, occupation and working time
#  migr_emi2           - Emigration by age and sex
#FILES BROKEN BY NUTS 2:
#  demo_r_d2jan        - Population on 1 January by age, sex and NUTS 2 region
#  ilc_li41            - At-risk-of-poverty rate by NUTS 2 regions 
#  ilc_lvhl21          - People living in households with very low work intensity by NUTS 2 regions (population aged 0 to 59 years) 
#  ilc_lvho04n         - Average number of rooms per person by NUTS 2 region 
#  ilc_mddd21          - Severe material deprivation rate by NUTS 2 regions 
#  ilc_peps11          - People at risk of poverty or social exclusion by NUTS 2 regions 
#  jvs_a_nace2         - Job vacancy statistics by occupation, NUTS 2 regions and NACE Rev. 2 activity - annual data (2008-2015) 
#  lfst_r_lfe2eedu     - Employment by sex, age, educational attainment level and NUTS 2 regions (1000)
#  lfst_r_lfe2eftpt    - Employment by full-time/part-time, sex and NUTS 2 regions (1 000)
#  lfst_r_lfe2emp      - Employment by sex, age and NUTS 2 regions (1 000)
#  lfst_r_lfe2en2      - Employment by age, economic activity and NUTS 2 regions (NACE Rev. 2) - 1 000
#  lfst_r_lfp2acedu    - Economically active population by sex, age, educational attainment level and NUTS 2 regions (1 000)
#  lfst_r_lfu3pers     - Unemployment by sex, age and NUTS 2 regions (1 000)
#  nama_10r_2gdp       - Gross domestic product (GDP) at current market prices by NUTS 2 regions 
#  nama_10r_2gvagr     - Real growth rate of regional gross value added (GVA) at basic prices by NUTS 2 regions - Percentage change on previous year 
#  nama_10r_2hhinc     - Income of households by NUTS 2 regions 
#  trng_lfse_04        - Participation rate in education and training (last 4 weeks) by NUTS 2 regions
#VARIABLES:
#  Age: age={...}
#  Geo: geo={...}
#  Sex: sex={F; M; T}
#  Unit: unit={...}
#  Age definition: agedef={COMPLET=Completed years ;REACH=Reached during the year}
#  Employment indicator: indic_em={JOBVAC=Number of job vacancies; JOBOCC=Number of occupied jobs; JOBRATE=Job vacancy rate; CH_Y_Y=Job vacancy rate year on year change}
#  Structure of earnings indicator: indic_se={...}
#  ISCED 11: isced11={ED0-2; ED3_4; ED5-8}
#  ISCED-F 13: iscedf13={...}
#  ISCO 08: isco08={OC1-5=Non manual workers; OC1; OC2; OC3; OC4; OC5; OC6-8=Skilled manual workers; OC6; OC7; OC8; OC9=Elementary occupations}
#  Migration status: mgstatus={...}
#  National account indicator: na_item={B5N=Balance of priamry icomes; B6N=Disposable income}
#  NACE Rev.2: nace_r2={...}
#  Size classes in number of employees: sizeclas={TOTAL=Total; GE10=10 employees or more}
#  Working time: worktime={FT=Full-time; PT=Part-time}
#  Time: time={X1990,...,X2017}
files_nuts1 = c("earn_ses_hourly.tsv","educ_uoe_fine06.tsv",
  "educ_uoe_grad02.tsv","lfso_14leeow.tsv","migr_emi2.tsv")
files_nuts2 = c("demo_r_d2jan.tsv","ilc_li41.tsv","ilc_lvhl21.tsv",
  "ilc_lvho04n.tsv","ilc_mddd21.tsv","ilc_peps11.tsv","jvs_a_nace2.tsv",
  "lfst_r_lfe2eedu.tsv","lfst_r_lfe2eftpt.tsv","lfst_r_lfe2emp.tsv",
  "lfst_r_lfe2en2.tsv","lfst_r_lfu3pers.tsv","nama_10r_2gdp.tsv",
  "nama_10r_2gvagr.tsv","nama_10r_2hhinc.tsv","trng_lfse_04.tsv")  
geo_cat = c("AT","AT1","AT11","AT12","AT13","AT2","AT21","AT22","AT3","AT31",
  "AT32","AT33","AT34","BE","BE1","BE10","BE2","BE21","BE22","BE23","BE24",
  "BE25","BE3","BE31","BE32","BE33","BE34","BE35","BG","BG3","BG31","BG32",
  "BG33","BG34","BG4","BG41","BG42","CY","CY0","CY00","CZ","CZ0","CZ01","CZ02",
  "CZ03","CZ04","CZ05","CZ06","CZ07","CZ08","DE","DE1","DE11","DE12","DE13",
  "DE14","DE2","DE21","DE22","DE23","DE24","DE25","DE26","DE27","DE3","DE30",
  "DE4","DE40","DE5","DE50","DE6","DE60","DE7","DE71","DE72","DE73","DE8",
  "DE80","DE9","DE91","DE92","DE93","DE94","DEA","DEA1","DEA2","DEA3","DEA4",
  "DEA5","DEB","DEB1","DEB2","DEB3","DEC","DEC0","DED","DED2","DED4","DED5",
  "DEE","DEE0","DEF","DEF0","DEG","DEG0","DK","DK0","DK01","DK02","DK03","DK04",
  "DK05","EE","EE0","EE00","EL","EL3","EL30","EL4","EL41","EL42","EL43","EL5",
  "EL51","EL52","EL53","EL54","EL6","EL61","EL62","EL63","EL64","EL65","ES",
  "ES1","ES11","ES12","ES13","ES2","ES21","ES22","ES23","ES24","ES3","ES30",
  "ES4","ES41","ES42","ES43","ES5","ES51","ES52","ES53","ES6","ES61","ES62",
  "ES63","ES64","ES7","ES70","FI","FI1","FI19","FI1B","FI1C","FI1D","FI2",
  "FI20","FR","FR1","FR10","FR2","FR21","FR22","FR23","FR24","FR25","FR26",
  "FR3","FR30","FR4","FR41","FR42","FR43","FR5","FR51","FR52","FR53","FR6",
  "FR61","FR62","FR63","FR7","FR71","FR72","FR8","FR81","FR82","FR83","FRA",
  "FRA1","FRA2","FRA3","FRA4","FRA5","HR","HR0","HR03","HR04","HU","HU1","HU10",
  "HU2","HU21","HU22","HU23","HU3","HU31","HU32","HU33","IE","IE0","IE01",
  "IE02","IT","ITC","ITC1","ITC2","ITC3","ITC4","ITF","ITF1","ITF2","ITF3",
  "ITF4","ITF5","ITF6","ITG","ITG1","ITG2","ITH","ITH1","ITH2","ITH3","ITH4",
  "ITH5","ITI","ITI1","ITI2","ITI3","ITI4","LT","LT0","LT00","LU","LU0","LU00",
  "LV","LV0","LV00","MT","MT0","MT00","NL","NL1","NL11","NL12","NL13","NL2",
  "NL21","NL22","NL23","NL3","NL31","NL32","NL33","NL34","NL4","NL41","NL42",
  "PL","PL1","PL11","PL12","PL2","PL21","PL22","PL3","PL31","PL32","PL33",
  "PL34","PL4","PL41","PL42","PL43","PL5","PL51","PL52","PL6","PL61","PL62",
  "PL63","PT","PT1","PT11","PT15","PT16","PT17","PT18","PT2","PT20","PT3",
  "PT30","RO","RO1","RO11","RO12","RO2","RO21","RO22","RO3","RO31","RO32",
  "RO4","RO41","RO42","SE","SE1","SE11","SE12","SE2","SE21","SE22","SE23","SE3",
  "SE31","SE32","SE33","SI","SI0","SI03","SI04","SK","SK0","SK01","SK02","SK03",
  "SK04","UK","UKC","UKC1","UKC2","UKD","UKD1","UKD3","UKD4","UKD6","UKD7",
  "UKE","UKE1","UKE2","UKE3","UKE4","UKF","UKF1","UKF2","UKF3","UKG","UKG1",
  "UKG2","UKG3","UKH","UKH1","UKH2","UKH3","UKI","UKI3","UKI4","UKI5","UKI6",
  "UKI7","UKJ","UKJ1","UKJ2","UKJ3","UKJ4","UKK","UKK1","UKK2","UKK3","UKK4",
  "UKL","UKL1","UKL2","UKM","UKM2","UKM3","UKM5","UKM6","UKN","UKN0")
geo_cat0 = geo_cat[sapply(geo_cat,nchar)==2]
geo_cat1 = geo_cat[sapply(geo_cat,nchar)==3]
geo_cat2 = geo_cat[sapply(geo_cat,nchar)==4]
sex_cat = c("F","M","T")
rem_col = c("indic_se","isced11","na_item","sizeclas","worktime")

#### 2.1. DATA OVERVIEW
{
f1 = "../data/"
cnam1 = overview_data1(files_nuts1,f1=f1)
t1 = table(unlist(lapply(cnam1,colnames)))   #get all colnames
#print(t1)

f1 = "../data/"
cnam2 = overview_data1(files_nuts2,f1=f1)
t2 = table(unlist(lapply(cnam2,colnames)))   #get all colnames
#print(t2)

#overview_data2(names(t1)[1:12],cnam1)
#overview_data2(names(t2)[1:11],cnam2)

}
#### 2.2. READ DATA
{
library("rjson")

f1 = "../results/csv/"
f2 = "../results/json/"
f3 = "../results/dat1.RDS"
if(!file.exists(f3)){
    dat1 = read_data_v2(files_nuts1,geo_cat,sex_cat,rem_col)
    save(dat1,file=f3)
    for(i1 in files_nuts1){
        f4 = paste(f1,gsub(".tsv",".csv",i1),sep="")
        f5 = paste(f2,gsub(".tsv",".json",i1),sep="")
        write.table(dat1[[i1]],file=f4,sep=";",dec=".",row.names=FALSE)
        write(toJSON(as.list(dat1[[i1]])),f5)
    }
}else{
    load(file=f3)
}
#lapply(dat1,dim)
#lapply(dat1,colnames)
#lapply(dat1,summary)
#t3 = table(unlist(lapply(dat1,colnames)))
#print(t3)
#for(i1 in names(t3)[2:8]){
#    print(lapply(dat1,function(x){
#      if(i1 %in% colnames(x) & i1 !="geo"){table(x[,"geo"],x[,i1])}}))
#}

f1 = "../results/csv/"
f2 = "../results/json/"
f3 = "../results/dat2.RDS"
if(!file.exists(f3)){
    dat2 = read_data_v2(files_nuts2,geo_cat,sex_cat,rem_col)
    save(dat2,file=f3)
    for(i1 in files_nuts2){
        f4 = paste(f1,gsub(".tsv",".csv",i1),sep="")
        f5 = paste(f2,gsub(".tsv",".json",i1),sep="")
        write.table(dat2[[i1]],file=f4,sep=";",dec=".",row.names=FALSE)
        write(toJSON(as.list(dat2[[i1]])),f5)
    }
}else{
    load(file=f3)
}
#lapply(dat2,dim)
#lapply(dat2,colnames)
#lapply(dat2,summary)
#t4 = table(unlist(lapply(dat2,colnames)))
#print(t4)
#for(i1 in names(t4)[c(2:6,9)]){
#    print(lapply(dat2,function(x){
#      if(i1 %in% colnames(x) & i1 !="geo"){table(x[,"geo"],x[,i1])}}))
#}

}
}
## 3. EDIT DATA
{
#### 3.1. RESHAPE TO WIDE FORMAT
{
dat1w = reshape2wide_v2(dat1)
#lapply(dat1w,dim)
#lapply(dat1w,colnames)
#table(unlist(lapply(dat1w,colnames)))

dat2w = reshape2wide_v2(dat2)
#lapply(dat2w,dim)
#lapply(dat2w,colnames)
#table(unlist(lapply(dat2w,colnames)))

}
#### 3.2. PREPARE VARIABLES
# y1: Measure mismatch (skills demand vs skills supplied by education)
#educ_uoe_grad02     - Graduates by education level, programme orientation, sex and field of education
#jvs_a_nace2         - Job vacancy statistics by occupation, NUTS 2 regions and NACE Rev. 2 activity - annual data (2008-2015) 
# y2: Migration status of employees
#lfso_14leeow        - Employees by migration status, educational attainment level, occupation and working time
# y3: Emigration
#migr_emi2           - Emigration by age and sex
# y4: cluster between countries
# xs: Attractiveness of labour market NUTS2
#1 demo_r_d2jan      - Population on 1 January by age, sex and NUTS 2 region
#2 ilc_li41          - At-risk-of-poverty rate by NUTS 2 regions 
#3 ilc_lvhl21        - People living in households with very low work intensity by NUTS 2 regions (population aged 0 to 59 years) 
#4 ilc_lvho04n       - Average number of rooms per person by NUTS 2 region 
#5 ilc_mddd21        - Severe material deprivation rate by NUTS 2 regions 
#6 ilc_peps11        - People at risk of poverty or social exclusion by NUTS 2 regions 
#7 lfst_r_lfe2eedu   - Employment by sex, age, educational attainment level and NUTS 2 regions (1000)
#8 lfst_r_lfe2eftpt  - Employment by full-time/part-time, sex and NUTS 2 regions (1000)
#9 lfst_r_lfe2emp    - Employment by sex, age and NUTS 2 regions (1000)
#10 lfst_r_lfe2en2   - Employment by age, economic activity and NUTS 2 regions (NACE Rev. 2) - 1000
#11 lfst_r_lfu3pers  - Unemployment by sex, age and NUTS 2 regions (1000)
#12 nama_10r_2gdp     - Gross domestic product (GDP) at current market prices by NUTS 2 regions 
#13 nama_10r_2gvagr   - Real growth rate of regional gross value added (GVA) at basic prices by NUTS 2 regions - Percentage change on previous year 
#14 nama_10r_2hhinc   - Income of households by NUTS 2 regions 
#15 trng_lfse_04      - Participation rate in education and training (last 4 weeks) by NUTS 2 regions
# xs: Attractiveness of labour market NUTS1
#1 earn_ses_hourly   - Structure of earnings survey: hourly earnings
#2 educ_uoe_fine06   - Total public expenditure on education by education level and programme orientation - as % of GDP
{
###### 3.2.1 Predictors at NUTS2-level

######## 3.2.1.1. Merge datasets
s1 = dat2w[[1]][,"sex"] == "T"
x1 = dat2w[[1]][s1,c(3:9)]; colnames(x1) = c("geo","pop_Total","pop_Y0-14",
  "pop_Y15-24","pop_Y25-64","pop_Y65-74","pop_YGE75")
s2 = dat2w[[2]][,"sex"] == "T"
x2 = dat2w[[2]][s2,4:5]; colnames(x2) = c("geo","ARPR")
s3 = dat2w[[3]][,"sex"] == "T"
x3 = dat2w[[3]][s3,4:5]; colnames(x3) = c("geo","low_work")
s4 = dat2w[[4]][,"sex"] == "T"
x4 = dat2w[[4]][s4,4:5]; colnames(x4) = c("geo","rooms_pp")
s5 = dat2w[[5]][,"sex"] == "T"
x5 = dat2w[[5]][s5,4:5]; colnames(x5) = c("geo","mat_depriv")
s6 = dat2w[[6]][,"sex"] == "T"
x6 = dat2w[[6]][s6,4:5]; colnames(x6) = c("geo","ARPR_socexcl")
s7 = dat2w[[8]][,"sex"] == "T"
x7 = dat2w[[8]][s7,3:9]; colnames(x7) = c("geo","emp_Y15-24_ED0-2",
  "emp_Y25-64_ED0-2","emp_Y15-24_ED3-4","emp_Y25-64_ED3-4",
  "emp_Y15-24_ED5-8","emp_Y25-64_ED5-8")
s8 = dat2w[[9]][,"sex"] == "T"
x8 = dat2w[[9]][s8,4:6]; colnames(x8) = c("geo","emp_TF","emp_TP")
s9 = dat2w[[10]][,"sex"] == "T"
x9 = dat2w[[10]][s9,3:6]; colnames(x9) = c("geo","emp_Y15-24","emp_Y25-64",
  "emp_YGE65")
s10 = dat2w[[11]][,"sex"] == "T"
x10 = dat2w[[11]][s10,3:23]; colnames(x10) = c("geo","emp_Y15-24_NaceA",
  "emp_Y25-64_NaceA","emp_Y15-24_NaceB-E","emp_Y25-64_NaceB-E",
  "emp_Y15-24_NaceF","emp_Y25-64_NaceF","emp_Y15-24_NaceG-I",
  "emp_Y25-64_NaceG-I","emp_Y15-24_NaceJ","emp_Y25-64_NaceJ",
  "emp_Y15-24_NaceK","emp_Y25-64_NaceK","emp_Y15-24_NaceL",
  "emp_Y25-64_NaceL","emp_Y15-24_NaceM_N","emp_Y25-64_NaceM_N",
  "emp_Y15-24_NaceO-Q","emp_Y25-64_NaceO-Q","emp_Y15-24_NaceR-U",
  "emp_Y25-64_NaceR-U")
s11 = dat2w[[12]][,"sex"] == "T"
x11 = dat2w[[12]][s11,3:5]; colnames(x11) = c("geo","unemp_Y15-24",
  "unemp_YGE25")
s12 = dat2w[[13]][,"sex"] == "T"
x12 = dat2w[[13]][s12,4:5]; colnames(x12) = c("geo","GDP")
s13 = dat2w[[14]][,"sex"] == "T"
x13 = dat2w[[14]][s13,4:5]; colnames(x13) = c("geo","GVAgr")
s14 = dat2w[[15]][,"sex"] == "T"
x14 = dat2w[[15]][s14,4:5]; colnames(x14) = c("geo","disp_income")
s15 = dat2w[[16]][,"sex"] == "T"
x15 = dat2w[[16]][s15,4:5]; colnames(x15) = c("geo","training")
l1 = list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
x_all1 = Reduce(function(dtf1, dtf2){merge(dtf1,dtf2,by="geo",all=TRUE)},l1)
w_all1 = c(1,rep(1/5,5),1,1,1,1,1,rep(1/6,6),rep(1/2,2),rep(1/3,3),rep(1/20,20),
  rep(1/2,2),1,1,1,1)

######## 3.2.1.2. Normalize for population size accross regions
s1 = c(seq(13,17,2),                    #emp_Y15-24_ED[?]
       seq(24,42,2))                    #emp_Y15-24_Nace[?]
a1 = x_all1[,21]                        #emp_Y15-24
x_all1[,s1] = apply(x_all1[,s1],2,function(x){100*(x)/a1}) 

s2 = c(seq(14,18,2),                    #emp_Y25-64_ED[?]
       seq(25,43,2))                    #emp_Y25-64_Nace[?]
a2 = x_all1[,22]                        #emp_Y25-64
x_all1[,s2] = apply(x_all1[,s2],2,function(x){100*(x)/a2})

s3 = c(21,                              #emp_Y15-24
       44)                              #unemp_Y15-24
a3 = x_all1[,4]                         #pop_Y15-24
x_all1[,s3] = apply(x_all1[,s3],2,function(x){100*(x*1000)/a3})

s4 = c(22)                              #emp_Y25-64
a4 = x_all1[,5]                         #pop_Y25-64
x_all1[,s4] = 100*(x_all1[,s4]*1000)/a4

s5 = 45                                 #unemp_Y25-74
a5 = apply(x_all1[,5:6],1,sum)          #pop_Y25-74
x_all1[,s5] = 100*(x_all1[,s5]*1000)/a5

s6 = c(19:20)                           #emp_T[?]
a6 = apply(x_all1[,4:6],1,sum)          #pop_YGE15
x_all1[,s6] = apply(x_all1[,s6],2,function(x){100*(x*1000)/a6})

s7 = 23                                 #emp_YGE65
a7 = apply(x_all1[,6:7],1,sum)          #pop_YGE65
x_all1[,s7] = 100*(x_all1[,s7]*1000)/a7

s8 = c(3:7)                             #pop_Y[?]
a8 = x_all1[,2]                         #pop_Total
x_all1[,s8] = apply(x_all1[,s8],2,function(x){100*(x)/a8})

###### 3.2.2. Predictors at NUTS0-level
s1 = match(substr(x_all1[,"geo"],1,2),dat1w[[1]][,"geo"]) #map NUTS0 to x_all
x_all2 = dat1w[[1]][s1,5:31]; colnames(x_all2) = c("earn_OC1_NaceB-F",
  "earn_OC2_NaceB-F","earn_OC3_NaceB-F","earn_OC4_NaceB-F","earn_OC5_NaceB-F",
  "earn_OC6_NaceB-F","earn_OC7_NaceB-F","earn_OC8_NaceB-F","earn_OC9_NaceB-F",
  "earn_OC1_NaceG-N","earn_OC2_NaceG-N","earn_OC3_NaceG-N","earn_OC4_NaceG-N",
  "earn_OC5_NaceG-N","earn_OC6_NaceG-N","earn_OC7_NaceG-N","earn_OC8_NaceG-N",
  "earn_OC9_NaceG-N","earn_OC1_NaceP-S","earn_OC2_NaceP-S","earn_OC3_NaceP-S",
  "earn_OC4_NaceP-S","earn_OC5_NaceP-S","earn_OC6_NaceP-S","earn_OC7_NaceP-S",
  "earn_OC8_NaceP-S","earn_OC9_NaceP-S")
w_all2 = rep(1/27,27)

s2 = match(substr(x_all1[,"geo"],1,2),dat1w[[2]][,"geo"]) #map NUTS0 to x_all
x_all3 = dat1w[[2]][s2,5]
w_all3 = 1

###### 3.2.3. All predictors
x_all = cbind(x_all1,x_all2,x_all3)
colnames(x_all) = c(colnames(x_all)[-ncol(x_all)],"expend_ED5-8")
rownames(x_all) = x_all[,1]
x_all = x_all[,-1]
w_all = c(w_all1,w_all2,w_all3)

f1 = "../results/csv/lmktattract.csv"
f2 = "../results/json/lmktattract.json"
write.table(x_all,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(x_all)),f2)
x_all_nuts0 =  x_all[rownames(x_all) %in% geo_cat0,] #nrows = 28
x_all_nuts1 =  x_all[rownames(x_all) %in% geo_cat1,] #nrows = 98
x_all_nuts2 =  x_all[rownames(x_all) %in% geo_cat2,] #nrows = 276

###### 3.2.4. Response variable (Skills mismatch)

######## 3.2.4.1. Skills supply
s1 = dat1w[["educ_uoe_grad02.tsv"]][,"sex"] == "T"
tmp1 = as.matrix(dat1w[["educ_uoe_grad02.tsv"]][s1,5:15])
tmp2 = as.matrix(dat1w[["educ_uoe_grad02.tsv"]][s1,16,drop=FALSE])
rownames(tmp1) = dat1w[["educ_uoe_grad02.tsv"]][s1,"geo"]
rownames(tmp2) = dat1w[["educ_uoe_grad02.tsv"]][s1,"geo"]
y1_educ = matrix(NA,nrow=length(geo_cat0),ncol=ncol(tmp1))
rownames(y1_educ) = geo_cat0
colnames(y1_educ) = c("F00","F01","F02","F03","F04","F05","F06","F07","F08",
  "F09","F10")
for(i1 in geo_cat0){
    if(i1 %in% rownames(tmp1)){
        y1_educ[i1,] = 100*tmp1[i1,]/tmp2[i1,1]
    }
}
s1 = match(substr(rownames(x_all),1,2),rownames(y1_educ)) #map NUTS0 to x_all
y1_educ = y1_educ[s1,]
rownames(y1_educ) = rownames(x_all)

f1 = "../results/csv/lmktsupply.csv"
f2 = "../results/json/lmktsupply.json"
write.table(y1_educ,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(y1_educ)),f2)
y1_educ_nuts0 =  y1_educ[rownames(y1_educ) %in% geo_cat0,] #nrows = 28
y1_educ_nuts1 =  y1_educ[rownames(y1_educ) %in% geo_cat1,] #nrows = 98
y1_educ_nuts2 =  y1_educ[rownames(y1_educ) %in% geo_cat2,] #nrows = 276

######## 3.2.4.2. Skills demand
onlyNUTS0 = FALSE
map1 = read.csv("../data/ESCO_v3.csv",stringsAsFactors=FALSE)
tmp1 = dat2w[["jvs_a_nace2.tsv"]][,c(5:14,27:36,sapply(seq(49,225,11),
  function(x){x:(x+9)}))]                         #{"A",...,"S"} per {"OC0",...,"OC9"}
tmp2 = dat2w[["jvs_a_nace2.tsv"]][,c(16:25)]      #{"A-S"} per {"OC0",...,"OC9"}
tmp3 = dat2w[["jvs_a_nace2.tsv"]][,c(38:47)]      #{"B-S"} per {"OC0",...,"OC9"}
tmp4 = dat2w[["jvs_a_nace2.tsv"]][,26,drop=FALSE] #{"Total"}
rownames(tmp1) = rownames(tmp2) = rownames(tmp3) = rownames(tmp4) = 
  dat2w[["jvs_a_nace2.tsv"]][,"geo"]
colnames(tmp1) = gsub("_nace_","_",gsub("X2014_indicem_JOBTOT_isco08_","",
  colnames(tmp1)))
colnames(tmp2) = gsub("_nace_","_",gsub("X2014_indicem_JOBTOT_isco08_","",
  colnames(tmp2)))
colnames(tmp3) = gsub("_nace_","_",gsub("X2014_indicem_JOBTOT_isco08_","",
  colnames(tmp3)))
colnames(tmp4) = gsub("_nace_","_",gsub("X2014_indicem_JOBTOT_isco08_","",
  colnames(tmp4)))
if(!onlyNUTS0){
    tmp1["LV00",] = tmp1["LV0",]        #correct missing info in LV NUTS3
    tmp2["LV00",] = tmp2["LV0",]        #correct missing info in LV NUTS3
    tmp3["LV00",] = tmp3["LV0",]        #correct missing info in LV NUTS3
    tmp4["LV00",] = tmp4["LV0",]        #correct missing info in LV NUTS3
    y1_jobv = mapSCQO(tmp1,tmp2,tmp3,tmp4,map1[,3:4],geo_cat,usetot=FALSE)
}else{
    y1_jobv = mapSCQO(tmp1,tmp2,tmp3,tmp4,map1[,3:4],usetot=FALSE)
    s1 = match(substr(rownames(x_all),1,2),rownames(y1_jobv)) #map NUTS0 to x_all
    y1_jobv = y1_jobv[s1,]
    rownames(y1_jobv) = rownames(x_all)
}

f1 = "../results/csv/lmktdemand.csv"
f2 = "../results/json/lmktdemand.json"
write.table(y1_jobv,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(y1_jobv)),f2)
y1_jobv_nuts0 =  y1_jobv[rownames(y1_jobv) %in% geo_cat0,] #nrows = 28
y1_jobv_nuts1 =  y1_jobv[rownames(y1_jobv) %in% geo_cat1,] #nrows = 98
y1_jobv_nuts2 =  y1_jobv[rownames(y1_jobv) %in% geo_cat2,] #nrows = 276

######## 3.2.4.3. Skills mismatch
# Note: Distances
# 1. "euclidean", sqrt(sum((x_i - y_i)^2))
# 2. "maximum",   max(x_i - y_i)
# 3. "manhattan", sum(abs(x_i - y_i))
# 4. "minkowski", (sum((x_i - y_i)^p))^(1/p)
y1 = matrix(NA,nrow=length(geo_cat),ncol=1)
rownames(y1) = geo_cat
colnames(y1) = c("mismatch")
for(i1 in geo_cat){
    y1[i1,1] = dist(rbind(y1_jobv[i1,],y1_educ[i1,]),method="euclidean")
}

f1 = "../results/csv/lmkmismatch.csv"
f2 = "../results/json/lmktmismatch.json"
write.table(y1,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(y1)),f2)
y1_nuts0 =  y1[rownames(y1) %in% geo_cat0,,drop=FALSE] #nrows = 28
y1_nuts1 =  y1[rownames(y1) %in% geo_cat1,,drop=FALSE] #nrows = 98
y1_nuts2 =  y1[rownames(y1) %in% geo_cat2,,drop=FALSE] #nrows = 276

###### 3.2.5. Response variable (Mobility)
tmp1 = as.matrix(dat1w[["lfso_14leeow.tsv"]][,c(seq(8,44,4),45:47,48)])
tmp2 = as.matrix(dat1w[["lfso_14leeow.tsv"]][,c(seq(52,88,4),89:91,92)])
rownames(tmp1) = dat1w[["lfso_14leeow.tsv"]][,"geo"]
rownames(tmp2) = dat1w[["lfso_14leeow.tsv"]][,"geo"]
y2 = matrix(NA,nrow=length(geo_cat0),ncol=ncol(tmp1))
rownames(y2) = geo_cat0
colnames(y2) = paste("lmktm_",c(paste("OC",0:9,sep=""),"ED0-2","ED3_4","ED5-8",
  "Total"),sep="")
for(i1 in geo_cat0){
    if(i1 %in% rownames(tmp1)){
        y2[i1,] = 100*tmp1[i1,]/tmp2[i1,]
    }
}
s1 = match(substr(rownames(x_all),1,2),rownames(y2)) #map NUTS0 to x_all
y2 = y2[s1,]
rownames(y2) = rownames(x_all)

f1 = "../results/csv/lmktmobil.csv"
f2 = "../results/json/lmktmobil.json"
write.table(y2,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(y2)),f2)
y2_nuts0 =  y2[rownames(y2) %in% geo_cat0,] #nrows = 28
y2_nuts1 =  y2[rownames(y2) %in% geo_cat1,] #nrows = 98
y2_nuts2 =  y2[rownames(y2) %in% geo_cat2,] #nrows = 276

###### 3.2.6. Response variable (Emigration rate)
s1 = dat1w[["migr_emi2.tsv"]][,"sex"] == "T"
tmp1 = as.matrix(dat1w[["migr_emi2.tsv"]][s1,c(4,6:7)])
tmp2 = as.matrix(dat1w[["migr_emi2.tsv"]][s1,c(10,12:13)])
rownames(tmp1) = dat1w[["migr_emi2.tsv"]][s1,"geo"]
rownames(tmp2) = dat1w[["migr_emi2.tsv"]][s1,"geo"]
y3 = matrix(NA,nrow=length(geo_cat0),ncol=ncol(tmp1))
rownames(y3) = geo_cat0
colnames(y3) = paste("migrt_",c("Total","Y15-24","Y25-64"),sep="")
for(i1 in geo_cat0){
    if(i1 %in% rownames(tmp1)){
        s1 = dat2w[["demo_r_d2jan.tsv"]][,"geo"] == i1 & 
          dat2w[["demo_r_d2jan.tsv"]][,"sex"] == "T"
        v1 = as.matrix(dat2w[["demo_r_d2jan.tsv"]][s1,c(4,6:7)])
        y3[i1,1] = 100*tmp1[i1,1]/v1[1]     #Normalize for pop_Total
        y3[i1,2] = 100*tmp1[i1,2]/v1[2]     #Normalize for pop_Y15-24
        y3[i1,3] = 100*tmp1[i1,3]/v1[3]     #Normalize for pop_Y25-64
        isna = any(is.na(y3[i1,2:3]) | y3[i1,2:3] == 0)
        if(isna){
            y3[i1,2] = 100*tmp2[i1,2]/v1[2] #Normalize for pop_Y15-24
            y3[i1,3] = 100*tmp2[i1,3]/v1[3] #Normalize for pop_Y25-64
        }
    }
}
s1 = match(substr(rownames(x_all),1,2),rownames(y3)) #map NUTS0 to x_all
y3 = y3[s1,]
rownames(y3) = rownames(x_all)

f1 = "../results/csv/migrate.csv"
f2 = "../results/json/migrate.json"
write.table(y3,file=f1,sep=";",dec=".",row.names=TRUE,col.names=NA)
write(toJSON(as.list(y3)),f2)
y3_nuts0 =  y3[rownames(y3) %in% geo_cat0,] #nrows = 28
y3_nuts1 =  y3[rownames(y3) %in% geo_cat1,] #nrows = 98
y3_nuts2 =  y3[rownames(y3) %in% geo_cat2,] #nrows = 276

}
}
## 4. ANALYSE DATA
{
col0 = c(rbind(hsv(h=seq(0,1,len=8),s=1/3,v=5/5)[1:7],
               hsv(h=seq(0,1,len=8),s=2/3,v=5/5)[1:7],
               hsv(h=seq(0,1,len=8),s=3/3,v=5/5)[1:7],
               hsv(h=seq(0,1,len=8),s=3/3,v=4/5)[1:7]))
names(col0) = geo_cat0                  #colours for EU28
col1 = col0[match(substring(geo_cat1,1,2),names(col0))]
names(col1) = geo_cat1
col2 = col0[match(substring(geo_cat2,1,2),names(col0))]
names(col2) = geo_cat2
xunits = c("NR",rep("PC_POP",6),"PC_YLE60","AVG",rep("PC_POP",2),
  rep(c("PC_EMP_Y15-24","PC_EMP_Y25-64"),3),rep("PC_YGE15",2),"PC_POP_Y15-24",
  "PC_POP_Y25-64","PC_POP_YGE65",rep(c("PC_EMP_Y15-24","PC_EMP_Y25-64"),10),
  "PC_YGE15-LE24","PC_YGE25-LE74","PPS_HAB","PCH_PRE",
  "PPCS_HAB","PC_YGE25-LE64",rep("MN_PPS",27),"PC_GDP")
names(xunits) = colnames(x_all)
yunits = c("DIFF",rep("PC_EMP",14),rep("PC_POP",3))
xprop = rep(FALSE,length=ncol(x_all))
names(xprop) = colnames(x_all)
xprop[grep("^PC",xunits)] = TRUE
x_all_min = rep(NA,length=ncol(x_all))
x_all_max = rep(NA,length=ncol(x_all))
names(x_all_min) = colnames(x_all)
names(x_all_max) = colnames(x_all)
x_all_min[1] = min(x_all[,1],na.rm=TRUE)         #pop_Total
x_all_max[1] = max(x_all[,1],na.rm=TRUE)         #pop_Total
x_all_min[9] = min(x_all[,9],na.rm=TRUE)         #rooms_pp
x_all_max[9] = max(x_all[,9],na.rm=TRUE)         #rooms_pp
x_all_min[45] = min(x_all[,45],na.rm=TRUE)       #GDP
x_all_max[45] = max(x_all[,45],na.rm=TRUE)       #GDP
x_all_min[47] = min(x_all[,47],na.rm=TRUE)       #disp_income
x_all_max[47] = max(x_all[,47],na.rm=TRUE)       #disp_income
x_all_min[49:75] = min(x_all[,49:75],na.rm=TRUE) #earn_???_???
x_all_max[49:75] = max(x_all[,49:75],na.rm=TRUE) #earn_???_???

#### 4.1. NUTS0-LEVEL (COUNTRY-LEVEL)
{
y_all_nuts0 = cbind(y1_nuts0,y2_nuts0,y3_nuts0)
w_all_nuts0 = w_all

###### 4.1.1. Summary statistics
f1 = "../results/data_overview/nuts0/descr_x_nuts0.csv"
describe_data_2(x_all_nuts0,f1)
f1 = "../results/data_overview/nuts0/descr_y_nuts0.csv"
describe_data_2(y_all_nuts0,f1)

###### 4.1.2. Scatterplots
f1 = "../results/data_overview/nuts0/descr_x_nuts0"
plot_data(x_all_nuts0,col1=col0,units1=xunits,f1=f1,rm_ylab=FALSE,rm_main=FALSE)

f1 = "../results/data_overview/nuts0/descr_y_nuts0"
plot_data(y_all_nuts0,col1=col0,units1=yunits,f1=f1,rm_ylab=FALSE,rm_main=FALSE)

###### 4.1.3. Barplots
f1 = "../results/data_overview/nuts0/lmkt_mismatch_nuts0"
plot_barplot(t(y1_educ_nuts0)[-1,],t(y1_jobv_nuts0)[-1,],ylab1="Skills supply",
  ylab2="Skills demand",main="NUTS0",f1=f1)
  
###### 4.1.4. Distances between regions
a1 = t(apply(scale(x_all_nuts0),1,function(x){x*w_all_nuts0}))
x_dist_nuts0 = dist(a1,method="euclidean")

###### 4.1.5. Correlation

######## 4.1.5.1. Correlation between predictors
f1 = "../results/data_overview/nuts0/cor_x_nuts0"
x_cor_nuts0 = calc_cor(x_all_nuts0,f1=f1)
plot_levelplots_cor_v2(x_cor_nuts0,reord=FALSE,rm_xlab=TRUE,rm_ylab=TRUE,f1=f1,
  main="NUTS0")

######## 4.1.5.2 Correlation between response and predictors
f1 = "../results/data_overview/nuts0/cor_xy_nuts0"
xy_cor_nuts0 = calc_cor_v2(y_all_nuts0,x_all_nuts0,f1=f1)
plot_levelplots_cor_v2(xy_cor_nuts0,stat="P",reord=FALSE,mirror=FALSE,
  rm_xlab=FALSE,rm_ylab=FALSE,f1=f1,main="NUTS0")

###### 4.1.6. Classification analysis
f1 = "../results/data_overview/nuts0/net_nuts0"
plot_sna_v2(x_dist_nuts0,col1=col0,edge_thr=0.65,f1=f1,
  main="Labour market attractiveness")
analyse_net(x_dist_nuts0,f1=f1)
f1 = "../results/data_overview/nuts0/pam_nuts0"
plot_pam_v2(x_dist_nuts0,col1=col0,f1=f1,main="Labour market attractiveness")
kbest = 10
f1 = paste("../results/data_overview/nuts0/pam_nuts0_grps",kbest,sep="")
plot_pam_v2(x_dist_nuts0,kbest=kbest,col1=col0,f1=f1,
  main="Labour market attractiveness")

###### 4.1.7. Study classification analysis
f1 = "../results/data_overview/nuts0/pam_nuts0_grps10_3.csv"
a1 = read.csv(f1)
clst1 = a1[,"cluster"]
clst1 = factor(clst1,levels=sort(unique(clst1)))
names(clst1) = a1[,"X"]
y4_nuts0 = clst1 = clst1[rownames(y_all_nuts0)]
y_all_nuts0 = cbind(data.frame(y_all_nuts0),y4_nuts0)
colnames(y_all_nuts0) = c(colnames(y_all_nuts0)[-ncol(y_all_nuts0)],"EU_groups")

######## 4.1.7.1. Descriptive statistics

########## 4.1.7.1.1. Summary statistics on clusters
f1 = "../results/analyse_classification/nuts0/summs_nuts0"
for(i1 in levels(clst1)){
    s1 = clst1 == i1
    f2 = paste(f1,"_clst",i1,".csv",sep="")
    describe_data_2(x_all_nuts0[s1,],f1=f2)
}

########## 4.1.7.1.2. Summary statistics on variables
f1 = "../results/analyse_classification/nuts0/summs_nuts0_vars"
f2 = paste(f1,"_vb.csv",sep="")
f3 = paste(f1,"_vw.csv",sep="")
f4 = paste(f1,"_sep.csv",sep="")
f5 = paste(f1,"_cb.csv",sep="")
f6 = paste(f1,"_cw.csv",sep="")
f7 = paste(f1,"_cor.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5) |
  !file.exists(f6) | !file.exists(f7)){
    describe_vars_v2(x_all_nuts0,clst1,f1)
}

########## 4.1.7.1.3. T-tests on variables
f1 = "../results/analyse_classification/nuts0/ttest_nuts0"
calc_pairwise_ttest_v2(x_all_nuts0,clst1,stat="pval",thr=0.05,f1=f1)

########## 4.1.7.1.4. Boxplots
f1 = "../results/analyse_classification/nuts0/boxplot_nuts0"
boxplot_data_v2(x_all_nuts0,clst1,units1=xunits,f1=f1,k1=8,k2=2,rm_ylab=FALSE,
  rm_main=FALSE)

######## 4.1.7.2. Modeling
f1 = "../results/analyse_classification/nuts0/datsel_nuts0"
f2 = "../results/analyse_classification/nuts0/datred_nuts0.log"
f3 = paste(f1,".log",sep="")
ina = apply(x_all_nuts0,2,function(x){any(is.na(x))})
x_complete = x_all_nuts0[,!ina]
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=FALSE,type="output")
x_complete = reduce_predictors_v2(x_complete,clst1,thr1=0.90,thr2=NULL,
  thr3=0.00,thr4=Inf,maxsize=30,f1=f2)
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=TRUE,type="output")
f4 = paste(f1,".csv",sep="")
write.table(t(apply(x_complete,2,summary)),file=f4,sep=",",dec=".",
  row.names=TRUE,col.names=NA)
x_scaled = x_complete
s1 = match(colnames(x_scaled),colnames(x_all))
iprop = xprop[s1]
x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,function(x){x/100}))
for(i1 in 1:ncol(x_scaled)){
    if(!iprop[i1]){
        xmin = x_all_min[s1][i1]
        xmax = x_all_max[s1][i1]
        x_scaled[,i1] = minmax(x_scaled[,i1],xmin,xmax)
    }
}
f5 = paste(f1,"_scaled.csv",sep="")
write.table(t(apply(x_scaled,2,summary)),file=f5,sep=",",dec=".",
  row.names=TRUE,col.names=NA)

########## 4.1.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
if(length(levels(clst1)) <= 2){
    f1 = "../results/analyse_classification/nuts0/modsel_binom_nuts0"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_binom(x_complete,clst1,thr1=1.00,thr2=0.00,thr3=0.00,
          thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts0/fit_binom_nuts0"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_binomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS0")
    }
}else{
    f1 = "../results/analyse_classification/nuts0/modsel_multinom_nuts0"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_multinom(x_scaled,clst1,thr1=1.00,thr2=0.00,
          thr3=0.00,thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts0/fit_multinom_nuts0"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_multinomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS0")
    }
}

###### 4.1.8. Enrichment analysis
for(i1 in c(4,5,10)){
    thr1 = (i1 - 1)/i1
    f1 = paste("../results/analyse_classification/nuts0/ora_ge",round(thr1*100),
      "_nuts0",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts0,prob=thr1,dirct="ge",
          stat="pval",f1=f2)
    }
    thr2 = 1/i1
    f1 = paste("../results/analyse_classification/nuts0/ora_le",round(thr2*100),
      "_nuts0",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts0,prob=thr2,dirct="le",
          stat="pval",f1=f2)
    }
}

###### 4.1.9. Weighted correlation network analysis (WCNA)
f1 = "../results/wcna/nuts0/datsel_nuts0.log"
thr = 0.50
srow = apply(x_all_nuts0,1,function(x){round(sum(is.na(x))/length(x),3)})
trow = table(srow)
#print(trow)
#plot(trow,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_all_nuts0[srow < thr,]
y_reduced = y_all_nuts0[srow < thr,]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=FALSE,type="output")  
thr = 0.15
scol = apply(x_reduced,2,function(x){round(sum(is.na(x))/length(x),3)})
tcol = table(scol)
#print(tcol)
#plot(tcol,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_reduced[,scol < thr]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=TRUE,type="output")

f1 = "../results/wcna/nuts0/wcna_nuts0"
f2 = paste(f1,"_1.csv",sep="")
f3 = paste(f1,"_2.csv",sep="")
f4 = paste(f1,"_subjInfo.csv",sep="")
f5 = paste(f1,"_varInfo.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5)){
    analyse_WCNA_v2(x_reduced,y_reduced,mins1=2,spl1=4,meth1="hybrid",
      cut1=0.2,pwr=NULL,pwrrg=0.90,f1=f1)
}

###### 4.1.10. Multivariate linear [General LM] and non-linear [ANN] regression
# Note: No optimization for network size! Too few points for successful bootstrap
# optimization leads to very different results between "model selection" and "fit
# model"!
library("caret")

f1 = "../results/regression/nuts0/datsel_nuts0"
f2 = "../results/regression/nuts0/datred_nuts0"
f3 = "../results/regression/nuts0/modsel_lm_nuts0"
f4 = "../results/regression/nuts0/fit_lm_nuts0"
f5 = "../results/regression/nuts0/modsel_nnet_nuts0"
f6 = "../results/regression/nuts0/fit_nnet_nuts0"
opt_nnet = FALSE
nn_size = 1
nn_decay = 1e-4
ina = apply(x_all_nuts0,2,function(x){any(is.na(x))})
x_complete = x_all_nuts0[,!ina]
for(i1 in colnames(y_all_nuts0)){
    y = y_all_nuts0[,i1,drop=FALSE]
    if(is.null(levels(y[,1]))){
        f7 = paste(f1,"_",i1,".log",sep="")
        f8 = paste(f1,"_",i1,".csv",sep="")
        f9 = paste(f2,"_",i1,".log",sep="")
        capture.output(
          dim(x_complete),
          rownames(x_complete),
          colnames(x_complete),
          summary(x_complete),
          file=f7,append=FALSE,type="output")
        ina = is.na(y[,1])
        y = y[!ina,,drop=FALSE]
        x = x_complete[!ina,]
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        x = reduce_predictors_v2(x,y[,1],thr1=0.90,thr2=NULL,thr3=0.00,thr4=Inf,
          maxsize=30,f1=f9)
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        write.table(t(apply(x,2,summary)),file=f8,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f3,"_",i1,sep="")
        f8 = paste(f7,"_mod1.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,"_1.csv",sep="")
        f11 = paste(f7,"_2.csv",sep="")
        if(!file.exists(f9) | !file.exists(f10) | !file.exists(f11)){
            mod1 = mod_select_lm(x,y,thr1=1.00,thr2=0.00,thr3=0.00,thr4=Inf,
              maxsize=NULL,maxs1="hard",confs=100,level=1,thrs2=200000,f1=f7)
            save(mod1,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f4,"_",i1,sep="")
        f8 = paste(f7,".log",sep="")
        f9 = paste(f7,".csv",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_lm(mod1$formula,mod1$data,f1=f7,main="NUTS0")
        }
        f7 = paste(f1,"_",i1,"_scaled.csv",sep="")
        x_scaled = x
        s1 = match(colnames(x_scaled),colnames(x_all))
        iprop = xprop[s1]
        x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,
          function(x){x/100}))
        for(i2 in 1:ncol(x_scaled)){
            if(!iprop[i2]){
                xmin = x_all_min[s1][i2]
                xmax = x_all_max[s1][i2]
                x_scaled[,i2] = minmax(x_scaled[,i2],xmin,xmax)
            }
        }
        write.table(t(apply(x_scaled,2,summary)),file=f7,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f5,"_",i1,sep="")
        if(opt_nnet){
            f8 = paste(f7,"_optnnet.log",sep="")
            vars = all.vars(mod1$formula)
            dat1 = x_scaled
            colnames(dat1) = sapply(colnames(dat1),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            colnames(y) = sapply(colnames(y),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            dat1 = as.data.frame(cbind(y,dat1[,vars[-1]]))
            grid1 = expand.grid(size=c(1,3*(1:10)),decay=1e-4)
            train1 = train(mod1$formula,dat1,method="nnet",abstol=1e-6,
              reltol=1e-12,linout=TRUE,trace=FALSE,tuneGrid=grid1,metric="RMSE",
              trControl=trainControl(method="boot",number=100))
            nn_size = train1$bestTune[1,1]
            nn_decay = train1$bestTune[1,2]
            capture.output(train1,file=f8,append=FALSE,type="output")
        }
        f8 = paste(f7,"_mod2.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,".csv",sep="")
        if(!file.exists(f9) | !file.exists(f10)){
            mod2 = mod_select_nnet_v2(x_scaled,y,thr1=1.00,thr2=0.00,thr3=0.00,
              thr4=Inf,maxsize=NULL,maxs1="hard",nn_size=nn_size,
              nn_decay=nn_decay,nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,
              confs=100,level=1,thrs2=200000,f1=f7)
            save(mod2,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f6,"_",i1,sep="")
        f8 = paste(f7,"_1.log",sep="")
        f9 = paste(f7,"_2.log",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_nnet_v2(mod2$formula,mod2$data,nn_size=nn_size,nn_decay=1e-4,
              nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,main="NUTS0",f1=f7)
        }
    }
}

}
#### 4.2. NUTS1-LEVEL
{
rm_x_NUTS0 = TRUE                       #remove x-vars at NUTS0-level
rm_y_NUTS0 = TRUE                       #remove y-vars at NUTS0-level
if(rm_x_NUTS0){
    x_all_nuts1 = x_all_nuts1[,-c(49:76)]
    w_all_nuts1 = w_all[-c(49:76)]
}else{
    w_all_nuts1 = w_all
}
if(rm_y_NUTS0){
    y_all_nuts1 = y1_nuts1
}else{
    y_all_nuts1 = cbind(y1_nuts1,y2_nuts1,y3_nuts1)
}

###### 4.2.1. Summary statistics
f1 = "../results/data_overview/nuts1/descr_x_nuts1.csv"
describe_data_2(x_all_nuts1,f1)
f1 = "../results/data_overview/nuts1/descr_y_nuts1.csv"
describe_data_2(y_all_nuts1,f1)

###### 4.2.2. Scatterplots
f1 = "../results/data_overview/nuts1/descr_x_nuts1"
legend1 = col0; names(legend1) = geo_cat0
plot_data(x_all_nuts1,col1=col1,units1=xunits,legend1=legend1,f1=f1,
  rm_ylab=FALSE,rm_main=FALSE)

f1 = "../results/data_overview/nuts1/descr_y_nuts1"
legend1 = col0; names(legend1) = geo_cat0
plot_data(y_all_nuts1,col1=col1,units1=yunits,legend1=legend1,f1=f1,
  rm_ylab=FALSE,rm_main=FALSE)

###### 4.2.3. Barplots
f1 = "../results/data_overview/nuts1/lmkt_mismatch_nuts1"
plot_barplot(t(y1_educ_nuts1)[-1,],t(y1_jobv_nuts1)[-1,],border=NA,
  ylab1="Skills supply",ylab2="Skills demand",main="NUTS1",cex_names=0.8,f1=f1)

###### 4.2.4. Distances between regions
a1 = t(apply(scale(x_all_nuts1),1,function(x){x*w_all_nuts1}))
x_dist_nuts1 = dist(a1,method="euclidean")

###### 4.2.5. Correlation

######## 4.2.5.1. Correlation between predictors
f1 = "../results/data_overview/nuts1/cor_x_nuts1"
x_cor_nuts1 = calc_cor(x_all_nuts1,f1=f1)
plot_levelplots_cor_v2(x_cor_nuts1,reord=FALSE,rm_xlab=FALSE,rm_ylab=FALSE,f1=f1,
  main="NUTS1")

######## 4.2.5.2 Correlation between response and predictors
f1 = "../results/data_overview/nuts1/cor_xy_nuts1"
xy_cor_nuts1 = calc_cor_v2(y_all_nuts1,x_all_nuts1,f1=f1)
plot_levelplots_cor_v2(xy_cor_nuts1,stat="padj",reord=FALSE,mirror=FALSE,
  rm_xlab=FALSE,rm_ylab=FALSE,f1=f1,main="NUTS1")

###### 4.2.6. Classification analysis
f1 = "../results/data_overview/nuts1/net_nuts1"
plot_sna_v2(x_dist_nuts1,col1=col1,edge_thr=0.8,label_cex=0.5,f1=f1,
  main="Labour market attractiveness")
analyse_net(x_dist_nuts1,f1=f1)
f1 = "../results/data_overview/nuts1/pam_nuts1"
plot_pam_v2(x_dist_nuts1,col1=col1,f1=f1,main="Labour market attractiveness")
kbest = 21
f1 = paste("../results/data_overview/nuts1/pam_nuts1_grps",kbest,sep="")
plot_pam_v2(x_dist_nuts1,kbest=kbest,col1=col1,f1=f1,
  main="Labour market attractiveness")

###### 4.2.7. Study classification analysis
f1 = "../results/data_overview/nuts1/pam_nuts1_grps21_3.csv"
a1 = read.csv(f1)
clst1 = a1[,"cluster"]
clst1 = factor(clst1,levels=sort(unique(clst1)))
names(clst1) = a1[,"X"]
y4_nuts1 = clst1 = clst1[rownames(y_all_nuts1)]
y_all_nuts1 = cbind(data.frame(y_all_nuts1),y4_nuts1)
colnames(y_all_nuts1) = c(colnames(y_all_nuts1)[-ncol(y_all_nuts1)],"EU_groups")

######## 4.2.7.1. Descriptive statistics

########## 4.2.7.1.1. Summary statistics on clusters
f1 = "../results/analyse_classification/nuts1/summs_nuts1"
for(i1 in levels(clst1)){
    s1 = clst1 == i1
    f2 = paste(f1,"_clst",i1,".csv",sep="")
    describe_data_2(x_all_nuts1[s1,],f1=f2)
}

########## 4.2.7.1.2. Summary statistics on variables
f1 = "../results/analyse_classification/nuts1/summs_nuts1_vars"
f2 = paste(f1,"_vb.csv",sep="")
f3 = paste(f1,"_vw.csv",sep="")
f4 = paste(f1,"_sep.csv",sep="")
f5 = paste(f1,"_cb.csv",sep="")
f6 = paste(f1,"_cw.csv",sep="")
f7 = paste(f1,"_cor.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5) |
  !file.exists(f6) | !file.exists(f7)){
    describe_vars_v2(x_all_nuts1,clst1,f1)
}

########## 4.2.7.1.3. T-tests on variables
f1 = "../results/analyse_classification/nuts1/ttest_nuts1"
calc_pairwise_ttest_v2(x_all_nuts1,clst1,stat="pval",thr=0.05,f1=f1)

########## 4.2.7.1.4. Boxplots
f1 = "../results/analyse_classification/nuts1/boxplot_nuts1"
boxplot_data_v2(x_all_nuts1,clst1,units1=xunits,f1=f1,k1=8,k2=1,rm_ylab=FALSE,
  rm_main=FALSE)

######## 4.2.7.2. Modeling
f1 = "../results/analyse_classification/nuts1/datsel_nuts1"
f2 = "../results/analyse_classification/nuts1/datred_nuts1.log"
f3 = paste(f1,".log",sep="")
x_complete = cbind(x_all_nuts1,clst1)
x_complete = remove_na(x_complete,init=2,plotit=FALSE,verbose=FALSE)
clst2 = x_complete[,"clst1"]
names(clst2) = rownames(x_complete)
clst2 = factor(clst2,levels=sort(unique(as.numeric(as.character(clst2))))) #in case some level was lost
x_complete = x_complete[,-match("clst1",colnames(x_complete))]
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=FALSE,type="output")
x_complete = reduce_predictors_v2(x_complete,clst2,thr1=0.90,thr2=NULL,
  thr3=0.00,thr4=Inf,maxsize=30,f1=f2)
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=TRUE,type="output")
f4 = paste(f1,".csv",sep="")
write.table(t(apply(x_complete,2,summary)),file=f4,sep=",",dec=".",
  row.names=TRUE,col.names=NA)
x_scaled = x_complete
s1 = match(colnames(x_scaled),colnames(x_all))
iprop = xprop[s1]
x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,function(x){x/100}))
for(i1 in 1:ncol(x_scaled)){
    if(!iprop[i1]){
        xmin = x_all_min[s1][i1]
        xmax = x_all_max[s1][i1]
        x_scaled[,i1] = minmax(x_scaled[,i1],xmin,xmax)
    }
}
f5 = paste(f1,"_scaled.csv",sep="")
write.table(t(apply(x_scaled,2,summary)),file=f5,sep=",",dec=".",
  row.names=TRUE,col.names=NA)
  
########## 4.2.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
if(length(levels(clst2)) <= 2){
    f1 = "../results/analyse_classification/nuts1/modsel_binom_nuts1"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_binom(x_complete,clst2,thr1=1.00,thr2=0.00,thr3=0.00,
          thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts1/fit_binom_nuts1"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_binomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS1")
    }
}else{
    f1 = "../results/analyse_classification/nuts1/modsel_multinom_nuts1"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_multinom(x_scaled,clst2,thr1=1.00,thr2=0.00,
          thr3=0.00,thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts1/fit_multinom_nuts1"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_multinomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS1")
    }
}

###### 4.2.8. Enrichment analysis
for(i1 in c(4,5,10)){
    thr1 = (i1 - 1)/i1
    f1 = paste("../results/analyse_classification/nuts1/ora_ge",round(thr1*100),
      "_nuts1",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts1,prob=thr1,dirct="ge",
          stat="pval",f1=f2)
    }
    thr2 = 1/i1
    f1 = paste("../results/analyse_classification/nuts1/ora_le",round(thr2*100),
      "_nuts1",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts1,prob=thr2,dirct="le",
          stat="pval",f1=f2)
    }
}

###### 4.2.9. Weighted correlation network analysis (WCNA)
f1 = "../results/wcna/nuts1/datsel_nuts1.log"
thr = 0.50
srow = apply(x_all_nuts1,1,function(x){round(sum(is.na(x))/length(x),3)})
trow = table(srow)
#print(trow)
#plot(trow,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_all_nuts1[srow < thr,]
y_reduced = y_all_nuts1[srow < thr,]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=FALSE,type="output")  
thr = 0.15
scol = apply(x_reduced,2,function(x){round(sum(is.na(x))/length(x),3)})
tcol = table(scol)
#print(tcol)
#plot(tcol,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_reduced[,scol < thr]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=TRUE,type="output")

f1 = "../results/wcna/nuts1/wcna_nuts1"
f2 = paste(f1,"_1.csv",sep="")
f3 = paste(f1,"_2.csv",sep="")
f4 = paste(f1,"_subjInfo.csv",sep="")
f5 = paste(f1,"_varInfo.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5)){
    analyse_WCNA_v2(x_reduced,y_reduced,mins1=2,spl1=4,meth1="hybrid",
      cut1=0.2,pwr=NULL,pwrrg=0.90,f1=f1)
}

###### 4.2.10. Multivariate linear [General LM] and non-linear [ANN] regression
library("caret")

f1 = "../results/regression/nuts1/datsel_nuts1"
f2 = "../results/regression/nuts1/datred_nuts1"
f3 = "../results/regression/nuts1/modsel_lm_nuts1"
f4 = "../results/regression/nuts1/fit_lm_nuts1"
f5 = "../results/regression/nuts1/modsel_nnet_nuts1"
f6 = "../results/regression/nuts1/fit_nnet_nuts1"
opt_nnet = TRUE
x_complete = remove_na(x_all_nuts1,init=2,plotit=FALSE,verbose=FALSE)
ina1 = !rownames(x_all_nuts1) %in% rownames(x_complete)
ina2 = !colnames(x_all_nuts1) %in% colnames(x_complete)
for(i1 in colnames(y_all_nuts1)){
    y = y_all_nuts1[,i1,drop=FALSE]
    if(is.null(levels(y[,1]))){
        f7 = paste(f1,"_",i1,".log",sep="")
        f8 = paste(f1,"_",i1,".csv",sep="")
        f9 = paste(f2,"_",i1,".log",sep="")
        capture.output(
          dim(x_complete),
          rownames(x_complete),
          colnames(x_complete),
          summary(x_complete),
          file=f7,append=FALSE,type="output")
        ina3 = is.na(y[,1]) 
        y = y[!(ina1 | ina3),,drop=FALSE]
        x = x_all_nuts1[!(ina1 | ina3),!ina2]
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        x = reduce_predictors_v2(x,y[,1],thr1=0.90,thr2=NULL,thr3=0.00,thr4=Inf,
          maxsize=30,f1=f9)
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        write.table(t(apply(x,2,summary)),file=f8,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f3,"_",i1,sep="")
        f8 = paste(f7,"_mod1.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,"_1.csv",sep="")
        f11 = paste(f7,"_2.csv",sep="")
        if(!file.exists(f9) | !file.exists(f10) | !file.exists(f11)){
            mod1 = mod_select_lm(x,y,thr1=1.00,thr2=0.00,thr3=0.00,thr4=Inf,
              maxsize=NULL,maxs1="hard",confs=100,level=1,thrs2=200000,f1=f7)
            save(mod1,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f4,"_",i1,sep="")
        f8 = paste(f7,".log",sep="")
        f9 = paste(f7,".csv",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_lm(mod1$formula,mod1$data,f1=f7,main="NUTS1")
        }
        f7 = paste(f1,"_",i1,"_scaled.csv",sep="")
        x_scaled = x
        s1 = match(colnames(x_scaled),colnames(x_all))
        iprop = xprop[s1]
        x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,
          function(x){x/100}))
        for(i2 in 1:ncol(x_scaled)){
            if(!iprop[i2]){
                xmin = x_all_min[s1][i2]
                xmax = x_all_max[s1][i2]
                x_scaled[,i2] = minmax(x_scaled[,i2],xmin,xmax)
            }
        }
        write.table(t(apply(x_scaled,2,summary)),file=f7,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f5,"_",i1,sep="")
        if(opt_nnet){
            f8 = paste(f7,"_optnnet.log",sep="")
            vars = all.vars(mod1$formula)
            dat1 = x_scaled
            colnames(dat1) = sapply(colnames(dat1),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            colnames(y) = sapply(colnames(y),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            dat1 = as.data.frame(cbind(y,dat1[,vars[-1]]))
            grid1 = expand.grid(size=c(1,3*(1:10)),decay=1e-4)
            train1 = train(mod1$formula,dat1,method="nnet",abstol=1e-6,
              reltol=1e-12,linout=TRUE,trace=FALSE,tuneGrid=grid1,metric="RMSE",
              trControl=trainControl(method="boot",number=100))
            nn_size = train1$bestTune[1,1]
            nn_decay = train1$bestTune[1,2]
            capture.output(train1,file=f8,append=FALSE,type="output")
        }
        f8 = paste(f7,"_mod2.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,".csv",sep="")
        if(!file.exists(f9) | !file.exists(f10)){
            mod2 = mod_select_nnet_v2(x_scaled,y,thr1=1.00,thr2=0.00,thr3=0.00,
              thr4=Inf,maxsize=NULL,maxs1="hard",nn_size=nn_size,
              nn_decay=nn_decay,nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,
              confs=100,level=1,thrs2=200000,f1=f7)
            save(mod2,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f6,"_",i1,sep="")
        f8 = paste(f7,"_1.log",sep="")
        f9 = paste(f7,"_2.log",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_nnet_v2(mod2$formula,mod2$data,nn_size=NULL,nn_decay=1e-4,
              nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,main="NUTS1",f1=f7)
        }
    }
}

}
#### 4.3. NUTS2-LEVEL
{
rm_x_NUTS0 = TRUE                       #remove x-vars at NUTS0-level
rm_y_NUTS0 = TRUE                       #remove y-vars at NUTS0-level
if(rm_x_NUTS0){
    x_all_nuts2 = x_all_nuts2[,-c(49:76)]
    w_all_nuts2 = w_all[-c(49:76)]
}else{
    w_all_nuts2 = w_all
}
if(rm_y_NUTS0){
    y_all_nuts2 = y1_nuts2
}else{
    y_all_nuts2 = cbind(y1_nuts2,y2_nuts2,y3_nuts2)
}

###### 4.3.1. Summary statistics
f1 = "../results/data_overview/nuts2/descr_x_nuts2.csv"
describe_data_2(x_all_nuts2,f1)
f1 = "../results/data_overview/nuts2/descr_y_nuts2.csv"
describe_data_2(y_all_nuts2,f1)

###### 4.3.2. Scatterplots
f1 = "../results/data_overview/nuts2/descr_x_nuts2"
legend1 = col0; names(legend1) = geo_cat0
plot_data(x_all_nuts2,col1=col2,units1=xunits,legend1=legend1,f1=f1,
  rm_ylab=FALSE,rm_main=FALSE)

f1 = "../results/data_overview/nuts2/descr_y_nuts2"
legend1 = col0; names(legend1) = geo_cat0
plot_data(y_all_nuts2[,1,drop=FALSE],col1=col2,units1=yunits,legend1=legend1,
  f1=f1,rm_ylab=FALSE,rm_main=FALSE)

###### 4.3.3. Barplots
f1 = "../results/data_overview/nuts2/lmkt_mismatch_nuts2"
plot_barplot(t(y1_educ_nuts2)[-1,],t(y1_jobv_nuts2)[-1,],border=NA,
  ylab1="Skills supply",ylab2="Skills demand",main="NUTS2",cex_names=0.3,f1=f1)

###### 4.3.4. Distances between regions
a1 = t(apply(scale(x_all_nuts2),1,function(x){x*w_all_nuts2}))
x_dist_nuts2 = dist(a1,method="euclidean")

###### 4.3.5. Correlation

######## 4.3.5.1. Correlation between predictors
f1 = "../results/data_overview/nuts2/cor_x_nuts2"
x_cor_nuts2 = calc_cor(x_all_nuts2,f1=f1)
plot_levelplots_cor_v2(x_cor_nuts2,reord=FALSE,rm_xlab=FALSE,rm_ylab=FALSE,
  f1=f1,main="NUTS2")

######## 4.3.5.2 Correlation between response and predictors
f1 = "../results/data_overview/nuts2/cor_xy_nuts2"
xy_cor_nuts2 = calc_cor_v2(y_all_nuts2,x_all_nuts2,f1=f1)
plot_levelplots_cor_v2(xy_cor_nuts2,stat="padj",reord=FALSE,mirror=FALSE,
  rm_xlab=FALSE,rm_ylab=FALSE,f1=f1,main="NUTS2")

###### 4.3.6. Classification analysis
f1 = "../results/data_overview/nuts2/net_nuts2"
plot_sna_v2(x_dist_nuts2,col1=col2,edge_thr=0.95,label_cex=0.3,f1=f1,
  main="Labour market attractiveness")
analyse_net(x_dist_nuts2,f1=f1)
f1 = "../results/data_overview/nuts2/pam_nuts2"
plot_pam_v2(x_dist_nuts2,col1=col2,border=NA,f1=f1,
  main="Labour market attractiveness")
kbest = 25
f1 = paste("../results/data_overview/nuts2/pam_nuts2_grps",kbest,sep="")
plot_pam_v2(x_dist_nuts2,kbest=kbest,col1=col2,border=NA,f1=f1,
  main="Labour market attractiveness")

###### 4.3.7. Study classification analysis
f1 = "../results/data_overview/nuts2/pam_nuts2_grps25_3.csv"
a1 = read.csv(f1)
clst1 = a1[,"cluster"]
clst1 = factor(clst1,levels=sort(unique(clst1)))
names(clst1) = a1[,"X"]
y4_nuts2 = clst1 = clst1[rownames(y_all_nuts2)]
y_all_nuts2 = cbind(data.frame(y_all_nuts2),y4_nuts2)
colnames(y_all_nuts2) = c(colnames(y_all_nuts2)[-ncol(y_all_nuts2)],"EU_groups")

######## 4.3.7.1. Descriptive statistics

########## 4.3.7.1.1. Summary statistics on clusters
f1 = "../results/analyse_classification/nuts2/summs_nuts2"
for(i1 in levels(clst1)){
    s1 = clst1 == i1
    f2 = paste(f1,"_clst",i1,".csv",sep="")
    describe_data_2(x_all_nuts2[s1,],f1=f2)
}

########## 4.3.7.1.2. Summary statistics on variables
f1 = "../results/analyse_classification/nuts2/summs_nuts2_vars"
f2 = paste(f1,"_vb.csv",sep="")
f3 = paste(f1,"_vw.csv",sep="")
f4 = paste(f1,"_sep.csv",sep="")
f5 = paste(f1,"_cb.csv",sep="")
f6 = paste(f1,"_cw.csv",sep="")
f7 = paste(f1,"_cor.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5) |
  !file.exists(f6) | !file.exists(f7)){
    describe_vars_v2(x_all_nuts2,clst1,f1)
}

########## 4.3.7.1.3. T-tests on variables
f1 = "../results/analyse_classification/nuts2/ttest_nuts2"
calc_pairwise_ttest_v2(x_all_nuts2,clst1,stat="pval",thr=0.05,f1=f1)

########## 4.3.7.1.4. Boxplots
f1 = "../results/analyse_classification/nuts2/boxplot_nuts2"
boxplot_data_v2(x_all_nuts2,clst1,units1=xunits,f1=f1,k1=8,k2=1,rm_ylab=FALSE,
  rm_main=FALSE)

######## 4.3.7.2. Modeling
f1 = "../results/analyse_classification/nuts2/datsel_nuts2"
f2 = "../results/analyse_classification/nuts2/datred_nuts2.log"
f3 = paste(f1,".log",sep="")
x_complete = cbind(x_all_nuts2,clst1)
x_complete = remove_na(x_complete,init=2,plotit=FALSE,verbose=FALSE)
clst2 = x_complete[,"clst1"]
names(clst2) = rownames(x_complete)
clst2 = factor(clst2,levels=sort(unique(as.numeric(as.character(clst2))))) #in case some level was lost
x_complete = x_complete[,-match("clst1",colnames(x_complete))]
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=FALSE,type="output")
x_complete = reduce_predictors_v2(x_complete,clst2,thr1=0.90,thr2=NULL,thr3=0.00,
  thr4=Inf,maxsize=30,f1=f2)
capture.output(
  dim(x_complete),
  rownames(x_complete),
  colnames(x_complete),
  summary(x_complete),
  file=f3,append=TRUE,type="output")
f4 = paste(f1,".csv",sep="")
write.table(t(apply(x_complete,2,summary)),file=f4,sep=",",dec=".",
  row.names=TRUE,col.names=NA)
x_scaled = x_complete
s1 = match(colnames(x_scaled),colnames(x_all))
iprop = xprop[s1]
x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,function(x){x/100}))
for(i1 in 1:ncol(x_scaled)){
    if(!iprop[i1]){
        xmin = x_all_min[s1][i1]
        xmax = x_all_max[s1][i1]
        x_scaled[,i1] = minmax(x_scaled[,i1],xmin,xmax)
    }
}
f5 = paste(f1,"_scaled.csv",sep="")
write.table(t(apply(x_scaled,2,summary)),file=f5,sep=",",dec=".",
  row.names=TRUE,col.names=NA)
 
########## 4.3.7.2.1. Binomial [GLM] and Multinomial [ANN] regression
if(length(levels(clst2)) <= 2){
    f1 = "../results/analyse_classification/nuts2/modsel_binom_nuts2"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_binom(x_complete,clst2,thr1=1.00,thr2=0.00,thr3=0.00,
          thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts2/fit_binom_nuts2"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_binomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS2")
    }
}else{
    f1 = "../results/analyse_classification/nuts2/modsel_multinom_nuts2"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4)){
        mod1 = mod_select_multinom(x_scaled,clst2,thr1=1.00,thr2=0.00,
          thr3=0.00,thr4=Inf,maxsize=NULL,maxs1="hard",confs=100,level=1,f1=f1)
    }
    f1 = "../results/analyse_classification/nuts2/fit_multinom_nuts2"
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    f5 = paste(f1,"_3.csv",sep="")
    if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) |
      !file.exists(f4)){
        fit_multinomial_v2(mod1$formula,mod1$data,f1=f1,main="NUTS2")
    }
}

###### 4.3.8. Enrichment analysis
for(i1 in c(4,5,10)){
    thr1 = (i1 - 1)/i1
    f1 = paste("../results/analyse_classification/nuts2/ora_ge",round(thr1*100),
      "_nuts2",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts2,prob=thr1,dirct="ge",
          stat="pval",f1=f2)
    }
    thr2 = 1/i1
    f1 = paste("../results/analyse_classification/nuts2/ora_le",round(thr2*100),
      "_nuts2",sep="")
    for(i2 in levels(clst1)){
        f2 = paste(f1,"_grp",i2,sep="")
        nuts_uni = names(clst1)
        nuts_sel = nuts_uni[clst1 == i2]
        performORA_v2(nuts_sel,nuts_uni,x_all_nuts2,prob=thr2,dirct="le",
          stat="pval",f1=f2)
    }
}

###### 4.3.9. Weighted correlation network analysis (WCNA)
f1 = "../results/wcna/nuts2/datsel_nuts2.log"
thr = 0.50
srow = apply(x_all_nuts2,1,function(x){round(sum(is.na(x))/length(x),3)})
trow = table(srow)
#print(trow)
#plot(trow,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_all_nuts2[srow < thr,]
y_reduced = y_all_nuts2[srow < thr,]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=FALSE,type="output")  
thr = 0.15
scol = apply(x_reduced,2,function(x){round(sum(is.na(x))/length(x),3)})
tcol = table(scol)
#print(tcol)
#plot(tcol,ylab="# of NAs")
#abline(v=thr,col="red",lty=2)
x_reduced = x_reduced[,scol < thr]
capture.output(
  dim(x_reduced),
  rownames(x_reduced),
  colnames(x_reduced),
  summary(x_reduced),
  file=f1,append=TRUE,type="output")

f1 = "../results/wcna/nuts2/wcna_nuts2"
f2 = paste(f1,"_1.csv",sep="")
f3 = paste(f1,"_2.csv",sep="")
f4 = paste(f1,"_subjInfo.csv",sep="")
f5 = paste(f1,"_varInfo.csv",sep="")
if(!file.exists(f2) | !file.exists(f3) | !file.exists(f4) | !file.exists(f5)){
    analyse_WCNA_v2(x_reduced,y_reduced,mins1=2,spl1=4,meth1="hybrid",
      cut1=0.2,pwr=NULL,pwrrg=0.90,f1=f1)
}

###### 4.3.10. Multivariate linear [General LM] and non-linear [ANN] regression
library("caret")

f1 = "../results/regression/nuts2/datsel_nuts2"
f2 = "../results/regression/nuts2/datred_nuts2"
f3 = "../results/regression/nuts2/modsel_lm_nuts2"
f4 = "../results/regression/nuts2/fit_lm_nuts2"
f5 = "../results/regression/nuts2/modsel_nnet_nuts2"
f6 = "../results/regression/nuts2/fit_nnet_nuts2"
opt_nnet = TRUE
x_complete = remove_na(x_all_nuts2,init=2,plotit=FALSE,verbose=FALSE)
ina1 = !rownames(x_all_nuts2) %in% rownames(x_complete)
ina2 = !colnames(x_all_nuts2) %in% colnames(x_complete)
for(i1 in colnames(y_all_nuts2)){
    y = y_all_nuts2[,i1,drop=FALSE]
    if(is.null(levels(y[,1]))){
        f7 = paste(f1,"_",i1,".log",sep="")
        f8 = paste(f1,"_",i1,".csv",sep="")
        f9 = paste(f2,"_",i1,".log",sep="")
        capture.output(
          dim(x_complete),
          rownames(x_complete),
          colnames(x_complete),
          summary(x_complete),
          file=f7,append=FALSE,type="output")
        ina3 = is.na(y[,1]) 
        y = y[!(ina1 | ina3),,drop=FALSE]
        x = x_all_nuts2[!(ina1 | ina3),!ina2]
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        x = reduce_predictors_v2(x,y[,1],thr1=0.90,thr2=NULL,thr3=0.00,thr4=Inf,
          maxsize=30,f1=f9)
        capture.output(
          dim(x),
          rownames(x),
          colnames(x),
          summary(x),
          file=f7,append=TRUE,type="output")
        write.table(t(apply(x,2,summary)),file=f8,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f3,"_",i1,sep="")
        f8 = paste(f7,"_mod1.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,"_1.csv",sep="")
        f11 = paste(f7,"_2.csv",sep="")
        if(!file.exists(f9) | !file.exists(f10) | !file.exists(f11)){
            mod1 = mod_select_lm(x,y,thr1=1.00,thr2=0.00,thr3=0.00,thr4=Inf,
              maxsize=NULL,maxs1="hard",confs=100,level=1,thrs2=200000,f1=f7)
            save(mod1,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f4,"_",i1,sep="")
        f8 = paste(f7,".log",sep="")
        f9 = paste(f7,".csv",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_lm(mod1$formula,mod1$data,f1=f7,main="NUTS2")
        }
        f7 = paste(f1,"_",i1,"_scaled.csv",sep="")
        x_scaled = x
        s1 = match(colnames(x_scaled),colnames(x_all))
        iprop = xprop[s1]
        x_scaled[,iprop] = as.data.frame(apply(x_scaled[,iprop],2,
          function(x){x/100}))
        for(i2 in 1:ncol(x_scaled)){
            if(!iprop[i2]){
                xmin = x_all_min[s1][i2]
                xmax = x_all_max[s1][i2]
                x_scaled[,i2] = minmax(x_scaled[,i2],xmin,xmax)
            }
        }
        write.table(t(apply(x_scaled,2,summary)),file=f7,sep=",",dec=".",
          row.names=TRUE,col.names=NA)
        f7 = paste(f5,"_",i1,sep="")
        if(opt_nnet){
            f8 = paste(f7,"_optnnet.log",sep="")
            vars = all.vars(mod1$formula)
            dat1 = x_scaled
            colnames(dat1) = sapply(colnames(dat1),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            colnames(y) = sapply(colnames(y),function(x){gsub("[-|.]",
              "_",x)})                  #remove symbols used in "formulas"
            dat1 = as.data.frame(cbind(y,dat1[,vars[-1]]))
            grid1 = expand.grid(size=c(1,3*(1:10)),decay=1e-4)
            train1 = train(mod1$formula,dat1,method="nnet",abstol=1e-6,
              reltol=1e-12,linout=TRUE,trace=FALSE,tuneGrid=grid1,metric="RMSE",
              trControl=trainControl(method="boot",number=100))
            nn_size = train1$bestTune[1,1]
            nn_decay = train1$bestTune[1,2]
            capture.output(train1,file=f8,append=FALSE,type="output")
        }
        f8 = paste(f7,"_mod2.RDS",sep="")
        f9 = paste(f7,".log",sep="")
        f10 = paste(f7,".csv",sep="")
        if(!file.exists(f9) | !file.exists(f10)){
            mod2 = mod_select_nnet_v2(x_scaled,y,thr1=1.00,thr2=0.00,thr3=0.00,
              thr4=Inf,maxsize=NULL,maxs1="hard",nn_size=nn_size,
              nn_decay=nn_decay,nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,
              confs=100,level=1,thrs2=200000,f1=f7)
            save(mod2,file=f8)
        }else{
            load(file=f8)
        }
        f7 = paste(f6,"_",i1,sep="")
        f8 = paste(f7,"_1.log",sep="")
        f9 = paste(f7,"_2.log",sep="")
        if(!file.exists(f8) | !file.exists(f9)){
            fit_nnet_v2(mod2$formula,mod2$data,nn_size=NULL,nn_decay=1e-4,
              nn_maxit=1000,nn_abstol=1e-6,nn_reltol=1e-12,main="NUTS2",f1=f7)
        }
    }
}

}
