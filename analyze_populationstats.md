Collect info about population of study
========================================================


```r
setwd("Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/experimentsRadiologypaper-revision/Tree-based-RF/ensemble-Treebased-RF")

library("RSQLite")
```

```
## Loading required package: DBI
```

```r

read_data <- function(subdata, ids) {
    sqlite <- dbDriver("SQLite")
    conn <- dbConnect(sqlite, "stage1localData.db")
    
    # 2) all T1W features
    lesionsQuery <- dbGetQuery(conn, "SELECT *\n           FROM  stage1features\n           INNER JOIN lesion ON (stage1features.lesion_id = lesion.lesion_id)\n           INNER JOIN f_dynamic ON (stage1features.lesion_id = f_dynamic.lesion_id)\n           INNER JOIN f_morphology ON (stage1features.lesion_id = f_morphology.lesion_id)\n           INNER JOIN f_texture ON (stage1features.lesion_id = f_texture.lesion_id)")
    
    # prune entries and extract feature subsets corresponds to 5 entries
    # lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
    lesionfields = names(lesionsQuery)
    lesioninfo = lesionsQuery[c(1, 2, 150, 151)]
    stage1features = lesionsQuery[c(3:103, 124:127)]
    dynfeatures = lesionsQuery[c(154:187)]
    morphofeatures = lesionsQuery[c(190:208)]
    texfeatures = lesionsQuery[c(211:234)]
    
    # combine all features
    allfeatures = cbind(lesioninfo[c(2, 3)], stage1features, dynfeatures, morphofeatures, 
        texfeatures)
    
    if (subdata == "stage2") {
        # organized the data by subdata
        allfeatures = allfeatures[ids, ]
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "NC", "C")
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "NC", "C")
        allfeatures = rbind(M, N)
    }
    if (subdata == "stage1") {
        # organized the data by subdata
        allfeatures = allfeatures[ids, ]
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "mass", "mass")
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "nonmass", "nonmass")
        allfeatures = data.frame(rbind(M, N))
    }
    if (subdata == "oneshot") {
        # organized the data by subdata
        allfeatures = allfeatures[ids, ]
        M <- subset(allfeatures, lesion_label == "massB" | lesion_label == "massM")
        M$lesion_label <- ifelse(M$lesion_label == "massB", "NC", "C")
        N <- subset(allfeatures, lesion_label == "nonmassB" | lesion_label == 
            "nonmassM")
        N$lesion_label <- ifelse(N$lesion_label == "nonmassB", "NC", "C")
        allfeatures = data.frame(rbind(M, N))
    }
    # procees data
    allfeatures$lesion_label <- as.factor(allfeatures$lesion_label)
    allfeatures$peakCr_inside <- as.factor(allfeatures$peakCr_inside)
    allfeatures$peakVr_inside <- as.factor(allfeatures$peakVr_inside)
    allfeatures$peakCr_countor <- as.factor(allfeatures$peakCr_countor)
    allfeatures$peakVr_countor <- as.factor(allfeatures$peakVr_countor)
    allfeatures$k_Max_Margin_Grad <- as.factor(allfeatures$k_Max_Margin_Grad)
    allfeatures$max_RGH_mean_k <- as.factor(allfeatures$max_RGH_mean_k)
    allfeatures$max_RGH_var_k <- as.factor(allfeatures$max_RGH_var_k)
    
    # 2) all T1W features
    lesionsQueryinfo <- dbGetQuery(conn, "SELECT *\n           FROM  lesion\n           INNER JOIN stage1features ON (stage1features.lesion_id = lesion.lesion_id)")
    
    # prune entries and extract feature subsets corresponds to 5 entries
    # lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
    lesionsfields = names(lesionsQueryinfo)
    lesionsinfo = lesionsQueryinfo[c(1:24)]
    
    # 2) all T1W features
    lesionsQuerymass <- dbGetQuery(conn, "SELECT *\n           FROM  lesion\n           INNER JOIN mass_lesion ON (mass_lesion.lesion_id = lesion.lesion_id)\n           INNER JOIN stage1features ON (stage1features.lesion_id = lesion.lesion_id)")
    
    # prune entries and extract mass info
    lesionsmass = names(lesionsQuerymass)
    lesionsmassinfo = lesionsQuerymass[c(1:33)]
    
    # 2) all T1W features
    lesionsQuerynonmass <- dbGetQuery(conn, "SELECT *\n           FROM  lesion\n           INNER JOIN nonmass_lesion ON (nonmass_lesion.lesion_id = lesion.lesion_id)\n           INNER JOIN stage1features ON (stage1features.lesion_id = lesion.lesion_id)")
    
    # prune entries and extract nonmass
    lesionsmass = names(lesionsQuerynonmass)
    lesionsnonmassinfo = lesionsQuerynonmass[c(1:33)]
    
    output <- list(features = allfeatures, info = lesionsinfo, mass = lesionsmassinfo, 
        nonmass = lesionsnonmassinfo)
    return(output)
}
```


Start by loading data

```r
all = read_data(subdata = "multi", ids = 1:409)
lesioninfo = all$info
# install.packages('reshape')
library(reshape)
```

```
## Loading required package: plyr
## 
## Attaching package: 'reshape'
## 
## The following object(s) are masked from 'package:plyr':
## 
##     rename, round_any
```

```r

# Number of patients:
nopatients = cast(lesioninfo, ~cad_pt_no_txt)
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```r
print(length(nopatients))
```

```
## [1] 240
```

```r

# histopathology source
proc_source = cast(lesioninfo, ~proc_proc_source_int)
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```r
summary(proc_source)
```

```
##    value         NA       Radiology  
##  (all):1   Min.   :14   Min.   :332  
##            1st Qu.:14   1st Qu.:332  
##            Median :14   Median :332  
##            Mean   :14   Mean   :332  
##            3rd Qu.:14   3rd Qu.:332  
##            Max.   :14   Max.   :332  
##  Surgical/Operating Rm (includes 'Sentinel Lymph Node Biopsy')
##  Min.   :63                                                   
##  1st Qu.:63                                                   
##  Median :63                                                   
##  Mean   :63                                                   
##  3rd Qu.:63                                                   
##  Max.   :63
```

```r
# by guidance type
print(cast(lesioninfo, proc_proc_tp_int ~ proc_proc_source_int))
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##         proc_proc_tp_int NA Radiology
## 1     Core Needle Biopsy  0       221
## 2 Fine Needle Aspiration  0        17
## 3                    N/A  0         2
## 4                     NA 14         0
## 5 Vacuum Assisted Biopsy  0        92
##   Surgical/Operating Rm (includes 'Sentinel Lymph Node Biopsy')
## 1                                                             0
## 2                                                             0
## 3                                                            59
## 4                                                             4
## 5                                                             0
```

```r
# by guidance localization
print(cast(lesioninfo, proc_proc_tp_int ~ proc_proc_guid_int))
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##         proc_proc_tp_int MRI N/A NA None Stereo  US
## 1     Core Needle Biopsy  11  17  0    5     37 151
## 2 Fine Needle Aspiration   0   2  0    0      2  13
## 3                    N/A   0  61  0    0      0   0
## 4                     NA   0   0 18    0      0   0
## 5 Vacuum Assisted Biopsy  91   0  0    0      1   0
```

```r

# subset by mass and non-masses
print(cast(lesioninfo, ~lesion_label))
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB massM nonmassB nonmassM
## 1 (all)   136   144       65       64
```

```r

# subset by pathology
patho = summary(as.factor(lesioninfo$lesion_diagnosis))
for (k in 1:length(patho)) {
    print(cast(lesioninfo, ~lesion_label, subset = lesion_diagnosis == names(patho[k])))
}
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value          massM       nonmassM
## 1 (all) Adenocarcinoma Adenocarcinoma
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value    massB
## 1 (all) Adenosis
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value         massB
## 1 (all) AtypicalCells
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     5       11
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     5        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value             massB
## 1 (all) AtypicalPapilloma
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value              massB
## 1 (all) BenignbyAssumption
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value nonmassB
## 1 (all)        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB nonmassM
## 1 (all)    25       17        1
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value nonmassB
## 1 (all)        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     1        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value        nonmassB
## 1 (all) ColumnarChanges
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                  massB
## 1 (all) ComplexFibroepithelial
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value         nonmassB
## 1 (all) ComplexPapillary
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB
## 1 (all)     6
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value              nonmassB
## 1 (all) DifuseStromalFibrosis
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value             massB
## 1 (all) DuctalHyperplasia
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB
## 1 (all)     2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     1        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)    10        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                   massB
## 1 (all) DYSTROPHICCALCIFICATION
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)    31        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)    26       12
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value           massB
## 1 (all) Fibroepithelial
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     3        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value nonmassB
## 1 (all) Hematoma
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value       massB
## 1 (all) Hyperplasia
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massM nonmassM
## 1 (all)    30       26
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value         massB
## 1 (all) InsituLobular
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                    massM                 nonmassM
## 1 (all) INSITUPAPILLARYCARCINOMA INSITUPAPILLARYCARCINOMA
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massM nonmassM
## 1 (all)    95       29
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massM nonmassM
## 1 (all)    15        6
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                    massM
## 1 (all) InvasiveLobularCarcinoma
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB nonmassB
## 1 (all)     1        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                massM
## 1 (all) MetaplasticCarcinoma
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                  massB
## 1 (all) Papillary(focalAtypia)
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value          massB
## 1 (all) PhyllodesTumor
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value                    massB
## 1 (all) ProliferativeFibrocystic
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value      massB
## 1 (all) RadialScar
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value              massB
## 1 (all) SclerosingAdenosis
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value              massB
## 1 (all) SCLEROSINGADENOSIS
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value     massB
## 1 (all) Sclerosis
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
```

```
##   value             massB          nonmassB
## 1 (all) SclerosisAdenosis SclerosisAdenosis
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value nonmassB
## 1 (all)        2
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   value massB
## 1 (all)     2
```

```r
print(patho)
```

```
##             Adenocarcinoma                   Adenosis 
##                          2                          1 
##              AtypicalCells  AtypicalDuctalHyperplasia 
##                          1                         16 
## AtypicalLobularHyperplasia          AtypicalPapilloma 
##                          7                          1 
##         BenignbyAssumption           BenignbyFollowUp 
##                          1                          2 
##               BenignTissue ColumnarAlterationwoAtypia 
##                         43                          2 
##        ColumnarCellChanges            ColumnarChanges 
##                          3                          1 
##     ComplexFibroepithelial           ComplexPapillary 
##                          1                          1 
##                       Cyst      DifuseStromalFibrosis 
##                          6                          1 
##          DuctalHyperplasia   DuctalHyperplasiaWAtypia 
##                          1                          2 
##  DuctalHyperplasiaWoAtypia              DuctPapilloma 
##                          3                         12 
##    DYSTROPHICCALCIFICATION               Fibroadenoma 
##                          1                         33 
##                Fibrocystic            Fibroepithelial 
##                         38                          1 
##                   Fibrosis                   Hematoma 
##                          5                          1 
##                Hyperplasia               InsituDuctal 
##                          1                         56 
##              InsituLobular   INSITUPAPILLARYCARCINOMA 
##                          1                          2 
##             InvasiveDuctal            InvasiveLobular 
##                        124                         21 
##   InvasiveLobularCarcinoma         LobularHyperplasia 
##                          1                          3 
##       MetaplasticCarcinoma     Papillary(focalAtypia) 
##                          1                          1 
##             PhyllodesTumor   ProliferativeFibrocystic 
##                          1                          1 
##                 RadialScar         SclerosingAdenosis 
##                          1                          1 
##         SCLEROSINGADENOSIS                  Sclerosis 
##                          1                          1 
##          SclerosisAdenosis            StromalFibrosis 
##                          2                          2 
##     UsualDuctalHyperplasia 
##                          2
```

```r


# histopathology source
proc_source_type = cast(lesioninfo, proc_proc_guid_int ~ proc_proc_tp_int ~ 
    proc_proc_source_int)
```

```
## Using lesion_diagnosis as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```r
summary(proc_source)
```

```
##    value         NA       Radiology  
##  (all):1   Min.   :14   Min.   :332  
##            1st Qu.:14   1st Qu.:332  
##            Median :14   Median :332  
##            Mean   :14   Mean   :332  
##            3rd Qu.:14   3rd Qu.:332  
##            Max.   :14   Max.   :332  
##  Surgical/Operating Rm (includes 'Sentinel Lymph Node Biopsy')
##  Min.   :63                                                   
##  1st Qu.:63                                                   
##  Median :63                                                   
##  Mean   :63                                                   
##  3rd Qu.:63                                                   
##  Max.   :63
```

```r

```


Analyze BIRADS descriptors
=========

```r
# curve types
massinfo = all$mass
nonmassinfo = all$nonmass
print(cast(massinfo, lesion_label ~ find_curve_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label I Ia II III None Other
## 1        massB 0 16 27  20   57    16
## 2        massM 1  7 13  45   62    16
```

```r
print(cast(nonmassinfo, lesion_label ~ find_curve_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label Ia Ib II III None Other
## 1     nonmassB 13  0  8   1   33    10
## 2     nonmassM  1  2  7   2   44     8
```

```r

# initial enhancement
print(cast(massinfo, lesion_label ~ find_mri_dce_init_enh_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label M Moderate to marked N/A Rapid Slow Slow to medium
## 1        massB 0                 15  67    50    2              2
## 2        massM 1                  4  79    55    4              1
```

```r

# Morphology
print(cast(massinfo, lesion_label ~ find_mammo_n_mri_mass_shape_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label Irregular Lobular N/A Oval R Reinform Round
## 1        massB        10      22  75   17 0        1    11
## 2        massM        27      16  77    9 1        1    13
```

```r
print(cast(massinfo, lesion_label ~ find_mri_mass_margin_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label I Irregular N/A Smooth Spiculated
## 1        massB 0         5 115      9          7
## 2        massM 1        15  91      8         29
```

```r

# Nonmass
print(cast(nonmassinfo, lesion_label ~ find_mri_nonmass_int_enh_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label Clumped Heterogeneous Homogeneous N/A Stippled or punctate
## 1     nonmassB      13             6           1  45                    0
## 2     nonmassM      10             4           0  47                    3
```

```r
print(cast(nonmassinfo, lesion_label ~ find_mri_nonmass_dist_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label Ductal Focal Linear Multiple regions N/A Regional Segmental
## 1     nonmassB      1    15     24                0  22        1         2
## 2     nonmassM      3     9     19                1  17        3        12
```

```r

print(cast(nonmassinfo, lesion_label ~ find_mri_nonmass_dist_int + find_mri_nonmass_int_enh_int))
```

```
## Using lesion_id.2 as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##   lesion_label Ductal_Clumped Ductal_Heterogeneous Focal_Clumped
## 1     nonmassB              0                    1             1
## 2     nonmassM              2                    1             0
##   Focal_Heterogeneous Focal_N/A Linear_Clumped Linear_Heterogeneous
## 1                   3        11              6                    1
## 2                   2         7              6                    0
##   Linear_N/A Multiple regions_N/A N/A_Clumped N/A_Heterogeneous N/A_N/A
## 1         17                    0           5                 1      16
## 2         13                    1           0                 0      16
##   N/A_Stippled or punctate Regional_Clumped Regional_N/A Segmental_Clumped
## 1                        0                1            0                 0
## 2                        1                0            3                 2
##   Segmental_Heterogeneous Segmental_Homogeneous Segmental_N/A
## 1                       0                     1             1
## 2                       1                     0             7
##   Segmental_Stippled or punctate
## 1                              0
## 2                              2
```

```r

#
```


Query Biomatrix via ODBC library
====

```r
library(RODBC)
con <- odbcConnect("PostgreSQL35W")
sqlTables(con, "tbl_pt_exam")
```

```
##                TABLE_QUALIFIER TABLE_OWNER                      TABLE_NAME
## 1  biomatrixdb_raccess_mri_cad      public                      tbl_option
## 2  biomatrixdb_raccess_mri_cad      public               tbl_option_detail
## 3  biomatrixdb_raccess_mri_cad      public             tbl_pt_demographics
## 4  biomatrixdb_raccess_mri_cad      public                     tbl_pt_exam
## 5  biomatrixdb_raccess_mri_cad      public             tbl_pt_exam_finding
## 6  biomatrixdb_raccess_mri_cad      public tbl_pt_exam_finding_lesion_link
## 7  biomatrixdb_raccess_mri_cad      public              tbl_pt_exam_lesion
## 8  biomatrixdb_raccess_mri_cad      public         tbl_pt_exam_lesion_link
## 9  biomatrixdb_raccess_mri_cad      public           tbl_pt_mri_cad_can_tp
## 10 biomatrixdb_raccess_mri_cad      public           tbl_pt_mri_cad_record
## 11 biomatrixdb_raccess_mri_cad      public               tbl_pt_mri_series
## 12 biomatrixdb_raccess_mri_cad      public        tbl_pt_mri_series_4_test
## 13 biomatrixdb_raccess_mri_cad      public      tbl_pt_path_exam_find_link
## 14 biomatrixdb_raccess_mri_cad      public                tbl_pt_pathology
## 15 biomatrixdb_raccess_mri_cad      public    tbl_pt_pathology_lesion_link
## 16 biomatrixdb_raccess_mri_cad      public                tbl_pt_procedure
## 17 biomatrixdb_raccess_mri_cad      public                    tbl_usr_code
##    TABLE_TYPE REMARKS
## 1       TABLE        
## 2       TABLE        
## 3        VIEW        
## 4        VIEW        
## 5        VIEW        
## 6        VIEW        
## 7        VIEW        
## 8        VIEW        
## 9        VIEW        
## 10       VIEW        
## 11      TABLE        
## 12      TABLE        
## 13       VIEW        
## 14       VIEW        
## 15       VIEW        
## 16       VIEW        
## 17      TABLE
```

```r

lesioninfo_reasons = data.frame()
for (k in 1:length(lesioninfo$cad_pt_no_txt)) {
    print(k)
    CADid = lesioninfo$cad_pt_no_txt[k]
    datetime = lesioninfo$exam_dt_datetime[k]
    date = strsplit(datetime, " ")[[1]][1]
    query <- paste0("SELECT * FROM public.tbl_pt_mri_cad_record, public.tbl_pt_exam WHERE public.tbl_pt_exam.pt_id = public.tbl_pt_mri_cad_record.pt_id and cad_pt_no_txt = '", 
        CADid, "' and exam_dt_datetime = '", date, "'")
    res <- sqlQuery(con, query)
    lesioninfo_reasons = rbind(lesioninfo_reasons, res)
}
```

```
## [1] 1
## [1] 2
## [1] 3
## [1] 4
## [1] 5
## [1] 6
## [1] 7
## [1] 8
## [1] 9
## [1] 10
## [1] 11
## [1] 12
## [1] 13
## [1] 14
## [1] 15
## [1] 16
## [1] 17
## [1] 18
## [1] 19
## [1] 20
## [1] 21
## [1] 22
## [1] 23
## [1] 24
## [1] 25
## [1] 26
## [1] 27
## [1] 28
## [1] 29
## [1] 30
## [1] 31
## [1] 32
## [1] 33
## [1] 34
## [1] 35
## [1] 36
## [1] 37
## [1] 38
## [1] 39
## [1] 40
## [1] 41
## [1] 42
## [1] 43
## [1] 44
## [1] 45
## [1] 46
## [1] 47
## [1] 48
## [1] 49
## [1] 50
## [1] 51
## [1] 52
## [1] 53
## [1] 54
## [1] 55
## [1] 56
## [1] 57
## [1] 58
## [1] 59
## [1] 60
## [1] 61
## [1] 62
## [1] 63
## [1] 64
## [1] 65
## [1] 66
## [1] 67
## [1] 68
## [1] 69
## [1] 70
## [1] 71
## [1] 72
## [1] 73
## [1] 74
## [1] 75
## [1] 76
## [1] 77
## [1] 78
## [1] 79
## [1] 80
## [1] 81
## [1] 82
## [1] 83
## [1] 84
## [1] 85
## [1] 86
## [1] 87
## [1] 88
## [1] 89
## [1] 90
## [1] 91
## [1] 92
## [1] 93
## [1] 94
## [1] 95
## [1] 96
## [1] 97
## [1] 98
## [1] 99
## [1] 100
## [1] 101
## [1] 102
## [1] 103
## [1] 104
## [1] 105
## [1] 106
## [1] 107
## [1] 108
## [1] 109
## [1] 110
## [1] 111
## [1] 112
## [1] 113
## [1] 114
## [1] 115
## [1] 116
## [1] 117
## [1] 118
## [1] 119
## [1] 120
## [1] 121
## [1] 122
## [1] 123
## [1] 124
## [1] 125
## [1] 126
## [1] 127
## [1] 128
## [1] 129
## [1] 130
## [1] 131
## [1] 132
## [1] 133
## [1] 134
## [1] 135
## [1] 136
## [1] 137
## [1] 138
## [1] 139
## [1] 140
## [1] 141
## [1] 142
## [1] 143
## [1] 144
## [1] 145
## [1] 146
## [1] 147
## [1] 148
## [1] 149
## [1] 150
## [1] 151
## [1] 152
## [1] 153
## [1] 154
## [1] 155
## [1] 156
## [1] 157
## [1] 158
## [1] 159
## [1] 160
## [1] 161
## [1] 162
## [1] 163
## [1] 164
## [1] 165
## [1] 166
## [1] 167
## [1] 168
## [1] 169
## [1] 170
## [1] 171
## [1] 172
## [1] 173
## [1] 174
## [1] 175
## [1] 176
## [1] 177
## [1] 178
## [1] 179
## [1] 180
## [1] 181
## [1] 182
## [1] 183
## [1] 184
## [1] 185
## [1] 186
## [1] 187
## [1] 188
## [1] 189
## [1] 190
## [1] 191
## [1] 192
## [1] 193
## [1] 194
## [1] 195
## [1] 196
## [1] 197
## [1] 198
## [1] 199
## [1] 200
## [1] 201
## [1] 202
## [1] 203
## [1] 204
## [1] 205
## [1] 206
## [1] 207
## [1] 208
## [1] 209
## [1] 210
## [1] 211
## [1] 212
## [1] 213
## [1] 214
## [1] 215
## [1] 216
## [1] 217
## [1] 218
## [1] 219
## [1] 220
## [1] 221
## [1] 222
## [1] 223
## [1] 224
## [1] 225
## [1] 226
## [1] 227
## [1] 228
## [1] 229
## [1] 230
## [1] 231
## [1] 232
## [1] 233
## [1] 234
## [1] 235
## [1] 236
## [1] 237
## [1] 238
## [1] 239
## [1] 240
## [1] 241
## [1] 242
## [1] 243
## [1] 244
## [1] 245
## [1] 246
## [1] 247
## [1] 248
## [1] 249
## [1] 250
## [1] 251
## [1] 252
## [1] 253
## [1] 254
## [1] 255
## [1] 256
## [1] 257
## [1] 258
## [1] 259
## [1] 260
## [1] 261
## [1] 262
## [1] 263
## [1] 264
## [1] 265
## [1] 266
## [1] 267
## [1] 268
## [1] 269
## [1] 270
## [1] 271
## [1] 272
## [1] 273
## [1] 274
## [1] 275
## [1] 276
## [1] 277
## [1] 278
## [1] 279
## [1] 280
## [1] 281
## [1] 282
## [1] 283
## [1] 284
## [1] 285
## [1] 286
## [1] 287
## [1] 288
## [1] 289
## [1] 290
## [1] 291
## [1] 292
## [1] 293
## [1] 294
## [1] 295
## [1] 296
## [1] 297
## [1] 298
## [1] 299
## [1] 300
## [1] 301
## [1] 302
## [1] 303
## [1] 304
## [1] 305
## [1] 306
## [1] 307
## [1] 308
## [1] 309
## [1] 310
## [1] 311
## [1] 312
## [1] 313
## [1] 314
## [1] 315
## [1] 316
## [1] 317
## [1] 318
## [1] 319
## [1] 320
## [1] 321
## [1] 322
## [1] 323
## [1] 324
## [1] 325
## [1] 326
## [1] 327
## [1] 328
## [1] 329
## [1] 330
## [1] 331
## [1] 332
## [1] 333
## [1] 334
## [1] 335
## [1] 336
## [1] 337
## [1] 338
## [1] 339
## [1] 340
## [1] 341
## [1] 342
## [1] 343
## [1] 344
## [1] 345
## [1] 346
## [1] 347
## [1] 348
## [1] 349
## [1] 350
## [1] 351
## [1] 352
## [1] 353
## [1] 354
## [1] 355
## [1] 356
## [1] 357
## [1] 358
## [1] 359
## [1] 360
## [1] 361
## [1] 362
## [1] 363
## [1] 364
## [1] 365
## [1] 366
## [1] 367
## [1] 368
## [1] 369
## [1] 370
## [1] 371
## [1] 372
## [1] 373
## [1] 374
## [1] 375
## [1] 376
## [1] 377
## [1] 378
## [1] 379
## [1] 380
## [1] 381
## [1] 382
## [1] 383
## [1] 384
## [1] 385
## [1] 386
## [1] 387
## [1] 388
## [1] 389
## [1] 390
## [1] 391
## [1] 392
## [1] 393
## [1] 394
## [1] 395
## [1] 396
## [1] 397
## [1] 398
## [1] 399
## [1] 400
## [1] 401
## [1] 402
## [1] 403
## [1] 404
## [1] 405
## [1] 406
## [1] 407
## [1] 408
## [1] 409
```

```r

summary(lesioninfo_reasons)
```

```
##  pt_mri_cad_record_id     pt_id      cad_pt_no_txt 
##  Min.   :   9         Min.   :   3   Min.   :   2  
##  1st Qu.: 526         1st Qu.: 830   1st Qu.: 667  
##  Median : 606         Median :1367   Median : 811  
##  Mean   : 613         Mean   :1696   Mean   :2208  
##  3rd Qu.: 667         3rd Qu.:2278   3rd Qu.:6014  
##  Max.   :1242         Max.   :5842   Max.   :7151  
##                                                    
##  latest_mutation_status_int specific_mutation_txt menopause_datetime  
##  BRCA1    : 47              Length:428            Min.   :1998-01-01  
##  BRCA2    : 29              Class :character      1st Qu.:2000-12-15  
##  High Risk:130              Mode  :character      Median :2008-03-22  
##  Other    :183                                    Mean   :2005-07-09  
##  NA's     : 39                                    3rd Qu.:2010-01-01  
##                                                   Max.   :2010-01-01  
##                                                   NA's   :414         
##      menopause_tp_int   pt_exam_id      pt_id.1           exam_tp_int 
##  Oophorectomy:  6     Min.   : 298   Min.   :   3   MRI         :408  
##  NA's        :422     1st Qu.:3717   1st Qu.: 830   Mammographic: 20  
##                       Median :3837   Median :1367                     
##                       Mean   :4037   Mean   :1696                     
##                       3rd Qu.:3916   3rd Qu.:2278                     
##                       Max.   :8237   Max.   :5842                     
##                                                                       
##  exam_dt_datetime     previous_exam_when_int previous_exam_where_txt
##  Min.   :2001-08-10   Length:428             Length:428             
##  1st Qu.:2009-06-02   Class :character       Class :character       
##  Median :2010-03-01   Mode  :character       Mode  :character       
##  Mean   :2010-02-09                                                 
##  3rd Qu.:2011-02-22                                                 
##  Max.   :2014-02-01                                                 
##                                                                     
##  previous_exam_reason_int      side_int   exam_tp_mammo_int
##  Length:428               Bilateral:422   Mode:logical     
##  Class :character         Left     :  1   NA's:428         
##  Mode  :character         Right    :  5                    
##                                                            
##                                                            
##                                                            
##                                                            
##  sty_indicator_high_risk_yn sty_indicator_high_risk_brca_1_yn
##  Min.   :0.000              Min.   :0.000                    
##  1st Qu.:0.000              1st Qu.:0.000                    
##  Median :0.000              Median :0.000                    
##  Mean   :0.486              Mean   :0.035                    
##  3rd Qu.:1.000              3rd Qu.:0.000                    
##  Max.   :1.000              Max.   :1.000                    
##                                                              
##  sty_indicator_high_risk_brca_2_yn sty_indicator_high_risk_brca_1_or_2_yn
##  Min.   :0.0000                    Min.   :0.000                         
##  1st Qu.:0.0000                    1st Qu.:0.000                         
##  Median :0.0000                    Median :0.000                         
##  Mean   :0.0374                    Mean   :0.014                         
##  3rd Qu.:0.0000                    3rd Qu.:0.000                         
##  Max.   :1.0000                    Max.   :1.000                         
##                                                                          
##  sty_indicator_high_risk_at_yn sty_indicator_high_risk_other_gene_yn
##  Min.   :0                     Min.   :0.0000                       
##  1st Qu.:0                     1st Qu.:0.0000                       
##  Median :0                     Median :0.0000                       
##  Mean   :0                     Mean   :0.0023                       
##  3rd Qu.:0                     3rd Qu.:0.0000                       
##  Max.   :0                     Max.   :1.0000                       
##                                                                     
##  sty_indicator_high_risk_prior_high_risk_marker_yn
##  Min.   :0.0000                                   
##  1st Qu.:0.0000                                   
##  Median :0.0000                                   
##  Mean   :0.0748                                   
##  3rd Qu.:0.0000                                   
##  Max.   :1.0000                                   
##                                                   
##  sty_indicator_high_risk_prior_personal_can_hist_yn
##  Min.   :0.000                                     
##  1st Qu.:0.000                                     
##  Median :0.000                                     
##  Mean   :0.142                                     
##  3rd Qu.:0.000                                     
##  Max.   :1.000                                     
##                                                    
##  sty_indicator_high_risk_hist_of_mantle_rad_yn
##  Min.   :0                                    
##  1st Qu.:0                                    
##  Median :0                                    
##  Mean   :0                                    
##  3rd Qu.:0                                    
##  Max.   :0                                    
##                                               
##  sty_indicator_high_risk_fam_hist_yn
##  Min.   :0.000                      
##  1st Qu.:0.000                      
##  Median :0.000                      
##  Mean   :0.147                      
##  3rd Qu.:0.000                      
##  Max.   :1.000                      
##                                     
##  sty_indicator_pre_operative_extent_of_dis_yn
##  Min.   :0.000                               
##  1st Qu.:0.000                               
##  Median :0.000                               
##  Mean   :0.121                               
##  3rd Qu.:0.000                               
##  Max.   :1.000                               
##                                              
##  sty_indicator_post_operative_margin_yn sty_indicator_pre_neoadj_trtmnt_yn
##  Min.   :0                              Min.   :0.0000                    
##  1st Qu.:0                              1st Qu.:0.0000                    
##  Median :0                              Median :0.0000                    
##  Mean   :0                              Mean   :0.0047                    
##  3rd Qu.:0                              3rd Qu.:0.0000                    
##  Max.   :0                              Max.   :1.0000                    
##                                                                           
##  sty_indicator_post_neoadj_trtmnt_yn sty_indicator_prob_solv_diff_img_yn
##  Min.   :0                           Min.   :0.0000                     
##  1st Qu.:0                           1st Qu.:0.0000                     
##  Median :0                           Median :0.0000                     
##  Mean   :0                           Mean   :0.0397                     
##  3rd Qu.:0                           3rd Qu.:0.0000                     
##  Max.   :0                           Max.   :1.0000                     
##                                                                         
##  sty_indicator_scar_vs_recurr_yn sty_indicator_axil_meta_lymph_node_yn
##  Min.   :0                       Min.   :0                            
##  1st Qu.:0                       1st Qu.:0                            
##  Median :0                       Median :0                            
##  Mean   :0                       Mean   :0                            
##  3rd Qu.:0                       3rd Qu.:0                            
##  Max.   :0                       Max.   :0                            
##                                                                       
##  sty_indicator_non_axil_meta_lymph_node_yn
##  Min.   :0                                
##  1st Qu.:0                                
##  Median :0                                
##  Mean   :0                                
##  3rd Qu.:0                                
##  Max.   :0                                
##                                           
##  sty_indicator_folup_recommend_yn sty_indicator_mri_guided_bio_yn
##  Min.   :0.0000                   Min.   :0                      
##  1st Qu.:0.0000                   1st Qu.:0                      
##  Median :0.0000                   Median :0                      
##  Mean   :0.0397                   Mean   :0                      
##  3rd Qu.:0.0000                   3rd Qu.:0                      
##  Max.   :1.0000                   Max.   :0                      
##                                                                  
##  sty_indicator_folup_post_mri_guided_bio_yn
##  Min.   :0.000                             
##  1st Qu.:0.000                             
##  Median :0.000                             
##  Mean   :0.007                             
##  3rd Qu.:0.000                             
##  Max.   :1.000                             
##                                            
##  sty_indicator_prior_2_prophy_mast_yn a_number_txt      
##  Min.   :0.0000                       Length:428        
##  1st Qu.:0.0000                       Class :character  
##  Median :0.0000                       Mode  :character  
##  Mean   :0.0023                                         
##  3rd Qu.:0.0000                                         
##  Max.   :1.0000                                         
##                                                         
##  us_was_done_1_mn_prior_2_mri_yn us_was_done_1_mn_after_mri_yn
##  N/A:345                         N/A:345                      
##  No : 77                         No : 56                      
##  Yes:  6                         Yes: 27                      
##                                                               
##                                                               
##                                                               
##                                                               
##  currently_taking_tamoxifen_or_raloxifene_yn exam_tp_target_int
##  N/A:374                                     Breast:428        
##  No : 53                                                       
##  Yes:  1                                                       
##                                                                
##                                                                
##                                                                
##                                                                
##  exam_img_dicom_txt other_info_curr_menstrual_status_int
##  Length:428         Length:428                          
##  Class :character   Class :character                    
##  Mode  :character   Mode  :character                    
##                                                         
##                                                         
##                                                         
##                                                         
##  other_info_last_period_datetime sty_indicator_add_eval_as_folup_yn
##  Min.   :2007-01-01              Min.   :0.000                     
##  1st Qu.:2009-08-26              1st Qu.:0.000                     
##  Median :2010-04-01              Median :0.000                     
##  Mean   :2010-05-17              Mean   :0.264                     
##  3rd Qu.:2011-04-28              3rd Qu.:1.000                     
##  Max.   :2013-01-10              Max.   :1.000                     
##  NA's   :335                                                       
##  radiology_rpt_generated_yn
##  Min.   :0.000             
##  1st Qu.:1.000             
##  Median :1.000             
##  Mean   :0.953             
##  3rd Qu.:1.000             
##  Max.   :1.000             
##                            
##                                                                                                                                                                                                                                                                 original_report_txt 
##  BILATERAL BREAST MRI\r\n\r\nTECHNIQUE: Imaging was performed using a 1.5 T GE magnet with a dedicated breast coil. Localizer, sagittal T2 weighted FSE (TE: 88, TR:3000) fat saturated images, sagittal T1 weighted 3D FSPGR (TE:4.2, TR:8.2) bilateral images, sim          :  8  
##  ?<Report title=""""><Sessions>\r\nDate: Jun 17, 2010\r\n\r\nBILATERAL BREAST MRI\r\n\r\nCLINICAL INDICATION: Left breast Ca 1 o'clock. Outside MRI\r\nimages describe three separate lesions around mass and left lower\r\nouter left breast enhancement significance unknown:  5  
##  ?<Report title=""""><Sessions>\r\nDate: May 05, 2011\r\nBILATERAL BREAST MRI\r\n\r\nCLINICAL INDICATION:  31 yo female. Bilateral stable masses.  Papilloma with atypia in the right breast excised in 2009. Treated for adenoma and ADH in 2000.\r\n\r\nRecent mammogram s  :  5  
##  ?<Report title=""""><Sessions>\r\n&lt;Report title=""""&gt;&lt;Sessions&gt;\r\nBILATERAL BREAST MRI\r\n\r\nCLINICAL INDICATION: Bilateral or needle biopsies. The right\r\nbreast with a diagnosis of atypical ductal epithelial hyperplasia\r\nwith microcalcifications an  :  5  
##  ?<Report title=""""><Sessions>\r\nBILATERAL BREAST MRI\r\n\r\nCLINICAL HISTORY: High risk of screening study. BSO 2005. Left\r\nlumpectomy (DCIS) February 2006.\r\n\r\nTECHNIQUE: Imaging was performed with a dedicated breast coil on\r\na 1.5 T GE magnet. Localizer imag:  4  
##  (Other)                                                                                                                                                                                                                                                                      :381  
##  NA's                                                                                                                                                                                                                                                                         : 20  
##                                                                                                                                                                                                    comment_txt       
##  Left breast Ca 1 o'clock. Outside MRI\r\nimages describe three separate lesions around mass and left lower\r\nouter left breast enhancement significance unknown\r\n                                          :  5  
##  31 yo female. Bilateral stable masses.  Papilloma with atypia in the right breast excised in 2009. Treated for adenoma and ADH in 2000.                                                                       :  5  
##  Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.\r\nLeft breast severe atypical intraductal proliferative lesion.:  5  
##  BRCA 1 variant. Known multiple fibroadenomas.                                                                                                                                                                 :  5  
##  High risk of screening study. BSO 2005. Left lumpectomy (DCIS) February 2006.                                                                                                                                 :  4  
##  (Other)                                                                                                                                                                                                       :320  
##  NA's                                                                                                                                                                                                          : 84  
##  sty_indicator_folup_after_pre_exam_yn
##  Min.   :0.0000                       
##  1st Qu.:0.0000                       
##  Median :0.0000                       
##  Mean   :0.0537                       
##  3rd Qu.:0.0000                       
##  Max.   :1.0000                       
##                                       
##  sty_indicator_rout_screening_obsp_yn sty_indicator_mammo_guided_bio_yn
##  Min.   :0.0000                       Min.   :0                        
##  1st Qu.:0.0000                       1st Qu.:0                        
##  Median :0.0000                       Median :0                        
##  Mean   :0.0678                       Mean   :0                        
##  3rd Qu.:0.0000                       3rd Qu.:0                        
##  Max.   :1.0000                       Max.   :0                        
##                                                                        
##  sty_indicator_us_guided_bio_yn            mri_cad_status_txt
##  Min.   :0                      Malignant           :134     
##  1st Qu.:0                      Benign by assumption: 27     
##  Median :0                      Unknown             :173     
##  Mean   :0                      Benign by pathology : 74     
##  3rd Qu.:0                      NA's                : 20     
##  Max.   :0                                                   
## 
```

```r

# by guidance localization
print(cast(lesioninfo_reasons, latest_mutation_status_int + sty_indicator_high_risk_yn ~ 
    sty_indicator_high_risk_brca_1_yn))
```

```
## Using mri_cad_status_txt as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##    latest_mutation_status_int sty_indicator_high_risk_yn   0  1
## 1                       BRCA1                          0  10  0
## 2                       BRCA1                          1  22 15
## 3                       BRCA2                          0   8  0
## 4                       BRCA2                          1  21  0
## 5                   High Risk                          0  56  0
## 6                   High Risk                          1  74  0
## 7                       Other                          0 129  0
## 8                       Other                          1  54  0
## 9                        <NA>                          0  17  0
## 10                       <NA>                          1  22  0
```

```r

print(cast(lesioninfo_reasons, cad_pt_no_txt + sty_indicator_high_risk_yn + 
    sty_indicator_high_risk_brca_1_or_2_yn + sty_indicator_high_risk_other_gene_yn + 
    sty_indicator_high_risk_prior_high_risk_marker_yn + sty_indicator_high_risk_prior_personal_can_hist_yn + 
    sty_indicator_high_risk_hist_of_mantle_rad_yn + sty_indicator_high_risk_fam_hist_yn ~ 
    latest_mutation_status_int))
```

```
## Using mri_cad_status_txt as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##     cad_pt_no_txt sty_indicator_high_risk_yn
## 1               2                          1
## 2              16                          1
## 3              18                          1
## 4              25                          0
## 5              27                          0
## 6              59                          0
## 7              93                          1
## 8             103                          0
## 9             111                          0
## 10            114                          1
## 11            114                          1
## 12            121                          1
## 13            122                          0
## 14            123                          1
## 15            130                          1
## 16            130                          1
## 17            132                          0
## 18            133                          1
## 19            133                          1
## 20            166                          1
## 21            177                          1
## 22            180                          1
## 23            186                          0
## 24            186                          1
## 25            190                          0
## 26            196                          0
## 27            197                          1
## 28            198                          1
## 29            220                          1
## 30            229                          0
## 31            232                          0
## 32            246                          1
## 33            252                          0
## 34            259                          1
## 35            261                          1
## 36            271                          1
## 37            272                          1
## 38            276                          0
## 39            280                          1
## 40            282                          1
## 41            293                          1
## 42            299                          1
## 43            311                          0
## 44            331                          1
## 45            340                          1
## 46            409                          1
## 47            420                          1
## 48            426                          0
## 49            455                          1
## 50            456                          1
## 51            462                          1
## 52            465                          1
## 53            503                          1
## 54            513                          1
## 55            547                          1
## 56            553                          1
## 57            556                          1
## 58            559                          1
## 59            576                          1
## 60            578                          1
## 61            580                          1
## 62            606                          1
## 63            635                          0
## 64            651                          0
## 65            657                          0
## 66            663                          1
## 67            664                          1
## 68            666                          1
## 69            667                          1
## 70            667                          1
## 71            668                          0
## 72            672                          0
## 73            673                          1
## 74            679                          1
## 75            681                          0
## 76            683                          0
## 77            684                          1
## 78            685                          1
## 79            687                          0
## 80            689                          1
## 81            690                          1
## 82            691                          0
## 83            700                          0
## 84            705                          0
## 85            706                          1
## 86            707                          0
## 87            710                          0
## 88            713                          0
## 89            714                          0
## 90            718                          0
## 91            721                          0
## 92            722                          0
## 93            724                          0
## 94            726                          1
## 95            727                          1
## 96            728                          1
## 97            729                          1
## 98            730                          0
## 99            731                          0
## 100           735                          0
## 101           736                          0
## 102           742                          0
## 103           744                          0
## 104           745                          0
## 105           747                          0
## 106           752                          0
## 107           755                          0
## 108           757                          1
## 109           758                          1
## 110           760                          0
## 111           764                          0
## 112           765                          0
## 113           771                          0
## 114           775                          0
## 115           775                          1
## 116           776                          0
## 117           778                          0
## 118           779                          0
## 119           781                          0
## 120           782                          1
## 121           783                          0
## 122           789                          1
## 123           790                          1
## 124           791                          1
## 125           792                          1
## 126           793                          0
## 127           793                          1
## 128           795                          1
## 129           796                          1
## 130           799                          0
## 131           802                          0
## 132           803                          0
## 133           805                          1
## 134           807                          0
## 135           809                          0
## 136           810                          0
## 137           812                          1
## 138           813                          0
## 139           814                          0
## 140           815                          0
## 141           817                          0
## 142           818                          0
## 143           827                          0
## 144           829                          0
## 145           830                          0
## 146           831                          0
## 147           837                          0
## 148           839                          1
## 149           843                          0
## 150           843                          1
## 151           845                          0
## 152           846                          0
## 153           847                          0
## 154           850                          0
## 155           851                          0
## 156           852                          0
## 157           853                          0
## 158           853                          1
## 159           855                          1
## 160           856                          1
## 161           857                          0
## 162           860                          0
## 163           861                          0
## 164           862                          1
## 165           863                          1
## 166           865                          0
## 167           867                          0
## 168           870                          0
## 169           870                          1
## 170           871                          1
## 171           873                          1
## 172           875                          0
## 173           877                          1
## 174           880                          1
## 175           881                          0
## 176           883                          0
## 177           884                          1
## 178           887                          1
## 179           888                          0
## 180           896                          0
## 181           900                          0
## 182           904                          1
## 183           913                          0
## 184           918                          1
## 185           920                          0
## 186           921                          0
## 187           937                          1
## 188           950                          0
## 189           956                          1
## 190          6001                          1
## 191          6004                          1
## 192          6005                          1
## 193          6008                          1
## 194          6014                          1
## 195          6015                          0
## 196          6017                          1
## 197          6018                          0
## 198          6019                          1
## 199          6020                          1
## 200          6021                          0
## 201          6022                          1
## 202          6023                          1
## 203          6024                          0
## 204          6025                          1
## 205          6026                          1
## 206          6027                          1
## 207          6029                          0
## 208          6029                          1
## 209          6032                          0
## 210          6034                          1
## 211          6035                          0
## 212          6036                          1
## 213          6037                          1
## 214          6038                          1
## 215          6039                          0
## 216          6040                          0
## 217          6041                          0
## 218          6041                          1
## 219          6042                          1
## 220          6043                          1
## 221          6044                          0
## 222          6045                          0
## 223          6046                          0
## 224          6047                          0
## 225          6048                          0
## 226          6050                          1
## 227          6051                          1
## 228          6052                          1
## 229          6054                          0
## 230          7008                          0
## 231          7011                          1
## 232          7018                          0
## 233          7024                          0
## 234          7029                          0
## 235          7030                          0
## 236          7043                          1
## 237          7045                          1
## 238          7054                          1
## 239          7066                          1
## 240          7076                          1
## 241          7077                          1
## 242          7085                          1
## 243          7086                          1
## 244          7088                          0
## 245          7094                          1
## 246          7096                          0
## 247          7097                          1
## 248          7104                          0
## 249          7110                          0
## 250          7127                          1
## 251          7151                          1
##     sty_indicator_high_risk_brca_1_or_2_yn
## 1                                        0
## 2                                        0
## 3                                        0
## 4                                        0
## 5                                        0
## 6                                        0
## 7                                        0
## 8                                        0
## 9                                        0
## 10                                       0
## 11                                       1
## 12                                       0
## 13                                       0
## 14                                       0
## 15                                       0
## 16                                       0
## 17                                       0
## 18                                       0
## 19                                       0
## 20                                       0
## 21                                       0
## 22                                       0
## 23                                       0
## 24                                       0
## 25                                       0
## 26                                       0
## 27                                       0
## 28                                       0
## 29                                       0
## 30                                       0
## 31                                       0
## 32                                       0
## 33                                       0
## 34                                       0
## 35                                       0
## 36                                       0
## 37                                       0
## 38                                       0
## 39                                       0
## 40                                       0
## 41                                       0
## 42                                       0
## 43                                       0
## 44                                       0
## 45                                       1
## 46                                       0
## 47                                       0
## 48                                       0
## 49                                       0
## 50                                       0
## 51                                       0
## 52                                       0
## 53                                       0
## 54                                       0
## 55                                       0
## 56                                       0
## 57                                       0
## 58                                       0
## 59                                       0
## 60                                       0
## 61                                       0
## 62                                       0
## 63                                       0
## 64                                       0
## 65                                       0
## 66                                       0
## 67                                       0
## 68                                       0
## 69                                       0
## 70                                       0
## 71                                       0
## 72                                       0
## 73                                       0
## 74                                       0
## 75                                       0
## 76                                       0
## 77                                       0
## 78                                       0
## 79                                       0
## 80                                       0
## 81                                       0
## 82                                       0
## 83                                       0
## 84                                       0
## 85                                       0
## 86                                       0
## 87                                       0
## 88                                       0
## 89                                       0
## 90                                       0
## 91                                       0
## 92                                       0
## 93                                       0
## 94                                       0
## 95                                       0
## 96                                       0
## 97                                       0
## 98                                       0
## 99                                       0
## 100                                      0
## 101                                      0
## 102                                      0
## 103                                      0
## 104                                      0
## 105                                      0
## 106                                      0
## 107                                      0
## 108                                      0
## 109                                      0
## 110                                      0
## 111                                      0
## 112                                      0
## 113                                      0
## 114                                      0
## 115                                      0
## 116                                      0
## 117                                      0
## 118                                      0
## 119                                      0
## 120                                      0
## 121                                      0
## 122                                      0
## 123                                      0
## 124                                      0
## 125                                      0
## 126                                      0
## 127                                      0
## 128                                      0
## 129                                      0
## 130                                      0
## 131                                      0
## 132                                      0
## 133                                      0
## 134                                      0
## 135                                      0
## 136                                      0
## 137                                      0
## 138                                      0
## 139                                      0
## 140                                      0
## 141                                      0
## 142                                      0
## 143                                      0
## 144                                      0
## 145                                      0
## 146                                      0
## 147                                      0
## 148                                      0
## 149                                      0
## 150                                      0
## 151                                      0
## 152                                      0
## 153                                      0
## 154                                      0
## 155                                      0
## 156                                      0
## 157                                      0
## 158                                      0
## 159                                      0
## 160                                      0
## 161                                      0
## 162                                      0
## 163                                      0
## 164                                      0
## 165                                      0
## 166                                      0
## 167                                      0
## 168                                      0
## 169                                      0
## 170                                      0
## 171                                      0
## 172                                      0
## 173                                      0
## 174                                      0
## 175                                      0
## 176                                      0
## 177                                      0
## 178                                      0
## 179                                      0
## 180                                      0
## 181                                      0
## 182                                      0
## 183                                      0
## 184                                      0
## 185                                      0
## 186                                      0
## 187                                      0
## 188                                      0
## 189                                      0
## 190                                      0
## 191                                      0
## 192                                      0
## 193                                      0
## 194                                      0
## 195                                      0
## 196                                      0
## 197                                      0
## 198                                      0
## 199                                      0
## 200                                      0
## 201                                      0
## 202                                      0
## 203                                      0
## 204                                      0
## 205                                      0
## 206                                      0
## 207                                      0
## 208                                      0
## 209                                      0
## 210                                      0
## 211                                      0
## 212                                      0
## 213                                      0
## 214                                      0
## 215                                      0
## 216                                      0
## 217                                      0
## 218                                      0
## 219                                      0
## 220                                      0
## 221                                      0
## 222                                      0
## 223                                      0
## 224                                      0
## 225                                      0
## 226                                      1
## 227                                      0
## 228                                      0
## 229                                      0
## 230                                      0
## 231                                      0
## 232                                      0
## 233                                      0
## 234                                      0
## 235                                      0
## 236                                      0
## 237                                      0
## 238                                      0
## 239                                      0
## 240                                      0
## 241                                      1
## 242                                      0
## 243                                      0
## 244                                      0
## 245                                      0
## 246                                      0
## 247                                      0
## 248                                      0
## 249                                      0
## 250                                      0
## 251                                      0
##     sty_indicator_high_risk_other_gene_yn
## 1                                       0
## 2                                       0
## 3                                       0
## 4                                       0
## 5                                       0
## 6                                       0
## 7                                       0
## 8                                       0
## 9                                       0
## 10                                      0
## 11                                      0
## 12                                      0
## 13                                      0
## 14                                      0
## 15                                      0
## 16                                      0
## 17                                      0
## 18                                      0
## 19                                      0
## 20                                      1
## 21                                      0
## 22                                      0
## 23                                      0
## 24                                      0
## 25                                      0
## 26                                      0
## 27                                      0
## 28                                      0
## 29                                      0
## 30                                      0
## 31                                      0
## 32                                      0
## 33                                      0
## 34                                      0
## 35                                      0
## 36                                      0
## 37                                      0
## 38                                      0
## 39                                      0
## 40                                      0
## 41                                      0
## 42                                      0
## 43                                      0
## 44                                      0
## 45                                      0
## 46                                      0
## 47                                      0
## 48                                      0
## 49                                      0
## 50                                      0
## 51                                      0
## 52                                      0
## 53                                      0
## 54                                      0
## 55                                      0
## 56                                      0
## 57                                      0
## 58                                      0
## 59                                      0
## 60                                      0
## 61                                      0
## 62                                      0
## 63                                      0
## 64                                      0
## 65                                      0
## 66                                      0
## 67                                      0
## 68                                      0
## 69                                      0
## 70                                      0
## 71                                      0
## 72                                      0
## 73                                      0
## 74                                      0
## 75                                      0
## 76                                      0
## 77                                      0
## 78                                      0
## 79                                      0
## 80                                      0
## 81                                      0
## 82                                      0
## 83                                      0
## 84                                      0
## 85                                      0
## 86                                      0
## 87                                      0
## 88                                      0
## 89                                      0
## 90                                      0
## 91                                      0
## 92                                      0
## 93                                      0
## 94                                      0
## 95                                      0
## 96                                      0
## 97                                      0
## 98                                      0
## 99                                      0
## 100                                     0
## 101                                     0
## 102                                     0
## 103                                     0
## 104                                     0
## 105                                     0
## 106                                     0
## 107                                     0
## 108                                     0
## 109                                     0
## 110                                     0
## 111                                     0
## 112                                     0
## 113                                     0
## 114                                     0
## 115                                     0
## 116                                     0
## 117                                     0
## 118                                     0
## 119                                     0
## 120                                     0
## 121                                     0
## 122                                     0
## 123                                     0
## 124                                     0
## 125                                     0
## 126                                     0
## 127                                     0
## 128                                     0
## 129                                     0
## 130                                     0
## 131                                     0
## 132                                     0
## 133                                     0
## 134                                     0
## 135                                     0
## 136                                     0
## 137                                     0
## 138                                     0
## 139                                     0
## 140                                     0
## 141                                     0
## 142                                     0
## 143                                     0
## 144                                     0
## 145                                     0
## 146                                     0
## 147                                     0
## 148                                     0
## 149                                     0
## 150                                     0
## 151                                     0
## 152                                     0
## 153                                     0
## 154                                     0
## 155                                     0
## 156                                     0
## 157                                     0
## 158                                     0
## 159                                     0
## 160                                     0
## 161                                     0
## 162                                     0
## 163                                     0
## 164                                     0
## 165                                     0
## 166                                     0
## 167                                     0
## 168                                     0
## 169                                     0
## 170                                     0
## 171                                     0
## 172                                     0
## 173                                     0
## 174                                     0
## 175                                     0
## 176                                     0
## 177                                     0
## 178                                     0
## 179                                     0
## 180                                     0
## 181                                     0
## 182                                     0
## 183                                     0
## 184                                     0
## 185                                     0
## 186                                     0
## 187                                     0
## 188                                     0
## 189                                     0
## 190                                     0
## 191                                     0
## 192                                     0
## 193                                     0
## 194                                     0
## 195                                     0
## 196                                     0
## 197                                     0
## 198                                     0
## 199                                     0
## 200                                     0
## 201                                     0
## 202                                     0
## 203                                     0
## 204                                     0
## 205                                     0
## 206                                     0
## 207                                     0
## 208                                     0
## 209                                     0
## 210                                     0
## 211                                     0
## 212                                     0
## 213                                     0
## 214                                     0
## 215                                     0
## 216                                     0
## 217                                     0
## 218                                     0
## 219                                     0
## 220                                     0
## 221                                     0
## 222                                     0
## 223                                     0
## 224                                     0
## 225                                     0
## 226                                     0
## 227                                     0
## 228                                     0
## 229                                     0
## 230                                     0
## 231                                     0
## 232                                     0
## 233                                     0
## 234                                     0
## 235                                     0
## 236                                     0
## 237                                     0
## 238                                     0
## 239                                     0
## 240                                     0
## 241                                     0
## 242                                     0
## 243                                     0
## 244                                     0
## 245                                     0
## 246                                     0
## 247                                     0
## 248                                     0
## 249                                     0
## 250                                     0
## 251                                     0
##     sty_indicator_high_risk_prior_high_risk_marker_yn
## 1                                                   0
## 2                                                   0
## 3                                                   0
## 4                                                   0
## 5                                                   0
## 6                                                   0
## 7                                                   0
## 8                                                   0
## 9                                                   0
## 10                                                  0
## 11                                                  0
## 12                                                  0
## 13                                                  0
## 14                                                  0
## 15                                                  0
## 16                                                  0
## 17                                                  0
## 18                                                  0
## 19                                                  0
## 20                                                  0
## 21                                                  1
## 22                                                  0
## 23                                                  0
## 24                                                  1
## 25                                                  0
## 26                                                  0
## 27                                                  0
## 28                                                  0
## 29                                                  0
## 30                                                  0
## 31                                                  0
## 32                                                  0
## 33                                                  0
## 34                                                  0
## 35                                                  0
## 36                                                  0
## 37                                                  0
## 38                                                  0
## 39                                                  0
## 40                                                  0
## 41                                                  0
## 42                                                  0
## 43                                                  0
## 44                                                  0
## 45                                                  0
## 46                                                  0
## 47                                                  1
## 48                                                  0
## 49                                                  0
## 50                                                  0
## 51                                                  0
## 52                                                  1
## 53                                                  0
## 54                                                  0
## 55                                                  0
## 56                                                  0
## 57                                                  0
## 58                                                  0
## 59                                                  0
## 60                                                  0
## 61                                                  0
## 62                                                  0
## 63                                                  0
## 64                                                  0
## 65                                                  0
## 66                                                  0
## 67                                                  0
## 68                                                  0
## 69                                                  0
## 70                                                  0
## 71                                                  0
## 72                                                  0
## 73                                                  0
## 74                                                  0
## 75                                                  0
## 76                                                  0
## 77                                                  1
## 78                                                  0
## 79                                                  0
## 80                                                  1
## 81                                                  0
## 82                                                  0
## 83                                                  0
## 84                                                  0
## 85                                                  1
## 86                                                  0
## 87                                                  0
## 88                                                  0
## 89                                                  0
## 90                                                  0
## 91                                                  0
## 92                                                  0
## 93                                                  0
## 94                                                  0
## 95                                                  0
## 96                                                  0
## 97                                                  1
## 98                                                  0
## 99                                                  0
## 100                                                 0
## 101                                                 0
## 102                                                 0
## 103                                                 0
## 104                                                 0
## 105                                                 0
## 106                                                 0
## 107                                                 0
## 108                                                 1
## 109                                                 1
## 110                                                 0
## 111                                                 0
## 112                                                 0
## 113                                                 0
## 114                                                 0
## 115                                                 0
## 116                                                 0
## 117                                                 0
## 118                                                 0
## 119                                                 0
## 120                                                 0
## 121                                                 0
## 122                                                 0
## 123                                                 0
## 124                                                 1
## 125                                                 0
## 126                                                 0
## 127                                                 0
## 128                                                 0
## 129                                                 1
## 130                                                 0
## 131                                                 0
## 132                                                 0
## 133                                                 0
## 134                                                 0
## 135                                                 0
## 136                                                 0
## 137                                                 0
## 138                                                 0
## 139                                                 0
## 140                                                 0
## 141                                                 0
## 142                                                 0
## 143                                                 0
## 144                                                 0
## 145                                                 0
## 146                                                 0
## 147                                                 0
## 148                                                 0
## 149                                                 0
## 150                                                 0
## 151                                                 0
## 152                                                 0
## 153                                                 0
## 154                                                 0
## 155                                                 0
## 156                                                 0
## 157                                                 0
## 158                                                 0
## 159                                                 1
## 160                                                 0
## 161                                                 0
## 162                                                 0
## 163                                                 0
## 164                                                 1
## 165                                                 0
## 166                                                 0
## 167                                                 0
## 168                                                 0
## 169                                                 0
## 170                                                 1
## 171                                                 1
## 172                                                 0
## 173                                                 0
## 174                                                 1
## 175                                                 0
## 176                                                 0
## 177                                                 0
## 178                                                 0
## 179                                                 0
## 180                                                 0
## 181                                                 0
## 182                                                 0
## 183                                                 0
## 184                                                 0
## 185                                                 0
## 186                                                 0
## 187                                                 0
## 188                                                 0
## 189                                                 0
## 190                                                 0
## 191                                                 0
## 192                                                 0
## 193                                                 0
## 194                                                 0
## 195                                                 0
## 196                                                 0
## 197                                                 0
## 198                                                 1
## 199                                                 0
## 200                                                 0
## 201                                                 0
## 202                                                 0
## 203                                                 0
## 204                                                 0
## 205                                                 0
## 206                                                 0
## 207                                                 0
## 208                                                 0
## 209                                                 0
## 210                                                 0
## 211                                                 0
## 212                                                 0
## 213                                                 0
## 214                                                 0
## 215                                                 0
## 216                                                 0
## 217                                                 0
## 218                                                 0
## 219                                                 0
## 220                                                 0
## 221                                                 0
## 222                                                 0
## 223                                                 0
## 224                                                 0
## 225                                                 0
## 226                                                 0
## 227                                                 0
## 228                                                 0
## 229                                                 0
## 230                                                 0
## 231                                                 0
## 232                                                 0
## 233                                                 0
## 234                                                 0
## 235                                                 0
## 236                                                 0
## 237                                                 0
## 238                                                 0
## 239                                                 0
## 240                                                 0
## 241                                                 0
## 242                                                 0
## 243                                                 0
## 244                                                 0
## 245                                                 0
## 246                                                 0
## 247                                                 0
## 248                                                 0
## 249                                                 0
## 250                                                 0
## 251                                                 0
##     sty_indicator_high_risk_prior_personal_can_hist_yn
## 1                                                    0
## 2                                                    0
## 3                                                    1
## 4                                                    0
## 5                                                    0
## 6                                                    0
## 7                                                    0
## 8                                                    0
## 9                                                    0
## 10                                                   0
## 11                                                   0
## 12                                                   0
## 13                                                   0
## 14                                                   0
## 15                                                   0
## 16                                                   0
## 17                                                   0
## 18                                                   0
## 19                                                   0
## 20                                                   0
## 21                                                   0
## 22                                                   0
## 23                                                   0
## 24                                                   0
## 25                                                   0
## 26                                                   0
## 27                                                   0
## 28                                                   0
## 29                                                   0
## 30                                                   0
## 31                                                   0
## 32                                                   0
## 33                                                   0
## 34                                                   0
## 35                                                   0
## 36                                                   0
## 37                                                   0
## 38                                                   0
## 39                                                   0
## 40                                                   0
## 41                                                   0
## 42                                                   1
## 43                                                   0
## 44                                                   0
## 45                                                   0
## 46                                                   1
## 47                                                   0
## 48                                                   0
## 49                                                   0
## 50                                                   0
## 51                                                   0
## 52                                                   0
## 53                                                   0
## 54                                                   0
## 55                                                   0
## 56                                                   0
## 57                                                   0
## 58                                                   0
## 59                                                   0
## 60                                                   0
## 61                                                   0
## 62                                                   0
## 63                                                   0
## 64                                                   0
## 65                                                   0
## 66                                                   1
## 67                                                   0
## 68                                                   0
## 69                                                   0
## 70                                                   0
## 71                                                   0
## 72                                                   0
## 73                                                   0
## 74                                                   1
## 75                                                   0
## 76                                                   0
## 77                                                   0
## 78                                                   0
## 79                                                   0
## 80                                                   0
## 81                                                   0
## 82                                                   0
## 83                                                   0
## 84                                                   0
## 85                                                   0
## 86                                                   0
## 87                                                   0
## 88                                                   0
## 89                                                   0
## 90                                                   0
## 91                                                   0
## 92                                                   0
## 93                                                   0
## 94                                                   1
## 95                                                   0
## 96                                                   1
## 97                                                   0
## 98                                                   0
## 99                                                   0
## 100                                                  0
## 101                                                  0
## 102                                                  0
## 103                                                  0
## 104                                                  0
## 105                                                  0
## 106                                                  0
## 107                                                  0
## 108                                                  0
## 109                                                  0
## 110                                                  0
## 111                                                  0
## 112                                                  0
## 113                                                  0
## 114                                                  0
## 115                                                  0
## 116                                                  0
## 117                                                  0
## 118                                                  0
## 119                                                  0
## 120                                                  0
## 121                                                  0
## 122                                                  0
## 123                                                  0
## 124                                                  0
## 125                                                  0
## 126                                                  0
## 127                                                  0
## 128                                                  0
## 129                                                  0
## 130                                                  0
## 131                                                  0
## 132                                                  0
## 133                                                  0
## 134                                                  0
## 135                                                  0
## 136                                                  0
## 137                                                  1
## 138                                                  0
## 139                                                  0
## 140                                                  0
## 141                                                  0
## 142                                                  0
## 143                                                  0
## 144                                                  0
## 145                                                  0
## 146                                                  0
## 147                                                  0
## 148                                                  0
## 149                                                  0
## 150                                                  0
## 151                                                  0
## 152                                                  0
## 153                                                  0
## 154                                                  0
## 155                                                  0
## 156                                                  0
## 157                                                  0
## 158                                                  0
## 159                                                  0
## 160                                                  0
## 161                                                  0
## 162                                                  0
## 163                                                  0
## 164                                                  0
## 165                                                  0
## 166                                                  0
## 167                                                  0
## 168                                                  0
## 169                                                  1
## 170                                                  0
## 171                                                  0
## 172                                                  0
## 173                                                  0
## 174                                                  0
## 175                                                  0
## 176                                                  0
## 177                                                  1
## 178                                                  0
## 179                                                  0
## 180                                                  0
## 181                                                  0
## 182                                                  0
## 183                                                  0
## 184                                                  0
## 185                                                  0
## 186                                                  0
## 187                                                  0
## 188                                                  0
## 189                                                  0
## 190                                                  1
## 191                                                  1
## 192                                                  1
## 193                                                  1
## 194                                                  1
## 195                                                  0
## 196                                                  1
## 197                                                  0
## 198                                                  0
## 199                                                  1
## 200                                                  0
## 201                                                  1
## 202                                                  0
## 203                                                  0
## 204                                                  1
## 205                                                  1
## 206                                                  1
## 207                                                  0
## 208                                                  1
## 209                                                  0
## 210                                                  1
## 211                                                  0
## 212                                                  1
## 213                                                  1
## 214                                                  1
## 215                                                  0
## 216                                                  0
## 217                                                  0
## 218                                                  1
## 219                                                  0
## 220                                                  1
## 221                                                  0
## 222                                                  0
## 223                                                  0
## 224                                                  0
## 225                                                  0
## 226                                                  1
## 227                                                  1
## 228                                                  1
## 229                                                  0
## 230                                                  0
## 231                                                  0
## 232                                                  0
## 233                                                  0
## 234                                                  0
## 235                                                  0
## 236                                                  0
## 237                                                  0
## 238                                                  0
## 239                                                  0
## 240                                                  0
## 241                                                  1
## 242                                                  0
## 243                                                  0
## 244                                                  0
## 245                                                  0
## 246                                                  0
## 247                                                  0
## 248                                                  0
## 249                                                  0
## 250                                                  0
## 251                                                  0
##     sty_indicator_high_risk_hist_of_mantle_rad_yn
## 1                                               0
## 2                                               0
## 3                                               0
## 4                                               0
## 5                                               0
## 6                                               0
## 7                                               0
## 8                                               0
## 9                                               0
## 10                                              0
## 11                                              0
## 12                                              0
## 13                                              0
## 14                                              0
## 15                                              0
## 16                                              0
## 17                                              0
## 18                                              0
## 19                                              0
## 20                                              0
## 21                                              0
## 22                                              0
## 23                                              0
## 24                                              0
## 25                                              0
## 26                                              0
## 27                                              0
## 28                                              0
## 29                                              0
## 30                                              0
## 31                                              0
## 32                                              0
## 33                                              0
## 34                                              0
## 35                                              0
## 36                                              0
## 37                                              0
## 38                                              0
## 39                                              0
## 40                                              0
## 41                                              0
## 42                                              0
## 43                                              0
## 44                                              0
## 45                                              0
## 46                                              0
## 47                                              0
## 48                                              0
## 49                                              0
## 50                                              0
## 51                                              0
## 52                                              0
## 53                                              0
## 54                                              0
## 55                                              0
## 56                                              0
## 57                                              0
## 58                                              0
## 59                                              0
## 60                                              0
## 61                                              0
## 62                                              0
## 63                                              0
## 64                                              0
## 65                                              0
## 66                                              0
## 67                                              0
## 68                                              0
## 69                                              0
## 70                                              0
## 71                                              0
## 72                                              0
## 73                                              0
## 74                                              0
## 75                                              0
## 76                                              0
## 77                                              0
## 78                                              0
## 79                                              0
## 80                                              0
## 81                                              0
## 82                                              0
## 83                                              0
## 84                                              0
## 85                                              0
## 86                                              0
## 87                                              0
## 88                                              0
## 89                                              0
## 90                                              0
## 91                                              0
## 92                                              0
## 93                                              0
## 94                                              0
## 95                                              0
## 96                                              0
## 97                                              0
## 98                                              0
## 99                                              0
## 100                                             0
## 101                                             0
## 102                                             0
## 103                                             0
## 104                                             0
## 105                                             0
## 106                                             0
## 107                                             0
## 108                                             0
## 109                                             0
## 110                                             0
## 111                                             0
## 112                                             0
## 113                                             0
## 114                                             0
## 115                                             0
## 116                                             0
## 117                                             0
## 118                                             0
## 119                                             0
## 120                                             0
## 121                                             0
## 122                                             0
## 123                                             0
## 124                                             0
## 125                                             0
## 126                                             0
## 127                                             0
## 128                                             0
## 129                                             0
## 130                                             0
## 131                                             0
## 132                                             0
## 133                                             0
## 134                                             0
## 135                                             0
## 136                                             0
## 137                                             0
## 138                                             0
## 139                                             0
## 140                                             0
## 141                                             0
## 142                                             0
## 143                                             0
## 144                                             0
## 145                                             0
## 146                                             0
## 147                                             0
## 148                                             0
## 149                                             0
## 150                                             0
## 151                                             0
## 152                                             0
## 153                                             0
## 154                                             0
## 155                                             0
## 156                                             0
## 157                                             0
## 158                                             0
## 159                                             0
## 160                                             0
## 161                                             0
## 162                                             0
## 163                                             0
## 164                                             0
## 165                                             0
## 166                                             0
## 167                                             0
## 168                                             0
## 169                                             0
## 170                                             0
## 171                                             0
## 172                                             0
## 173                                             0
## 174                                             0
## 175                                             0
## 176                                             0
## 177                                             0
## 178                                             0
## 179                                             0
## 180                                             0
## 181                                             0
## 182                                             0
## 183                                             0
## 184                                             0
## 185                                             0
## 186                                             0
## 187                                             0
## 188                                             0
## 189                                             0
## 190                                             0
## 191                                             0
## 192                                             0
## 193                                             0
## 194                                             0
## 195                                             0
## 196                                             0
## 197                                             0
## 198                                             0
## 199                                             0
## 200                                             0
## 201                                             0
## 202                                             0
## 203                                             0
## 204                                             0
## 205                                             0
## 206                                             0
## 207                                             0
## 208                                             0
## 209                                             0
## 210                                             0
## 211                                             0
## 212                                             0
## 213                                             0
## 214                                             0
## 215                                             0
## 216                                             0
## 217                                             0
## 218                                             0
## 219                                             0
## 220                                             0
## 221                                             0
## 222                                             0
## 223                                             0
## 224                                             0
## 225                                             0
## 226                                             0
## 227                                             0
## 228                                             0
## 229                                             0
## 230                                             0
## 231                                             0
## 232                                             0
## 233                                             0
## 234                                             0
## 235                                             0
## 236                                             0
## 237                                             0
## 238                                             0
## 239                                             0
## 240                                             0
## 241                                             0
## 242                                             0
## 243                                             0
## 244                                             0
## 245                                             0
## 246                                             0
## 247                                             0
## 248                                             0
## 249                                             0
## 250                                             0
## 251                                             0
##     sty_indicator_high_risk_fam_hist_yn BRCA1 BRCA2 High Risk Other NA
## 1                                     0     0     1         0     0  0
## 2                                     0     0     0         1     0  0
## 3                                     0     2     0         0     0  0
## 4                                     0     0     0         0     1  0
## 5                                     0     0     0         0     1  0
## 6                                     0     0     0         1     0  0
## 7                                     1     0     2         0     0  0
## 8                                     0     0     0         0     1  0
## 9                                     0     0     0         0     1  0
## 10                                    0     2     0         0     0  0
## 11                                    0     1     0         0     0  0
## 12                                    0     0     0         1     0  0
## 13                                    0     0     0         0     0  1
## 14                                    1     0     0         1     0  0
## 15                                    0     0     0         1     0  0
## 16                                    1     0     0         2     0  0
## 17                                    0     0     0         2     0  0
## 18                                    0     1     0         0     0  0
## 19                                    1     1     0         0     0  0
## 20                                    0     0     0         1     0  0
## 21                                    0     0     0         0     2  0
## 22                                    0     1     0         0     0  0
## 23                                    0     0     0         1     0  0
## 24                                    0     0     0         1     0  0
## 25                                    0     0     3         0     0  0
## 26                                    0     0     0         4     0  0
## 27                                    0     3     0         0     0  0
## 28                                    0     0     2         0     0  0
## 29                                    1     0     0         2     0  0
## 30                                    0     0     0         0     2  0
## 31                                    0     0     0         0     4  0
## 32                                    1     0     0         2     0  0
## 33                                    0     0     0         0     2  0
## 34                                    0     0     0         1     0  0
## 35                                    0     0     0         1     0  0
## 36                                    0     1     0         0     0  0
## 37                                    0     0     0         1     0  0
## 38                                    0     3     0         0     0  0
## 39                                    1     0     2         0     0  0
## 40                                    1     1     0         0     0  0
## 41                                    0     0     0         1     0  0
## 42                                    0     0     0         0     0  2
## 43                                    0     0     0         0     1  0
## 44                                    1     0     1         0     0  0
## 45                                    0     1     0         0     0  0
## 46                                    0     0     1         0     0  0
## 47                                    1     0     0         1     0  0
## 48                                    0     0     0         3     0  0
## 49                                    0     2     0         0     0  0
## 50                                    0     0     1         0     0  0
## 51                                    0     5     0         0     0  0
## 52                                    1     1     0         0     0  0
## 53                                    0     1     0         0     0  0
## 54                                    0     0     1         0     0  0
## 55                                    0     0     2         0     0  0
## 56                                    0     2     0         0     0  0
## 57                                    0     4     0         0     0  0
## 58                                    0     0     1         0     0  0
## 59                                    0     1     0         0     0  0
## 60                                    1     0     0         1     0  0
## 61                                    0     1     0         0     0  0
## 62                                    1     1     0         0     0  0
## 63                                    0     1     0         0     0  0
## 64                                    0     0     1         0     0  0
## 65                                    0     1     0         0     0  0
## 66                                    0     1     0         0     0  0
## 67                                    0     0     1         0     0  0
## 68                                    0     1     0         0     0  0
## 69                                    0     2     0         0     0  0
## 70                                    1     1     0         0     0  0
## 71                                    0     0     1         0     0  0
## 72                                    0     0     0         1     0  0
## 73                                    1     0     0         1     0  0
## 74                                    0     0     0         0     1  0
## 75                                    0     0     0         0     1  0
## 76                                    0     0     0         2     0  0
## 77                                    1     0     0         1     0  0
## 78                                    1     0     0         1     0  0
## 79                                    0     0     0         4     0  0
## 80                                    0     0     0         0     0  1
## 81                                    0     0     3         0     0  0
## 82                                    0     0     0         2     0  0
## 83                                    0     0     0         0     1  0
## 84                                    0     0     0         0     2  0
## 85                                    0     0     0         0     1  0
## 86                                    0     0     0         0     1  0
## 87                                    0     0     0         0     3  0
## 88                                    0     0     0         0     3  0
## 89                                    0     0     0         0     4  0
## 90                                    0     0     0         0     1  0
## 91                                    0     0     0         0     1  0
## 92                                    0     0     0         0     3  0
## 93                                    0     0     0         0     1  0
## 94                                    0     0     0         0     1  0
## 95                                    0     0     0         1     0  0
## 96                                    0     0     0         5     0  0
## 97                                    1     0     0         1     0  0
## 98                                    0     0     0         1     0  0
## 99                                    0     0     0         0     1  0
## 100                                   0     0     0         0     1  0
## 101                                   0     0     0         0     1  0
## 102                                   0     0     0         0     0  2
## 103                                   0     0     0         0     1  0
## 104                                   0     0     0         1     0  0
## 105                                   0     0     0         1     0  0
## 106                                   0     0     0         0     3  0
## 107                                   0     0     0         0     2  0
## 108                                   0     0     0         0     1  0
## 109                                   0     0     0         0     2  0
## 110                                   0     0     0         0     1  0
## 111                                   0     0     0         0     2  0
## 112                                   0     0     0         0     1  0
## 113                                   0     0     0         2     0  0
## 114                                   0     0     0         3     0  0
## 115                                   1     0     0         2     0  0
## 116                                   0     0     0         0     4  0
## 117                                   0     0     0         0     1  0
## 118                                   0     0     0         0     2  0
## 119                                   0     0     0         0     2  0
## 120                                   0     0     0         0     1  0
## 121                                   0     0     0         0     1  0
## 122                                   1     0     0         1     0  0
## 123                                   1     0     0         1     0  0
## 124                                   0     0     0         2     0  0
## 125                                   1     0     0         2     0  0
## 126                                   0     0     0         1     0  0
## 127                                   0     0     0         1     0  0
## 128                                   1     0     0         1     0  0
## 129                                   0     0     0         0     1  0
## 130                                   0     0     0         1     0  0
## 131                                   0     0     0         0     1  0
## 132                                   0     0     0         0     2  0
## 133                                   1     0     0         1     0  0
## 134                                   0     0     0         0     1  0
## 135                                   0     0     0         0     2  0
## 136                                   0     0     0         1     0  0
## 137                                   0     0     0         0     1  0
## 138                                   0     0     0         0     2  0
## 139                                   0     0     0         0     1  0
## 140                                   0     0     0         0     1  0
## 141                                   0     0     0         0     4  0
## 142                                   0     0     0         1     0  0
## 143                                   0     0     0         0     3  0
## 144                                   0     0     0         0     1  0
## 145                                   0     0     0         3     0  0
## 146                                   0     0     0         1     0  0
## 147                                   0     5     0         0     0  0
## 148                                   1     0     0         2     0  0
## 149                                   0     0     0         1     0  0
## 150                                   1     0     0         2     0  0
## 151                                   0     0     0         0     1  0
## 152                                   0     0     0         0     1  0
## 153                                   0     0     0         1     0  0
## 154                                   0     0     0         0     3  0
## 155                                   0     0     0         0     2  0
## 156                                   0     0     0         1     0  0
## 157                                   0     0     0         2     0  0
## 158                                   1     0     0         1     0  0
## 159                                   1     0     0         2     0  0
## 160                                   1     0     0         3     0  0
## 161                                   0     0     0         2     0  0
## 162                                   0     0     0         0     2  0
## 163                                   0     0     0         1     0  0
## 164                                   0     0     0         1     0  0
## 165                                   1     0     0         2     0  0
## 166                                   0     0     0         0     6  0
## 167                                   0     0     0         0     2  0
## 168                                   0     0     0         0     1  0
## 169                                   0     0     0         0     1  0
## 170                                   1     0     0         3     0  0
## 171                                   0     0     0         2     0  0
## 172                                   0     0     0         0     0  1
## 173                                   0     0     0         1     0  0
## 174                                   0     0     0         0     8  0
## 175                                   0     0     0         0     1  0
## 176                                   0     0     0         0     1  0
## 177                                   0     0     0         0     2  0
## 178                                   0     0     0         1     0  0
## 179                                   0     0     0         0     1  0
## 180                                   0     0     0         1     0  0
## 181                                   0     0     0         0     1  0
## 182                                   1     0     0         1     0  0
## 183                                   0     0     0         1     0  0
## 184                                   0     0     0         1     0  0
## 185                                   0     0     0         1     0  0
## 186                                   0     0     0         0     1  0
## 187                                   0     0     0         1     0  0
## 188                                   0     0     0         0     2  0
## 189                                   0     0     0         1     0  0
## 190                                   0     0     0         0     4  0
## 191                                   0     0     0         0     2  0
## 192                                   0     0     0         0     4  0
## 193                                   0     0     0         0     1  0
## 194                                   0     0     0         0     1  0
## 195                                   0     0     0         0     2  0
## 196                                   0     0     0         0     2  0
## 197                                   0     0     0         0     1  0
## 198                                   0     0     0         0     1  0
## 199                                   0     0     0         1     0  0
## 200                                   0     0     0         1     0  0
## 201                                   0     0     0         0     2  0
## 202                                   0     0     0         2     0  0
## 203                                   0     0     0         0     2  0
## 204                                   0     0     0         0     3  0
## 205                                   0     0     0         1     0  0
## 206                                   0     0     0         0     2  0
## 207                                   0     0     0         0     1  0
## 208                                   0     0     0         0     2  0
## 209                                   0     0     0         1     0  0
## 210                                   1     0     3         0     0  0
## 211                                   0     0     3         0     0  0
## 212                                   0     0     0         1     0  0
## 213                                   0     0     0         2     0  0
## 214                                   0     0     0         0     3  0
## 215                                   0     0     0         3     0  0
## 216                                   0     0     0         0     2  0
## 217                                   0     0     0         2     0  0
## 218                                   0     0     0         2     0  0
## 219                                   1     0     0         2     0  0
## 220                                   0     0     0         0     1  0
## 221                                   0     0     0         0     4  0
## 222                                   0     0     0         0     2  0
## 223                                   0     0     0         0     4  0
## 224                                   0     0     0         2     0  0
## 225                                   0     0     0         0     1  0
## 226                                   0     0     0         0     1  0
## 227                                   0     0     0         0     1  0
## 228                                   0     0     0         0     2  0
## 229                                   0     0     0         0    10  0
## 230                                   0     0     0         0     0  2
## 231                                   1     0     0         0     0  1
## 232                                   0     0     0         0     0  2
## 233                                   0     0     0         0     0  1
## 234                                   0     0     0         0     0  1
## 235                                   0     0     0         0     0  2
## 236                                   0     0     0         0     0  1
## 237                                   0     0     0         0     0  1
## 238                                   0     0     0         0     0  1
## 239                                   0     0     0         0     0  2
## 240                                   1     0     0         0     0  2
## 241                                   1     0     0         0     0  3
## 242                                   0     0     0         0     0  1
## 243                                   1     0     0         0     0  1
## 244                                   0     0     0         0     0  1
## 245                                   0     0     0         0     0  1
## 246                                   0     0     0         0     0  2
## 247                                   1     0     0         0     0  3
## 248                                   0     0     0         0     0  1
## 249                                   0     0     0         0     0  1
## 250                                   0     0     0         0     0  1
## 251                                   1     0     0         0     0  1
```

```r


print(cast(lesioninfo_reasons, cad_pt_no_txt + sty_indicator_pre_operative_extent_of_dis_yn + 
    sty_indicator_post_operative_margin_yn + sty_indicator_pre_neoadj_trtmnt_yn + 
    sty_indicator_post_neoadj_trtmnt_yn + sty_indicator_prob_solv_diff_img_yn + 
    sty_indicator_scar_vs_recurr_yn + sty_indicator_folup_recommend_yn + sty_indicator_add_eval_as_folup_yn + 
    sty_indicator_folup_after_pre_exam_yn + sty_indicator_rout_screening_obsp_yn ~ 
    latest_mutation_status_int))
```

```
## Using mri_cad_status_txt as value column.  Use the value argument to cast to override this choice
## Aggregation requires fun.aggregate: length used as default
```

```
##     cad_pt_no_txt sty_indicator_pre_operative_extent_of_dis_yn
## 1               2                                            0
## 2              16                                            0
## 3              18                                            0
## 4              25                                            0
## 5              27                                            0
## 6              59                                            0
## 7              93                                            0
## 8             103                                            0
## 9             111                                            0
## 10            114                                            0
## 11            114                                            0
## 12            121                                            0
## 13            122                                            0
## 14            123                                            0
## 15            130                                            0
## 16            130                                            0
## 17            132                                            0
## 18            133                                            0
## 19            166                                            0
## 20            177                                            0
## 21            180                                            0
## 22            186                                            0
## 23            190                                            0
## 24            196                                            0
## 25            197                                            0
## 26            198                                            0
## 27            220                                            0
## 28            229                                            0
## 29            232                                            0
## 30            232                                            0
## 31            246                                            0
## 32            252                                            0
## 33            259                                            0
## 34            261                                            0
## 35            271                                            0
## 36            272                                            0
## 37            276                                            0
## 38            280                                            0
## 39            282                                            0
## 40            293                                            0
## 41            299                                            0
## 42            311                                            0
## 43            331                                            0
## 44            340                                            0
## 45            409                                            0
## 46            420                                            0
## 47            426                                            0
## 48            455                                            0
## 49            456                                            0
## 50            462                                            0
## 51            465                                            0
## 52            503                                            0
## 53            513                                            0
## 54            547                                            0
## 55            553                                            0
## 56            556                                            0
## 57            559                                            0
## 58            576                                            0
## 59            578                                            0
## 60            580                                            0
## 61            606                                            0
## 62            635                                            0
## 63            651                                            0
## 64            657                                            0
## 65            663                                            0
## 66            664                                            0
## 67            666                                            0
## 68            667                                            0
## 69            668                                            0
## 70            672                                            0
## 71            673                                            0
## 72            679                                            0
## 73            681                                            0
## 74            683                                            0
## 75            684                                            0
## 76            685                                            0
## 77            687                                            0
## 78            689                                            0
## 79            690                                            0
## 80            690                                            0
## 81            691                                            0
## 82            700                                            1
## 83            705                                            1
## 84            706                                            0
## 85            707                                            0
## 86            710                                            0
## 87            713                                            0
## 88            714                                            0
## 89            718                                            0
## 90            721                                            1
## 91            722                                            0
## 92            724                                            0
## 93            726                                            0
## 94            727                                            0
## 95            728                                            0
## 96            729                                            0
## 97            730                                            0
## 98            731                                            0
## 99            735                                            0
## 100           736                                            0
## 101           742                                            0
## 102           744                                            0
## 103           745                                            0
## 104           747                                            0
## 105           752                                            0
## 106           755                                            0
## 107           757                                            1
## 108           758                                            0
## 109           760                                            0
## 110           764                                            0
## 111           765                                            0
## 112           771                                            1
## 113           775                                            0
## 114           775                                            0
## 115           776                                            0
## 116           778                                            0
## 117           779                                            1
## 118           781                                            0
## 119           782                                            0
## 120           783                                            0
## 121           789                                            0
## 122           790                                            0
## 123           791                                            0
## 124           792                                            0
## 125           793                                            0
## 126           793                                            0
## 127           795                                            0
## 128           796                                            0
## 129           799                                            0
## 130           802                                            0
## 131           803                                            0
## 132           805                                            0
## 133           807                                            1
## 134           809                                            0
## 135           810                                            0
## 136           812                                            1
## 137           813                                            0
## 138           814                                            1
## 139           815                                            0
## 140           817                                            0
## 141           817                                            0
## 142           818                                            0
## 143           827                                            0
## 144           829                                            0
## 145           830                                            1
## 146           831                                            1
## 147           837                                            0
## 148           837                                            0
## 149           839                                            0
## 150           843                                            0
## 151           845                                            0
## 152           846                                            0
## 153           847                                            0
## 154           850                                            0
## 155           851                                            0
## 156           852                                            0
## 157           853                                            0
## 158           853                                            0
## 159           855                                            0
## 160           856                                            0
## 161           857                                            0
## 162           857                                            0
## 163           860                                            0
## 164           861                                            0
## 165           862                                            0
## 166           863                                            0
## 167           865                                            0
## 168           865                                            0
## 169           867                                            0
## 170           867                                            0
## 171           870                                            0
## 172           871                                            0
## 173           873                                            1
## 174           875                                            0
## 175           877                                            0
## 176           880                                            0
## 177           881                                            1
## 178           883                                            1
## 179           884                                            0
## 180           887                                            0
## 181           888                                            0
## 182           896                                            0
## 183           900                                            0
## 184           904                                            0
## 185           913                                            0
## 186           918                                            0
## 187           920                                            0
## 188           921                                            0
## 189           937                                            0
## 190           950                                            0
## 191           956                                            0
## 192          6001                                            1
## 193          6004                                            0
## 194          6005                                            1
## 195          6008                                            0
## 196          6014                                            1
## 197          6015                                            1
## 198          6017                                            1
## 199          6018                                            0
## 200          6019                                            0
## 201          6020                                            0
## 202          6021                                            0
## 203          6022                                            0
## 204          6023                                            0
## 205          6024                                            0
## 206          6025                                            0
## 207          6026                                            1
## 208          6027                                            1
## 209          6029                                            0
## 210          6032                                            0
## 211          6034                                            0
## 212          6035                                            0
## 213          6036                                            1
## 214          6037                                            0
## 215          6038                                            1
## 216          6039                                            1
## 217          6040                                            0
## 218          6041                                            0
## 219          6041                                            1
## 220          6042                                            0
## 221          6043                                            0
## 222          6044                                            0
## 223          6044                                            1
## 224          6045                                            1
## 225          6046                                            0
## 226          6047                                            0
## 227          6048                                            0
## 228          6050                                            0
## 229          6051                                            1
## 230          6052                                            0
## 231          6054                                            0
## 232          6054                                            0
## 233          7008                                            1
## 234          7011                                            0
## 235          7018                                            0
## 236          7024                                            0
## 237          7029                                            0
## 238          7030                                            0
## 239          7043                                            0
## 240          7045                                            0
## 241          7054                                            0
## 242          7066                                            0
## 243          7076                                            0
## 244          7077                                            0
## 245          7085                                            0
## 246          7086                                            0
## 247          7088                                            0
## 248          7094                                            0
## 249          7096                                            0
## 250          7097                                            0
## 251          7104                                            0
## 252          7110                                            0
## 253          7127                                            0
## 254          7151                                            0
##     sty_indicator_post_operative_margin_yn
## 1                                        0
## 2                                        0
## 3                                        0
## 4                                        0
## 5                                        0
## 6                                        0
## 7                                        0
## 8                                        0
## 9                                        0
## 10                                       0
## 11                                       0
## 12                                       0
## 13                                       0
## 14                                       0
## 15                                       0
## 16                                       0
## 17                                       0
## 18                                       0
## 19                                       0
## 20                                       0
## 21                                       0
## 22                                       0
## 23                                       0
## 24                                       0
## 25                                       0
## 26                                       0
## 27                                       0
## 28                                       0
## 29                                       0
## 30                                       0
## 31                                       0
## 32                                       0
## 33                                       0
## 34                                       0
## 35                                       0
## 36                                       0
## 37                                       0
## 38                                       0
## 39                                       0
## 40                                       0
## 41                                       0
## 42                                       0
## 43                                       0
## 44                                       0
## 45                                       0
## 46                                       0
## 47                                       0
## 48                                       0
## 49                                       0
## 50                                       0
## 51                                       0
## 52                                       0
## 53                                       0
## 54                                       0
## 55                                       0
## 56                                       0
## 57                                       0
## 58                                       0
## 59                                       0
## 60                                       0
## 61                                       0
## 62                                       0
## 63                                       0
## 64                                       0
## 65                                       0
## 66                                       0
## 67                                       0
## 68                                       0
## 69                                       0
## 70                                       0
## 71                                       0
## 72                                       0
## 73                                       0
## 74                                       0
## 75                                       0
## 76                                       0
## 77                                       0
## 78                                       0
## 79                                       0
## 80                                       0
## 81                                       0
## 82                                       0
## 83                                       0
## 84                                       0
## 85                                       0
## 86                                       0
## 87                                       0
## 88                                       0
## 89                                       0
## 90                                       0
## 91                                       0
## 92                                       0
## 93                                       0
## 94                                       0
## 95                                       0
## 96                                       0
## 97                                       0
## 98                                       0
## 99                                       0
## 100                                      0
## 101                                      0
## 102                                      0
## 103                                      0
## 104                                      0
## 105                                      0
## 106                                      0
## 107                                      0
## 108                                      0
## 109                                      0
## 110                                      0
## 111                                      0
## 112                                      0
## 113                                      0
## 114                                      0
## 115                                      0
## 116                                      0
## 117                                      0
## 118                                      0
## 119                                      0
## 120                                      0
## 121                                      0
## 122                                      0
## 123                                      0
## 124                                      0
## 125                                      0
## 126                                      0
## 127                                      0
## 128                                      0
## 129                                      0
## 130                                      0
## 131                                      0
## 132                                      0
## 133                                      0
## 134                                      0
## 135                                      0
## 136                                      0
## 137                                      0
## 138                                      0
## 139                                      0
## 140                                      0
## 141                                      0
## 142                                      0
## 143                                      0
## 144                                      0
## 145                                      0
## 146                                      0
## 147                                      0
## 148                                      0
## 149                                      0
## 150                                      0
## 151                                      0
## 152                                      0
## 153                                      0
## 154                                      0
## 155                                      0
## 156                                      0
## 157                                      0
## 158                                      0
## 159                                      0
## 160                                      0
## 161                                      0
## 162                                      0
## 163                                      0
## 164                                      0
## 165                                      0
## 166                                      0
## 167                                      0
## 168                                      0
## 169                                      0
## 170                                      0
## 171                                      0
## 172                                      0
## 173                                      0
## 174                                      0
## 175                                      0
## 176                                      0
## 177                                      0
## 178                                      0
## 179                                      0
## 180                                      0
## 181                                      0
## 182                                      0
## 183                                      0
## 184                                      0
## 185                                      0
## 186                                      0
## 187                                      0
## 188                                      0
## 189                                      0
## 190                                      0
## 191                                      0
## 192                                      0
## 193                                      0
## 194                                      0
## 195                                      0
## 196                                      0
## 197                                      0
## 198                                      0
## 199                                      0
## 200                                      0
## 201                                      0
## 202                                      0
## 203                                      0
## 204                                      0
## 205                                      0
## 206                                      0
## 207                                      0
## 208                                      0
## 209                                      0
## 210                                      0
## 211                                      0
## 212                                      0
## 213                                      0
## 214                                      0
## 215                                      0
## 216                                      0
## 217                                      0
## 218                                      0
## 219                                      0
## 220                                      0
## 221                                      0
## 222                                      0
## 223                                      0
## 224                                      0
## 225                                      0
## 226                                      0
## 227                                      0
## 228                                      0
## 229                                      0
## 230                                      0
## 231                                      0
## 232                                      0
## 233                                      0
## 234                                      0
## 235                                      0
## 236                                      0
## 237                                      0
## 238                                      0
## 239                                      0
## 240                                      0
## 241                                      0
## 242                                      0
## 243                                      0
## 244                                      0
## 245                                      0
## 246                                      0
## 247                                      0
## 248                                      0
## 249                                      0
## 250                                      0
## 251                                      0
## 252                                      0
## 253                                      0
## 254                                      0
##     sty_indicator_pre_neoadj_trtmnt_yn sty_indicator_post_neoadj_trtmnt_yn
## 1                                    0                                   0
## 2                                    0                                   0
## 3                                    0                                   0
## 4                                    0                                   0
## 5                                    0                                   0
## 6                                    0                                   0
## 7                                    0                                   0
## 8                                    0                                   0
## 9                                    0                                   0
## 10                                   0                                   0
## 11                                   0                                   0
## 12                                   0                                   0
## 13                                   0                                   0
## 14                                   0                                   0
## 15                                   0                                   0
## 16                                   0                                   0
## 17                                   0                                   0
## 18                                   0                                   0
## 19                                   0                                   0
## 20                                   0                                   0
## 21                                   0                                   0
## 22                                   0                                   0
## 23                                   0                                   0
## 24                                   0                                   0
## 25                                   0                                   0
## 26                                   0                                   0
## 27                                   0                                   0
## 28                                   0                                   0
## 29                                   0                                   0
## 30                                   0                                   0
## 31                                   0                                   0
## 32                                   0                                   0
## 33                                   0                                   0
## 34                                   0                                   0
## 35                                   0                                   0
## 36                                   0                                   0
## 37                                   0                                   0
## 38                                   0                                   0
## 39                                   0                                   0
## 40                                   0                                   0
## 41                                   0                                   0
## 42                                   0                                   0
## 43                                   0                                   0
## 44                                   0                                   0
## 45                                   0                                   0
## 46                                   0                                   0
## 47                                   0                                   0
## 48                                   0                                   0
## 49                                   0                                   0
## 50                                   0                                   0
## 51                                   0                                   0
## 52                                   0                                   0
## 53                                   0                                   0
## 54                                   0                                   0
## 55                                   0                                   0
## 56                                   0                                   0
## 57                                   0                                   0
## 58                                   0                                   0
## 59                                   0                                   0
## 60                                   0                                   0
## 61                                   0                                   0
## 62                                   0                                   0
## 63                                   0                                   0
## 64                                   0                                   0
## 65                                   0                                   0
## 66                                   0                                   0
## 67                                   0                                   0
## 68                                   0                                   0
## 69                                   0                                   0
## 70                                   0                                   0
## 71                                   0                                   0
## 72                                   0                                   0
## 73                                   0                                   0
## 74                                   0                                   0
## 75                                   0                                   0
## 76                                   0                                   0
## 77                                   0                                   0
## 78                                   0                                   0
## 79                                   0                                   0
## 80                                   0                                   0
## 81                                   0                                   0
## 82                                   0                                   0
## 83                                   0                                   0
## 84                                   0                                   0
## 85                                   0                                   0
## 86                                   0                                   0
## 87                                   0                                   0
## 88                                   0                                   0
## 89                                   0                                   0
## 90                                   0                                   0
## 91                                   0                                   0
## 92                                   0                                   0
## 93                                   0                                   0
## 94                                   0                                   0
## 95                                   0                                   0
## 96                                   0                                   0
## 97                                   0                                   0
## 98                                   0                                   0
## 99                                   0                                   0
## 100                                  0                                   0
## 101                                  0                                   0
## 102                                  0                                   0
## 103                                  0                                   0
## 104                                  0                                   0
## 105                                  0                                   0
## 106                                  0                                   0
## 107                                  0                                   0
## 108                                  0                                   0
## 109                                  0                                   0
## 110                                  0                                   0
## 111                                  0                                   0
## 112                                  0                                   0
## 113                                  0                                   0
## 114                                  0                                   0
## 115                                  0                                   0
## 116                                  0                                   0
## 117                                  0                                   0
## 118                                  0                                   0
## 119                                  0                                   0
## 120                                  0                                   0
## 121                                  0                                   0
## 122                                  0                                   0
## 123                                  0                                   0
## 124                                  0                                   0
## 125                                  0                                   0
## 126                                  0                                   0
## 127                                  0                                   0
## 128                                  0                                   0
## 129                                  0                                   0
## 130                                  0                                   0
## 131                                  0                                   0
## 132                                  0                                   0
## 133                                  0                                   0
## 134                                  0                                   0
## 135                                  0                                   0
## 136                                  0                                   0
## 137                                  0                                   0
## 138                                  0                                   0
## 139                                  0                                   0
## 140                                  0                                   0
## 141                                  0                                   0
## 142                                  0                                   0
## 143                                  0                                   0
## 144                                  0                                   0
## 145                                  0                                   0
## 146                                  0                                   0
## 147                                  0                                   0
## 148                                  0                                   0
## 149                                  0                                   0
## 150                                  0                                   0
## 151                                  0                                   0
## 152                                  0                                   0
## 153                                  0                                   0
## 154                                  0                                   0
## 155                                  0                                   0
## 156                                  0                                   0
## 157                                  0                                   0
## 158                                  0                                   0
## 159                                  0                                   0
## 160                                  0                                   0
## 161                                  0                                   0
## 162                                  0                                   0
## 163                                  0                                   0
## 164                                  0                                   0
## 165                                  0                                   0
## 166                                  0                                   0
## 167                                  0                                   0
## 168                                  0                                   0
## 169                                  0                                   0
## 170                                  0                                   0
## 171                                  0                                   0
## 172                                  0                                   0
## 173                                  0                                   0
## 174                                  0                                   0
## 175                                  0                                   0
## 176                                  0                                   0
## 177                                  0                                   0
## 178                                  0                                   0
## 179                                  0                                   0
## 180                                  0                                   0
## 181                                  0                                   0
## 182                                  0                                   0
## 183                                  0                                   0
## 184                                  0                                   0
## 185                                  0                                   0
## 186                                  0                                   0
## 187                                  0                                   0
## 188                                  0                                   0
## 189                                  0                                   0
## 190                                  0                                   0
## 191                                  0                                   0
## 192                                  0                                   0
## 193                                  0                                   0
## 194                                  0                                   0
## 195                                  0                                   0
## 196                                  0                                   0
## 197                                  0                                   0
## 198                                  0                                   0
## 199                                  0                                   0
## 200                                  0                                   0
## 201                                  0                                   0
## 202                                  0                                   0
## 203                                  0                                   0
## 204                                  0                                   0
## 205                                  0                                   0
## 206                                  0                                   0
## 207                                  0                                   0
## 208                                  0                                   0
## 209                                  0                                   0
## 210                                  0                                   0
## 211                                  0                                   0
## 212                                  0                                   0
## 213                                  0                                   0
## 214                                  1                                   0
## 215                                  0                                   0
## 216                                  0                                   0
## 217                                  0                                   0
## 218                                  0                                   0
## 219                                  0                                   0
## 220                                  0                                   0
## 221                                  0                                   0
## 222                                  0                                   0
## 223                                  0                                   0
## 224                                  0                                   0
## 225                                  0                                   0
## 226                                  0                                   0
## 227                                  0                                   0
## 228                                  0                                   0
## 229                                  0                                   0
## 230                                  0                                   0
## 231                                  0                                   0
## 232                                  0                                   0
## 233                                  0                                   0
## 234                                  0                                   0
## 235                                  0                                   0
## 236                                  0                                   0
## 237                                  0                                   0
## 238                                  0                                   0
## 239                                  0                                   0
## 240                                  0                                   0
## 241                                  0                                   0
## 242                                  0                                   0
## 243                                  0                                   0
## 244                                  0                                   0
## 245                                  0                                   0
## 246                                  0                                   0
## 247                                  0                                   0
## 248                                  0                                   0
## 249                                  0                                   0
## 250                                  0                                   0
## 251                                  0                                   0
## 252                                  0                                   0
## 253                                  0                                   0
## 254                                  0                                   0
##     sty_indicator_prob_solv_diff_img_yn sty_indicator_scar_vs_recurr_yn
## 1                                     0                               0
## 2                                     0                               0
## 3                                     0                               0
## 4                                     0                               0
## 5                                     1                               0
## 6                                     0                               0
## 7                                     0                               0
## 8                                     0                               0
## 9                                     0                               0
## 10                                    0                               0
## 11                                    0                               0
## 12                                    1                               0
## 13                                    0                               0
## 14                                    0                               0
## 15                                    0                               0
## 16                                    0                               0
## 17                                    0                               0
## 18                                    0                               0
## 19                                    0                               0
## 20                                    0                               0
## 21                                    0                               0
## 22                                    0                               0
## 23                                    0                               0
## 24                                    1                               0
## 25                                    0                               0
## 26                                    0                               0
## 27                                    0                               0
## 28                                    0                               0
## 29                                    0                               0
## 30                                    0                               0
## 31                                    0                               0
## 32                                    0                               0
## 33                                    0                               0
## 34                                    0                               0
## 35                                    0                               0
## 36                                    0                               0
## 37                                    0                               0
## 38                                    0                               0
## 39                                    0                               0
## 40                                    0                               0
## 41                                    0                               0
## 42                                    0                               0
## 43                                    0                               0
## 44                                    0                               0
## 45                                    0                               0
## 46                                    0                               0
## 47                                    0                               0
## 48                                    0                               0
## 49                                    0                               0
## 50                                    0                               0
## 51                                    0                               0
## 52                                    0                               0
## 53                                    0                               0
## 54                                    0                               0
## 55                                    0                               0
## 56                                    0                               0
## 57                                    0                               0
## 58                                    0                               0
## 59                                    0                               0
## 60                                    0                               0
## 61                                    0                               0
## 62                                    0                               0
## 63                                    0                               0
## 64                                    0                               0
## 65                                    0                               0
## 66                                    0                               0
## 67                                    0                               0
## 68                                    0                               0
## 69                                    0                               0
## 70                                    0                               0
## 71                                    0                               0
## 72                                    0                               0
## 73                                    0                               0
## 74                                    0                               0
## 75                                    0                               0
## 76                                    0                               0
## 77                                    0                               0
## 78                                    0                               0
## 79                                    0                               0
## 80                                    0                               0
## 81                                    0                               0
## 82                                    0                               0
## 83                                    0                               0
## 84                                    0                               0
## 85                                    0                               0
## 86                                    0                               0
## 87                                    0                               0
## 88                                    0                               0
## 89                                    0                               0
## 90                                    0                               0
## 91                                    0                               0
## 92                                    1                               0
## 93                                    0                               0
## 94                                    0                               0
## 95                                    0                               0
## 96                                    0                               0
## 97                                    0                               0
## 98                                    0                               0
## 99                                    0                               0
## 100                                   0                               0
## 101                                   0                               0
## 102                                   0                               0
## 103                                   0                               0
## 104                                   0                               0
## 105                                   0                               0
## 106                                   0                               0
## 107                                   0                               0
## 108                                   0                               0
## 109                                   0                               0
## 110                                   0                               0
## 111                                   0                               0
## 112                                   0                               0
## 113                                   0                               0
## 114                                   1                               0
## 115                                   0                               0
## 116                                   0                               0
## 117                                   0                               0
## 118                                   0                               0
## 119                                   0                               0
## 120                                   0                               0
## 121                                   0                               0
## 122                                   0                               0
## 123                                   0                               0
## 124                                   1                               0
## 125                                   0                               0
## 126                                   0                               0
## 127                                   0                               0
## 128                                   0                               0
## 129                                   0                               0
## 130                                   1                               0
## 131                                   0                               0
## 132                                   1                               0
## 133                                   0                               0
## 134                                   0                               0
## 135                                   0                               0
## 136                                   0                               0
## 137                                   0                               0
## 138                                   0                               0
## 139                                   0                               0
## 140                                   0                               0
## 141                                   0                               0
## 142                                   0                               0
## 143                                   0                               0
## 144                                   0                               0
## 145                                   0                               0
## 146                                   0                               0
## 147                                   0                               0
## 148                                   0                               0
## 149                                   0                               0
## 150                                   0                               0
## 151                                   0                               0
## 152                                   0                               0
## 153                                   0                               0
## 154                                   0                               0
## 155                                   0                               0
## 156                                   0                               0
## 157                                   0                               0
## 158                                   0                               0
## 159                                   0                               0
## 160                                   0                               0
## 161                                   0                               0
## 162                                   0                               0
## 163                                   0                               0
## 164                                   0                               0
## 165                                   0                               0
## 166                                   1                               0
## 167                                   0                               0
## 168                                   0                               0
## 169                                   0                               0
## 170                                   0                               0
## 171                                   0                               0
## 172                                   0                               0
## 173                                   0                               0
## 174                                   0                               0
## 175                                   0                               0
## 176                                   0                               0
## 177                                   0                               0
## 178                                   0                               0
## 179                                   0                               0
## 180                                   0                               0
## 181                                   0                               0
## 182                                   0                               0
## 183                                   0                               0
## 184                                   1                               0
## 185                                   0                               0
## 186                                   0                               0
## 187                                   1                               0
## 188                                   0                               0
## 189                                   0                               0
## 190                                   0                               0
## 191                                   0                               0
## 192                                   0                               0
## 193                                   0                               0
## 194                                   0                               0
## 195                                   0                               0
## 196                                   0                               0
## 197                                   0                               0
## 198                                   0                               0
## 199                                   0                               0
## 200                                   0                               0
## 201                                   0                               0
## 202                                   0                               0
## 203                                   0                               0
## 204                                   0                               0
## 205                                   0                               0
## 206                                   0                               0
## 207                                   0                               0
## 208                                   0                               0
## 209                                   0                               0
## 210                                   0                               0
## 211                                   0                               0
## 212                                   0                               0
## 213                                   0                               0
## 214                                   0                               0
## 215                                   0                               0
## 216                                   0                               0
## 217                                   0                               0
## 218                                   0                               0
## 219                                   0                               0
## 220                                   0                               0
## 221                                   0                               0
## 222                                   0                               0
## 223                                   0                               0
## 224                                   0                               0
## 225                                   0                               0
## 226                                   0                               0
## 227                                   0                               0
## 228                                   0                               0
## 229                                   0                               0
## 230                                   0                               0
## 231                                   0                               0
## 232                                   0                               0
## 233                                   0                               0
## 234                                   0                               0
## 235                                   0                               0
## 236                                   0                               0
## 237                                   0                               0
## 238                                   0                               0
## 239                                   0                               0
## 240                                   0                               0
## 241                                   0                               0
## 242                                   0                               0
## 243                                   0                               0
## 244                                   0                               0
## 245                                   0                               0
## 246                                   0                               0
## 247                                   0                               0
## 248                                   0                               0
## 249                                   0                               0
## 250                                   0                               0
## 251                                   0                               0
## 252                                   0                               0
## 253                                   0                               0
## 254                                   0                               0
##     sty_indicator_folup_recommend_yn sty_indicator_add_eval_as_folup_yn
## 1                                  0                                  0
## 2                                  0                                  0
## 3                                  0                                  0
## 4                                  0                                  0
## 5                                  0                                  1
## 6                                  0                                  0
## 7                                  0                                  0
## 8                                  0                                  1
## 9                                  0                                  1
## 10                                 0                                  0
## 11                                 0                                  1
## 12                                 0                                  0
## 13                                 1                                  0
## 14                                 0                                  0
## 15                                 0                                  0
## 16                                 0                                  0
## 17                                 0                                  1
## 18                                 0                                  0
## 19                                 0                                  0
## 20                                 0                                  0
## 21                                 0                                  0
## 22                                 0                                  0
## 23                                 1                                  0
## 24                                 1                                  0
## 25                                 0                                  0
## 26                                 0                                  0
## 27                                 0                                  0
## 28                                 0                                  1
## 29                                 0                                  0
## 30                                 0                                  1
## 31                                 0                                  0
## 32                                 0                                  0
## 33                                 0                                  0
## 34                                 0                                  0
## 35                                 0                                  0
## 36                                 0                                  1
## 37                                 0                                  0
## 38                                 0                                  0
## 39                                 0                                  0
## 40                                 0                                  0
## 41                                 0                                  0
## 42                                 1                                  0
## 43                                 0                                  0
## 44                                 0                                  0
## 45                                 0                                  1
## 46                                 0                                  0
## 47                                 0                                  0
## 48                                 0                                  1
## 49                                 0                                  0
## 50                                 0                                  0
## 51                                 0                                  0
## 52                                 0                                  0
## 53                                 0                                  0
## 54                                 0                                  1
## 55                                 0                                  0
## 56                                 0                                  0
## 57                                 0                                  0
## 58                                 0                                  0
## 59                                 0                                  1
## 60                                 0                                  0
## 61                                 0                                  0
## 62                                 0                                  0
## 63                                 0                                  0
## 64                                 0                                  0
## 65                                 0                                  0
## 66                                 0                                  0
## 67                                 0                                  1
## 68                                 0                                  0
## 69                                 0                                  0
## 70                                 0                                  1
## 71                                 0                                  0
## 72                                 0                                  1
## 73                                 0                                  1
## 74                                 0                                  1
## 75                                 0                                  0
## 76                                 0                                  1
## 77                                 0                                  0
## 78                                 0                                  0
## 79                                 0                                  0
## 80                                 0                                  0
## 81                                 0                                  0
## 82                                 0                                  0
## 83                                 0                                  0
## 84                                 0                                  0
## 85                                 0                                  1
## 86                                 0                                  1
## 87                                 0                                  1
## 88                                 0                                  0
## 89                                 0                                  1
## 90                                 0                                  0
## 91                                 0                                  0
## 92                                 0                                  0
## 93                                 0                                  0
## 94                                 0                                  0
## 95                                 0                                  0
## 96                                 1                                  0
## 97                                 0                                  0
## 98                                 0                                  1
## 99                                 0                                  1
## 100                                0                                  1
## 101                                0                                  1
## 102                                0                                  1
## 103                                0                                  1
## 104                                0                                  1
## 105                                0                                  1
## 106                                0                                  0
## 107                                0                                  0
## 108                                0                                  0
## 109                                0                                  1
## 110                                0                                  0
## 111                                0                                  1
## 112                                0                                  0
## 113                                0                                  1
## 114                                0                                  0
## 115                                0                                  1
## 116                                0                                  0
## 117                                0                                  0
## 118                                0                                  1
## 119                                0                                  0
## 120                                0                                  0
## 121                                0                                  0
## 122                                0                                  0
## 123                                0                                  0
## 124                                0                                  1
## 125                                0                                  0
## 126                                0                                  0
## 127                                0                                  1
## 128                                0                                  0
## 129                                0                                  1
## 130                                0                                  0
## 131                                0                                  1
## 132                                0                                  1
## 133                                0                                  0
## 134                                0                                  1
## 135                                0                                  0
## 136                                0                                  0
## 137                                0                                  1
## 138                                0                                  0
## 139                                0                                  0
## 140                                0                                  0
## 141                                0                                  1
## 142                                0                                  1
## 143                                0                                  1
## 144                                0                                  0
## 145                                0                                  0
## 146                                0                                  0
## 147                                0                                  0
## 148                                1                                  1
## 149                                0                                  0
## 150                                0                                  0
## 151                                0                                  0
## 152                                0                                  0
## 153                                0                                  1
## 154                                0                                  1
## 155                                0                                  1
## 156                                0                                  1
## 157                                0                                  0
## 158                                0                                  1
## 159                                0                                  0
## 160                                0                                  0
## 161                                0                                  0
## 162                                0                                  1
## 163                                0                                  1
## 164                                0                                  1
## 165                                0                                  0
## 166                                0                                  0
## 167                                0                                  0
## 168                                0                                  1
## 169                                0                                  0
## 170                                0                                  1
## 171                                0                                  0
## 172                                0                                  0
## 173                                0                                  0
## 174                                0                                  0
## 175                                0                                  0
## 176                                0                                  0
## 177                                0                                  0
## 178                                0                                  0
## 179                                0                                  0
## 180                                0                                  0
## 181                                0                                  1
## 182                                1                                  0
## 183                                0                                  1
## 184                                0                                  0
## 185                                0                                  0
## 186                                0                                  0
## 187                                1                                  0
## 188                                1                                  0
## 189                                1                                  0
## 190                                0                                  1
## 191                                1                                  0
## 192                                0                                  0
## 193                                0                                  0
## 194                                0                                  0
## 195                                0                                  0
## 196                                0                                  0
## 197                                0                                  0
## 198                                0                                  0
## 199                                0                                  0
## 200                                0                                  0
## 201                                0                                  0
## 202                                0                                  0
## 203                                0                                  0
## 204                                0                                  0
## 205                                0                                  1
## 206                                0                                  0
## 207                                0                                  0
## 208                                0                                  0
## 209                                0                                  0
## 210                                0                                  0
## 211                                0                                  0
## 212                                0                                  1
## 213                                0                                  0
## 214                                0                                  0
## 215                                0                                  0
## 216                                0                                  0
## 217                                0                                  0
## 218                                0                                  0
## 219                                0                                  0
## 220                                0                                  0
## 221                                0                                  0
## 222                                0                                  0
## 223                                0                                  0
## 224                                0                                  0
## 225                                0                                  0
## 226                                0                                  1
## 227                                0                                  1
## 228                                0                                  0
## 229                                0                                  0
## 230                                0                                  0
## 231                                0                                  0
## 232                                0                                  1
## 233                                0                                  0
## 234                                0                                  0
## 235                                0                                  1
## 236                                0                                  1
## 237                                0                                  0
## 238                                0                                  0
## 239                                0                                  0
## 240                                0                                  0
## 241                                0                                  0
## 242                                0                                  0
## 243                                0                                  0
## 244                                0                                  0
## 245                                0                                  0
## 246                                0                                  0
## 247                                1                                  0
## 248                                0                                  0
## 249                                0                                  1
## 250                                0                                  0
## 251                                0                                  1
## 252                                0                                  1
## 253                                0                                  0
## 254                                0                                  0
##     sty_indicator_folup_after_pre_exam_yn
## 1                                       0
## 2                                       0
## 3                                       0
## 4                                       0
## 5                                       0
## 6                                       0
## 7                                       0
## 8                                       0
## 9                                       0
## 10                                      0
## 11                                      0
## 12                                      0
## 13                                      0
## 14                                      0
## 15                                      0
## 16                                      0
## 17                                      0
## 18                                      0
## 19                                      0
## 20                                      0
## 21                                      0
## 22                                      0
## 23                                      0
## 24                                      0
## 25                                      0
## 26                                      0
## 27                                      0
## 28                                      0
## 29                                      0
## 30                                      0
## 31                                      0
## 32                                      1
## 33                                      0
## 34                                      0
## 35                                      0
## 36                                      0
## 37                                      0
## 38                                      0
## 39                                      0
## 40                                      0
## 41                                      0
## 42                                      0
## 43                                      0
## 44                                      0
## 45                                      0
## 46                                      0
## 47                                      0
## 48                                      0
## 49                                      0
## 50                                      0
## 51                                      0
## 52                                      0
## 53                                      0
## 54                                      0
## 55                                      0
## 56                                      0
## 57                                      0
## 58                                      0
## 59                                      0
## 60                                      0
## 61                                      0
## 62                                      0
## 63                                      1
## 64                                      0
## 65                                      0
## 66                                      0
## 67                                      0
## 68                                      0
## 69                                      0
## 70                                      0
## 71                                      0
## 72                                      0
## 73                                      0
## 74                                      0
## 75                                      0
## 76                                      0
## 77                                      0
## 78                                      1
## 79                                      0
## 80                                      1
## 81                                      0
## 82                                      0
## 83                                      0
## 84                                      0
## 85                                      0
## 86                                      0
## 87                                      0
## 88                                      1
## 89                                      0
## 90                                      0
## 91                                      1
## 92                                      0
## 93                                      0
## 94                                      0
## 95                                      0
## 96                                      0
## 97                                      1
## 98                                      0
## 99                                      0
## 100                                     0
## 101                                     0
## 102                                     0
## 103                                     0
## 104                                     0
## 105                                     0
## 106                                     0
## 107                                     0
## 108                                     0
## 109                                     0
## 110                                     0
## 111                                     0
## 112                                     0
## 113                                     0
## 114                                     0
## 115                                     0
## 116                                     1
## 117                                     0
## 118                                     0
## 119                                     0
## 120                                     1
## 121                                     0
## 122                                     0
## 123                                     0
## 124                                     0
## 125                                     0
## 126                                     0
## 127                                     0
## 128                                     0
## 129                                     0
## 130                                     0
## 131                                     0
## 132                                     0
## 133                                     0
## 134                                     0
## 135                                     1
## 136                                     0
## 137                                     0
## 138                                     0
## 139                                     0
## 140                                     0
## 141                                     0
## 142                                     0
## 143                                     0
## 144                                     0
## 145                                     0
## 146                                     0
## 147                                     0
## 148                                     0
## 149                                     0
## 150                                     0
## 151                                     0
## 152                                     0
## 153                                     0
## 154                                     0
## 155                                     0
## 156                                     0
## 157                                     1
## 158                                     0
## 159                                     0
## 160                                     0
## 161                                     0
## 162                                     0
## 163                                     0
## 164                                     0
## 165                                     0
## 166                                     0
## 167                                     0
## 168                                     0
## 169                                     0
## 170                                     0
## 171                                     0
## 172                                     0
## 173                                     0
## 174                                     1
## 175                                     0
## 176                                     0
## 177                                     0
## 178                                     0
## 179                                     0
## 180                                     0
## 181                                     0
## 182                                     1
## 183                                     0
## 184                                     1
## 185                                     0
## 186                                     0
## 187                                     0
## 188                                     0
## 189                                     0
## 190                                     0
## 191                                     0
## 192                                     0
## 193                                     0
## 194                                     0
## 195                                     0
## 196                                     0
## 197                                     0
## 198                                     0
## 199                                     0
## 200                                     0
## 201                                     0
## 202                                     0
## 203                                     0
## 204                                     0
## 205                                     0
## 206                                     0
## 207                                     0
## 208                                     0
## 209                                     0
## 210                                     0
## 211                                     0
## 212                                     0
## 213                                     0
## 214                                     0
## 215                                     0
## 216                                     0
## 217                                     0
## 218                                     0
## 219                                     0
## 220                                     0
## 221                                     0
## 222                                     0
## 223                                     0
## 224                                     0
## 225                                     0
## 226                                     0
## 227                                     0
## 228                                     0
## 229                                     0
## 230                                     0
## 231                                     0
## 232                                     0
## 233                                     0
## 234                                     0
## 235                                     0
## 236                                     0
## 237                                     0
## 238                                     0
## 239                                     0
## 240                                     1
## 241                                     0
## 242                                     0
## 243                                     0
## 244                                     0
## 245                                     0
## 246                                     0
## 247                                     0
## 248                                     0
## 249                                     0
## 250                                     0
## 251                                     0
## 252                                     0
## 253                                     1
## 254                                     0
##     sty_indicator_rout_screening_obsp_yn BRCA1 BRCA2 High Risk Other NA
## 1                                      0     0     1         0     0  0
## 2                                      1     0     0         1     0  0
## 3                                      0     2     0         0     0  0
## 4                                      0     0     0         0     1  0
## 5                                      0     0     0         0     1  0
## 6                                      0     0     0         1     0  0
## 7                                      1     0     2         0     0  0
## 8                                      0     0     0         0     1  0
## 9                                      0     0     0         0     1  0
## 10                                     0     1     0         0     0  0
## 11                                     0     2     0         0     0  0
## 12                                     0     0     0         1     0  0
## 13                                     0     0     0         0     0  1
## 14                                     1     0     0         1     0  0
## 15                                     0     0     0         2     0  0
## 16                                     1     0     0         1     0  0
## 17                                     0     0     0         2     0  0
## 18                                     0     2     0         0     0  0
## 19                                     0     0     0         1     0  0
## 20                                     0     0     0         0     2  0
## 21                                     0     1     0         0     0  0
## 22                                     0     0     0         2     0  0
## 23                                     0     0     3         0     0  0
## 24                                     0     0     0         4     0  0
## 25                                     0     3     0         0     0  0
## 26                                     0     0     2         0     0  0
## 27                                     0     0     0         2     0  0
## 28                                     0     0     0         0     2  0
## 29                                     0     0     0         0     2  0
## 30                                     0     0     0         0     2  0
## 31                                     0     0     0         2     0  0
## 32                                     0     0     0         0     2  0
## 33                                     1     0     0         1     0  0
## 34                                     0     0     0         1     0  0
## 35                                     0     1     0         0     0  0
## 36                                     0     0     0         1     0  0
## 37                                     1     3     0         0     0  0
## 38                                     0     0     2         0     0  0
## 39                                     0     1     0         0     0  0
## 40                                     1     0     0         1     0  0
## 41                                     0     0     0         0     0  2
## 42                                     0     0     0         0     1  0
## 43                                     0     0     1         0     0  0
## 44                                     0     1     0         0     0  0
## 45                                     0     0     1         0     0  0
## 46                                     0     0     0         1     0  0
## 47                                     1     0     0         3     0  0
## 48                                     0     2     0         0     0  0
## 49                                     0     0     1         0     0  0
## 50                                     0     5     0         0     0  0
## 51                                     0     1     0         0     0  0
## 52                                     0     1     0         0     0  0
## 53                                     0     0     1         0     0  0
## 54                                     0     0     2         0     0  0
## 55                                     0     2     0         0     0  0
## 56                                     0     4     0         0     0  0
## 57                                     0     0     1         0     0  0
## 58                                     1     1     0         0     0  0
## 59                                     0     0     0         1     0  0
## 60                                     0     1     0         0     0  0
## 61                                     0     1     0         0     0  0
## 62                                     1     1     0         0     0  0
## 63                                     0     0     1         0     0  0
## 64                                     1     1     0         0     0  0
## 65                                     0     1     0         0     0  0
## 66                                     1     0     1         0     0  0
## 67                                     0     1     0         0     0  0
## 68                                     0     3     0         0     0  0
## 69                                     1     0     1         0     0  0
## 70                                     0     0     0         1     0  0
## 71                                     0     0     0         1     0  0
## 72                                     0     0     0         0     1  0
## 73                                     0     0     0         0     1  0
## 74                                     0     0     0         2     0  0
## 75                                     0     0     0         1     0  0
## 76                                     0     0     0         1     0  0
## 77                                     0     0     0         4     0  0
## 78                                     0     0     0         0     0  1
## 79                                     0     0     2         0     0  0
## 80                                     0     0     1         0     0  0
## 81                                     0     0     0         2     0  0
## 82                                     0     0     0         0     1  0
## 83                                     0     0     0         0     2  0
## 84                                     0     0     0         0     1  0
## 85                                     0     0     0         0     1  0
## 86                                     0     0     0         0     3  0
## 87                                     0     0     0         0     3  0
## 88                                     0     0     0         0     4  0
## 89                                     0     0     0         0     1  0
## 90                                     0     0     0         0     1  0
## 91                                     0     0     0         0     3  0
## 92                                     0     0     0         0     1  0
## 93                                     0     0     0         0     1  0
## 94                                     0     0     0         1     0  0
## 95                                     0     0     0         5     0  0
## 96                                     0     0     0         1     0  0
## 97                                     0     0     0         1     0  0
## 98                                     0     0     0         0     1  0
## 99                                     0     0     0         0     1  0
## 100                                    0     0     0         0     1  0
## 101                                    0     0     0         0     0  2
## 102                                    0     0     0         0     1  0
## 103                                    0     0     0         1     0  0
## 104                                    0     0     0         1     0  0
## 105                                    0     0     0         0     3  0
## 106                                    0     0     0         0     2  0
## 107                                    0     0     0         0     1  0
## 108                                    0     0     0         0     2  0
## 109                                    0     0     0         0     1  0
## 110                                    0     0     0         0     2  0
## 111                                    0     0     0         0     1  0
## 112                                    0     0     0         2     0  0
## 113                                    0     0     0         3     0  0
## 114                                    0     0     0         2     0  0
## 115                                    0     0     0         0     4  0
## 116                                    0     0     0         0     1  0
## 117                                    0     0     0         0     2  0
## 118                                    0     0     0         0     2  0
## 119                                    0     0     0         0     1  0
## 120                                    0     0     0         0     1  0
## 121                                    0     0     0         1     0  0
## 122                                    0     0     0         1     0  0
## 123                                    0     0     0         2     0  0
## 124                                    0     0     0         2     0  0
## 125                                    0     0     0         1     0  0
## 126                                    1     0     0         1     0  0
## 127                                    0     0     0         1     0  0
## 128                                    0     0     0         0     1  0
## 129                                    0     0     0         1     0  0
## 130                                    0     0     0         0     1  0
## 131                                    0     0     0         0     2  0
## 132                                    0     0     0         1     0  0
## 133                                    0     0     0         0     1  0
## 134                                    0     0     0         0     2  0
## 135                                    0     0     0         1     0  0
## 136                                    0     0     0         0     1  0
## 137                                    0     0     0         0     2  0
## 138                                    0     0     0         0     1  0
## 139                                    0     0     0         0     1  0
## 140                                    0     0     0         0     2  0
## 141                                    0     0     0         0     2  0
## 142                                    0     0     0         1     0  0
## 143                                    0     0     0         0     3  0
## 144                                    0     0     0         0     1  0
## 145                                    0     0     0         3     0  0
## 146                                    0     0     0         1     0  0
## 147                                    0     4     0         0     0  0
## 148                                    0     1     0         0     0  0
## 149                                    0     0     0         2     0  0
## 150                                    0     0     0         3     0  0
## 151                                    0     0     0         0     1  0
## 152                                    0     0     0         0     1  0
## 153                                    0     0     0         1     0  0
## 154                                    0     0     0         0     3  0
## 155                                    0     0     0         0     2  0
## 156                                    0     0     0         1     0  0
## 157                                    0     0     0         2     0  0
## 158                                    0     0     0         1     0  0
## 159                                    0     0     0         2     0  0
## 160                                    0     0     0         3     0  0
## 161                                    0     0     0         1     0  0
## 162                                    0     0     0         1     0  0
## 163                                    0     0     0         0     2  0
## 164                                    0     0     0         1     0  0
## 165                                    0     0     0         1     0  0
## 166                                    0     0     0         2     0  0
## 167                                    0     0     0         0     3  0
## 168                                    0     0     0         0     3  0
## 169                                    0     0     0         0     1  0
## 170                                    0     0     0         0     1  0
## 171                                    0     0     0         0     2  0
## 172                                    0     0     0         3     0  0
## 173                                    0     0     0         2     0  0
## 174                                    0     0     0         0     0  1
## 175                                    0     0     0         1     0  0
## 176                                    0     0     0         0     8  0
## 177                                    0     0     0         0     1  0
## 178                                    0     0     0         0     1  0
## 179                                    0     0     0         0     2  0
## 180                                    0     0     0         1     0  0
## 181                                    0     0     0         0     1  0
## 182                                    0     0     0         1     0  0
## 183                                    0     0     0         0     1  0
## 184                                    0     0     0         1     0  0
## 185                                    1     0     0         1     0  0
## 186                                    1     0     0         1     0  0
## 187                                    0     0     0         1     0  0
## 188                                    0     0     0         0     1  0
## 189                                    0     0     0         1     0  0
## 190                                    0     0     0         0     2  0
## 191                                    0     0     0         1     0  0
## 192                                    0     0     0         0     4  0
## 193                                    0     0     0         0     2  0
## 194                                    0     0     0         0     4  0
## 195                                    0     0     0         0     1  0
## 196                                    0     0     0         0     1  0
## 197                                    0     0     0         0     2  0
## 198                                    0     0     0         0     2  0
## 199                                    0     0     0         0     1  0
## 200                                    0     0     0         0     1  0
## 201                                    0     0     0         1     0  0
## 202                                    0     0     0         1     0  0
## 203                                    0     0     0         0     2  0
## 204                                    0     0     0         2     0  0
## 205                                    0     0     0         0     2  0
## 206                                    0     0     0         0     3  0
## 207                                    0     0     0         1     0  0
## 208                                    0     0     0         0     2  0
## 209                                    0     0     0         0     3  0
## 210                                    0     0     0         1     0  0
## 211                                    0     0     3         0     0  0
## 212                                    0     0     3         0     0  0
## 213                                    0     0     0         1     0  0
## 214                                    0     0     0         2     0  0
## 215                                    0     0     0         0     3  0
## 216                                    0     0     0         3     0  0
## 217                                    0     0     0         0     2  0
## 218                                    0     0     0         2     0  0
## 219                                    0     0     0         2     0  0
## 220                                    0     0     0         2     0  0
## 221                                    0     0     0         0     1  0
## 222                                    0     0     0         0     2  0
## 223                                    0     0     0         0     2  0
## 224                                    0     0     0         0     2  0
## 225                                    0     0     0         0     4  0
## 226                                    0     0     0         2     0  0
## 227                                    0     0     0         0     1  0
## 228                                    1     0     0         0     1  0
## 229                                    0     0     0         0     1  0
## 230                                    0     0     0         0     2  0
## 231                                    0     0     0         0     5  0
## 232                                    0     0     0         0     5  0
## 233                                    0     0     0         0     0  2
## 234                                    0     0     0         0     0  1
## 235                                    0     0     0         0     0  2
## 236                                    0     0     0         0     0  1
## 237                                    1     0     0         0     0  1
## 238                                    1     0     0         0     0  2
## 239                                    0     0     0         0     0  1
## 240                                    0     0     0         0     0  1
## 241                                    0     0     0         0     0  1
## 242                                    0     0     0         0     0  2
## 243                                    1     0     0         0     0  2
## 244                                    0     0     0         0     0  3
## 245                                    1     0     0         0     0  1
## 246                                    0     0     0         0     0  1
## 247                                    0     0     0         0     0  1
## 248                                    0     0     0         0     0  1
## 249                                    0     0     0         0     0  2
## 250                                    0     0     0         0     0  3
## 251                                    0     0     0         0     0  1
## 252                                    0     0     0         0     0  1
## 253                                    0     0     0         0     0  1
## 254                                    1     0     0         0     0  1
```

```r

save.image("Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/experimentsRadiologypaper-revision/Tree-based-RF/ensemble-Treebased-RF/biomatrix_results_analyze_populationstats.RData")
```



