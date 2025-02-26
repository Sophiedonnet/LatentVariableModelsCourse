
# Meadow vegetation data --------------------------------------------------

# Data from 
# Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects,
# by Lepvs, de Bell, Vsmilauer and Dole (2011)

library(tidyverse)
library(traitor)
data(ohrazeni)
# Get only numerical data and remove outliers (huge seed weight)
meadow_data <- select_if(ohrazeni$traitdata, is.numeric) %>% 
  filter(seedweight < 3) %>% 
  rownames_to_column(var = "Species") %>% 
  mutate(Species = str_to_sentence(Species)) %>% 
  rename_all(str_to_title) %>% 
  rename(SLA = Sla, LDMC = Ldmc, "Seed weight" = "Seedweight")

write.table(meadow_data,
            "Bohemia_vegetation.csv",
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

# Phytoplankton data ------------------------------------------------------

# Best mixture has only one component...

# Table extracted from the supplementary material of 
# "A Trait-Based Clustering for Phytoplankton Biomass Modeling and Prediction" 
# from Crispin M. Mutshinda1, Zoe V. Finkel2, Claire E. Widdicombe and Andrew J. Irwin
# available at http://www.mdpi.com/1424-2818/12/8/295/s1 (last page)

library(tidyverse) # For data manipulation
betaPresAbs=structure(.Data=c(-1.7804, -3.3985, -0.9553, -3.3017, -4.6143, -0.913, -3.9706, 0.2541, -3.9386, -0.8467, -3.3948, -0.8613, -5.216, 1.2042, -1.5271, -1.0645, -
                                3.0575, -5.0847, -3.8777, -3.5276, -2.6594, -3.4337, -4.0988, -2.3227, -1.761, -1.9429,
                              -1.7882, -4.4765, -0.1549, -2.7757, -2.8001, -0.254, -5.6402, -2.9366, -5.2134, -2.3857,
                              -3.3872, -0.4646, -3.2497, -2.7591, -1.4039, 0.9388, -1.0245, -2.3736, -1.7159, -4.024,
                              -5.1915, -0.2813, -5.6714, -1.2517, -4.2411, -3.1646, -3.5243, -1.3897, -0.9836, -
                                0.2551, -4.282, 1.0114, 2.3744, 0.2159, 2.6663, 2.7818, 4.6187, 1.8405, 1.4967, 3.1075,
                              3.4975, 1.9648, -1.6284, 3.3773, -0.6049, -2.164, 3.8184, -0.3598, 0.7395, -1.5339,
                              0.9599, -3.3392, -3.5744, 0.6017, -2.5201, 3.1398, -0.7818, 0.7403, -1.2601, -1.6011, -
                                3.005, -1.9843, 3.7159, -3.4616, 2.4238, -0.3936, 0.456, -3.6889, -2.6073, -1.6753, -
                                0.3039, 1.6336, -4.1608, -1.1272, -4.6249, -0.6035, -3.5953, -4.1946, -1.7987, 1.0409,
                              -4.7419, -0.3685, -0.552, 0.3437, -0.0851, 1.8895, -1.2578, -2.045, 1.9628, 2.6537, -
                                1.9244, -3.5454, -3.673, 0.1594, -3.2286, -0.6222, -3.2884, -2.4181, -0.2293, 2.3696,
                              1.5711, -1.7666, -2.2756, 0.1544, -3.2906, 0.9742, 3.3775, 1.4281, 4.8336, 0.4363,
                              3.9934, 1.0711, 0.331, 1.662, 3.3096, -0.6149, -1.2966, 1.7187, 3.5392, 2.7204, 4.0562,
                              0.3868, 1.1353, 0.4253, 2.1873, -0.8954, 2.3027, 1.817, 0.4272, -0.1727, -0.3549, -
                                0.435, 0.3848, -0.453, 1.2762, 0.3788, -0.7983, 1.4031, 1.5708, 0.9315, 0.7879, -0.0226,
                              4.4986, -2.108, 0.4411, -0.7583, 4.163, 0.6375, 1.4684, -1.0362, -0.1055, -0.1408,
                              1.8273, -0.5575, -2.1447, 1.2632, 1.6627, 1.819, 0.1754, 0.0251, 1.0962, -1.7865, -
                                2.9985, 2.3285, -2.1464, 1.7101, 2.2795, -1.7103, -0.911, -2.6259, 0.1163, 1.9046,
                              1.5747, 2.5197, -4.7137, 0.0149, 1.6448, 1.7346, -1.2939, 1.161, -2.8731, 1.0056, -
                                0.2309, 3.2278, 1.1206, 0.2882, 1.0474, 0.0662, 3.0801, 0.0707, -2.8912, 1.7626, -
                                1.1496, -1.6043, -1.1104, -0.5768, 4.8278, -5.0819, 0.3675, 5.2247, -6.42, 4.1182, -
                                7.3498, -1.6833, 5.4448, -2.1967, -1.3697, -2.5, 1.6821, -4.6621, -3.2013, -0.4862, -
                                4.4853, -1.9304, -2.0543, -3.7109, -0.0487, -3.8126, 1.3937, -1.9217, -2.4488, 0.2698,
                              -0.2659, -0.9579, 1.7902, -3.6302, -1.7792, -3.5458, -3.6109, -2.3181, 1.6374, 0.5882,
                              2.5113, 1.515, 3.1349, 0.0195, 0.9058, 2.1182, -0.3746, -2.8558, 0.3578, -2.5913, -
                                2.9484, -4.1801, -1.3984, -1.5187, -4.3928, 3.8878, -2.319, -2.6252, 0.3831, 3.8862, -
                                3.6841, -1.6164, -1.4201, 1.2517, -0.0614, -3.8044, -0.159, -2.2885, -1.5097, -2.0067,-0.0886, -2.7, -2.2762, -2.4815, -1.672, 0.0901, -1.0598, -2.3694, 1.0326, -0.0317, -
                                0.0865, -0.9288, -0.9044, -1.2203, -0.9371, -0.8732, -2.4214, -2.1623, -1.3172, -2.1758,
                              0.3228, 2.5705, -0.4122, -0.554, 0.3453, -2.1714, 0.3394, 1.4702, 1.6039, 1.3111, -
                                1.4815, -0.4735, -2.376, -1.3784, 0.2193, -0.8861, -2.097, 0.6051, -2.2315, 0.6377, -
                                2.2749, 0.1472, 0.5074, -1.9099, 1.4767, 1.3341, -1.8107, 2.2175, -1.0907, -1.8925, -
                                2.378, 1.8482, -0.4944, 1.0401, -0.6884, -2.659, 1.332, 1.0247, -0.6522, -2.8289, 2.3857,
                              1.5086, -0.029, -1.8034, 1.4775, 1.614, -1.2942, 2.2604, 2.4554, 2.2381, 2.2747, -
                                1.2557, -1.3189, -0.3369, 2.0862, 0.6612, 0.2964, 1.5578, 0.9364, -0.0879, 0.6671,
                              0.3726, -2.7634, -1.3228, -2.836, -2.9884, 0.9977, -0.6052, 2.2188, -1.8147, -3.1825, -
                                0.2542, 0.0371, -2.7301, -1.3749, -3.5535, -2.1568, 0.2325, 4.2499, 0.366, -0.5492, -
                                1.5655, 1.3613, -2.3621, 0.2167, -3.3304, -2.0645, -1.1004, -0.7491, 2.0104, -3.5807, -
                                3.6936, -3.1603, -2.8523, -4.6696, -1.1794, 1.937, -1.4851, -0.8839, 1.1368, -0.5223, -
                                0.3726, 2.6465, 0.8076, -0.775, -2.3085, -0.8969, 1.4361, -1.951, 1.2919, 3.7678, 3.7227,
                              2.0886, 1.0931, 3.0334, 0.8949, 1.2544, 0.3863, -1.1673, 1.3015, 2.608, 0.3413, 1.4305,
                              0.0825, 0.8561, -0.5248, 1.7528, -0.4091, 1.0045, -1.8824, -1.3967, -0.575, 1.3541, -
                                1.8486, 0.8975, -0.8737), .Dim = c(74, 6),
                      .Dimnames = list(c("N.clost", "G.delic", "P-n.delic", "P.sulcata", "M.membr", "Pleurosigma", "C.pelagica", "P-n.seriata", "D.cabro",
                                         "L.annulata", "Thalass.10μm", "Ch.debilis", "N.distans", "P.alata.5μm", "L.danicus",
                                         "C.danicus", "E.zodiacus", "T.nitzsch", "G.striata", "G.flaccida", "N.sigmoidea",
                                         "D.fragil", "R.tesselata", "R.setigera.5μm", "C.densus", "Navicula.sp.",
                                         "C.criophilum", "D.brightwel", "S.costatum", "R.imbric.5μm", "R.imbric.15μm",
                                         "L.minimus", "T.rotula", "C.affinis", "O.mobil", "C.decipiens", "Pennate.50μm",
                                         "R.imbric.10μm", "P.stelligera", "P.plancton", "Thalass.20μm", "Thalass.4μm",
                                         "C.simplex", "R.stylif", "B.paradoxa", "R.setigera.25μm", "P-n.pungens", "C.socialis",
                                         "T.punctigera", "Small.Pennate", "P.truncata", "Pennate.30μm", "L.mediterr", "P.alata",
                                         "C.radiatus", "P.pandurif", "D.pumila", "C.fusus", "C.horridum", "C.lineatum",
                                         "C.tripos", "D.acuminata", "K.mikimotoi", "G.spinifera", "Gymnod.sp", "G.cf.pygmaeum",
                                         "M.perforatus", "Micranthodinium.sp.", "P.balticum", "P.micans", "Ptrum.minimum",
                                         "P.triestinum", "S.trochoidea", "Scrip.sp.cyst"), 
                                       c("PAR", "Temperature", "Salinity",
                                         "Nitrogen", "Silicate", "Phosphate"))) |>
  as.data.frame() |> 
  rename(Irradiance = PAR) 

# Set proper output path before running

write.table(betaPresAbs,
            "Phytoplankton.csv",
            sep = ",",
            col.names = TRUE,
            row.names = TRUE)
