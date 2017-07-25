chd_OR = c(1.15,1.34,1.20,1.19,1.40,1.40)
alz_hetero_17 = c(0.3150501672, 0.4856020311, 0.4246593641, 0.2474167449,0.3535061795,0.4564337825,
                  0.4387864207,0.4817572099,0.2689378758,0.4495394473,0.3088490843,0.0987604958,0.4439470871,
                  0.14939228,0.1115641711,0.01945369753,0.3287890938,0.4933253237,0.3134787846,0.2221405665,
                  0.3962420601,0.4237764108,0.03914899293,0.4762381524,0.2108796931,0.1228491397,0.4440208668,
                  0.3643421229,0.08186077644,0.1133136887)
alz_homo_17 = c(0.07785953177, 0.3441693047, 0.1860806839, 0.04255931638, 0.1042450296, 0.2467106311, 0.2143811815,
                0.3265593561, 0.05263861055, 0.2265385129, 0.07238337121, 0.005597760896, 0.2241448423, 0.01409109122,
                0.007620320856,	0.0002664890073, 0.08740978348,	0.3936056601, 0.0751572748, 0.03180117584,
                0.1478456773, 0.1857448516, 0.00133386688,	0.3161126685,	0.02671964922,	0.009603841537,
                0.2178972713,	0.1131709276,	0.00421686747,	0.006094294133)
alz_logOR_17 = c(-0.15, -0.25, 0.21, 0.16, 0.09, 0.11, 0.28, 0.18, -0.09, 0.3, 0.08, 0.14, -0.23, 0.1, -0.11,
                 -0.23, 0.09, -0.11, -0.07, 0.09, 0.07, -0.08, 0.18, -0.09, -0.07, -0.09, -0.2, -0.06, 0.1, 0.08)
chd_logOR = c(0.06069784035,	0.1271047984,	0.08092190762,	0.07554696139,	0.1461280357,	0.1461280357)
chd_hetero = c(0.4573400251,	0.07320099256,	0.1652281134,	0.3064918851,	0.05778894472,	0.1171003717)
chd_homo = c(0.2321204517,	0.002481389578,	0.8224414303,	0.07553058677,	0.003768844221,	0.006815365551)


#returns a simulated genotype using heterozygote and homozygote rates in the target population. 
simulate_genotype_2 <- function(hetero_count,homo_count) {
  SNPs_in_model = length(hetero_count)
  genotype = array(0,SNPs_in_model)
  for(SNP in 1:SNPs_in_model) {
    zygote_distribution = c((1-hetero_count[SNP]-homo_count[SNP]),hetero_count[SNP],homo_count[SNP])
    genotype[SNP] = sample(0:2,1,prob=zygote_distribution)
  }
  return(genotype)
}
#returns a single possible progeny genotype from 2 parental genotypes
simulate_progeny <- function(mother_genotype,father_genotype) {
  SNPs_in_model = length(mother_genotype)
  progeny_genotype = array(0, SNPs_in_model)
  for(SNP in 1:SNPs_in_model) {
    if(mother_genotype[SNP] == 0 && father_genotype[SNP] == 0) {
      progeny_genotype[SNP] <- 0
    } else if (mother_genotype[SNP] == 2 && father_genotype [SNP] == 2) {
      progeny_genotype[SNP] <- 2
    } else if (mother_genotype[SNP] == 1 && father_genotype [SNP] == 2) {
      progeny_genotype[SNP] <- sample(1:2, 1)
    } else if (mother_genotype[SNP] == 2 && father_genotype [SNP] == 1) {
      progeny_genotype[SNP] <- sample(1:2, 1)
    } else if (mother_genotype[SNP] == 0 && father_genotype [SNP] == 1) {
      progeny_genotype[SNP] <- sample(0:1, 1)
    } else if (mother_genotype[SNP] == 1 && father_genotype [SNP] == 0) {
      progeny_genotype[SNP] <- sample(0:1, 1)
    } else if (mother_genotype[SNP] == 2 && father_genotype [SNP] == 0) {
      progeny_genotype[SNP] <- 1
    } else if (mother_genotype[SNP] == 0 && father_genotype [SNP] == 2) {
      progeny_genotype[SNP] <- 1
    } else if (mother_genotype[SNP] == 1 && father_genotype [SNP] == 1) {
      progeny_genotype[SNP] <- sample(0:2, 1, prob= c(0.25,0.50,0.25))
    }
  }
  return(progeny_genotype)
}
#calcualtes GRS by adding up logOdds multiplied by allele frequency at each SNP.
calculate_GRS_from_logOR <-function(genotype,log_OR=alz_logOR_17) {
  log_GRS = sum(genotype*log_OR)
  return(log_GRS)
}
#simulates a mother & father and outputs a vector of GRS scores of progeny embryos
model_all_embryo_GRS_2 <-function(mother_genotype,father_genotype,log_OR=alz_logOR_17,
                                  num_embryos=15,hetero_count=alz_hetero_17, homo_count = alz_homo_17) {
  GRS_scores = array(0,num_embryos)
  for(embryo in 1:num_embryos) {
    embryo_genotype = simulate_progeny(mother_genotype, father_genotype)
    GRS_scores[embryo] = calculate_GRS_from_logOR(embryo_genotype,log_OR)
  }
  return(GRS_scores)
}
#returns number of risk alleles in a given genotype.
count_risk_alleles <- function(genotype) {
  return(sum(genotype))
}
count_no_zero_risk <- function(genotype) {
  return(sum(genotype !=0))
}
#returns a vector of the GRS of the x-healthiest embryo available to a sample_size of simulated parental couples.
find_healthiest_embryos_2 <-function(log_OR, homo_count, hetero_count, num_embryos,sample_size, embryo_cutoff) {
  chosen_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    GRS_scores = model_all_embryo_GRS_2(log_OR,num_embryos,hetero_count,homo_count)
    chosen_embryo_GRS[i] = mean(sort(GRS_scores)[1:embryo_cutoff])
  }
  return(chosen_embryo_GRS)
}
#returns a vector of mean GRS of embryos of a sample_size of simulated parental couples.
find_mean_GRS_2 <-function(log_OR, homo_count, hetero_count, num_embryos,sample_size) {
  mean_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    GRS_scores = model_all_embryo_GRS_2(log_OR,num_embryos,hetero_count, homo_count)
    mean_embryo_GRS[i] = mean(GRS_scores)
  }
  return(mean_embryo_GRS)
}
#returns a vector of the safest embryos GRS
find_healthy_embryo_GRS <-function(log_OR, homo_count, hetero_count, num_embryos, sample_size, embryo_cutoff = 3) {
  mean_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    GRS_scores = model_all_embryo_GRS_2(log_OR,num_embryos,hetero_count, homo_count)
    mean_embryo_GRS[i] = mean(sort(GRS_scores)[1:embryo_cutoff])
  }
  return(mean_embryo_GRS)
}
#find GRS improved for high risk
find_GRS_delta_improved_high_risk <-function(log_OR=alz_logOR_17, homo_count=alz_homo_17, hetero_count=alz_hetero_17, 
                                             num_embryos = 15, sample_size=100, embryo_cutoff = 3, high_risk_GRS_cutoff=0.74) {
  embryo_GRS_delta = array(0,sample_size)
  for(i in 1:sample_size) {
    high_risk_genotype = model_high_risk_genotype(log_OR, homo_count, hetero_count, high_risk_GRS_cutoff)
    normal_genotype = simulate_genotype_2(hetero_count, homo_count)
    GRS_scores = model_all_embryo_GRS_2(high_risk_genotype, normal_genotype, log_OR,num_embryos,hetero_count, homo_count)
    embryo_GRS_delta[i] = mean(GRS_scores)-mean(sort(GRS_scores)[1:embryo_cutoff])
  }
  return(exp(embryo_GRS_delta))
}
# find the GRS distribution expected given a log_OR, homo_count, hetero alele count, and sample size.
distribution_of_GRS_not_normalized = function(log_OR=alz_logOR_17, homo_count=alz_homo_17,
                                              hetero_count = alz_homo_17, sample_size=1000){
  GRS_array = array(0,sample_size)
  for(i in 1:sample_size){
    genotype = simulate_genotype_2(hetero_count, homo_count)
    GRS_array[i] = calculate_GRS_from_logOR(genotype, log_OR)
  }
  return(GRS_array)
}
model_high_risk_genotype <-function(log_OR = alz_logOR_17,homo_count = alz_homo_17,hetero_count = alz_hetero_17,
                                    high_risk_GRS_cutoff = 0.74) {
  high_risk_genotype = simulate_genotype_2(hetero_count, homo_count)
  high_risk_GRS = calculate_GRS_from_logOR(high_risk_genotype, log_OR)
  while(high_risk_GRS < high_risk_GRS_cutoff) {
    high_risk_genotype = simulate_genotype_2(hetero_count, homo_count)
    high_risk_GRS = calculate_GRS_from_logOR(high_risk_genotype, log_OR)
  }
  return(high_risk_genotype)
}
find_GRS_improved_risk_stratification <-function(log_OR=alz_logOR_17, homo_count=alz_homo_17, hetero_count=alz_hetero_17, 
                                                 num_embryos = 15, sample_size=100, embryo_cutoff = 15, high_risk_GRS_cutoff=1) {
  mean_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    high_risk_genotype = model_high_risk_genotype(log_OR, homo_count, hetero_count, high_risk_GRS_cutoff)
    normal_risk_genotype = simulate_genotype_2(hetero_count, homo_count)
    GRS_scores = model_all_embryo_GRS_2(high_risk_genotype, high_risk_genotype, log_OR,num_embryos,hetero_count, 
                                        homo_count)
    mean_embryo_GRS[i] = mean(sort(GRS_scores)[1:embryo_cutoff])
  }
  return(mean_embryo_GRS)
}



############################








#running models
hist(find_GRS_delta_improved_high_risk(sample_size = 1000, high_risk_GRS_cutoff = 1),breaks = seq(1,3,0.1),xlim=c(1,3), col = "blue")
hist(find_GRS_delta_improved_high_risk(sample_size = 1000, high_risk_GRS_cutoff = -2),breaks = seq(1,3,0.1),xlim=c(1,3),
     add = T)
