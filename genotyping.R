alz_eaf = c(0.25,0.20,0.41,0.49,0.59,0.28,0.27,0.63,0.71,0.66,0.37,0.62,0.32,0.6,0.64,0.96,0.09,0.78,0.19,0.92)
alz_OR = c(2.08,1.18,1.22,1.08,1.08,1.11,1.10,1.08,1.10,1.11,1.10,1.16,1.08,1.11,1.15,1.30,1.14,1.10,1.15,1.14)
alz_eaf_no_ApoE4 = c(0.20,0.41,0.49,0.59,0.28,0.27,0.63,0.71,0.66,0.37,0.62,0.32,0.6,0.64,0.96,0.09,0.78,0.19,0.92)
alz_OR_no_ApoE4 = c(1.18,1.22,1.08,1.08,1.11,1.10,1.08,1.10,1.11,1.10,1.16,1.08,1.11,1.15,1.30,1.14,1.10,1.15,1.14)
chd_eaf = c(0.33,0.08,0.80,0.39,0.08,0.11)
chd_OR = c(1.15,1.34,1.20,1.19,1.40,1.40)

#returns a simulated genotype using effect allele frequency data.
simulate_genotype <- function(SNP_eaf){
  SNPs_in_model = length(SNP_eaf)
  genotype = array(0,SNPs_in_model)
  for(SNP in 1:SNPs_in_model) {
    odds_of_allele = SNP_eaf[SNP] #p^2+2pq+q^2=1 so eaf=p^2+p(1-p)
    heterozygote_freq = 2*odds_of_allele*(1-odds_of_allele)
    homozygote_freq = heterozygote_freq^2
    no_allele = (1-odds_of_allele)^2
    zygote_distribution = c(no_allele,heterozygote_freq,homozygote_freq) #this isn't the most accurate way to go from eaf to hterozygote quantity
    genotype[SNP] = sample(0:2, 1,prob=zygote_distribution)
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
#calculates GRS score using alzheimers paper genotype info
calculate_GRS <-function(genotype,SNP_weights) {
  log_GRS=0
  for(SNP in 1:length(genotype)) {
    if(!(genotype[SNP] == 0))
    log_GRS = log_GRS+log(SNP_weights[SNP])
  }
  GRS = exp(log_GRS)
  return(GRS)
}
#outputs vector of GRS scores from maternal and paternal genotype information.
model_embryo_GRS <- function(mother_genotype,father_genotype, SNP_weights, num_embryos = 20) {
  GRS_scores = array(0,num_embryos)
  for(embryo in 1:num_embryos) {
    embryo_genotype = simulate_progeny(mother_genotype, father_genotype)
    GRS_scores[embryo] = calculate_GRS(embryo_genotype,SNP_weights)
  }
  return(GRS_scores)
}
#simulates a mother & father and outputs a vector of GRS scores of progeny embryos
model_all_embryo_GRS <-function(SNP_weights,num_embryos,SNP_eaf) {
  mother_genotype = simulate_genotype(SNP_eaf)
  father_genotype = simulate_genotype(SNP_eaf)
  GRS_scores = model_embryo_GRS(mother_genotype, father_genotype, SNP_weights, num_embryos)
  return(GRS_scores)
}
#returns number of risk alleles in a given genotype.
count_risk_alleles <- function(genotype) {
  return(sum(genotype))
}
count_no-zero_risk <- function(genotype) {
  return(sum(genotype !=0))
}
#returns a vector of the distribution of risk alleles in simulated embryos from a simulated couple.
model_embryo_risk_alleles <-function(num_embryos,SNP_eaf) {
  risk_allele_count = array(0,num_embryos)
  mother_genotype = simulate_genotype(SNP_eaf)
  father_genotype = simulate_genotype(SNP_eaf)
  for(embryo in 1:num_embryos) {
    embryo_genotype = simulate_progeny(mother_genotype, father_genotype)
    risk_allele_count[embryo] = count_risk_alleles(embryo_genotype)
  }
  return(risk_allele_count)
}
#returns a vector of the GRS of the x-healthiest embryo available to a sample_size of simulated parental couples.
find_healthiest_embryos <-function(SNP_weights, SNP_eaf, num_embryos,sample_size, embryo_cutoff) {
  chosen_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    GRS_scores = model_all_embryo_GRS(SNP_weights,num_embryos,SNP_eaf)
    chosen_embryo_GRS[i] = sort(GRS_scores)[embryo_cutoff]
  }
  return(chosen_embryo_GRS)
}
#returns a vector of mean GRS of embryos of a sample_size of simulated parental couples.
find_mean_GRS <-function(SNP_weights, SNP_eaf, num_embryos,sample_size) {
  mean_embryo_GRS = array(0,sample_size)
  for(i in 1:sample_size) {
    GRS_scores = model_all_embryo_GRS(SNP_weights,num_embryos,SNP_eaf)
    mean_embryo_GRS[i] = mean(GRS_scores)
  }
  return(mean_embryo_GRS)
}
#finds the number of alleles in the 3rd healthiest baby in a sample_size of simulated parental couples.
find_low_risk_allele_embryos <-function(SNP_eaf, num_embryos,sample_size, embryo_cutoff = 3) {
  chosen_embryo_risk_allele_distribution = array(0,sample_size)
  for(i in 1:sample_size) {
    risk_allele_distribution = model_embryo_risk_alleles(num_embryos,SNP_eaf)
    chosen_embryo_risk_allele_distribution[i] = sort(risk_allele_distribution)[embryo_cutoff]
  }
  return(chosen_embryo_risk_allele_distribution)
}
#finds the mean nuber of risk alleles in a sample_size of siulated parental couples.
find_mean_risk_allele_embryos <-function(SNP_eaf, num_embryos, sample_size) {
  mean_embryo_risk_allele_num = array(0,sample_size)
  for(i in 1:sample_size) {
    risk_allele_distribution = model_embryo_risk_alleles(num_embryos,SNP_eaf)
    mean_embryo_risk_allele_num[i] = mean(risk_allele_distribution)
  }
  return(mean_embryo_risk_allele_num)
}
#finds the average number of risk alleles reduced when choosing 3rd best embryo
find_risk_alleles_avoided <-function(SNP_eaf, num_embryos, sample_size, embryo_cutoff = 3) {
  risk_alleles_avoided = array(0,sample_size)
  for(i in 1:sample_size) {
    risk_allele_distribution = model_embryo_risk_alleles(num_embryos,SNP_eaf)
    risk_alleles_avoided[i] = (mean(risk_allele_distribution)-mean(sort(risk_allele_distribution)[1:embryo_cutoff]))
  }
  return(risk_alleles_avoided)
}

############################
#running models

all_embryo_GRS = model_all_embryo_GRS(alz_OR,1000,alz_eaf)
hist(all_embryo_GRS,xlim= c(1,20), breaks=1000)

#mean_embryo_GRS = find_mean_GRS(SNP_weights = alz_OR, SNP_eaf = alz_eaf, num_embryos = 20, sample_size = 1000)
#pop_mean_GRS = mean(mean_embryo_GRS)
#mean_embryo_GRS= mean_embryo_GRS_2-pop_mean_GRS
#hist(mean_embryo_GRS, ylim = c(0,100), breaks = 20, col = "red")

#chosen_embryo_GRS = find_healthiest_embryos(SNP_weights = alz_OR, SNP_eaf = alz_eaf, num_embryos = 20, sample_size = 100, embryo_cutoff = 3)
#chosen_embryo_GRS = chosen_embryo_GRS - pop_mean_GRS
#hist(chosen_embryo_GRS, add = T, breaks = 20, col = "blue")

#mean_embryo_risk_alleles = find_mean_risk_allele_embryos(alz_eaf,num_embryos = 20,sample_size = 10000)
#hist(mean_embryo_risk_alleles,xlim=c(0,20),breaks=20)

#low_risk_allele_num_embryos = find_low_risk_allele_embryos(alz_eaf,num_embryos = 20, sample_size = 67000, embryo_cutoff = 3)
#hist(low_risk_allele_num_embryos,xlim=c(0,20),breaks=20)

#alz_risk_alleles_avoided = find_risk_alleles_avoided(alz_eaf, num_embryos = 20, sample_size = 1000, embryo_cutoff = 3)
#hist(risk_alleles_avoided,xlim=c(0,8), col="blue")

#chd_risk_alleles_avoided = find_risk_alleles_avoided(chd_eaf, num_embryos = 20, sample_size = 1000, embryo_cutoff = 3)
#hist(chd_risk_alleles_avoided,xlim=c(0,8), col="red")

