############ ICD-10 extraction of first diagnosis PD
{
  library(data.table)
  first_diagnosis<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv",select=c('eid','53-0.0',
                                                                                                              '42030-0.0','42031-0.0','42032-0.0','42033-0.0'),
                         header = T,sep=",")
  first_diagnosis=data.frame(first_diagnosis)
  diagnosis_parkinson<-data.frame(id=first_diagnosis$eid,diagnosis_date=as.Date(NA),source=NA,time_interval=NA,attending_date=first_diagnosis$X53.0.0)
  
  diagnosis_parkinson$diagnosis_date=pmin(first_diagnosis$X42030.0.0,first_diagnosis$X42032.0.0,na.rm = TRUE)
  
  index_0=which(first_diagnosis$X42030.0.0==diagnosis_parkinson$diagnosis_date)
  diagnosis_parkinson$source[index_0]=first_diagnosis$X42031.0.0[index_0]
  diagnosis_parkinson$time_interval=diagnosis_parkinson$diagnosis_date-diagnosis_parkinson$attending_date
  rm(first_diagnosis,index_0)
  
  index_nonNan=which(is.na(diagnosis_parkinson$diagnosis_date)==FALSE)
  diagnosis_parkinson=diagnosis_parkinson[index_nonNan,]
  rm(index_nonNan)
  
  
  ##### max date of icd record (use to determine the censored date)
  censoring_date=max(diagnosis_parkinson$diagnosis_date,na.rm = TRUE)
  
  
  ######## lost of follow-up
  library(data.table)
  datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv", 
               select=c('eid','190-0.0','191-0.0','53-0.0'),header = T,sep=",")
  lost_follow_up=data.frame(datas)
  colnames(lost_follow_up)=c('id','lost_reason','lost_date','attending_date')
  rm(datas)
  index_noNaN<-which(is.na(lost_follow_up$lost_date)==FALSE)
  lost_follow_up<-lost_follow_up[index_noNaN,]
  lost_follow_up$lost_time_interval<-lost_follow_up$lost_date-lost_follow_up$attending_date
  rm(index_noNaN)
  
  ###death information
  library(data.table)
  datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv", select=
                 c('eid','40000-0.0','53-0.0'), header = T,sep=",")
  colnames(datas)=c('id','death_date','attending_date')
  index_death<-which(!is.na(datas$death_date)==TRUE)
  death_information<-datas[index_death,]
  rm(datas)
  rm(index_death)
  death_information$time_interval=death_information$death_date-death_information$attending_date
  
  
  ################ censoring data, days between the date attending the assessment center and the date recording information
  library(data.table)
  datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv", select=c('eid','53-0.0'),header = T,sep=",")
  datas=data.frame(datas)
  data_censoring=data.frame(id=datas$eid,attending_date=datas$X53.0.0)
  data_censoring$censoring=rep(as.IDate(censoring_date),nrow(data_censoring))-data_censoring$attending_date
  rm(censoring_date)
  rm(datas)
  
  #incorporate lost information into the data_censoring table
  data_censoring$lost=NA
  data_censoring$lost_reason=NA
  lost_index=match(lost_follow_up$id,data_censoring$id)
  data_censoring$lost[lost_index]=lost_follow_up$lost_time_interval
  data_censoring$lost_reason[lost_index]=lost_follow_up$lost_reason
  rm(lost_index,lost_follow_up)
  
  #incorporate death information into the data_censoring table
  data_censoring$death=NA
  death_index=match(death_information$id,data_censoring$id)
  data_censoring$death[death_index]=death_information$time_interval
  rm(death_index)
  rm(death_information)
  
  #incorporate information of diagnosis into the data_censoring table
  data_censoring$diagnosis=NA
  data_censoring$source=NA
  diagnosis_index=match(diagnosis_parkinson$id,data_censoring$id)
  data_censoring$diagnosis[diagnosis_index]=diagnosis_parkinson$time_interval
  data_censoring$source[diagnosis_index]=diagnosis_parkinson$source
  rm(diagnosis_index,diagnosis_parkinson)
  
  
  ######calculate censoring days and label, can modify when need to record the participants with diagnosis before enrollment
  #diagnosis
  data_censoring$label=0
  index_diagnosis=which(is.na(data_censoring$diagnosis)==FALSE)
  data_censoring$label[index_diagnosis]=1
  data_censoring$censoring[index_diagnosis]=data_censoring$diagnosis[index_diagnosis]
  rm(index_diagnosis)
  
  #lost follow-up
  index_lost_wo_diagnosis=which(is.na(data_censoring$lost)==FALSE & data_censoring$label==0)#being lost after have a diagnosis
  data_censoring$censoring[index_lost_wo_diagnosis]=data_censoring$lost[index_lost_wo_diagnosis]
  rm(index_lost_wo_diagnosis)
  
  #death
  index_death_wo_diagnosis=which(is.na(data_censoring$death)==FALSE & data_censoring$label==0)
  data_censoring$censoring[index_death_wo_diagnosis]=data_censoring$death[index_death_wo_diagnosis]
  data_censoring$label[index_death_wo_diagnosis]=0
  rm(index_death_wo_diagnosis)
}

######################pain locations
{
  library(data.table)
  datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv",header = T,sep=",")
  datas=data.frame(datas)
  
  temp_pain=datas$X6159.0.0
  index_no_pain=which(temp_pain==-7)
  index_not_answer=which(temp_pain==-3)
  index_pain=which(temp_pain>0)
  temp_pain[index_no_pain]=0
  temp_pain[index_pain]=NA
  temp_pain[index_not_answer]=NA
  rm(index_no_pain,index_pain,index_not_answer)
  
  pain_duration=data.frame(headache=datas$X3799.0.0,
                           facial=datas$X4067.0.0,
                           neck_shoulder=datas$X3404.0.0,
                           back=datas$X3571.0.0,
                           stomach=datas$X3741.0.0,
                           hip_pain=datas$X3414.0.0,
                           knee=datas$X3773.0.0,
                           body=datas$X2956.0.0)
  
  for (i in 1:ncol(pain_duration)){
    chronic_pain=which((pain_duration[,i]==1))
    short_pain=which((pain_duration[,i]==0))
    not_answer=which((pain_duration[,i]<0))
    index_no_pain=which(temp_pain==0)
    pain_duration[chronic_pain,i]=2
    pain_duration[short_pain,i]=1
    pain_duration[not_answer,i]=NA
    pain_duration[index_no_pain,i]=0
    rm(chronic_pain,short_pain,not_answer,index_no_pain)
  }
  rm(i,temp_pain)
}

############################covariate variables######
{
  datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/Pain_parkinson.csv",header = T,sep=",")
  datas=data.frame(datas)
  age=datas$X21003.0.0
  gender=datas$X31.0.0
  # genetic_sex=datas$X22001.0.0
  # index_unmatch=which(genetic_sex!=gender)
  # gender[index_unmatch]=NA
  # rm(index_unmatch,genetic_sex)
  
  #ethnicity
  ethnicity=cbind(datas$X21000.0.0,datas$X21000.1.0,datas$X21000.2.0)
  index_nan=which(ethnicity[,1]==-1|ethnicity[,1]==-3)
  ethnicity[index_nan,1]=pmax(ethnicity[index_nan,1],ethnicity[index_nan,2],ethnicity[index_nan,3],na.rm = TRUE)
  race=ethnicity[,1]
  rm(ethnicity,index_nan)
  index_white=which((race %in% c(1,1001,1002,1003))==TRUE)
  index_nonwhite=which((race %in% c(1,1001,1002,1003))==FALSE)
  index_not_answer=which(race==-1|race==-3)
  race[index_white]=1
  race[index_nonwhite]=0
  race[index_not_answer]=NA
  rm(index_white,index_nonwhite,index_not_answer)
  ethnicity=race
  rm(race)
  
  sites=datas$X54.0.0
  
  deprivation=datas$X22189.0.0
  tertile_point=quantile(deprivation, probs = seq(0, 1, 1/3),na.rm = TRUE)
  index_low=which(deprivation<tertile_point[2])
  index_middle=which(deprivation>=tertile_point[2]&deprivation<tertile_point[3])
  index_high=which(deprivation>=tertile_point[3])
  deprivation[index_low]=0
  deprivation[index_middle]=1
  deprivation[index_high]=2
  rm(index_low,index_middle,index_high,tertile_point)
  
  
  whr=datas$X48.0.0/datas$X49.0.0
  bmi=datas$X21001.0.0

  
  income=datas$X738.0.0
  index_not_answe=which(income==-3)
  index_non=which(income==-1)
  index_low=which(income<=3&income>=1)
  index_medium=which(income==4)
  index_high=which(income==5)
  income[index_non]=0
  income[index_low]=1
  income[index_medium]=2
  income[index_high]=3
  income[index_not_answe]=0
  rm(index_non,index_low,index_medium,index_high,index_not_answe)
  
  
  smoking=datas$X20116.0.0
  index_non=which(is.na(smoking)==TRUE)
  smoking[index_non]=NA
  index_no_answer=which(smoking==-3)
  smoking[index_no_answer]=NA
  index_current=which(smoking==2)
  smoking[index_current]=1
  rm(index_non,index_no_answer,index_current)
  
  alcohol=datas$X20117.0.0
  index_no_answer=which(alcohol<0)
  alcohol[index_no_answer]=NA
  index_current=which(alcohol==2)
  alcohol[index_current]=1
  rm(index_no_answer,index_current)
  
  education=cbind(datas$X6138.0.0,datas$X6138.0.1,datas$X6138.0.2,datas$X6138.0.3,datas$X6138.0.4,datas$X6138.0.5)
  for (i in 1:6){
    index_no_answer=which(education[,i]==-3)
    education[index_no_answer,i]=NA
    index_none=which(education[,i]==-7)
    education[index_none,i]=7
    
    index_20=which(education[,i]==1)
    index_13=which(education[,i]==2)
    index_10=which(education[,i]==3|education[,i]==4)
    index_19=which(education[,i]==5)
    index_15=which(education[,i]==6)
    index_7=which(education[,i]==7)
    education[index_20,i]=20
    education[index_13,i]=13
    education[index_10,i]=10
    education[index_19,i]=19
    education[index_15,i]=15
    education[index_7,i]=7
    rm(index_20,index_13,index_10,index_19,index_15,index_7,index_none,index_no_answer)
  }
  rm(i)
  university_education=pmax(education[,1],education[,2],education[,3],education[,4],education[,5],education[,6],na.rm = TRUE)
  index_no_education=which(university_education==7)
  index_no_university=which(university_education<20&university_education!=7)
  index_university=which(university_education==20)
  university_education[index_no_education]=0
  university_education[index_no_university]=1
  university_education[index_university]=2
  rm(index_no_university,index_university,education,index_no_education)
  
  
  #sleep
  sleep=datas$X1160.0.0
  index_no_answer=which(sleep<0)
  index_excess=which(sleep>9)
  index_lack=which(sleep<7)
  index_normal=which(sleep>=7&sleep<=9)
  sleep[index_no_answer]=NA
  sleep[index_normal]=0
  sleep[index_excess]=1
  sleep[index_lack]=1
  rm(index_no_answer,index_normal,index_excess,index_lack)
  
  #watching tv
  sedentary=datas$X1070.0.0
  index_na=which(sedentary==-1|sedentary==-3)
  index_low=which((sedentary<4&sedentary>=0)|sedentary==-10)
  index_high=which(sedentary>=4)
  sedentary[index_na]=NA
  sedentary[index_low]=0
  sedentary[index_high]=1
  rm(index_na,index_low,index_high)
  
  #diabetes
  diabetes=datas$X2443.0.0
  glucose=datas$X30740.0.0
  index_diabetes=which(diabetes==1|glucose>=5.6)
  index_nan=which(diabetes==-1|diabetes==-3)
  diabetes[index_nan]=NA
  rm(index_nan)
  index_na=which((is.na(diabetes)==TRUE)&(is.na(glucose)==TRUE))
  diabetes=rep(0,nrow(datas))
  diabetes[index_diabetes]=1
  diabetes[index_na]=NA
  rm(glucose,index_diabetes,index_na)
  
  #obesity
  obesity=rep(NA,nrow(datas))
  index_underweight=which(bmi<18.5)
  index_normal=which(bmi<25&bmi>=18.5)
  index_overweight=which(bmi>=25&bmi<30)
  index_obese1=which(bmi>=30&bmi<40)
  index_obese2=which(bmi>=40)
  obesity[index_underweight]=0
  obesity[index_normal]=0
  obesity[index_overweight]=1
  obesity[index_obese1]=2
  obesity[index_obese2]=3#extreme obesity
  rm(index_underweight,index_normal,index_overweight,index_obese1,index_obese2)
  
  #cancer history
  cancer=datas$X2453.0.0
  index_not_answer=which(cancer<0)
  cancer[index_not_answer]=NA
  rm(index_not_answer)
  
  #vascular history
  vascular=cbind(datas$X6150.0.0,datas$X6150.0.1,datas$X6150.0.2,datas$X6150.0.3)
  for (i in 1:4){
    index_no_answer=which(vascular[,i]==-3)
    vascular[index_no_answer,i]=NA
    rm(index_no_answer)
  }
  rm(i)
  vascular_history=pmin(vascular[,1],vascular[,2],vascular[,3],vascular[,4],na.rm = TRUE)
  index_no_vascular=which(vascular_history==-7)
  index_vascular=which(vascular_history>0)
  vascular_history[index_no_vascular]=0
  vascular_history[index_vascular]=1
  rm(vascular,index_vascular,index_no_vascular)

  #prs
  pd_prs=datas$X26260.0.0
  tertile_point=quantile(pd_prs, probs = seq(0, 1, 1/3),na.rm = TRUE)
  index_low=which(pd_prs<tertile_point[2])
  index_middle=which(pd_prs>=tertile_point[2]&pd_prs<tertile_point[3])
  index_high=which(pd_prs>=tertile_point[3])
  pd_prs[index_low]=0
  pd_prs[index_middle]=1
  pd_prs[index_high]=2
  rm(index_low,index_middle,index_high,tertile_point)
  
  comp_1=datas$X22009.0.1
  comp_2=datas$X22009.0.2
  comp_3=datas$X22009.0.3
  comp_4=datas$X22009.0.4
  comp_5=datas$X22009.0.5
  comp_6=datas$X22009.0.6
  comp_7=datas$X22009.0.7
  comp_8=datas$X22009.0.8
  comp_9=datas$X22009.0.9
  comp_10=datas$X22009.0.10
  
  ###############pain medication###############
  medication=cbind(datas$X6154.0.0,datas$X6154.0.1,datas$X6154.0.2,datas$X6154.0.3,datas$X6154.0.4,datas$X6154.0.5)
  for (i in 1:6){
    index_no_answer=which(medication[,i]==-3|medication[,i]==-1)
    medication[index_no_answer,i]=NA
    index_pain_relief=which(medication[,i]==1|medication[,i]==2|medication[,i]==3)
    medication[index_pain_relief,i]=1
    index_other_medication=which(medication[,i]==4|medication[,i]==5|medication[,i]==6)
    medication[index_other_medication,i]=2
    index_no_medication=which(medication[,i]==-7)
    medication[index_no_medication,i]=0
    rm(index_no_answer,index_pain_relief,index_other_medication,index_no_medication)
  }
  rm(i)
  pain_medication_max=pmax(medication[,1],medication[,2],medication[,3],medication[,4],medication[,5],medication[,6],na.rm = TRUE)
  pain_medication_min=pmin(medication[,1],medication[,2],medication[,3],medication[,4],medication[,5],medication[,6],na.rm = TRUE)
  index_no_medication=which(pain_medication_max==0)
  index_pain_medication=which(pain_medication_max==1)
  index_other_medication=which(pain_medication_min==2)
  index_2_medication=which(pain_medication_min==1&pain_medication_max==2)
  pain_medication=rep(NA,nrow(datas))
  pain_medication[index_no_medication]=0
  pain_medication[index_pain_medication]=1
  pain_medication[index_other_medication]=0
  pain_medication[index_2_medication]=1
  rm(pain_medication_max,pain_medication_min,index_no_medication,index_pain_medication,index_other_medication,index_2_medication)
  rm(datas,medication)
}

##############################combine data###########################################
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,
                age=age,gender=gender,sites=sites,ethnicity=ethnicity,
                university_education=university_education,deprivation=deprivation,income=income,
                sedentary=sedentary,smoking=smoking,alcohol=alcohol,obesity=obesity,prs=pd_prs,
                comp_1=comp_1,comp_2=comp_2,comp_3=comp_3,comp_4=comp_4,comp_5=comp_5,
                comp_6=comp_6,comp_7=comp_7,comp_8=comp_8,comp_9=comp_9,comp_10=comp_10)

data$ethnicity=factor(data$ethnicity)
data$university_education=factor(data$university_education)
data$gender=factor(data$gender)
data$income=factor(data$income)
data$alcohol=factor(data$alcohol)
data$smoking=factor(data$smoking)
data$obesity=factor(data$obesity)
data$deprivation=factor(data$deprivation)
data$sedentary=factor(data$sedentary)
data$prs=factor(data$prs)

rm(age, gender,alcohol,bmi,deprivation,diabetes,ethnicity,income,sleep,sedentary,
   sites,smoking,university_education,whr,data_censoring,obesity,cancer,vascular_history,pain_medication,
   pd_prs,comp_1,comp_2,comp_3,comp_4,comp_5,comp_6,comp_7,comp_8,comp_9,comp_10)
####################################################################################

##############################combine data, extended covariates###########################################
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,
                age=age,gender=gender,sites=sites,ethnicity=ethnicity,
                university_education=university_education,deprivation=deprivation,income=income,
                sedentary=sedentary,smoking=smoking,alcohol=alcohol,obesity=obesity,
                vascular_history=vascular_history,cancer=cancer,diabetes=diabetes,
                sleep=sleep,prs=pd_prs,medication=pain_medication,
                comp_1=comp_1,comp_2=comp_2,comp_3=comp_3,comp_4=comp_4,comp_5=comp_5,
                comp_6=comp_6,comp_7=comp_7,comp_8=comp_8,comp_9=comp_9,comp_10=comp_10)
rm(age, gender,alcohol,bmi,deprivation,diabetes,ethnicity,income,sleep,sedentary,
   sites,smoking,university_education,whr,data_censoring,obesity,cancer,vascular_history,pain_medication,
   pd_prs,comp_1=comp_1,comp_2=comp_2,comp_3=comp_3,comp_4=comp_4,comp_5=comp_5,
   comp_6=comp_6,comp_7=comp_7,comp_8=comp_8,comp_9=comp_9,comp_10=comp_10)


data$ethnicity=factor(data$ethnicity)
data$university_education=factor(data$university_education)
data$gender=factor(data$gender)
data$income=factor(data$income)
data$alcohol=factor(data$alcohol)
data$smoking=factor(data$smoking)
data$obesity=factor(data$obesity)
data$deprivation=factor(data$deprivation)
data$sedentary=factor(data$sedentary)
data$vascular_history=factor(data$vascular_history)
data$cancer=factor(data$cancer)
data$diabetes=factor(data$diabetes)
data$sleep=factor(data$sleep)
data$medication=factor(data$medication)
data$prs=factor(data$prs)

####################################################################################


##calculate the number of chronic pain sites
data$chronic_pain=rowSums(pain_duration[,1:7]==2,na.rm = TRUE)#the number of chronic chronic_pain sites (may also have acute pain)
index_more_5=which(data$chronic_pain>=5)#merge sites>=5 as 5
data$chronic_pain[index_more_5]=5
rm(index_more_5)
pain_sum=rowSums(pain_duration[,1:7],na.rm = TRUE)
index_short=which(data$chronic_pain==0&pain_sum!=0)#only acute pain
data$chronic_pain[index_short]=0.5
rm(index_short,pain_sum)
index_na=which(rowSums(is.na(pain_duration))==8)### NA values, NOo value for any of the 8 sites
data$chronic_pain[index_na]=NA
rm(index_na)

##calculate the number of acute pain sites
data$acute_pain=rowSums(pain_duration[,1:7]==1,na.rm = TRUE)
index_more_5=which(data$acute_pain>=5)#merge sites>=5 as 5
data$acute_pain[index_more_5]=5
rm(index_more_5)
pain_sum=rowSums(pain_duration[,1:7],na.rm = TRUE)
index_chronic=which(data$acute_pain==0&pain_sum!=0)#only chronic pain
data$acute_pain[index_chronic]=0.5
rm(index_chronic,pain_sum)
index_na=which(rowSums(is.na(pain_duration))==8)### NA values, NOo value for any of the 8 sites
data$acute_pain[index_na]=NA
rm(index_na)

######################################pain all over the body
index_body_pain=which(pain_duration$body>0)
data=data[-index_body_pain,]##remove pain all over the body
pain_duration=pain_duration[-index_body_pain,1:7]
rm(index_body_pain)
######################################do not run this code when running the analysis for each single pain site

############# retain complete data
number_na<-rowSums(is.na(data))
index_nonNA<-which(number_na==0)
datas<-data[index_nonNA,]
pain_duration=pain_duration[index_nonNA,]
#inflammation_biomarker=inflammation_biomarker[index_nonNA,]
rm(number_na,index_nonNA)

#can exclude participants younger than 40 years
index_included=which(datas$age>=40)
datas=datas[index_included,]
pain_duration=pain_duration[index_included,]
#inflammation_biomarker=inflammation_biomarker[index_included,]
datas$sites=factor(datas$sites)
rm(index_included)

index_2year_lanmark=which(!(datas$time<365*2&datas$status==1))
datas=datas[index_2year_lanmark,]
pain_duration=pain_duration[index_2year_lanmark,]
#inflammation_biomarker=inflammation_biomarker[index_2year_lanmark,]
rm(index_2year_lanmark)

##########pain category, single-site pain, multiple-sizte pain, short-term pain
datas$pain=datas$chronic_pain
index_multiple=which(datas$pain>=2)
datas$pain[index_multiple]=2
rm(index_multiple)

datas$pain=factor(datas$pain)
datas$pain = relevel(datas$pain, ref = "0")



##############################################models###############################################
library("survival")
library("survminer")
#extended covariates
hr <- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+
              vascular_history+cancer+diabetes+medication+sleep+
              prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
summary(hr)
test.zph <- cox.zph(hr,'rank')
test.zph

#covariates
hr <- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+
              prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
summary(hr)


###########################################HR for each pain site #######################################################
library("survival")
library("survminer")
result=data.frame(matrix(nrow=2*ncol(pain_duration),ncol=8))
row.names(result)=c(paste(colnames(pain_duration),'1'),paste(colnames(pain_duration),'2'))
colnames(result)=c('HR','lower','upper','CI','P_value','Participants','Events','beta')

for (index_pain in c(1:ncol(pain_duration))) {
  datas$pain_site=pain_duration[,index_pain]
  datas$pain_site=factor(datas$pain_site)
  hr <- coxph(Surv(time, status) ~ pain_site+age+gender+ethnicity+deprivation+income+university_education+
                smoking+alcohol+sedentary+obesity+
                #vascular_history+cancer+diabetes+medication+sleep+
                prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
  
  ##pain_1, short-term pain
  coeff<-summary(hr)
  result$HR[index_pain]=round(coeff$coefficients['pain_site1','exp(coef)'],2)
  result$P_value[index_pain]=signif(coeff$coefficients['pain_site1','Pr(>|z|)'],3)
  lower=round(coeff$conf.int['pain_site1','lower .95'],2)
  upper=round(coeff$conf.int['pain_site1','upper .95'],2)
  result$lower[index_pain]=lower
  result$upper[index_pain]=upper
  result$CI[index_pain]=paste(result$HR[index_pain],' (',lower,'-',upper,')',sep='')
  result$Participants[index_pain]=length(which(datas$pain_site==1))
  result$Events[index_pain]=length(which(datas$pain_site==1&(datas$status==1)))
  result$beta[index_pain]=signif(coeff$coefficients['pain_site1','coef'],3)
  rm(lower,upper)
  print(c(index_pain,length(which(datas$pain_site==1)),length(which(datas$pain_site==1&(datas$status==1)))))
  ##pain_2, chronic pain
  result$HR[index_pain+7]=round(coeff$coefficients['pain_site2','exp(coef)'],2)
  result$P_value[index_pain+7]=signif(coeff$coefficients['pain_site2','Pr(>|z|)'],3)
  lower=round(coeff$conf.int['pain_site2','lower .95'],2)
  upper=round(coeff$conf.int['pain_site2','upper .95'],2)
  result$lower[index_pain+7]=lower
  result$upper[index_pain+7]=upper
  result$CI[index_pain+7]=paste(result$HR[index_pain+7],' (',lower,'-',upper,')',sep='')
  result$Participants[index_pain+7]=length(which(datas$pain_site==2))
  result$Events[index_pain+7]=length(which(datas$pain_site==2&(datas$status==1)))
  result$beta[index_pain+7]=signif(coeff$coefficients['pain_site2','coef'],3)
  
  print(c(index_pain,length(which(datas$pain_site==2)),length(which(datas$pain_site==2&(datas$status==1)))))
  rm(hr,lower,upper,coeff)
}
rm(index_pain)
result=result[c(1,8,2,9,3,10,4,11,5,12,6,13,7,14),]
PF_events=length(which(datas$pain_site==0&datas$status==1))
for (i in c(1:nrow(result))) {
  result$paf[i]=round(result$Events[i]/(PF_events+result$Events[i])*(as.numeric(result$HR[i])-1)/as.numeric(result$HR[i]),4)*100
}
rm(PF_events)
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_pain_locations.csv')


####################################HR 0/0.5/1/2,write results##################################################
hr <- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+
              smoking+alcohol+sedentary+obesity+
              prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,
              data = datas)
coeff<-summary(hr)
result=data.frame(HR=round(coeff$coefficients[,2],2),lower=round(coeff$conf.int[,3],2),upper=round(coeff$conf.int[,4],2),p_value=signif(coeff$coefficients[,5],3))
for (i in 1:nrow(result)){
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
rm(i)
{
  result<-rbind(c(1, 1, 1, 1, '1.00 (Ref)'),
               result[1:4,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[5,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[6,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[7:8,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[9:11,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[12:13,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[14,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[15,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[16,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[17:19,],c(1, 1, 1, NA, '1.00 (Ref)'),
               result[20:21,])
  row.names(result)[c(1,6,8,10,13,17,20,22,24,26,30)]=c('Pain free','Female','Ethnic minorities','Deprivation_low','Income_Unknown',
                                                              'No education','Never_smoking','Never drinking','0-4 hours/day','Healthy weight','Low genetic risk')
  result$Participants=NA
  result$Events=NA
}#pain
unique_pain=c(0,0.5,1,2)
for (i in 1:4){
  result$Participants[i]=length(which(datas$pain==unique_pain[i]))
  result$Events[i]=length(which(datas$pain==unique_pain[i]&(datas$status==1)))
}
rm(unique_pain)
{
  #age
    result$Participants[5]=length(datas$age)
    result$Events[5]=length(which(datas$status==1))
  #gender
  for (i in 6:7){
    result$Participants[i]=length(which(datas$gender==i-6))
    result$Events[i]=length(which(datas$gender==(i-6)&datas$status==1))
  }
  #ethnicity
  for (i in 8:9){
    result$Participants[i]=length(which(datas$ethnicity==i-8))
    result$Events[i]=length(which(datas$ethnicity==(i-8)&datas$status==1))
  }
  #deprivation
  for (i in 10:12){
    result$Participants[i]=length(which(datas$deprivation==i-10))
    result$Events[i]=length(which(datas$deprivation==(i-10)&datas$status==1))
  }
  #income
  for (i in 13:16){
    result$Participants[i]=length(which(datas$income==i-13))
    result$Events[i]=length(which(datas$income==(i-13)&datas$status==1))
  }
  #education
  for (i in 17:19){
    result$Participants[i]=length(which(datas$university_education==i-17))
    result$Events[i]=length(which(datas$university_education==(i-17)&datas$status==1))
  }
  #smoking
  for (i in 20:21){
    result$Participants[i]=length(which(datas$smoking==i-20))
    result$Events[i]=length(which(datas$smoking==(i-20)&datas$status==1))
  }
  #alcohol
  for (i in 22:23){
    result$Participants[i]=length(which(datas$alcohol==i-22))
    result$Events[i]=length(which(datas$alcohol==(i-22)&datas$status==1))
  }
  #sedentary
  for (i in 24:25){
    result$Participants[i]=length(which(datas$sedentary==i-24))
    result$Events[i]=length(which(datas$sedentary==(i-24)&datas$status==1))
  }
  #Obesity
  for (i in 26:29){
    result$Participants[i]=length(which(datas$obesity==i-26))
    result$Events[i]=length(which(datas$obesity==(i-26)&datas$status==1))
  }
  #genetic risk
  for (i in 30:32){
    result$Participants[i]=length(which(datas$prs==i-30))
    result$Events[i]=length(which(datas$obesity==(i-30)&datas$status==1))
  }
  rm(coeff,i)
}
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status.csv')
paf_short_pain=round(result$Events[2]/(result$Events[1]+result$Events[2])*(as.numeric(result$HR[2])-1)/as.numeric(result$HR[2]),4)*100
paf_single_pain=round(result$Events[3]/(result$Events[1]+result$Events[3])*(as.numeric(result$HR[3])-1)/as.numeric(result$HR[3]),4)*100
paf_multiple_pain=round(result$Events[4]/(result$Events[1]+result$Events[4])*(as.numeric(result$HR[4])-1)/as.numeric(result$HR[4]),4)*100


####################################nonlinear between No and HR####################
library(survival)
library(rms)
library(dplyr)
library(Gmisc)
library(splines)
library(Greg)
library(tidyverse)
#0,1,2,3,4,5,6,7
index=which(datas$pain==0.5)
datas=datas[-index,]
rm(index)
dd <- datadist(datas) #为后续程序设定数据环境
options(datadist='dd') #为后续程序设定数据环境

fit<- cph(Surv(time, status) ~ rcs(pain,3) +age+gender+ethnicity+deprivation+income+university_education+
            smoking+alcohol+sedentary+obesity+
            prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
anova(fit)#test the significance of nolinear regression

plotHR(fit, term="pain", xlab="No of acute pain sites",se=TRUE,
       polygon_ci=TRUE,
       col.term = "#0070b9",
       col.se = "#DEEBF7BB",
       lwd.term=2,
       lty.term = 1,
       ylog=FALSE,
       alpha=0.05,#alpha level, 95%
       cex=1,
       axes = TRUE,
       plot.bty="l", xlim=c(0,7),ylim=c(0,2.5),rug="density")




##########################################0/1/2,subgroup analysis by all covariates#########################################################
result<-data.frame(id=seq(1:37),
                   HR_single=NA,P_single=NA,CI_single=NA,lower_single=NA,upper_single=NA,
                   HR_multiple=NA,P_multiple=NA,CI_multiple=NA,lower_multiple=NA,upper_multiple=NA)
row.names(result)=c('Sex','Females','Males','Race','Ethnic minorities','White','Deprivation','Lower','Middle','Higher',
                    'Income','Unknown','Low (<£51999)','Middle (£52k-£100k)','High (>£100k)','Education','No education','Less than college','Above College',
                    'Smoking status','Never','Ever','Alcohol intake status','Never_0','Ever_1',
                    'Sedentary behavior','0-4 hours/day','>4 hours/day',
                    'BMI category','Healthy weight','Overweight','Obese','Severe Obese','Genetic risk','Low','Intermediate','High')
index_short=which(datas$pain==0.5)
datas=datas[-index_short,]
rm(index_short)
datas$pain=as.numeric(as.character(datas$pain))

for (pain_status in c(1:2)) {
  
  temp_data <- subset(datas, pain==0|pain==pain_status)
  temp_data$pain=factor(temp_data$pain)

  ##sex
  fit1<- coxph(Surv(time, status) ~ pain:gender+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1,2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+1,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+1,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+1,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+1,5*pain_status]=lower
    result[i+1,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##race
  fit1<- coxph(Surv(time, status) ~ pain:ethnicity+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1,2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+4,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+4,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+4,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+4,5*pain_status]=lower
    result[i+4,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##deprivation
  fit1<- coxph(Surv(time, status) ~ pain:deprivation+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+7,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+7,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+7,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+7,5*pain_status]=lower
    result[i+7,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##income
  fit1<- coxph(Surv(time, status) ~ pain:income+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:4)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+11,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-4+i,'exp(coef)']),2)
    result[i+11,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-4+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-4+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-4+i,'upper .95'],2)
    result[i+11,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+11,5*pain_status]=lower
    result[i+11,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##education
  fit1<- coxph(Surv(time, status) ~ pain:university_education+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+16,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+16,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+16,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+16,5*pain_status]=lower
    result[i+16,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##smoking
  fit1<- coxph(Surv(time, status) ~ pain:smoking+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+20,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+20,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+20,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+20,5*pain_status]=lower
    result[i+20,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##alcohol
  fit1<- coxph(Surv(time, status) ~ pain:alcohol+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+23,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+23,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+23,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+23,5*pain_status]=lower
    result[i+23,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##sedentary
  fit1<- coxph(Surv(time, status) ~ pain:sedentary+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+26,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+26,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+26,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+26,5*pain_status]=lower
    result[i+26,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##bmi category
  fit1<- coxph(Surv(time, status) ~ pain:obesity+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:4)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+29,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-4+i,'exp(coef)']),2)
    result[i+29,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-4+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-4+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-4+i,'upper .95'],2)
    result[i+29,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+29,5*pain_status]=lower
    result[i+29,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##genetic risk
  fit1<- coxph(Surv(time, status) ~ pain:prs+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+34,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+34,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+34,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+34,5*pain_status]=lower
    result[i+34,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  
  rm(temp_data)
}
result$CI_single=paste(result$HR_single,result$CI_single,sep=' ')
result$CI_multiple=paste(result$HR_multiple,result$CI_multiple,sep=' ')
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/subgroup_hr_status_chronic.csv')


###############################################interaction at the whole group level#########################################
index=which(datas$pain==0.5)
datas=datas[-index,]
rm(index)
fit0<- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit1<- coxph(Surv(time, status) ~ pain*age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit2<- coxph(Surv(time, status) ~ pain*gender+age+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit3<- coxph(Surv(time, status) ~ pain*ethnicity+age+gender+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit4<- coxph(Surv(time, status) ~ pain*deprivation+age+gender+ethnicity+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit5<- coxph(Surv(time, status) ~ pain*income+age+gender+ethnicity+deprivation+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit6<- coxph(Surv(time, status) ~ pain*university_education+age+gender+ethnicity+deprivation+income+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit7<- coxph(Surv(time, status) ~ pain*smoking+age+gender+ethnicity+deprivation+income+university_education+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit8<- coxph(Surv(time, status) ~ pain*alcohol+age+gender+ethnicity+deprivation+income+university_education+smoking+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit9<- coxph(Surv(time, status) ~ pain*sedentary+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit10<-coxph(Surv(time, status) ~ pain*obesity+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit11<-coxph(Surv(time, status) ~ pain*prs+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
fit_models=list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
rm(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
P_interaction=data.frame(matrix(ncol = 1, nrow = 11))
row.names(P_interaction)=c('Age','Sex','Race','Deprivation','Income','University_education','Smoking','Alcohol','Sedentary','Obesity','Genetic risk')
colnames(P_interaction)=c('P_interaction')
for (i in c(1:11)) {
  P_interaction$P_interaction[i]=signif(anova(fit0,fit_models[[i]],test="Chisq")$`Pr(>|Chi|)`[2],3)
}
rm(i,fit0,fit_models)
#write.csv(P_interaction,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/Interaction_pain_status_chronic.csv')





##################################interaction for each pain type####################################################
P_interaction=data.frame(matrix(ncol = ncol(pain_duration), nrow = 11))
row.names(P_interaction)=c('Age','Sex','Race','Deprivation','Income','University_education','Smoking','Alcohol','Sedentary','Obesity','Genetic risk')
colnames(P_interaction)=colnames(pain_duration)

for (index_pain in c(1:ncol(pain_duration))) {
  datas$pain_site=pain_duration[,index_pain]
  datas$pain_site=factor(datas$pain_site)
  
  fit0<- coxph(Surv(time, status) ~ pain_site+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit1<- coxph(Surv(time, status) ~ pain_site*age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit2<- coxph(Surv(time, status) ~ pain_site*gender+age+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit3<- coxph(Surv(time, status) ~ pain_site*ethnicity+age+gender+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit4<- coxph(Surv(time, status) ~ pain_site*deprivation+age+gender+ethnicity+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit5<- coxph(Surv(time, status) ~ pain_site*income+age+gender+ethnicity+deprivation+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit6<- coxph(Surv(time, status) ~ pain_site*university_education+age+gender+ethnicity+deprivation+income+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit7<- coxph(Surv(time, status) ~ pain_site*smoking+age+gender+ethnicity+deprivation+income+university_education+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit8<- coxph(Surv(time, status) ~ pain_site*alcohol+age+gender+ethnicity+deprivation+income+university_education+smoking+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit9<- coxph(Surv(time, status) ~ pain_site*sedentary+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit10<-coxph(Surv(time, status) ~ pain_site*obesity+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit11<-coxph(Surv(time, status) ~ pain_site*prs+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = datas)
  fit_models=list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
  rm(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11)
  
  for (i in c(1:11)) {
    P_interaction[i,index_pain]=signif(anova(fit0,fit_models[[i]],test="Chisq")$`Pr(>|Chi|)`[2],3)
  }
  rm(i,fit0,fit_models)
}
#write.csv(P_interaction,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/interaction_pain_type.csv')
###################################################################################################################
##########################################subgroup analysis by each covariate for each pain type######################################
#################need to manually specify the pain type in codes
result<-data.frame(id=seq(1:37),
                   HR_single=NA,P_single=NA,CI_single=NA,lower_single=NA,upper_single=NA,
                   HR_multiple=NA,P_multiple=NA,CI_multiple=NA,lower_multiple=NA,upper_multiple=NA)
row.names(result)=c('Sex','Females','Males','Race','Ethnic minorities','White','Deprivation','Lower','Middle','Higher',
                    'Income','Unknown','Low (<£51999)','Middle (£52k-£100k)','High (>£100k)','Education','No education','Less than college','Above College',
                    'Smoking status','Never','Ever','Alcohol intake status','Never_0','Ever_1',
                    'Sedentary behavior','0-4 hours/day','>4 hours/day',
                    'BMI category','Healthy weight','Overweight','Obese','Severe Obese','Genetic risk','Low','Intermediate','High')

for (pain_status in c(1:2)) {#acute/chronic
  datas$pain=pain_duration$facial##########pain type
  datas$pain=factor(datas$pain)
  
  temp_data <- subset(datas, pain==0|pain==pain_status)
  temp_data$pain=factor(temp_data$pain)
{
  ##sex
  fit1<- coxph(Surv(time, status) ~ pain:gender+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1,2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+1,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+1,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+1,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+1,5*pain_status]=lower
    result[i+1,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##race
  fit1<- coxph(Surv(time, status) ~ pain:ethnicity+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1,2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+4,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+4,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+4,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+4,5*pain_status]=lower
    result[i+4,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##deprivation
  fit1<- coxph(Surv(time, status) ~ pain:deprivation+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+7,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+7,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+7,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+7,5*pain_status]=lower
    result[i+7,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##income
  fit1<- coxph(Surv(time, status) ~ pain:income+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:4)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+11,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-4+i,'exp(coef)']),2)
    result[i+11,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-4+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-4+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-4+i,'upper .95'],2)
    result[i+11,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+11,5*pain_status]=lower
    result[i+11,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##education
  fit1<- coxph(Surv(time, status) ~ pain:university_education+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+16,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+16,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+16,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+16,5*pain_status]=lower
    result[i+16,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##smoking
  fit1<- coxph(Surv(time, status) ~ pain:smoking+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+20,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+20,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+20,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+20,5*pain_status]=lower
    result[i+20,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##alcohol
  fit1<- coxph(Surv(time, status) ~ pain:alcohol+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+23,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+23,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+23,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+23,5*pain_status]=lower
    result[i+23,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##sedentary
  fit1<- coxph(Surv(time, status) ~ pain:sedentary+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:2)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+26,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
    result[i+26,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
    result[i+26,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+26,5*pain_status]=lower
    result[i+26,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##bmi category
  fit1<- coxph(Surv(time, status) ~ pain:obesity+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:4)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+29,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-4+i,'exp(coef)']),2)
    result[i+29,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-4+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-4+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-4+i,'upper .95'],2)
    result[i+29,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+29,5*pain_status]=lower
    result[i+29,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
  ##genetic risk
  fit1<- coxph(Surv(time, status) ~ pain:prs+age+gender+ethnicity+deprivation+income+university_education+
                 smoking+alcohol+sedentary+obesity+
                 prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,data = temp_data)
  for (i in c(1:3)) {
    nrow_summary<-nrow(summary(fit1)$conf.int)
    result[i+34,5*pain_status-3]=round((summary(fit1)$conf.int[nrow_summary-3+i,'exp(coef)']),2)
    result[i+34,5*pain_status-2]=signif((summary(fit1)$coefficients[nrow_summary-3+i,'Pr(>|z|)']),3)
    lower=round(summary(fit1)$conf.int[nrow_summary-3+i,'lower .95'],2)
    upper=round(summary(fit1)$conf.int[nrow_summary-3+i,'upper .95'],2)
    result[i+34,5*pain_status-1]=paste('(',lower,'-',upper,')',sep ='')
    result[i+34,5*pain_status]=lower
    result[i+34,5*pain_status+1]=upper
    rm(lower,upper,nrow_summary)
  }
  rm(i,fit1)
}
  rm(temp_data)
}
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/subgroup_hr_facial.csv')
#####################################################################################################################################

####################################HR 0/0.5/1/2,extended covariates##################################################
hr <- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+
              smoking+alcohol+sedentary+obesity+
              prs+
              vascular_history+cancer+diabetes+medication+sleep+
              comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10,
            data = datas)
coeff<-summary(hr)
result=data.frame(HR=round(coeff$coefficients[,2],2),lower=round(coeff$conf.int[,3],2),upper=round(coeff$conf.int[,4],2),p_value=signif(coeff$coefficients[,5],3))
for (i in 1:nrow(result)){
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
rm(i)
{
  result<-rbind(c(1, 1, 1, 1, '1.00 (Ref)'),
                result[1:4,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[5,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[6,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[7:8,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[9:11,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[12:13,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[14,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[15,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[16,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[17:19,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[20:21,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[22,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[23,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[24,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[25,],c(1, 1, 1, NA, '1.00 (Ref)'),
                result[26,])
  row.names(result)[c(1,6,8,10,13,17,20,22,24,26,30,33,35,37,39,41)]=c('Pain free','Female','Ethnic minorities','Deprivation_low','Income_Unknown',
                                                        'No education','Never_smoking','Never drinking','0-4 hours/day','Healthy weight','Low genetic risk',
                                                        'Vascular_history_No','cancer_history_No','Diabetes_No','Medication_No','Sleep_unhealth')
  result$Participants=NA
  result$Events=NA
}#pain
unique_pain=c(0,0.5,1,2)
for (i in 1:4){
  result$Participants[i]=length(which(datas$pain==unique_pain[i]))
  result$Events[i]=length(which(datas$pain==unique_pain[i]&(datas$status==1)))
}
rm(unique_pain)
{
  #age
  result$Participants[5]=length(datas$age)
  result$Events[5]=length(which(datas$status==1))
  #gender
  for (i in 6:7){
    result$Participants[i]=length(which(datas$gender==i-6))
    result$Events[i]=length(which(datas$gender==(i-6)&datas$status==1))
  }
  #ethnicity
  for (i in 8:9){
    result$Participants[i]=length(which(datas$ethnicity==i-8))
    result$Events[i]=length(which(datas$ethnicity==(i-8)&datas$status==1))
  }
  #deprivation
  for (i in 10:12){
    result$Participants[i]=length(which(datas$deprivation==i-10))
    result$Events[i]=length(which(datas$deprivation==(i-10)&datas$status==1))
  }
  #income
  for (i in 13:16){
    result$Participants[i]=length(which(datas$income==i-13))
    result$Events[i]=length(which(datas$income==(i-13)&datas$status==1))
  }
  #education
  for (i in 17:19){
    result$Participants[i]=length(which(datas$university_education==i-17))
    result$Events[i]=length(which(datas$university_education==(i-17)&datas$status==1))
  }
  #smoking
  for (i in 20:21){
    result$Participants[i]=length(which(datas$smoking==i-20))
    result$Events[i]=length(which(datas$smoking==(i-20)&datas$status==1))
  }
  #alcohol
  for (i in 22:23){
    result$Participants[i]=length(which(datas$alcohol==i-22))
    result$Events[i]=length(which(datas$alcohol==(i-22)&datas$status==1))
  }
  #sedentary
  for (i in 24:25){
    result$Participants[i]=length(which(datas$sedentary==i-24))
    result$Events[i]=length(which(datas$sedentary==(i-24)&datas$status==1))
  }
  #Obesity
  for (i in 26:29){
    result$Participants[i]=length(which(datas$obesity==i-26))
    result$Events[i]=length(which(datas$obesity==(i-26)&datas$status==1))
  }
  #genetic risk
  for (i in 30:32){
    result$Participants[i]=length(which(datas$prs==i-30))
    result$Events[i]=length(which(datas$obesity==(i-30)&datas$status==1))
  }
  #vascular history
  for (i in 33:34){
    result$Participants[i]=length(which(datas$vascular_history==i-33))
    result$Events[i]=length(which(datas$vascular_history==(i-33)&datas$status==1))
  }
  #cancer history
  for (i in 35:36){
    result$Participants[i]=length(which(datas$cancer==i-35))
    result$Events[i]=length(which(datas$cancer==(i-35)&datas$status==1))
  }
  #diabetes history
  for (i in 37:38){
    result$Participants[i]=length(which(datas$diabetes==i-37))
    result$Events[i]=length(which(datas$diabetes==(i-37)&datas$status==1))
  }
  #medication
  for (i in 39:40){
    result$Participants[i]=length(which(datas$medication==i-39))
    result$Events[i]=length(which(datas$medication==(i-39)&datas$status==1))
  }
  #sleep
  for (i in 41:42){
    result$Participants[i]=length(which(datas$vascular_history==i-41))
    result$Events[i]=length(which(datas$vascular_history==(i-41)&datas$status==1))
  }
  rm(coeff,i)
}
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_extended_covariates.csv')

###############################################interaction: extended covariates#########################################
fit0<- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit1<- coxph(Surv(time, status) ~ pain*age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit2<- coxph(Surv(time, status) ~ pain*gender+age+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit3<- coxph(Surv(time, status) ~ pain*ethnicity+age+gender+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit4<- coxph(Surv(time, status) ~ pain*deprivation+age+gender+ethnicity+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit5<- coxph(Surv(time, status) ~ pain*income+age+gender+ethnicity+deprivation+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit6<- coxph(Surv(time, status) ~ pain*university_education+age+gender+ethnicity+deprivation+income+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit7<- coxph(Surv(time, status) ~ pain*smoking+age+gender+ethnicity+deprivation+income+university_education+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit8<- coxph(Surv(time, status) ~ pain*alcohol+age+gender+ethnicity+deprivation+income+university_education+smoking+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit9<- coxph(Surv(time, status) ~ pain*sedentary+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit10<-coxph(Surv(time, status) ~ pain*obesity+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit11<-coxph(Surv(time, status) ~ pain*prs+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication+sleep,data = datas)
fit12<-coxph(Surv(time, status) ~ pain*vascular_history+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+cancer+diabetes+medication+sleep,data = datas)
fit13<-coxph(Surv(time, status) ~ pain*cancer+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+diabetes+medication+sleep,data = datas)
fit14<-coxph(Surv(time, status) ~ pain*diabetes+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+medication+sleep,data = datas)
fit15<-coxph(Surv(time, status) ~ pain*medication+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+sleep,data = datas)
fit16<-coxph(Surv(time, status) ~ pain*sleep+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10+vascular_history+cancer+diabetes+medication,data = datas)
fit_models=list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14,fit15,fit16)
rm(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14,fit15,fit16)
P_interaction=data.frame(matrix(ncol = 1, nrow = 16))
row.names(P_interaction)=c('Age','Sex','Race','Deprivation','Income','University_education','Smoking','Alcohol','Sedentary','Obesity','Genetic risk','Vascular history','Cancer history','Diabetes history','Medication','Sleep')
colnames(P_interaction)=c('P_interaction')
for (i in c(1:16)) {
  P_interaction$P_interaction[i]=signif(anova(fit0,fit_models[[i]],test="Chisq")$`Pr(>|Chi|)`[2],3)
}
rm(i,fit0,fit_models)


####################################Pain status and PRS group#####################
index=which(datas$pain!=0.5)
datas=datas[index,]
rm(index)

index_0_0=which(datas$prs==0&datas$pain==0)
index_0_1=which(datas$prs==0&datas$pain==1)
index_0_2=which(datas$prs==0&datas$pain==2)
index_1_0=which(datas$prs==1&datas$pain==0)
index_1_1=which(datas$prs==1&datas$pain==1)
index_1_2=which(datas$prs==1&datas$pain==2)
index_2_0=which(datas$prs==2&datas$pain==0)
index_2_1=which(datas$prs==2&datas$pain==1)
index_2_2=which(datas$prs==2&datas$pain==2)
datas$prs_pain=NA
datas$prs_pain[index_0_0]='00'
datas$prs_pain[index_0_1]='01'
datas$prs_pain[index_0_2]='02'
datas$prs_pain[index_1_0]='10'
datas$prs_pain[index_1_1]='11'
datas$prs_pain[index_1_2]='12'
datas$prs_pain[index_2_0]='20'
datas$prs_pain[index_2_1]='21'
datas$prs_pain[index_2_2]='22'
datas$prs_pain=factor(datas$prs_pain)
model <- coxph(Surv(time, status) ~ prs_pain+age+gender+ethnicity+deprivation+income+university_education+
              smoking+alcohol+sedentary+obesity+
              comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
# model <- coxph(Surv(time, status) ~ prs_pain+age+gender+ethnicity+deprivation+income+university_education+
#               smoking+alcohol+sedentary+obesity+
#               vascular_history+cancer+diabetes+medication+sleep+
#               comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)

coeff=summary(model)
rm(index_0_0,index_0_1,index_0_2,index_1_0,index_1_1,index_1_2,index_2_0,index_2_1,index_2_2,model)

result=data.frame(matrix(nrow=9,ncol=7))
row.names(result)=c('Low_PF','Low_SCP','Low_MCP','Intermediate_PF','Intermediate_SCP','Intermediate_MCP','High_PF','High_SCP','High_MCP')
colnames(result)=c('HR','lower','upper','CI','P_value','Participants','Events')

for (i in 2:9) {
  result$HR[i]=round(coeff$coefficients[i-1,2],2)
  result$P_value[i]=signif(coeff$coefficients[i-1,5],2)
  result$lower[i]=round(coeff$conf.int[i-1,3],2)
  result$upper[i]=round(coeff$conf.int[i-1,4],2)
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
i=1
for (genetic_risk in 1:3){
  for (pain_status in 1:3) {
    result$Participants[i]=length(which(datas$prs==genetic_risk-1&datas$pain==pain_status-1))
    result$Events[i]=length(which(datas$prs==genetic_risk-1&datas$pain==pain_status-1&datas$status==1))
    i=i+1
  }
}
result[1,1:5]=c(1, 1, 1, 1, '1.00 (Ref)')
rm(genetic_risk,pain_status,i,coeff)
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_prs_chronic.csv')

########################R, additive interaction###########################################
library(interactionR)
result<-data.frame(id=seq(1:2),reri_mp=NA,ap_mp=NA,s_mp=NA,reri_sp=NA,ap_sp=NA,s_sp=NA)
colnames(result)=c('index','reri_sp','ap_sp','s_sp','reri_mp','ap_mp','s_mp')

#metabolic_syndrome
for (pain_index in 1:2) {
  datas1=subset(datas,(pain==0|pain==pain_index)&(prs==0|prs==2))
  datas1$pain=factor(datas1$pain)
  datas1$prs=factor(datas1$prs)
  model <- coxph(Surv(time, status) ~ pain*prs+
                   age+gender+ethnicity+deprivation+income+university_education+
                   smoking+alcohol+sedentary+obesity+
                   comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas1)
  
  summary(model)
  out<-interactionR(model, exposure_names = c('pain','prs'),ci.type = "delta",ci.level = 0.95,em=F)
  outcome<-out$dframe
  if(pain_index==1){
    result$reri_sp[1]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
    result$ap_sp[1]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
    result$s_sp[1]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
  }else{
    result$reri_mp[1]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
    result$ap_mp[1]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
    result$s_mp[1]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
  }
  rm(model,outcome,datas1,out)
}

#########################################demographic information for pain/no pain#####################
#age in years
index_no_pain=which(datas$pain==0)
index_pain=which(datas$pain!=0)
datas$pain[index_no_pain]=1
datas$pain[index_pain]=2
rm(index_no_pain,index_pain)

summary_information=data.frame(matrix(ncol = 3, nrow = 39))
colnames(summary_information)=c('Whole population','Pain free','With pain')
row.names(summary_information)=c('Total N','Age at baseline, years (mean, std)','Sex','Females','Males','Race','Ethnic minorities','White','Deprivation','Lower','Middle','Higher',
                                 'Income','Unknown','Low (<£51999)','Middle (£52k-£100k)','High (>£100k)','Education','No education','Less than college','Above College',
                                 'Smoking status','Never','Ever','Alcohol intake status','Never_0','Ever_1',
                                 'Sedentary behavior','0-4 hours/day','>4 hours/day','BMI category','Normal weight','Overweight','Obese','Severe Obese',
                                 'Genetic risk','Low','Intermediate','High')
for (pain_status in c(10,1,2)) {
  if (pain_status==10){
    index_pain=c(1:nrow(datas))#10=whole samples
    column_index=1
  }else{
    index_pain=which(datas$pain==pain_status)
    column_index=pain_status+1
  }
  temp_datas=datas[index_pain,]
  rm(index_pain)
  #N
  summary_information[1,column_index]=nrow(temp_datas)
  #age,age is numeric
  summary_information[2,column_index]=paste(round(mean(temp_datas$age),2),' (',round(sd(temp_datas$age),2),')',sep='')
  #gender
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$gender==i-1))
    summary_information[i+3,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #ethnicity
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$ethnicity==i-1))
    summary_information[i+6,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #deprivation
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$deprivation==i-1))
    summary_information[i+9,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #income
  for (i in c(1:4)) {
    n_temp=length(which(temp_datas$income==i-1))
    summary_information[i+13,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #education
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$university_education==i-1))
    summary_information[i+18,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #smoking
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$smoking==i-1))
    summary_information[i+22,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #alcohol
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$alcohol==i-1))
    summary_information[i+25,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #sedentary
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$sedentary==i-1))
    summary_information[i+28,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #bmi_category
  for (i in c(1:4)) {
    n_temp=length(which(temp_datas$obesity==i-1))
    summary_information[i+31,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #Genetic risk
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$prs==i-1))
    summary_information[i+36,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  rm(temp_datas)
  rm(i)
}
rm(column_index,pain_status)
#write.csv(summary_information,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/demographics_status_pain_no_pain.csv')


#########################################demographic information for 0/0.5/1/2#####################
#age in years
datas$pain=as.numeric(as.character(datas$pain))
index_short=which(datas$pain==0.5)
index_single=which(datas$pain==1)
index_multiple=which(datas$pain==2)
datas$pain[index_short]=1
datas$pain[index_single]=2
datas$pain[index_multiple]=3
rm(index_short,index_single,index_multiple)

summary_information=data.frame(matrix(ncol = 5, nrow = 39))
colnames(summary_information)=c('Whole population','Pain free','Short-term pain','Single chronic pain','Multiple chronic pain')
row.names(summary_information)=c('Total N','Age at baseline, years (mean, std)','Sex','Females','Males','Race','Ethnic minorities','White','Deprivation','Lower','Middle','Higher',
                                 'Income','Unknown','Low (<£51999)','Middle (£52k-£100k)','High (>£100k)','Education','No education','Less than college','Above College',
                                 'Smoking status','Never','Ever','Alcohol intake status','Never_0','Ever_1',
                                 'Sedentary behavior','0-4 hours/day','>4 hours/day','BMI category','Normal weight','Overweight','Obese','Severe Obese',
                                 'Genetic risk','Low','Intermediate','High')
for (pain_status in c(10,0,1,2,3)) {
  if (pain_status==10){
    index_pain=c(1:nrow(datas))#10=whole samples
    column_index=1
  }else{
    index_pain=which(datas$pain==pain_status)
    column_index=pain_status+2
  }
  temp_datas=datas[index_pain,]
  rm(index_pain)
  #N
  summary_information[1,column_index]=nrow(temp_datas)
  #age,age is numeric
  summary_information[2,column_index]=paste(round(mean(temp_datas$age),2),' (',round(sd(temp_datas$age),2),')',sep='')
  #gender
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$gender==i-1))
    summary_information[i+3,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #ethnicity
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$ethnicity==i-1))
    summary_information[i+6,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #deprivation
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$deprivation==i-1))
    summary_information[i+9,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #income
  for (i in c(1:4)) {
    n_temp=length(which(temp_datas$income==i-1))
    summary_information[i+13,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #education
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$university_education==i-1))
    summary_information[i+18,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #smoking
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$smoking==i-1))
    summary_information[i+22,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #alcohol
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$alcohol==i-1))
    summary_information[i+25,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #sedentary
  for (i in c(1:2)) {
    n_temp=length(which(temp_datas$sedentary==i-1))
    summary_information[i+28,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #bmi_category
  for (i in c(1:4)) {
    n_temp=length(which(temp_datas$obesity==i-1))
    summary_information[i+31,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  #Genetic risk
  for (i in c(1:3)) {
    n_temp=length(which(temp_datas$prs==i-1))
    summary_information[i+36,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
  }
  rm(n_temp)
  rm(temp_datas)
  rm(i)
}
rm(column_index,pain_status)
#write.csv(summary_information,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/demographics_status_chronic.csv')






##############################################weighted pain score-mutual adjustment##############################################
library("survival")
library("survminer")
new_data<-data.frame(matrix(0,nrow=nrow(datas),ncol=ncol(pain_duration)*2))
colnames(new_data)<-c(paste(colnames(pain_duration),'_acute',sep =''),paste(colnames(pain_duration),'_chronic',sep =''))

for (i in c(1:ncol(pain_duration))) {
  index_acute=which(pain_duration[,i]==1)
  new_data[index_acute,i]=1
  new_data[,i]=factor(new_data[,i])
  index_chronic=which(pain_duration[,i]==2)
  new_data[index_chronic,i+ncol(pain_duration)]=1
  new_data[,i+ncol(pain_duration)]=factor(new_data[,i+ncol(pain_duration)])
  rm(index_acute,index_chronic)
}
rm(i)

new_data=cbind(datas[,1:24],new_data)#16=original status+time+covariates
hr <- coxph(Surv(time, status) ~ headache_acute+headache_chronic+facial_acute+facial_chronic+neck_shoulder_acute+neck_shoulder_chronic+
              back_acute+back_chronic+stomach_acute+stomach_chronic+hip_pain_acute+hip_pain_chronic+knee_acute+knee_chronic+
              age+gender+ethnicity+deprivation+income+university_education+
              smoking+alcohol+sedentary+obesity+
              prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = new_data)
coeff<-summary(hr)$coefficients
rm(new_data,hr)

###############construct weighted pain scores#####################
new_data<-data.frame(matrix(0,nrow=nrow(datas),ncol=ncol(pain_duration)*2))
colnames(new_data)<-c(paste(colnames(pain_duration),'_acute',sep =''),paste(colnames(pain_duration),'_chronic',sep =''))

for (i in c(1:ncol(pain_duration))) {
  index_acute=which(pain_duration[,i]==1)
  new_data[index_acute,i]=1
  index_chronic=which(pain_duration[,i]==2)
  new_data[index_chronic,i+ncol(pain_duration)]=1
  rm(index_acute,index_chronic)
}
rm(i)
beta<-coeff[1:(ncol(pain_duration)*2),1]
pain_score<-(new_data$headache_acute*beta[1]+new_data$headache_chronic*beta[2]+new_data$facial_acute*beta[3]+new_data$facial_chronic*beta[4]+
               new_data$neck_shoulder_acute*beta[5]+new_data$neck_shoulder_chronic*beta[6]+new_data$back_acute*beta[7]+new_data$back_chronic*beta[8]+
               new_data$stomach_acute*beta[9]+new_data$stomach_chronic*beta[10]+new_data$hip_pain_acute*beta[11]+new_data$hip_pain_chronic*beta[12]+
               new_data$knee_acute*beta[13]+new_data$knee_chronic*beta[14])/sum(beta)
rm(beta,new_data)
datas$pain_score=pain_score
rm(pain_score,coeff)


#############################construct weighted pain scores/indvidual beta###########################
beta=result$beta
significance=result$P_value
for (pain_index in 1:ncol(pain_duration)) {
  index_acute=which(pain_duration[,pain_index]==1)
  index_chronic=which(pain_duration[,pain_index]==2)
  pain_duration[index_acute,pain_index]=beta[pain_index*2-1]*(significance[pain_index*2-1]<0.05/7)
  pain_duration[index_chronic,pain_index]=beta[pain_index*2]*(significance[pain_index*2]<0.05/7)
  print(c(beta[pain_index*2-1],(significance[pain_index*2-1]<0.05/7),beta[pain_index*2],(significance[pain_index*2]<0.05/7)))
  rm(index_acute,index_chronic)
}
datas$pain_score=rowSums(pain_duration,na.rm=TRUE)
rm(beta,significance,pain_index)

########tertiles of weighted scores########
datas$weighted_pain_status=0
scores_excluding_0=datas$pain_score[which(datas$pain_score!=0)]
tertile_point=quantile(scores_excluding_0, probs = seq(0, 1, 1/3),na.rm = TRUE)
rm(scores_excluding_0)
index_low=which(datas$pain_score<tertile_point[2]&datas$pain_score>=tertile_point[1])
index_middle=which(datas$pain_score>=tertile_point[2]&datas$pain_score<=tertile_point[3])
index_high=which(datas$pain_score>tertile_point[3])
datas$weighted_pain_status[index_low]=1
datas$weighted_pain_status[index_middle]=2
datas$weighted_pain_status[index_high]=3
index_pain_free=which(datas$pain_score==0)
datas$weighted_pain_status[index_pain_free]=0
rm(index_low,index_middle,index_high,tertile_point,index_pain_free)
datas$weighted_pain_status=factor(datas$weighted_pain_status)

ordinal_model <- coxph(Surv(time, status) ~ weighted_pain_status+
                         age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+
                         prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
summary(ordinal_model)

#########write ordinal models to csv#########
coeff<-summary(ordinal_model)
result=data.frame(HR=round(coeff$coefficients[,2],2),lower=round(coeff$conf.int[,3],2),upper=round(coeff$conf.int[,4],2),p_value=signif(coeff$coefficients[,5],3))
for (i in 1:nrow(result)){
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
result=result[1:3,]
rm(i)
result<-rbind(c(1, 1, 1, 1, '1.00 (Ref)'),result[1:3,])
row.names(result)=c('Pain free','Lower tertile','Middle tertile','Higher tertile')
result$Participants=NA
result$Events=NA

for (i in 1:4){
  result$Participants[i]=length(which(datas$weighted_pain_status==i-1))
  result$Events[i]=length(which(datas$weighted_pain_status==i-1&(datas$status==1)))
}
rm(coeff,i)
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_composite_tertile.csv')
####################################nonlinear between No and HR####################
library(survival)
library(rms)
library(dplyr)
library(Gmisc)
library(splines)
library(Greg)
library(tidyverse)
dd <- datadist(datas)
options(datadist='dd')

fit<- cph(Surv(time, status) ~ rcs(pain_score,3) +age+gender+ethnicity+deprivation+income+university_education+
            smoking+alcohol+sedentary+obesity+
            prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
anova(fit)#test the significance of nolinear regression
plotHR(fit, term="pain_score", xlab="Composite pain scores",se=TRUE,
       polygon_ci=TRUE,
       col.term = "#0070b9",
       col.se = "#DEEBF7BB",
       lwd.term=2,
       lty.term = 1,
       ylog=FALSE,
       alpha=0.05,#alpha level, 95%
       cex=1,
       axes = TRUE,
       plot.bty="l", xlim=c(0,2.1),ylim=c(0,3.5),rug="density")



##################################################plot hr for 8 types ##################################################
library(forestplot)
library(dplyr)
library(ggthemes)

rm(list=ls())
#result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_pain_locations.csv',header=T)
#result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_pain_locations_extended_covariates.csv',header=T)
result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_pain_locations_exclude_dementia_2_year.csv',header=T)
#insert a new row to better show the range of x-lab
result1=rbind(result1,result1[1,])
result1$HR[nrow(result1)]=1
result1$lower[nrow(result1)]=0.5
result1$upper[nrow(result1)]=2.5
row.names(result1)[nrow(result1)]='Null'
result_HR <- tibble(HR  = result1$HR, pain_type = result1$X,p_value=as.character(result1$P_value),ci=result1$CI,
                    lower = result1$lower,upper = result1$upper,participants=as.character(result1$Participants),
                    events=as.character(result1$Events))
#rm(result)
header <- tibble(pain_type='pain_type',p_value='P value',ci='HR (95% CI)',participants='Participants',events='Events')
combined_data <- bind_rows(header,result_HR)
rm(header,result_HR)
combined_data %>% 
  forestplot(labeltext = c(pain_type, ci, participants,events,p_value),
             mean=HR,
             graph.pos =3,boxsize=0.35,zero=1,
             is.summary = FALSE,clip = c(0.2,2.5), xlog = FALSE,graphwidth = unit(.28,"npc"),
             col = fpColors(box = "#3c5a66",line = "black",summary = "grey"),
             mar=unit(rep(1.0, times = 4), "cm"),
             txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.6), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
             xlab = "HR (95% CI)",lwd.xaxis = 1.0)





###############################plot interaction for each pain type########################################################

result<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/interaction_pain_type.csv',header=T)
row.names(result)=result$X
result1=result[,2:8]
result1=-log10(result1)
colnames(result1)=c('Headache','facial pain','Neck pain','Back pain','Stomach pain','Hip pain','Knee pain')

library(ggcorrplot)
library(RColorBrewer)
mycolors<-(brewer.pal(9, "Greens"))[1:9]

ggcorrplot(t(result1), method = "circle",type = "full", outline.color = "white",show.diag = TRUE,tl.cex = 8,
           lab='TRUE',lab_size=3,lab_col = 'black',
           ggtheme = ggplot2::theme_classic)+ scale_fill_gradientn(limit = c(0,3.5),colours=mycolors)


#################################plot hr for subgroup#############################################################
library(forestplot)
library(dplyr)
library(ggthemes)

mydata<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/subgroup_hr_status_chronic.csv',header=T)
mydata=rbind(mydata,mydata[1,])
mydata$HR_single[nrow(mydata)]=0.5
mydata$lower_single[nrow(mydata)]=2
mydata$upper_single[nrow(mydata)]=2
mydata$HR_multiple[nrow(mydata)]=2
mydata$lower_multiple[nrow(mydata)]=2
mydata$upper_multiple[nrow(mydata)]=2
mydata$X[nrow(mydata)]='Null'
result1=mydata[,c(1,3,6,7)]
result1=mydata[,c(1,8,11,12)]
colnames(result1)=c('subgroup','HR','lower','upper')

result_HR <- tibble(HR  = result1$HR, subgroup = result1$subgroup,lower = result1$lower,upper = result1$upper)
#rm(result)
header <- tibble(subgroup='subgroup')
combined_data <- bind_rows(header,result_HR)
rm(header,result_HR)
combined_data %>% 
  forestplot(labeltext = c(subgroup),
             mean=HR,
             graph.pos =2,boxsize=0.38,zero=1,
             is.summary = FALSE,clip = c(0.5,2.5), xlog = FALSE,graphwidth = unit(.32,"npc"),
             col = fpColors(box = "#3c5a66",line = "black",summary = "grey"),
             mar=unit(rep(1.0, times = 4), "cm"),
             txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.6), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
             xlab = "HR (95% CI)",lwd.xaxis = 1.0)
##################################################plot hr for prs*pain ##################################################
library(forestplot)
library(dplyr)
library(ggthemes)

rm(list=ls())
#result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status-chronic.csv',header=T)
#result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_acute.csv',header=T)
result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_acute_excluding_dementia.csv',header=T)
result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_acute_extended_covariates.csv',header=T)
result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_excluding_dementia.csv',header=T)
result1=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_frailty_status_extended_covariates.csv',header=T)
#insert a new row to better show the range of x-lab
result1=result1[c(1,3,4),]
result1=rbind(result1,result1[1,])
result1$HR[nrow(result1)]=1.5
result1$lower[nrow(result1)]=1.5
result1$upper[nrow(result1)]=1.5
row.names(result1)[nrow(result1)]='Null'
result_HR <- tibble(HR  = result1$HR, pain_type = result1$X,p_value=as.character(result1$p_value),ci=result1$CI,
                    lower = result1$lower,upper = result1$upper,participants=as.character(result1$Participants),
                    events=as.character(result1$Events))
header <- tibble(pain_type='pain_type',p_value='P value',ci='HR (95% CI)',participants='Total N',events='Events')
combined_data <- bind_rows(header,result_HR)
rm(header,result_HR)
combined_data %>% 
  forestplot(labeltext = c(pain_type, participants,events,ci,p_value),
             mean=HR,
             graph.pos =4,boxsize=0.35,zero=1,
             is.summary = FALSE,clip = c(0,1.5), xlog = FALSE,graphwidth = unit(.25,"npc"),
             col = fpColors(box = "#3c5a66",line = "black",summary = "grey"),
             mar=unit(rep(1.0, times = 4), "cm"),
             txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.6), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
             xlab = "HR (95% CI)",lwd.xaxis = 1.0)


library(ggplot2)
result<-read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_prs_chronic.csv',header=T)
result<-read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/datas/HR_prs_acute.csv',header=T)
#behavior_order=order((result$HR))
ggplot(result, aes(x=X, y=HR)) +
  geom_pointrange(aes(ymin=lower, ymax=upper))+
  geom_point(data=result, mapping=aes(x=X, y=HR,fill=X), size=2.6, shape=21,stroke = 0.9)+
  scale_x_discrete(limits=c(result$X))+ylim(0.5,3.8)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.85))

##########################pain medication#####################
index=which(datas$pain==1|datas$pain==2)
datas1=datas[index,]
hr <- coxph(Surv(time, status) ~ medication+age+gender+ethnicity+deprivation+income+university_education+
              smoking+alcohol+sedentary+vascular_history+cancer+obesity+diabetes, data = datas1)
summary(hr)

index_pain_yes=which(datas1$pain!=0&datas1$medication==1)
index_pain_no=which(datas1$pain!=0&datas1$medication==0)


# index_no_medication_no_pain=which(datas$pain==0)
# index_no_medication_pain=which(datas$pain!=0&datas$medication==0)
# index_medication_pain=which(datas$pain!=0&datas$medication==1)
# datas$medication[index_no_medication_no_pain]=0
# datas$medication[index_no_medication_pain]=0
# datas$medication[index_medication_pain]=1
# datas$medication=factor(datas$medication)
# hr <- coxph(Surv(time, status) ~ pain+age+gender+ethnicity+deprivation+income+university_education+smoking+alcohol+sedentary+obesity+
#               vascular_history+cancer+diabetes+medication+sleep+
#               prs+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10, data = datas)
# summary(hr)




###################################################mdd pgc no ukb and 23andMe###################################################
for (i in c(1:9)) {
  exposure=exposures[[i]]
  write.csv(exposure,'/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/TwoSampleMR/exposure.csv')
  exposure<-system.file('exposure.csv',package = "TwoSampleMR")
  exposure_exp_dat<-read_exposure_data(filename=exposure,sep=",",snp_col = "rsid",beta_col = "beta",
                                       se_col = "se",
                                       effect_allele_col = "a_0",
                                       other_allele_col = "a_1",
                                       eaf_col = "af",
                                       pval_col = "p")
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  exposure_dat_clumped<-clump_data(exposure_exp_dat,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 1,clump_p2 = 1,pop="EUR")
  rm(exposure,exposure_exp_dat)
  
  # outcome<-fread('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-PD/GCST90275127.tsv',header=T)
  # c<-merge(exposure_dat_clumped,outcome,by.x = "SNP",by.y = "rs_id")
  # write.csv(c,'/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/TwoSampleMR/outcome.csv')
  # rm(c)
  # outcome<-system.file('outcome.csv',package = "TwoSampleMR")
  # outcome_dat<-read_outcome_data(snps=exposure_dat_clumped$SNP,
  #                                filename=outcome,sep=",",snp_col = "SNP",beta_col = "beta",
  #                                se_col = "standard_error",
  #                                effect_allele_col = "effect_allele",
  #                                other_allele_col = "other_allele",
  #                                pval_col = "p_value")
  # dat<-harmonise_data(exposure_dat =  exposure_dat_clumped, outcome_dat=outcome_dat)
  # retained_index=which(dat$mr_keep==TRUE)
  # dat<-dat[retained_index,]
  # 
  outcome_data<-extract_outcome_data(snps=exposure_dat_clumped$SNP,outcomes='ieu-b-7')
  data<-harmonise_data(exposure_dat_clumped,outcome_data)
  
  retained_index=which(data$mr_keep==TRUE)
  dat<-data[retained_index,]
  
  rm(exposure_dat_clumped,outcome_dat,outcome,retained_index)
  result=mr(dat,method_list = c("mr_egger_regression_bootstrap","mr_ivw_fe","mr_ivw_mre","mr_simple_median","mr_weighted_median"))
  results$p_mr[(5*i-4):(5*i)]=signif(result$pval,3)
  
  results$nsnp[(5*i)]=result$nsnp[1]
  
  or=generate_odds_ratios(result)
  results$or_mr[(5*i-4):(5*i)]=round(or$or,2)
  
  results$lower_mr[(5*i-4):(5*i)]=round(or$or_lci95,2)
  results$upper_mr[(5*i-4):(5*i)]=round(or$or_uci95,2)
  
  pleitropy<-mr_pleiotropy_test(dat)
  results$p_pleiotropy[5*i]=signif(pleitropy$pval,3)
  results$egger_intercept[5*i]=signif(pleitropy$egger_intercept,3)
  
  hetegeneity<-mr_heterogeneity(dat)
  results$p_hetegeneity[5*i]=signif(hetegeneity$Q_pval[1],3)
  results$p_hetegeneity[5*i-1]=signif(hetegeneity$Q_pval[2],3)
  results$q_hetegeneity[5*i]=hetegeneity$Q[1]
  results$q_hetegeneity[5*i-1]=hetegeneity$Q[2]
  rm(or,pleitropy,hetegeneity,result,dat)
}



