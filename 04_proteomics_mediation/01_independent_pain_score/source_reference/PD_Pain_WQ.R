############ ICD-10 extraction of first diagnosis PD
data_path <- "PD_data.csv"
.libPaths(c(normalizePath("R_libs"), .libPaths()))
expand_multi_select <- function(x) {
  raw_values <- as.character(x)
  raw_values[is.na(raw_values)] <- ""
  split_values <- strsplit(raw_values, "\\|")
  max_len <- max(lengths(split_values))
  expanded <- matrix(NA_real_, nrow=length(raw_values), ncol=max_len)
  for (i in seq_along(split_values)) {
    values <- suppressWarnings(as.numeric(split_values[[i]]))
    if (length(values)>0) {
      expanded[i,seq_along(values)] <- values
    }
  }
  as.data.frame(expanded)
}

get_multi_select <- function(dat, field_id) {
  cols <- grep(paste0("^X",field_id,"\\.0\\."), names(dat), value=TRUE)
  if (length(cols)==1) {
    return(expand_multi_select(dat[[cols]]))
  }
  as.data.frame(dat[,cols,drop=FALSE])
}
{
  library(data.table)
  first_diagnosis<-fread(data_path,select=c('eid','53-0.0',
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
  datas<-fread(data_path, 
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
  datas<-fread(data_path, select=
                 c('eid','40000-0.0','53-0.0'), header = T,sep=",")
  colnames(datas)=c('id','death_date','attending_date')
  index_death<-which(!is.na(datas$death_date)==TRUE)
  death_information<-datas[index_death,]
  rm(datas)
  rm(index_death)
  death_information$time_interval=death_information$death_date-death_information$attending_date
  
  
  ################ censoring data, days between the date attending the assessment center and the date recording information
  library(data.table)
  datas<-fread(data_path, select=c('eid','53-0.0'),header = T,sep=",")
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
  datas<-fread(data_path,header = T,sep=",")
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
  datas<-fread(data_path,header = T,sep=",")
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
  
  education=get_multi_select(datas,"6138")
  for (i in seq_len(ncol(education))){
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
  university_education=do.call(pmax,c(education,list(na.rm = TRUE)))
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
  vascular=get_multi_select(datas,"6150")
  for (i in seq_len(ncol(vascular))){
    index_no_answer=which(vascular[,i]==-3)
    vascular[index_no_answer,i]=NA
    rm(index_no_answer)
  }
  rm(i)
  vascular_history=do.call(pmin,c(vascular,list(na.rm = TRUE)))
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
  medication=get_multi_select(datas,"6154")
  for (i in seq_len(ncol(medication))){
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
  pain_medication_max=do.call(pmax,c(medication,list(na.rm = TRUE)))
  pain_medication_min=do.call(pmin,c(medication,list(na.rm = TRUE)))
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
  participant_eid=datas$eid
  rm(datas,medication)
}

##############################combine data###########################################
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,
                age=age,gender=gender,sites=sites,ethnicity=ethnicity,
                university_education=university_education,deprivation=deprivation,income=income,
                sedentary=sedentary,smoking=smoking,alcohol=alcohol,obesity=obesity,prs=pd_prs,
                comp_1=comp_1,comp_2=comp_2,comp_3=comp_3,comp_4=comp_4,comp_5=comp_5,
                comp_6=comp_6,comp_7=comp_7,comp_8=comp_8,comp_9=comp_9,comp_10=comp_10,
                eid=participant_eid)

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
   pd_prs,comp_1,comp_2,comp_3,comp_4,comp_5,comp_6,comp_7,comp_8,comp_9,comp_10,participant_eid)
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


##############################################weighted pain score-mutual adjustment##############################################
library("survival")
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
              prs, data = new_data)
#+comp_1+comp_2+comp_3+comp_4+comp_5+comp_6+comp_7+comp_8+comp_9+comp_10

coeff<-summary(hr)$coefficients
index=which(coeff[1:14,'coef']<0)
coeff[index,'coef']=0
rm(new_data,hr,index)

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
write.csv(datas[,c("eid","pain_score","weighted_pain_status")],"composite_pain_scores.csv",row.names = FALSE)

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
write.csv(result,"HR_composite_tertile.csv")
