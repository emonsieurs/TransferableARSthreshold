#####################################################################
#####################################################################
# ************************************************************** ####
#                       ALGORITHM FOR THE                        ####
#               AUTOMATIC CALCULATION OF RAINFALL 		     ####	
#	         THRESHOLDS FOR LANDSLIDE OCCURRENCE               ####
#             BASED ON AN ANTECEDENT RAINFALL INDEX              ####
#                 AND LANDSLIDE SUSCEPTIBILITY                   ####
#         	        - without bootstrapping - 			     ####
#                                                                ####
#  AUTHOR: ELISE MONSIEURS                                       #### 
#  CONTACT: elise.monsieurs@africamuseum.be 			     ####
#  DATE CREATED: September 2019					     ####		 
#  LICENSE: this code is licensed under GPL(version 2 or later)  ####
#  CITATION: Monsieurs, E., Dewitte, O., Depicker, A.            ####
#  Demoulin, A. 2019. Towards a transferable antecedent          ####
#  rainfall – susceptibility threshold approach for              ####
#  landsliding. Water							     ####
# ************************************************************** ####
#  This script was prepared using                         	     ####
#  R Core Team (2017). R: A language and environment for         ####
#  statistical computing. R Foundation for Statistical           ####
#  Computing, Vienna,Austria.                                    ####
#  URL http://www.R-project.org/                                 ####
#  R version: R-3.4.3                                            ####
# ************************************************************** ####
#  The script requires the following input data (stored in       ####
#  one folder named 'INPUT'):                                    ####
#  1) LS_Inventory.csv, including following fields:              ####
#	   - ID: ID landslides					           ####
#	   - Date: Date landslide occurrence (dd/mm/yyyy)	     ####
#  2) LS_Susceptibility.csv, including following fields:         ####
#	   - ID: ID landslides					           ####
#	   - Susceptibility: Landslide susceptibility [0,1]	     #### 
#  3) SRE_Daily.csv                                              ####
#	   - Date: dd/mm/yyyy						     ####	
#	   - ID landslides: daily rainfall(mm) for each landslide  ####
#	in the inventory (column header is the landslide ID)       ####
#	This script was prepared with daily satellite rainfall     ####
#	estimates (SRE) from TMPA 3B42 RT                          ####
# ************************************************************** ####
# No Packages are required to be installed.                      ####
# ************************************************************** ####
#####################################################################
#####################################################################


#INPUT SETTINGS*******************************************************************************************************
	#Set the working directory where the INPUT folder is located	
		wdMAIN<-"C:/Documents/AR-S_Threshold/"
		wdIN<-paste(wdMAIN,"INPUT/",sep="")
		NDirOUT<-paste(wdMAIN,"OUTPUT_NoBootstrap",sep="")
		dir.create(NDirOUT)
		wdOUT<-paste(NDirOUT,"/",sep="")

	#Set desired threshold probability of exceedance level
	#Here, thresholds at the 5% and 10% exceedance probability levels are calibrated 
		TPE<-c(0.05,0.1) 
		
	#Set a, b, n parameters for antecedent rainfall calculation
		n<-42 #days
		n_min<-n-1
		n.char<-paste(n,"Days",sep=" ")
		aa<--1.2
		aa.char<-as.character(abs(aa))
		b<-1.2

	#Set ARmin
		ARmin<-5

	#Read data landslide inventory 
		filename<-paste(wdIN,"LS_Inventory.csv",sep="")
		LSinv<-read.csv(filename, header=TRUE)
		TotalLS<-dim(LSinv)[1]

	#Read data Satellite Rainfall Estimates 
		SREall<-read.csv(paste(wdIN,"SRE_Daily.csv",sep="")) 
		SREcol<-colnames(SREall)

	#Read Susceptibility data  
		SuscINPUT<-read.csv(paste(wdIN,"LS_Susceptibility.csv",sep=""))  

	#Set Breaks for 10 S classes 
		BreaksN<-11

	#Identify axis labels
		xls<-"LS Susceptibility"
		xlr<-"Antecedent Rain (mm)"

	#Set vectors for threshold results for multiple probability levels of exceedance
	beta<-vector(mode = "numeric", length = length(TPE))
	alpha<-vector(mode = "numeric", length = length(TPE))
	R_Sqr<-vector(mode = "numeric", length = length(TPE))
	P_CSA<-vector(mode = "numeric", length = length(TPE))
	P_CSB<-vector(mode = "numeric", length = length(TPE))
	FNR<-vector(mode = "numeric", length = length(TPE))
	Threshold_Results<-data.frame()


#***************************************  A. DATA PREPARATION  ************************************************
#A.1 Extract AR values for landslide reported day ± 1 

	#A.1.1 Calculate AR time series
	#Initiate dataframe for AR values of landslide dates ± 1 days
		LS_AntR<-data.frame(matrix(ncol = 3, nrow = TotalLS*3)) 
		names(LS_AntR)<-c("ID","Date","AntR")

	#Initiate row values for LS_AntR dataframe
		p<-1
		d<-2
		a<-3

	#Start loop through landslides for AR time series calculation
	for (i in 1:TotalLS){
	
	#LS details
	ID<-as.character(LSinv$ID[i],stringsAsFactors = FALSE)
	Date_LS<-as.Date(as.character(LSinv$Date[i],stringsAsFactors = FALSE),format="%d/%m/%Y")

	#Read LS event complete daily time series 
	indSRE_LS<-which(SREcol==ID)
	SREls<-SREall[,c(1,indSRE_LS)] 
	RI<-SREls[,2] 
	lSREls<-length(RI) 

	#Initialize new time series for calculated Antecedent rainfall values
	lAR<-list()

		#Start loop time series
		for (j in n:lSREls){
		R0<-RI[j]  
		AR<-R0

		#Calculate weighted rain 
		for (k in 1:n_min){AR<-AR+(RI[j-k]*(exp((aa*k)/(RI[j-k])^b)))}  
		lAR[j]<-AR

		#End loop time series
		}
	
	#Add dates to new rainfall time series
	Dates<-as.Date(as.character(SREls[n:lSREls,1],stringsAsFactors=FALSE),format="%d/%m/%Y")
	Date_ARls<-data.frame(Date=Dates, AntecR=unlist(lAR[n:lSREls]),stringsAsFactors=FALSE)
	names(Date_ARls)<-c("Date",ID)	

	#A.1.2 Extract AR of Day landslide (D), one day Prior to D (P), and one day After D(A)
	Ind_LS_D<-which(Date_ARls$Date==Date_LS)
	Ind_LS_P<-Ind_LS_D-1
	Ind_LS_A<-Ind_LS_D+1
	
	AntecR_D<-Date_ARls[Ind_LS_D,2]
	AntecR_P<-Date_ARls[Ind_LS_P,2]
	AntecR_A<-Date_ARls[Ind_LS_A,2]
	
	ID_D<-paste(ID,"_D",sep="")
	ID_P<-paste(ID,"_P",sep="")
	ID_A<-paste(ID,"_A",sep="")
	
	Date_LS_P<-as.Date(Date_ARls[Ind_LS_P,1],format = "%d/%m/%Y")
	Date_LS_A<-as.Date(Date_ARls[Ind_LS_A,1],format = "%d/%m/%Y")
	
	#Save dates and antecedent rainfall relevant for the threshold calibration
	LS_AntR$ID[p]<-ID_P
	LS_AntR$Date[p]<-as.character(Date_LS_P)
	LS_AntR$AntR[p]<-AntecR_P
	
	LS_AntR$ID[d]<-ID_D
	LS_AntR$Date[d]<-as.character(as.Date(Date_LS,format = "%d/%m/%Y"))
	LS_AntR$AntR[d]<-AntecR_D
	
	LS_AntR$ID[a]<-ID_A
	LS_AntR$Date[a]<-as.character(Date_LS_A)
	LS_AntR$AntR[a]<-AntecR_A
	
	p<-p+3
	d<-d+3
	a<-a+3

	#End loop landslides
	}

	#A.1.3 Discard data with AR<ARmin from data set
	LS_AntR<-LS_AntR[LS_AntR$AntR>=ARmin,] 
	q<-dim(LS_AntR)[1]


#A.2 Weigh data for event date uncertainty

	#A.2.1 Attribute weight to data
	WeightAttr<-data.frame(matrix(0,dim(LS_AntR)[1],1))
	for (j in 1:dim(LS_AntR)[1]) {
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="D") {WeightAttr[j,1]<-as.numeric(2/3)}
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="P") {WeightAttr[j,1]<-as.numeric(1/6)}
	if (substr(LS_AntR$ID[j], nchar(LS_AntR$ID[j]), nchar(LS_AntR$ID[j]))=="A") {WeightAttr[j,1]<-as.numeric(1/6)}
	}
	LS_AntR$Weight<-WeightAttr[,1]


	# Add landslide susceptibility data
	names(SuscINPUT)<-c("ID2","Susceptibility")
	#Extract original ID data from LS_AntR
	ID<-list()
	for (i in 1:dim(LS_AntR)[1]){
	ID[i]<-substr(LS_AntR$ID[i], 1,nchar(LS_AntR$ID[i])-2)
	}
	ID<-unlist(ID)
	LS_AntR$ID2<-ID
	#Merge landslide susceptibility data with LS_AntR
	LS_AntR_S<-merge(LS_AntR,SuscINPUT,by="ID2")
	
	#A.2.2 Create expanded inventory according to weights
	LS_AntR_S$SampleFreq<-0
	for (i in 1:dim(LS_AntR_S)[1]){
	if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="D"){LS_AntR_S$SampleFreq[i]<-4}
	if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="P"){LS_AntR_S$SampleFreq[i]<-1}
	if(substr(LS_AntR_S$ID[i],nchar(LS_AntR_S$ID[i]),nchar(LS_AntR_S$ID[i]))=="A"){LS_AntR_S$SampleFreq[i]<-1}
	}
	Susc_AR<- LS_AntR_S[rep(row.names(LS_AntR_S), LS_AntR_S$SampleFreq), 1:length(colnames(LS_AntR_S))-1]
	r<-dim(Susc_AR)[1]

#**************************************  B. THRESHOLD CALIBRATION  ********************************************

	#Start loop probability of exceedance levels
	for (FL in 1:length(TPE)){
	LowestAR<-data.frame()

#B.1 Calculate number of data to be selected per S class 
	tc<-(2*TPE[FL]*r)/(BreaksN-1) 

#B.2 Select lowest AR data per S class

	#B.2.1 Group data according to S classes
	#Set Breaks for 10 logarithmic equidistant classes 
	breaks<-exp(seq(log(min(Susc_AR$Susceptibility)), log(max(Susc_AR$Susceptibility)), length.out=BreaksN))
	S1 <- Susc_AR[which(Susc_AR$Susceptibility<=breaks[2]),]
	S2 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[2] & Susc_AR$Susceptibility<=breaks[3]),]
	S3 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[3] & Susc_AR$Susceptibility<=breaks[4]),]
	S4 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[4] & Susc_AR$Susceptibility<=breaks[5]),]
	S5 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[5] & Susc_AR$Susceptibility<=breaks[6]),]
	S6 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[6] & Susc_AR$Susceptibility<=breaks[7]),]
	S7 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[7] & Susc_AR$Susceptibility<=breaks[8]),]
	S8 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[8] & Susc_AR$Susceptibility<=breaks[9]),]
	S9 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[9] & Susc_AR$Susceptibility<=breaks[10]),]
	S10 <- Susc_AR[which(Susc_AR$Susceptibility>breaks[10]),]

	#B.2.2 Select data with lowest AR per S class until tc is reached 
	ItLar<-1
	S1_S<-S1[order(S1$AntR),]
	while((ItLar-1)<tc){
	if(dim(S1_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S1_S$ID2[ItLar],ID=S1_S$ID[ItLar],Date=S1_S$Date[ItLar], AntR=S1_S$AntR[ItLar],
				Weight=S1_S$Weight[ItLar],Susceptibility=S1_S$Susceptibility[ItLar],Sclass="1"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S2_S<-S2[order(S2$AntR),]
	while((ItLar-1)<tc){
	if(dim(S2_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S2_S$ID2[ItLar],ID=S2_S$ID[ItLar],Date=S2_S$Date[ItLar], AntR=S2_S$AntR[ItLar],
				Weight=S2_S$Weight[ItLar],Susceptibility=S2_S$Susceptibility[ItLar],Sclass="2"))
	ItLar<-ItLar+1
	}

	ItLar<-1	
	S3_S<-S3[order(S3$AntR),]
	while((ItLar-1)<tc){
	if(dim(S3_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S3_S$ID2[ItLar],ID=S3_S$ID[ItLar],Date=S3_S$Date[ItLar], AntR=S3_S$AntR[ItLar],
				Weight=S3_S$Weight[ItLar],Susceptibility=S3_S$Susceptibility[ItLar],Sclass="3"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S4_S<-S4[order(S4$AntR),]
	while((ItLar-1)<tc){
	if(dim(S4_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S4_S$ID2[ItLar],ID=S4_S$ID[ItLar],Date=S4_S$Date[ItLar],AntR=S4_S$AntR[ItLar],
				Weight=S4_S$Weight[ItLar],Susceptibility=S4_S$Susceptibility[ItLar],Sclass="4"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S5_S<-S5[order(S5$AntR),]
	while((ItLar-1)<tc){
	if(dim(S5_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S5_S$ID2[ItLar],ID=S5_S$ID[ItLar],Date=S5_S$Date[ItLar], AntR=S5_S$AntR[ItLar],
				Weight=S5_S$Weight[ItLar],Susceptibility=S5_S$Susceptibility[ItLar],Sclass="5"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S6_S<-S6[order(S6$AntR),]
	while((ItLar-1)<tc){
	if(dim(S6_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S6_S$ID2[ItLar],ID=S6_S$ID[ItLar],Date=S6_S$Date[ItLar], AntR=S6_S$AntR[ItLar],
				Weight=S6_S$Weight[ItLar],Susceptibility=S6_S$Susceptibility[ItLar],Sclass="6"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S7_S<-S7[order(S7$AntR),]
	while((ItLar-1)<tc){
	if(dim(S7_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S7_S$ID2[ItLar],ID=S7_S$ID[ItLar],Date=S7_S$Date[ItLar], AntR=S7_S$AntR[ItLar],
				Weight=S7_S$Weight[ItLar],Susceptibility=S7_S$Susceptibility[ItLar],Sclass="7"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S8_S<-S8[order(S8$AntR),]
	while((ItLar-1)<tc){
	if(dim(S8_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S8_S$ID2[ItLar],ID=S8_S$ID[ItLar],Date=S8_S$Date[ItLar], AntR=S8_S$AntR[ItLar],
				Weight=S8_S$Weight[ItLar],Susceptibility=S8_S$Susceptibility[ItLar],Sclass="8"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S9_S<-S9[order(S9$AntR),]
	while((ItLar-1)<tc){
	if(dim(S9_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S9_S$ID2[ItLar],ID=S9_S$ID[ItLar],Date=S9_S$Date[ItLar], AntR=S9_S$AntR[ItLar],
				Weight=S9_S$Weight[ItLar],Susceptibility=S9_S$Susceptibility[ItLar],Sclass="9"))
	ItLar<-ItLar+1
	}

	ItLar<-1
	S10_S<-S10[order(S10$AntR),]
	while((ItLar-1)<tc){
	if(dim(S10_S)[1]<ItLar){break}
	LowestAR<-rbind(LowestAR, data.frame(ID2=S10_S$ID2[ItLar],ID=S10_S$ID[ItLar],Date=S10_S$Date[ItLar], AntR=S10_S$AntR[ItLar],
				Weight=S10_S$Weight[ItLar],Susceptibility=S10_S$Susceptibility[ItLar],Sclass="10"))
	ItLar<-ItLar+1
	}

	write.csv(LowestAR,paste(wdOUT,"LowestAR_TPE",TPE[FL]*100,"Perc.csv",sep=""),row.names=FALSE)

#B.3 Linear regression with log transformed AR-S data

	AR<-LowestAR$AntR
	Susc<-LowestAR$Susceptibility

	#POWER LAW RELATION
	log.AR<-log10(AR)
      log.susc<-log10(Susc)
      fit.straight.line <- lm(log.AR~ log.susc) 
	
	a.straight<-coef(fit.straight.line)[1]
      b.straight<-coef(fit.straight.line)[2]

      #Calculation of the AR threshold 
	alpha[FL]<-round(10^a.straight,1)
	beta[FL]<-round(b.straight,5) 

	R_Sqr[FL]<-summary(fit.straight.line)$r.squared
	Sig_Alpha<-summary(fit.straight.line)$coeff[1,4]<0.05 
	Sig_Beta<-summary(fit.straight.line)$coeff[2,4]<0.05
	ifelse(Sig_Alpha==TRUE,P_CSA[FL]<- "Alpha is significant at the 95% level", P_CSA[FL]<- "Alpha is not significant at the 95% level")
	ifelse(Sig_Beta==TRUE,P_CSB[FL]<- "Beta is significant at the 95% level", P_CSA[FL]<- "Beta is not significant at the 95% level")


#**************************************** C. THRESHOLD EVALUATION *********************************************

#Calculate False Negative Rate 

	FNsum<-0
	for(fn in 1:dim(Susc_AR)[1]){
		AR<-Susc_AR$AntR[fn]
		susc<-Susc_AR$Susceptibility[fn]
		Thr.AR<-alpha[FL]*susc^beta[FL]
			if(Susc_AR$AntR[fn]<Thr.AR){FNsum<-FNsum+1} 
		}

	FNR[FL]<-round(FNsum/r,2)


#Summerize results threshold calibration
	Threshold_Results<-rbind(Threshold_Results,data.frame(a=aa.char,b=b,n=n.char,TPE=TPE[FL],FNR=FNR[FL],
								Alpha=round(alpha[FL],3),AlphaSign=P_CSA[FL],
								Beta=round(beta[FL],3),BetaSign=P_CSB[FL], 
                                                Rsquare=round(R_Sqr[FL],2),BreakN=(BreaksN-1),
								p=TotalLS,q=q,r=r,tc=tc))



#End loop probability of exceedance levels
}


#*************************************** WRITE/PLOT RESULTS ***************************************************

#Write results threshold calibration
	filename2<-paste(wdOUT,"AR-S_Threshold_Results_",Sys.Date(),".csv",sep="") 
	write.csv(Threshold_Results, file=filename2,quote=F,row.names=F)

#Plot thresholds
	Susc<-Susc_AR$Susceptibility
	l_do<-length(Susc)
	x_vector<-seq.int(from=min(Susc), to=max(Susc),length.out=l_do)

	filename3<-paste(wdOUT,"AR-S_Threshold.pdf",sep="")  
	pdf(file=filename3,paper = "special", width =30,height =21,pointsize=30,colormodel="srgb",
	useDingbats=F,fillOddEven=T,version="1.7")

	Prob_66<-Susc_AR[which(Susc_AR$Weight>0.4),]
	Prob_17<-Susc_AR[which(Susc_AR$Weight<0.4),]

	par(mgp=c(2.2,0.8,0),mar=c(5,4,5,2))

	plot(Prob_66$Susceptibility,Prob_66$AntR,type="p",yaxs="r",xaxs="r",log="xy",xlab=xls,ylab=xlr,pch=21,col="gray7",bg="blue",cex.lab=1.35,cex=1.2)
	points(Prob_17$Susceptibility,Prob_17$AntR,type="p",yaxs="r",xaxs="r",pch=21,col="gray7",bg="blue",cex=0.7)

	#Draw threshold lines
      lines(x_vector,alpha[1]*x_vector^beta[1],col="darkgreen")
      lines(x_vector,alpha[2]*x_vector^beta[2],col="darkred")           

	#Threshold labels
      label_1 = bquote("T" ~ .(format(paste(TPE[1]*100,"%",sep="")), digits = 1)~":"  ~  italic(AR)==~ .(format(alpha[1], digits = 2,nsmall = 1)) ~ 
                        italic(S)^ .(format(beta[1], digits = 3,nsmall = 1))~ " [R²="~ .(format(R_Sqr[1],digits=2,nsmall=2)) ~"]")

      label_2 = bquote("T" ~ .(format(paste(TPE[2]*100,"%",sep="")), digits = 1)~":" ~  italic(AR)==~ .(format(alpha[2], digits = 2,nsmall = 1)) ~ 
                        italic(S)^ .(format(beta[2], digits = 3,nsmall = 1))~ " [R²="~ .(format(R_Sqr[2],digits=2,nsmall=2)) ~"]")

	#Add labels in margins
	mtext(label_1,side=3, line=0 , lty=1, bty="n", col="darkgreen",adj=1,cex=1.5)
	mtext(label_2,side=3, line=1.3 , lty=1, bty="n", col="darkred",adj=1,cex=1.5)
	mtext(paste("Ndata = ",dim(Susc_AR)[1],sep=""),side=3, line=0 , lty=1, bty="n", col="black",adj=0,cex=1.3)  

	# Highlight subset data
	FL_1<-read.csv(paste(wdOUT,"LowestAR_TPE",TPE[1]*100,"Perc.csv",sep=""))
	FL_2<-read.csv(paste(wdOUT,"LowestAR_TPE",TPE[2]*100,"Perc.csv",sep=""))
	FL_1_66<-FL_1[which(FL_1$Weight>0.4),]
	FL_1_17<-FL_1[which(FL_1$Weight<0.4),]
	FL_2_66<-FL_2[which(FL_2$Weight>0.4),]
	FL_2_17<-FL_2[which(FL_2$Weight<0.4),] 

	points(FL_2_66$Susceptibility,FL_2_66$AntR,type="p",yaxs="r",xaxs="r",pch=21,col="gray7",bg="darkred",cex=1.2)
	points(FL_2_17$Susceptibility,FL_2_17$AntR,type="p",yaxs="r",xaxs="r",pch=21,col="gray7",bg="darkred",cex=0.7)
	points(FL_1_66$Susceptibility,FL_1_66$AntR,type="p",yaxs="r",xaxs="r",pch=21,col="gray7",bg="darkgreen",cex=1.2)
	points(FL_1_17$Susceptibility,FL_1_17$AntR,type="p",yaxs="r",xaxs="r",pch=21,col="gray7",bg="darkgreen",cex=0.7)

	# Add S classes
	abline(v=breaks,lty=3,yaxs="r",lwd=4.5,col="darkgrey")

	# Add legend
	legend("bottomleft",c("Day LS event","Day ± 1 LS event"),pt.bg="white",col="gray7",lwd=NA,pch=21,cex=1,pt.cex=c(1.2,0.7),pt.lwd=1,y.intersp = 1.25,bg="white")	
   
	# Plot parameter settings
      par(mfrow=c(1,1),mar=c(1, 1,1, 1))
      plot.new()
	text(0,0.99, "Antecedent rainfall", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1)
	text(0,0.94, "Function: AR=sum((e^(-a*t/r^b))*r)", lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.91, paste("a = ",aa.char,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
      text(0,0.88, paste("b = ",b,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.85, paste("Accumulation period = ",n.char,sep=""), lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)

	text(0,0.75, "Significance threshold parameters", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1)
	text(0,0.70, paste("T",TPE[1]*100,": ",P_CSA[1],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.67, paste("T",TPE[1]*100,": ",P_CSB[1],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.64, paste("T",TPE[2]*100,": ",P_CSA[2],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.61, paste("T",TPE[2]*100,": ",P_CSB[2],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)

	text(0,0.51, "Threshold evaluation", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1)
	text(0,0.46, paste("Obtained FNR for TPE level set at ",TPE[1], ": ",FNR[1],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.43, paste("Obtained FNR for TPE level set at ",TPE[2], ": ",FNR[2],sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)

	text(0,0.33, "Determination coefficient R²", lwd=2, lty=1, bty="n", col="blue",adj=0,cex=1)
	text(0,0.28, paste("R² of least-square fit of subset for TPE level set at ",TPE[1], ": ",round(R_Sqr[1],2),sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)
	text(0,0.25, paste("R² of least-square fit of subset for TPE level set at ",TPE[2], ": ",round(R_Sqr[2],2),sep=""),lwd=2, lty=1, bty="n", col="black",adj=0,cex=0.8)

	
dev.off()
