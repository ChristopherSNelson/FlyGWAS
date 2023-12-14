""" Chris Nelson 2013 FlyGWAS.py
This script is intended to analyze DGRP drosophila GWAS data. It takes in a genotype file in vcf format, and a phenotype file, and outputs pvalues by genotype marker.
Only intended for two alles or missing data at a given genotype marker.
VCF file format:

##fileformat=VCFv4.1
##source=DGRP2
##reference=dm3
##INFO=<ID=REFCOUNT,Number=1,Type=Integer,Description="Reference Allele Count">
##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description="Alternative Allele Count">
##INFO=<ID=ANNOTATION,Number=1,Type=String,Description="Variant effect">
##INFO=<ID=REFLINE,Number=1,Type=String,Description="Line numbers homozygous for the reference allele">
##INFO=<ID=ALTLINE,Number=1,Type=String,Description="Line numbers homozygous for the alternative allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	line_21	line_26	line_28	line_31	line_32	line_38	line_40	line_41	line_42	line_45	line_48	line_49	line_57	line_59	line_69	line_73	line_75	line_83	line_85	line_88	line_91	line_93	line_100	line_101	line_105	line_109	line_129	line_136	line_138	line_142	line_149	line_153	line_158	line_161	line_176	line_177	line_181	line_189	line_195	line_208	line_217	line_223	line_227	line_228	line_229	line_233	line_235	line_237	line_239	line_256	line_280	line_287	line_301	line_303	line_304	line_306	line_307	line_309	line_310	line_313	line_315	line_317	line_318	line_319	line_320	line_321	line_324	line_325	line_332	line_335	line_336	line_338	line_340	line_348	line_350	line_352	line_354	line_355	line_356	line_357	line_358	line_359	line_360	line_361	line_362	line_365	line_367	line_370	line_371	line_373	line_374	line_375	line_377	line_379	line_380	line_381	line_382	line_383	line_385	line_386	line_390	line_391	line_392	line_395	line_397	line_399	line_405	line_406	line_409	line_426	line_427	line_437	line_439	line_440	line_441	line_443	line_461	line_486	line_491	line_492	line_502	line_505	line_508	line_509	line_513	line_517	line_528	line_530	line_531	line_535	line_551	line_555	line_559	line_563	line_566	line_584	line_589	line_595	line_596	line_627	line_630	line_634	line_639	line_642	line_646	line_703	line_705	line_707	line_712	line_714	line_716	line_721	line_727	line_730	line_732	line_737	line_738	line_748	line_757	line_761	line_765	line_774	line_776	line_783	line_786	line_787	line_790	line_796	line_799	line_801	line_802	line_804	line_805	line_808	line_810	line_812	line_818	line_819	line_820	line_821	line_822	line_832	line_837	line_843	line_849	line_850	line_852	line_853	line_855	line_857	line_859	line_861	line_879	line_882	line_884	line_887	line_890	line_892	line_894	line_897	line_900	line_907	line_908	line_911	line_913
chr2L	5002	2L_5002_SNP	G	T	999	PASS	REFCOUNT=127;ALTCOUNT=1;ANNOTATION=|||;REFLINE=21,28,31,38,40,41,45,48,49,69,73,75,83,85,88,91,93,100,109,129,136,153,158,161,176,189,223,227,228,229,235,237,280,303,307,310,315,317,318,319,320,321,324,332,335,336,338,348,350,352,354,355,358,359,361,362,365,367,371,373,375,377,380,381,382,385,390,391,392,395,397,405,426,427,440,441,443,492,505,509,517,530,551,559,563,589,595,596,627,630,634,639,707,712,714,716,721,727,730,737,761,774,783,790,796,799,805,808,812,818,820,821,822,849,850,852,853,857,861,884,887,892,897,900,907,911,913;ALTLINE=837	GT	0/0	./.	0/0	0/0	./.	0/0	0/0	0/0	./.	0/0	0/0	0/0	./.	./.	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	./.	./.	0/0	0/0	0/0	./.	./.	./.	0/0	0/0	0/0	0/0	./.	./.	0/0	./.	./.	./.	0/0	0/0	0/0	0/0	./.	0/0	0/0	./.	./.	0/0	./.	./.	0/0	./.	./.	0/0	./.	0/0	./.	0/0	0/0	0/0	0/0	0/0	0/0	0/0	./.	0/0	0/0	0/0	0/0	./.	0/0	0/0	0/0	0/0	0/0	./.	./.	0/0	0/0	./.	0/0	0/0	./.	0/0	./.	0/0	0/0	./.	0/0	0/0	./.	0/0	0/0	0/0	./.	0/0	./.	0/0	0/0	0/0	0/0	0/0	./.	0/0	./.	./.	0/0	0/0	./.	./.	0/0	0/0	0/0	./.	./.	./.	0/0	./.	0/0	./.	0/0	./.	0/0	./.	0/0	./.	./.	0/0	./.	0/0	0/0	./.	./.	0/0	0/0	0/0	0/0	0/0	0/0	0/0	1/1	./.	./.	./.	0/0	0/0	0/0	0/0	0/0	0/0	0/0	./.	0/0	./.	./.	./.	0/0	./.	0/0	./.	0/0	./.	./.	0/0	0/0	0/0	./.	./.	./.	0/0	0/0	./.	0/0	0/0	./.	0/0	0/0	0/0	./.	1/1	./.	0/0	0/0	0/0	0/0	./.	0/0	./.	0/0	./.	./.	0/0	0/0	./.	0/0	./.	0/0	0/0	0/0	./.	0/0	0/0
chr2L	5040	2L_5040_SNP	G	A	999	PASS	REFCOUNT=1;ALTCOUNT=118;ANNOTATION=|||;REFLINE=176;ALTLINE=21,28,31,38,40,41,45,49,69,73,83,85,93,109,129,136,153,158,161,189,223,228,229,233,235,280,301,303,307,310,315,317,318,319,320,321,324,332,335,338,348,350,352,354,355,357,358,359,360,361,362,365,367,371,373,375,377,380,381,382,390,391,392,395,397,399,405,426,427,440,441,443,486,492,505,517,530,559,589,595,596,634,639,646,707,712,714,716,721,727,737,761,774,783,790,796,799,805,808,812,818,820,821,822,837,850,852,853,857,861,884,887,892,897,900,908,911,913	GT	1/1	./.	1/1	1/1	./.	1/1	1/1	1/1	./.	1/1	./.	1/1	./.	./.	1/1	1/1	./.	1/1	1/1	./.	./.	1/1	./.	./.	./.	1/1	1/1	1/1	./.	./.	./.	1/1	1/1	1/1	0/0	./.	./.	1/1	./.	./.	./.	1/1	./.	1/1	1/1	1/1	1/1	./.	./.	./.	1/1	./.	1/1	1/1	./.	./.	1/1	./.	1/1	./.	1/1	1/1	1/1	1/1	1/1	1/1	1/1	./.	1/1	1/1	./.	1/1	./.	1/1	1/1	1/1	1/1	1/1	./.	1/1	1/1	1/1	1/1	1/1	1/1	1/1	1/1	./.	1/1	1/1	./.	1/1	1/1	./.	1/1	1/1	1/1	./.	./.	./.	1/1	1/1	1/1	1/1	1/1	1/1	1/1	./.	./.	1/1	1/1	./.	./.	1/1	1/1	1/1	./.	1/1	./.	1/1	./.	1/1	./.	./.	./.	1/1	./.	1/1	./.	./.	./.	./.	1/1	./.	./.	./.	1/1	1/1	1/1	./.	./.	1/1	1/1	./.	1/1	./.	./.	1/1	1/1	1/1	1/1	1/1	1/1	./.	./.	1/1	./.	./.	./.	1/1	./.	1/1	./.	1/1	./.	./.	1/1	1/1	1/1	./.	./.	./.	1/1	1/1	./.	1/1	1/1	./.	1/1	1/1	1/1	./.	1/1	./.	./.	1/1	1/1	1/1	./.	1/1	./.	1/1	./.	./.	1/1	1/1	./.	1/1	./.	1/1	1/1	./.	1/1	1/1	1/1

I'm typically working with input from FlyGWAS_KapahiOnlyFilterMAFandMissing.py, which has ~20 fewer columns than the above example.

usage: python FlyGWAS.py

see also FlyGWAS_KapahiOnlyFilterMAFandMissing.py for filtering missing genotype data and minor allele frequency (MAF).


In the output file 10th column=DR 11th column=AL 12th column= ratio pvalue.

###
4-14-2014
In this latest version I'm taking in lifespan phenotype and date the lifespan started and feeding into a multiple linear regression. I added a second loop doing the same analysis except 
throwing out date info and doing a Wilcoxon test, as previous versions in March and earlier had done.

###
8-28-14
In this version I'm reading from lifespan data and pulling out different percentiles of survival. Watch DRphen definition to make sure you're pulling from the right data columns in the phenotype data.

###
~9-16-14
converted to pandas handling, so the phenotype parsing simplified and the data structure changed

###
9-24-14
In this version I figured out that interaction terms of categorical variables require consideration for encoding. 
So I've changed diet(in the input file) and genotype(in the genotype parsing in theis script) to -1/+1 integers. 
The default 0 vs 1 "dummy" categorical encoding might have let interesting interactions slip through the cracks.
I should probably force the parsing of the phenotype file to be (-1/1) in case this convention is forgotten.

###
9-25-14
In this version I am using phenotype data that's been rescaled and recentered so mean=0 and stdev=1. We're calling this Daniel Promislow's idea.
As such there is no longer any week covariate or any mention of weeks in this script itself.
I'm keeping effects coding (-1/1) for diet(defined in the phenotype input file) and genotype.

###
1-8-15
Adding output to a directory for tidiness.
Adding MAF cutoff of variable stringency.

####
3-24-15
converting to chisquared test of MAF of longevity alleles. Kind of ala Barzilai/Suh


"""
import sys
from os import system
import datetime
#from statsmodels.regression.linear_model import OLS #this version of OLS supports R-style formulas for model definition.
import numpy as np
from time import time
from pandas import read_csv
from pandas import merge
from pandas import DataFrame
from numpy import nan
print "done importing"
start=time()
from scipy.stats import fisher_exact #3-24-15 function importing for guts of the test.
from scipy.stats import chi2_contingency #3-24-15 function importing for guts of the test.
MAFcutoff=0.25
DRShortClimbingThreshold=2.5#number of days used to decide between cases and controls
ALShortClimbingThreshold=2#picking thresholds based on top 15% of distribution. kinda arbitrary.

OutputDirectoryName="PERM_Climbing50pct_caseControl_MAF"+str(MAFcutoff)+"_"+str(datetime.date.today())
system("mkdir "+OutputDirectoryName)

#get phenotypes and covariatefrom phenotype file and put them in a dataframe.
for PhenoFileName in ["PERM_50percentClimbing_1.csv","PERM_50percentClimbing_2.csv","PERM_50percentClimbing_3.csv"]:
	PhenotypeDataFrame=read_csv(PhenoFileName)# Use (PhenoFileName,sep="\t") for tabbed files. Expects "strain", "Phenotype", "diet", and "week" headers. Column order shouldn't matter.
	#Convert position in the table to strain name with a list.
	#DGRP205LineList=['DGRP-21','DGRP-26','DGRP-28','DGRP-31','DGRP-32','DGRP-38','DGRP-40','DGRP-41','DGRP-42','DGRP-45','DGRP-48','DGRP-49','DGRP-57','DGRP-59','DGRP-69','DGRP-73','DGRP-75','DGRP-83','DGRP-85','DGRP-88','DGRP-91','DGRP-93','DGRP-100','DGRP-101','DGRP-105','DGRP-109','DGRP-129','DGRP-136','DGRP-138','DGRP-142','DGRP-149','DGRP-153','DGRP-158','DGRP-161','DGRP-176','DGRP-177','DGRP-181','DGRP-189','DGRP-195','DGRP-208','DGRP-217','DGRP-223','DGRP-227','DGRP-228','DGRP-229','DGRP-233','DGRP-235','DGRP-237','DGRP-239','DGRP-256','DGRP-280','DGRP-287','DGRP-301','DGRP-303','DGRP-304','DGRP-306','DGRP-307','DGRP-309','DGRP-310','DGRP-313','DGRP-315','DGRP-317','DGRP-318','DGRP-319','DGRP-320','DGRP-321','DGRP-324','DGRP-325','DGRP-332','DGRP-335','DGRP-336','DGRP-338','DGRP-340','DGRP-348','DGRP-350','DGRP-352','DGRP-354','DGRP-355','DGRP-356','DGRP-357','DGRP-358','DGRP-359','DGRP-360','DGRP-361','DGRP-362','DGRP-365','DGRP-367','DGRP-370','DGRP-371','DGRP-373','DGRP-374','DGRP-375','DGRP-377','DGRP-379','DGRP-380','DGRP-381','DGRP-382','DGRP-383','DGRP-385','DGRP-386','DGRP-390','DGRP-391','DGRP-392','DGRP-395','DGRP-397','DGRP-399','DGRP-405','DGRP-406','DGRP-409','DGRP-426','DGRP-427','DGRP-437','DGRP-439','DGRP-440','DGRP-441','DGRP-443','DGRP-461','DGRP-486','DGRP-491','DGRP-492','DGRP-502','DGRP-505','DGRP-508','DGRP-509','DGRP-513','DGRP-517','DGRP-528','DGRP-530','DGRP-531','DGRP-535','DGRP-551','DGRP-555','DGRP-559','DGRP-563','DGRP-566','DGRP-584','DGRP-589','DGRP-595','DGRP-596','DGRP-627','DGRP-630','DGRP-634','DGRP-639','DGRP-642','DGRP-646','DGRP-703','DGRP-705','DGRP-707','DGRP-712','DGRP-714','DGRP-716','DGRP-721','DGRP-727','DGRP-730','DGRP-732','DGRP-737','DGRP-738','DGRP-748','DGRP-757','DGRP-761','DGRP-765','DGRP-774','DGRP-776','DGRP-783','DGRP-786','DGRP-787','DGRP-790','DGRP-796','DGRP-799','DGRP-801','DGRP-802','DGRP-804','DGRP-805','DGRP-808','DGRP-810','DGRP-812','DGRP-818','DGRP-819','DGRP-820','DGRP-821','DGRP-822','DGRP-832','DGRP-837','DGRP-843','DGRP-849','DGRP-850','DGRP-852','DGRP-853','DGRP-855','DGRP-857','DGRP-859','DGRP-861','DGRP-879','DGRP-882','DGRP-884','DGRP-887','DGRP-890','DGRP-892','DGRP-894','DGRP-897','DGRP-900','DGRP-907','DGRP-908','DGRP-911','DGRP-913']# should be of length 205 (strain names).
	DGRP184LineList=[28122,28123,28124,28125,29651,28126,28127,28128,29652,28129,28130,28131,28132,28134,28274,28135,28136,28137,28138,28139,28140,28141,28142,28143,28144,28145,28146,28147,28148,28149,28150,28151,28152,28153,25174,28154,28155,28156,28157,29653,28159,28275,28160,28161,28162,28164,28165,25175,25176,25177,37525,25179,28166,28276,25180,25181,28167,28168,29654,29655,25182,28170,28171,25183,28172,28173,28174,28176,28177,28178,25184,25185,28179,25186,28180,25187,25445,28181,28182,28183,28184,28185,25188,28186,25189,25190,28188,28189,28190,28191,28192,25191,28194,25192,29656,29657,28278,28196,25193,25194,29658,28197,28198,28199,28200,25195,28202,28203,28204,28205,28206,29659,25197,29660,28207,28208,25198,28211,28212,28213,28215,25199,28216,28217,28218,25744,25200,25201,25745,28219,28220,28221,25202,25203,28222,28223,28224,28226,28227,25204,25205,28229,28230,25206,28231,28232,28233,25207,28234,28235,28236,28237,28238,28239,28240,28241,28242,25208,28243,28244,28245,28246,28247,28248,28249,25209,28250,28251,28252,25210,28253,28254,28255,28256,28279,28257,28258,28259,28260,28261,28262,28263,28264,28265]
	GenotypeDict={'genotype':[],'strain':DGRP184LineList}#initialize the genotype dictionary with strain names and a empty genotype list.
	resultsFileList=[]
	ChrList=["4"]#,"X","3R","3L","2R","2L"]#
	#ModelFormula='medianSurvival~diet*genotype+C(DateStarted)'#R-style formula definition. Moving it out of the loop for a tiny bit of efficiency.#C(week)
	FileNameStem="GWASoutput_caseControl_"+PhenoFileName[:-3]+"_"+str(datetime.date.today())+"_MAF"+str(MAFcutoff)+".txt"

	for ChrS in ChrList:#Now I'm done with the phenotype parsing, and I'm going to take chromosomes and run the stats on them in a series of large loops. Each chromosome takes around 25 minutes as of Jan 8th 2014.
		print "starting chromosome "+ChrS
		prunedGenotypeFileName="dgrp2_"+ChrS+".vcf_KapahiLinesOnlyVcf_filtered_missingData0.3_and_MAF0.05_.txt"# #vcfSnippet3_KapahiLinesOnlyVcf__filtered_missingData0.3_and_MAF0.05_.txt or dgrp2_KapahiLinesOnlyVcf_filtered_missingData0.3_and_MAF0.05_.txt  dgrp2_KapahiLinesOnlyVcf__filtered_missingData0.3_and_MAF0.05_Nov252013
		prunedGenotypeFile=open(prunedGenotypeFileName,"r")#reopening filtered genotype file #for each genotype line split phenotypes
		GWASoutputFileName=OutputDirectoryName+"/"+prunedGenotypeFileName[:8]+"_"+FileNameStem #building up the output file name
		resultsFileList.append(GWASoutputFileName)
		GWASoutputFile=open(GWASoutputFileName,"w")
		print"Splitting phenotypes, and performing stat test."
		for line in prunedGenotypeFile:
			tempFields=line.strip().split("\t")
			#print "tempFields length should be ~193 for the real data lines, longer tempfields cause file parsing errors."# it  was in the kapahi only vcf file filtering - line splits were wrong
			if "##" in line:
				#print "skipping doublehash"
				continue
			elif line[0]=="#":# find the most relevant header line explaining the marker info.Should be something like "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
					markerInfo=tempFields[2:8]#default
					markerInfo="\t".join(markerInfo[:3]+markerInfo[-1].split("|")[:3])#just trimming down and parsing the marker info.
					GWASoutputFile.write(markerinfo+"\n" )# writing the output file's header
					print "Wrote to header"
					continue		    
			
			markerInfo=tempFields[2:8]
			markerInfo="\t".join(markerInfo[:3]+markerInfo[-1].split("|")[:3])#just trimming down and parsing the marker info.
			
			genoList=[]#erase old genotypes from previous marker
			GenotypeDict['genotype']=[]#erase old genotypes from previous marker
			for field in tempFields[8:]:# if you're not NA, then add your phenotype to group 1 or group 0
				if field=="0/0":#find the 0 allele positions
					genoList.append(-1)#thought has gone into the choice of -1 vs +1 for integer encoding here. Too much to summarize in a short note, but it effects the construction and interpretation of the interaction term.
				elif field=="1/1":#find the 1 allele positions
					genoList.append(1)#thought has gone into the choice of -1 vs +1 for integer encoding here. Too much to summarize in a short note, but it effects the construction and interpretation of the interaction term.
				elif field=="./.":#find NA positions
					genoList.append(nan)
			#print len(genoList)#for Kapahi lines files this is 184, could go up to 205 if you include all lines.
			GenotypeDict['genotype']=genoList				
			geno_df=DataFrame(GenotypeDict)#converting to a pandas dataframe so we can merge it with the phenotypes and covariates dataframe.
			
			MergedDataFrame=merge(PhenotypeDataFrame,geno_df,on='strain',how='inner')# Merging the phenotype and genotype dataframes. "inner" means the intersection of the strain names (makes a smaller table of variable size)
			#print MergedDataFrame.dropna()# Sometimes this is empty if there's a type mismatch between the two strain columns.
			if len(MergedDataFrame)==0:
				sys.exit("Hit an empty merged dataframe! In ChrS = "+ChrS+" markerinfo="+str(markerInfo))#temporary exit while I debug the empty merged dataframe problem
			
			#11-11-14 Adding this bit so you can alter the MAF cutoff.
			CountPhen0=len(MergedDataFrame[(MergedDataFrame.genotype==-1)].genotype)#Phenotype missing data check is superfluous, should be already taken care of. It'd be ...&(MergedDataFrame.Phenotype!=nan) inside the brackets.
			CountPhen1=len(MergedDataFrame[(MergedDataFrame.genotype==1)].genotype)
			MinorCount=min(CountPhen0,CountPhen1)
			MajorCount=max(CountPhen0,CountPhen1)
			#print MinorCount,MajorCount
			MAF=float(MinorCount)/float(CountPhen0+CountPhen1)
			#print MAF
			if MAF<MAFcutoff:
				continue# Skips the stats test and doesn't return a result for low MAF markers. #3-24-15 might need to get more stringent here to avoid low n producing weird effects
				
			#3-24-15 assembling case-control counts based on phenotype thresholds set near the beginning of the script.
			#Fill out the two-by-two contingency table below (for each diet separately)
#                                            alleles
# 			            -1 reference allele     +1 alternative allele   
# 		long lifespan  |         #ofStrains     |        #ofStrains       |
# 		short lifespan |         #ofStrains     |        #ofStrains       |
# 			
# 			
#
			RefTable=MergedDataFrame[(MergedDataFrame.genotype==-1)]
			AltTable=MergedDataFrame[(MergedDataFrame.genotype==1)]
			DRrefList=list(RefTable[(RefTable.Diet==-1)].WeekBelow50pct) #splitting up into diets
			DRaltList=list(AltTable[(AltTable.Diet==-1)].WeekBelow50pct)
			ALrefList=list(RefTable[(RefTable.Diet==1)].WeekBelow50pct) 
			ALaltList=list(AltTable[(AltTable.Diet==1)].WeekBelow50pct)
			
			#print len(DRrefList)#early snps I'm looking at break down into length80 and length65 lists
			#print len(ALrefList)
			#print len(DRaltList)
			#print len(ALaltList)
			
			ALrefLongClimbingCount=0
			ALaltLongClimbingCount=0
			ALrefShortClimbingCount=0
			ALaltShortClimbingCount=0
			DRrefLongClimbingCount=0 
			DRaltLongClimbingCount=0
			DRrefShortClimbingCount=0
			DRaltShortClimbingCount=0
			
			for phenotype in ALrefList:
				if phenotype<=ALShortClimbingThreshold:
					ALrefShortClimbingCount+=1
				elif phenotype>ALShortClimbingThreshold:
					ALrefLongClimbingCount+=1
					
			for phenotype in ALaltList:
				if phenotype<=ALShortClimbingThreshold:
					ALaltShortClimbingCount+=1
				elif phenotype>ALShortClimbingThreshold:
					ALaltLongClimbingCount+=1
					
			for phenotype in DRrefList:
				if phenotype<=DRShortClimbingThreshold:
					DRrefShortClimbingCount+=1
				elif phenotype>DRShortClimbingThreshold:
					DRrefLongClimbingCount+=1
			
			for phenotype in DRaltList:
				if phenotype<=DRShortClimbingThreshold:
					DRaltShortClimbingCount+=1
				elif phenotype>DRShortClimbingThreshold:
					DRaltLongClimbingCount+=1		
		
# 			print ALrefShortLivedCount
# 			print ALaltShortLivedCount
# 			print ALrefShortLivedCount
# 			print ALaltShortLivedCount
# 			print DRrefShortLivedCount 
# 			print DRaltShortLivedCount
# 			print DRrefShortLivedCount
# 			print DRaltShortLivedCount
	
			ALcontingencyTable=[[ALrefShortClimbingCount,ALaltShortClimbingCount],[ALrefLongClimbingCount,ALaltLongClimbingCount]]
			DRcontingencyTable=[[DRrefShortClimbingCount,DRaltShortClimbingCount],[DRrefLongClimbingCount,DRaltLongClimbingCount]]
			AL_Pval =fisher_exact( ALcontingencyTable )[1]
			DR_Pval =fisher_exact( DRcontingencyTable )[1]
			AL_Pval_chi2 =chi2_contingency( ALcontingencyTable )[1]
			DR_Pval_chi2 =chi2_contingency( DRcontingencyTable )[1]			
			
			#print ALchiSqPval
			#print DRchiSqPval
			#print chi2_contingency( ALcontingencyTable )
			#print chi2_contingency( DRcontingencyTable )
			 
			GWASoutputFile.write(markerInfo+"\t"+str(DR_Pval_chi2)+"\t"+str(AL_Pval_chi2)+"\t"+str(DR_Pval)+"\t"+str(AL_Pval)+"\n")
			#sys.exit("debugging, stopping here: "+str(markerInfo))
			#MergedDataFrame.__delitem__('genotype')#this line deletes the genotype column from the dataframe and does it in place. Alternative is df=df.drop("geno",1)
		print "Done. Wrote to "+GWASoutputFileName
		prunedGenotypeFile.close()
		GWASoutputFile.close()
		print"FlyGWAS.py finishing Chromosome "+ChrS
		runlength=(time()-start)/60
		print "Run took "+str(runlength)+" minutes, so far."
	print "number of files to concatenate: "+str(len(resultsFileList))#typically 6 chromosomal segments
	#This next bit takes all the output files for the different chromosomes and crams them together.
	catOutputFileName="dgrp2_allChr_"+FileNameStem
	catShellCommand=" ".join(['cat']+resultsFileList+["> "+OutputDirectoryName+"/"+catOutputFileName])#Building up the command to send to cat shell command.
	print catShellCommand
	system(catShellCommand)#concatenates chromosome results files
	#Now I'm going to take the concatenated result file and get only the lines with low pvalues. The purpose of this is to reduce the size of the results file.
	thresholdPval=0.01#Adjust this for pval stringency.
	InFile=OutputDirectoryName+"/"+catOutputFileName
	InF=open(InFile,"r")
	OutFile=InFile[:-4]+"_"+str(thresholdPval)+"_grep.txt"
	OutF=open(OutFile,"w")
	for line in InF:
		[DR_Pval_chi2, AL_Pval_chi2, DR_Pval,AL_Pval]=line.split()[-4:]
		if float(DR_Pval)<thresholdPval:
			OutF.write(line)
		elif float(AL_Pval)<thresholdPval:
			OutF.write(line)
		elif float(DR_Pval_chi2)<thresholdPval:
			OutF.write(line)
		elif float(AL_Pval_chi2)<thresholdPval:
			OutF.write(line)
	InF.close()
	OutF.close()
	print "Done grepping the results. Output file is "+OutFile

print "FlyGWAS.py finished."
