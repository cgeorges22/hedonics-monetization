//  Hedonics
//  v2.12 CG (with github copilot in vscode) May 2026   for monetization 
//                add artificial and effective quality -- search by effective utilty 
//                (under development)
//  v2.11 CG      June 2019 for macro gentrification -- experiments with one time or gradual chage in single parameter
//	(v3 CG and students -- not included here -- developed Sept 2016-2017 for monte carlo experiments with various parameters) 
//	v2.10 CG		  June 2016	experiments with single changes in params such as eta during run
//	v2.9 CG				June 2016	allow bias toward own market for immiation at restart
//	v2.8 CG				May 2016	migration of firms across markets
//	v2.7 CG				May 2016	R&D by market as well as avg R&D performance
//	v2.6 CG				May 2016	track profit by market, markup in input.txt
//	v2.5 CG				May 2016	segment consumers and firms into two markets by class
//	v2.4 CG 			May 2016	account for wageBill, salaryBill, and production and OH labor separately
//  v2.3 CG				May 2016 	break out OH labor and salary -- prelim to two classes of consumer and firm
//	v2.2 CG				April 2016	add profit accounting recording
//	v2.1 CG				Nov 2015	adjust quality updating
//  v2.0 CG				Oct 2015	is v1.9 with bools removed from input.txt -- will work on all platforms (hpc, mac, gambs2)
//  v1.9 Philip Ewing	June 2015   break out parameter values to input.txt file
//       Chris Georges              two versions -- this one for mac and hpc, other for gambs2
//  v1.8 Chris Georges  June 2011   multiplicative quality increments 
//  v1.7 Chris Georges  Oct 2010    adds labor overhead, can suppress intermed goods
//                                  restart accounting fixed 6-3-11
//  v1.6 Chris Georges  July 2010   adds quasi discrete choice option for R&D choice 
//  v1.6 Chris Georges  March 2010  adds endogenous innovation -- adoption of R&D with labor cost
//  v1.5 Chris Georges  July 2009	adds loop over random seed and calculates means and variances of output
//                                  and corrects memory leak in consumer search aug 2 2009
//  v1.4 Chris Georges  June 2009 	adds complementarities
//  v1.3 Chris Georges  June 2009 	adds indivisibilities
//  v1.2 Chris Georges  June 2009   add quality and inflation adjusted output to v1.1    
//  v1.1 Chris Georges  Dec 2008    fixes product shares in v1.0
//  v1.0 Chris Georges  July 2008	alters CES Firm Dynamics model: firms1 devc++ v2.6 (see history below)
//	   purpose is to add hedonic product characterization and explicit 
//         model of product innovation

//  Firm Dynamics
//  v2.6 Chris Georges Jan 2008 
//  v2.0 Chris Georges March 2007
//  v1.0 Dane Johnson  Dec 2006
	
   
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "hedonics-firm-v2.12.h"
#include "hedonics-utilities-v2.12.h"
#include "time.h"
#include <cassert> //7/4/09
#include <cstring>

//for mac and hpc removed bools from input.txt
//but also removed the intermediate "saved_" versions of the input.txt variables and just directly input them in main()

//Global variables

int firmNum;
int firmNum1; //number of fims in low consumer segment (caters to production workers)
int firmNum2; //number of firms in high consumer segment (caters to OH labor)
int rounds;
int startOutput1; // start recording output
int outputFrequency; // record each xth round, 10 is every 10th round, 1 is every round
int printRoundsFrequency; // during simulation on std output 
int startRandDInvest; // updating choice of investment status
int startProductInnovation; // actual updating of hedonic product qualities
int startConsumerSearch; // rep consumer searches for new product mix
int startVolatilityMeasure; // for mean and volatility stats -- volatility as average absolute percent change from previous round
int volatilityMeasureFrequency; // this is freqency for sampling between startVolatilityMeasure and rounds
int startTechShocks;
int endTechShocks;
int endRandDInvest;
int endProductInnovation;
int endConsumerSearch;
int numSearchTrials;
int numMutationPairs;
double probInnovate;
double probUpdateZeroHQElement;
double probUpTick;
double probDownTick;
double probEffShock;
double probIntegerConstrained;
double probComplements;
double complementIntConstraintFactor;
double elasDem;
double markup; //5/16/16
double overHeadGReq;
double overHeadLReq;
double overHeadRandDReq;
double salary;
double wage;
double intPriceDiscount;
double probRandomRestart;
double incrementFactor;
double cesHedonicsElastGoods;
double cesHedonicsElastChars;
double riGain;
double rpGain;
double intensityOfChoice;
double probRandDInvestMutation;
double weightRandDDC; //5/28/16 weight placed on R&D vs. own market's profit in deciding R&D investment status
double avTechEff;
double strtTech;
double gdpDeflator;
double nomGDP;
double realGDP;
double wageBill;
double salaryBill; //5/12/16
double income; //5/13/16
int numRandD;
int numRandD1; //5/29/16 market 1
int numRandD2; //5/29/16 market 2
double medRecentRandDInvestment;  
double avgRecentProfitAboveMedRandD; 
double avgRecentProfitBelowMedRandD;
double avgProfit; 
double avgProfitType1; //market 1 caters to PL consumption //5/25/16
double avgProfitType2; //market 2 caters to OHL consumption
double avgRecentProfitType1;
double avgRecentProfitType2; 
double probSwitchMarkets; //5/29/16 at restart to higher profit market
double probRandomSwitchMarkets; //6/8/16 plus some purely random to keep n1 n2 from getting stuck at zero
int ownMarketBias; //6/10/16 //immitation - number of search trials to try to find firm to imitate in own market at restart
double testParam; //5/29/16
double artificialQualityDepreciationRate; // 2026-04-25 copilot: probability of decrementing positive artificial quality by 1 per round
int qualityConversionThreshold; // 2026-04-25 copilot: hedonic quality threshold above which conversion to artificial quality occurs
double qualityConversionRate; // 2026-04-25 copilot: probability of moving 1 unit from hedonic to artificial quality per round
bool doPostRunConsumerOptimization; // 2026-05-10 copilot: run consumer optimization routine at end of simulation
int optimRounds; // 2026-05-10 copilot: number of additional rounds for consumer optimization search

void getInput();
void applyArtificialHedonicShockToType2Subset(firm firms[], const int firmNum);
//copilot: 4/25/26 added function for depreciation of artificial quality and conversion from hedonic to artificial quality -- called in main loop below -- parameters added to input.txt
void doDepreciationAndConversion(firm firms[], const int firmNum, double depreciation_rate, int conversion_threshold, double conversion_rate, std::ofstream &testoutput, int round, bool anotation);

//copilot: 4/17/26 function for exogenous artificial quality shock to subset of type 2 firms in round 2000
//can change to type 1 firms as well as probability per hedonic element as desired  
//turn on/off on line 308 below 
void applyArtificialHedonicShockToType2Subset(firm firms[], const int firmNum) {
  int i, j;
  int type2Count = 0;
  int selectedType2Count = 0;

  for(i = 0; i < firmNum; i++) {
    if(firms[i].getType() == 2) {
      type2Count++;
    }
  }

  for(i = 0; i < firmNum && selectedType2Count < (type2Count / 2); i++) {
    if(firms[i].getType() != 2) {
      continue;
    }

    for(j = 0; j < numHedonicElements; j++) {
      if(randNumUni() < 0.2) { // 20 % chance to apply shock to each hedonic element
        firms[i].setArtificialQualityElementJ(j, firms[i].getArtificialQualityElementJ(j) + (rand() % 21)); // e.g rand() % 6 is rand integer between 0 and 5
      }
    }
    selectedType2Count++;
  }

  printf("applied artificial hedonic shock to %d of %d type-2 firms in round 2000\n", selectedType2Count, type2Count);
}

// Main Function
int main(int argc, char **argv) {
  
  //output files
  std::ofstream output1("data1.txt");		//output timm series of aggregate data
  std::ofstream output2("data2.txt");		//output last round firm distributions of size, productivity, etc.
  std::ofstream output3("data3.txt");		//output time series for six representative firms
  std::ofstream output4("data4.txt");		//output last round firm dist of hedonic quality vectors and final demand shares
  std::ofstream output4a("data4a.txt");		//copilot 2026-04-26: output first-round firm dist of hedonic quality vectors and final demand shares
  std::ofstream output4opt("data4opt.txt");	//copilot 2026-05-10: output consumer optimization results
  std::ofstream output5("data5.txt");		//output mean and volatility of Yc for each randSeed //7/6/09	
  std::ofstream testoutput("test.txt");	//track firm 6 if anotation switch is on
  
  //for looping over random seeds
  int randSeed; //7/5/09
  int randSeedStart = 10; //12, 14
  int randSeedEnd = 10; //12, 14
  //double mean[randSeedEnd - randSeedStart + 1]; //7/6/09 12/9/16 causing seg fault, move below getInput
  //double vol[randSeedEnd - randSeedStart + 1]; //7/6/09 12/9/16 causing seg fault, move below getInput
  
  //more declarations
  int i,j, round, restarts;
  double totOutput, totUtility, totUtilityPL, totUtilityOH, totEffectiveUtility, totEffectiveUtilityPL, totEffectiveUtilityOH, lastOutput, secs, mutScale, totProdLEmployment, totOHLEmployment; // 2026-04-15 copilot: track standard and effective utility separately.
  double finalWageBill, finalSalaryBill;
  bool commonShocks, imitation, anotation, altPriceIndex, cesHedonics, debugging;
  bool endogInnovation, discChoice, intermediateGoods, multiplicativeMutation;
  bool twoClasses; //5/12/16 alow two consumer/firm classes/submarkets
  bool endogN1N2; //5/29/16 endogenous switching of firms to high profit market
  bool artificialHedonicShockEnabled;
  time_t startTime, endTime;

  //initializaation
  // default values (may be overridden by input.txt in getInput)
  doPostRunConsumerOptimization = false; // disabled by default
  optimRounds = 1000; // default optimization rounds
  getInput(); // ints and reals from input.txt -- set bools below here
  
  double *mean = new double[randSeedEnd - randSeedStart + 1]; //7/6/09 12/9/16 causing seg fault, move below getInput
  double *vol = new double[randSeedEnd - randSeedStart + 1]; //7/6/09 12/9/16 causing seg fault, move below getInput
    
  cesHedonics = true;
  discChoice = true;
  endogInnovation = true;
  imitation = true; // for restarts -- if false, reset A_i and theta_i to 1.0 at restart
  multiplicativeMutation = true; // for hedonic innovation, if false then mutation is additive (characteristics incremented by +-1), if true increments are nearest integer proportional //6-9-11
  intermediateGoods = false; //10/6/10 //if this is false will bypass doIntermediateGoods, should also set overHeadG=0 above and set inputFactor=0 in firm.h manually -- though not necessary
  altPriceIndex = false; //if false, real GDP and gdpDeflators are standard (w/based year period 0). if true then use alt index that uses nom shares as weights in gdp deflator -- this removes hysteresis in y from idiosyncratic productivity shocks
  commonShocks = false; // heterogeneous vs. common shocks version
  anotation = false; // blow by blow description of activity of firm 4 to file test.txt
  debugging = false; // prints sequence of procedures //8/2/09
  twoClasses = true; // 5/12/16
  endogN1N2 = true; // 5/29/16 endogenous switching of firms to higher profit market
    
  if(twoClasses) {mutScale = (double) 1.0 / (double) firmNum;} //really should vary with relative firmNums ***
  else {mutScale = (double) 1.0 / (double) firmNum;} //1.0 //for consumer search -- upper limit of share mutations -- should be <1 and on order of 1/firmNum -- 2/firmNum if twoClasses
  firmNum1 = 0; firmNum2 = 0;
  
  //firmNum1 = ((int) firmNum/2); //5/11/16 for segmented consumer goods market by class *** not true if not 50/50 in initialization or if endog updating
  //firmNum2 = ((int) firmNum/2); //5/11/16 for segmented consumer goods market by class *** not true if not 50/50 in initialization or if endog updaing
  
  
  //for recording mean and volatility of output each run
  for(int i=0;i<randSeedEnd-randSeedStart+1;i++) {mean[i]=0;}
  for(int i=0;i<randSeedEnd-randSeedStart+1;i++) {vol[i]=0;}
  
  //top loop: loop over random seeds for multiple runs
  for(randSeed = randSeedStart; randSeed <= randSeedEnd; randSeed ++) { //7/5/09
    
    //set random seed for the current run
    srand(randSeed); 

    //removed all of phil's setting of input.txt values from "saved_" versions to true versions 7-11-15
     
    time(&startTime);
    
    //establish the array of firms -- allocate memory in heap  
    firm *firms;
    firms = new firm[firmNum];
    //firm firms[firmNum];  //this would place array in the stack 
    
    //set random seed for the current run //moved above 10-6-10
    //srand(randSeed); 
    
    //print some basic info
    printf("randSeed = %d \n", randSeed); //7/5/09
    printf("number of firms is %d \n", firmNum);  
    printf("size of firm is %lu B; size of population is %f MB \n", sizeof(firm), (double) firmNum * (double) sizeof(firm) /1000000 ); //%d changed to %ul (unsigned long) 6-14-15
    printf("intensity of choics is %f \n", intensityOfChoice); //7-12-15
    printf("weightRandDDC is %f \n", weightRandDDC); //5/28/16
    printf("endConsumerSearch is %d \n", endConsumerSearch); //5/28/16
    printf("testParam is %f \n", testParam); //5/28/16
	  printf("probRandDInvestMutation is %f \n", probRandDInvestMutation); //5/28/16
	  printf("probSwitchMarkets to higher profit market at restart is %f \n", probSwitchMarkets); //5/29/16
	  printf("probRandomSwitchMarkets at restart is %f \n", probRandomSwitchMarkets); //6/8/16
	  //printf("salary and markup are %f and %f \n", salary, markup); 6-17-19
	
	
	
    //testoutput << "number of rounds is " << rounds << "\n";
    //testoutput << "number of firms is " << firmNum << "\n";
    //testoutput << "size of firm is " << sizeof(firm) << " B; size of population is " << sizeof(firm)/1000000 << "\n";  
    
    //inititalize firms //moved to utilities 10-4-10
    //initializeFirms(firms, firmNum, round, testoutput, anotation, probComplements, probInnovate, probIntegerConstrained, intPriceDiscount);
    initializeFirms(firms, firmNum, twoClasses, testoutput, anotation, probComplements, probInnovate, probIntegerConstrained, intPriceDiscount, wage, salary, markup); // round removed 6-14-15 //markup added 5/26/16
      
    
    //output data on number of products with int constraints and complements
    //printf("number of products is %d \n", firmNum);
    int numWithIntConstraints = 0;
    for(i=0;i<firmNum;i++) {if(firms[i].getIntegerConstrained()) {++numWithIntConstraints;}}
    printf("number of products with integer constraints is %d \n", numWithIntConstraints);
    int numWithComplements = 0;
    for(i=0;i<firmNum;i++) {
      if(firms[i].getHasComplements()) {
		  ++numWithComplements;
		  //printf("does firm %d have complementary products? %d\n",i,firms[i].getHasComplements());
		  //printf("firm %d has %d complementary products\n", i, firms[i].getNumComplements());	
		  //printf("firm %d 's first complementary product is %d \n", i, firms[i].getComplementaryProductK(0)); //7/3/09
		  //printf("complementary hedonic element 1 for firm %d is %d \n", i, firms[i].getComplementaryHedonicElementKJ(0,0));
      }
    }
    printf("number of products with complements at start of simulation is %d \n", numWithComplements);
    
    //output data on first firm with complements
    i=0;
    while(! firms[i].getHasComplements() && i < firmNum) {i++;} 
    if(i<firmNum){
      printf("firm %d is the first firm with complements \n",i);
      printf("firm %d has %d complements \n",i,firms[i].getNumComplements());
      printf("firm %d 's first complementary product is %d \n", i, firms[i].getComplementaryProductK(0));
      //printf("at start of simulation firm %d's first complementary hedonic quality vector is (",i);
      //for(j=0;j<numHedonicElements -1;j++) {printf("%d,", firms[i].getComplementaryHedonicElementKJ(0,j));}
      printf("at start of simulation, firm %d's final demand share is %f \n", i,firms[i].getFinalDemandShare());  
    }
    int qqq = i; //3/19/10
    
    for(i=0;i<firmNum;i++) { //5/24/16
  		if(firms[i].getType()==1){firmNum1++;}
  		else {firmNum2++;}
    }
  	printf("The numbers of firms of types 1 and 2 are %d and %d \n", firmNum1, firmNum2); 
  
    
    //printf("the type of firm 2 is %d \n", firms[1].getType()); //5/11/16
    //printf("the type of firm 8 is %d \n", firms[7].getType()); //5/11/16
	  //printf("the wage is %f \n", wage); 
	  //printf("firm 6's wage is %f \n", firms[5].getWage());
	  //printf("the type of firm 12 is %d \n", firms[11].getType()); //5/11/16
	  printf("at start of simulation, firm 5 has markup %f \n", firms[4].getMarkup());//5/26/16 testing
    printf("at start of simulation, firm 5 has price %f \n", firms[4].getPrice());//5/26/16 testing
    printf("at start of simulation, firm 5 has wage %f \n", firms[4].getWage());//5/26/16 testing
    
	  restarts = 0;
    round = 0;
    // 2026-04-17 copilot: added switch to turn on artificial hedonic shock to subset of type 2 firms in round 2000 for testing.
    artificialHedonicShockEnabled = false; 

    //activity loop: loop over rounds in a run
    for(i = 0; i < rounds; i++) {
      round = i;
      if(anotation) {testoutput << "round " << round + 1 << ":" << "\n";}
      // 2026-04-15 copilot: reset per-round utility tracking before consumer search.
      totUtility = 0.0;
      totUtilityPL = 0.0;
      totUtilityOH = 0.0;
      totEffectiveUtility = 0.0;
      totEffectiveUtilityPL = 0.0;
      totEffectiveUtilityOH = 0.0;
      wageBill = 0.0;
      numRandD=0; //7/9/10
      numRandD1=0; //5/29/16
      numRandD2=0; //5/29/16
      medRecentRandDInvestment=0.0;  //4/6/16
	    avgRecentProfitAboveMedRandD=0.0; //4/6/16
      avgRecentProfitBelowMedRandD=0.0; //4/6/16
      avgProfit=0.0; //4/6/16
      avgProfitType1 = 0.0; //5/25/16
      avgProfitType2 = 0.0; 
      avgRecentProfitType1 = 0.0;
      avgRecentProfitType2 = 0.0;

      // 2026-04-17 copilot: apply artificial hedonic shock to subset of type 1 firms in round 2000.
      // if switch above on line 306 is set to true
      if(artificialHedonicShockEnabled && round == 2000) {
        applyArtificialHedonicShockToType2Subset(firms, firmNum);
      }
      
      //reset current profits to 0
      for(j = 0; j<firmNum; j++) {firms[j].setProfit(0.0);}
      
      /*
      //experiment -- e.g., change markup from 2 to 3 after round 2000 //6/10/16
      //or draw out over longer period //6-17-19
      //if change, may need to update this AND next chunk of code 
      //if((round + 1) == 2010) { //2-20-26 2000 for std experiments , other for testing specific firms
        //markup = 3.0; // 5, 10, 20 all work as well even though very extreem //2-25-26
        //salary = 5.0; //2-20-26  
        //markup = markup + 0.0004; //0.0004 good
      if(round >= 3000 && round < 4000)	{ //2-20-26
      //if(round >= 2000) { //2-20-26
	  	//markup = 3.0; // 5, 10, 20 all work as well even though very extreem //2-25-26
      //salary = 5.0; //2-20-26  
      markup = markup + 0.0004; //0.0004 good 
	  	//markup = markup + 0.0002;
      //markup = markup + 0.001; //6-17-19 2-20-26
	  	//salary = 2.0; //2-20-26
		  //salary = salary + 0.001; //2-20-26
      for(j=0; j<firmNum; j++){
		  firms[j].setMarkup(markup); //2-20-26
		  //firms[j].setSalary(salary); //2-20-26
		  //if(firms[j].getType()==2) { //i'm now inclined to drop this reset // 2-23-26
		  //	firms[j].setExptdFinalDemand(salary * firms[j].getExptdFinalDemand());
		  //}
		  //firms[72].setMarkup(1.5); //2-20-26 // note that this goes back to default when firm restarts
      } 
      }
      */
      // /*
        if((round + 1) % 1000 == 0){ //3-6-26 testing
                printf("in round %d, firm 5 has markup %f \n", round + 1, firms[4].getMarkup());//3/6/26 testing
                printf("in round %d, firm 5 has price %f \n", round + 1, firms[4].getPrice());//3/6/26 testing
                printf("in round %d, firm 5 has profit %f \n", round + 1, firms[4].getProfit());//3/6/26 testing
        }
      // */
      /*
      //if change experiment above, may need to updat this AND previous chunk of code
      //but i'm now inclined to drop this adjustment // 2-23-26
      //if(round >= 3000) {
      if(round == 3000) { //2/21/26
      //if(round >= 2000) { //2-20-26
      for(j=0; j<firmNum; j++){
        markup = 3.0; //2.4
        //salary = 2.0; //1.0 1.2
		  firms[j].setMarkup(markup);
		  //firms[j].setSalary(salary);
		  if(firms[j].getType()==2) {
		  	firms[j].setExptdFinalDemand(salary * firms[j].getExptdFinalDemand());
		  }
		  }   
      }
      */
        
      //R and D investment choice and accounting
      if(debugging){printf("update RandD investment decisions for round %d \n", round + 1);} //3/20/10
      if(endogInnovation && round >= startRandDInvest && round <= endRandDInvest) { 
	  	if(discChoice){updateRandDInvestmentStatusDC(firms, firmNum, intensityOfChoice, weightRandDDC, probRandDInvestMutation, intermediateGoods, overHeadLReq, overHeadRandDReq, salary, intPriceDiscount, medRecentRandDInvestment, avgRecentProfitAboveMedRandD, avgRecentProfitBelowMedRandD, avgRecentProfitType1, avgRecentProfitType2, testoutput, round, anotation);} //7/20/10 quasi discrete choice
	  	else {updateRandDInvestmentStatus(firms, firmNum, intensityOfChoice, probRandDInvestMutation, intermediateGoods, overHeadLReq, overHeadRandDReq, salary, intPriceDiscount, testoutput, round, anotation);} //3/22/10 individual social learning
      } //modified 8/30/10
      doOHandRandDAccounting(firms, firmNum, riGain,  overHeadLReq, overHeadRandDReq, salary, wageBill, salaryBill, numRandD, numRandD1, numRandD2, testoutput, round, anotation); //8/31/10 //moved out of updateR&DInvestmentStatus utility
      
      //depreciation of artificial quality and conversion from hedonic to artificial quality
      if(debugging){printf("do depreciation and conversion for round %d \n", round + 1);} //2026-04-25
      doDepreciationAndConversion(firms, firmNum, artificialQualityDepreciationRate, qualityConversionThreshold, qualityConversionRate, testoutput, round, anotation); //2026-04-25 copilot: turn on depreciation of artificial quality and conversion from hedonic to artificial quality
      //comment out lne above to turn off depreciation and conversion for testing

      //innovation
      if(debugging){printf("do product and process innovation for round %d \n", round + 1);} //8/2/09 
      if(commonShocks) {updateTechEffCommonShock(firms, firmNum, testoutput, round, avTechEff, strtTech, anotation);} //needs to be fleshed out
      else {
	  	if(round >= startTechShocks && round <= endTechShocks) { //6/02/09
	  		updateTechEff(firms, firmNum, testoutput, round, avTechEff, strtTech, anotation, probEffShock);
	  	}  //6/02/09
      	if(!endogInnovation && round >= startProductInnovation && round <= endProductInnovation) { //6/02/09
	  		updateHedonicQual(firms, firmNum, testoutput, round, probInnovate, probUpTick, probDownTick, anotation);
	  	} //6/02/09
	  	if(endogInnovation && round >= startProductInnovation && round <= endProductInnovation) { //3/20/10
	  		doEndogProductInnovation(firms, firmNum, testoutput, round, probInnovate, probUpdateZeroHQElement, probUpTick, probDownTick, multiplicativeMutation, incrementFactor, anotation); //3/22/10
	  	} //3/20/10 // multiplicativeMutation added 6-9-11
      }
      
      //prices, intermediate goods, and production
      if(debugging){printf("set prices for round %d \n", round + 1);} //8/2/09 
      updatePrice(firms, firmNum, /*inputFactor,*/ testoutput, anotation, intPriceDiscount); 
      if(debugging){printf("order produce and sell intermediate goods for round %d \n", round + 1);} //8/2/09 
      //doProduction(firms, wageBill, firmNum, totDeficit, defMax,  overHeadG, testoutput, round, anotation, intPriceDiscount); //split up 10-6-10
      doIntermediateGoods(firms, wageBill, firmNum, overHeadGReq, testoutput, round, anotation, intPriceDiscount, intermediateGoods); //10-6-10
      if(debugging){printf("produce consumption goods for round %d \n", round + 1);} //8/2/09 
      doProduction(firms, wageBill, firmNum, overHeadGReq, testoutput, round, anotation, intPriceDiscount, intermediateGoods); //10-6-10
      
      //consumer search and sales
      if(round >= startConsumerSearch && round <= endConsumerSearch) { //6/02/09  
		  if(debugging){printf("do consumer search for round %d \n", round + 1);} //8/2/09
		  if(twoClasses){
			totUtilityPL = doConsumerSearchClass1(firms, wageBill, firmNum, testoutput, round, numMutationPairs, numSearchTrials, mutScale, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor); 
			totUtilityOH = doConsumerSearchClass2(firms, salaryBill, firmNum, testoutput, round, numMutationPairs, numSearchTrials, mutScale, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor); 
      // 2026-04-15 copilot: compute effective utilities from the same final shares.
      double *sharesPL = new double[firmNum];
      double *sharesOH = new double[firmNum];
      for(j=0;j<firmNum;j++) {
        if(firms[j].getType()==1) {sharesPL[j] = firms[j].getFinalDemandShare();}
        else {sharesPL[j] = 0.0;}
        if(firms[j].getType()==2) {sharesOH[j] = firms[j].getFinalDemandShare();}
        else {sharesOH[j] = 0.0;}
      }
      totEffectiveUtilityPL = calcEffectiveUtility(firms, wageBill, firmNum, sharesPL, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
      totEffectiveUtilityOH = calcEffectiveUtility(firms, salaryBill, firmNum, sharesOH, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
      delete[] sharesPL;
      sharesPL = NULL;
      delete[] sharesOH;
      sharesOH = NULL;
		}
    else{
      totUtility = doConsumerSearch(firms, wageBill, salaryBill, firmNum, testoutput, round, numMutationPairs, numSearchTrials, mutScale, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
    // 2026-04-15 copilot: compute effective utility from the final single-class shares.
      double *sharesAll = new double[firmNum];
      for(j=0;j<firmNum;j++) {sharesAll[j] = firms[j].getFinalDemandShare();}
      totEffectiveUtility = calcEffectiveUtility(firms, wageBill + salaryBill, firmNum, sharesAll, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
      delete[] sharesAll;
      sharesAll = NULL;
    } //6/02/09 also return totUtility is new  
	  }

      if(debugging){printf("do final sales for round %d \n", round + 1);} //8/2/09
      totOutput = doFinalSales(firms, twoClasses, wageBill, salaryBill, firmNum, firmNum1, firmNum2, rpGain, testoutput, round, 
			       anotation, altPriceIndex, nomGDP, realGDP, gdpDeflator, avgProfit, avgProfitType1, avgProfitType2); //6/8/09 add nomGDP, realGDP, gdpDeflator 
	  
	  //calc employment
	  totProdLEmployment = wageBill/wage; //5/13/16
	  totOHLEmployment = salaryBill/salary; //5/13/16
	  if(round == 1000) {
		printf("total output and R&D are %f and %d \n", totOutput, numRandD);
		printf("total employment of production labor and OH labor are %f and %f \n", totProdLEmployment , totOHLEmployment);
    printf("total wageBill and salaryBill are %f and %f \n", wageBill, salaryBill);
    printf("total Utility of PL and OHL are %f and %f \n", totUtilityPL, totUtilityOH);
    printf("salary and markup are %f and %f \n", salary, markup); //6-17-19
  	}
      
    //calc mean and volatility of output
    if(debugging){printf("calc volatility for round %d \n", round + 1);} //8/2/09
    if(round == startVolatilityMeasure - volatilityMeasureFrequency) {lastOutput = totOutput;}			
    if(round>=startVolatilityMeasure && round % volatilityMeasureFrequency == 0){	//7/6/09
		  mean[randSeed - randSeedStart] += totOutput;
		  vol[randSeed - randSeedStart] += fabs(log(totOutput) - log(lastOutput))*100;
		  //if(round % 100 == 0 ) printf("abs percent change Yc in round %d is %f \n", round, /*fabs(log(totOutput) - log(lastOutput))*100*/ totOutput - lastOutput);
		  lastOutput = totOutput; //7/6/09
      } 
    if(round==rounds-1){ //7/6/09
		mean[randSeed - randSeedStart] = (mean[randSeed - randSeedStart]/((double)rounds - (double) startVolatilityMeasure))*((double) volatilityMeasureFrequency);
		vol[randSeed - randSeedStart] = (vol[randSeed - randSeedStart]/((double)rounds - (double) startVolatilityMeasure))*((double) volatilityMeasureFrequency);
		//printf("volatility of Tc is %f \n", vol[randSeed-randSeedStart] );
		//printf("mean of Yc is %f \n", mean[randSeed-randSeedStart]);
    } 
      
    //printf("RandD overall and in markets 1 and 2 are %d %d and %d \n", numRandD, numRandD1, numRandD2 );
      
    //output main time series to output1 (for first randSeed only) (and before restarts)
    if(round >= startOutput1 && round % outputFrequency == 0 && randSeed == randSeedStart){  //6/02/09 added totUtility, removed avQual 6/8/09 add nomGDP realGDP gdpdeflator
		if(debugging){printf("output1 for round %d \n", round + 1);} //8/2/09 
		output1 << round << " " << totOutput << " " << nomGDP << " " 
			<< realGDP << " " << gdpDeflator << " " << firmNum1 << " " << firmNum2 << " " 
			<< totUtility << " " << totUtilityPL << " " << totUtilityOH  << " "
			<< avTechEff << " " << strtTech << " " << wageBill << " " 
			<< salaryBill << " " << totProdLEmployment << " " << totOHLEmployment << " "
			<< restarts << " " << numRandD << " " << numRandD1 << " " << numRandD2 << " "
			<< avgProfit << " " << avgRecentProfitAboveMedRandD << " " << avgRecentProfitBelowMedRandD <<  " " 
			<< avgRecentProfitType1 << " " << avgRecentProfitType2 << "\n";
      }
      //output sixfirm time series to output3 (for first randSeed only)
      if(round >= startOutput1 && round % outputFrequency == 0 && randSeed == randSeedStart){	
		if(debugging){printf("output3 for round %d \n", round + 1);} //8/2/09
		output3 << round << " " << (firms[3]).getFinalProduction() << " " << (firms[3]).getExptdFinalDemand()/(firms[3]).getPrice() << " " << (firms[3]).getCapital() << " ";
		output3 << (firms[4]).getFinalProduction() << " " << (firms[4]).getExptdFinalDemand()/(firms[4]).getPrice() << " " << (firms[4]).getCapital() << " ";
		output3 << (firms[5]).getFinalProduction() << " " << (firms[5]).getExptdFinalDemand()/(firms[5]).getPrice() << " " << (firms[5]).getCapital() << " ";
		output3 << (firms[6]).getFinalProduction() << " " << (firms[6]).getExptdFinalDemand()/(firms[6]).getPrice() << " " << (firms[6]).getCapital() << " ";
		output3 << (firms[7]).getFinalProduction() << " " << (firms[7]).getExptdFinalDemand()/(firms[7]).getPrice() << " " << (firms[7]).getCapital() << " ";
		output3 << (firms[8]).getFinalProduction() << " " << (firms[8]).getExptdFinalDemand()/(firms[8]).getPrice() << " " << (firms[8]).getCapital() << "\n";
      }
      //testoutput (for first randSeed only)
      if(round == 3 && randSeed == randSeedStart) {
		if(debugging){printf("testoutput for round %d \n", round + 1);} //8/2/09 
		testoutput <<  "  final production: " << (firms[0]).getFinalProduction() << " " << (firms[1]).getFinalProduction() << " " 
		   << (firms[2]).getFinalProduction() << " " <<  (firms[3]).getFinalProduction() << " " 
		   << (firms[4]).getFinalProduction() << " " <<  (firms[5]).getFinalProduction() << " " 
		   << (firms[6]).getFinalProduction() << " " <<  (firms[7]).getFinalProduction() << " " 
		   << (firms[8]).getFinalProduction() << " "  << (firms[9]).getFinalProduction() << "\n" ;
      }
      
      //restarts 
      if(debugging){printf("do restarts for round %d \n", round + 1);} //8/2/09 
      if(commonShocks) {restarts = doRestartsCommonShock(firms, firmNum, strtTech, wageBill, overHeadGReq, intPriceDiscount, imitation, probRandomRestart, wage, salary, markup);}
      else {
	  	if (endogN1N2) {
		  restarts = doRestartsHedonicsEndogN1N2(firms, firmNum, firmNum1, firmNum2, strtTech, wageBill, salaryBill, overHeadGReq, intPriceDiscount, imitation, probRandomRestart, probIntegerConstrained, probComplements, wage, salary, markup, avgRecentProfitType1, avgRecentProfitType2, probSwitchMarkets, probRandomSwitchMarkets, ownMarketBias);
		}
      	else {
		  restarts = doRestartsHedonics(firms, firmNum, strtTech, wageBill, salaryBill, overHeadGReq, intPriceDiscount, imitation, probRandomRestart, probIntegerConstrained, probComplements, wage, salary, markup);
		}
      }

      finalWageBill = wageBill;
      finalSalaryBill = salaryBill;

      //copilot 2026-04-26: output first-round firm dist of hedonic quality vectors and final demand shares (for first randSeed only)
      if(round == 0 && randSeed == randSeedStart){
        if(twoClasses) {
          output4a << "effectiveUtilityPL " << totEffectiveUtilityPL << "\n";
          output4a << "effectiveUtilityOH " << totEffectiveUtilityOH << "\n";
        }
        else {
          output4a << "effectiveUtility " << totEffectiveUtility << "\n";
        }
        for(j = 0; j < firmNum; j++) {
          for(int k = 0; k < numHedonicElements; k++) {output4a << (firms[j]).getHedonicQualityElementJ(k) << " ";}
          output4a << "\n";
          for(int k = 0; k < numHedonicElements; k++) {output4a << (firms[j]).getArtificialQualityElementJ(k) << " ";}
          output4a << "\n";
          for(int k = 0; k < numHedonicElements; k++) {output4a << (firms[j]).getEffectiveQualityElementJ(k) << " ";}
          output4a << "  " << (firms[j]).getFinalDemandShare() << "\n";
          output4a << " " << (firms[j]).getType() << "\n";
        }
      }

      wageBill = 0.0; //fixed from == 8/2/09
      salaryBill = 0.0; // 5/13/16

      
      //if((i+1) >= 2000){  //3-7-26 testing -- print round by round -- set printRoundsFrequency to 1 or 2 for this
      if ((i+1) % printRoundsFrequency == 0) printf("Completed round %d\n", i + 1); 
      if ((i+1) % printRoundsFrequency == 0) printf("salary and markup are %f and %f \n", salary, markup); //6-17-19
      if ((i+1) % printRoundsFrequency == 0) printf("avgRecentProfitType1 and avgRecentProfitType2 are %f and %f \n ", avgRecentProfitType1, avgRecentProfitType2);
	
      //if ((i+1) % printRoundsFrequency == 0) printf("markup of firm 73 is %f\n", firms[72].getMarkup()); //3-6-26 testing
      //if ((i+1) % printRoundsFrequency == 0) printf("price of firm 73 is %f\n", firms[72].getPrice()); //3-6-26 testing
      //if ((i+1) % printRoundsFrequency == 0) printf("profit of firm 73 is %f\n", firms[72].getProfit()); //3-6-26 testing
      //if ((i+1) % printRoundsFrequency == 0) printf("share of firm 73 %f\n", firms[72].getFinalDemandShare()); //3-6-26 testing
      //if ((i+1) % printRoundsFrequency == 0) printf("type of firm 73 is %d\n", firms[72].getType()); //3-7-26 testing
      //} //3-7-26 testing
      

    }//end of activity loop
    
    //output last round firm distributions to output2 (for first randSeed only)
    if(randSeed == randSeedStart){
      if(debugging){printf("output2 for final round \n");} //8/2/09
      for(i = 0; i < firmNum; i++) {
	    output2 << (firms[i]).getCapital() << " ";
	    if((firms[i]).getExptdFinalDemand()/(firms[i]).getPrice() >= (firms[i]).getFinalProduction() ) { //output quantity sold 
	    output2 << (firms[i]).getFinalProduction() << " ";
	    output2 << (firms[i]).getFinalProduction() * (firms[i]).getPrice() << " "; //add nominal sales //7/8/09
	}
	else {
	  output2 << (firms[i]).getExptdFinalDemand()/(firms[i]).getPrice() << " ";
	  output2 << (firms[i]).getExptdFinalDemand() << " "; //add nominal sales //7/8/09
	}
	output2 << (firms[i]).getFinalProduction() << " "; //ouput quantity produced
	output2 << (firms[i]).getExptdFinalDemand()/(firms[i]).getPrice() << " "; //output quantity demanded 
	output2 <<  (firms[i]).getTechEff() << " " <<  /*(firms[i]).getQual() <<*/ "\n";
      } 
    }
    
    //output last round firm dist of hedonic quality vectors and final demand shares (for first randSeed only) (and after restarts)
    if(randSeed == randSeedStart){
      if(debugging){printf("output4 for last round \n");} //8/2/09
      // 2026-04-15 copilot: write effective utility values before the firm-level arrays.
      if(twoClasses) {
    output4 << "effectiveUtilityPL " << totEffectiveUtilityPL << "\n";
    output4 << "effectiveUtilityOH " << totEffectiveUtilityOH << "\n";
    }
    else {
    output4 << "effectiveUtility " << totEffectiveUtility << "\n";
    }
    for(i = 0; i < firmNum; i++) {
		// 2026-04-15 copilot: output standard and effective hedonic arrays back-to-back for each firm.
		for(j = 0; j < numHedonicElements; j++) {output4 << (firms[i]).getHedonicQualityElementJ(j) << " ";}
    output4 << "\n";
    for(j = 0; j < numHedonicElements; j++) {output4 << (firms[i]).getArtificialQualityElementJ(j) << " ";}
    output4 << "\n"; //cg 2026-04-24 add artificial hedonic array
    for(j = 0; j < numHedonicElements; j++) {output4 << (firms[i]).getEffectiveQualityElementJ(j) << " ";}
    output4 << "  " << (firms[i]).getFinalDemandShare() << "\n";
		output4 << " " << (firms[i]).getType() << "\n"; //5/18/16
	//printf("at end of simulation, firm %d has innovation status %d \n",i+1,firms[i].getRandDInvestment());//3/24/10 testing
      }
    }
    
    //print end of simulation Yc stats, and demand share for first firm with complementscd
    if(qqq < firmNum) {printf("at end of simulation, firm %d's final demand share is %f \n", qqq,firms[qqq].getFinalDemandShare());}
    printf("mean of Yc is %f \n", mean[randSeed-randSeedStart]); //7/6/09
    printf("volatility of Yc is %f \n", vol[randSeed-randSeedStart]); //7/6/09	
    printf("final utility of prodution and overhead workers are %f and %f \n", totUtilityPL, totUtilityOH);
    printf("final numbers of firms in markets 1 and 2 are %d and %d \n ", firmNum1, firmNum2);
    
    //for all randSeeds, output mean and volatility of Yc to output5 -- volatility is average absolute % change from volatilityMeasureFrequency periods ago
    output5 << randSeed << " " << mean[randSeed-randSeedStart] << " " << vol[randSeed-randSeedStart] << "\n"; //7/6/09
    
    time(&endTime);
    secs = difftime(endTime, startTime);
    printf("Done Running Simulation. \n" ); 
    printf("Simulation took %f minutes. \n", secs/60);
    
    //testoutput (for first randSeed only)
    if(randSeed == randSeedStart){testoutput << "simulation took " << secs/60 << " minutes " << "\n";}
    
    //Run consumer post run optimization routine if enabled (for first randSeed only) //copilot 2026-05-10
    if(doPostRunConsumerOptimization && randSeed == randSeedStart) {
      if(debugging){printf("Running consumer optimization routine \n");}

      //if post-run optimization enabled, indicate this 
      printf("Running consumer optimization routine for %d rounds per class \n", optimRounds);
      
      // Output initial state (end of regular simulation)
      output4opt << "INITIAL STATE (End of Round " << rounds - 1 << ")\n";
      if(twoClasses) {
        output4opt << "PL_utility " << totUtilityPL << " PL_effective_utility " << totEffectiveUtilityPL << " PL_income " << finalWageBill << "\n";
        output4opt << "OH_utility " << totUtilityOH << " OH_effective_utility " << totEffectiveUtilityOH << " OH_income " << finalSalaryBill << "\n";
      } else {
        output4opt << "Utility " << totUtility << " Effective_utility " << totEffectiveUtility << " Income " << (finalWageBill + finalSalaryBill) << "\n";
      }
      output4opt << "Firm_shares: ";
      for(i=0; i<firmNum; i++) {output4opt << firms[i].getFinalDemandShare() << " ";}
      output4opt << "\n";
      
      // Run optimization for each consumer class
      double optUtilityPL = 0.0, optUtilityOH = 0.0;
      double optEffUtilityPL = 0.0, optEffUtilityOH = 0.0;
      
      if(twoClasses) {
        optUtilityPL = doPostRunConsumerOptimizationClass1(firms, finalWageBill, firmNum, optimRounds, output4opt, numMutationPairs, numSearchTrials, mutScale, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
        optUtilityOH = doPostRunConsumerOptimizationClass2(firms, finalSalaryBill, firmNum, optimRounds, output4opt, numMutationPairs, numSearchTrials, mutScale, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
        
        // Calculate effective utilities at optimized shares
        double *sharesPL = new double[firmNum];
        double *sharesOH = new double[firmNum];
        for(j=0;j<firmNum;j++) {
          if(firms[j].getType()==1) {sharesPL[j] = firms[j].getFinalDemandShare();}
          else {sharesPL[j] = 0.0;}
          if(firms[j].getType()==2) {sharesOH[j] = firms[j].getFinalDemandShare();}
          else {sharesOH[j] = 0.0;}
        }
        optEffUtilityPL = calcEffectiveUtility(firms, finalWageBill, firmNum, sharesPL, testoutput, rounds, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
        optEffUtilityOH = calcEffectiveUtility(firms, finalSalaryBill, firmNum, sharesOH, testoutput, rounds, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);
        delete[] sharesPL;
        delete[] sharesOH;
      }
      
      // Output post run optimized state
      output4opt << "AFTER POST-RUN OPTIMIZATION (After " << optimRounds << " optimization rounds per class)\n";
      if(twoClasses) {
        output4opt << "PL_utility " << optUtilityPL << " PL_effective_utility " << optEffUtilityPL << " PL_income " << finalWageBill << "\n";
        output4opt << "OH_utility " << optUtilityOH << " OH_effective_utility " << optEffUtilityOH << " OH_income " << finalSalaryBill << "\n";
      } else {
        output4opt << "Utility " << (optUtilityPL + optUtilityOH) << " Effective_utility " << (optEffUtilityPL + optEffUtilityOH) << " Income " << (finalWageBill + finalSalaryBill) << "\n";
      }
      output4opt << "Firm_shares: ";
      for(i=0; i<firmNum; i++) {output4opt << firms[i].getFinalDemandShare() << " ";}
      output4opt << "\n";
      
      // Output final firm qualities
      output4opt << "FINAL FIRM QUALITIES AFTER POST-RUN OPTIMIZATION\n";
      for(i = 0; i < firmNum; i++) {
        for(j = 0; j < numHedonicElements; j++) {output4opt << (firms[i]).getHedonicQualityElementJ(j) << " ";}
        output4opt << "\n";
        for(j = 0; j < numHedonicElements; j++) {output4opt << (firms[i]).getArtificialQualityElementJ(j) << " ";}
        output4opt << "\n";
        for(j = 0; j < numHedonicElements; j++) {output4opt << (firms[i]).getEffectiveQualityElementJ(j) << " ";}
        output4opt << "\n";
      }
      
    }
    // End of post-run optimization routine - record time taken in total //added 5/22/26
    time(&endTime);
    secs = difftime(endTime, startTime);
    printf("Done Running Post Run Optimization. \n" ); 
    printf("Full program took %f minutes. \n", secs/60);


    delete[] firms; //7/6/09 
    firms = NULL; //7/6/09
    
  } //end for randSeed loop
  
 delete[] mean; //12/9/16
 delete[] vol;  //12/9/16
  
    //close output files at end of all runs
  output1.close();
  output2.close();
  output3.close();
  output4.close(); //7/5/09
  output5.close(); //7/6/09
  output4a.close(); //copilot 2026-04-26
  output4opt.close(); //copilot 2026-05-10
  testoutput.close();
  
  //system("PAUSE"); //only for windows
  
  
  return 0;
  
} //end main

void toBool(bool & var, char * val);

void getInput() {
  // Strings to store line by line, variable name, start, and final value
  char line[256];
  char* variable;
  char* value;
  
  // Uses the file "input.txt"
  std::fstream input;
  input.open("input.txt", std::fstream::in);
  
  // Go through the input file line by line, tokenize the line into three strings
  input.getline(line,256);
  while (strcmp(line, "") != 0) {
    
    variable = strtok(line, " ");
    // start = strtok(NULL, " ");
    value = strtok(NULL, " ");
    
    // cout << variable << " " << start << " " << last << endl;
    // outfile <<  variable << " " << start << " " << last << endl;
    
    // Find and set the given variable
    
     if (strcmp(variable, "probEffShock") == 0) {
     probEffShock = atof(value);
     }
     else if (strcmp(variable, "probIntegerConstrained") == 0) {
     probIntegerConstrained = atof(value);
     }
     else if (strcmp(variable, "probComplements") == 0) {
     probComplements = atof(value);
     }
     else if (strcmp(variable, "complementIntConstraintFactor") == 0) { //6/6/12
     complementIntConstraintFactor = atof(value);
     }
     else if (strcmp(variable, "numSearchTrials") == 0) {  //6/6/12
     numSearchTrials = atoi(value);
     }
     else if (strcmp(variable, "numMutationPairs") == 0) {
     numMutationPairs = atoi(value);
     }
     else if (strcmp(variable, "probInnovate") == 0) {
     probInnovate = atof(value);
     }
     else if (strcmp(variable, "probUpdateZeroHQElement") == 0) {
     probUpdateZeroHQElement = atof(value);
     } 
     else if (strcmp(variable, "probUpTick") == 0) {
     probUpTick = atof(value);
     }
     else if (strcmp(variable, "probDownTick") == 0) {
     probDownTick = atof(value);
     }
     else if (strcmp(variable, "markup") == 0) {
     markup = atof(value);
     }
     else if (strcmp(variable, "overHeadGReq") == 0) {
     overHeadGReq = atof(value);
     }
     else if (strcmp(variable, "overHeadLReq") == 0) {
     overHeadLReq = atof(value);
     }
     else if (strcmp(variable, "overHeadRandDReq") == 0) {
     overHeadRandDReq = atof(value);
     }
     else if (strcmp(variable, "salary") == 0) {
     salary = atof(value);
	 }
	 else if (strcmp(variable, "wage") == 0) {
     wage = atof(value);
	 }
     else if (strcmp(variable, "intPriceDiscount") == 0) {
     intPriceDiscount = atof(value);
     }
     else if (strcmp(variable, "probRandomRestart") == 0) {
     probRandomRestart = atof(value);
     }
     else if (strcmp(variable, "incrementFactor") == 0) {
     incrementFactor = atof(value);
     }
     else if (strcmp(variable, "cesHedonicsElastGoods") == 0) {
     cesHedonicsElastGoods = atof(value);
     }
     else if (strcmp(variable, "cesHedonicsElastChars") == 0) {
         cesHedonicsElastChars = atof(value);
     }
     else if (strcmp(variable, "riGain") == 0) {
         riGain = atof(value);
     }
     else if (strcmp(variable, "rpGain") == 0) {
         rpGain = atof(value);
     }
     else if (strcmp(variable, "intensityOfChoice") == 0) {
         intensityOfChoice = atof(value);
     }
     else if (strcmp(variable, "probRandDInvestMutation") == 0) {
         probRandDInvestMutation = atof(value);
     }
     else if (strcmp(variable, "avTechEff") == 0) {
         avTechEff = atof(value);
     }
     else if (strcmp(variable,"strtTech") == 0) {
         strtTech = atof(value);
     }
     /*
     else if (strcmp(variable, "gdpDeflator") == 0) {
         gdpDeflator = atof(value);
     }
     else if (strcmp(variable, "nomGDP") == 0) {
         nomGDP = atof(value);
     }
     else if (strcmp(variable, "realGDP") == 0) {
         realGDP = atof(value);
     }
     else if (strcmp(variable, "wageBill") == 0) {
         wageBill = atof(value);
     }
     else if (strcmp(variable, "salaryBill") == 0) {
         salaryBill = atof(value);
     }
     else if (strcmp(variable, "numRandD") == 0) {
         numRandD = atoi(value);
     }
     */
     else if (strcmp(variable, "firmNum") == 0) {
         firmNum = atoi(value);
     }
     else if (strcmp(variable, "rounds") == 0) {
         rounds = atoi(value);
     }
     else if (strcmp(variable, "startOutput1") == 0) {
         startOutput1 = atoi(value);
     }
     else if (strcmp(variable, "outputFrequency") == 0) {
         outputFrequency = atoi(value);
     }
     else if (strcmp(variable, "printRoundsFrequency") == 0) {
         printRoundsFrequency = atoi(value);
     }
     else if (strcmp(variable, "startRandDInvest") == 0) {
         startRandDInvest = atoi(value);
     }
     else if (strcmp(variable, "startProductInnovation") == 0) {
         startProductInnovation = atoi(value);
     }
     else if (strcmp(variable, "startConsumerSearch") == 0) {
         startConsumerSearch = atoi(value);
     }
     else if (strcmp(variable, "startVolatilityMeasure") == 0) {
         startVolatilityMeasure = atoi(value);
     }
     else if (strcmp(variable, "volatilityMeasureFrequency") == 0) {
         volatilityMeasureFrequency = atoi(value);
     }
     else if (strcmp(variable, "startTechShocks") == 0) {
         startTechShocks = atoi(value);
     }
     else if (strcmp(variable, "endTechShocks") == 0) {
         endTechShocks = atoi(value);
     }
     else if (strcmp(variable, "endRandDInvest") == 0) {
         endRandDInvest = atoi(value);
     }
     else if (strcmp(variable, "endProductInnovation") == 0) {
         endProductInnovation = atoi(value);
     }
	 else if (strcmp(variable, "weightRandDDC") == 0) {
     	 weightRandDDC = atof(value);  
	 }
     else if (strcmp(variable, "endConsumerSearch") == 0) {
         endConsumerSearch = atoi(value);
     }
     else if (strcmp(variable, "probSwitchMarkets") == 0) {
         probSwitchMarkets = atof(value);
     }
     else if (strcmp(variable, "probRandomSwitchMarkets") == 0) {
         probRandomSwitchMarkets = atof(value);
     }
     else if (strcmp(variable, "ownMarketBias") == 0) {
         ownMarketBias = atoi(value);
     }
     else if (strcmp(variable, "testParam") == 0) {
         testParam = atof(value);
     }
     else if (strcmp(variable, "artificialQualityDepreciationRate") == 0) {
         artificialQualityDepreciationRate = atof(value);
     }
     else if (strcmp(variable, "qualityConversionThreshold") == 0) {
         qualityConversionThreshold = atoi(value);
     }
     else if (strcmp(variable, "qualityConversionRate") == 0) {
         qualityConversionRate = atof(value);
     }
     else if (strcmp(variable, "doPostRunConsumerOptimization") == 0) {
      // Accept numeric or textual true values; handle possible CRLF or trailing chars
      doPostRunConsumerOptimization = (atoi(value) == 1 || strcmp(value, "true") == 0 || strcmp(value, "True") == 0);
     }
     else if (strcmp(variable, "optimRounds") == 0) {
         optimRounds = atoi(value);
     }
     
    else {
      std::cout << "Variable name: " << variable << " not recognized. No value set." << std::endl;
    }  		
    input.getline(line, 256);
  }
  input.close();
}

void toBool(bool & var, char * val) {
  if (strcmp(val, "true") == 0)
    var = true;
  else
    var = false;
}
