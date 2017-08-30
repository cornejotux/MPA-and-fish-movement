#include "fish.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm> 
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <vector>

#if MPI_ON  == 1
       #include <mpi.h>
#endif
#define test false

//Every function is from fish.cpp

// How to call gridtest 
//
// $> ./gridtest -MPA 20 -SST -Int -Cur -f 20 -tga -Name 1

int main(int argc, char** argv)
{
    //printf("----> %i %s \n", argc, args[1]);]
    if (argc == 1)
    {
        printf("\n -------------------------------------\n | You NEED to pass a MPA size!!     |\n | e.g $$ ./gristest -mpa 25 -tac 10 |\n -------------------------------------\n\n");
        
        exit(1);
    }
    int inputTAC, inputMPA, centerX, centerY;
    int lengthX = 0;
    int lengthY = 0;
    int YearOfTAC = 150;
    int YearOfMPA = 150;
    int LagOfManagement = 10;
    int resTime = 0;
    int UseMultiSST = 0;
    int nEggs = 10;
    float F_Mortality = 0.0022; // correspond to 1.23 [1/y]
    float f2;
    char _F[] = "NoF";
    char manag[] = "NoMPA";
    char _managX[] = "NaN";
    char _managY[] = "NaN";
    double MPAsize;
    double sMPA = 0;
    double minMPAx = 0;
    double maxMPAx = 0;
    double minMPAy = 0;
    double maxMPAy = 0;
    double boatsC = 3;
    bool USE_Interactions = false;
    bool USETEMPFIELD = false;
    bool USECURRENTS;
    bool mpaFlag = false;
    bool isMPA = false;
    bool tacFlag = false ;
    bool tga = false;
    char sst[] = "NoSST";
    char cur[] = "NoCur";
    char interact[] = "NaN";
    char RunName[] = "NoName";

    if (argc > 1)
    {
    for (int i=1; i<(argc); i++)
    {
        if (strncmp(argv[i], "-lag", 4) == 0 or strncmp(argv[i], "-LAG", 4) == 0)
	{
	    LagOfManagement = (double)strtoul(argv[i+1], 0, 0);
	}
        else if (strncmp(argv[i], "-TAC", 4) == 0 or strncmp(argv[i], "-tac", 4) == 0)
        {
            inputTAC = (int)strtoul(argv[i+1], 0, 0);
            printf("%s, %i \n",  argv[i],  inputTAC);
            sprintf(manag, "TAC");
            tacFlag = true;
	    YearOfTAC = YearOfFishing + LagOfManagement;
        } 
        else if (strncmp(argv[i], "-MPA",4) == 0 or strncmp(argv[i], "-mpa",4) == 0)
        {
            inputMPA=(int)strtoul(argv[i+1], 0, 0); 
            printf("%s, %i \n",  argv[i],  inputMPA);
            sprintf(manag, "MPA");
            MPAsize = inputMPA;
            sMPA = MPAsize /2;
            minMPAx = 0;   // min is 0: This is left in the image
            maxMPAx = 100; // max is 83: This is right in the image
            minMPAy = 50-sMPA; // Min is 0: This is the down in the image
            maxMPAy = 50+sMPA;
            mpaFlag = true;
            isMPA = true;
            YearOfMPA = YearOfFishing + LagOfManagement;
	    printf("Year of MPA starts %i\n\n", YearOfMPA); 
       }
        else if (strncmp(argv[i], "-x",2) == 0 or strncmp(argv[i], "-X",2) == 0)
        {
            centerX=(int)strtoul(argv[i+1], 0, 0);
            lengthX=(int)strtoul(argv[i+2], 0, 0);
            printf("%s, %i \n",  argv[i],  centerX);
            minMPAx = centerX - lengthX;   // min is 0: This is left in the image
            maxMPAx = centerX + lengthX; // max is 83: This is right in the image
            sprintf(_managX, "%.f-%.f", minMPAx, maxMPAx);
            sprintf(manag, "MPA");
            isMPA = true;
	    YearOfMPA = YearOfFishing + LagOfManagement;
        }
        else if (strncmp(argv[i], "-y",2) == 0 or strncmp(argv[i], "-y",2) == 0)
        {
            centerY=(int)strtoul(argv[i+1], 0, 0);
            lengthY=(int)strtoul(argv[i+2], 0, 0);
            printf("%s, %i \n",  argv[i],  centerY);
            minMPAy = centerY - lengthY; // Min is 0: This is the down in the image
            maxMPAy = centerY + lengthY;
            sprintf(_managY, "%.f-%.f", minMPAy, maxMPAy);
        }
        else if (strncmp(argv[i], "-Int",4) == 0 or strncmp(argv[i], "-int",4) == 0)
        {
        USE_Interactions = true;
        sprintf(interact, "Int");
        }
        else if (strncmp(argv[i], "-noInt",6) == 0 or strncmp(argv[i], "-NoInt",6) == 0)
        {
        USE_Interactions = false;
        sprintf(interact, "NoInt");
        }
        else if (strncmp(argv[i], "-SST",4) == 0 or strncmp(argv[i], "-sst",4) == 0)
        {
        USETEMPFIELD = true;
        sprintf(sst, "SST");
        }
        else if ((strncmp(argv[i], "-NoSST",6) == 0) or (strncmp(argv[i], "-noSST",6) == 0))
        {
        USETEMPFIELD = false;  
        sprintf(sst, "NoSST");
        }
        else if (strncmp(argv[i], "-Cur",4) == 0 or strncmp(argv[i], "-cur",4) == 0)
        {
        USECURRENTS = true;
        sprintf(cur, "Cur");
        }
        else if (strncmp(argv[i], "-NoCur",6) == 0 or strncmp(argv[i], "-noCur",6) == 0)
        {
        USECURRENTS = false;
        sprintf(cur, "NoCur");
        }
        else if (strncmp(argv[i], "-f", 2) == 0 or strncmp(argv[i], "-F", 2) == 0)
        {
            double f = (double)strtoul(argv[i+1], 0, 0);
            float f_ = f;
            f2= f/100;
            F_Mortality = pow(f2, 0.002739726)-1;
            sprintf(_F, "%4.2f", (f2-1));
            printf("Year rate: %s, and as Day Rate %6.6f  \n", _F, F_Mortality);
        }
        else if (strncmp(argv[i], "-tga", 4) == 0)
        {
            tga = true;
        }
	else if (strncmp(argv[i], "-C", 4) == 0) 
	{
            boatsC = (double)strtoul(argv[i+1], 0, 0);
	}
       	else if (strncmp(argv[i], "-egg", 4) == 0)
	{
		nEggs = (double)strtoul(argv[i+1], 0, 0);
	}
        // This is to automatically change the names of the output files, 
        // the idea with this is to be able to run a Montecarlo analysis
        else if (strncmp(argv[i], "-Name", 5) == 0)
        {
            char _RunName_ = (char)strtoul(argv[i+1], 0, 0);
            sprintf(RunName, "%i", _RunName_);

            printf("File Begins with: %s_  \n", RunName);
        }
	else if (strncmp(argv[i], "-MSST", 5) == 0)
	{
		UseMultiSST = 1;
	}
	if (isMPA and tacFlag)
	{
            sprintf(manag, "TAC-MPA");
	}        
    }
    }
    
    if (strncmp(interact, "NaN",3) == 0 or strncmp(sst, "NaN",3) ==0 or strncmp(cur, "NaN",3) == 0)
    {
        printf("\n ERROR: \n Int or SST or Cur not defined!!! \n %s %s %s \n \n", interact, sst, cur);
        exit(1);
    }
    
    // UnComment this for the MPA loop
    //for (double MPAsize = 4; MPAsize < 25; MPAsize = MPAsize + 2)
    //for (int iA = 0; iA < sizeof(MPAs)/sizeof(*MPAs); ++iA)
    
    // UnComment this for the TACpp loop
    
    //for (int iA = 0; iA < sizeof(TACs)/sizeof(*TACs); ++iA)
            // old version: for (double TACpp = 0.03; TACpp < 0.15; TACpp = TACpp + 0.05)
    
    // UnComment this for the TAC msy
    // This value is irrelevant!!
    //double TACpp = 9999;
    
    {
    //#if TACloop
    //double TACpp = TACs[iA]; // this will only work when TACLOOP is true
    double TACpp = inputTAC;
    //double MPAsize = MPAs[iA]; //This will pass the size from the array to the real code variable
    //#endif 
    printf("Year of MPA: %i", YearOfMPA);
    printf("MPA: Left Bottom Corner [x,y]: %2.0f - %2.0f \n ", minMPAx, minMPAy);
    printf("MPA: Right Top Corner  [x,y]: %2.0f - %2.0f \n ", maxMPAx, maxMPAy);
    printf("%s_%s \n", _managX, _managY);
    printf("%s %s %s \n", interact, sst, cur);
    
            seedRandom();  

            //Let's create a World:
            World W;
            W.TACpp = 0;
            W.YearlyAccumulatedLandings = 0;
            W.howManyFishToRepToday = 0;
            if (tacFlag)
            {
                W.TACpp = TACpp/100;        
            }
            W.initSize(7,6, (GRIDSIZEX-1)/7.0, (GRIDSIZEY-1)/6.0 );	// 7 oceans by 
                                    // 6 oceans, size of each ocean NB "...-1"
            //W.initFishFromFile("2013initialposwithangles.dat");
            //W.initFishFromFile("FishInitialLong.dat");
            W.initFishFromFile(IniFishFile);
            W.randomFish(FISHPEROCEAN, DEFAULT_SPEED_LOWER_BOUND, 
                    DEFAULT_SPEED_UPPER_BOUND, DEFAULT_SELFWEIGHT);
 
            /*puts fish into the oceans defines the density, speed distribution and 
             * selfweight*/
            //W.alignedFish(FISHPEROCEAN, DEFAULT_SPEED_LOWER_BOUND, 
            //          DEFAULT_SPEED_UPPER_BOUND, 0.5, -0.866, DEFAULT_SELFWEIGHT);	
            /*puts fish into the oceans defines the density, speed distribution, 
             * direction and selfweight*/	

            W.isTorus = false;

            W.initPictureSerial();	//Allocates memory for the picture of the world
            W.allocateGrid(GRIDSIZEX, GRIDSIZEY);	//Allocates the grid
            W.zeroGrid( GRIDSIZEX, GRIDSIZEY ); //Puts zero into all values in the grid
            W.emptyLandOfFish(); //Remove all the fish that is over land.
            W.translateFtoGrid( GRIDSIZEX, GRIDSIZEY, F_Mortality ); //Create 0 and 1 to define the MPA, where 0 is the MPA

            #if KAI_CURRENTS
                    W.translateStraumurToGrid("KaiCurrents/KaiCurrent_1999_305.dat",
                            GRIDSIZEX, GRIDSIZEY);	//Same thing for the straumur
            #endif
            if (USECURRENTS)
            {
                W.translateStraumurToGrid("CurrentsFields.dat",GRIDSIZEX, GRIDSIZEY);
		printf("Currents files loaded!");
            } else
            {
                 W.translateStraumurToGrid("CurrentsNO.dat",GRIDSIZEX, GRIDSIZEY);    
            }
            W.findNearestGridPoint( GRIDSIZEX, GRIDSIZEY );	
                    //Every fish finds its nearest grid point
                    //The World calls a function (of the same name) of every ocean
                    //The Oceans then loop through the fish (in the ocean)
                    //The Fish then in turn call their member functions and find the
                    //      nearest grid point

            for (int i=0; i<7; i++)
                    for (int j=0; j<6; j++)
                            {
                            W.oceans[i][j].removeMarkedFish();
                            for (int k=0; k<W.oceans[i][j].numberOfFish; k++)
                                    {
                                    if (TYPEOFREP == 0)
                                    {
                                        W.oceans[i][j].fish[k].dayOfReproduction = probOfReproduction(TYPEOFREP);
                                    }
                                    else if (TYPEOFREP == 1)
                                    {
                                        W.oceans[i][j].fish[k].dayOfReproduction = probOfReproduction(TYPEOFREP, MU, STD);
                                    }
                                    //Here I assign the age of the fish at the start of the simulation
                                    W.oceans[i][j].fish[k].age = rand()%(STARTAGE) + STARTAGE/2;
                                    //printf("%i ", W.oceans[i][j].fish[k].age);
                                    }
                            }
            W.totalNumberOfFish();
            W.initializeFishID();

            #if DEB_ON
                    //W.initializeDEB(0.83,1.0,0.15,0.022,MATURE,0.0,0.0); 
                    //(L/Lm, [E]/[Em], uR, er, MaturityLevel, daysSinceMature, Time)
                    W.initializeDEB(0.83275,0.66647,0.279,0.12126,MATURE,9.0,0.0); 
                    //(L/Lm, [E]/[Em], uR, er, MaturityLevel, daysSinceMature, Time)
            #else
                    W.initializeMaturityLevels(MATURE);
            #endif

            printf("Total number of fish is now %d \n", W.howManyFish);

            int iters = ITERS;
            int numTimeStepsPerDay = (int) (1/TIMESTEP);
            printf("Time step %9.6f \n", TIMESTEP);
            printf("Number of time steps per day: %d \n", numTimeStepsPerDay);
            int numDays=0;

            for(int i=0; i<ITERS; i++)
            {   
                int year, day, dayYear;
                year = i / numTimeStepsPerDay / 365;
                day  = i / numTimeStepsPerDay;
                bool mpaFlag = false;

                if (UseMultiSST == 0)
                    {
                    if (i == 0) 
                            {
                            printf("No Multi SST -->> Ago \n");
                            W.initializePrefTempRanges(Summer_TC_IMMATURE, 
                                    Summer_TH_IMMATURE, Summer_TC_MATURE, 
                                    Summer_TH_MATURE);// Summer
                            if (USETEMPFIELD)
                            {
                                    W.translateTemperatureToGrid("SST/Box/y_2010_m_08.dat", GRIDSIZEX, GRIDSIZEY ); 
                                    //W.translateTemperatureToGrid("SST/Box/y_2010_m_08_17Box.dat", GRIDSIZEX, GRIDSIZEY ); 
                            } else 
                            {
                                    //W.translateTemperatureToGrid("SST/NoGradient.dat", GRIDSIZEX, GRIDSIZEY ); 
                                    W.translateTemperatureToGrid("SST/TempFieldRectangle_noGradient.dat",GRIDSIZEX, GRIDSIZEY );  
                            }
                            //W.translateTemperatureToGrid("SST/TempFieldRectangle_noGradient.dat", GRIDSIZEX, GRIDSIZEY ); 
                            }
                    }
                if (UseMultiSST == 1)
                    {
                    int iter = i - (year * 365 * numTimeStepsPerDay);
                    if (i == 0 or iter == dJan ) 
                         {
                         printf("Jan \n");
                         W.initializePrefTempRanges(Summer_TC_IMMATURE, 
                                 Summer_TH_IMMATURE, Summer_TC_MATURE, 
                                 Summer_TH_MATURE);// Summer
                         W.translateTemperatureToGrid("SST/Box/y_2009_m_01.dat", 
                                 GRIDSIZEX, GRIDSIZEY );  
                         }
                    if (iter == dFeb) 
                        {
                        printf("Feb \n");
                        W.initializePrefTempRanges(Summer_TC_IMMATURE, 
                                Summer_TH_IMMATURE, Summer_TC_MATURE, 
                                Summer_TH_MATURE);// Summer
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_02.dat", 
                                GRIDSIZEX, GRIDSIZEY ); 
                        }
                    //Fall
                    if (iter == dMar) 
                        {
                        printf("Mar \n");
                        W.initializePrefTempRanges(Fall_TC_IMMATURE, 
                                Fall_TH_IMMATURE, Fall_TC_MATURE, Fall_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_03.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    if (iter == dApr) 
                        {
                        printf("Apr \n");
                        W.initializePrefTempRanges(Fall_TC_IMMATURE, 
                                Fall_TH_IMMATURE, Fall_TC_MATURE, Fall_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_04.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    if (iter == dMay) 
                        {
                        printf("May \n");
                        W.initializePrefTempRanges(Fall_TC_IMMATURE, 
                                Fall_TH_IMMATURE, Fall_TC_MATURE, Fall_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_05.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    // Winter
                    if (iter == dJun) 
                        {
                        printf("Jun \n");
                        W.initializePrefTempRanges(Winter_TC_IMMATURE, 
                            Winter_TH_IMMATURE, Winter_TC_MATURE, Winter_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_06.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    if (iter == dJul) 
                        {
                        printf("Jul \n");
                        W.initializePrefTempRanges(Winter_TC_IMMATURE, 
                            Winter_TH_IMMATURE, Winter_TC_MATURE, Winter_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_07.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    if (iter == dAgo) 
                        {
                        printf("Ago \n");
                        W.initializePrefTempRanges(Winter_TC_IMMATURE, 
                            Winter_TH_IMMATURE, Winter_TC_MATURE, Winter_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_08.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    // Spring
                    if (iter == dSep) 
                        {
                        printf("Sep \n");
                        W.initializePrefTempRanges(Spring_TC_IMMATURE, 
                            Spring_TH_IMMATURE, Spring_TC_MATURE, Spring_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_09.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        } 
                    if (iter == dOct) 
                        {
                        printf("Oct \n");
                        W.initializePrefTempRanges(Spring_TC_IMMATURE, 
                            Spring_TH_IMMATURE, Spring_TC_MATURE, Spring_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_10.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    if (iter == dNov) 
                        {
                        printf("Nov \n");
                        W.initializePrefTempRanges(Spring_TC_IMMATURE, 
                            Spring_TH_IMMATURE, Spring_TC_MATURE, Spring_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_11.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    // Summer       
                    if (iter == dDec) 
                        {
                        printf("Dec \n");
                        W.initializePrefTempRanges(Summer_TC_IMMATURE, 
                            Summer_TH_IMMATURE, Summer_TC_MATURE, Summer_TH_MATURE);
                        W.translateTemperatureToGrid("SST/Box/y_2009_m_12.dat", 
                                GRIDSIZEX, GRIDSIZEY );  
                        }
                    }

                    //W.setReproduction(i, true); // Here I call the function to change the reproduction status of the fish

                double toKillToday;
                if (i == 1)
                {
                    W.emptyLandOfFish();
                    W.nOfAdults = W.GetDensities(GRIDSIZEX, GRIDSIZEY);
                    //toKillToday = W.adjustLarvaeMortality(0, W.nOfAdults, 300000); //W.howManyFish - W.nOfAdults);
                    toKillToday = 0.01;
                     // To define potential value to start...
                    printf("To kill today %2.8f\n", toKillToday);
                    printf("Accumulated Landings %10.0f \n", W.howManyFishToRepToday);
                }

           // Defining the different temperature preferences for each season 
       
            if (i%100 == 0)
                          {
                                W.totalNumberOfFish();
                                W.nOfAdults = W.GetDensities(GRIDSIZEX, GRIDSIZEY);
                                if (mpaFlag)
                                {
                                    printf("[%2.0f, %i, %i] #Fish: %d, #Adults: %d, #Rep: %6.0f, AccLand: %i, TAC is: %2.0f \n", MPAsize, day, year, W.howManyFish, W.nOfAdults, W.howManyFishToRepToday, W.YearlyAccumulatedLandings, W.TAC);
                                } else if (tacFlag)
                                {
                                    printf("[%2.3f, %i, %i] #Fish: %d, #Adults: %d, #Rep: %d, AccLand: %i, TAC is: %2.0f \n", W.TACpp, day, year, W.howManyFish, W.nOfAdults, W.howManyFishToRepToday, W.YearlyAccumulatedLandings, W.TAC);
                                }
                            }

             if ( i%numTimeStepsPerDay==0 ) //We have a new day                                  
                    {
                        //printf("%i \n", W.howManyFishToRepToday);
                        W.totalToRepToday(i*TIMESTEP);
                       if (i/numTimeStepsPerDay % (365/1) == 0  or i/numTimeStepsPerDay == 365)
                        {
                            if (W.howManyFish - W.nOfAdults > 0)
                            {
                                toKillToday = W.adjustLarvaeMortality(0, W.nOfAdults, W.howManyFish - W.nOfAdults);
                            }
                            else (toKillToday = 0.01);
                            //printf("%i, %i, %2.8f\n", W.howManyFish, W.nOfAdults, toKillToday);
                        }
                        W.FishAging(W.howManyFish, toKillToday, W.howManyFish - W.nOfAdults);

                        double mpaCenter;
                        
                        //This make MPA, it is defined only once!!
                        if (i%(365*20*(YearOfMPA-1))==0 and mpaFlag == false and lengthX == 0 and lengthY == 0) //This means that the MPA will be set at year before it is implemented!!
                                {
                                    mpaCenter = W.mpaDef(MOVMPA, MPAsize); // This calls the function to get the new center for the MPA!!
                                    minMPAy = mpaCenter - sMPA;
                                    maxMPAy = mpaCenter + sMPA;
                                    printf("%.f, %.f, %.f <<-- \n", minMPAy, mpaCenter, maxMPAy);
                                    mpaFlag = true;
                                }
                        
                        if (i >= 100) 
                        {
                            #if USEREPRODUCTION 
                                    //printf("%i  ",W.oceans[4][4].numberOfFish);
                                    W.LastUsedIDFish = W.addFish(W.howManyFish, i*TIMESTEP, nEggs);
                                    //printf("%i  \n",W.oceans[4][4].numberOfFish);
                            #endif

                            W.GetDensities(GRIDSIZEX, GRIDSIZEY);
                            W.BoatsDist(i, GRIDSIZEX, GRIDSIZEY, minMPAx, minMPAy, maxMPAx, maxMPAy, YearOfMPA, boatsC);

                            #if WRITEGRID // This is to write the grid with the detailed data!!		    	
                                if (i%(365*20) == 0)
                                {   
                                    char file[255];
                                    file[0] = 0;
                                    
                                    if (lengthX == 0 and lengthY == 0){
                                        if (mpaFlag && !tacFlag)
                                        {
                                                sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_%.f_%s_GridData_%i.txt", manag, interact, sst, cur, RunName,  MPAsize, _F, 100000000 + year * 100000 + i/numTimeStepsPerDay);
                                                printf("mpaFlag \n");

                                        } else if (!mpaFlag && tacFlag)
                                        {
                                                sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_%2.2f_%s_GridData_%i.txt", manag, interact, sst, cur, RunName, W.TACpp, _F, 100000000 + year * 100000 + i/numTimeStepsPerDay);
                                                printf("tacFlag \n");
						
                                        } else if (mpaFlag && tacFlag)
					{
                                                sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_%2.2f_%s_%s_GridData_%i.txt", manag, interact, sst, cur, RunName, MPAsize, W.TACpp, _F, 100000000 + year * 100000 + i/numTimeStepsPerDay);
						printf("Print to file %s \n", file);
					}
                                        else if (lengthX == 0 and lengthY == 0){
                                                sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_%.f_%s_GridData_%i.txt", manag, interact, sst, cur, RunName, MPAsize, _F,  100000000 + year * 100000 + i/numTimeStepsPerDay);
                                                printf("lengthX and lengthY == 0 \n");
                                        }
                                    }  else if (lengthX != 0 and lengthY !=0) {
	                                //sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_x_%s_y_%s_%s_GridData_%i.txt", manag, interact, sst, cur, RunName, _managX, _managY, _F,  100000000 + year * 100000 + i/numTimeStepsPerDay);
                                        sprintf(file, "./GridOutPut_%s_%s_%s_%s/%s_x_%s_y_%s_C_%2.1f_e_%i_tac_%2.2f_%s_GridData_%i.txt", manag, interact, sst, cur, RunName, _managX, _managY, boatsC, nEggs, W.TACpp, _F, 100000000 + year * 100000 + i/numTimeStepsPerDay);
                                        printf("!=0 \n");
                                    } 
                                    
                                    
                                    W.WriteXYDensity(file, GRIDSIZEX, GRIDSIZEY, sst, cur, interact, manag, _managX, _managY);
                                }
                            #endif

                            if (year >= YearOfFishing)
                            {
                                if (i%(365*20)==0)  
                                {
                                    W.YearlyAccumulatedLandings = 0;
                                    printf("%d\n", W.YearlyAccumulatedLandings);
                                }
                                W.TAC = W.checkTAC(i, W.howManyFish, tacFlag);
                                if (year < YearOfTAC or tacFlag == false)
                                {
                                    W.FishingMortality(i, GRIDSIZEX, GRIDSIZEY);
                                }
                                else if (tacFlag == true and W.TAC >= W.YearlyAccumulatedLandings)
                                {
                                    W.FishingMortality(i, GRIDSIZEX, GRIDSIZEY);
                                }
                            }

                            

                        }

                    //printf("Last used ID in a fish is %i (%i, %i)\n", W.oceans[1][2].fish[0].age, 1, 2);

                            numDays++;
                            char tFilename[255];
                            tFilename[0]=0;

                            #if KAI_CURRENTS
                                    char CFilename[255];
                                    CFilename[0] = 0;
                            #endif

                            int numDays = (int) i/numTimeStepsPerDay;

                            #if KAI_TEMPERATURE
                            printf("i: %d numDays %d ", i, numDays);	
                            if (i< 1220) //Still in 1999
                               {
                                    sprintf( tFilename, 
                                            "KaiTemperature/KaiTemp_1999_%d.dat", 305+numDays );
                                    #if KAI_CURRENTS
                                            sprintf( CFilename, "KaiCurrents/KaiCurrent_1999_%d.dat", 305+numDays );
                                    #endif
                            }
                            else //Made it to 2000
                            {

                                    sprintf( tFilename, "KaiTemperature/KaiTemp_2000_%d.dat", numDays-60 );
                                    #if KAI_CURRENTS
                                            sprintf( CFilename, "KaiCurrents/KaiCurrent_2000_%d.dat", numDays-60 );
                                    #endif
                            }
                            #endif	

                            //W.translateTemperatureToGrid("KaiTemp_1999_nov.dat", GRIDSIZEX, GRIDSIZEY );	//Reads temperature info into grid
                            //W.translateTemperatureToGrid(tFilename, GRIDSIZEX, GRIDSIZEY );	//Reads temperature info into grid
                            #if KAI_CURRENTS
                                    W.translateStraumurToGrid(CFilename,GRIDSIZEX, GRIDSIZEY);	//Same thing for the straumur
                            #endif

                    } 

                    W.resetInteractionCounters();	/* ... to zero */

                    W.communicateGhostFish();
                    if (USE_Interactions)
                    {
                            W.interact();	/* just calculate the interaction vector */
                    }
                    
                    W.removeGhostFish();

                    W.factorInTemperatureGradients(NORM_TEMPERATURE_SENSITIVITY, GRIDSIZEX, GRIDSIZEY);	/* calculate temperature gradient and
                                                                                                                    adjust directional heading accordingly */
                    #if ADDNOISE
                            W.addNoiseToDirectionAngle(AMPLITUDEOFNOISE);
                    #endif

                    #if DEB_ON
                            W.solveDEB(TIMESTEP/DEB_RK_TIMESTEPS_PER_TIMESTEP, i*TIMESTEP, (i+1)*TIMESTEP);
                    #endif

                    W.move();
                    W.translateByStraumur(TIMESTEP, GRIDSIZEX, GRIDSIZEY);
                    W.transportFish(); /* put the fish into the correct ocean */           
                    
		    if ((year >= 20 & year < 31))
		    {                    
                    if( i%(CHOOSEEVERYNTHFRAME*(365/365)) == 0 & tga == true) //365
                            {
                            char s[255];
                            s[0] = 0;
                            char sa[255];
                            sa[0] = 0;
                            int a = 1;
			    int mmpa = MPAsize;
                            if (isMPA == true and tacFlag == false)
                            {
                                sprintf( s, "./fig/%s_%i_out%d.tga", RunName, mmpa, i/numTimeStepsPerDay+100 );
                                
                            }
                            if (isMPA == false and tacFlag == true)
                                {
                               int percent = W.TACpp * 100;
                               sprintf( s, "./fig/%s_%i_out%d.tga", RunName, percent, i/numTimeStepsPerDay+100 );
                                }
			    if (isMPA == true and tacFlag == true)
				{
				int percent = W.TACpp * 100;
                                sprintf( s, "./fig/%s_%i_out%d.tga", RunName, percent, i/numTimeStepsPerDay+100 );
				}
                            if (lengthX != 0 and lengthY !=0) 
                            {
                               sprintf(s, "./fig/%s_%2.0f_out%d.tga", RunName, MPAsize, i/numTimeStepsPerDay+100);
                               sprintf(sa, "./fig/Adults_%s_%2.0f_out%d.tga", RunName, MPAsize, i/numTimeStepsPerDay+100);
                            }                                          
                            W.drawPicture(s, GRIDSIZEX, GRIDSIZEY);
                            W.drawPictureAdults(sa, GRIDSIZEX, GRIDSIZEY);
                            }
                    }                      
                    if (W.howManyFish == 0) {break;}
                    
                    if (resTime == 1 and i%(numTimeStepsPerDay*365/365)==0 and year >= 15 and year <= 19)
                    {
                        char s[255];
                        s[0]=0;
			char temp[255];
			temp[0]=0;
			float size;
			size = (maxMPAx-minMPAx)*(maxMPAy-minMPAy)/40;
			sprintf(temp, "%.f", size);
                        sprintf(s, "./resTime/%s_mpaSize-%s_%s_%s_%s_%i.csv", RunName, temp, interact, sst, cur,  1000000+day); 
                        W.WriteGridFishIDs(s, day, GRIDSIZEX, GRIDSIZEY, minMPAx, maxMPAx, minMPAy, maxMPAy);
                        if (year == 2)
                        {
                            exit(1);
                        }
                    }
                    if (i%(numTimeStepsPerDay) == 0 & 1 == 0 & year >=0) //This is just to stop this execution;
                    {
			if (year>30){break;}
			printf("#Year: %i, #Day %i  ", year, day);
                        W.writeRparam(year, day, W.howManyFish, W.nOfAdults, W.Rparam(W.howManyFish, W.nOfAdults));
                    }                    
            }
            printf("Total number of fish is now %d \n", W.howManyFish);
            W.deallocateGrid();
    }
    return 0;
}
