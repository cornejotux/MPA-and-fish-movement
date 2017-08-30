#pragma once

using namespace std;

/* this is the header file to make the .tga images work */
#include "tga.h"

#if MPI_ON  == 1
        #include "mpi.h"
#endif

#define ITERS (365/TIMESTEP)*70 
#define TIMESTEP 0.05
#define ITERPERDAY 20
#define CHOOSEEVERYNTHFRAME 20  //1 image each # of iterations, this will depends on how many 
                                //time steps we 

/*	These determine how much weight the particles place on the
direction determined by the interaction vs. the temperature gradients */
#define ALPHA 0.995  // This only works when DEB_ON 
//alpha is used to determine the acceleration
#define BETA 0.005
#define SSTmin 16
#define SSTmax 20
#define DEFAULT_INTERACTIONWEIGHT	(1-BETA)
#define DEFAULT_TEMPERATUREWEIGHT	BETA

//From what file the fish is read!
//#define IniFishFile "FishInitialLong.dat"
#define IniFishFile "FishInitialMedium.dat"

#define MAINDIR "PICS_2010_08"

#define KAI_CURRENTS 0
#define KAI_TEMPERATURE 0

/*	This defines the amount of weight a particle places on its own directional heading
	NB: updateVelocity might produce NaN if selfweight == 0 since it divides by it (plus possibly more) */
#define DEFAULT_SELFWEIGHT 1.0

// These define the radii of the zones of interaction
#define DEFAULT_RADIUS_OF_REPULSION     0.02
#define DEFAULT_RADIUS_OF_ORIENTATION	0.10
#define DEFAULT_RADIUS_OF_ATTRACTION	0.10
	/* NOTE: Radius of Attraction must be at least as large as radius of orientation
	or the code will not work the way you expect (the code sorts fish so that they ignore
	all fish not inside their radius of attraction!) */

/*	These are pretty self-explanatory */
#define DEFAULT_SPEED_LOWER_BOUND  0.42 
	//about 12 km
#define DEFAULT_SPEED_NO_ROE       5 
    //NB In km/day!
#define MAX_SPEED_INCREASE         10
    // NB In km/day
#define DEFAULT_SPEED_UPPER_BOUND  .5
#define EGGSPEED        0.002  //Speed that the eggs will move until AGEOFSELFMOVIIMENT

#define AMPLITUDEOFNOISE 0.0

//Define some population dynamics parameters
#define TOTALCC          30000 // This is the total number of adults that the system can sustain.
#define OCEANCC          TOTALCC/(6*7)
#define LARVAECC         200000
#define OCEANLARVAECC    LARVAECC/(6*7)
#define NOFEGGSPERFEMALE 10
///// Parameters to set the day of reproduction
#define TYPEOFREP        1 // Type of rep = 0 means uniform distribution trough the year 
                           // 1 mean normal distribution with mu and std
#define MU               250
#define STD              20 
///////////////////////////////

#define USEREPRODUCTION  true 
#define USEMORTALITY     true
#define USEFISHING       true
#define WRITEGRID        true //This is to write the Densities and catches at each x,y

#define MOVMPA           0     // 0 is not moving and 1 define a moving block.
                               // more option will be added...


#define NORM_TEMPERATURE_SENSITIVITY 0.00001
#define DECISION_TOL pow(10,-5)
#define SMALL_NUMBER .5
#define DEFAULT_PADDING (DEFAULT_RADIUS_OF_ATTRACTION + SMALL_NUMBER)

/* These give us the preferred temperature range, between TOOCOLD and TOOHOT for each maturity level
   See World::initializePrefTempRanges */

///// Here I define the if temp preferences are seasonal or the same for the whole year
#define SeasonPref 0

#if SeasonPref == 0
        #define Summer_TC_IMMATURE SSTmin //17.5
        #define Summer_TH_IMMATURE SSTmax // 23
        #define Summer_TC_MATURE   SSTmin //17.5
        #define Summer_TH_MATURE   SSTmax //23

        #define Fall_TC_IMMATURE   Summer_TC_IMMATURE
        #define Fall_TH_IMMATURE   Summer_TH_IMMATURE
        #define Fall_TC_MATURE     Summer_TC_MATURE
        #define Fall_TH_MATURE     Summer_TH_MATURE 

        #define Winter_TC_IMMATURE Summer_TC_IMMATURE
        #define Winter_TH_IMMATURE Summer_TH_IMMATURE
        #define Winter_TC_MATURE   Summer_TC_MATURE
        #define Winter_TH_MATURE   Summer_TH_MATURE

        #define Spring_TC_IMMATURE Summer_TC_IMMATURE
        #define Spring_TH_IMMATURE Summer_TC_IMMATURE
        #define Spring_TC_MATURE   Summer_TC_MATURE
        #define Spring_TH_MATURE   Summer_TH_MATURE
#else
        #define Summer_TC_IMMATURE 17.5
        #define Summer_TH_IMMATURE 23.0
        #define Summer_TC_MATURE   17.5
        #define Summer_TH_MATURE   23.0

        #define Fall_TC_IMMATURE   16.5
        #define Fall_TH_IMMATURE   20.0
        #define Fall_TC_MATURE     16.5
        #define Fall_TH_MATURE     20.0 

        #define Winter_TC_IMMATURE 16.5
        #define Winter_TH_IMMATURE 17.5
        #define Winter_TC_MATURE   16.5
        #define Winter_TH_MATURE   17.5

        #define Spring_TC_IMMATURE 16.0
        #define Spring_TH_IMMATURE 19.0
        #define Spring_TC_MATURE   16.0
        #define Spring_TH_MATURE   19.0
#endif



// Here are the parameters defined for the fish!
/* This number has to correspond to the number of maturity levels, enumerated below */
#define NUMMATURITYLEVELS 2
#define AGEMATURE    365 //This is the age (in days!) when a fish mature and can reproduce
#define AGEOLD       365 * 4
#define AGESELFMOVE  365*2/3 //Define when the fish start swimming.
#define STARTAGE     365*2  // This is the age at what the fish from the file will start. Asummes all the fish are are adults of 1 year
#define Z_JUVENILES  0.002  
#define Z_ADULTS     0.00035 
#define Z_OLD        0.0015 

#define totOfBoats   5000
#define r0           1.2    // Intrinsic rate of increase

#define BaseTAC      0 // This is the base number that the fishery can take any year, it is not affected by the number of fish available

///////////////////////////////////
// Define the MPA implementation //

#define YearOfFishing  20 // Year then the fishing start

///////////////////////////////////


/*	This is how many fish are initially uniformly distributed inside each ocean */
#define FISHPEROCEAN 10 

#define DOUBLE_DETERMINING_RED_DENSITY 75

#define DEB_ON  0 // originally 1

/*	Allows us to print every nth fish out to a file; prints a fish's info out if fish.ID is 0 mod thisNumber */
#define CHOOSE_EVERY_NTH_FISH_TO_PRINT_DEB 2000
#define DEB_RK_TIMESTEPS_PER_TIMESTEP 1 
//	For Runge-Kutta method in solveDEB

/* A bunch of DEB constants which need to be determined in some reasonable way */
#define DEB_GAMMA 0.20
#define DEB_NU    0.02
#define DEB_KJ    0.001
#define DEB_KAPPA 0.4
#define DEB_J_EAm 0.23
    //mmol d-1 cm-20
#define DEB_mu_E  500.0
    //J mmol-3 chemical potential
#define DEB_M_V   4.4
    //mmol cm-3
#define DEB_y_VE  0.8
    //- yield of structure

#define DEB_TA 9100.0 //Done by BE
	/* (K) Found by BE */
#define DEB_T1 279.65
	/* (K) i.e. 6.5 degrees */
#define DEB_LM 0.6086956521739 //delta*max(physical)Length 0.161*17.5
		/* NB delta times the maximum physical length, gives the max volumetric length */

#define DEB_EM (DEB_mu_E*DEB_J_EAm/DEB_NU)
    //5860 is the maximum from Anthony et al. 2000
#define DEB_E_G (DEB_mu_E*DEB_M_V/DEB_y_VE)
    //2800 is the value from the anchovies
#define DEB_U2E (DEB_EM*DEB_LM*DEB_LM*DEB_LM)
#define DEB_DV 1 //correct
#define DEB_RHO_E 39300
	// J g^{-1} NB see Anthony et al (2000)
#define DEB_RHO_ROE 10000
// J g^{-1} (This is where we take into account the increasing water content of roe) Need from Matis perhaps
#define DEB_G (DEB_E_G/(DEB_KAPPA*DEB_EM))
#define DEB_UHP (3930.4/DEB_U2E)

#define DEB_ROE_PERCENTAGE_MARK 8.0
    /* Used to determine "Maturity", i.e. when water content of roe starts increasing,
     and also when speed starts speeodelding up. */
#define DEB_ROE_PERCENTAGE_MAX 25.0
    /* After this has been reached, the */
#define DEB_TIME2MAX 30.0
    // Days - water content increases by 20% in that time
#define DEB_ROE_MAX_WATERINCREASE 0.20


/*	These determine how the picture is drawn, i.e. if TOOCOLD and TOOHOT are drawn,
	if the land is solid, if the outline of the land is there */
#define DRAWTEMP false		
#define DRAW_LAND_OUTLINE 0
#define DRAW_LAND_SOLIDLY 1

#define USEGRID true

#define ADDNOISE false

/*	These tell us how large the grid that we need to define is */
#define GRIDSIZEX 50 
#define GRIDSIZEY 100

/*	These all affect the size of the final .tga output */
#define PICTURESCALING 10
#define PICTURESIZEX   100
#define PICTURESIZEY   100

/* These four bounds help to initialize the fish correctly and are passed into initializeFishPositions_Rectangle() */
#define FISHPOSITIONLOWERX 0 
#define FISHPOSITIONUPPERX 50
#define FISHPOSITIONLOWERY 0 
#define FISHPOSITIONUPPERY 100

/* This switches the parallel algorithm on and off.  On is 1, off is 0 */
#define MPI_ON 0  // <--- Here is defined to run in parallel...

#define TRIANGULAR_GRADIENT 0

// initial iteration corresponding for each month
// specially important if seasonal SST and/or currents are used.

#define dJan 1
#define dFeb (1/TIMESTEP)*31+1
#define dMar (1/TIMESTEP)*(31+28)+1
#define dApr (1/TIMESTEP)*(31+30+28)+1
#define dMay (1/TIMESTEP)*(31*2+30+28)+1
#define dJun (1/TIMESTEP)*(31*3+30+28)+1
#define dJul (1/TIMESTEP)*(31*3+30*2+28)+1
#define dAgo (1/TIMESTEP)*(31*4+30*2+28)+1
#define dSep (1/TIMESTEP)*(31*5+30*2+28)+1
#define dOct (1/TIMESTEP)*(31*5+30*3+28)+1
#define dNov (1/TIMESTEP)*(31*6+30*3+28)+1
#define dDec (1/TIMESTEP)*(31*6+30*4+28)+1



/* Fish live in the Ocean class */
class Fish;
/* GridPoint stores info on temp and current for a point */
class GridPoint;
/* Ocean is a rectangular piece of the "World", so that the fish don't need to interact with all fish in the world */
class Ocean;
/* World is an array of Oceans */
class World;

/* Processor has to do with the parallel algorithm */
class Processor;


/* Defines tags MIGRATING and GHOST to tell which fish are migrating and which fish are ghost fish
 (for when we pass fish between processors) */
enum {
	MIGRATING =1,
	GHOST =2
};


/* We now define integer constants which correspond to the maturity levels of the fish.
   The numbers should be used as indices in arrays - i.e. 0,1,2,...
   See also #define numMaturityLevels above */
enum {
	IMMATURE=0,
	MATURE=1
};


/*	A FishRecord is a struct that will be used for passing messages about fish.  
The goal is to pass less total data in the message to
improve performance. */

typedef struct FishRecord {
	double x;
	double y;
	double cosPhi;
	double sinPhi;
	float speed;
	double selfWeight;
} FishRecord;

/*	This defines a vector in R^2 for use for the straumur data */
typedef struct Vector {
        float x, y;
        } Vector;
       
typedef struct rVector {
        float all, adults;
        } rVector;

/*	This defines a vector in Z^2 for use for the gridpoint data */
typedef struct IntegerVector {
int x, y;
} IntegerVector;




/* */
typedef struct ExchangeData {
	int noExchange;
	int x, y;
	int targetThread;
} ExchangeData;



/*	These two functions are not member functions of any class  */
void seedRandom();

double randomVariable();

int probOfReproduction(int type) ; // This function will assign the probability of 
int probOfReproduction(int type, double mu, double std);//reproduction to happen in a any given day of the year!
                            //type is what kind of probability function will be used




class Fish {
public:
	Fish();

/*	FISH CLASS MEMBER VARIABLES */

	Ocean* myOcean; // This is a pointer to the ocean where the fish currently resides

	bool removeMe;
	bool isGhost;
        bool isAdult;     //If they are adult, they can follow the temperature an can reproduce!
        bool repThisYear; //If the fish reproduced this year
        float pp;         //Probability to to be used in the mortality functions 

        
	/*	Each fish has an ID and a current position x,y
		and a velocity made up of a unit-length vector: (cosPhi,sinPhi)
		and a speed (the magnitude of the velocity) */
	int ID;
        int age; //This is used to track the age: Identify when can reproduce, recruit to the fishery
                // follow the temperature gradient, and increase the natural mortality rate as they get old!
        int dayOfReproduction; //This is the day when the fish if going to reproduce
	double x;
	double y;
	//Vector oldPosition;

	/*	This is all for using the grid of temperature and straumur */
	IntegerVector nearestGridPoint;	// this is to find the info to use from the grid of straumur and temperature
	int quadrant;	// options here are 0,1,2,3, where they are labelled like the quadrants except minus 1
					// we use this to keep track of which temperature gradient to use in World member functions getTemperatureGradient()
					// and getTemperatureGradient2()

	//Used as the escape region:
	int region;

	/*	These define the weight that each fish places on the interaction and temperature
			Could change with maturity, perhaps to allow different reaction during feeding */
	double interactionWeight;
	double temperatureWeight;


	/*	These member variables control the velocity of the fish;
		phi is the direction angle, speed is the magnitude of the velocity,
		selfWeight is the weight that the fish gives to its own directional heading */
	double cosPhi;
	double sinPhi;
	double speed;
	double selfWeight;
	/*	When the velocity updates, we need temporary storage for the old velocity information so that
		as fish update their current velocities, their old velocities stick around so other fish can update using them.*/
	double oldCos;
	double oldSin;
	double oldSpeed;

	/* DEB variables */
	double l;		// scaled volumetric length [-]
	double e;		// scaled energy reserves (fat) [-]
	double uR;	// scaled reproduction energy buffer [-]
	double er;	// scaled energy contained in roe
	
	int maturityLevel;
	double daysSinceMature;//
	double roePercentage;

/*	FISH CLASS MEMBER FUCTIONS */

	void print() const;
	void printToFile() const;
        //This function is to change the reproduction status of the fish
        void setReproduction(int i, bool Reproduced);
   
	/*	this function finds the nearest gridpoint on the grid that stores straumur and temperature */
	void findNearestGridPoint(int gridSizeX, int gridSizeY );

	/*	this function finds the quadrant of the nearest grid point where the fish is located */
	void findQuadrant();

	/*	Creating the neighbor and temperature weights */
	void updateWeights(double nWeight, double tWeight, double sWeight);

	/*	These are the radii of interaction which determine how fish react to neighbors*/
	double radiusOfRepulsion, radiusOfOrientation, radiusOfAttraction;

	/*	makeRadii sets the radii of interaction to default values*/
	void makeRadii();

	/*	zero out position and velocity*/
	void zero();

	/*	these functions allow us to randomly initialize fish*/
	void randomPosition(double displacementX, double displacementY, double amplitudeX, double amplitudeY);
	void randomDirection();
	void randomSpeed(double speedLowerBound, double speedUpperBound);

	/*	these are some accessory functions for the velocity*/
	void setPosition(double inx, double iny);
	void setDirection(double angle);
	void setVelocity(double inSpeed, double inCosPhi, double inSinPhi, double inSelfWeight);
	void resetDirection(double theta, double noiseAmp);
	void setVelocity(double vx, double vy, double inSelfWeight);
	void setSpeed(double inSpeed);
	void resetSpeed(double inSpeed);

	/*	copyVelocity sets the "old" velocity variables equal to the current ones */
	void copyVelocity();

	/*	distanceToSquared simply computes the square of the distance between 'this' fish and the fish F */
	double distanceToSquared(const Fish& F);

	
	/*	updateVelocity is where the interaction happens between this fish and its neighbors
		the list of other fish in the ocean is iterated and fish in the zone of attraction are considered */
	void updateVelocity();

	/*	addNoiseToDirectionAngle adds noise of the specified amplitude to the direction angle of the fish */
	void addNoiseToDirectionAngle(double noiseAmplitude);


	/*	This function takes the temperature and weighs it into the direction obtained
		previously from the neighbors--then normalized */
	//void weighInTemperature(int x, int y, int quadrant);

	/*	move according to current velocity*/
	void move(double timestep);



	/*	emptyLandOfFish() removes fish from the land (if the nearest gridPoint has temperature above 998) */
	void emptyLandOfFish();
	/*	this removes fish if it's outside the specified bounds */
	void initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary);
	/*	initializes fish inside the 8 ellipsoids specified below */
	void initializeFishPositions(Vector center1, double majax1, double minAx1, Vector center2, double majax2, double minAx2,
							Vector center3, double majax3, double minAx3, Vector center4, double majax4, double minAx4,
							Vector center5, double majax5, double minAx5, Vector center6, double majax6, double minAx6,
							Vector center7, double majax7, double minAx7, Vector center8, double majax8, double minAx8);


	/*	translate according to the current at the nearest node */
//	void translateByStraumur(double timestep);  This function was moved to the World class.

//private:
	/*	escapeRegion returns myOcean->escapeRegion(*this) so look at that method's comment to see what this does*/
	int escapeRegion();

	/*	dumps relevant data about a fish into the fishrecord data structure*/
	void makeRecord(FishRecord& R) const;

	/*	does the inverse of makeRecord*/
	void mimmicRecord(const FishRecord& R);

	void initializeMaturityLevels(int initialMaturityLevel);

	/*	DEB functions */
	void initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime);	// This function initializes the variables for the DEB
	void solveDEB(double h, double startTime, double endTime);		// This function solves the coupled DEB ODEs using 4th order Runge Kutta.
									// It saves a copy of the current length, fat, roe so they won't be overwritten while needed
									// At the end, it writes the new values into F.length, etc.
	double dl(double t, double length, double energy, double nu);	// Defines dl/dt for the volumetric length
	double de(double t, double length, double energy, double nu);	// Defines de/dt for the energy reserves
	double duR(double t, double length, double energy, double nu, double kJ);	// Defines duR/dt for the reproduction energy reserves
	double der(double t, double repBuffer, double roeEn, double gamma);			// Defines der/dt for the energy of the roe content
	
	/* Depends on roe maturity of the fish, see function for details*/
	double getPrefSpeed();

	void printInitialBiologyToFile();		// Begins the .m file with a variable name, etc.
	void printBiologyHistoryToFile(double time);		/*	Prints the length, energy reserves, reproduction buffer, roe energy, time 
												into an array of numbers which matlab can read */              
}; //class Fish



class GridPoint {
public:
	GridPoint();
	~GridPoint();

/*	GRIDPOINT CLASS MEMBER VARIABLES */
	double temperature;
	bool   tempCorrectionFlag;      //Used for DEB
	double temperatureCorrection; //Used for DEB calculations, default 1 in case DEB not in use
	double tempCorrected_nu;      //Used for DEB calculations, NB calculated in getTemperatureGradient when flags have not been set
	double tempCorrected_kJ;      //Used for DEB calculations, NB calculated in getTemperatureGradient when flags have not been set
	double tempCorrected_gamma;   //Used for DEB calculations, NB calculated in getTemperatureGradient when flags have not been set
	bool*  flagQuadrant; /* This tells you if you've computed the temperature gradient for a particular quadrant using method with triangles*/
	bool   gradientFlag; /* This tells you if you've computed temperatureGradientAtCenter */
	Vector*  temperatureGradientAtCenter; /* This stores the temperature gradient at a gridPoint made using Sven's method */
	Vector** temperatureGradient; /* This stores temperature gradients for the four quadrants and all maturity levels */
	double*  temperatureGradAtCenterFactor;
	double** temperatureGradFactor;

	Vector straumur;  /* This means "current" in Icelandic! */
        
        int    density; //Keep the density of all the fish 
        int    AdultDensity; //This is the density of fish > of AGEMATURE!
        double FM; //:These are the values of fishing mortality for each XY point
        int    FCatches; //Track of the accumulated amount of fish caught in each YX pixel
        int    AnnualCatches; //Annual catches in each grid point;
        float  FMpp;       // Number of potential fish to be captured based on the fishing mortality and density of adults.
        float  bp;         // B'i in Hilborn model
        float  boats;      // Number of boats to allocate in the each pixel of the grid.
        float  BoatsXFMpp; // Number of potential fish to be captured based on the fishing mortality and number of boats in each pixel.
        double    age0;       // Age class of 0 to 1 years
        double    age1;       // Age class of  1 to 2
        double    age2;       // Age class of  2 to 3
        double    age3;       // Age class of  4 to 4
        double    age4;       // Age class of  4 to 5
        double    age5;       // Age class of  5 and up
        
        
/*	GRIDPOINT CLASS MEMBER FUNCTIONS */
	void zero();
};




class Ocean {
public:
	Ocean(int inNumberOfFish);
	Ocean(const Ocean& O);
	Ocean();
	~Ocean();

/*	OCEAN CLASS MEMBER VARIABLES */

	/*	the difference between numberOfFish and sizeOfFish is: sizeOfFish is the size of the allocated block
		in memory (in fish, not in bytes), and sizeOfFish is the */
	int sizeOfFish;
	int numberOfFish;
        int numberOfFishToRep; //Store the # of fish that can reproduce
	int interactionCounter;
        int nJuveniles;
        int nAdults   ;
        int nOld      ;

	int RunOnThread;

	World* myWorld;
	Fish* fish;

private:
	/* tempPadding is a variable for width of the strip around the boundary 
         * of the ocean that gets sent to other oceans as "ghost" information.
	It gets set immediately before the call to qsort that sorts by escape 
         * region, so that subsequent calls to escapeRegion account for the 
         * padding or not, depending on whether tempPadding is 0 or... something 
         * positive.
	*/
	mutable double tempPadding;

public:

	/*bounds for the ocean, i.e. this is the rectangle that the ocean really
         * uses.  Everything outside this rectangle gets sent to another ocean*/
	double left, top, right, bottom;

/*	OCEAN CLASS MEMBER FUNCTIONS */

	int countFish(int fishInOcean); //counts number of fish in ocean
	void print() const;
	void printGhostFish() const;

	/*	clear sets numberOfFish to 0 and reallocates the block to */
	void clear();

	/* reallocate looks a sizeOfFish and numberOfFish and makes a decision 
         * about whether to reallocate fish.  If it does reallocate it does so 
         * with a call to realloc so if the memory happens to exist lying at the
         * end, realloc just expands the current block (b/c that's what realloc 
         * does)*/
	void reallocate();

	/*setSize deletes the current fish array and replaces it with a brand 
         * new one of size insize*/
	void setSize(int insize);

	/*	adds a fish to the ocean which mimmicks the record R*/
	void add(const FishRecord& R, bool inIsGhost);

	/* add adds a fish to the list adjusting the size somehow */
	void add(const Fish& F);

	/*	addAsGhost adds a fish and on the way sets the flag isGhost*/
	void addAsGhost(const Fish& F);

	/*	sets all fish in the array to position zero, velocity zero*/
	void zero();

	/*	sorts the fish in order of x-coordinate to make testing faster */
	void sortByx();

	/*	determines the first fish in the list whose x-coordinate is 
         * greater than or equal to M;
	it returns the index of the first fish satisfying this condition */
	int binarySearchxCoordinate(double M);

	/*helpful function for placing all fish randomly in a rectangle with 
         * random velocities*/
	void initializeFish(double displacementX, double displacementY, double amplitudeX, double amplitudeY, double speedLowerBound, double speedUpperBound, double inSelfWeight);

	/*	Assigns ID numbers to each fish in its fish array */
	void initializeFishID(int currentFishID);

	/*	This allows us to reset the speed midway through a simulation */
	void resetSpeed(double inSpeed);

	/* Sets initial speed*/
	void setSpeed(double inSpeed);

	/* Takes in an angle, amplitude of noise, and a bounding box and calls the function F.resetDirection() for fish within the bounding box */
	void resetDirection(double theta, double noiseAmp, double xMin, double xMax, double yMin, double yMax);


	/*copyVelocities calls copyVelocity for all fish (copyVelocity sets the "old" velocity variables equal to the current ones)*/
	void copyVelocities();

	/*	interact calls updateVelocity for all fish*/
	void interact();

	/*	addNoiseToDirectionAngle calls addNoiseToDirectionAngle for every fish in the ocean */
	void addNoiseToDirectionAngle(double noiseAmplitude);


	/*	move calls move for all fish*/
	void move(double timestep = TIMESTEP);
        
        //Update the reproduction status
        void setReproduction(int i, bool Reproduced);


	/* These functions are for the measures for swarming */
	Vector sumOfDirectionalHeadings(); //Computes the sum of the directional headings over all fish in the ocean
	Vector sumOfCoordinates();  // Computes the sum of the x- and y-coordinates for all fish in the ocean


	void findNearestGridPoint( int gridSizeX, int gridSizeY );
	void emptyLandOfFish();
	void initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary);
	void initializeFishPositions(Vector center1, double majax1, double minAx1, Vector center2, double majax2, double minAx2,
						Vector center3, double majax3, double minAx3, Vector center4, double majax4, double minAx4,
						Vector center5, double majax5, double minAx5, Vector center6, double majax6, double minAx6,
						Vector center7, double majax7, double minAx7, Vector center8, double majax8, double minAx8);

	/*	setBounds is an accessor */
	void setBounds(double inLeft, double inRight, double inBottom, double inTop);

	/*	inBounds returns whether or not the fish is really in the bounds of this ocean*/
	bool inBounds(const Fish& F) const;

	/*	finds all fish not in bounds and sends them away*/
	void sendFishOut();

	/*	removeMarkedFish runs through all the fish and removes them from the list if they have their removeMe flag set, and  */
	void removeMarkedFish();

	/*	removeGhostFish removes all fish which are ghosts (I hope this function is not actually necessary)*/
	void removeGhostFish();


//private:
	/*	escapeRegion assigns a number to each region around the current ocean.  Hopefully this picture explains the numbering:
		7 6 5
		8 0 4   <-- the 0 is the current ocean and the others are the oceans around it.
		1 2 3
	*/
	int escapeRegion(const Fish& F) const;

	/*	determines the first fish in the list of fish sorted by escaping ocean who escapes into ocean bigger than or equal to M+1;
		If there is no such fish, it returns a hack*/
	int binarySearchEscapeRegion(int M);

	/*	sorts the fish by escape region.  Takes an argument inPadding*/
	void sortByEscapeRegion(double inPadding);

public:
	//calls F.setDirection for all fish in ocean with identical parameters
	void setDirection(double angle);

	void makeDensityField(double** densityField, int sizeOfDesnityFieldX, int sizeOfDesnityFieldY );


	/*	randomFish initializes fish and makes all the fish have a random position and speed (according to the argument) and stuff*/
	void randomFish(double speedLowerBound, double speedUpperBound, double inSelfWeight);

	/*	alignedFish initializes fish and makes all the fish's directional headings equal */
	void alignedFish(double speedLowerBound, double speedUpperBound, double cosPhi, double sinPhi, double inSelfWeight);

	void initializeMaturityLevels(int initialMaturityLevel);

	/*	DEB functions */
	void initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime);
	void solveDEB(double h, double startTime, double endTime);		// This function solves the coupled DEB ODEs using 4th order Runge Kutta.

	FishRecord* fishRecordList;
	int constructList();
	void destroyList();
        
        //Adding functions to the main model
        void getOld(int TotalNumberFish, float toKillToday, float tJuv);
        void NaturalMortality(int type, int NaturalMortality, float tJuv, int nJuveniles, int nAdults, int nOld, float toKillToday);
        int numberOfFishToRepToday(int day);
        int isBorn(int iniID, int day, int nEggs);

        
        
}; //class Ocean



class World {
public:
	World(const World& W);
	World(int inSizeX, int inSizeY);
	World();
	~World();

/*	WORLD CLASS MEMBER VARIABLES */

	int headNode;
	int numberOfProcessors;

/*	WORLD CLASS MEMBER FUNCTIONS */

	void initFishFromFile(const char* filename); /*Reads in initial coordinates of fish from a file. In (x,y) format, one pair per line */

	void allocateOceans(int inSizeX, int inSizeY);
	void allocateProcessors(int inSize); // should only be called on the head node

	void deallocateOceans();
	void deallocateProcessors(); // should only be called on the head node
        
	/*	allocates and deallocates the grid */
	void allocateGrid(int x, int y);
	void deallocateGrid();
	void allocatePrefTempRanges();
	void deallocatePrefTempRanges();
	void zeroGrid(int xDimension, int yDimension);
	void findNearestGridPoint( int gridSizeX, int gridSizeY );


	/*   given a filename, an array name, and dimensions n and m, this function translates a list of numbers from a file into an n by m matrix */
	
        void translateTemperatureToGrid(const char* filename, int xDimension, int yDimension);
	void resetTemperatureGradientFlags();
	void translateStraumurToGrid(const char* filename, int xDimension, int yDimension);

	void initSize(int inSizeX, int inSizeY, double inOWidth, double inOHeight);

	void initMultiThread();

/*	WORLD CLASS MEMBER VARIABLES */

	int rank;
        
	bool isTorus;
	int sizeX, sizeY;
	double oWidth, oHeight;
	Ocean** oceans;
	Processor* processors;
	double padding;
	ExchangeData exchangeData;


    /* These variables are used in the measure functions */
    int howManyFish; //Stores the total number of fish  in the world
    int howManyFishToRep; //Stores the number of fish able to reproduce.
    int howManyFishToRepToday; //Stores the number of fish able to reproduce today.
    int nOfAdults; //Keep track of the adults...
    int LastUsedIDFish; //Last ID used assigned to a fish
    int AccumulatedLandings; //Self explanatory....
    int YearlyAccumulatedLandings; //This is reset to 0 each year!
    float MaxDensAdults;       //Used to calculate the number of boats per pixel.
    float totalBoatPP;
    float MaxNBoats; // Keep track of the pixel with more boats
    float TAC;
    float TACpp;
    
    Vector averageDirection; //Stores the average directional heading over all the fish in the world
    Vector centerOfMass; //Stores the center of mass of the fish in the world

    	/*	These are the results of the measurements defined for swarming solutions */
	double alignmentFactor;
	double averageSpeed;
	double polarity;
	double bouncingFactor;


     /* These are all for drawing pictures */
	int totalFishFromAllProcessors;
	int sizeOfImageX;
	int sizeOfImageY;
	int sizeOfDensityFieldX;
	int sizeOfDensityFieldY;
	TGAImage picture;

	double* densityFieldData;
	double** densityField;


    	/*	This is all useful for the creation of the big grid which holds info
		about the temperature and the currents */
	int sizeGridsX;
	int sizeGridsY;
	GridPoint** grid;
	
	/* Temperature preference is [tooCold,tooHot] for each maturity level of the fish 
		Initialized in World::initializeTempPrefRanges()
		See #defined temperature ranges above */
	double* tooCold;
	double* tooHot;


/*	WORLD CLASS MEMBER FUNCTIONS */

    /*  An initialization function puts randomly-located random-velocity fish into each ocean*/
	void randomFish(int inNumberOfFish, double speedLowerBound, double speedUpperBound, double inSelfWeight);

	/*  An initialization function puts randomly-located random-speed, directionally aligned fish into each ocean*/
	void alignedFish(int inNumberOfFish, double speedLowerBound, double speedUpperBound, double cosPhi, double sinPhi, double inSelfWeight);


	/*	these two functions get the info from GridPoints that's needed to update the fish's directions */
	double getTemperature(IntegerVector nearestGridPoint);
	Vector getTemperatureGradient( double normTempSens, IntegerVector nearestGridPoint, int quadrant, int fishMaturityLevel);
	Vector getTemperatureGradient2(double normTempSens, IntegerVector nearestGridPoint, int fishMaturityLevel);


	/*	this function computes the temperature part of the comfort function for grid[x][y] given upper and
		lower comfort bounds on the temperature */
	double R(int x, int y, double T1, double T2);
	void move();
        
        //Function to update the reproduction status
        void setReproduction(int i, bool Reproduced);
        
	void emptyLandOfFish();
	void initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary);
	void initializeFishPositions(Vector center1, double majax1, double minAx1, Vector center2, double majax2, double minAx2,
						Vector center3, double majax3, double minAx3, Vector center4, double majax4, double minAx4,
						Vector center5, double majax5, double minAx5, Vector center6, double majax6, double minAx6,
						Vector center7, double majax7, double minAx7, Vector center8, double majax8, double minAx8);

	void translateByStraumur(double timestep, int gridSizeX, int gridSizeY);
	void factorInTemperatureGradients(double normTempSens, int gridSizeX, int gridSizeY);
	//void initializePrefTempRanges();
        void initializePrefTempRanges(float TOOCOLD_IMMATURE, float TOOHOT_IMMATURE, float TOOCOLD_MATURE, float TOOHOT_MATURE);

	/*  sendFishTo copies fish from the array fishArray in the range of indices from placeholderStart to placeholderNextStart-1
        It is a local function i.e. does not call any mpi code*/
	int sendFishTo(int targetX, int targetY, Fish* fishArray, int placeholderStart, int placeholderNextStart, bool inIsGhost = false);
	int sendFishTo(Ocean& targetOcean, Fish* fishArray, int placeholderStart, int placeholderStartNext, bool inIsGhost = false);

	//calls O.setDirection for each ocean using identical parameters. used to easily align all fishes direction/speed
	void setDirection(double angle);

	void resetInteractionCounters();
	void transportFish();
	void sendFishOut();
	void communicateGhostFish();
	void removeGhostFish();

     //Calls interact for each Ocean
	void interact();

	/*	This allows us to reset the speed midway through a simulation */
	void resetSpeed(double inSpeed);

	/*sets initial speed for  different  oceans*/
	void setSpeed(double inSpeed);

	/*	This calls O.resetDirection, which takes in an angle, amplitude of noise, and a bounding box and calls the function F.resetDirection() for fish within the bounding box */
	void resetDirection(double theta, double noiseAmp, double xMin, double xMax, double yMin, double yMax);

	/*	addNoiseToDirecitonAngle calls addNoiseToDirectionAngle for every ocean */
	void addNoiseToDirectionAngle(double noiseAmplitude);


     /* Functions to do with the drawing of fish density into a file:*/

	void allocatePicture(double scaling);
	void allocatePicture(int width, int height);
	void allocateDensityField(double scaling);
	void allocateDensityField(int width, int height);
	void deallocateDensityField();

	void initPictureSerial();
	void initPictureHeadNode();
	void initPictureClientNode();

	void densityPointillism(Fish& fish);
	void drawPicture(const char* filename, int gridSizeX, int gridSizeY);
        void drawPictureAdults(const char* filename, int gridSizeX, int gridSizeY);
	void writePictureToFile(const char* filename);

	void copyDensityFieldIntoPicture(int offsetX, int offsetY, int thread);
	void sendTotalNumberOfFish();
	void getTotalFishFromAllProcessors();
	void sendDensities();
	void receiveDensityForOcean(const Ocean& O, int x, int y);
	void getPicture();
	void clearDensity();

	/*  Functions to do with creating a measure on our model */

	void findAverageDirectionalHeading();  //Computes the average directional heading over all the fish in the world
	void findAlignmentFactor();

	void findCenterOfMass();
	//void findAverageAngularVelocity();
	//void findPolarity();
	//void findBounceFactor();

	void initializeMaturityLevels(int initialMaturityLevel);

	/*	DEB functions */
	void initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime);
	void solveDEB(double h, double startTime, double endTime);		// This function solves the coupled DEB ODEs using 4th order Runge Kutta.

	/*  Basic utility functions */

	void totalNumberOfFish(); //This function updates the member variable howManyFish; it's called at the end of W.transportFish()
	void initializeFishID(); // Assigns a positive integer to each fish (from being -1) - works in serial for now
	void printNumbers() const;
	void printWhoLivesWhere() const;
	int tagForTask(int srcI, int srcJ, int flag, int iscount );
	void clearShadowOceans();


	/*  For dynamic load balancing.*/

	void printInteractionCounter() const;
	int totalInteractionCounter() const;
	int biggestLoadedProcessor();
	Processor smallestLoadedNeighbor(Processor& processorHigh);
	bool findOptimalOceanForExchange(Processor& processorHigh, Processor& processorLow, int& outX, int& outY);

	int runHeadNode();
	void sendLoadToHead();
	void getAllOceanLoads();
	int updateHeadNode();
	void updateMyWorld();
	void MPI_SendFish(int, int, int, int, int);
	void MPI_RecvFish(int, int, int, int, int);
	void constructExchangeData(Processor processorHigh, Processor processorLow);

        int addFish(int TotalFishCount, int day, int nEggs);
        void totalToRepToday(int day);
        void FishAging(int NaturalMortality, float toKillToday, float tJuv);
        void translateFtoGrid (int xDimension, int yDimension, float F_Mortality);
        int GetDensities(int gridSizeX, int gridSizeY);
        void FishingMortality(int i, int gridSizeX, int gridSizeY);
        
        void WriteXYDensity(char file[255], int gridSizeX, int gridSizeY, char sst[20], 
                        char cur[20], char interact[20], char manag[20], char _managX[20], char _managY[20]);
        void WriteGridFishIDs(char file[255], int day, int gridSizeX, int gridSizeY, double minMPAx, double maxMPAx, double minMPAy, double maxMPAy);
        void IniFishMovement();
        void printFish(int ID);
        void BoatsDist(int i, int gridSizeX, int gridSizeY, double minMPAx, 
                double minMPAy, double maxMPAx, double maxMPAy, int YearOfMPA, double boatsC);
        void BoatsDistNCatches(int i, int gridSizeX, int gridSizeY, double minMPAx, 
                double minMPAy, double maxMPAx, double maxMPAy, int YearOfMPA, float NCatches, double boatsC);
        
        double checkTAC(int i, int nOfAdults, bool tacFlag);
        float adjustLarvaeMortality(int type, int nOfAdults, int nOfJuv);
        double mpaDef(int mpaDef, int mpaSize);
        rVector Rparam(int HowManyFish, int nOfAdults);
        void writeRparam(int year, int day, int howManyFish, int nOfAdults, rVector R);
        

}; //class World



class Processor {
public:
	int thread;
	int load;
};


