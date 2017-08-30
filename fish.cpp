#include "fish.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tga.h"
#include <strings.h>
#include <sys/time.h>
#include <fstream>

#if MPI_ON  == 1
        #include "mpi.h"
#endif

const double pi = M_PI;  


/* This seeds the function srandom() so that it's different every run (\mu a.e.) */
void seedRandom()
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	srandom(tv.tv_usec);
}

/* This returns a randomVariable with value between -1 and 1 */
double randomVariable()
        {
	/*  the maximum value of random() is (2**31)-1, so we scale down to something in the range from -1 to 1
		and then multiply by amplitude to get a random number of the right size */
	double rand = (double)random()/((0x40000000)-1)-1; //NB 0x40000000 is hex for 2**30
	return rand;
        }

//This function describe how the probability of reproduction is assigned for any given day
int probOfReproduction(int type)
        {
        int dayOfRep;
        if (type == 0)
                {
                dayOfRep = rand()%366;
                }
        else printf("Not yet implemented!!");
        return dayOfRep;
        }

int probOfReproduction(int type, double mu, double std)
        {
        int dayOfRep;
        if (type == 1)
        {
                double u =(double)(random() %100000 + 1)/100000; //for precision
                double v =(double)(random() %100000 + 1)/100000; //for precision
                double x = sqrt(-2*log(u))*cos(2*pi*v);   //or sin(2*pi*v)
                double y = x * std + mu;
                dayOfRep = y;
        }
        else printf("Not yet implemented!!");

        return dayOfRep;
        }

//This is a recursive function to make sure the day stay between 1 and 365
    int JulDay(int day) // 
        {
        if (day > 365) {day = JulDay(day - 365);}
        return day;
        }

/* Instantiates a Fish: */
Fish::Fish()
{
	ID = -1; //Should be assigned a positive integer later
	removeMe = false;
	isGhost = false;
	quadrant = 0;
	updateWeights(DEFAULT_INTERACTIONWEIGHT, DEFAULT_TEMPERATUREWEIGHT, DEFAULT_SELFWEIGHT);
	makeRadii();
	zero();
}

void Fish::print() const
{
	printf( "p:(%5.3f, %5.3f) d:(%5.3f %5.3f) s:%5.3f  :  %d\n", x,y, cosPhi, sinPhi, speed, quadrant );
}

/* NB has not been tested */
void Fish::printToFile() const
{
	FILE* file;
	char positionInfo[20];
	positionInfo[0]=0;
	sprintf( positionInfo, "%8.6f %8.6f %d\n", x, y, (int)removeMe );
	char filename[255];
	filename[0] = 0;
	sprintf( filename, "positionInfo_fish%d.dat", ID );
	file = fopen(filename, "a+b");
	if( file==NULL )
	{
		printf( "file %s failed to open for writing.\n", filename );
	}
	fwrite(positionInfo, 1, sizeof(positionInfo), file);
	fclose(file);

	printf("fish.%d is calling printToFile() \n", ID);
}



/*	this function finds the nearest gridpoint on the grid which stores the straumur and temperature */
void Fish::findNearestGridPoint(int gridSizeX, int gridSizeY)
{
 	nearestGridPoint.x = (int)(x + 0.5);
	nearestGridPoint.y = (int)(y+0.5);

	//printf("(F.x,F.y) = (%f, %f) BUT (nearestGridPoint.x, nearestGridPoint.y) = (%d, %d) \n", (double)x, (double)y, nearestGridPoint.x, nearestGridPoint.y);

	if( nearestGridPoint.x < 0 ) nearestGridPoint.x=0;
	if( nearestGridPoint.x > gridSizeX-1 ) nearestGridPoint.x=gridSizeX-1;

	if( nearestGridPoint.y < 0 ) nearestGridPoint.y=0;
	if( nearestGridPoint.y > gridSizeY-1 ) nearestGridPoint.y=gridSizeY-1;
}

/*	this function finds the quadrant of the nearest grid point where the fish is located
	Doesn't yet work for isTorus=true */
void Fish::findQuadrant()
{
	if (y > nearestGridPoint.y)
	{
		if ( x > nearestGridPoint.x)
			quadrant = 0;
		else
			quadrant = 1;
	}
	else
	{
		if ( x < nearestGridPoint.x)
			quadrant = 2;
		else
			quadrant = 3;
	}
}


void Fish::updateWeights(double iWeight, double tWeight, double sWeight)
{
	interactionWeight = iWeight;
	temperatureWeight = tWeight;
	selfWeight = sWeight;
}


void Fish::makeRadii()
{
	/*  For the moment I just want them to swim in random directions*/
	radiusOfRepulsion  =  DEFAULT_RADIUS_OF_REPULSION;
	radiusOfOrientation = DEFAULT_RADIUS_OF_ORIENTATION;
	radiusOfAttraction =  DEFAULT_RADIUS_OF_ATTRACTION;
}



void Fish::zero()
{
	x=y=speed=cosPhi=sinPhi=0.0;
	nearestGridPoint.x=nearestGridPoint.y=0;
}


/* Randomizes the position of a Fish in terms of global World coordinates */
void Fish::randomPosition(double displacementX, double displacementY, double amplitudeX, double amplitudeY)
{
	/*	The center of the Ocean that the Fish lives in is(displacementX,displacementY)
		The distances from the center to the boundaries of the Ocean are amplitudeX and amplitudeY*/
	x = amplitudeX*randomVariable() + displacementX;
	y = amplitudeY*randomVariable() + displacementY;
}

void Fish::randomDirection()
{
	double angle = pi*(randomVariable()+1);	/* Uniform [0,2*pi] */
	cosPhi = cos(angle);
	sinPhi = sin(angle);
}

void Fish::randomSpeed(double lowerBound, double upperBound)
{ //Gives a uniform distribution in [lowerBound, upperBound]
//	speed = amplitude*((double)random()/(2*(0x40000000)-1));
	speed = lowerBound + (upperBound-lowerBound)*fabs(randomVariable());
}


void Fish::setPosition(double inx, double iny)
{
	x=inx;
	y=iny;
}

void Fish::setDirection(double angle)
{
        cosPhi = cos(angle*(pi/180));
	sinPhi = sin(angle*(pi/180));
}

void Fish::setSpeed(double inSpeed)
{
	speed=inSpeed;
}

void Fish::setVelocity(double inSpeed, double inCosPhi, double inSinPhi, double inSelfWeight)
{
	speed = inSpeed;
	cosPhi = inCosPhi;
	sinPhi = inSinPhi;
	selfWeight = inSelfWeight;
}

void Fish::resetSpeed(double inSpeed)
{
	//if(x>68) //This corresponds to 13.5 West
	printf("Fish.%d.speed is resetting! From Feb through July 2008, this only was reset if F.x > 68 \n", ID);
	{	speed=inSpeed;
	}
}

/* This function resets the direction of a fish to theta with a random perturbation of amplitude noiseAmp */
void Fish::resetDirection(double theta, double noiseAmp)
{
	double direction = theta+noiseAmp*randomVariable();
	cosPhi = cos(direction);
	sinPhi = sin(direction);
}


void Fish::setVelocity(double vx, double vy, double inSelfWeight)
{
    speed = sqrt(vx*vx+vy*vy);
    if(speed < DECISION_TOL)
        cosPhi = sinPhi = 0.0;
    else
    {
        cosPhi = vx/speed;
        sinPhi = vy/speed;
    }
    selfWeight = inSelfWeight;
}



void Fish::copyVelocity()
{
	oldCos = cosPhi;
	oldSin = sinPhi;
	oldSpeed = speed;
}


/* Calculates the distance squared from this Fish to Fish F */
double Fish::distanceToSquared(const Fish& F)
{
	return (F.x-x)*(F.x-x) + (F.y-y)*(F.y-y);
}

double Fish::getPrefSpeed()
{
	//Fish keeps track of its roe percentage: F.roePercentage
	//NB The numbers here are set manually!
	double prefSpeed;
	if(roePercentage < DEB_ROE_PERCENTAGE_MARK) //About 8% roe
		prefSpeed = DEFAULT_SPEED_NO_ROE; // km/day 
	else if(roePercentage < DEB_ROE_PERCENTAGE_MAX) //Abotu 25%
		prefSpeed = DEFAULT_SPEED_NO_ROE + MAX_SPEED_INCREASE*(roePercentage-DEB_ROE_PERCENTAGE_MARK)/(DEB_ROE_PERCENTAGE_MAX - DEB_ROE_PERCENTAGE_MARK); // km/day
	else
		prefSpeed = DEFAULT_SPEED_NO_ROE + MAX_SPEED_INCREASE; // km/day
		
	prefSpeed /= 12.0; //Done to switch to grid units
	return prefSpeed;
}


void Fish::updateVelocity()
{
  	if( !isGhost and age >= AGESELFMOVE)
 	{

	double averageCos = selfWeight*oldCos;
	double averageSin = selfWeight*oldSin;
	double averageSpeed = oldSpeed;
	double radiusOfAttractionSquared = radiusOfAttraction*radiusOfAttraction;
	double neighborCounter = selfWeight;

	int orientationCounter = 1;

	for(int i=myOcean->binarySearchxCoordinate(x-radiusOfAttraction); i < myOcean->numberOfFish; i++)
	{
            
		Fish& neighbor(myOcean->fish[i]);

		if(&neighbor==this)
			continue;
                
                if (neighbor.age < AGESELFMOVE) //Avoid the juveniles to affect the movement of the adults!!)
                    break; 
                
		if(neighbor.x > x+radiusOfAttraction)
			break;

		if(neighbor.y > y+radiusOfAttraction || neighbor.y < y - radiusOfAttraction)
			continue;

		double proximitySquared = distanceToSquared(neighbor);
		if (proximitySquared < radiusOfAttractionSquared)	// if its in the radius of attraction we'll consider it
		{
			double proximity = sqrt(proximitySquared);	// first take the distance to the neighbor fish
			neighborCounter = neighborCounter+1;					// count it as a neighbor to be used in averaging later
			if(proximity < radiusOfOrientation)		// determine if its in the radius of orientation
			{
				averageSpeed += neighbor.oldSpeed;
				orientationCounter++;			// only things in the zone of orientation affect the speed, so they get counted separately
				if(proximity < radiusOfRepulsion)
				{					// if its in the radius of repulsion, we add to the average velocity in the direction pointed away from the neighbor...
					if( proximity > DECISION_TOL)
					{
						averageCos += (x-neighbor.x)/proximity;
						averageSin += (y-neighbor.y)/proximity;
					}
					continue;			// ... and that's all we do for this neighbor
				}
				averageCos += neighbor.oldCos;          // if it's in the zone of orientation, but not the zone of repulsion, then we add the velocity of the neighbor to the average velocity...
				averageSin += neighbor.oldSin;
				continue;						// ... and that's that
			}

			averageCos += (neighbor.x-x)/proximity;         //  if it's outside the zone of orientation, but inside the zone of attraction, we add to the average velocity in the direction pointed toward from the neighbor
			averageSin += (neighbor.y-y)/proximity;
		}
	}

	averageCos /= neighborCounter;
	averageSin /= neighborCounter;

	double normOfAveDirection = sqrt(averageCos*averageCos + averageSin*averageSin);

	if(normOfAveDirection < DECISION_TOL)
	{
		//printf("uh oh, at least one inf \n");
		cosPhi = averageCos;
		sinPhi = averageSin;
	}
	else
	{
		cosPhi = averageCos/normOfAveDirection; // (cosPhi, sinPhi) should be a unit vector, so we divide by the norm.
		sinPhi = averageSin/normOfAveDirection;
	}
	//Calculate the speed of neighboring particles:
	averageSpeed /= orientationCounter; //nb orientation counter is at least 1
	
	//Finally, set the found speed:
	#if DEB_ON
		double prefSpeed;
		prefSpeed = getPrefSpeed();
		speed = (1-ALPHA)*averageSpeed + ALPHA*prefSpeed;
	#else
		//Just average the interactions:
		speed = averageSpeed;
	#endif

	myOcean->interactionCounter += int(neighborCounter - selfWeight);
       
	} // Ending the isGhost check
	
} // Ending Fish::updateVelocity



/*	This function adds noise of the given amplitude to the direction angle of the fish */
void Fish::addNoiseToDirectionAngle(double noiseAmplitude)
{
	double R = noiseAmplitude*randomVariable();
	double cosR = cos(R);
	double sinR = sin(R);

	sinPhi = sinPhi*cosR + cosPhi*sinR;
	cosPhi = cosPhi*cosR - sinPhi*sinR;
}

void Fish::move(double timestep)
        {
	x += cosPhi*speed*timestep;
	y += sinPhi*speed*timestep;
        }

void Fish::emptyLandOfFish()
{
	if( myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].temperature > 998)
	{
		removeMe = true;
	}

if(myOcean->myWorld->getTemperature(nearestGridPoint)>998)
	removeMe=true;
}


void Fish::initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary)
{
	if( myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].temperature < 998)
		if(nearestGridPoint.x < leftBoundary || nearestGridPoint.x > rightBoundary
			|| nearestGridPoint.y < lowerBoundary || nearestGridPoint.y > upperBoundary)
				removeMe = true;
}

void Fish::initializeFishPositions(Vector center1, double majAx1, double minAx1, Vector center2, double majAx2, double minAx2,
						Vector center3, double majAx3, double minAx3, Vector center4, double majAx4, double minAx4,
						Vector center5, double majAx5, double minAx5, Vector center6, double majAx6, double minAx6,
						Vector center7, double majAx7, double minAx7, Vector center8, double majAx8, double minAx8)
{
	/* removeMe is set false if the fish are in "safe zones" */
	if(	minAx1*minAx1*(x-center1.x)*(x-center1.x) + majAx1*majAx1*(y-center1.y)*(y-center1.y) < minAx1*minAx1*majAx1*majAx1 ||
		minAx2*minAx2*(x-center2.x)*(x-center2.x) + majAx2*majAx2*(y-center2.y)*(y-center2.y) < minAx2*minAx2*majAx2*majAx2 ||
		minAx3*minAx3*(x-center3.x)*(x-center3.x) + majAx3*majAx3*(y-center3.y)*(y-center3.y) < minAx3*minAx3*majAx3*majAx3 ||
		minAx4*minAx4*(x-center4.x)*(x-center4.x) + majAx4*majAx4*(y-center4.y)*(y-center4.y) < minAx4*minAx4*majAx4*majAx4 ||
		minAx5*minAx5*(x-center5.x)*(x-center5.x) + majAx5*majAx5*(y-center5.y)*(y-center5.y) < minAx5*minAx5*majAx5*majAx5 ||
		minAx6*minAx6*(x-center6.x)*(x-center6.x) + majAx6*majAx6*(y-center6.y)*(y-center6.y) < minAx6*minAx6*majAx6*majAx6 ||
		minAx7*minAx7*(x-center7.x)*(x-center7.x) + majAx7*majAx7*(y-center7.y)*(y-center7.y) < minAx7*minAx7*majAx7*majAx7 ||
		minAx8*minAx8*(x-center8.x)*(x-center8.x) + majAx8*majAx8*(y-center8.y)*(y-center8.y) < minAx8*minAx8*majAx8*majAx8 )
			removeMe = false;
	else
	{
			removeMe = true;
	}
}

int Fish::escapeRegion()
{
	return myOcean->escapeRegion(*this);
}


void Fish::mimmicRecord(const FishRecord& R)
{
	x = R.x;
	y = R.y;
	cosPhi = R.cosPhi;
	sinPhi = R.sinPhi;
	speed = R.speed;
	selfWeight = R.selfWeight;
}

void Fish::makeRecord(FishRecord& R) const
{
	R.x = x;
	R.y = y;
	R.cosPhi = cosPhi;
	R.sinPhi = sinPhi;
	R.speed = speed;
	R.selfWeight = selfWeight;
}


void Fish::initializeMaturityLevels(int initialMaturityLevel)
{
	maturityLevel = initialMaturityLevel;
}


/*	DEB functions: */

/*	Initiates the variables for the DEB model */
void Fish::initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime)
{
	l = l0;
	e = e0;
	uR = uR0;
	er = er0;
	maturityLevel = initialMaturityLevel;
	daysSinceMature = initialDaysSinceMature;
	
	double Wroe = (DEB_U2E/DEB_RHO_ROE) * er;
	double Lm3 = DEB_LM*DEB_LM*DEB_LM;
	double l3 = l*l*l;   
	double V = Lm3*l3;
	double Weight = V*DEB_DV + (DEB_U2E/DEB_RHO_E)*(l3*e + uR + er) + Wroe;
		
	roePercentage = 100*Wroe/Weight;


	if(ID%CHOOSE_EVERY_NTH_FISH_TO_PRINT_DEB == 0)
	{ 
		printInitialBiologyToFile();
	}
}


/*	This function solves the coupled DEB ODEs using 4th order Runge Kutta with time step h.
	Functions when h<(endTime-startTime)
	It saves a copy of the current length, energy reserves, reproduction buffer and roe energy so they won't be overwritten while needed.
	At the end of the function, it updates F.weight, etc */
void Fish::solveDEB(double h, double startTime, double endTime)
{
	//printf("ID of fish = %d \n",ID);

	double length = l;
	double energy = e;
	double repBuffer = uR;
	double roeEn = er;
	double t = startTime;

	double k11,k12,k13,k14;
	double k21,k22,k23,k24;
	double k31,k32,k33,k34;
	double k41,k42,k43,k44;

	bool flaggy = myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].tempCorrectionFlag;
	if (!flaggy)
	{	printf("Holy shit, %d %9.6f %9.6f %9.6f\n",ID,x,y,endTime);
		printf("Temp in nearest grid point: %9.6f\n",myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].temperature);
		exit(0);
	}
		
	double nu    = myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].tempCorrected_nu;
	double kJ    = myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].tempCorrected_kJ;
	double gamma = myOcean->myWorld->grid[nearestGridPoint.x][nearestGridPoint.y].tempCorrected_gamma;

	while (t+h <= endTime)
	{
		/* Here comes 4th order Runge-Kutta.  Summary:
			y_prime = f(t,y),  y(t_0) = y_0
			Then:
			k_1 = f(t_n, y_n)
			k_2 = f(t_n + h/2, y_n + (h/2)*k1)
			k_3 = f(t_n + h/2, y_n + (h/2)*k2)
			k_4 = f(t_n + h, y_n + h*k_3)

			y_{n+1} = y_n + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)
			t+{n+1} = t_n + h
		*/

		k11 = dl( t, length, energy,nu);
		k12 = de( t, length, energy,nu);
		k13 = duR(t, length, energy,nu,kJ);
		k14 = der(t,                 repBuffer, roeEn,gamma);

		k21 = dl( t+h/2, length+(h/2)*k11, energy+(h/2)*k12,nu);
		k22 = de( t+h/2, length+(h/2)*k11, energy+(h/2)*k12,nu);
		k23 = duR(t+h/2, length+(h/2)*k11, energy+(h/2)*k12,nu,kJ);
		k24 = der(t+h/2,                                     repBuffer+(h/2)*k13, roeEn+(h/2)*k14,gamma);

		k31 = dl( t+h/2, length+(h/2)*k21, energy+(h/2)*k22,nu);
		k32 = de( t+h/2, length+(h/2)*k21, energy+(h/2)*k22,nu);
		k33 = duR(t+h/2, length+(h/2)*k21, energy+(h/2)*k22,nu,kJ);
		k34 = der(t+h/2,                                    repBuffer+(h/2)*k23, roeEn+(h/2)*k24,gamma);

		k41 = dl( t+h, length+h*k31, energy+h*k32,nu);
		k42 = de( t+h, length+h*k31, energy+h*k32,nu);
		k43 = duR(t+h, length+h*k31, energy+h*k32,nu,kJ);
		k44 = der(t+h,                            repBuffer+h*k33, roeEn+h*k34,gamma);

		length +=		(h/6)*(k11 + 2*k21 + 2*k31 + k41);
		energy +=		(h/6)*(k12 + 2*k22 + 2*k32 + k42);
		repBuffer +=	(h/6)*(k13 + 2*k23 + 2*k33 + k43);
		roeEn +=		(h/6)*(k14 + 2*k24 + 2*k34 + k44);
		t += h;
	}

	// Now it's ok to overwrite these variables, since the functions which need the old
	//info have already been computed
	l = length;
	e = energy;
	uR = repBuffer;
	er = roeEn;
	
	if (maturityLevel == IMMATURE)
	{// Check whether the roe percentage has exceeded 8% -
		double Wroe = (DEB_U2E/DEB_RHO_ROE) * er;
		double Lm3 = DEB_LM*DEB_LM*DEB_LM;
		double l3 = l*l*l;   
		double V = Lm3*l3;
		double Weight = V*DEB_DV + (DEB_U2E/DEB_RHO_E)*(l3*e + uR + er) + Wroe;
		
		roePercentage = 100*Wroe/Weight;

		// NB The roe percentage is 100*Wroe/Weight;
		if( roePercentage > DEB_ROE_PERCENTAGE_MARK) // here denoting a mature fish
		{
			maturityLevel = MATURE;
		}
	}
	else
	{// We need to account for a higher water content of the roe
		double rho_r;
		//double time2max = DEB_TIME2MAX; //days
		//double rhoMaxIncrease = DEB_ROE_MAX_WATERINCREASE;// percentage
		daysSinceMature += endTime-startTime;

		if ( DEB_TIME2MAX <= daysSinceMature)
		{//Take to be the maximum
			rho_r = DEB_RHO_ROE/( 1 + DEB_ROE_MAX_WATERINCREASE );
		}
		else
		{//We increase by some value:
			rho_r = DEB_RHO_ROE/( 1 + DEB_ROE_MAX_WATERINCREASE*(daysSinceMature/DEB_TIME2MAX) );
		}
		//And now, calculate the roe percentage:
		double Wroe = (DEB_U2E/rho_r) * er;
		//double Lm3 = DEB_LM*DEB_LM*DEB_LM;
		double l3 = l*l*l;   
		//double V = Lm3*l3;
        double V = (DEB_LM*DEB_LM*DEB_LM)*l3;
		double Weight = V*DEB_DV + (DEB_U2E/DEB_RHO_E)*(l3*e + uR + er) + Wroe;
		roePercentage = 100*Wroe/Weight;
	}

	// This might slow down the code a bit
	if(ID%CHOOSE_EVERY_NTH_FISH_TO_PRINT_DEB == 0 && ((int)(endTime*100))%100 == 0) // writes out every day
	{
		printBiologyHistoryToFile(endTime);
	}
}//void Fish::solveDEB


//NB nu needs to be temperature dependent, should be done when called!
double Fish::dl(double t, double length, double energy, double nu)
{
	double l_prime = nu/(3*DEB_LM) * (energy - length)/(energy + DEB_G);
	if (l_prime<0)
		return 0;
	else
		return l_prime;
}

//NB nu needs to be temperature dependent, should be done when called!
double Fish::de(double t, double length, double energy, double nu)
{
	double f;//=0.5*(cos(t/365.0*2*pi)+1.0)/2.0+0.5;
	double daysWithFood = 30.0;
	if (t<daysWithFood)
		f = 0.5*( 1.0 - t/daysWithFood );	// This is food being eaten or food density or something
	else
		f=0;

	double e_prime = nu/(length*DEB_LM) * (f - energy);

	return e_prime;
}

//NB nu and kJ need to be temperature dependent, should be done when called!
/*	This includes the temperature of the nearest gridPoint as the temperature */
double Fish::duR(double t, double length, double energy, double nu, double kJ)
        {
	double uR_prime = nu/DEB_LM * (1-DEB_KAPPA)* energy*length*length * (length + DEB_G) / (energy + DEB_G) - kJ*DEB_UHP;
	return uR_prime;
        }

//NB nu needs to be temperature dependent, should be done when called!
double Fish::der(double t, double repBuffer, double roeEn, double gamma)
{	
	double er_prime = gamma * (repBuffer - roeEn) * roeEn;
	return er_prime;
}

// Begins the .m file with a variable name, etc.
void	Fish::printInitialBiologyToFile()
{
	FILE* file;

	char biologyData[90];
	biologyData[0]=0;
	sprintf( biologyData, "%9.6f %9.6f %10.6f %10.6f %10.6f %10.6f %11.6f %d %11.6f\n", l, e, uR, er, roePercentage, speed, 0.0, maturityLevel, daysSinceMature );

	char filename[255];
	filename[0] = 0;
	sprintf( filename, "DEBoutput/biologicalHistory_fish%d.dat", ID );

	char initialComment[202]; // length of following string
	initialComment[0]=0;
	sprintf( initialComment, "%% This file provides biological history for fish with ID given in file name\n%% columns record, from left to right: \n%% l, e, u_R, e_r, roePercentage, speed, time in days, maturityLevel, daysSinceMature\n\n");

	file = fopen(filename, "wb");
	if( file==NULL )
	{
		printf( "file %s failed to open for writing.\n", filename );
	}
	fwrite(initialComment, 1, sizeof(initialComment), file);
	fwrite(biologyData, 1, sizeof(biologyData), file);
	fclose(file);

	//printf("fish.%d is calling printInitialBiologyToFile() \n", ID);
}

// Prints the weight, fat, roe, time into an array of numbers which matlab can read
void Fish::printBiologyHistoryToFile(double time)
{
	FILE* file;

	char biologyData[90];
	biologyData[0]=0;
	sprintf( biologyData, "%9.6f %9.6f %10.6f %10.6f %10.6f %10.6f %11.6f %d %11.6f\n", l, e, uR, er, roePercentage, speed, time, maturityLevel, daysSinceMature );

	char filename[255];
	filename[0] = 0;
	sprintf( filename, "%s/Beta_%f/SST_%d_%d/DEBoutput/biologicalHistory_fish%d.dat", MAINDIR, BETA, SSTmin, SSTmax, ID );

	file = fopen(filename, "a+b");
	if( file==NULL )
	{
		printf( "file %s failed to open for writing.\n", filename );
	}
	fwrite(biologyData, 1, sizeof(biologyData), file);
	fclose(file);

	//printf("fish.%d is calling printBiologyHistoryToFile() \n", ID);
}


// --------------------------------------------------//
/* Here are all the functions that I added to *
 * the class Fish                                   */
// --------------------------------------------------//
// This function create a new fish when the prpulation reproduce


// This function allows me to update the reproductive status of each fish! 
void Fish::setReproduction(int k, bool Reproduced)
        {
        if (age >= AGEMATURE) //365/TIMESTEP)
                {
                isAdult = true;
                }
        }


/* En of Fish class*/

GridPoint::GridPoint()
{
	flagQuadrant = new bool[4];
	temperatureGradientAtCenter = new Vector[NUMMATURITYLEVELS];
	temperatureGradAtCenterFactor = new double[NUMMATURITYLEVELS];
	temperatureGradient = new Vector*[4];
	temperatureGradFactor = new double*[4];
	for (int i=0; i<4; i++)
	{
		temperatureGradient[i]   = new Vector[NUMMATURITYLEVELS];
		temperatureGradFactor[i] = new double[NUMMATURITYLEVELS];
	}
	zero();}


GridPoint::~GridPoint()
{
	delete[] flagQuadrant;
	//These guys depend on the maturity levels:
	delete[] temperatureGradientAtCenter;
	delete[] temperatureGradAtCenterFactor;
	//The rest is [quadrants][maturity levels]
	for (int j=0; j<4; j++)
	{//First, take care of the maturity levels:
		delete[] temperatureGradient[j];
		delete[] temperatureGradFactor[j];
	}
	//And then the quadrants:
	delete[] temperatureGradient;
	delete[] temperatureGradFactor;
}


void GridPoint::zero()
{
	temperature = 0.0;
        density = 0; 
        AdultDensity = 0; 
        FM = 0;   
        FMpp = 0; //Number of fish that can potentially be catch!
        FCatches = 0; //Number of fish caught by the fishery.
        AnnualCatches = 0; 
        boats = totOfBoats/(GRIDSIZEX * GRIDSIZEY);      // Number of boats to allocate in the each pixel of the grid.
        BoatsXFMpp = 0;
        bp = 0;
        
	//Used if DEB is on:
	tempCorrectionFlag = false;
	temperatureCorrection = 1.0;
	tempCorrected_nu = 0.0;
	tempCorrected_kJ = 0.0;
	tempCorrected_gamma = 0.0;

	flagQuadrant[0] = false;
	flagQuadrant[1] = false;
	flagQuadrant[2] = false;
	flagQuadrant[3] = false;

	gradientFlag = false;

	for (int j=0; j<NUMMATURITYLEVELS; j++)
	{
		//Note the quadrants:
		temperatureGradient[0][j].x = temperatureGradient[0][j].y = 0;
		temperatureGradient[1][j].x = temperatureGradient[1][j].y = 0;
		temperatureGradient[2][j].x = temperatureGradient[2][j].y = 0;
		temperatureGradient[3][j].x = temperatureGradient[3][j].y = 0;

		temperatureGradFactor[0][j] = 0;
		temperatureGradFactor[1][j] = 0;
		temperatureGradFactor[2][j] = 0;
		temperatureGradFactor[3][j] = 0;
		
		//Here there are no quadrants:
		temperatureGradientAtCenter[j].x = 0;
		temperatureGradientAtCenter[j].y = 0;
		temperatureGradAtCenterFactor[j] = 0;	
	}

	straumur.x = straumur.y = 0;
}


Ocean::Ocean(const Ocean& O)
{
	printf( "SNAP!  You called Ocean's copy constructor.\n" );
	exit(0);
}

Ocean::Ocean()
{
	RunOnThread = 0;
	fishRecordList = NULL;
	numberOfFish = 0;
        int NumberOfFishToRep = 0;
	sizeOfFish = 32;
	fish = (Fish*)malloc( sizeof(Fish)*sizeOfFish );
}

Ocean::~Ocean()
{
	free(fishRecordList);
	free(fish);
}

void Ocean::add(const FishRecord& R, bool inIsGhost)
{
	Fish F;
	F.mimmicRecord(R);
	if( inIsGhost )
		addAsGhost(F);
	else
		add(F);
}

void Ocean::add(const Fish& F)
{
	reallocate();
	fish[numberOfFish] = F;
	fish[numberOfFish].myOcean = this;
	numberOfFish++;        
}

void Ocean::addAsGhost(const Fish& F)
{
	reallocate();

	fish[numberOfFish] = F;
	fish[numberOfFish].isGhost = true;
	fish[numberOfFish].myOcean = this;
	numberOfFish++;
}

//  Here sizeOfFish is the size of the memory block used for storing the fish
//  It gets doubled when not big enough and halved when huge
void Ocean::reallocate()
{
	if( numberOfFish<16 )
		return; //(if it's that small don't even bother.)

	if( numberOfFish<sizeOfFish/2 )
	{
		sizeOfFish/=2;
		fish = (Fish*)realloc( fish, sizeof(Fish)*sizeOfFish );
	}
	else if( numberOfFish >= sizeOfFish )
	{
		sizeOfFish*=2;
		fish = (Fish*)realloc( fish, sizeof(Fish)*sizeOfFish );
	}
}

void Ocean::print() const
{
	for(int i=0; i<numberOfFish; i++)
	{
		Fish& F(fish[i]);
		printf( "%3.d ", i );
		F.print();
	}
	printf( "\n\n" );
}

void Ocean::printGhostFish() const
{
  for(int i=0; i<numberOfFish; i++)
    {
      Fish& F(fish[i]);
      if( F.isGhost )
	{
	  printf( "%3.d ", i );
	  F.print();
	}
    }
  printf( "\n\n" );
}

void Ocean::clear()
{
	numberOfFish = 0;
	sizeOfFish = 32;
	free(fish);
	fish = (Fish*)malloc(sizeof(Fish)*sizeOfFish); //heh
}

void Ocean::setSize(int inNumberOfFish)
{
	clear();
	Fish F;
	for (int i=0; i < inNumberOfFish; i++)
		add(F); /*	every call to add copies F, so this will actually add numberOfFish brand new fish to the ocean
					each of which is identical to F*/
}

void Ocean::setDirection(double angle)
{
	for(int i=0; i<numberOfFish; i++)
		{
			fish[i].setDirection(angle);
		}
}

Ocean::Ocean(int inNumberOfFish)
{
	RunOnThread = 0;
	fishRecordList = NULL;
	sizeOfFish = 32;
	fish = (Fish*)malloc(sizeof(Fish)*sizeOfFish);
	setSize(inNumberOfFish);
}

void Ocean::zero()
{
	for(int i=0; i<numberOfFish; i++)
		fish[i].zero();
}

int compare(const void * a, const void * b)
{
	double t = ((Fish*)a) -> x - ((Fish*)b) -> x;
	if(t<0)
		return -1;
	else if(t==0)
		return 0;
	else
		return 1;
}

void Ocean::sortByx()
{
	qsort(fish, numberOfFish, sizeof(Fish), compare);
}

int Ocean::binarySearchxCoordinate(double M)
{
	int lower = 0;
	int upper = numberOfFish-1;
	while(upper-lower>1)
	{
		int middle = (upper+lower)/2;
		if( fish[middle].x <= M)
			lower = middle;
		else
			upper = middle;
	}
	if(fish[lower].x < M)
		return upper;
	else
		return lower;
}

void Ocean::initializeFish(double displacementX, double displacementY, double amplitudeX, double amplitudeY, double speedLowerBound, double speedUpperBound, double inSelfWeight)
{
	for(int i=0; i<numberOfFish; i++)
	{
		fish[i].randomSpeed(speedLowerBound, speedUpperBound);
		fish[i].selfWeight = inSelfWeight;
	}
}

/*	Assigns ID numbers to each fish in its fish array */
void Ocean::initializeFishID(int currentFishID)
{
	for(int i=0; i<numberOfFish; i++)
	{
		fish[i].ID = currentFishID+i;
		//printf("fish[%d].ID = %d \n", i, fish[i].ID);
	}
}

/*	This allows us to reset the speed midway through a simulation */
void Ocean::resetSpeed(double inSpeed)
{
	for(int i=0; i<numberOfFish; i++)
	{
		fish[i].resetSpeed(inSpeed);
	}
}

void Ocean::setSpeed(double inSpeed)
{
	for(int i=0; i<numberOfFish; i++)
		{
			fish[i].setSpeed(inSpeed);
		}
}

/* Takes in an angle, amplitude of noise, and a bounding box and calls the function F.resetDirection() for fish within the bounding box */
void Ocean::resetDirection(double theta, double noiseAmp, double xMin, double xMax, double yMin, double yMax)
{
	for(int i=0; i<numberOfFish; i++)
	{
		double x = fish[i].x;
		double y = fish[i].y;
		if(x >= xMin && x <= xMax && y >= yMin && y <= yMax)
			fish[i].resetDirection(theta, noiseAmp);
	}
}

void Ocean::copyVelocities()
{
	for (int i=0; i < numberOfFish; i++)
	{
		fish[i].copyVelocity();
	}
}

void Ocean::interact()
{
	sortByx();
	copyVelocities();
	for(int i=0; i<numberOfFish; i++)
	{
            if (fish[i].age >= AGESELFMOVE)
                {
		Fish& current(fish[i]);  /* This does the same exact thing as Fish& current=fish[i];*/
		current.updateVelocity();
                }
	}
}

void Ocean::addNoiseToDirectionAngle(double noiseAmplitude)
{
	for(int i=0; i<numberOfFish; i++)
	{
		Fish& current(fish[i]);  /* This does the same exact thing as Fish& current=fish[i]; */
		current.addNoiseToDirectionAngle(noiseAmplitude);
	}
}

void Ocean::move(double timestep)
{
	for(int i=0; i<numberOfFish; i++)
		fish[i].move(timestep);
}

/* This function computes the sum of the directional headings over all fish in the ocean*/
Vector Ocean::sumOfDirectionalHeadings()
{
	Vector sum;
	sum.x = 0.0;
	sum.y = 0.0;
	for(int i=0; i<numberOfFish; i++)
	{
		sum.x += fish[i].cosPhi;
		sum.y += fish[i].sinPhi;
	}

	return sum;
}

/* Computes the sum of the x- and y-coordinates for all fish in the ocean */
Vector Ocean::sumOfCoordinates()
{
	Vector sum;
	sum.x = 0.0;
	sum.y = 0.0;
	for(int i=0; i<numberOfFish; i++)
	{
		sum.x += fish[i].x;
		sum.y += fish[i].y;
	}
	return sum;
}

void Ocean::findNearestGridPoint( int gridSizeX, int gridSizeY )
{
	for(int i=0; i<numberOfFish; i++)
		fish[i].findNearestGridPoint( gridSizeX, gridSizeY );

}

void Ocean::emptyLandOfFish()
{
	for(int i=0; i<numberOfFish; i++)
		fish[i].emptyLandOfFish();
}

void Ocean::initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary)
{
	for(int i=0; i<numberOfFish; i++)
	{
		fish[i].initializeFishPositions_Rectangle( leftBoundary, rightBoundary, upperBoundary, lowerBoundary);
	}
}

void Ocean::initializeFishPositions(Vector center1, double majAx1, double minAx1, Vector center2, double majAx2, double minAx2,
					Vector center3, double majAx3, double minAx3, Vector center4, double majAx4, double minAx4,
					Vector center5, double majAx5, double minAx5, Vector center6, double majAx6, double minAx6,
					Vector center7, double majAx7, double minAx7, Vector center8, double majAx8, double minAx8)
{
	for(int i=0; i<numberOfFish; i++)
	{
		fish[i].initializeFishPositions(center1, majAx1, minAx1, center2, majAx2, minAx2,
						center3,  majAx3, minAx3, center4, majAx4, minAx4,
						center5,	majAx5, minAx5, center6, majAx6, minAx6,
						center7,  majAx7, minAx7, center8, majAx8, minAx8);
	}
}

//	Checks if fish is in bounds of the ocean
bool Ocean::inBounds(const Fish& F) const
{
	return( F.x<right && F.x>=left && F.y<top && F.y>=bottom );
}

//	These are the boundaries of the ocean itself
void Ocean::setBounds(double inLeft, double inRight, double inBottom, double inTop)
{
	left = inLeft;
	top = inTop;
	right = inRight;
	bottom = inBottom;
}

//	Should take care of fish outside the current ocean and get rid of them
void Ocean::sendFishOut()
{
	for(int i=0; i<numberOfFish; i++)
	{
		Fish& F(fish[i]);
		if(!inBounds(F))
		{
			F.removeMe = true;
		}
	}

	removeMarkedFish();
}

void Ocean::removeMarkedFish()
{

	/*	okay, the way this works: j and i are indices in the fish array.  i iterates directly through the array, fish i is copied into position j unless it's
		marked for removal (unless F.removeMe is true) in which case it doesn't get copied into position j and j doesn't increment*/
	int j=0;

	for(int i=0; i<numberOfFish; i++)
	{
		Fish& F(fish[i]);
		if(F.removeMe)
			continue;

		if(j<i)
			fish[j] = F;

		j++;
	}

	numberOfFish = j;
	reallocate();
}

void Ocean::removeGhostFish()
{
	/* okay, the way this works: j and i are indices in the fish array.  
         * i iterates directly through the array, fish i is copied into position 
         * j unless it's marked for removal (unless F.removeMe is true) in which
         * case it doesn't get copied into position j and j doesn't increment*/
	int j=0;
        int k = 0;

	for(int i=0; i<numberOfFish; i++)
	{
		Fish& F(fish[i]);
		if(F.isGhost)
			continue;

		if(j<i)
			fish[j] = F;
		j++;
	}

	numberOfFish = j;
	reallocate();
}

/*	This assigns a number to each region around the current ocean.  Hopefully this picture explains the numbering:
		765
		804   <-- the 0 is the current ocean and the others are the oceans around the current ocean.
		123
*/
int Ocean::escapeRegion(const Fish& F) const
{
	if(F.isGhost)
		return 9; // an experimental hack to automatically separate ghosts from the list

	int table[] = {1,2,3,8,0,4,7,6,5};

	int r = (right-tempPadding<F.x)-(F.x<left+tempPadding) + 3*((top-tempPadding<F.y)-(F.y<bottom+tempPadding)); //wa-CHA!!!
	/*	r is a number from -4 to 4 corresponding to the region that the fish is in.
		To get the ordering described above, add 4 and then lookup the new number in the table.*/
	return table[r+4];
}

int Ocean::binarySearchEscapeRegion(int M)
{
	int lower = -1;
	int upper = numberOfFish;

	while(upper-lower>1)
	{
		int middle = (upper+lower)/2;
		if(fish[middle].region <= M)
			lower = middle;
		else
			upper = middle;
	}

	return upper;
}

/*	function for use by Ocean::sortByEscapeRegion in call to qsort*/
int compareEscapeRegion(const void * a, const void * b)
{
	int ai = ((Fish*)a)->region,
		bi = ((Fish*)b)->region;

	return ai-bi;
}

void Ocean::sortByEscapeRegion(double inPadding)
{
	tempPadding = inPadding;

	for(int i=0; i<numberOfFish; i++) /*First memoize the escapeRegion*/
		fish[i].region = fish[i].escapeRegion();

	/*	then call the sort.  The sort uses the memoized escapeRegion, i.e. the member variable region.*/
	qsort(fish, numberOfFish, sizeof(Fish), compareEscapeRegion);
}

void Ocean::makeDensityField(double** densityField, int sizeOfDensityFieldX, int sizeOfDensityFieldY )
{
	for(int i=0; i<numberOfFish; i++)
	{
		Fish &F(fish[i]);

		int x = (int)(sizeOfDensityFieldX * (F.x - left)/((double)(right-left)));
		int y = (int)(sizeOfDensityFieldY * (F.y - bottom)/((double)(top-bottom)));
		if(x<0) x=0;
		if(y<0) y=0;
		if(x>sizeOfDensityFieldX-1) x = sizeOfDensityFieldX-1;
		if(y>sizeOfDensityFieldY-1) y = sizeOfDensityFieldY-1;

		densityField[x][y] += 1.0;
	}
}

/*	randomFish initializes fish and makes all the fish have a random position and speed (according to the argument) and stuff*/
void Ocean::randomFish(double speedLowerBound, double speedUpperBound, double inSelfWeight)
{
	initializeFish( 0.5*(left+right), 0.5*(bottom+top), 0.5*(right-left), 0.5*(top-bottom), speedLowerBound, speedUpperBound, inSelfWeight );
}


/*	alignedFish initializes fish and makes all the fish's directional headings equal */
void Ocean::alignedFish(double speedLowerBound, double speedUpperBound, double cosPhi, double sinPhi, double inSelfWeight)
{
	initializeFish( 0.5*(left+right), 0.5*(bottom+top), 0.5*(right-left), 0.5*(top-bottom), speedLowerBound, speedUpperBound, inSelfWeight );
	for( int i=0; i<numberOfFish; i++)
	{
		fish[i].cosPhi = cosPhi;
		fish[i].sinPhi = sinPhi;
	}
}

void Ocean::initializeMaturityLevels(int initialMaturityLevel)
{
	for( int i=0; i<numberOfFish; i++)
		fish[i].initializeMaturityLevels(initialMaturityLevel);
}


/*	Initiates the variables for the DEB model */

void Ocean::initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime)
{
	for( int i=0; i<numberOfFish; i++)
		fish[i].initializeDEB(l0, e0, uR0, er0, initialMaturityLevel, initialDaysSinceMature, currentTime);
}

void Ocean::solveDEB(double h, double startTime, double endTime)		// This function solves the coupled DEB ODEs using 4th order Runge Kutta.
{
	for( int i=0; i<numberOfFish; i++)
		fish[i].solveDEB(h, startTime, endTime);
}


int Ocean::constructList()
{
	if(fishRecordList==NULL)
		fishRecordList = (FishRecord*)malloc(sizeof(FishRecord)*numberOfFish);
	else
		fishRecordList = (FishRecord*)realloc(fishRecordList, sizeof(FishRecord)*numberOfFish);

	for(int i=0; i<numberOfFish; i++)
	  fish[i].makeRecord(fishRecordList[i]);

	return numberOfFish;
}


void Ocean::destroyList()
{
	if(fishRecordList)
	{
		free(fishRecordList);
		fishRecordList = NULL;
	}
}

// Ocean
// --------------------------------------------------//
/* Here are all the functions that I added to        *
 * the class Ocean                                   */
// --------------------------------------------------//

//Jorge: Set the reproductive status of the fish!
void Ocean::setReproduction(int k, bool Reproduced)
	{
	for(int i=0; i<numberOfFish; i++)
		fish[i].setReproduction(k,true);
	}

int Ocean::numberOfFishToRepToday(int day)
        {
        int RepToday = 0; 
        for (int z=0; z<numberOfFish; z++)
                {
                int temp = fish[z].dayOfReproduction;
                if (temp == day & fish[z].age >= AGEMATURE)
                        {
                        RepToday++;
                        }
                }
        return RepToday;
        }


void Ocean::getOld(int TotalNumberFish, float toKillToday, float tJuv)
{
    //srand (time(NULL)); // new seed for every set of fish
    nJuveniles = 0;
    nAdults    = 0;
    nOld       = 0;
   
    if (numberOfFish > 0)
    {
        for (int z=0; z<numberOfFish; z++)
            {
                fish[z].age ++;
                //printf("Age: %i\n", fish[z].age);
                if      (fish[z].age < AGEMATURE) 
                {
                    nJuveniles++;
                }
                else if (fish[z].age >= AGEMATURE and fish[z].age < AGEOLD) 
                {
                    nAdults++;
                }
                else    //(fish[z].age >= AGEOLD) 
                {
                    nOld++;
                }
            }

       // printf("In Get Old [%i, %i, %i]: ", nJuveniles, nAdults, nOld);
        #if USEMORTALITY
                NaturalMortality(1, TotalNumberFish, tJuv, nJuveniles, nAdults, nOld, toKillToday);
        #endif
    }                    
}

void Ocean::NaturalMortality(int type, int TotalNumberFish, float tJuv, int nJuveniles, int nAdults, int nOld, float toKillToday)
{
        //Natural mortality is evaluated
    if (type == 0)
    {

        float juvDie =  (exp(Z_JUVENILES)); //Number of fish that will die at this day
        if (nJuveniles == 0) {juvDie = 0;}
        float aduDie =  (exp(Z_ADULTS    * (nAdults + nOld))) / nAdults ;
        if (nAdults == 0) {aduDie = 0;}
        float oldDie =  (exp(Z_OLD       * (nAdults + nOld))) / nOld ;
        if (nOld == 0) {oldDie = 0;}
        //printf("[%.3f, %.3f, %.3f] \n", juvDie, aduDie, oldDie);
        
        //This is how the sample of fish is selected to be removed
        
        if (juvDie > 0 or aduDie>0 or oldDie > 0)
        {
            for (int z = 0; z<numberOfFish; z++)
            {
               fish[z].pp = (float) rand()/RAND_MAX;
               
               if (fish[z].age < AGEMATURE and fish[z].pp < juvDie)
                {
                        fish[z].removeMe = true;
                    }
               if (fish[z].age >= AGEMATURE and fish[z].age < AGEOLD and fish[z].pp < aduDie)
                {
                        fish[z].removeMe = true;
                    }
               if (fish[z].age >= AGEOLD and fish[z].pp < oldDie)
                {
                        fish[z].removeMe = true;
                    }

                }
            //printf("%lu, %i, %lu, %g, %s\n", juvDie, fish[z].age, aduDie, exp(Z_JUVENILES * (nAdults + nOld + nJuveniles)), (fish[z].removeMe)?"true":"false");
        }
        
    }
    else if (type == 1)
    {
        //printf("%2.2f \n", toKillToday);
        float oCC = OCEANCC;
        float ccFactor = (nAdults + nOld)/oCC;
        float lCC = OCEANLARVAECC;
        float lccFactor = 0;
        float ppToDie_J;
        if (toKillToday <= 0) 
                {ppToDie_J = 0;}
        else {ppToDie_J = toKillToday;} 
        float ppToDie_A = Z_ADULTS     * ccFactor; 
        float ppToDie_O = Z_OLD        * ccFactor; 
        int o = 0;
        for (int z = 0; z < numberOfFish; z++)
        {
            Fish &F(fish[z]);
            float temp2 = (rand()%1000);
            float temp = temp2 /  1000;
            
            if (F.age < AGEMATURE and temp < ppToDie_J and ppToDie_J != 0) //and temp != 0) 
            {
                F.removeMe = true;
            }
            if (F.age >= AGEMATURE and F.age < AGEOLD and temp < ppToDie_A and ppToDie_A != 0 and temp != 0)
            {
                F.removeMe = true;
            }
            if (F.age >= AGEOLD and temp < ppToDie_O and ppToDie_O != 0 and temp != 0)
            {
                F.removeMe = true;
            }
        }
    }
    removeMarkedFish();
}

int Ocean::isBorn(int iniID, int day, int nEggs)
{ 
    //Ocean& O;
    int nOfFemales;
    day = JulDay(day);
    nOfFemales = numberOfFishToRepToday(day);
    int newFish = nOfFemales * nEggs;
    int trackID = iniID;

    if (newFish > 0)
    {
        int CountOfFemales = 1;  
        for (int i=0; i<numberOfFish; i++)
        {   
            
            if (fish[i].age >= AGEMATURE and fish[i].dayOfReproduction == day)
            {
                double x = 0, y=0;
                double lon = 0, lat=0;
                for (int z=0; z<nEggs; z++)
                {
                    x = rand()%50;
                    y = rand()%50;
                    lon = fish[i].x + x/100;
                    lat = fish[i].y + y/100;
                    
                    Fish F;
                    F = fish[i];
                    F.age = 0;
                    F.ID = trackID;
                    F.setDirection(rand()%360);
                    float temp = rand()%500;
                    F.speed = temp/10000;
                    F.setPosition(lon, lat);
                    F.selfWeight = .1;
                    if (TYPEOFREP == 0)
                    {
                        F.dayOfReproduction = probOfReproduction(TYPEOFREP);
                    }
                    else if (TYPEOFREP == 1)
                    {
                        F.dayOfReproduction = probOfReproduction(TYPEOFREP, MU, STD);
                    }
                    trackID++;
                    add(F);
                }
                CountOfFemales++;
            }
            if (CountOfFemales == nOfFemales) {break;}
        }
    }
    return trackID;
    
}


/* Start of World class */
void World::copyDensityFieldIntoPicture(int offsetX, int offsetY, int thread)
{
	double totalFish = totalFishFromAllProcessors;

	for(int j=0; j<sizeOfDensityFieldY; j++)
	for(int i=0; i<sizeOfDensityFieldX; i++)
	{
		double red = densityField[i][j]/((double)totalFish/(double)(sizeOfImageX*sizeOfImageY));
		double blue = (1-red) - 0.1*((double)thread);
		double green = densityField[i][j]/((double)totalFish/(double)(sizeOfImageX*sizeOfImageY))/2-1;

		picture.setPixel(i+offsetX, j+offsetY,  red, green, blue);
	}
}

void World::sendTotalNumberOfFish()
{
	int t = howManyFish;

	#if MPI_ON
		MPI_Send(&t, 1, MPI_DOUBLE, 0, 20000, MPI_COMM_WORLD);
	#endif
}


void World::getTotalFishFromAllProcessors()
{
	totalFishFromAllProcessors = 0;

	for(int p=1; p<numberOfProcessors; p++)
	{
		int t;
		#if MPI_ON
			MPI_Status status;
			MPI_Recv(&t, 1, MPI_DOUBLE, p, 20000, MPI_COMM_WORLD, &status);
		#endif

		totalFishFromAllProcessors += t;
	}
}

void World::sendDensities()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
        Ocean& O = oceans[x][y];
		if(O.RunOnThread == rank)
		{
		  clearDensity();
			O.makeDensityField(densityField, sizeOfDensityFieldX, sizeOfDensityFieldY);

			#if MPI_ON
				MPI_Send(densityFieldData, sizeOfDensityFieldX*sizeOfDensityFieldY, MPI_DOUBLE, 0, 10000+x+sizeX*y, MPI_COMM_WORLD);
			#endif
		}
	}
}

void World::receiveDensityForOcean(const Ocean& O, int x, int y)
{
	#if MPI_ON
		MPI_Status status;
		MPI_Recv(densityFieldData, sizeOfDensityFieldX*sizeOfDensityFieldY, MPI_DOUBLE, O.RunOnThread, 10000+x+sizeX*y, MPI_COMM_WORLD, &status);
	#endif
}

void World::getPicture()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O = oceans[x][y];

		clearDensity();

		if(O.RunOnThread == 0)
		{
			/*If the ocean is mine just draw its own density field*/
			O.makeDensityField(densityField, sizeOfDensityFieldX, sizeOfDensityFieldY);
		}
		else
		{
			/*If the ocean is somebody elses, then we receive the message from it*/
			receiveDensityForOcean(O, x,y);
		}
		copyDensityFieldIntoPicture(x*sizeOfDensityFieldX, y*sizeOfDensityFieldY, O.RunOnThread);
	}
}


void World::clearDensity()
{
        bzero(densityFieldData, sizeof(double)*sizeOfDensityFieldX*sizeOfDensityFieldY);
}


/* FUNCTIONS SPECIFICALLY BUILT FOR THE MEASURES PERTAINING TO SWARMING */

/* Computes the average directional heading over all the fish in the world */
void World::findAverageDirectionalHeading()
{
	Vector temporaryDirection;
	Vector sumDirectionalHeading;
	sumDirectionalHeading.x = 0.0;
	sumDirectionalHeading.y = 0.0;
	double norm;

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		temporaryDirection = oceans[x][y].sumOfDirectionalHeadings();
		sumDirectionalHeading.x += temporaryDirection.x;
		sumDirectionalHeading.y += temporaryDirection.y;
	}

	temporaryDirection.x = (1/double(howManyFish))*sumDirectionalHeading.x;
	temporaryDirection.y = (1/double(howManyFish))*sumDirectionalHeading.y;

	printf( "temporaryDirection is (%f, %f) \n", temporaryDirection.x, temporaryDirection.y);

	norm = sqrt( temporaryDirection.x*temporaryDirection.x + temporaryDirection.y*temporaryDirection.y );
	if( norm > 0.000001)
	{
		averageDirection.x = (1/norm)*temporaryDirection.x;
		averageDirection.y = (1/norm)*temporaryDirection.y;
	}
	else
		averageDirection = temporaryDirection;
}

void World::findAlignmentFactor()
{
	findAverageDirectionalHeading();
	alignmentFactor = 0.0;
	averageSpeed = 0.0;

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		for(int i=0; i<oceans[x][y].numberOfFish; i++)
		{
			double cosPhi = oceans[x][y].fish[i].cosPhi;
			double sinPhi = oceans[x][y].fish[i].sinPhi;

			if( oceans[x][y].fish[i].speed < 0)
			printf( "oceans[%d][%d].fish[%d].speed = %f!!!!!!!!!!!!!!!!!!!!!!!!! \n", x, y, i, oceans[x][y].fish[i].speed);

			averageSpeed += oceans[x][y].fish[i].speed;

			alignmentFactor += fabs(cosPhi*averageDirection.y - sinPhi*averageDirection.x);
//			printf( "Adding %f to the alignment Factor; cosPhi = %f, sinPhi = %f, averageDirection.x = %f, averageDirection.y = %f \n", alignmentFactor, cosPhi, sinPhi, averageDirection.x, averageDirection.y);
				/* This is the area of the parallelogram between the unit direction vector and the unit average direction vector */
		}
	}

	alignmentFactor = alignmentFactor/double(howManyFish);  // Taking into account how many fish there are in the world

	alignmentFactor = 1/alignmentFactor;  // Trying to make it so that the alignment is very low when the fish point different ways

	averageSpeed = averageSpeed/double(howManyFish);

	printf( "alignmentFactor = %f \n averageSpeed = %f \n",alignmentFactor, averageSpeed);

}

/* Finds the center of mass for the world */
void World::findCenterOfMass()
{
	Vector localSumOfCoordinates;
	Vector sumOfCoordinates;
	sumOfCoordinates.x = 0;
	sumOfCoordinates.y = 0;

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		localSumOfCoordinates = oceans[x][y].sumOfCoordinates();
		sumOfCoordinates.x += localSumOfCoordinates.x;
		sumOfCoordinates.y += localSumOfCoordinates.y;
	}

	centerOfMass.x = (1/double(howManyFish))*(sumOfCoordinates.x);
	centerOfMass.y = (1/double(howManyFish))*(sumOfCoordinates.y);
}

void World::initializeMaturityLevels(int initialMaturityLevel)
{
	for(int y=0; y<sizeY; y++)
		for(int x=0; x<sizeX; x++)
			oceans[x][y].initializeMaturityLevels(initialMaturityLevel);
}



/*	DEB functions */

/*	Initiates the variables for the DEB model */
void World::initializeDEB(double l0, double e0, double uR0, double er0, int initialMaturityLevel, double initialDaysSinceMature, double currentTime)
{
	for(int y=0; y<sizeY; y++)
		for(int x=0; x<sizeX; x++)
			oceans[x][y].initializeDEB(l0, e0, uR0, er0, initialMaturityLevel, initialDaysSinceMature, currentTime);
}

void World::solveDEB(double h, double startTime, double endTime)		// This function solves the coupled DEB ODEs using 4th order Runge Kutta.
{
	for(int y=0; y<sizeY; y++)
		for(int x=0; x<sizeX; x++)
		{
			//printf("solveDEB for oceans[%d][%d], number of fish in ocean = %d \n", x,y,oceans[x][y].numberOfFish);
			oceans[x][y].solveDEB(h,startTime,endTime);
		}
}


World::World(const World& W)
{
	printf( "WOAH!  You shouldn't call World's copy constructor.  Trust me.\n" );
	exit(0);
}


World::World(int inSizeX, int inSizeY)
{
        printf("In World(int inSizeX, int inSizeY)\n");
	headNode = 0;
	rank = 0;

	isTorus = false;
	padding = DEFAULT_PADDING;
 	densityField = NULL;
	processors = NULL;
	allocateOceans(inSizeX,inSizeY);
}


World::World()
{
        printf("In World()\n");
	headNode = 0;

	rank = 0;

	isTorus = false;
	padding = DEFAULT_PADDING;
	densityField = NULL;
	processors = NULL;
	oceans = NULL;
}

void World::allocatePrefTempRanges()
{
	tooCold = new double[NUMMATURITYLEVELS];
	tooHot = new double[NUMMATURITYLEVELS];
}

void World::allocateOceans(int inSizeX, int inSizeY)
{
	sizeX = inSizeX;
	sizeY = inSizeY;

	oceans = new Ocean*[sizeX];
	for(int x=0; x<sizeX; x++)
	{
		oceans[x] = new Ocean[sizeY];
	}

	for(int x=0; x<sizeX; x++)
	for(int y=0; y<sizeY; y++)
		oceans[x][y].myWorld = this;
}

void World::allocateProcessors(int inSize)
{
	processors = new Processor[inSize];
	numberOfProcessors = inSize;
}

World::~World()
{
        deallocateDensityField();
	deallocateOceans();
	deallocateProcessors();
}


void World::deallocatePrefTempRanges()
{
	delete[] tooCold;
	delete[] tooHot;
}



void World::deallocateOceans()
{
	if( oceans )
	{
		for(int x=0; x<sizeX; x++)
			delete[] oceans[x];

		delete[] oceans;

		oceans = NULL;
	}
}

void World::deallocateProcessors()
{
	if( processors )
	{
		delete[] processors;
		processors = NULL;
	}
}


void World::initSize(int inSizeX, int inSizeY, double inOWidth, double inOHeight)
{
	if( oceans != NULL )
	{
		printf( "World::initSize should not be called on a World constructed with the default constructor.\n" );
		exit(0);
	}

	allocateOceans(inSizeX, inSizeY);

	oWidth = inOWidth;
	oHeight = inOHeight;

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);
                if (x == 0 && y ==0)
                {
                    O.setBounds( x*oWidth, (x+2)*oWidth,  y*oHeight, (y+2)*oHeight );
                } else 
                {
		O.setBounds( x*oWidth, (x+1)*oWidth,  y*oHeight, (y+1)*oHeight );
                }
	}
}


void World::initMultiThread()
{
        #if MPI_ON
                /*INSERT MPI CODE HERE*/
                MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcessors );
                MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        #else
                numberOfProcessors = 1;
                rank = 0;
        #endif
        printf("\n This is the rank %d of a total of %d processors \n", rank, numberOfProcessors );
        
        allocateProcessors(numberOfProcessors);
        for(int p =0; p<numberOfProcessors; p++)
                {
                processors[p].thread = p; ///LATER WE WILL WRITE CODE THAT ASSUMES THAT PROCESSOR[P].THREAD = P
                processors[p].load = 0;
                }
}


void World::setDirection(double angle)
        {
	for(int y=0; y<sizeY; y++)
				for(int x=0; x<sizeX; x++)
					oceans[x][y].setDirection(270);
        }

int World::sendFishTo(int targetX, int targetY, Fish* fishArray, int placeholderStart, int placeholderNextStart, bool inIsGhost)
{
	int fishSent = 0;

	if( isTorus )
	{
		double offsetX=0, offsetY=0;

		if( targetX >= sizeX )
		{
			targetX = 0;
			offsetX = -oWidth*sizeX;
		}

		if( targetX < 0 )
		{
			targetX = sizeX-1;
			offsetX = oWidth*sizeX;
		}

		if( targetY >= sizeY )
		{
			targetY = 0;
			offsetY = -oHeight*sizeY;
		}

		if( targetY < 0 )
		{
			targetY = sizeY-1;
			offsetY = oHeight*sizeY;
		}

		Ocean& targetOcean(oceans[targetX][targetY]);

		for(int i = placeholderStart; i < placeholderNextStart; i++)
		{
			Fish F(fishArray[i]);
			F.x += offsetX;
			F.y += offsetY;
			if( inIsGhost )
				targetOcean.addAsGhost(F);
			else
				targetOcean.add(F);

			fishSent++;
		}
	}
	else
	{
		if( targetX >= sizeX || targetX < 0 || targetY >= sizeY || targetY < 0 )
			return 0;

		fishSent = sendFishTo(oceans[targetX][targetY], fishArray, placeholderStart, placeholderNextStart, inIsGhost);
	}


	return fishSent;
}



int World::sendFishTo(Ocean& targetOcean, Fish* fishArray, int placeholderStart, int placeholderNextStart, bool inIsGhost)
{
	if( inIsGhost )
	{
		for(int i = placeholderStart; i < placeholderNextStart; i++)
			targetOcean.addAsGhost(fishArray[i]);
	}
	else
	{
		for(int i = placeholderStart; i < placeholderNextStart; i++)
			targetOcean.add(fishArray[i]);
	}

	return placeholderNextStart-placeholderStart;
}

void World::sendFishOut()
{
	//Loop through oceans:
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].sendFishOut();
	}
}

void World::resetInteractionCounters()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].interactionCounter = 0;
	}
}



void World::transportFish()
{
	int relativeXPosition[] = {0,-1,0,1,1,1,0,-1,-1};	// Defining the relative x-position of the various neighboring oceans
	int relativeYPosition[] = {0,-1,-1,-1,0,1,1,1,0};	// Defining the relative y-position of the various neighboring oceans

	//	Iterate through the oceans:
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);

		O.sortByEscapeRegion(0); //Here 0 is the padding


		int placeholderStart = 0;
		int placeholderNextStart = O.binarySearchEscapeRegion(0);

		int fishsent = 0;

		for(int i=1; i<9; i++)
		{
			placeholderStart = placeholderNextStart;
			placeholderNextStart = O.binarySearchEscapeRegion(i);

			if (placeholderStart == placeholderNextStart)
				continue;

			int targetX = x + relativeXPosition[i];
			int targetY = y + relativeYPosition[i];

			fishsent += sendFishTo(targetX, targetY, O.fish, placeholderStart, placeholderNextStart);
		}

		O.numberOfFish = O.binarySearchEscapeRegion(0);
	}
	totalNumberOfFish();
}




void World::communicateGhostFish()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);

		O.sortByEscapeRegion(padding);

		int placeholders[10];

		for(int i=0; i<10; i++)
			placeholders[i] = O.binarySearchEscapeRegion(i);

		/*	7 6 5
			8 0 4
			1 2 3 */

		sendFishTo(x-1, y-1,	O.fish,		placeholders[1-1], placeholders[1], true );
		sendFishTo(x, y-1,		O.fish,		placeholders[1-1], placeholders[3], true );
		sendFishTo(x+1, y-1,	O.fish,		placeholders[3-1], placeholders[3], true );
		sendFishTo(x+1, y,		O.fish,		placeholders[3-1], placeholders[5], true );
		sendFishTo(x+1, y+1,	O.fish,		placeholders[5-1], placeholders[5], true );
		sendFishTo(x, y+1,		O.fish,		placeholders[5-1], placeholders[7], true );
		sendFishTo(x-1, y+1,	O.fish,		placeholders[7-1], placeholders[7], true );

		sendFishTo(x-1, y,		O.fish,		placeholders[7-1], placeholders[8], true );
		sendFishTo(x-1, y,		O.fish,		placeholders[1-1], placeholders[1], true );	/*It's beautiful*/
	}
}


void World::removeGhostFish()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].removeGhostFish();
	}
}


void World::interact()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
	  oceans[x][y].interactionCounter = 0;
	  oceans[x][y].interact();
	}
}

/*	This allows us to reset the speed midway through a simulation */
void World::resetSpeed(double inSpeed)
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
	  oceans[x][y].resetSpeed(inSpeed);
	}
}

void World::setSpeed(double inSpeed)
{

	for(int y=0; y<sizeY; y++)
		for(int x=0; x<sizeX; x++)
		{
		  oceans[x][y].setSpeed(.21);
		}

}

/*	This calls O.resetDirection, which takes in an angle, amplitude of noise, and a bounding box and calls the function F.resetDirection() for fish within the bounding box */
void World::resetDirection(double theta, double noiseAmp, double xMin, double xMax, double yMin, double yMax)
{
printf("In World::resetDirection, which resets the direction of fish inside a box\n");
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
	  oceans[x][y].resetDirection(theta, noiseAmp, xMin, xMax, yMin, yMax);
	}
}


void World::addNoiseToDirectionAngle(double noiseAmplitude)
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].addNoiseToDirectionAngle(noiseAmplitude);
	}
}


void World::allocatePicture(double scaling)
{  
	allocatePicture( (int)(sizeX*scaling*oWidth + 0.5 ), (int)(sizeY*scaling*oHeight + 0.5 ) );
}

void World::allocatePicture(int width, int height)
{
    sizeOfImageX = width;
    sizeOfImageY = height;
    picture.resize(24, sizeOfImageX, sizeOfImageY);
}



void World::allocateDensityField(double scaling)
{
    allocateDensityField( (int)(sizeX*scaling*oWidth + 0.5 ), (int)(sizeY*scaling*oHeight + 0.5 ) );
}

void World::allocateDensityField(int width, int height)
{
	sizeOfDensityFieldX = width;
	sizeOfDensityFieldY = height;

	densityFieldData = new double[sizeOfDensityFieldX*sizeOfDensityFieldY];
	densityField = new double*[sizeOfDensityFieldX];

	for(int j=0; j < sizeOfDensityFieldX; j++)
	{
		densityField[j] = densityFieldData +j*sizeOfDensityFieldY;
	}
}

void World::deallocateDensityField()
{
    if(densityField)
        delete densityField;
    if(densityFieldData)
        delete densityFieldData;
}


void World::initPictureSerial()
{
    allocateDensityField(PICTURESCALING);
    allocatePicture(PICTURESCALING);
}


void World::initPictureHeadNode()
{
	allocateDensityField(PICTURESIZEX, PICTURESIZEY);
	allocatePicture(sizeX*PICTURESIZEX, sizeY*PICTURESIZEY);
}

void World::initPictureClientNode()
        {
      	allocateDensityField(100, 100);
        }


void World::densityPointillism(Fish& fish)
{
	int x = (int)(fish.x*(double)sizeOfDensityFieldX/(oWidth*(double)sizeX));
	int y = (int)(fish.y*(double)sizeOfDensityFieldY/(oHeight*(double)sizeY));

	densityField[x][y] += 1;
}


void World::drawPicture(const char* filename, int gridSizeX, int gridSizeY)
{
	for (int i=0; i < sizeOfDensityFieldX; i++)
	for (int j=0; j < sizeOfDensityFieldY; j++)
	{
		densityField[i][j] = 0.0;
	}

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
                {
                Ocean& O = oceans[x][y];
                for(int i=0; i < O.numberOfFish; i++)
                        densityPointillism(O.fish[i]);
                }

	int totalFish = howManyFish;

	for (int i=0; i < sizeOfImageX; i++)
//	for (int j=1; j < sizeOfImageY; j++)
	for (int j=0; j < sizeOfImageY; j++)
	{
		double interpolationConstant=0.0;
		double red = 0.0;
		double green = 0.0;
		double blue = 1.0;
		if(densityField[i][j]>0.0)
		{
			interpolationConstant=((double)densityField[i][j]-1.0)/DOUBLE_DETERMINING_RED_DENSITY;	//creates a range to show the fish density between yellow and red, for instance 25 fish make it fully red :D
			blue = 0;
			red = 1;
			green = 1-interpolationConstant;
		}
		picture.setPixel(i, j, red, green, blue);
	}

#if USEGRID
	for (int i=0; i < sizeOfImageX; i++)
	for (int j=0; j < sizeOfImageY; j++)
	{
		double worldWidth = sizeX;
		double worldHeight = sizeY;

		IntegerVector nearestestGridPoint;

		nearestestGridPoint.x = (int)((double)i/(double)PICTURESCALING + 0.5);
		nearestestGridPoint.y = (int)((double)j/(double)PICTURESCALING + 0.5);

		if( nearestestGridPoint.x < 0 ) nearestestGridPoint.x=0;
		if( nearestestGridPoint.x > gridSizeX-1 ) nearestestGridPoint.x=gridSizeX-1;

		if( nearestestGridPoint.y < 0 ) nearestestGridPoint.y=0;
		if( nearestestGridPoint.y > gridSizeY-1 ) nearestestGridPoint.y=gridSizeY-1;

		double temp = grid[nearestestGridPoint.x][nearestestGridPoint.y].temperature;
#if DRAWTEMP
		if( fabs(temp - TOOHOT[0]) <= .1 )
//			picture.setPixel(i, j-1, 1, 0, 0);
			picture.setPixel(i, j, 1, 0, 0);
		else if( fabs(temp - TOOCOLD[0]) <= .1 )
//			picture.setPixel(i, j-1, 0, 1, 0);
			picture.setPixel(i, j, 0, 1, 0);
#endif

#if DRAW_LAND_OUTLINE
printf("Beware of Greenland \n");
		if( temp > 900 )
		{
			if(  grid[nearestestGridPoint.x+1][nearestestGridPoint.y].temperature < 900 ||   grid[nearestestGridPoint.x-1][nearestestGridPoint.y].temperature < 900
			||   grid[nearestestGridPoint.x][nearestestGridPoint.y+1].temperature < 900 ||   grid[nearestestGridPoint.x][nearestestGridPoint.y-1].temperature < 900	
			|| grid[nearestestGridPoint.x+1][nearestestGridPoint.y+1].temperature < 900 || grid[nearestestGridPoint.x-1][nearestestGridPoint.y-1].temperature < 900
			|| grid[nearestestGridPoint.x+1][nearestestGridPoint.y-1].temperature < 900 || grid[nearestestGridPoint.x-1][nearestestGridPoint.y+1].temperature < 900 )
//				picture.setPixel(i, j-1, 0, 0, 0);
				picture.setPixel(i, j, 0, 0, 0);
		}
#endif

#if DRAW_LAND_SOLIDLY

		if( temp > 900 )
			picture.setPixel(i, j, 0, 0, 0);
#endif
	}
#endif

	picture.write(filename);

}

void World::drawPictureAdults(const char* filename, int gridSizeX, int gridSizeY)
{
	for (int i=0; i < sizeOfDensityFieldX; i++)
	for (int j=0; j < sizeOfDensityFieldY; j++)
	{
		densityField[i][j] = 0.0;
	}

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
                {
                Ocean& O = oceans[x][y];
                for(int i=0; i < O.numberOfFish; i++)
                    if (O.fish[i].age >=365)
                        densityPointillism(O.fish[i]);
                }

	int totalFish = howManyFish;

	for (int i=0; i < sizeOfImageX; i++)
//	for (int j=1; j < sizeOfImageY; j++)
	for (int j=0; j < sizeOfImageY; j++)
	{
		double interpolationConstant=0.0;
		double red = 0.0;
		double green = 0.0;
		double blue = 1.0;
		if(densityField[i][j]>0.0)
		{
			interpolationConstant=((double)densityField[i][j]-1.0)/DOUBLE_DETERMINING_RED_DENSITY;	//creates a range to show the fish density between yellow and red, for instance 25 fish make it fully red :D
			blue = 0;
			red = 1;
			green = 1-interpolationConstant;
		}
//		picture.setPixel(i, j-1, red, green, blue);
		picture.setPixel(i, j, red, green, blue);
	}
	picture.write(filename);

}

void World::writePictureToFile(const char* filename)
{
	picture.write(filename);
}

void World::initFishFromFile(const char* filename)
{
    FILE* fp = fopen(filename, "r");
	if(!fp)
	{
        printf("fish file %s didn't open\n",filename);
		exit(0);
	}
        fseek(fp, 0, SEEK_END);
	int fileSize = ftell(fp);
	printf("fileSize %d \n",fileSize);
	char* content = (char*)malloc(sizeof(char)*fileSize);
	fseek(fp, 0, SEEK_SET);
	fread(content, sizeof(char), fileSize, fp);

	fclose(fp);



        int numberOfEntries = fileSize;// number of characters in file, used to iterate through content array
        int howManyEntriesSoFar = 0;	// this keeps track of which entry we are getting data for
	int i=0;	// i keeps track of which character we are reading from the file
	int j=0;	// j keeps track of what character of an entry we are writing
        int k=0;

	char entry[2000];
	double x[60000];
	double y[60000];
	double angle[60000];
	double entryDouble;

	while(howManyEntriesSoFar <= numberOfEntries && i <= fileSize)
                {
		if(content[i] == ' ' && content[i+1] == ' '&& content[i+2] == ' '&& content[i+3] == ' ')
                        {
			entry[j] = 0;
			entryDouble = atof(entry);
                        y[k]=entryDouble;
                        i=i+3;
			j = 0;
                        }
		else if(content[i] == ' ' && content[i+1] == ' ' && content[i+2] == ' '&& content[i+3] != ' ')
                        {
			entry[j] = 0;
                        entryDouble = atof(entry);
			x[k]=entryDouble;
			i = i+2;
			j = 0;
        		}
        	else if(content[i] == '\n')//angle
				{
					entry[j] = 0;
					entryDouble = atof(entry);
					angle[k]=entryDouble;
					for(int countX=0;countX<sizeX; countX++) //loop determines what ocean the fish is in and places it in the correct array
					for(int countY=0;countY<sizeY;countY++)
					if(x[k]<oceans[countX][countY].right && x[k]>=oceans[countX][countY].left && y[k]<oceans[countX][countY].top && y[k]>=oceans[countX][countY].bottom)
        					{
						//creates a fish, sets its position to the coordinates that were read then places that fish into its ocean
						Fish F;
						Ocean& O = oceans[countX][countY];
						F.setDirection(angle[k]);
						F.setPosition(x[k],y[k]);
                                                F.setReproduction(0, true);
						O.add(F);
                                                }
					howManyEntriesSoFar++;
					k++;
					i = i+2;
					j = 0;
				}
		else
		{
			entry[j]=content[i];
			j++;
			i++;
		}

            }
	
	free(content);
	printf("Imported: %s",filename);
    printf("\n");

	return;

};
void World::randomFish(int inNumberOfFish, double speedLowerBound, double speedUpperBound, double inSelfWeight)
{
printf("In World::randomFish, which randomizes both speed and direction\n");
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		//oceans[x][y].setSize(inNumberOfFish); removed, now called by initFishFishFromFile()
		oceans[x][y].randomFish(speedLowerBound, speedUpperBound, inSelfWeight);
	}
}


/*  An initialization function puts randomly-located random-speed, directionally aligned fish into each ocean*/
void World::alignedFish(int inNumberOfFish, double speedLowerBound, double speedUpperBound, double cosPhi, double sinPhi, double inSelfWeight)
{
printf("In World::alignedFish, which randomizes both speed BUT fixes the direction\n");
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].alignedFish(speedLowerBound, speedUpperBound, cosPhi, sinPhi, inSelfWeight);
	}
}


void World::allocateGrid(int x, int y)
{
	//	Here x and y are the number of grid points,
	//	hence the range of the grid will be from 0 to x-1...
	sizeGridsX = x;
	sizeGridsY = y;

	grid = new GridPoint*[sizeGridsX];
	for (int i = 0; i < sizeGridsX; i++)
	{
		grid[i] = new GridPoint[sizeGridsY];
		//Allocating the stuff in a grid point is done in the GridPoint constructor
	}

	zeroGrid(x,y);

	allocatePrefTempRanges();
}


void World::deallocateGrid()
{
	for (int i=0; i<sizeGridsX; i++)
	{	//NB the stuff stored in a grid point is deleted in the GridPoint destructor ~GridPoint()
		//Delete the collumns of grid points:
		delete[] grid[i];
	}
	//At last, we delete the row of pointers to the grid point columns:
	delete[] grid;

	deallocatePrefTempRanges();
}

void World::zeroGrid(int xDimension, int yDimension)
{
	//	Here x and y are the number of grid points,
	//	hence the range of the grid will be from 0 to xDim-1...
	for (int i=0; i<xDimension; i++)
		for (int j = 0; j<yDimension; j++)
		{
			grid[i][j].zero();
		}
}


void World::findNearestGridPoint( int gridSizeX, int gridSizeY )
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].findNearestGridPoint( gridSizeX, gridSizeY );
	}
}


/* See #defined preferred temperature ranges in fish.h */
void World::initializePrefTempRanges(float TOOCOLD_IMMATURE, float TOOHOT_IMMATURE, 
        float TOOCOLD_MATURE, float TOOHOT_MATURE)
        {
        printf("  Temp Preferences TC: %f TH: %f \n", TOOCOLD_IMMATURE, TOOHOT_MATURE);
	tooCold[IMMATURE]=TOOCOLD_IMMATURE;
	tooHot[IMMATURE]=TOOHOT_IMMATURE;

	tooCold[MATURE]=TOOCOLD_MATURE;
	tooHot[MATURE]=TOOHOT_MATURE;
        }


/*   given a filename, and dimensions xDimension and yDimension, this function translates a list of numbers from a file into an n by m matrix */
void World::translateTemperatureToGrid(const char* filename, int xDimension, int yDimension)
{
/*	file lists temperatures by row left to right, from bottom to top */
	FILE* fp = fopen(filename, "r");
	if(!fp)
	{
        printf("temperature file %s didn't open\n",filename);
		exit(0);
	}
	fseek(fp, 0, SEEK_END);
	int fileSize = ftell(fp);
	char* content = (char*)malloc(sizeof(char)*fileSize);
	fseek(fp, 0, SEEK_SET);
	fread(content, sizeof(char), fileSize, fp);

	fclose(fp);
	int numberOfEntries = xDimension*yDimension;
	int howManyEntriesSoFar = 0;	// this keeps track of which entry we are getting data for
	int i=0;	// i keeps track of which character we are reading from the file
	int j=0;	// j keeps track of what character of an entry we are writing
	char entry[1000];
	double entryDouble;
	while( howManyEntriesSoFar <= numberOfEntries && i <= fileSize)
	{
		if (content[i] == ' ')
			i++;

		if(content[i] == '\n')
		{
			entry[j] = 0;
			int yIndex = (int) (howManyEntriesSoFar/xDimension);
			int xIndex = howManyEntriesSoFar - yIndex*xDimension;
			entryDouble = atof(entry);						// could also write: entryDouble = strtod(entry, NULL);
			grid[xIndex][yIndex].temperature = entryDouble;
			howManyEntriesSoFar++;
			i++;
			j = 0;
		}
		else
		{
			entry[j]=content[i];
			j++;
			i++;
		}
	}
	free(content);

	resetTemperatureGradientFlags();
	return;
}

void World::resetTemperatureGradientFlags()
{
	for(int i=0; i<sizeGridsX; i++)
		for(int j=0; j<sizeGridsY; j++)
		{
			grid[i][j].gradientFlag=false;
			grid[i][j].tempCorrectionFlag=false; //Used for DEB
			for(int k=0; k<4; k++)
				grid[i][j].flagQuadrant[k]=false;
		}
}

void World::translateStraumurToGrid(const char* filename, int xDimension, int yDimension)
{
/*	file lists straumur (the currents at the gridpoints) as coordinate pairs (x,y) by row left to right, from bottom to top */
	FILE* fp = fopen(filename, "r");
	if(!fp)
	{
		printf("straumur file didn't open\n");
		exit(0);
	}
	fseek(fp, 0, SEEK_END);
	int fileSize = ftell(fp);
	char* content = (char*)malloc(sizeof(char)*fileSize);
	fseek(fp, 0, SEEK_SET);
	fread(content, sizeof(char), fileSize, fp);
	fclose(fp);
	int numberOfEntries = xDimension*yDimension;
	int howManyEntriesSoFar = 0;	// this keeps track of which entry we are getting data for
	int i=0;	// i keeps track of which character we are reading from the file
	int j=0;	// j keeps track of what character of an entry we are writing
	char entry[1000];
	double entryDouble;
	while( howManyEntriesSoFar < numberOfEntries && i <= fileSize)
	{
		if(content[i] == '\n')
		{
			entry[j] = 0;

			int yIndex = (int) (howManyEntriesSoFar/xDimension);
			int xIndex = howManyEntriesSoFar - yIndex*xDimension;
			entryDouble = atof(entry);					
			grid[xIndex][yIndex].straumur.y = entryDouble;
			howManyEntriesSoFar++;
			i=i+2;
			j = 0;
		}
		else if(content[i] == ' ' && content[i+1] == ' ')
		{
			entry[j] = 0;

			int yIndex = (int) (howManyEntriesSoFar/xDimension);
			int xIndex = howManyEntriesSoFar - yIndex*xDimension;
			entryDouble = atof(entry);	
			grid[xIndex][yIndex].straumur.x = entryDouble;
			i = i+2;
			j = 0;
		}
		else
		{
			entry[j]=content[i];
			j++;
			i++;
		}
	}
	free(content);
	printf("Imported: %s",filename);
    printf("\n");

	return;
}//World::translateStraumurToGrid



/*	the following two functions get the info from GridPoints that's needed to update the fish's directions */
double World::getTemperature(IntegerVector nearestGridPoint)
{
	return grid[nearestGridPoint.x][nearestGridPoint.y].temperature;
}


//Uses quadrants in order to calculate gradients in triangles
Vector World::getTemperatureGradient(double normTempSens, IntegerVector nearestGridPoint, int quadrant, int fishMaturityLevel)
{
	int x = nearestGridPoint.x;
	int y = nearestGridPoint.y;
		
	//If not, then we have to calculate the temperature gradients and norms:
	GridPoint& G = grid[x][y];
	double norm;
	G.flagQuadrant[quadrant] =true; // So that we won't calculate it again unless it gets updated
	
	//Also the Arrhenius temperature correction (might be unnecessary)
	#if DEB_ON
		double currentTemp = G.temperature;
		if( currentTemp > 100)		// This is to correct for the land points being extremely hot
		{	currentTemp = DEB_T1;	// The current temperature should be in Kelvin
		}
		else
		{	currentTemp += 273.15; //Haha, I forgot this line at first					
		}
		G.temperatureCorrection = exp( (DEB_TA/DEB_T1) - (DEB_TA/currentTemp) );
		G.tempCorrected_nu    = G.temperatureCorrection*DEB_NU;
		G.tempCorrected_kJ    = G.temperatureCorrection*DEB_KJ;
		G.tempCorrected_gamma = G.temperatureCorrection*DEB_GAMMA;
	#endif

	if( quadrant == 0)
	{
		for (int k=0; k<NUMMATURITYLEVELS; k++)
		{
			
			G.temperatureGradient[quadrant][k].x = R( x+1, y, tooCold[k], tooHot[k]) - R(x, y, tooCold[k], tooHot[k]);
			G.temperatureGradient[quadrant][k].y = R( x, y+1, tooCold[k], tooHot[k]) - R(x, y, tooCold[k], tooHot[k]);
			
			norm = G.temperatureGradient[quadrant][k].x*G.temperatureGradient[quadrant][k].x+G.temperatureGradient[quadrant][k].y*G.temperatureGradient[quadrant][k].y;
			norm = sqrt(norm);

			// Normalize only if norm>normTempSens
			if(norm > normTempSens)
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/norm;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/norm;
				// Now, save the temperature gradient factor in quadrant 0 as 1:
				G.temperatureGradFactor[quadrant][k] = 1;
			}
			else
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/normTempSens;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/normTempSens;
				// Now, save the temperature gradient factor in quadrant 0:
				G.temperatureGradFactor[quadrant][k] = norm/normTempSens;
			}
		}
		//Finally, return the gradient according to the fish's maturity level:
		return G.temperatureGradient[quadrant][fishMaturityLevel];
	}
	if( quadrant == 1)
	{
		for (int k=0; k<NUMMATURITYLEVELS; k++)
		{
			//Have to return G.temperatureGradient[quadrant][fishMaturityLevel] 
		
			G.temperatureGradient[quadrant][k].x = R( x, y, tooCold[k], tooHot[k] ) - R(x-1, y, tooCold[k], tooHot[k]);
			G.temperatureGradient[quadrant][k].y=R( x, y+1, tooCold[k], tooHot[k] ) - R( x , y, tooCold[k], tooHot[k]);
			
			norm = G.temperatureGradient[quadrant][k].x*G.temperatureGradient[quadrant][k].x+G.temperatureGradient[quadrant][k].y*G.temperatureGradient[quadrant][k].y;
			norm = sqrt(norm);		
			
			// Normalize only if norm>normTempSens
			if(norm > 1)
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/norm;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/norm;
				// Now, save the temperature gradient factor in quadrant 1 as 1:
				G.temperatureGradFactor[quadrant][k] = 1; 
			}
			else
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/normTempSens;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/normTempSens;
				// Now, save the temperature gradient factor in quadrant 1:
				G.temperatureGradFactor[quadrant][k] = norm/normTempSens;		
			}
		}
		return G.temperatureGradient[quadrant][fishMaturityLevel];
	}
	if( quadrant == 2)
	{	
		for (int k=0; k<NUMMATURITYLEVELS; k++)
		{
			//Have to return G.temperatureGradient[quadrant][fishMaturityLevel] 
		
			G.temperatureGradient[quadrant][k].x = R( x, y, tooCold[k], tooHot[k] ) - R(x-1, y, tooCold[k], tooHot[k]);
			G.temperatureGradient[quadrant][k].y = R( x, y, tooCold[k], tooHot[k] ) - R(x, y-1, tooCold[k], tooHot[k]);
			
			norm = G.temperatureGradient[quadrant][k].x*G.temperatureGradient[quadrant][k].x+G.temperatureGradient[quadrant][k].y*G.temperatureGradient[quadrant][k].y;
			norm = sqrt(norm);

			// Normalize only if norm>normTempSens
			if(norm > 1)
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/norm;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/norm;
				// Now, save the temperature gradient factor in quadrant 2 as 1:
				G.temperatureGradFactor[quadrant][k] = 1; 
			}
			else
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/normTempSens;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/normTempSens;
				// Now, save the temperature gradient factor in quadrant  2:
				G.temperatureGradFactor[quadrant][k] = norm/normTempSens;	
			}
		}
		return G.temperatureGradient[quadrant][fishMaturityLevel];
	}
	if( quadrant == 3)
	{	
		for (int k=0; k<NUMMATURITYLEVELS; k++)
		{
			//Have to return G.temperatureGradient[quadrant][fishMaturityLevel] 
		
			G.temperatureGradient[quadrant][k].x = R( x+1, y, tooCold[k], tooHot[k] ) - R(x,  y,  tooCold[k], tooHot[k]);
			G.temperatureGradient[quadrant][k].y = R( x ,  y, tooCold[k], tooHot[k] ) - R(x, y-1, tooCold[k], tooHot[k]);
			
			norm = G.temperatureGradient[quadrant][k].x*G.temperatureGradient[quadrant][k].x+G.temperatureGradient[quadrant][k].y*G.temperatureGradient[quadrant][k].y;
			norm = sqrt(norm);

			// Normalize only if norm>normTempSens
			if(norm > 1)
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/norm;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/norm;
				// Now, save the temperature gradient factor in quadrant 3 as 1:
				G.temperatureGradFactor[quadrant][k] = 1; 
			}
			else
			{
				G.temperatureGradient[quadrant][k].x=G.temperatureGradient[quadrant][k].x/normTempSens;
				G.temperatureGradient[quadrant][k].y=G.temperatureGradient[quadrant][k].y/normTempSens;
				// Now, save the temperature gradient factor in quadrant 3:
				G.temperatureGradFactor[quadrant][k] = norm/normTempSens;
			}
		}
		return G.temperatureGradient[quadrant][fishMaturityLevel];
	}

	printf("Oops! Quadrant = %d doesn't exist! \n", quadrant);
} //World::getTemperatureGradient


//Sven's simplified version. The gradient is calculated as an average from the nodes adjacent to the nearest one
Vector World::getTemperatureGradient2(double normTempSens, IntegerVector nearestGridPoint, int fishMaturityLevel)
{
	int x = nearestGridPoint.x;
	int y = nearestGridPoint.y;
	
	//NB should depend on the maturity levels

	if(grid[x][y].gradientFlag)
	{
		return grid[x][y].temperatureGradientAtCenter[fishMaturityLevel];
	}
	else
	{
		//If not, then we have to calculate the temperature gradients and norms:
		GridPoint& G = grid[x][y];
		double norm;
		G.gradientFlag = true; //So that we won't calculate it again
		
		//Also the Arrhenius temperature correction:
		#if DEB_ON
			if (!G.tempCorrectionFlag) //Needs to be calculated for this temperature
			{
				G.tempCorrectionFlag = true; //Could also be calculated if temperature is within the preferred temp range
				double currentTemp = G.temperature;
				if( currentTemp > 25.0)		// This is to correct for the land points being extremely hot
				{	currentTemp = DEB_T1;	// The current temperature should be in Kelvin
					//printf("Oops, fish is close to land and it messes with the DEB\n");
				}
				else
				{
					currentTemp += 273.15; //Haha, at first I forgot this line
				}
				G.temperatureCorrection = exp( (DEB_TA/DEB_T1) - (DEB_TA/currentTemp) );
				G.tempCorrected_nu    = G.temperatureCorrection*DEB_NU;
				G.tempCorrected_kJ    = G.temperatureCorrection*DEB_KJ;
				G.tempCorrected_gamma = G.temperatureCorrection*DEB_GAMMA;
			}
		#endif
		
		//Now loop through the maturity levels:
		for (int k=0; k<NUMMATURITYLEVELS; k++)
		{
			//Have to return G.temperatureGradientAtCenter[fishMaturityLevel] 
			if( x == 0)
				G.temperatureGradientAtCenter[k].x = R(x+1, y, tooCold[k], tooHot[k]) - R(x, y, tooCold[k], tooHot[k]);
			else if( x == GRIDSIZEX-1)
				G.temperatureGradientAtCenter[k].x = R(x, y, tooCold[k], tooHot[k]) - R(x-1, y, tooCold[k], tooHot[k]);
			else
				G.temperatureGradientAtCenter[k].x = 0.5*(R(x+1,y, tooCold[k], tooHot[k]) - R(x-1, y, tooCold[k], tooHot[k]));

			if( y == 0)
				G.temperatureGradientAtCenter[k].y = R(x, y+1, tooCold[k], tooHot[k]) - R(x, y, tooCold[k], tooHot[k]);
			else if( y == GRIDSIZEY-1)
				G.temperatureGradientAtCenter[k].y = R(x, y, tooCold[k], tooHot[k]) - R(x, y-1, tooCold[k], tooHot[k]);
			else
				G.temperatureGradientAtCenter[k].y = 0.5*(R(x, y+1, tooCold[k], tooHot[k]) - R(x, y-1, tooCold[k], tooHot[k]));

			norm = G.temperatureGradientAtCenter[k].x*G.temperatureGradientAtCenter[k].x+G.temperatureGradientAtCenter[k].y*G.temperatureGradientAtCenter[k].y;
			norm = sqrt(norm);
			//printf("normReturnValue = %f\n", normReturnValue);
			
			// Normalize and return 1 as temperature gradient factor
			if(norm > normTempSens)
			{
			//	printf("dividing by normReturnValue!\n");
				G.temperatureGradientAtCenter[k].x=G.temperatureGradientAtCenter[k].x/norm;
				G.temperatureGradientAtCenter[k].y=G.temperatureGradientAtCenter[k].y/norm;
				grid[x][y].temperatureGradAtCenterFactor[k] = 1;
			}
			else //Don't normalize but divide by normTempSens.
			{	G.temperatureGradientAtCenter[k].x=G.temperatureGradientAtCenter[k].x/normTempSens;
				G.temperatureGradientAtCenter[k].y=G.temperatureGradientAtCenter[k].y/normTempSens;
				grid[x][y].temperatureGradAtCenterFactor[k] = norm/normTempSens;
			}
		
		
		return G.temperatureGradientAtCenter[fishMaturityLevel];
		
	}
	printf("We should never get to here--if so, getTemperatureGradient2() is broken!\n");
	exit (1);
} //World::getTemperatureGradient2



/*	this function computes the temperature part of the comfort function for grid[x][y] given upper and
	lower comfort bounds on the temperature */
// Modified  temperature function to make the gradual ///
double World::R(int x, int y, double T1, double T2)
{
	double T = grid[x][y].temperature;
	double R;
	if( T <= T1)
		R = -0.0001 * (T-T1)*(T-T1) ;  //
                //R = - (T - T1)*(T - T1)*(T - T1)*(T - T1);
	else if( T1 < T && T < T2)
		R = 0;
	else
		R = -0.0001 * (T - T2)*(T - T2); //
                //R = - ( T - T2)*(T - T2);
	return R;
}



void World::move()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].move();
	}
}


void World::emptyLandOfFish()
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].emptyLandOfFish();
		oceans[x][y].removeMarkedFish();
	}
}


void World::initializeFishPositions_Rectangle(double leftBoundary, double rightBoundary, double upperBoundary, double lowerBoundary)
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].initializeFishPositions_Rectangle( leftBoundary, rightBoundary, upperBoundary, lowerBoundary);
		oceans[x][y].removeMarkedFish();
	}
}

void World::initializeFishPositions(Vector center1, double majAx1, double minAx1, Vector center2, double majAx2, double minAx2,
					Vector center3, double majAx3, double minAx3, Vector center4, double majAx4, double minAx4,
					Vector center5, double majAx5, double minAx5, Vector center6, double majAx6, double minAx6,
					Vector center7, double majAx7, double minAx7, Vector center8, double majAx8, double minAx8)
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		oceans[x][y].initializeFishPositions(center1, majAx1, minAx1, center2, majAx2, minAx2,
						center3,  majAx3, minAx3, center4, majAx4, minAx4,
						center5,	majAx5, minAx5, center6, majAx6, minAx6,
						center7,  majAx7, minAx7, center8, majAx8, minAx8);
		oceans[x][y].removeMarkedFish();
	}
}



void World::translateByStraumur(double timestep, int gridSizeX, int gridSizeY)
{
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);
		for(int i=0; i<O.numberOfFish; i++)
		{
			Fish& F(O.fish[i]);
			F.findNearestGridPoint(gridSizeX, gridSizeY);
			F.x += timestep*grid[F.nearestGridPoint.x][F.nearestGridPoint.y].straumur.x;
			F.y += timestep*grid[F.nearestGridPoint.x][F.nearestGridPoint.y].straumur.y;
		}
	}
}


/* This function gets called from gridtest. It in turn calls getTemperatureGradient */
void World::factorInTemperatureGradients(double normTempSens, int gridSizeX, int gridSizeY)
{
	//Initially, tempWeight is set to the default value which is increased if deemed necessary.
	double tempWeight = DEFAULT_TEMPERATUREWEIGHT;
	Vector temperatureDirection;// to keep track of the information that getTemperatureGradient returns
	Vector newDirection;	// to not overwrite cosPhi and sinPhi until later
	double normNewDirection;			// the norm of the newDirection vector, to make (cosPhi, sinPhi) unit length

	//Iterates through the oceans:
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);
		//Iterates through the fish of the oceans:
		for(int i=0; i<O.numberOfFish; i++)
		{
			Fish& F(O.fish[i]);
			F.findNearestGridPoint(gridSizeX, gridSizeY);
			//printf("(F.x, F.y)=(%f, %f) and (F.nearestGridPoint.x, F.nearestGridPoint.y) = (%d, %d) \n", F.x, F.y, F.nearestGridPoint.x, F.nearestGridPoint.y);
		
			GridPoint& G(grid[F.nearestGridPoint.x][F.nearestGridPoint.y]
			if( G.temperature > tooCold[F.maturityLevel] && G.temperature < tooHot[F.maturityLevel])
			{
				#if DEB_ON
				if (!G.tempCorrectionFlag)
				{
					G.tempCorrectionFlag = true; //Also calculated when temperature gradients are found
					double currentTemp = G.temperature;
					if( currentTemp > 100.0)		// This is to correct for the land points being extremely hot
					{	currentTemp = DEB_T1;	// The current temperature should be in Kelvin
						//printf("Oops, fish is close to land and it messes with the DEB\n");
					}
					else
					{
						currentTemp += 273.15; //Haha, I forgot this line at first
					}
					G.temperatureCorrection = exp( (DEB_TA/DEB_T1) - (DEB_TA/currentTemp) );
					G.tempCorrected_nu    = G.temperatureCorrection*DEB_NU;
					G.tempCorrected_kJ    = G.temperatureCorrection*DEB_KJ;
					G.tempCorrected_gamma = G.temperatureCorrection*DEB_GAMMA;
				}
				#endif

				//The temperature therefore doesn't affect the directional heading of the fish and so we...
				continue; 
			}

                        #if TRIANGULAR_GRADIENT
			if( F.x > F.nearestGridPoint.x )
			{
				if( F.y > F.nearestGridPoint.y )
					F.quadrant = 0;
				else
					F.quadrant = 3;
			}
			else
			{
				if( F.y > F.nearestGridPoint.y)
					F.quadrant = 1;
				else
				{
					F.quadrant = 2;
				}
			}
			//Now, check whether outside boundary:
			if(F.nearestGridPoint.x==0 && (F.quadrant==1 || F.quadrant==2))
				{
					printf("Oops! F.nearestGridPointx==0 and F.quadrant = %d doesn't exist! \n", F.quadrant);
					printf("F.x=%f and F.y=%f, F.removeMe %d\n",F.x,F.y,F.removeMe);
					printf("F.myOcean->left: %f,top: %f,right: %f,bottom: %f\n",F.myOcean->left,F.myOcean->top,F.myOcean->right,F.myOcean->bottom);
					exit(1);
				}

				if(F.nearestGridPoint.y==0 && (F.quadrant==2 || F.quadrant==3))
				{
					printf("Oops! F.nearestGridPoint.y==0 and F.quadrant = %d doesn't exist! \n", F.quadrant);
					printf("F.x=%f and F.y=%f, F.removeMe %d\n",F.x,F.y,F.removeMe);
					printf("F.myOcean->left: %f,top: %f,right: %f,bottom: %f\n",F.myOcean->left,F.myOcean->top,F.myOcean->right,F.myOcean->bottom);
					exit(1);
				}

				if(F.nearestGridPoint.x==GRIDSIZEX-1 && (F.quadrant==0 || F.quadrant==3))
				{
					printf("Oops! F.nearestGridPoint.x==%d and F.quadrant = %d doesn't exist! \n", GRIDSIZEX-1, F.quadrant);
					printf("F.x=%f and F.y=%f, F.removeMe %d\n",F.x,F.y,F.removeMe);
					printf("F.myOcean->left: %f,top: %f,right: %f,bottom: %f\n",F.myOcean->left,F.myOcean->top,F.myOcean->right,F.myOcean->bottom);
					exit(1);
				}

				if(F.nearestGridPoint.y==GRIDSIZEY-1 && (F.quadrant==0 || F.quadrant==1))
				{
					printf("Oops! F.nearestGridPoint.y==%d and quadrant = %d doesn't exist! \n", GRIDSIZEY-1, F.quadrant);
					printf("F.x=%f and F.y=%f, F.removeMe %d\n",F.x,F.y,F.removeMe);
					printf("F.myOcean->left: %f,top: %f,right: %f,bottom: %f\n",F.myOcean->left,F.myOcean->top,F.myOcean->right,F.myOcean->bottom);
					exit(1);
				}


			// getTemperatureGradient2 is Sven's simplified idea, getTemperatureGradient is the one with the quadrants and triangles
			temperatureDirection = getTemperatureGradient(normTempSens, F.nearestGridPoint, F.quadrant, F.maturityLevel);
			
			/*if (!G.gradientFlag)
			{
			printf("Sheize, something isn't working in factorInTemperatueGradients \n");
			exit(0);
			}*/
			
			
			//Older version with LAND_TEMPWEIGHT:
			//Check whether fish is close to land or not:
			/*if(G.temperature > 100)
			{
				newDirection.x = (1-LAND_TEMPWEIGHT)*F.cosPhi + LAND_TEMPWEIGHT*temperatureDirection.x;
				newDirection.y = (1-LAND_TEMPWEIGHT)*F.sinPhi + LAND_TEMPWEIGHT*temperatureDirection.y;
			}*/
			
			bool loop = true;
			while(loop)
			{
				if(tempWeight>1)
				{
					//Let the new direction only depend on the temperature gradient:
					newDirection.x = temperatureDirection.x;
					newDirection.y = temperatureDirection.y;
					loop = false;
				}
				else
				{
			
					/* NB: temperatureDirection is normalized only if the norm is greater than normTempSens
					Note that the factor is imbedded in the temperatureGradientDirection (!)*/
					newDirection.x = (1-tempWeight*G.temperatureGradFactor[F.quadrant][F.maturityLevel])*F.cosPhi + tempWeight*temperatureDirection.x;
					newDirection.y = (1-tempWeight*G.temperatureGradFactor[F.quadrant][F.maturityLevel])*F.sinPhi + tempWeight*temperatureDirection.y;
			
					if(grid[(int)(F.x + TIMESTEP*F.speed*newDirection.x + 0.5)][(int)(F.y + TIMESTEP*F.speed*newDirection.y + 0.5)].temperature > 100.0)
					{//Things are too hot at the next location. Increase the tempWeigt:
						tempWeight += DEFAULT_TEMPERATUREWEIGHT;
					}
					else //we have found the new direction which leads to a place not near land.
	 					loop = false;
				}
			}

                        #endif

                        #if 1-TRIANGULAR_GRADIENT
			// getTemperatureGradient2 is Sven's simplified idea, getTemperatureGradient is the one with the quadrants and triangles
			temperatureDirection = getTemperatureGradient2(normTempSens, F.nearestGridPoint , F.maturityLevel);

			bool landFlag = false;
			bool loop = true;
			while(loop)
			{
				if(tempWeight>1)
				{
					//Let the new direction only depend on the temperature gradient:
					newDirection.x = temperatureDirection.x;
					newDirection.y = temperatureDirection.y;
                                        landFlag = true;
					loop = false;
				}
				else
				{
                                    landFlag = false;
					/* NB: temperatureDirection is normalized only if the norm is greater than normTempSens
					Note that the factor is imbedded in the temperatureGradientDirection (!)*/
					newDirection.x = (1-tempWeight*G.temperatureGradAtCenterFactor[F.maturityLevel])*F.cosPhi + tempWeight*temperatureDirection.x;
					newDirection.y = (1-tempWeight*G.temperatureGradAtCenterFactor[F.maturityLevel])*F.sinPhi + tempWeight*temperatureDirection.y;
			
					if(grid[(int)(F.x + TIMESTEP*F.speed*newDirection.x + 0.5)][(int)(F.y + TIMESTEP*F.speed*newDirection.y + 0.5)].temperature > 100.0)
					{//Things are too hot at the next location. Increase the tempWeigt:
						tempWeight += DEFAULT_TEMPERATUREWEIGHT;
					}
					else //we have found the new direction which leads to a place not near land.
						loop = false;
				}
			}

#endif

			// Here we have calculated normNewDirection either using triangle gradients or Sven's method
			normNewDirection = sqrt( newDirection.x*newDirection.x + newDirection.y*newDirection.y);
                        if (F.age >= AGESELFMOVE ) //Added here to avoid eggs and larvae react to temperature!
                        {
                            if( normNewDirection < DECISION_TOL )
                            {	//Don't normalize since vector is too small to make any decisions from
                                    F.cosPhi = newDirection.x;
                                    F.sinPhi = newDirection.y;
                                    //if(F.ID==0)
                                    //	printf("not dividing by normNewDirection! normNewDirection = %f \n", normNewDirection);
                            }
                            else
                            {
                                    F.cosPhi = newDirection.x/normNewDirection;
                                    F.sinPhi = newDirection.y/normNewDirection;
                            }
                            if(F.ID== - 5001) 
                                printf("%2.2f, %2.2f  ", newDirection.x, F.cosPhi);
                            //if(F.ID==0)
                            	//printf("At the end of factorInTemperatureGradients(), (cosPhi,sinPhi)=(%f,%f)\n",F.cosPhi,F.sinPhi);
                        }
       
                            float tempCosPhi = F.cosPhi;
                            float tempSinPhi = F.sinPhi;
                            if( normNewDirection < DECISION_TOL )
                            {	//Don't normalize since vector is too small to make any decisions from
                                    F.cosPhi = newDirection.x;
                                    F.sinPhi = newDirection.y;
                                    //if(F.ID==0)
                                    //	printf("not dividing by normNewDirection! normNewDirection = %f \n", normNewDirection);
                            }
                            else
                            {
                                    //if(F.ID==0)
                                    //	printf("dividing by normNewDirection.  normNewDirection = %f \n", normNewDirection);
                                    F.cosPhi = newDirection.x/normNewDirection;
                                    F.sinPhi = newDirection.y/normNewDirection;
                            }
                            if (F.cosPhi == 0) {F.cosPhi = tempCosPhi;}
                            if (F.sinPhi == 0) {F.sinPhi = tempSinPhi;}
                            if(F.ID== - 5001) 
                                printf("%2.2f, %2.2f  ", newDirection.x, F.cosPhi);
                            //if(F.ID==0)
                            	//printf("At the end of factorInTemperatureGradients(), (cosPhi,sinPhi)=(%f,%f)\n",F.cosPhi,F.sinPhi);
                        
                        
                        
		}
	}
} //World::factorInTemperatureGradients



void World::totalNumberOfFish()
{
	int R = 0;
        int Rep = 0;
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		R+=oceans[x][y].numberOfFish;
                //Number of fish that can reproduce
                Rep+= oceans[x][y].numberOfFishToRep;
	}
	howManyFish = R;
        howManyFishToRep = Rep;
}


// Assigns a positive integer to each fish (from being -1) - works in serial for now
void World::initializeFishID()
{
	int currentFishID = 0;
	for(int y=0; y<sizeY; y++)
		for(int x=0; x<sizeX; x++)
		{
			//printf("oceans[%d][%d].numberOfFish = %d: \n", x, y, oceans[x][y].numberOfFish);
			oceans[x][y].initializeFishID(currentFishID);
			currentFishID += oceans[x][y].numberOfFish;
		}

	if(currentFishID != howManyFish)
		printf("Error in assigning fish ID: World::initializeFishID \n");
	else
		printf("IDs successfully assigned. Current id %d \n", currentFishID);
        LastUsedIDFish = currentFishID;
}



void World::printWhoLivesWhere() const
{
	printf("Thread [%d]: Printing Who lives where: \n", rank);
	for (int y=sizeY-1; y >=0 ; y--)
	{
		printf("y = %d | ", y);

		for (int x=0; x < sizeX; x++)
		{
			printf("%d ",  oceans[x][y].RunOnThread );
		}
		printf("\n");
	}
}




void World::printNumbers() const
{
    printf("fish in world: %d\n", howManyFish );
    for(int y=sizeY-1; y>=0; y--)
    {
        for(int x=0; x<sizeX; x++)
        {
            int m = oceans[x][y].numberOfFish;
            printf( m?"%3.0d ":"  %d ", m );
        }

        printf( "\n" );
    }
    printf( "\n\n" );
}


void World::printInteractionCounter() const
{
    printf( "fish in world: %d\n", howManyFish );

    for(int y=sizeY-1; y>=0; y--)
    {
        for(int x=0; x<sizeX; x++)
        {
            int m = oceans[x][y].interactionCounter;
            printf("%d ", m );
        }

        printf( "\n" );
    }
    printf( "\n\n" );
}

int World::totalInteractionCounter() const
{
	int totalInteractions = 0;
	for(int y=sizeY-1; y>=0; y--)
        for(int x=0; x<sizeX; x++)
            totalInteractions += oceans[x][y].interactionCounter;
	return totalInteractions;
}

int World::biggestLoadedProcessor()
{
        int biggestSoFar = 0;

	for(int i=0; i<numberOfProcessors; i++)
	{
	  if(processors[i].load > processors[biggestSoFar].load)
	  biggestSoFar = i;
	}

	return biggestSoFar;
}


Processor World::smallestLoadedNeighbor(Processor& processorHigh)
{
	int neighborX;
	int neighborY;
	int minSoFar=processorHigh.load;
	int leastLoadSoFar = processorHigh.load;
	Processor lowestLoadedNeighborSoFar = processorHigh;

	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O = oceans[x][y];

		if( O.RunOnThread == processorHigh.thread )
		{
			for(int deltaX = -1; deltaX <= 1; deltaX++)
			for(int deltaY = -1; deltaY <= 1; deltaY++)
			{
				if ( !((deltaX + deltaY) == 1 || (deltaX + deltaY) == -1 ) )
					continue;

				neighborX = x+deltaX;
				neighborY = y+deltaY;

				if( (!isTorus) && (neighborX<0 || neighborX>=sizeX || neighborY<0 || neighborY>=sizeY) )
					continue;

				neighborX = (neighborX+sizeX)%sizeX;
				neighborY = (neighborY+sizeY)%sizeY;

				int neighborThread = oceans[neighborX][neighborY].RunOnThread;

				if(neighborThread == processorHigh.thread)
				  continue;

				int neighborLoad = processors[neighborThread].load;

				if ( neighborLoad < leastLoadSoFar)
				  {
					lowestLoadedNeighborSoFar = processors[neighborThread];
					leastLoadSoFar = neighborLoad;
				  }
			}
		}
	}

	//printf("lowestLoadedNeighborSoFar = %d \n",lowestLoadedNeighborSoFar);
	return lowestLoadedNeighborSoFar;
}


bool World::findOptimalOceanForExchange(Processor& processorHigh, Processor& processorLow, int& outX, int& outY)
{
	double badness=0;
	int neighborX;
	int neighborY;
	int neighborProcessLoad;
	int oceanCandidateX = -1;
	int oceanCandidateY = -1;
	double minSoFar = (double) processorHigh.load;
	int leastLoadSoFar = processorHigh.load;
	int myOceanCounter=0;

	Processor lowestLoadedNeighborSoFar = processorHigh;
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O = oceans[x][y];

		if( O.RunOnThread == processorHigh.thread )
		{
		  myOceanCounter++;
			for(int deltaX = -1; deltaX <= 1; deltaX++)
			for(int deltaY = -1; deltaY <= 1; deltaY++)
			{
				if ( !((deltaX + deltaY) == 1 || (deltaX + deltaY) == -1 ) )
					continue;

				neighborX = x+deltaX;
				neighborY = y+deltaY;

				if( (!isTorus) && (neighborX<0 || neighborX>=sizeX || neighborY<0 || neighborY>=sizeY) )
					continue;

				neighborX = (neighborX+sizeX)%sizeX;
				neighborY = (neighborY+sizeY)%sizeY;

				Ocean &neighborOcean(oceans[neighborX][neighborY]);
				if( neighborOcean.RunOnThread != processorLow.thread )
				  continue;

				neighborProcessLoad = processors[neighborOcean.RunOnThread].load;

				printf( "neighbor: %d %d has load: %d\n", neighborX, neighborY, neighborOcean.interactionCounter );

				badness = abs(neighborProcessLoad + O.interactionCounter - (processorHigh.load + neighborProcessLoad)/2);
				  printf("badness = %f \n",badness);
				if (badness < minSoFar)
				{
				        minSoFar = badness;
					oceanCandidateX = x;
					oceanCandidateY = y;
				}
			}
		}
	}


	if ( oceanCandidateX < 0 || oceanCandidateY < 0 || myOceanCounter<=1 )
		return true;

	outX = oceanCandidateX;
	outY = oceanCandidateY;

	return false;
}



void World::sendLoadToHead()
{
#if MPI_ON
	for(int i=0; i<sizeX; i++)
	for(int j=0; j<sizeY; j++)
	{
		Ocean &O(oceans[i][j]);
		if(O.RunOnThread==rank)
		{
		  MPI_Send( &(O.interactionCounter), 1, MPI_INT, headNode, i*sizeY+j, MPI_COMM_WORLD );
		}
	}
#endif
}


void World::getAllOceanLoads()
{
#if MPI_ON

  for(int p=0; p<numberOfProcessors; p++)
    {
      processors[p].load = 0;
    }

	MPI_Status status;

	for(int i=0; i<sizeX; i++)
	for(int j=0; j<sizeY; j++)
	{
	        Ocean &O(oceans[i][j]);
		if(O.RunOnThread!=rank)
		{
			MPI_Recv( &(O.interactionCounter), 1, MPI_INT, O.RunOnThread, i*sizeY+j, MPI_COMM_WORLD, &status );
		}
		processors[O.RunOnThread].load += O.interactionCounter;
	}


#endif
}


int World::updateHeadNode()
{
#if MPI_ON
	sendLoadToHead();

	MPI_Status status;
	MPI_Recv(&exchangeData, sizeof(ExchangeData), MPI_CHAR, headNode, 20001, MPI_COMM_WORLD, &status);

   return 1;
#endif
}


void World::updateMyWorld()
{
#if MPI_ON
  if( exchangeData.noExchange )
    return;

  if( oceans[exchangeData.x][exchangeData.y].RunOnThread == exchangeData.targetThread )
    return;     //this should actually never happen, I think



  if( oceans[exchangeData.x][exchangeData.y].RunOnThread == rank )
    {
      //I'm giving my ocean to another process.
      // rank here is rank if the current processor; i.e. the source of the message
      // thread is the destination of the message
      //void World::MPI_SendFish(int srcI, int srcJ, int rank, int thread, int  flag)

      MPI_SendFish( exchangeData.x, exchangeData.y, rank, exchangeData.targetThread, MIGRATING);
    }
  else if( exchangeData.targetThread == rank )
    {
      //I'm receiving a new ocean from somebody.
      // rank here is rank of the current thread (i.e. the destination of the message
      // thread is source of the message
      //void World::MPI_RecvFish(int targetI, int targetJ, int rank, int thread, int flag)
      MPI_RecvFish(exchangeData.x, exchangeData.y, rank, oceans[exchangeData.x][exchangeData.y].RunOnThread, MIGRATING);
    }

  oceans[exchangeData.x][exchangeData.y].RunOnThread = exchangeData.targetThread;
#endif
}



// A simple OceanMigrating function.
// This will find the heaviest ocean on processorHigh, and choose it
// evacuation from the current to the target processor (i.e. processorLow)
// probably we will need to construct a more sophisticated and intelligent
// method on how to exchange oceans...based on their neighbors and where they live,
// the source processor that ocean will be evacuated from and the dst processor to
// which ocean will migrate to.
void World::constructExchangeData(Processor processorHigh, Processor processorLow){
#if MPI_ON
        int x, y; // x and y are co-ordinate of potential ocean to be evacuated

	// first, find the biggest ocean on processorHigh

	exchangeData.noExchange = findOptimalOceanForExchange(processorHigh, processorLow, x, y);
	exchangeData.x = x;
	exchangeData.y = y;
	exchangeData.targetThread = processorLow.thread;

	return;
#endif
}



int World::runHeadNode()
{
#if MPI_ON
	//receive other processors' information

	//uses mpi to collect everybody's load
        getAllOceanLoads();

	//find biggest loaded processor
	Processor processorHigh;
	int numberOfProcessorWithBiggestLoad = biggestLoadedProcessor();
	processorHigh = processors[numberOfProcessorWithBiggestLoad];

	//find the smallest loaded neighbor to that processor
	Processor processorLow  = smallestLoadedNeighbor(processorHigh);


	//find optimal ocean for exchange and use that to construct the data structure exchangeData
	constructExchangeData(processorHigh, processorLow);

	for(int p=0; p<numberOfProcessors; p++)
	  {
	    if(p!=rank)
	      {
		MPI_Send(&exchangeData, sizeof(ExchangeData), MPI_CHAR, p, 20001, MPI_COMM_WORLD);
	      }
	  }
	//deallocate anything we allocated to receive

	return 0;
#endif
}


int World::tagForTask(int srcI, int srcJ, int flag, int iscount )
{
#if MPI_ON
	int oceanNumber = srcI + sizeX*srcJ;
	int numOceans   = sizeX*sizeY;
	int typeNumber  = flag*numOceans + iscount*numOceans*numOceans*3;

	return typeNumber+oceanNumber;
#endif
}


void World::clearShadowOceans()
{
#if MPI_ON
  for (int i=0; i < sizeX; i++)
    for(int j=0; j< sizeY; j++)
      if (oceans[i][j].RunOnThread != rank)
	oceans[i][j].clear();
#endif
}


// rank here is rank if the current processor; i.e. the source of the message
// thread is the destination of the message
void World::MPI_SendFish(int srcI, int srcJ, int rank, int thread, int  flag)
{
#if MPI_ON
	int count;
	#ifdef  WhoSendsWhat
		printf("Thread [%d]: Sending message from Ocean[%d,%d] to thread %d with tag %d & flag %d.\n", rank, srcI, srcJ, thread, srcI*sizeX + srcJ, flag);
	#endif

	count = oceans[srcI][srcJ].constructList();

	MPI_Send(&count, 1, MPI_INT, thread, tagForTask(srcI, srcJ, flag, 1),  MPI_COMM_WORLD);
	MPI_Send(oceans[srcI][srcJ].fishRecordList, count * sizeof(FishRecord), MPI_CHAR, thread, tagForTask(srcI, srcJ, flag, 0), MPI_COMM_WORLD);

	return;
#endif
}


// rank here is rank of the current thread (i.e. the destination of the message
// thread is source of the message
void World::MPI_RecvFish(int targetI, int targetJ, int rank, int thread, int flag)
{
#if MPI_ON
	int count;
	FishRecord* recvFishRecord;
	MPI_Status status;

	#ifdef WhoReceivesWhat
	printf("Thread [%d]: Receiving message from Ocean[%d,%d] on thread %d to thread %d with tag %d & flag %d.\n",
		 rank, targetI, targetJ, rank, thread,  targetI*sizeX + targetJ, flag);
	#endif

    	MPI_Recv(&count, 1, MPI_INT, thread,  tagForTask(targetI, targetJ, flag, 1), MPI_COMM_WORLD, &status);

	if ( count < 0) {
		printf("Thread [%d]: Error in receiving count from thread %d for ocean[%d,%d].\n", rank, thread, targetI, targetJ);
		MPI_Finalize();
		exit(0);
	}

	recvFishRecord = (FishRecord*) malloc ( sizeof(FishRecord)*count );

	if ((recvFishRecord == NULL) && ( count != 0 )){
		printf("Thread [%d]: Error in allocating receiving buffer from thread %d for ocean[%d,%d].\n", rank, thread, targetI, targetJ);
		MPI_Finalize();
		exit(0);
	}


	MPI_Recv(recvFishRecord, count * sizeof(FishRecord), MPI_CHAR, thread, tagForTask(targetI, targetJ, flag, 0), MPI_COMM_WORLD, &status);

	if( flag == MIGRATING ) // it is a migrating fish
		for (int k=0; k < count; k++)
			oceans[targetI][targetJ].add(recvFishRecord[k], 0);

	else if ( flag == GHOST ) // it is a ghost fish
		for (int k=0; k < count; k++)
			oceans[targetI][targetJ].add(recvFishRecord[k], 1);

	free(recvFishRecord);
#endif
}


// --------------------------------------------------//
/* Here are all the functions that I added to        *
 * the class World                                   */
// --------------------------------------------------//

// Function to update the reproductive status of the fish
void World::setReproduction(int k, bool Reproduced)
	{
	for (int y=0; y<sizeY; y++)
		for (int x=0; x<sizeX; x++)
			{
			oceans[x][y].setReproduction(k, true);
			}
	}



int World::addFish(int TotalFishCount, int day, int nEggs)
        {
        for (int X=0; X<7; X++)
            for (int Y=0; Y<6; Y++)
            {
                //if (TotalFishCount < TOTALCC)
                //{
                LastUsedIDFish = oceans[X][Y].isBorn(LastUsedIDFish, day, nEggs); 
                //}
	}
        //printf("Last used ID %i", LastUsedIDFish);
        return LastUsedIDFish;
}

void World::totalToRepToday(int day)
{
        day = JulDay(day);
        howManyFishToRepToday = 0;
	for(int y=0; y<6; y++)
                {
                for(int x=0; x<7; x++)
                        {
                        howManyFishToRepToday += oceans[x][y].numberOfFishToRepToday(day);
                        //printf("Ocean %i, %i have %i to reproduce today (%i) \n", x, y, howManyFishToRepToday, day);
                        }
                }
        
}

void World::FishAging(int NaturalMortality, float toKillToday, float tJuv)
{
    //printf("%i\n", toKillToday);
    nOfAdults = 0;
	for(int x=0; x<7; x++)
                for(int y=0; y<6; y++)
                {
                   // printf("[%i, %i, %i]: ", x, y, oceans[x][y].numberOfFish);
                    oceans[x][y].getOld(NaturalMortality, toKillToday, tJuv);
                    //printf("Adults: %i, Juv: %i, Old: %i, [%i, %i]\n", oceans[x][y].nAdults, oceans[x][y].nJuveniles, oceans[x][y].nOld, x, y);
                }
    
}


void World::translateFtoGrid (int xDimension, int yDimension, float F_Mortality)
{
	for (int xIndex=0; xIndex<xDimension; xIndex++)
        {
            for (int yIndex=0; yIndex<yDimension; yIndex++)
            {
                (grid[xIndex][yIndex].FM = F_Mortality);
            }    
        }
}


/* This function gets called from gridtest. It in turn calls getTemperatureGradient */
int World::GetDensities(int gridSizeX, int gridSizeY)
{
    for (int gx = 0; gx<gridSizeX; gx++)
    {
        for (int gy = 0; gy<gridSizeY; gy++)
        {
            GridPoint& G(grid[gx][gy]);
            G.AdultDensity = 0;
            G.density = 0;
            G.age0 = 0;
            G.age1 = 0;
            G.age2 = 0;
            G.age3 = 0;
            G.age4 = 0;
            G.age5 = 0;
            //printf("%i, %i, %i, %i \n", G.age0, G.age1, G.age2, G.age3);
        }
    }
            
	//Iterates through the oceans:
        nOfAdults = 0;
	for(int y=0; y<sizeY; y++)
	for(int x=0; x<sizeX; x++)
	{
		Ocean& O(oceans[x][y]);
		//Iterates through the fish of the oceans:
		for(int i=0; i<O.numberOfFish; i++)
		{
			Fish& F(O.fish[i]);
			F.findNearestGridPoint(gridSizeX, gridSizeY);
			GridPoint& G(grid[F.nearestGridPoint.x][F.nearestGridPoint.y]);
                        G.FCatches = 0;
			G.density++;
                        
                        if (F.age >= AGEMATURE)
                        {
                            G.AdultDensity++;
                            howManyFishToRep++;
                            nOfAdults++;
                        }
                        // Obtain the frequency for each age class
                        if (F.age < 365)
                        {
                            G.age0++;
                        }
                        else if (F.age >= 365 and F.age < 365*2)
                        {
                            G.age1++;
                        }
                        else if (F.age >= 365*2 and F.age < 365*3)
                        {
                            G.age2++;
                        }
                         else if (F.age >= 365*3 and F.age < 365*4)
                        {
                            G.age3++;
                        }
                         else if (F.age >= 365*4 and F.age < 365*5)
                        {
                            G.age4++;
                        }
                         else if (F.age >= 365*5)
                        {
                            G.age5++;
                        }       
		}       
	}
       MaxDensAdults = 0;
       MaxNBoats = 0;
       for (int j=0; j<gridSizeX; j++)
            for (int k=0; k<gridSizeY; k++)
            {    
                GridPoint& G(grid[j][k]);
                //nOfAdults += G.AdultDensity;

                if (G.AdultDensity > MaxDensAdults) 
                {
                   MaxDensAdults = G.AdultDensity;  
                }
                if (G.boats > MaxNBoats)
                {
                    MaxNBoats = G.boats;
                }
            }
       return(nOfAdults);
} //World::factorInTemperatureGradients


void World::FishingMortality(int i, int gridSizeX, int gridSizeY)
{
	//Iterates through the oceans:
	for(int y=0; y<sizeY; y++)
        {
            for(int x=0; x<sizeX; x++)
            {
                    Ocean& O(oceans[x][y]);
                    //Iterates through the fish of the oceans:

                    for(int z=0; z<O.numberOfFish; z++)
                    {
                            Fish& F(O.fish[z]);
                            
                            if (F.age>=AGEMATURE)
                            {
                                F.pp = (float)rand()/RAND_MAX;
                                F.findNearestGridPoint(gridSizeX, gridSizeY);
                                GridPoint& G(grid[F.nearestGridPoint.x][F.nearestGridPoint.y]);
                                if (G.BoatsXFMpp > 0)
                                {
                                    if (F.pp <= G.BoatsXFMpp)
                                    {
                                        F.removeMe = true;
                                        G.FCatches++;
                                    }
                                }
                            }
                    }
                    O.removeMarkedFish();
            }
        }
        int totalCatch = 0;
        int pp = 0;
        int temp = i * TIMESTEP;

        for (int j=0; j<gridSizeX; j++)
            for (int k=0; k<gridSizeY; k++)
            {
                totalCatch += grid[j][k].FCatches;
                YearlyAccumulatedLandings += grid[j][k].FCatches;
                grid[j][k].AnnualCatches  += grid[j][k].FCatches;
                grid[j][k].FCatches = 0;
                if (temp%365==0) {
                    grid[j][k].AnnualCatches = 0;
                }
            }
} //World::factorInTemperatureGradients


//Write info in the grid to a file...

void World::WriteXYDensity(char file[255], int gridSizeX, int gridSizeY, char sst[20], 
        char cur[20], char interact[20], char manag[20], char _managX[20], char _managY[20])
{
        std::ofstream myfile;
        myfile.open(file);
        // This will be the header
        myfile << "This file was created with: " << manag << ", " << interact << ", " << sst  << ", " << cur << "\n";
        myfile << "The MPA covers for X and Y (center and length): " << _managX << ", " << _managY << "\n";
        myfile << "x, y, temperature, totDens, aduDens, MaxDens, AnnualCatches, nBoats, CorrectedFM, Bp, TotalBoatPP, Class0, Class1, Class2, Class3, Class4, Class5+  \n";
	//Iterates through the oceans:
	for(int y=0; y<gridSizeY; y++)
	for(int x=0; x<gridSizeX; x++)
	{
		GridPoint& G(grid[x][y]);
                double temper  = G.temperature;
                int    catches = G.AnnualCatches;
                int    totDens = G.density;
                float  maxDens = MaxDensAdults;
                int    aduDens = G.AdultDensity;
                float  nBoats  = G.boats;
                float  CorrFM  = G.BoatsXFMpp;
                float  Bp      = G.bp;
                float  BoatPP  = totalBoatPP;
                double Class0    = G.age0;
                double Class1    = G.age1;
                double Class2    = G.age2;
                double Class3    = G.age3;
                double Class4    = G.age4;
                double Class5    = G.age5; 
                myfile << x << ", " << y << ", " << temper << ", " << totDens << ", " << aduDens 
                        << ", " << maxDens <<", " << catches<< ", " << nBoats 
                        << ", " << CorrFM <<", " << Bp << ", " << BoatPP << ", " 
                        << Class0 << ", " << Class1 << ", " << Class2 << ", "
                        << Class3 << ", " << Class4 << ", " << Class5 << "\n" ;
	}
        myfile.close();
        printf("File %s wrote!!\n", file);
} 

void World::WriteGridFishIDs(char file[255], int day, int gridSizeX, int gridSizeY, double minMPAx, double maxMPAx, double minMPAy, double maxMPAy)
{
        std::ofstream myfile;
        myfile.open(file);
        // This will be the header
        myfile << "This file was created to find the residence time for day: " << day << "\n";
        //myfile << "The MPA covers for X and Y (center and length): " << _managX << ", " << _managY << "\n";
        myfile << "Day, x, y, fishID, inMPA \n";
	//Iterates through the oceans:
	for(int y=0; y<sizeY; y++)
        for(int x=0; x<sizeX; x++)
            {
                Ocean& O(oceans[x][y]);
                //Iterates through the fish of the oceans:

                for(int z=0; z<O.numberOfFish; z++)
                {
                        Fish& F(O.fish[z]);
			int inMPA = 0;
                        if (F.age>=AGEMATURE)
                        {
                            F.findNearestGridPoint(gridSizeX, gridSizeY);
                            int x = F.nearestGridPoint.x;
                            int y = F.nearestGridPoint.y;
                            int iDs = F.ID;
			    if (x >= minMPAx & x <= maxMPAx & y >= minMPAy & y <= maxMPAy)
				{
				inMPA = 1;
				}
                            myfile << day << ", " << x << ", " << y << ", " << iDs << ", " << inMPA <<"\n" ;
                        }
                }
            }
        myfile.close();
        printf("ResTime File %s wrote!!\n", file);
}
       
void World::IniFishMovement()
{
    for (int x = 0; x < 7; x++)
        for (int y = 0; y < 6; y++)
        {
            Ocean& O(oceans[1][0]);
            for(int z=0; z<O.numberOfFish; z++)
            {
                Fish& F(O.fish[z]);
                F.setDirection(20);
                F.setSpeed(.3);
            }
        }
}

void World::printFish(int ID)
{
    for (int x=0; x<7; x++)
        for (int y=0; y<6; y++)
        {
            oceans[x][y].removeMarkedFish();
            for (int k=0; k<oceans[x][y].numberOfFish; k++)
                {
                    Fish& F(oceans[x][y].fish[k]);
                    if (F.ID == ID)
                    {
                        printf("%2.2f, %2.2f, %2.2f, %2.2f, %2.2f \n", F.x, F.y, F.cosPhi, F.sinPhi, F.speed);
                    }
                }
        }
}

void World::BoatsDist(int i, int gridSizeX, int gridSizeY, double minMPAx, 
        double minMPAy, double maxMPAx, double maxMPAy, int YearOfMPA, double boatsC)
{       
    int day = i / 20;
    int year = day / 365;
    //Loop to get the "fraction" of boats in each pixel
    totalBoatPP = 0;
    for (int j=0; j<gridSizeX; j++)
    {
        for (int k=0; k<gridSizeY; k++)
            {    
                GridPoint& G(grid[j][k]);
                G.bp = exp(-boatsC * (1-(G.AdultDensity/MaxDensAdults)));
                if (year >= YearOfMPA)
                {
                    if (j >= minMPAx & j <= maxMPAx & k >=  minMPAy & k <= maxMPAy) // Check if the pixel is inside the MPA
                    {
                       G.bp = 0;
                    }
                }
                if (j >= 40)
                {
                    G.bp = 0; //This is to avoid putting boats over land....
                }
                totalBoatPP = totalBoatPP + G.bp;
            }
    }
    // This loop transform the relative number of boats into a proportion of the total number in each pixel 
    for (int j=0; j<gridSizeX; j++)
    {
        for (int k=0; k<gridSizeY; k++)
        {
            GridPoint& G(grid[j][k]);
            G.boats = totOfBoats * G.bp / totalBoatPP;
            G.BoatsXFMpp = (G.FM * 45  * G.boats/MaxNBoats); //45 is the scaling factor that characterizes the likelihood 
                //of capture for a given density of fish and a given effort of fishing calculated as in White C, Costello C. 
                //Close the High Seas to Fishing? PLoS Biol. 2014;12: e1001826. doi:10.1371/journal.pbio.1001826
	    if (G.FM == 0) {G.BoatsXFMpp = 0;}
            if (G.boats == 0) {G.BoatsXFMpp = 0;}
           
             if (year >= YearOfMPA)
                {      
                 if (j >= minMPAx & j <= maxMPAx & k >=  minMPAy & k <= maxMPAy) // Check if the pixel is inside the MPA
                    {
                     G.BoatsXFMpp = 0;
                    }
                }
            //printf("[%i, %i] %2.3f  %2.4f %i %i %2.3f \n", j, k, G.boats, G.bp, G.AdultDensity, MaxDensAdults, 5000*G.bp);
        }
    }   
} //World::factorInTemperatureGradients

void World::BoatsDistNCatches(int i, int gridSizeX, int gridSizeY, double minMPAx, 
        double minMPAy, double maxMPAx, double maxMPAy, int YearOfMPA, float NCatches, double boatsC)
{       
    int day = i / 20;
    int year = day / 365;
    //Loop to get the "fraction" of boats in each pixel
    totalBoatPP = 0;
    for (int j=0; j<gridSizeX; j++)
    {
        for (int k=0; k<gridSizeY; k++)
            {    
                GridPoint& G(grid[j][k]);
                // Here I need to add a flag to check for the MPA!!!              
                G.bp = exp(-boatsC * (1-(G.AdultDensity/MaxDensAdults)));
                if (year >= YearOfMPA)
                {
                    if (j >= minMPAx & j <= maxMPAx & k >=  minMPAy & k <= maxMPAy) // Check if the pixel is inside the MPA
                    {
                       G.bp = 0;
                    }
                }
                if (j >= 40)
                {
                    G.bp = 0; //This is to avoid putting boats over land....
                }
                totalBoatPP = totalBoatPP + G.bp;
            }
    }
    // This loop transform the relative number of boats into a proportion of the total number in each pixel 
    for (int j=0; j<gridSizeX; j++)
    {
        for (int k=0; k<gridSizeY; k++)
        {
            GridPoint& G(grid[j][k]);
            G.boats = totOfBoats * G.bp / totalBoatPP;
            if (G.AdultDensity == 0)
            {
                G.BoatsXFMpp = 1;
            } else if (G.AdultDensity > 0)
            {
                G.BoatsXFMpp = (G.boats * NCatches) / G.AdultDensity;
            }
	    if (G.FM == 0) {G.BoatsXFMpp = 0;}
            if (G.boats == 0) {G.BoatsXFMpp = 0;}
           
             if (year >= YearOfMPA)
                {      
                 if (j >= minMPAx & j <= maxMPAx & k >=  minMPAy & k <= maxMPAy) // Check if the pixel is inside the MPA
                    {
                     G.BoatsXFMpp = 0;
                    }
                }
            //printf("[%i, %i] %2.3f  %2.4f %i %i %2.3f \n", j, k, G.boats, G.bp, G.AdultDensity, MaxDensAdults, 5000*G.bp);
        }
    }   
} //World::factorInTemperatureGradients


//This function is to check if in the previous iteration the TAC was reached.
double World::checkTAC(int i, int nOfAdults, bool tacFlag)
{
    int day = i / 20;
    int year = day / 365;
    
    if (day%365 == 0 | i == 120) 
    {
        #if USEmsy 
                TAC = r0 * (TOTALCC) / 4;
        #endif 
        if (tacFlag)
        {
                TAC = BaseTAC + TACpp * nOfAdults;  
        }
        printf("The TAC is now: %2.2f \n", TAC);
                
    }
    return(TAC);
}

// Function to change the mortality rate of the larvae to recruit the
// correct number of organisms corresponding to the Gordon-Schaefer model
 float World::adjustLarvaeMortality(int type, int nOfAdults, int nOfJuv)
 {
     
     //Every day we kill 1/365 of the larvae - recruits as predicted by
     // Gordon Schaefer model.
     if (type == 0) 
     {
         float k = TOTALCC;
         float ccFactor = nOfAdults / k;
         float nt = nOfAdults + nOfAdults * r0 * (1 - ccFactor);
         float recruit = nt - nOfAdults;
         //printf("%2.2f, %2.2f, %i - ", nt, recruit, nOfJuv);
         
         if (recruit < 0) {recruit = 0;}
         float tempRecToJuv = (recruit / nOfJuv);
         float tempDayBase = pow(300, -1); //365^(-1);
         
         float ppSurvive = pow(tempRecToJuv, tempDayBase);
         if (isnan(ppSurvive)) {ppSurvive = 0;}
         float ppToKill = 1 - ppSurvive ; /// (nOfJuv - recruit);
         if (ppToKill < 0){ppToKill = 0;}
         printf("Adjusted Larvae Mortality:  %2.4f \n", ppToKill);
         return(ppToKill);
     }
 }

 double World::mpaDef(int mpaDef, int mpaSize)
 {
    double MPArange = mpaSize / 2;
    double MPAcenter = 0;
    //No moving MPA    
    if (mpaDef == 0 or mpaSize == 100)
    {
        printf("No Moving MPA, defined at the center of size %i \n", mpaSize);
        MPAcenter = 50;
    }
    //Basic block moving MPA
    if (mpaDef == 1 and mpaSize < 100)
    {
        int distFish[GRIDSIZEY] = {0};
        int blockFish[GRIDSIZEY] = {0};
        int yLargest = 0;
        int largest  = 0;
        for(int y=0; y<GRIDSIZEY; y++)
            for(int x=0; x<GRIDSIZEX; x++)
            {
                    GridPoint& G(grid[x][y]);
                    distFish[y] += G.AdultDensity;
            } 
        for (int i=MPArange; i<(GRIDSIZEY-MPArange); i++)
        {
            for (int j = (i-MPArange); j < (i + MPArange); j ++)
            {
                blockFish[i] += distFish[j];
            }
        }
        for (int k = 0; k < GRIDSIZEY; k++)
        {
            if (yLargest <= blockFish[k])
            {
                yLargest = blockFish[k];
                largest  = k; 
            }
        }
        MPAcenter = largest;
        if ((MPAcenter + MPArange) >= 100)
        {
            MPAcenter = 100 - MPArange;
        }
        if ((MPAcenter - MPArange) < 1)
        {
            MPAcenter = 1 + MPArange;
        }
        printf("MPA block of size %i, center at: %i \n", mpaSize, largest);
    }
       
    return(MPAcenter);
 }
 
 // This functions calculates the Kuramoto parameter (R) for the world!
rVector World::Rparam(int HowManyFish, int nOfAdults)
{
    float sumCosPhi = 0;
    float sumSinPhi = 0;
    float reipsi1 = 0;
    float reipsi2 = 0;
    float psi = 0;
    rVector R;
    float sumCosPhiAdults = 0;
    float sumSinPhiAdults = 0;
    float reipsi1a = 0;
    float reipsi2a = 0;
    float psiA = 0;
    int count = 0;
    R.all = 0;
    R.adults = 0;
    for (int x=0; x<7; x++)
    {
        for (int y=0; y<6; y++)
        {
            for (int k=0; k<oceans[x][y].numberOfFish; k++)
            {
                Fish& F(oceans[x][y].fish[k]);
                sumCosPhi += F.cosPhi;
                sumSinPhi += F.sinPhi;
                if (F.age>=AGEMATURE)
                {
                    sumCosPhiAdults += F.cosPhi;
                    sumSinPhiAdults += F.sinPhi;
                }
            }   
        }
    }
    reipsi1 = sumCosPhi / HowManyFish;
    reipsi2 = sumSinPhi / HowManyFish;
    reipsi1a = sumCosPhiAdults / nOfAdults;
    reipsi2a = sumSinPhiAdults / nOfAdults;

    R.all = sqrt(reipsi1*reipsi1 + reipsi2*reipsi2);
    R.adults = sqrt(reipsi1a*reipsi1a + reipsi2a*reipsi2a);

    printf("#Fish: %d [%i] and Rall is:  %2.2f , Radults is: %2.2f \n", HowManyFish, nOfAdults, R.all, R.adults);
    return(R);
}

 void World::writeRparam(int year, int day, int howManyFish, int nOfAdults, rVector R)
 {
     if (day == 1)
     {
        ofstream myfile;
        myfile.open("./Rparam.csv");
        myfile << "R parameter for year: " << year << "\n";
        myfile << "Year, Day, n_all, n_adults, R_all, R_adults  \n";
        myfile << year << ", " << day << ", " << howManyFish << ", " << nOfAdults << ", " << R.all << ", " << R.adults << "\n";
        myfile.close();
     }
     else 
     {
         std::ofstream myfile("./Rparam.csv", std::ios_base::app | std::ios_base::out);
         myfile << year << ", " << day << ", "  << howManyFish << ", " << nOfAdults << ", " << R.all << ", " << R.adults <<  "\n";
         myfile.close();
     }
 }
        
