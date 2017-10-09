// Kinetic Slug Model.cpp : Defines the entry point for the console application.
// Rotation, translation, and everything!
// Simpler SlugModelRot getting rid of viscoelasticity and FV properties
// Simpler SlugModelNeuralInput putting neural recordings into the model
// doing scans for square waves across fitness, with fitness calculations for biting, swallow
// and rejection, without changing the biomechanics.
// XXII model with dynamic hinge.

#include <iostream>
#include <fstream>
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;

#include <cstdlib> // for exit function


FILE *valout;
FILE *neuralinputs;

/* including some files of my own in the rebulid */
#include <string.h>
#include "math.h"
#include "stdio.h"
#include <vector>
#include <sstream>

/* The constants */
#define eom 1 /*Equation of Motion, if this is one - using complex equation of motion (NEOM)
									if this is two - using simple equation of motion (SEOM)
									if this is three - using the new simple equation of motion (NSEOM)*/

#define inputcheck 1 /*if this equals 1 will let you input neuron values into printmodel (ignores GA) 
						if this does not equal one, code will run a GA */

#define protractfit -1 /*if this equals 1 fitcalc will give fitness for protraction
						if this does not equal one, fitness will be given at retraction only*/

// TATE
double valueTBD;

const double pi = 3.14159265359;
const double StepSize = .000008;  //.00001; size of the time step used to evaluate the Euler integrators in the model
const double timer = 0.01;  //printmodel will print to the text file every "timer" seconds - keeps text file from being too big

//Newton/Raphson constants
const double delta = 0.00001;
const double Error = 0.01;

//forces in mN, distances in meters, masses in grams
//odontophore mass and dimensions taken from hinge individual "Usurp"
const double Radius = .005; //radius of sphere odontophore in meters
const double r = 0.00125; //inner radius of the I1/I3 torus in m, remains constant in this model
const double Odonmass = 1.7; //mass of sphere in grams (labeled M in Guidebook)
const double I1I3mass = 1.3; //mass of torus in grams (labeled m in Guidebook)

//Free parameters
const double LI2opt = 0.0225; //.0825; //optimal length of I2
const double LI1I3opt = 0.0417;//0.04; //optimal length of I1/I3

//activation parameters, from Yu et al 1999
const double tauact = 2.45; //sec
const double beta = 0.496;//0.703 <- reported in paper, this was with a tendon, the tendonless model is used in the code
const double Ao = 0.165;
const double g = 1; //different from Yu et al because of a typo in Yu et al.

const double uprimebound = -2.0; //uprime was bounded in the paper at 0, to fit the data better on deactivation this was changed
                                  //Original was -0.3, seeing how this affects in vivo records.
//muscle parameters
const double FnotI2 = 240;//mN - maximum tension that I2 can produce, different from Yu to reproduce figure 2A
const double FnotI1I3 = 720; //mN - maximum tension that I1/I3 can produce

//visco-elastic hinge force parameters, from Sutton et al 2003 (hinge paper)
const double A1 = 0.617; //in mN/mm^2
const double A2 = 24.17; //in mN/mm
const double Xo = 13.27;
const double DX = 2.422;
const double B1 = -58.7; //in mN/mm^2
const double B2 = 12.6; //in mN/mm
const double B3 = -0.895; //in mN
const double B4 = 0.024;
const double C1 = 0.0997; //in mN/mm^2
const double C2 = 0.27;
const double C3 = -0.009;
const double D1 = -0.363;
const double D2 = 1.096; 
const double E1 = 1.69;
const double E2 = -8.19;
const double E3 = 11.03;
const double AF1 = 0.288;
const double AF2 = 6.78;
const double viscdamp = 4.0;  //  simple model it was 4.0;
const double xp = 3.4; //mm, point in protraction when hinge force is first measured
const double xr = 6.9; //mm, point in retraction when hinge force is last measured
const double passiveoffset = 9.8; //shifts the hinge forces in order to get a biologically accurate rest position
const double lengthshift = 0.0;  // 0.0; //-.0062;  the shift on the hinge LT in order to get good active forces.

const double MAXSEAWEEDFORCE = 100;  // the maximum force on the seaweed in mN.

const int seednum = 42;
const double RunDuration = 8.5; // 8.5; //seconds, length of time each individual is run

const int inputrows = 740;
const int NeuronNum = 4; //number of neurons in the network

const double anglestiffness = .0001; //maximum additional rotation speed for odontophore in deg/second

/* Variables used for opening input file*/
ifstream inputFile;
string times;
string positions;
string radiuss;
string angles;
string hingeFs;
string fitnesss;
string freqI2s;
string freqI1I3s;
string freqN3s;
string freqHinges;
string actI2s;
string actI1I3s;
string acthinges;
string seaweedforces;

/*Arrays of stored values from input file*/
double positionarray[1000];
double radiusarray[1000];
double anglearray[1000];
double hingeFarray[1000];
double fitnessarray[1000];
double freqI1I3array[1000];
double freqN3array[1000];
double freqHingearray[1000];
double actI2array[1000];
double actI1I3array[1000];
double acthingearray[1000];
double seaweedforcearray[1000];
double freqI2array[1000];
double timearray[1000];

/* Value of columns in the input file*/
int numberColumns = 0;

/* My function prototypes */
//calling the functions
	
double updateforces (double x, double ytop, double ybottom, double a, double ydot, double xdot, double oldx, double oldytop, double oldybottom, double olda, 
					 double freqI2, double freqI1I3, double *aprimeI2, double *aprimeI1I3, double *force1, double *force2, 
					 double *F1, double *HingeAct, double hingefrequency, double Oangle, double *hingeforce, double *externalforce);
  //updateforces takes the current position of the odontophore and I1/I3 and the frequency of stimulation of the muscles 
  //and calculates the forces on the odontophore

double musclelengthI2 (double a, double x, double ytop, double ybottom, double Oangle); 
	//takes the position of the odontophore and calculates the length of the I2 muscle
	
double fixI2length (double); 
	//adjusts I2 length to give results that are consistent with the biology

double musclelengthI1I3 (double ytop, double ybottom); 
	//takes the position of I1/I3 and calculates the length of I1I3

double musclevelocityI2 (double oldx, double oldytop, double oldybottom, double olda, double x, double ytop, double ybottom, double a, double Oangle); 
	//takes the current and previous positions of the odontophore and I1/I3 and calculates the velocity of I2 

double musclevelocityI1I3 (double ytop, double ybottom, double oldytop, double oldybottom); 
	//takes the change in position of I1/I3 and calculates the velocity of I1/I3

double activationI2 (double frequency, double *aprime); 
	//calculates the activation of I2 given the frequency of stimulation
	
double activationI1I3 (double frequency, double *aprime);
	//calculates the activation of I1/I3 given the frequency of stimulation
	
double activationN3(double frequency, double *act, double time);
	//calculates the activation of the shape change given the output frequency of N3

double secondactivationN3(double frequency, double *act, double time);
	//calculates the activation of the shape change given the output frequency of N3
	
double passive (double); 
	//calculates the passive forces of a muscle given its normalized length
	
double lengthtens (double); 
	//calculates the length tension of a muscle given its normalized length

double forcevelocity (double); 
	//calculates the force velocity of a muscle given its velocity

double calchingeforce (double, double, double *); 
	//calculates the passive hinge force given the current position and velocity of the odontophore

double calchingeforce2 (double, double);

double updateposition (double *x, double *ytop, double *ybottom, double *a, double *xdot, double *ydot, double *adot, double f1, double f2, double actN3, double *Txc, double *Bxc, double *output, double *Oangle, double hingeforce); 
	//takes the forces on the odontophore and finds the new position of the odontophore and I1/I3

double calcphi (double x, double ytop, double ybottom, double a, double *Tphi, double *Bphi, double Oangle);
	//calculates the angle of the straight piece of I2 with respect to the horizontal
	//this is used in the calculation of I2's force

double contactslope (double *, double *, double *, double *, double *, double Oangle);
	//calculates the slope of the tangent angle at the point of contact between the odontophore and
	//I1/I3 and updates the point of contact (xc, yc)

void fitcalc (double x, double a, double oldx, double olda, double & fitness);
	//updates an individual's fitness

double power(double, double);
	//takes the first number and raises it to the second number's power (positive powers only)

double topEllipseSlope(double a, double xc, double rotation, double & yc);

double topEllipseSlope2(double a, double xc, double rotation, double & yc);

double bottomEllipseSlope(double a, double xc, double rotation, double & yc);

double topI1I3slope(double x, double xc);
	/*slope of top I1I3 ring (bottom circle curve*/

double bottomI1I3slope(double x, double xc);
	/*slope of bottom I1I3 ring (top circle curve*/ 

double OdonAngle(double x);
	/*gives the odontophore's rotation given its position*/

double OdonAngle2(double x, double hingeforce, double radius, double oldodonangle);
    /* This code rotates the odontophore as a function of hinge force and current radius */
	
void xcCalc(double a, double x, double rotation, double & Txc, double & Bxc, double & Tyc, double & Byc, double & topslope, double & bottomslope);

double activehingeforce (double activation, double velocity, double length);

void updateinputs(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2, string behaviorType);

void updateinputsRejectionB(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

void updateinputsBite(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

void updateinputsSwallowB(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

void updateinputsSwallowA(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

void updateinputsRejectionA(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

void updateinputsSwallowPerturbed(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

bool openAndRead(string behaviorType);

int timeAdjuster(double timeArray[], double timeStamp);

void updateDynamicInputs(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2);

string returnFirstLine(string s);

int counter(string str);

int columnChecker();

int main(int argc, char* argv[])
{  
    
    /* Taking an argument, parsing it to a string*/
    std::string programName = argv[0];
    std::string first_argument;
    std::vector<std::string> all_args;
    if (argc > 1) {
        first_argument = argv[1];
        all_args.assign(argv + 1, argv + argc);
    }
    
    /* Opening and reading a file*/
    numberColumns = columnChecker();
    openAndRead(first_argument);
    if(openAndRead(first_argument) == true){

    /* Absent minded code variable definitions */
	const char* filename = "SlugOutput2.txt";
	const char* filename2 = "NeuralInputs.txt";

    //variables used to calculate the muscle forces and odontophore position - explained when initialized
    //variables used to calculate the muscle forces and odontophore position - explained when initialized
    double x, y, xdot, ydot, adot, xacc;
    double a, b;
    double odontophoreangle;
    double vol, force1, force2, time;
    double oldx, oldy, olda, oldytop, oldybottom;
    double aprimeI2, aprimeI1I3, aprimeN3, freqI2, freqI1I3, freqN3, actN3;
    double F1; //visco-elastic hinge force, F1 from Sutton et al 2002
    double hingeF = 0;
    double xctop, xcbottom, yctop, ycbottom, ytop, ybottom, topslope, bottomslope, rotation;
	
    //fitness calculation variables - explained when initialized
    double fitness;
    double SArea;
    double fitcounter1, fitcounter2, fitcheck;
    double fitpertime, oldfit, oldtime, totalfitness;
    double oldxdot;
    double aPrevious = 0;
    double printtimer;
    double eggtimer = 0;
    double actI3 = 0;
    double acthinge = 0;
    int i = 1;
    int I2inputcounter = 1;
    int j = 1;
    int k = 1;
    int I1I3inputcounter1, I1I3inputcounter2, I1I3inputcounter3, RNinputcounter, Hingeinputcounter;
    double printvar = 0;
    double printvar2 = 0; // used in update position to ou tput the effective slope and phi
		
    /* Absent minded code variable definitions */
    int columncounter = 0;
    int rowcounter = 0;
    ifstream indata;  //the input stream
    float num; //variable for input value
    float neuralvariables[12][inputrows];
    double I1I3metafreq1 = 0;
    double I1I3metafreq2 = 0;
    double I1I3metafreq3 = 0;
    double freqHinge = 0;
    double startloop1 = 1;
    double frequencyiterationtime = startloop1; //Code added to run loops of frequency code
    double endfrequencytime = 1.1;  //BLARF CHANGED TO ONLY CYCLE ONCE
    double frequencystep = .1;
	double frequencyiterationtime2 = 0;
    double endfrequencytime2 = 0.1;
    double frequencystep2 = .1;  //BLARF CHANGED TO ONLY CYCLE ONCE
    double checker = 1;
    double period = 0;
    double dummyvariable = 0;  //just something to hand to functions for outputs that dont matter
    ifstream inFile;
    valout = fopen(filename, "w"); //opens a text file, valout.txt
    double seaweedforce = 0;  //This is added seaweed force on the odontophore

    //********INITIALIZING THE VARIABLES********
		
    //variables for the physics - see geometry diagram in Valerie Snyder's notebook 4/24/02
    //Guidebook pg. 1 diagram
    a = Radius; //the major axis of the odontophore, the odontophore is intitialized as a sphere
    b = (5.0*sqrt(5.0)/sqrt(a * 1000.0))/1000.0; //the minor axis of the odontophore, equation keeps the odontophore isovolumetric
    x = -1 * a; //displacement of center of odontophore with respect to vertical line down the center of the I1/I3 torus
    y = 0.00375; //major radius of I1/I3 ring
    //r = sqrt(x*x + y*y) - Radius; //Pythagorean equality, r = thickness of I1/I3 torus (minor/cross-sectional radius)
    vol = 2 * pi *r*r*y; //volume of the I1/I3 torus
    xdot = 0; //velocity of odontophore
    ydot = 0; //velocity of contraction/expansion of I1/I3 torus
    adot = 0; //velocity of major axis of the odontophore
    xacc = 0; //acceleration of the odontophore
    //rotation = OdonAngle(x);
    odontophoreangle = 90;
    xctop = -0.0037;
    xcbottom = -0.0037;
		
    xcCalc(a, -1*x, odontophoreangle, xctop, xcbottom, yctop, ycbottom, topslope, bottomslope);
	
    ytop = yctop + sqrt(r*r - (-1*x - xctop)*(-1*x - xctop));
    ybottom = ycbottom - sqrt(r*r - (-1*x - xcbottom)*(-1*x - xcbottom)); //initializing contact point
		
    force1 = 0.000; //horizontal force acting on odontophore
    force2 = 0.000; //vertical force acting on the I1/I3 torus
	
    //oldx, oldy, and oldr keep track of the previous position of the individual in order to calculate the velocity of I2 and adot
    oldx = x;
    oldy = y;
    olda = a;
    oldytop = ytop;
    oldybottom = ybottom;
    
    //activation variables, from Yu et al 1999
    aprimeI2 = 0;
    aprimeI1I3 = 0;
    aprimeN3 = 0.01;
    actN3 = 0;

    //freqI2 and freqI1I3 are the frequency of stimulation of I2 and I1/I3 respectively
    freqI2 = 0;
    freqI1I3 = 0;
	freqN3 = 0;  //controls the length of the odontophore's major axis
		
    //component of the visco-elastic hinge force, F1 from Sutton et al 2002
    F1 = 0;
	
    //fitness calculation variables
    //fitness = 0.005 ; //0.1; //running total of the fitness an individual has earned
    fitness = 0.0;  //starting fitness at zero because we're not doing GA's right now.
    SArea = 0; //amount of surface area exposed when the individual reaches peak protraction
    fitcounter1 = 0; //counts consecutive negative velocities in order to determine when the protraction is finished
    fitcounter2 = 0; //counts consecutive positive velocities in order to determine when retraction is finished
    fitcheck = 0;  //keeps track of whether the individual is currently protracting (1) or retracting (-1)
    fitpertime = 0; //amount of fitness an individual recieves per second
    oldfit = 0; //used to store the amount of fitness at the last peak protraction to calculate fitpertime
    oldtime = 0;  //used to store the time when the last peak protraction occurred to calculate fitpertime
    totalfitness = 0; //fitness an individual would receive if the run were stopped (includes partial surface area or partial retraction)
    time = 0; //keeps track of the amount of time the model has run so model will stop after the RunDuration is up
    printtimer = 0;
    I2inputcounter = 0;
    I1I3inputcounter1 = 0;
    I1I3inputcounter2 = 0;
    I1I3inputcounter3 = 0;
    RNinputcounter = 0;
    Hingeinputcounter = 0;
		
		//neuralinputs = fopen(filename2, "r"); //sets up the text file for reading the input.

    fprintf(valout, "time	position	radius	angle	hingeF	fitness	freqI2	freqI1I3	freqN3	freqHinge	actI2	actI1I3	acthinge	fitness	seaweedforce\n");
		
    /* Reading the file for neural input */
    i = 0;
    j = 1;
		
        /*while (j<inputrows)
		{
		printf("Scan3 %i \n", j);

		i = j;
		printf("Scan4 %i \n", j);

		fscanf(neuralinputs, "%f	%f	%f	%f	%f	%f	%f	%f	%f	%f", &neuralvariables[1][i], &neuralvariables[2][i],&neuralvariables[3][i], &neuralvariables[4][i], &neuralvariables[5][i], &neuralvariables[6][i],&neuralvariables[7][i], &neuralvariables[8][i], &neuralvariables[9][i], &neuralvariables[10][i]);
		
		printf("Scan %i \n", j);
		j=j+1;
		printf("Scan2 %i \n", j);

		}*/

	/*indata.open("NeuralInputs.txt"); // opens the file

	 if(!indata) 
	  { // file couldn't be opened
		 cerr << "Error Spam: file could not be opened" << endl;
		 exit(1);
	  }
	 else
	 {
		 cout << "Found It!" << endl;
	 }  */
/*  BLARF COMMENT REMOVING READING INPUT
// Reading the data into the model
		indata >> num;

		while ( !indata.eof() ) 
		{ 
		// keep reading until end-of-file
		


		if (rowcounter > 11)
		{
			rowcounter = 0;
			columncounter++;
		}

		neuralvariables[rowcounter][columncounter] = num;
		
		//cout << "The next number is " << neuralvariables[rowcounter][columncounter] << endl;
		rowcounter = rowcounter +1;

	
		//cout << columncounter << endl;
		//cout << rowcounter << endl;

		indata >> num; // sets EOF flag if no value found
  	   }

   indata.close();
   cout << "End-of-file reached.." << endl;
			
		
BLARF COMMENT REMOVING READING INPUT END */
			
			//********ENDING INITIALIZATION********

/* Running the model here */
    printf("Dreaming the impossible dream \n" );
/*printf("%f	%f \n", neuralvariables[1][I2inputcounter], neuralvariables[2][I2inputcounter]);
printf("%f	%f \n", neuralvariables[1][1], neuralvariables[2][1]);

printf("%f	%f \n", neuralvariables[1][3], neuralvariables[2][3]);*/

/* WHILE LOOP REMOVED  while(frequencyiterationtime2 < endfrequencytime2) //a second loop to scan two dimensions
{
	
while(frequencyiterationtime < endfrequencytime)  //loop added to do cyclic frequency checks
{   */

//redo of initialization for the run
		
    //variables for the physics - see geometry diagram in Valerie Snyder's notebook 4/24/02
    //Guidebook pg. 1 diagram
    a = Radius; //the major axis of the odontophore, the odontophore is intitialized as a sphere
    b = (5.0*sqrt(5.0)/sqrt(a * 1000.0))/1000.0; //the minor axis of the odontophore, equation keeps the odontophore isovolumetric
    x = -1 * a; //displacement of center of odontophore with respect to vertical line down the center of the I1/I3 torus
    y = 0.00375; //major radius of I1/I3 ring
    //r = sqrt(x*x + y*y) - Radius; //Pythagorean equality, r = thickness of I1/I3 torus (minor/cross-sectional radius)
    vol = 2 * pi *r*r*y; //volume of the I1/I3 torus
    xdot = 0; //velocity of odontophore
    ydot = 0; //velocity of contraction/expansion of I1/I3 torus
    adot = 0; //velocity of major axis of the odontophore
    xacc = 0; //acceleration of the odontophore
    rotation = OdonAngle(x);
    xctop = -0.0037;
    xcbottom = -0.0037;
    xcCalc(a, -1*x, rotation, xctop, xcbottom, yctop, ycbottom, topslope, bottomslope);
    ytop = yctop + sqrt(r*r - (-1*x - xctop)*(-1*x - xctop));
    ybottom = ycbottom - sqrt(r*r - (-1*x - xcbottom)*(-1*x - xcbottom)); //initializing contact point
    force1 = 0.000; //horizontal force acting on odontophore
    force2 = 0.000; //vertical force acting on the I1/I3 torus
   
    //oldx, oldy, and oldr keep track of the previous position of the individual in order to calculate the velocity of I2 and adot
    oldx = x;
    oldy = y;
    olda = a;
    oldytop = ytop;
    oldybottom = ybottom;
		
    //activation variables, from Yu et al 1999
    aprimeI2 = 0;
    aprimeI1I3 = 0;
    aprimeN3 = 0.01;
    actN3 = 0;
		
    //freqI2 and freqI1I3 are the frequency of stimulation of I2 and I1/I3 respectively
    freqI2 = 0;
    freqI1I3 = 0;
    freqN3 = 0;  //controls the length of the odontophore's major axis
    freqHinge = 0;
		
    //component of the visco-elastic hinge force, F1 from Sutton et al 2002
    F1 = 0;
	
    //fitness calculation variables, making really low for bites!
    fitness = 0.0; //running total of the fitness an individual has earned, tube manipulation
    //fitness = x; //fitness for the maximum protraction trials of biting.

    time = 0;

        //TATE
        //printf("%f",freqI2array[0]);
        //printf("%f",freqHingearray[0]);
        //printf("%f",freqI1I3array[0]);
        //printf("%f",freqN3array[0]);
        
        
        
	while(time < RunDuration) //runs the individual until time is greater than RunDuration
		{
            
            
        //NeuronOutput is from 0 to 1, multip ly be 20 to get freq from 0 to 20
        freqI2 = 0;//circ.NeuronOutput(1) * 20;
        freqI1I3 = 0;//circ.NeuronOutput(2) * 20;
        freqN3 = 0;//circ.NeuronOutput(3) * 20;
        I1I3metafreq1 = 0;
        I1I3metafreq2 = 0;
        I1I3metafreq3 = 0;
        freqHinge = 0;
            
        //Update Neural variables and seaweed force based on the current time
        updateinputs(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime,  frequencyiterationtime2, first_argument);
            
            
        //BLARF looking at first protraction first
        //freqI1I3 = 0;
        //freqN3 = 0;
        //freqN3 = 15;

        eggtimer += StepSize;
				
        //actN3 = activationN3(freqN3, &aprimeN3, time);
        actN3 = secondactivationN3(freqN3, &aprimeN3, time);

        oldxdot = xdot; //for use in hinge fix
			
        //take frequencies and determine the forces
        //printvar2 = updateforces (x, ytop, ybottom, a, ydot, xdot, oldx, oldytop, oldybottom, olda, freqI2, freqI1I3, &aprimeI2, &aprimeI1I3, &force1, &force2, &F1, &acthinge, freqHinge);
        printvar = force1;

        //take frequencies and determine the forces
        printvar =  updateforces (x, ytop, ybottom, a, ydot, xdot, oldx, oldytop, oldybottom, olda, freqI2, freqI1I3, &aprimeI2, &aprimeI1I3, &force1, &force2, &F1, &acthinge, freqHinge, odontophoreangle, &hingeF, &seaweedforce);

        //oldx, oldy, and olda are set equal to x, y, and a respectively before x, y, and a are updated to the new position
        oldx = x;
        oldy = y;
        olda = a;
        oldytop = ytop;
        oldybottom = ybottom;
				
        //take forces and the frequency of N3 and determine odontophore position
        updateposition(&x, &ytop, &ybottom, &a, &xdot, &ydot, &adot, force1, force2, actN3, &xctop, &xcbottom, &dummyvariable, &odontophoreangle, hingeF);
								
        //resetting the visco-elastic force to zero to prevent accumulation of visco-elastic forces over large amounts of time
        //from overwelming the elastic hinge force
        if ((xdot/oldxdot) < 0)
        {
            F1 = 0;
        }
				
        eggtimer += StepSize;
        time += StepSize; //time advances one time step, run will stop when time becomes larger than the RunDuration
        printtimer += StepSize;

			
        if (printtimer > timer) //function will only print every printtimer seconds - keeps file from being too big
        {
            //print statement
            //fprintf(valout, "%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	\n",time, x, freqI2, freqI1I3, freqN3, actN3, a, fitness, aprimeI2 - Ao, aprimeI1I3 - Ao,  LI2opt*musclelengthI2 (a, x, y));
            //fprintf(valout, "%f	%f	%f	%f	%f	\n",x, printvar2, printvar, passive(printvar), passive(printvar)*2*pi*FnotI2);
            
            //fprintf(valout, "%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", x, a, OdonAngle(x), lengthtens(printvar), printvar2, force1 - calchingeforce(x, xdot, &F1), calchingeforce(x,xdot,&F1), lengthtens(musclelengthI1I3(ytop,ybottom)), printvar, force2 * printvar);
            
            //fprintf(valout, "%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	\n", time, x, ybottom, ytop, a, freqI2, freqI1I3, I1I3metafreq1, I1I3metafreq2, I1I3metafreq3, freqN3, freqHinge, -calchingeforce2(x, xdot),  -activehingeforce(acthinge, xdot, x), printvar2, dummyvariable, printvar, acthinge, aprimeI1I3, aprimeI2, OdonAngle(x), fitness);
						
            //This is the Fprint for the simulations in the example PDF's I sent Hillel.
            fprintf(valout, "%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f \n", time, x, a, odontophoreangle, hingeF, fitness, freqI2, freqI1I3, freqN3, freqHinge, aprimeI2, aprimeI1I3, acthinge, fitness, seaweedforce);

            //resets printtimer
            printtimer = 0.0;
        }

			
        if (time > 0) //fitness function ignores the initial jiggle into the equilibrium position
						   // because of the delay in the activation function, no muscle moves before time = 0.4
        {
            //updates the fitness / surface area / fitcheck
            fitcalc(x, a, oldx, olda, fitness); //commenting out the old fitness

            //for odontophore protraction magnitude
            /*if (x > fitness)
            {
                fitness = x;
            } */
            //fitcalcbite(x,a,oldx, olda, fitness);
        }
            // TATE TEST
            //if(time < 0.00008){
              //printf("%f \n", aprimeI2);
            //}
            //printf("%d \n", numberColumns);
    }

    printf("that's another simulation completed at %f	%f\n", frequencyiterationtime, frequencyiterationtime2);

		/* While loop Removal
    fprintf(valout, "%f	%f	%f \n", frequencyiterationtime, frequencyiterationtime2, fitness);
	printf("that's another simulation completed at %f	%f\n", frequencyiterationtime, frequencyiterationtime2);
    frequencyiterationtime = frequencyiterationtime + frequencystep;

	}

	frequencyiterationtime = startloop1;
    frequencyiterationtime2 = frequencyiterationtime2 + frequencystep2;

	} While loop Removal end */

    
    }
    else{
        printf("Please use a correct 'windows formatted text' file with either 6 or 15 columns \n");
    }
}



	/* The absent minded code that just runs after everything else works */
 
	

/* Putting in my function definitions */
double updateforces (double x, double ytop, double ybottom, double a, double ydot, double xdot, double oldx, double oldytop, double oldybottom, double olda, 
					 double freqI2, double freqI1I3, double *aprimeI2, double *aprimeI1I3, double *force1, double *force2, 
					 double *F1, double *HingeAct, double hingefrequency, double Oangle, double *hinge, double *externalforce)
	
	/*updateforces takes the new and old position of the odontophore (x, oldx), the new and old radius of the I1/I3 torus 
(y, oldy), the new and old major axes (a, olda), the velocity of the odontophore (xdot) velocity of the radius of the torus 
(ydot), and the frequency of muscle stimulation for I2 and I1/I3 (freqI2, freqI1I3).  The function uses these to calculate 
the new horizontal force on the odontophore (f1 in the model, force1 in function) and the new vertical force on the I1/I3 
torus (f2 in the model, force2 in the function).  The function first calculates the lengths of I2 and I1/I3 using 
musclelengthI2 and musclelengthI1I3.  I2's length is renormalized using the function fixI2length.  The velocities of the 
muscles are then calculated using the musclevelocity functions.  From the frequencies of stimulation of the muscles, the 
activation of I2 and I1/I3 (actI2, actI1I3) are calculated.  The passive force of each muscle is calculated by giving the 
muscle's length to the function passive and then multiplying the result by the maximum force of the muscle (FnotI2 or FnotI1I3).  
Using a Hill analysis the muscle tension = Fnot * act * FV * LT.  This equation is used in the function with the passive 
tension added.  Because the muscle tension is along the muscle in direction, some geometry needed to be done in order to 
determine the total horizontal force and total vertical force that was acting on the odontophore and I1/I3 torus.  In order 
to perform this calculation, the angle of the straight piece of I2 (phi) had to first be calculated.  These derivations were 
done on 5/14/02 in Val's first notebook.  They are also in the Guidebook pg. 12-13.  The visco-elastic hinge force was 
calculated using the function calchingeforce.  Since the hinge force was characterized from the horizontal direction it could 
be added directly to the horizontal force acting on the odontophore.*/
	
	{
		double lengthI2, lengthI1I3, velocityI2, velocityI1I3, actI2, actI1I3, tensionI2, tensionI1I3;
		double passiveI2, passiveI1I3;
		double acthinge;
		//double legx, legy, hyp; BLARF

		
		//double Tphi, Bphi, A, B; BLARF
		//double vertI2force, horI2force; BLARF
        
		double Tphi, Bphi;

		double effectivePhi;
		
		//calculate the muscle lengths, are normalized
		lengthI2 = musclelengthI2(a, x, ytop, ybottom, Oangle);   //this will update phi
		lengthI2 = fixI2length (lengthI2);   //rescaling of I2 length to stay on lengthtension curve
		


		lengthI1I3 = musclelengthI1I3 (ytop, ybottom);

		// cancelling output for now *output = lengthI2;
		
		//calculating muscle velocities, are normalized
		velocityI2 = musclevelocityI2 (oldx, oldytop, oldybottom, olda, x, ytop, ybottom, a, Oangle);
		velocityI1I3 = musclevelocityI1I3 (ytop, ybottom, oldytop, oldybottom);
		
		//calculating activations
		actI2 = activationI2 (freqI2, aprimeI2);
		actI1I3 = activationI1I3 (freqI1I3, aprimeI1I3);
		acthinge = activationI2 (hingefrequency, HingeAct);
		
		//calculating the passive muscle forces
		passiveI2 = FnotI2 * passive(lengthI2);
		passiveI1I3 = 0.0; //taking out passive I1I3 forces FnotI1I3 * passive(lengthI1I3);
		
		//calculating the muscle tensions
		tensionI2 = FnotI2 * lengthtens (lengthI2) * forcevelocity (velocityI2) * actI2 + 1.0*passiveI2; //BLARF 1.8 for stable bites
		tensionI1I3 = FnotI1I3 * lengthtens (lengthI1I3) * forcevelocity (velocityI1I3) * actI1I3 + passiveI1I3/2;
		
		//update phi
		calcphi (x, ytop, ybottom, a, &Tphi, &Bphi, Oangle);
		effectivePhi = (Tphi - Bphi)/2;
		
		//if(cos(effectivePhi) < 0)
		//effectivePhi = pi/2;
		
		//based on tensions, calculating force one and force two
		//geometry done 5/14/02  --  pg 11 and 12 in Val's first lab notebook
		//pg. 12-13 in Guidebook

		//only including the elastic component visco = calchingeforce (x, xdot, F1);
		*hinge = calchingeforce2(x,xdot) + activehingeforce(acthinge, xdot, x);
		//ActiveHingeForce Here
		

		//Putting in an active hinge with a basic Zajac McLT curve

		*force1 = 2 * pi * tensionI2 * cos(effectivePhi) + *hinge + *externalforce;   ///JEFF THIS IS WHERE THE EXTERNAL FORCE GOES!!!!  

		*force2 =  2 * tensionI1I3 * pi; //+ 2*pi*tensionI2*sin(effectivePhi) taking out I2component on I3 simplemodel;
		
        //valueTBD = aprimeI2; //TATE
        
        //printf("%f \n",actI2);
        
		return (lengthI2);
	}
	
	
double calcphi (double x, double ytop, double ybottom, double a, double *Tphi, double *Bphi, double Oangle)

	/*The function calcphi takes the position and shape of the odontophore (x, y, a) and calculates the angle of the straight 
section of I2 with respect to the horizontal.  The output is in the form of the length of the legs of the right triangle that
contains phi.  The function updateforces can then used these lengths to convert muscle tension to muscle force in the desired
direction.*/
	
	{
		double x1, Ty1, By1; //point of I2 contact at back of I1/I3
		double Tx2, Bx2, Ty2, By2; //point of contact of straight section of I2 on odontophore
		double Tx1prime, Bx1prime, Ty1prime, By1prime, Tx2prime, Bx2prime, Ty2prime, By2prime;
		double b = (0.005 * sqrt(0.005)/sqrt(a)); //minor axis of the odontophore

		//double rotation = OdonAngle(x);
		double rotation = Oangle;
		double backslope = tan((90-rotation)*pi/180);
		double backxpoint = -1 * a * a * backslope/(sqrt(b*b + a*a*backslope*backslope));
		double backypoint = b * sqrt(1 - backxpoint * backxpoint/(a * a));
		double furthestx = backxpoint*cos((-1*rotation)*pi/180) + backypoint*sin((-1*rotation)*pi/180);
		
		double angle = rotation*pi/180;
		
		if (x < -1*furthestx - r)  //if I2 is attached to the odontophore 
		{
			
			
			//calculate the point where I2 connects to I1/I3
			x1 = -x - r;
			Ty1 = ytop;
			By1 = ybottom;
			
			//rotate the two points
			Tx1prime = x1*cos(angle)+Ty1*sin(angle);
			Ty1prime = Ty1*cos(angle)-x1*sin(angle);
			
			Bx1prime = x1*cos(angle)+By1*sin(angle);
			By1prime = By1*cos(angle)-x1*sin(angle);
			
			//find the (x2, y2) points in the rotated case
			if(Ty1prime < 0 && Tx1prime > 0)
				Tx2prime = (a * a * b * b * Tx1prime + sqrt(a * a * a * a * Ty1prime * Ty1prime * (a * a * Ty1prime * Ty1prime - a * a * b * b + b * b * Tx1prime * Tx1prime)))/(b * b * Tx1prime * Tx1prime + a * a * Ty1prime * Ty1prime);
			else
				Tx2prime = (a * a * b * b * Tx1prime - sqrt(a * a * a * a * Ty1prime * Ty1prime * (a * a * Ty1prime * Ty1prime - a * a * b * b + b * b * Tx1prime * Tx1prime)))/(b * b * Tx1prime * Tx1prime + a * a * Ty1prime * Ty1prime);
			
			if(Ty1prime > 0 || Tx1prime > a || Tx1prime < -1*a)
				Ty2prime = b*sqrt(1-Tx2prime*Tx2prime/(a*a));
			else
				Ty2prime = -1*b*sqrt(1-Tx2prime*Tx2prime/(a*a));
			
			
			if(By1prime < 0 && Bx1prime > 0)
				Bx2prime = (a * a * b * b * Bx1prime + sqrt(a * a * a * a * By1prime * By1prime * (a * a * By1prime * By1prime - a * a * b * b + b * b * Bx1prime * Bx1prime)))/(b * b * Bx1prime * Bx1prime + a * a * By1prime * By1prime);
			else
				Bx2prime = (a * a * b * b * Bx1prime - sqrt(a * a * a * a * By1prime * By1prime * (a * a * By1prime * By1prime - a * a * b * b + b * b * Bx1prime * Bx1prime)))/(b * b * Bx1prime * Bx1prime + a * a * By1prime * By1prime);
			
			if(By1prime > 0 || Bx1prime > a || Bx1prime < -1*a)
				By2prime = b*sqrt(1-Bx2prime*Bx2prime/(a*a));
			else
				By2prime = -1*b*sqrt(1-Bx2prime*Bx2prime/(a*a));
			
			//rotate points back
			Tx2 = Tx2prime*cos(-1*angle) + Ty2prime*sin(-1*angle);
			Ty2 = Ty2prime*cos(-1*angle) - Tx2prime*sin(-1*angle);
			
			Bx2 = Bx2prime*cos(-1*angle) + By2prime*sin(-1*angle);
			By2 = By2prime*cos(-1*angle) - Bx2prime*sin(-1*angle);
			
			//calculate phi's
			*Tphi = atan2(Ty1 - Ty2, x1 - Tx2);
			*Bphi = atan2(By1 - By2, x1 - Bx2);
			
		}
		else //if I2 is not attached to the odontophore, phi is 90 degrees
		{
			*Tphi = pi/2;
			*Bphi = -pi/2;
		}		
		
		
		return 0;
	}		


double musclelengthI2 (double a, double x, double ytop, double ybottom, double Oangle)

	/* musclelengthI2 uses the shape of the odontophore (a), the positions of the odontophore (x) and the major radius of 
I1/I3 (y) to calculate the length of I2.  The function first determines if I2 has disengaged from the odontophore.  If I2 is 
still attached to the odontophore, the length of I2 is calculated using a weighted average between the spherical equation 
and a fit to the elliptical equation.  For the spherical approximation, the length of I2 is broken into two segments: a 
straight line from the contact point of I2 on the odontophore to I2's attachment point to I1/I3 (L), and a curved segment 
that wraps around the odontophore (Radius * gamma).  Full geometry for determining the length of I2 was done on 5/14/02 and 
is found on pg.10 of Val's 1st notebook and pg 8-9 of the Guidebook.  The elliptical approximation is a polynomial fit to the 
analytic solution of I2's length for a given position and major axis length.  For major axis lengths from 5 mm to 5.6 mm, the 
two approximations are combined in a weighted average.  For major axis values larger than 5.6 mm, the elliptical fit is used 
alone.  If I2 has disengaged from the odontophore its length is equal to twice the major radius of the I1/I3 torus (y).  The 
function returns the normalized length of I2 (divided by the optimum length of I2) */

	{
		double lengthI2;
		double b = 0.005*sqrt(0.005)/sqrt(a);
		
		double rotation =  Oangle; //OdonAngle(x);
		double backslope = tan((90-rotation)*pi/180);
		double backxpoint = -1 * a * a * backslope/(sqrt(b*b + a*a*backslope*backslope));
		double backypoint = b * sqrt(1 - backxpoint * backxpoint/(a * a));
		double furthestx = backxpoint*cos((-1*rotation)*pi/180) + backypoint*sin((-1*rotation)*pi/180);
		
		
		if (x < -1*furthestx - r)  //if I2 is attached to the odontophore 
		{  
			//lengthI2 = (0.037768 - 2.15703 * x - 13.311 * a + 407.588 * x * a + 278.757 * x * x + 2595.27 * a * a + 116443* power(x, 3) - 133018. * x * x * a - 80848.4 * x * a * a - 
  			//			194671 * power(a, 3) + 8.47687*power(10,6) * power(x, 4) - 2.22748*power(10,7) * power(x,3)*a + 1.02207*power(10,7)*x*x*a*a + 6.0973*power(10,6)*x*power(a, 3) + 5.73107*power(10,6)*power(a, 4))/LI2opt;
		
			lengthI2 = (0.0402964 - 1.47509 * x - 15.4796 * a + 119.129 * x * a + 347.797 * x * x + 3193.83 * a * a + 114949 * x * x * x - 153302 * x * x * a - 41919.9 * x * a * a - 
  							260505 * a * a * a + 8.69996 * power(10, 6) * power(x, 4) - 2.17551 * power(10,7) * x * x * x * a + 1.13701 * power(10,7) * x * x * a * a + 4.29927 * power(10, 6) * x * a * a * a + 8.27669 * power(10,6) * power(a, 4))/LI2opt;
				
		}
		else  //if I2 has pulled away from the odontophore (peak protraction)
		{
			lengthI2 = (ytop - ybottom)/LI2opt;
		}
		
		
		return (lengthI2); //returns normalized lengths
	}


double fixI2length (double lengthI2)

	/*The function fixI2length renormalizes the length of I2.  The length of I2 had to be adjusted since the range of 
I2's normalized lengths fell outside the length tension curve.  This means that I2 would lose its ability to contract 
with appreciable force when it became very short or very long (peak protraction / retraction).  I2 was thus unable to 
protract the odontophore far enough or fast enough.  Full derivation on pg. 9-10 of the Guidebook.*/

	{
		double fixedlength;
		
		fixedlength = 0.5875 * lengthI2 + 0.35;
		
		if (fixedlength < 0.3)
		{
			fixedlength = 0.3;
		}
		
		return (fixedlength);	
	}			
	

double musclelengthI1I3 (double ytop, double ybottom)

	/*musclelengthI1I3 uses the major radius of the I1/I3 torus (y) to calculate the length of the I1/I3 muscle.  It does this 
using the equation for the circumference of a circle with y as the radius.  The length is then divided by the optimum I1/I3
length and the function returns the normalized length.*/

	{
		
		double radius = (ybottom + ytop)/2;
		
		if(ybottom < 0)
			radius = (-1*ybottom + ytop)/2;
		
		return ((2 * radius * pi)/LI1I3opt); //returns normalized lengths
	}		


double musclevelocityI2 (double oldx, double oldytop, double oldybottom, double olda, double x, double ytop, double ybottom, double a, double Oangle)
	
	/*musclevelocityI2 - we were unable to derive a direct equation for calculating the velocity of the I2 muscle.  Instead
the function musclevelocityI2 takes the old x, y, a and the new x, y, and a and calculates the old length of I2 and the 
new length of I2.  The function subtracts these two to get the muscle velocity.  Since both the old and new lengths of I2
are normalized, the function returns a normalized velocity.*/

	{
		double oldlength, currentlength;
		
		oldlength = musclelengthI2 (olda, oldx, oldytop, oldybottom, Oangle);
		currentlength = musclelengthI2 (a, x, ytop, ybottom, Oangle);
			
		return ((currentlength - oldlength) / StepSize);//returns normalized velocity
	}


double musclevelocityI1I3 (double ytop, double ybottom, double oldytop, double oldybottom)

	/*musclevelocityI1I3 uses the velocity of y to calculate the velocity of I1/I3.  Since the length of I1/I3 equals 
2 * pi * y, the velocity of I1/I3 is the derivative of this function (2 * pi * ydot).  Since the result is divided by
the optimal length of I1/I3, the function returns a normalized velocity.*/

	{
		double ybot, oldybot;
		
		if(ybottom > 0)
			ybot = ybottom;
		else
			ybot = -1*ybottom;
			
		if(oldybottom > 0)
			oldybot = oldybottom;
		else
			oldybot = -1*oldybottom;
		
		
		double radius = (ybot + ytop)/2;
		double oldradius = (oldybot + oldytop)/2;
		double radiusdot = (radius - oldradius)/StepSize;
		
		return (2 * pi * radiusdot / LI1I3opt);//returns normalized velocity
	}
	

double activationI2 (double frequency, double *aprime)

	/*activationI2 uses the Hill type activation function from Yu et al 1999.  The function takes the frequency of 
stimulation and calculates the amount of activation.  Some changes had to be made with the original function.  
These are explained in the Guidebook, pg. 11.*/
	
	{
		double uprime, activation, daprime;
		
		uprime = 1.03 - 4.31 * exp(-0.198 * frequency);
		
		if (uprime < uprimebound)
		{
			uprime = uprimebound;
		}
		
		if (uprime > 1)
		{
			uprime = 1;
		}	

		
		daprime = 1/tauact * (uprime - (beta + (1-beta)*uprime)* *aprime);
		*aprime = *aprime + daprime * StepSize;
		
		if (*aprime < 0)
		{
			*aprime = 0;
		}	
		
		activation = g * (*aprime - Ao);
		if (activation <0)
		{
			activation = 0;
		}
		if (activation >1)
		{
			activation = 1;
		}	
		
		return activation;
	}
	

double activationI1I3 (double frequency, double *aprime)

	/*activationI1I3 is a modified version of activationI2 created to match Hui's data (6/5/02).  The function is the
same as activationI2 except where indicated.  Changes are explained in the Guidebook pg. 11.*/

	{  	
		double uprime, activation, daprime;
		
					//the 4.31 was replaced with 1.0 to reflect Hui's data
		uprime = 1.03 - 1.0 * exp(-0.198 * frequency);
		
		//Changing uprime in order to reflect a faster relaxation time of I1I3
		if (frequency<5.0)
			{							//the 6.31 was 3.31
				uprime = 1.03 - (1.0 - (6.31 * (frequency - 5)/5)) * exp(-0.198 * frequency);
			}
		
		if (uprime < uprimebound)
		{
			uprime = uprimebound;
		}
		
		if (uprime > 1)
		{
			uprime = 1;
		}	
		
		if (frequency < 5.0)  //added to give better results from the neural recordings
		{
			uprime = -1.0;
		}	
			
		
		daprime = 1/tauact * (uprime - (beta + (1-beta)*uprime)* *aprime);
		*aprime = *aprime + daprime * StepSize;
		
		if (*aprime < 0)
		{
			*aprime = 0;
		}	
		
		activation = g * (*aprime - Ao);
		
		if (activation <0)
		{
			activation = 0;
		}
		if (activation >1)
		{
			activation = 1;
		}	
		
		return activation;
	}


double activationN3(double frequency, double *act, double time)
	{
		double A, B, actprime = 0;
		
		A = frequency * 0.386779 - 2.93698;
		
		B = (A * time + log(1/ *act - 1))/A;
		
		actprime = (A * exp(-1*A*(time - B)))/pow((1 + exp(-1*A*(time - B))),2);
		
		*act = *act + actprime*StepSize;
		
		if(*act > 0.9789)
			*act = 0.9789;
		if(*act < 0.01)
			*act = 0.01;
			
		return (*act*1.0321-0.010321);
	
	}
double secondactivationN3(double frequency, double *act, double time)
	{
		double A, B, actprime = 0;
		
		//A = frequency * 0.386779 - 2.93698;
		A = 0.6 * frequency - 7.0;  //for october 2008 simulations
		//A = 4.8 * frequency - 56.0;
		A = 0.9 * frequency - 14.0;

		//B = (A * time + log(1/ *act - 1))/A;
		
		//actprime = (A * exp(-1*A*(time - B)))/pow((1 + exp(-1*A*(time - B))),2);

		B = (log(1/ *act - 1))/A;
		
		actprime = (A * exp(-1*A*( - B)))/pow((1 + exp(-1*A*(- B))),2);
		
		*act = *act + actprime*StepSize;
		
		if(*act > 0.9789)
			*act = 0.9789;

		
		if(*act < 0.01)
			*act = 0.01;
		

		return (*act*1.0321-0.010321);
	
	}	

double passive (double musclength)

	/*passive uses the equation for the amount of passive tension given in Yu et al 1999.  The function takes muscle length
and calculates the passive tension.  Since the passive tension cannot be negative, the function returns zero if the passive
tension is less that zero.*/

	{
		double passtension;
		passtension = -0.41 + 1.04 * exp((musclength - 1.18) / 0.28);
		
		if (passtension < 0)
		{
			return 0;
		}
		else
		{
			return passtension; 
		}	
	}		


double lengthtens (double length)

	/*lengthtens uses the Hill type length tension function from Yu et al 1999.  The function uses the length of a muscle
to calculate the amount of length tension.  This can be calculated as though a tendon was present or without a tendon.*/
	
	{
		double output;
		
		//output = 1.81 - 10.8 * length + 19.2 * length * length - 9.2 * length * length * length; complexmod//length tension without a tendon
		output = (-5.27*length*length + 10.54*length-4.27); //length tension with a tendon simpmodel
		
		if (output > 0.0) //lengthtension property cannot be < 0 (function returns 0 if it is)
		{
			return (output);
		}
		else
		{
			return (0.0);
		}
	}


double forcevelocity (double musclevelocity)

	/*forcevelocity uses the Hill type force velocity function from Yu et al 1999.  The function uses the velocity of a
muscle to calculate the force velocity term.  The velocity that is handed to the function must have its sign reversed
because in the function from the paper assumes that lengthening is negative and shortening is positive.*/

	{
		double musclevel;
		musclevel = -1 * musclevelocity; //changing the sign because the forcevelocity 
										//function assumes that lengthening is negative
		/*if (musclevel < 0) //lengthening
		{
			return (1.0 + 0.61 / (1.0 - 0.038 / musclevel));
		}
		else //shortening
		{
			return (1.0 / (1.0 + 10.8 * musclevel));
		}*/
		return 0.33;  //Simplemodel
	}


double calchingeforce (double xdisplacement, double	xvelocity, double *F1)

	/*calchingforce uses the function for the hinge force from Sutton et al 2002.  Explanation of passiveoffset: without
an offset the hinge forces on the odontophore were too weak to retract the odontophore back to a point where x became 
negative again.  The biology tells us that the slug is able to retract its odontophore.  To model this, the passiveoffset
parameter was added to the calchingeforce function so that the passive force curve was shifted to allow the hinge forces 
to be stronger.  The hinge function is not completely contained within calchingeforce.  The hinge function also has to keep 
track of the velocity of the odontophore (xdot).  If the odontophore changes direction (xdot switches signs) then the value 
of F1 has to be reset to zero so that the visco-elastic force does not overwhelm the elastic hinge force over large amounts 
of time.  This operation is performed in the model code rather than in calchingeforce.
	This model accepts displacements in M, and velocities in M/s, and then does
	the conversion to mm and mm/s internally, and returns all the forces in mN
	this is what we'll work with until we get it all working, and then will fix
	it so the conversions are external*/

	{
		double Fo, So, S1, Damp1, F1dot, x, xdot;
		
		x = (xdisplacement) * 1000; //convert to mm
		xdot = xvelocity * 1000; //convert to mm
		
		//check whether we are using the protraction or retraction curve
	if (xdot >= 0.0)	//protraction curve  
		{
			if ((x + passiveoffset) < xp) //if odontophore is not protracted past xp
			{
				So = 0;  //So is the spring stiffness of the Kelvin element
			}		
			else  //if the odontophore has protracted past xp 
			{
				So = A2 + (A1-A2)/(1+exp(((x + passiveoffset)-Xo)/DX));		
			}
		}
			
		else //retraction curve
		{
			if ((x + passiveoffset) < xr)
			{
				So = 0;
			}
			else
			{
				So = B1 + B2 * (x + passiveoffset) + B3 * (x + passiveoffset) * (x + passiveoffset) + B4 * (x + passiveoffset) * (x + passiveoffset) * (x + passiveoffset); ;
			}
		}  
		
		Fo = So * (x + passiveoffset); //force from the spring in the Kelvin element
		
		if (Fo < 0)
		{
			Fo = 0;
		}
		if (xdot >= 0.0)  //protraction curve
			{
				S1 = C1 * (x + passiveoffset)*(x + passiveoffset) + C2*(x + passiveoffset) + C3; //spring stiffness of the Maxwell element

				Damp1 = E1*(x + passiveoffset)*(x + passiveoffset) + E2 * (x + passiveoffset) + E3; //dashpot viscocity of the Maxwell element
				
			if (S1 < 0)
				{
					S1 = 0;
				}
			}
		else  //retraction curve
			{
				 S1 = -1.5*(D1 * (x + passiveoffset) + D2);
				
				 Damp1 = 3.0 * (AF1 * (x + passiveoffset) + AF2);
				
				/*if(S1>0)
				{
					S1 = 0;
				}*/
			}
			
			//equation of motion for Maxwell element
			F1dot = xdot * S1 - (S1 / Damp1) * *F1;
			*F1 = *F1 + F1dot * StepSize;
			
			if (xdot < 0) //prevents the negative spring constant from causing trouble
			{
				if (*F1 < (-6.5594 - 857.325 * xdisplacement))
				{
					*F1 = -6.5594 - 857.325 * xdisplacement;
				}
			}
			
			//Last minute checks of what the force is to prevent 
			//odontophores from flinging to infinity and BEYOND!
			if (Fo < 0)
			{
				Fo = 0;
			}
			if ((Fo + *F1) < 0)
			{
				*F1 =  -Fo;
			}
			
		/* Making viscoelastic force = 0  Simpmodel */
		*F1 = 0;
        
		return (-1 * Fo - *F1 - viscdamp * xdot);

	}
	
double updateposition (double *x, double *ytop, double *ybottom, double *a, double *xdot, double *ydot, double *adot, double f1, double f2, double actN3, double *Txc, double *Bxc, double *output, double *Oangle, double hingeforce)

	/*updateposition takes the previous positions of the odontophore (x), the major radius of I1/I3 (y), the major axis of the 
odontophore (a), the velocity of the odontophore (xdot), the velocity of I1/I3 (ydot), the horizontal force on the odontophore
(force1 in model, f1 in updateposition), and the vertical force on the I1/I3 (force2 in model, f2 in updateposition).  Three 
different equations can be used (different equations are selected by changing the value of eom at the beginning of the code.
The code first updates the value of the major axis using the frequency of N3 and a linear transform.  It then calculates adot
and addot using a numerical approximation.  Then the value of the contact point between the odontophore and I1/I3 (xc, yc) and
the slope at that contact point (slope) is updated.  xddot is then calculated using one of the equations of motion.
	The first equation is the complex equation of motion for a shape changing odontophore.  The full derivation is shown in 
the Guidebook pg. 21-33.  The equation used is a fit of the full analytic equation which was too complex to use in the model 
code.  In this equation, all lengths are converted to mm and forces are converted to microN.  The equation takes the shape and 
position of the odontophore (a, x) and calculates the coefficients of the xdot^2, adot^2, addot, f1, f2, and adot*xdot terms 
in the equation of motion.  Then using the velocity of the odontophore (xdot), the velocity and acceleration of the major axis 
(adot, addot) and the horizontal and vertical forces on the odontophore (f1, f2) the equation calculates the acceleration of 
the odontophore.
	The second equation is the simple equation of motion (SEOM).  The simple equation of motion is based on the assumption 
that the model is quasi-static and therefore, the contact constraint, velocity terms, and mass of the I1/I3 ring do not 
contribute a significant amount to the acceleration of the odontophore.
	The third equation is for an almost quasi-static model.  The derivation of this equation is given in the Guidebook 
pg. 33-34.  Once xddot is calculated, the function then integrates xddot using an Euler integrator to find the new velocity 
of the odontophore (xdot).  xdot is integrated using an Euler integrator to get the new position of the odontophore (x). The 
major radius of I1/I3 (y) is calculated using the contact point between I1/I3 and the odontophore (derivation in Guidebook 
pg. 33).  All of the updated variables are given to the function as pointers and are thus updated in the model function as 
they are changed in updateposition.  The function does not need to return anything.*/
	
	{
		//eom variables - explained where initialized
		//double xc, yc, slope; blarf
		double xddot, oldy;
		double Tyc, Byc, topslope, bottomslope;
		double effectiveslope;
		
		*Oangle = OdonAngle2(*x, hingeforce, *a, *Oangle);
			
		oldy = *ytop; //oldy is updated before y is changed (used to find ydot numerically)
		
		//xcCalc(*a, -1* *x, OdonAngle(*x), *Txc, *Bxc, Tyc, Byc, topslope, bottomslope);
        xcCalc(*a, -1* *x, *Oangle, *Txc, *Bxc, Tyc, Byc, topslope, bottomslope);
		//contactslope (a, x, &yc, &xc, &slope); //initializing contact point and slope at that point
		
		effectiveslope = (topslope - bottomslope)/2;
		*output = effectiveslope; // output back in removed output simplemodel
		
		//slightly more complex simple eom (NSEOM)
		xddot =(f1 + f2*effectiveslope)/(Odonmass+I1I3mass*effectiveslope*effectiveslope);

		//*a = freqN3 * 0.00015 + 0.005; //major axis length is updated from freqN3
		*a = 0.003*actN3+0.005; //removing ability to change a, simplemodel this was original function above was commented out	
		
		//BLARF  making aspect ratio maximum at 0.007
		/*if (*a > 0.006)
		{
			*a = .006;
		} */

		//now we integrate to get x
		*xdot = *xdot + xddot * StepSize; //get xdot using integration of xddot

		//BLARF putting in a maximum retraction velocity to try to get stability
		if (*xdot < -.1)
		{
			*xdot = -0.1;
		}
		*x = *x + *xdot * StepSize; //get x with integration


		
		//use geometry to find y and a numerical approximation for ydot
		*ytop = Tyc + sqrt(r * r - (-1 * *x - *Txc) * (-1 * *x - *Txc)); 
		*ybottom = Byc - sqrt(r*r - (-1* *x - *Bxc)*(-1* *x - *Bxc));
		
		*ydot = (*ytop - oldy)/StepSize;
		
		return 0;
	}


double contactslope (double *a, double *x, double *yc, double *xc, double *slope, double Oangle)

	/*Updates the contact point between the odontophore and the ellipse (xc, yc) and calculates the slope of the tangent
to that point.*/
	
	{
		double b;
		double amm,xmm;
		
		amm = *a * 1000;  //converts a to mm
		xmm = *x * -1000; //converts x to mm
		
		
		b = (5.0 * sqrt(5.0)/ sqrt(*a * 1000.0))/1000.0; //calculates the minor axis of the odontophore
		
		//calculates xc using the polynomial fit (xc derivation.nb and Guidebook pg. 22-23)
		*xc = (0.0350629 - 0.0198081* amm + 0.00438099* amm*amm - 
  				0.00045905* power(amm,3) + 0.0000191278 *power(amm,4) + 0.74469*xmm - 
 			 	0.030564 *amm* xmm + 0.0102451* power(amm,2)* xmm - 0.00036768 *power(amm,3)*xmm + 
 			 	0.000153275* xmm*xmm - 0.0000534595 *amm *xmm*xmm + 0.0000045673*power(amm,2)*power(xmm,2) + 
 			 	0.002039379 *power(xmm,3) - 0.00041629 *amm*power(xmm,3) - 0.0000000139376555* power(xmm,4))/1000;
		
		*yc = b * sqrt(1 - *xc * *xc/(*a * *a)); //yc found using the equation for an ellipse
		
		if (*yc < 0)  //check to make sure yc doesn't go negative
		{
			*yc = -1 * *yc;
		}	
		
		//slope at the contact point found using the derivative of the ellipse equation
		*slope = (-1 * b * *xc)/(*a * *a * sqrt(1 - (*xc * *xc)/(*a * *a))); 
		
		
		return 0;	
	}	
	// putting in fitcalc for the simulations for the model model sweeps

void fitcalc (double x, double a, double oldx, double olda, double & fitness)
{
	double xdist, ydist;
	if (a>.007)  //I've put in a threshhold on the fitness and seeing what it does to rejection.
	{
	// old function no Y fitness += (a*cos(OdonAngle(x)*pi/180)+x - olda*cos(OdonAngle(oldx)*pi/180)-oldx)*((a-0.005)/.003);

		xdist =  (a*cos(OdonAngle(x)*pi/180)+x - olda*cos(OdonAngle(oldx)*pi/180)-oldx)*((a-0.005)/.003);
	    ydist =  (a*sin(OdonAngle(x)*pi/180) - olda*sin(OdonAngle(oldx)*pi/180))*((a-0.005)/.003);

		fitness += -ydist + xdist;
	}
	
	//fitness += (x - oldx)*(a/0.003 - 5.0/3.0);
}

void fitcalcbite (double x, double a, double oldx, double olda, double & fitness)
	//updates an individual's fitness
{
	if (x>fitness)
	{
		fitness = x;
	}
}

void fitcalcswallow (double x, double a, double oldx, double olda, double & fitness)
{
	
}
	//updates an individual's fitness

void fitcalcreject (double x, double a, double oldx, double olda, double & fitness)
{
	
}
	//updates an individual's fitness
			

	
double power(double num, double ex)

	/* the function power takes the value of num and raises it to the "ex" power.  This is used to shorten the polynomial fit
equations.*/

	{
		int counter;
		double output;
		
		counter = 1;
		output = num;
		
		
		while(counter < ex)
		{
			output = output*num;
			
			counter = counter+1;
		}
		
		return output;
	}		


double topEllipseSlope(double a, double xc, double rotation, double & yc)
{
	double theta = rotation*pi/180;
	double b = 0.005*sqrt(0.005)/sqrt(a);
	double part1 = a*a*cos(theta)*xc;
	double part2 = power(a,4)*b*b*cos(theta)*cos(theta)*sin(theta)*sin(theta) - a*a*b*b*xc*xc*sin(theta)*sin(theta)+a*a*power(b,4)*power(sin(theta),4);
	double part3 = a*a*cos(theta)*cos(theta)+b*b*sin(theta)*sin(theta);
	double xcHoriz = (part1+sqrt(part2))/part3;
	double horizslope = -b*xcHoriz/(a*a*sqrt(1-xcHoriz*xcHoriz/(a*a)));
	
	yc = -1*sin(-1*theta)*xcHoriz + cos(-1*theta)*b*sqrt(1-xcHoriz*xcHoriz/(a*a));
	
	return ((-sin(-1*theta) + horizslope * cos(-1*theta))/(cos(-1*theta) + horizslope*sin(-1*theta)));
}

double topEllipseSlope2(double a, double xc, double rotation, double & yc)
{
	double theta = rotation*pi/180;
	double b = 0.005*sqrt(0.005)/sqrt(a);
	double part1 = a*a*cos(theta)*xc;
	double part2 = power(a,4)*b*b*cos(theta)*cos(theta)*sin(theta)*sin(theta) - a*a*b*b*xc*xc*sin(theta)*sin(theta)+a*a*power(b,4)*power(sin(theta),4);
	double part3 = a*a*cos(theta)*cos(theta)+b*b*sin(theta)*sin(theta);
	double xcHoriz = (part1+sqrt(part2))/part3;
	double horizslope = b*xcHoriz/(a*a*sqrt(1-xcHoriz*xcHoriz/(a*a)));
	
	yc = -1*sin(-1*theta)*xcHoriz + cos(-1*theta)*-1*b*sqrt(1-xcHoriz*xcHoriz/(a*a));
	
	return ((-sin(-1*theta) + horizslope * cos(-1*theta))/(cos(-1*theta) + horizslope*sin(-1*theta)));
}

double bottomEllipseSlope(double a, double xc, double rotation, double & yc)
{
	double theta = rotation*pi/180;
	double b = 0.005*sqrt(0.005)/sqrt(a);
	double part1 = a*a*cos(theta)*xc;
	double part2 = power(a,4)*b*b*cos(theta)*cos(theta)*sin(theta)*sin(theta) - a*a*b*b*xc*xc*sin(theta)*sin(theta)+a*a*power(b,4)*power(sin(theta),4);
	double part3 = a*a*cos(theta)*cos(theta)+b*b*sin(theta)*sin(theta);
	double xcHoriz = (part1 - sqrt(part2))/part3;
	double horizslope = b*xcHoriz/(a*a*sqrt(1-xcHoriz*xcHoriz/(a*a)));
	
	yc = -1*sin(-1*theta)*xcHoriz + cos(-1*theta)*-1*b*sqrt(1-xcHoriz*xcHoriz/(a*a));
	
	return ((-sin(-1*theta) + horizslope * cos(-1*theta))/(cos(-1*theta) + horizslope*sin(-1*theta)));
}

double topI1I3slope(double x, double xc)
/*slope of top I1I3 ring (bottom circle curve*/
{
	return ((xc - x)/sqrt(r*r - (xc - x)*(xc - x)));
}

double bottomI1I3slope(double x, double xc)
/*slope of bottom I1I3 ring (top circle curve)*/

{
	return (-1*(xc - x)/sqrt(r*r - (xc - x)*(xc - x)));
}


void xcCalc(double a, double x, double rotation, double & Txc, double & Bxc, double & Tyc, double & Byc, double & topslope, double & bottomslope)
{
	double TES, BES, TI1I3S, BI1I3S, n, ndelta, Tslope;
	double b, odonslope, frontxpoint, frontypoint, furthestx;
	double oldTxc = Txc, oldBxc = Bxc;
	int counter = 0;
	
	b = 0.005*sqrt(0.005)/sqrt(a);
	odonslope = tan((90 - rotation)*pi/180);
	frontxpoint = a*a*odonslope/sqrt(b*b+a*a*odonslope*odonslope);
	frontypoint = -1*b*sqrt(1-frontxpoint*frontxpoint/(a*a));
	furthestx = frontxpoint*cos((-1*rotation)*pi/180) + frontypoint*sin((-1*rotation)*pi/180);	

	
	if(Bxc < (-r + x))
		Bxc = x - 0.00125 + 0.0000101;
	if(Bxc > (r + x))
		Bxc = x + r - 0.0000101;
	if(Txc < (x - r))
		Txc = x - r + 0.0000101;
	if(Txc > (r + x))
		Txc = x + r - 0.0000101;

	if(Txc > furthestx)
		Txc = (furthestx + x - r)/2;
	if(Txc < -1*furthestx)
		Txc = (-1*furthestx + x + r)/2;
	if(Bxc > furthestx)
		Bxc = (furthestx + x - r)/2;
	if(Bxc < -1*furthestx)
		Bxc = (-1*furthestx + x + r)/2; 
	
	for(;;)
	{
		if(rotation > 45 && (Txc > (a * cos(rotation*pi/180))))
			TES = topEllipseSlope2(a, Txc, rotation, Tyc);
		else
			TES = topEllipseSlope(a, Txc, rotation, Tyc);
			
		TI1I3S = topI1I3slope(x, Txc);
		
		n = TES - TI1I3S;
		
		topslope = TES;
		
		if((n < 0 && -1*n < Error) || (n > 0 && n < Error))
			break;
		
		counter = counter + 1;

		if(counter > 10)
		{
			Txc = oldTxc;
			
			if(rotation > 45 && (Txc > (a * cos(rotation*pi/180))))
				topslope = topEllipseSlope2(a, Txc, rotation, Tyc);
			else
				topslope = topEllipseSlope(a, Txc, rotation, Tyc);
			
			break;
		}
			
		if(rotation > 45 && ((Txc+delta) > (a * cos(rotation*pi/180))))
			Tslope = topEllipseSlope2(a, Txc+delta, rotation, Tyc);
		else
			Tslope = topEllipseSlope(a, Txc+delta, rotation, Tyc);
		
		ndelta = Tslope - topI1I3slope(x, Txc + delta);
		
		if((Txc - n*delta/(ndelta - n)) > x - 0.00125)
			Txc -= n*delta/(ndelta - n);
		else
			Txc = x - r + 0.0000001;
	}
	
	counter = 0;

	for(;;)
	{
		BES = bottomEllipseSlope(a, Bxc, rotation, Byc);
		BI1I3S = bottomI1I3slope(x, Bxc);
		
		n = BES - BI1I3S;
		
		bottomslope = BES;
		
		if((n < 0 && -1*n < Error) || (n > 0 && n < Error))
			break;
		
		counter=counter+1;

		if(counter > 10)
		{	
			Bxc = oldBxc;
			bottomslope = bottomEllipseSlope(a, Bxc, rotation, Byc);
			break;		
		}
		
		ndelta = bottomEllipseSlope(a, Bxc + delta, rotation, Byc) - bottomI1I3slope(x, Bxc + delta);
		
		if((Bxc - (n*delta)/(ndelta - n)) > (-0.00125 + x))
			Bxc -= n*delta/(ndelta - n);
		else
			Bxc = x - r + 0.0000001;
	}
	
}

/*double OdonAngle(double x)
{
	if(x < -0.00518916)
		return 90;
	else
	{
		if(x < -0.000349894)
			return -11905*x+28.223;
		else
		{
			if(x < -0.000146141)
				return -99083*x-2.2801;
			else
				return 12.2;
		} Val's old function
	}*/

double OdonAngle(double x)  //my new function for rejection, based on Hui's paper, Figure 5E,F
{
	if(x < -0.0011)
		return 90;
	else
	{
		if ( x>= .0051)
			return 0;
			
		else
		{
			if ((90 - 18000*(x+.0011)) > 0)
			{
				return (90 - 18000*(x+.0011));
			}
		
			 else
			{
				return 0;
			}
			;  //new rot function
			//return -14516*x + 74.03;  //old rot function
		}
	}

	
	
	/*if(x < -0.00517083)
		return 90;
	else
	{
		if(x < 0.000345961)
			return -10443*x+36.001;
		else
		{
			if(x < 0.00057778)
				return -86914*x+62.457;
			else
				return 12.2;
		}
	}*/
}
double calchingeforce2 (double xdisplacement, double xvelocity)
/* calchingeforce2 only includes the elastic anti-protraction component of the hinge force */

	{
		double Fo, So, x, xdot;
		
		x = (xdisplacement) * 1000; //convert to mm
		xdot = xvelocity*1000;

		/* Klugy way of only doing the protraction curve */
		/*if(xvelocity>0)
		{
		xdot = xvelocity * 1000;
		}
		else
		{
			xdot = -xvelocity*1000;
		}*/
		//convert to mm
		
		//check whether we are using the protraction or retraction curve
	if (xdot >= 0.0)	//protraction curve  
		{
			if ((x + passiveoffset) < xp) //if odontophore is not protracted past xp
			{
				So = 0;  //So is the spring stiffness of the Kelvin element
			}		
			else  //if the odontophore has protracted past xp 
			{
				So = A2 + (A1-A2)/(1+exp(((x + passiveoffset)-Xo)/DX));		
			}
		}
			
		else //retraction curve
		{
			if ((x + passiveoffset) < xp) //if odontophore is not protracted past xp
			{
				So = 0;  //So is the spring stiffness of the Kelvin element
			}		
			else  //if the odontophore has protracted past xp 
			{
				So = .95*(A2 + (A1-A2)/(1+exp(((x + passiveoffset)-Xo)/DX)));		
			}


			/*if ((x + passiveoffset) < xr)
			{
				So = 0;
			}
			else
			{
				So = B1 + B2 * (x + passiveoffset) + B3 * (x + passiveoffset) * (x + passiveoffset) + B4 * (x + passiveoffset) * (x + passiveoffset) * (x + passiveoffset); ;
			}*/
		}  
		
		Fo = So * (x + passiveoffset) + viscdamp * xdot; //force from the spring in the Kelvin element putting in some damping
		
		if (Fo < 0)
		{
			Fo = 0;
		}
			
		

		return (-1 * Fo  - viscdamp * xdot);
	}

double activehingeforce (double activation, double velocity, double length)
{
    
	double HingeLT;
	double FoHinge =  -600;  // -800;

	HingeLT = -.024*((length+lengthshift)*1000)*((length+lengthshift)*1000) - .04*((length+lengthshift)*1000) +1;
    

	if ((HingeLT > 0) && (activation > Ao))
		return FoHinge*HingeLT*(activation-Ao);
	else
		return 0;
}
double OdonAngle2(double x, double hingeforce, double radius, double oldodonangle)
{
	double output;
	double forceadjustment, rotationadjustment, totaladjustment, maxangleoffset;
	double equilibriumangle;

	if(x < -0.0011)
		output = 90;
	else
	{
		if ( x>= .0051)
			output = 0;
			
		else
		{
			if ((90 - 18000*(x+.0011)) > 0)
			{
				output = (90 - 18000*(x+.0011));
			}
			else
			{
				output = 0;
			}
			;  //new rot function
			//return -14516*x + 74.03;  //old rot function
		}
	}

	equilibriumangle = output;

	//maxangleoffset = output;
    maxangleoffset = 90 - output;

	//putting in the dynamic offset as a function of hinge force (which is negative) and radius
    if (hingeforce < -75)
	{
		rotationadjustment = (radius-.005)/ .003;
		forceadjustment = (-1*hingeforce - 75)/75;
		if (forceadjustment > 1.0)
		{
			forceadjustment = 1.0;
		}
		totaladjustment = 2*(maxangleoffset * rotationadjustment*forceadjustment);

		// putting in a maximum rotation speed to prevent massive rotation. Max speed, 
		equilibriumangle = output - totaladjustment;
		if(equilibriumangle <0)
		{
			equilibriumangle=0;
		}
	} 

	output = oldodonangle + anglestiffness*(equilibriumangle - oldodonangle);
	
	return output;
}

/* Changes the neural variables and seaweed force variable in place using references based on the current time
 */

void updateinputs (double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2, string behaviorType)
{
    
    if(behaviorType == "RejectionB")
    {
        updateinputsRejectionB(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else if(behaviorType == "Bite")
    {
        updateinputsBite(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else if(behaviorType == "SwallowB")
    {
        updateinputsSwallowB(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);

    }
    else if(behaviorType == "SwallowA")
    {
        updateinputsSwallowA(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else if(behaviorType == "RejectionA")
    {
        updateinputsRejectionA(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else if(behaviorType == "SwallowPerturbed")
    {
        updateinputsSwallowPerturbed(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else if (behaviorType == "Dynamic")
    {
        updateDynamicInputs(time, freqI2, freqHinge, freqI1I3, freqN3, seaweedforce, a, frequencyiterationtime, frequencyiterationtime2);
    }
    else
    {
        cerr << "Behavior Type Not Recognized"  << endl;
        exit(1);
    }
}

void updateinputsRejectionB(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    //Rejection B sqaure wave inputs
    if (time > 0.38)
     {
     freqI2 = 20;
     }
     
     if (time > 2.75)
     {
     freqI2 = 0;
     }
     
     if (time > 1.01)
     {
     freqHinge = 20;  //BLARF was 20
     }
     
     if (time> 7.56)
     {
     freqHinge = 0;
     }
     
     if (time> 3.56)
     {
     freqI1I3 = 20;   //BLARF was 35
     }
     
     if (time > 7.69)
     {
     freqI1I3 = 0;
     }
     
     if (time > .8)
     {
     freqN3 = 30;  //BLARF was 20
     }
     
     if (time > 2.5)
     {
     freqN3 = 0;
     }
}

void updateinputsBite(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    // Bite square wave inputs
    if (time > 0.0)
     {
     freqI2 = 20;
     }
     
     if (time > 3.41)
     {
     freqI2 = 0;
     }
     
     if (time > 2.36)
     {
     freqHinge = 20;  //BLARF was 20
     }
     
     if (time> 6.80)
     {
     freqHinge = 0;
     }
     
     if (time> 2.21)
     {
     freqI1I3 = 20;   //BLARF was 35
     }
     
     if (time > 6.56)
     {
     freqI1I3 = 0;
     }
     
     if (time > 2.15)
     {
     freqN3 = 30;  //BLARF was 20
     }
     
     if (time > 4.85)
     {
     freqN3 = 0;
     }
}

void updateinputsSwallowB(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    //Proposed swallow B code
    if (time > 0.0)
     {
     freqI2 = 20;
     }
     
     if (time > 2.05)
     {
     freqI2 = 0;
     }
     
     if (time > 2.48)
     {
     freqHinge = 20;
     }
     
     if (time> 6.62)
     {
     freqHinge = 0;
     }
     
     if (time> 2.33)
     {
     freqI1I3 = 20;
     }
     
     if (time > 6.62)
     {
     freqI1I3 = 0;
     }
     
     if (time > frequencyiterationtime) //(time > 1.4)
     {
     freqN3 = 30;
     }
     
     if (time > frequencyiterationtime + 2.85) //(time > 4.25)
     {
     freqN3 = 0;
     }  //End of proposed Swallow B
}

void updateinputsSwallowA(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    // Beginning of proposed swallow A
    if (time > 0.0)
     {
     freqI2 = 20;
     }
     
     if (time > 1.51)
     {
     freqI2 = 0;
     }
     
     if (time > 1.8)
     {
     freqHinge = 20;
     }
     
     if (time> 5.8)
     {
     freqHinge = 00;
     }
     
     if (time> 1.95)
     {
     freqI1I3 = 20;
     }
     
     if (time > 6.3)
     {
     freqI1I3 = 0;
     }
     
     if (time > 1.0)
     {
     freqN3 = 30;
     }
     
     if (time > 3.9)
     {
     freqN3 = 0;
     }    //End of proposed swallow A
}

void updateinputsRejectionA(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    //this is rejection A
    if (time >.6)
     {
     freqN3 = 30;
     }
     
     if (time > 1.7)
     {
     freqN3 = 0;
     }
     
     if (time > .26)
     {
     freqI2 = 20;
     }
     
     if (time > 1.9)
     {
     freqI2 = 0;
     }
     
     if (time>2.5)
     {
     freqI1I3 =  20;
     }
     if (time>6.6)
     {
     freqI1I3 = 0;
     }
     
     if (time>1.2)
     {
     freqHinge = 20;
     }
     if (time> 6.0)
     {
     freqHinge = 0;
     }
}

void updateinputsSwallowPerturbed(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    //this is is Hillel's perturbed swallow
    // TODO come back to this:
    {
        freqI2 = 20;
    }
    
    if (time > 2.05)
    {
        freqI2 = 0;
    }
    
    if (time > 2.48)
    {
        freqHinge = 20;
    }
    
    if (time> 6.62)
    {
        freqHinge = 0;
    }
    
    if (time > 3.5)  //time>2.55
    {
        freqI1I3 = 20;
    }
    
    if (time > 6.62)
    {
        freqI1I3 = 0;
    }
    
    if (time > 1.4) //(time > 1.4)
    {
        freqN3 = 30;
    }
    
    if (time > 3.0) //(time > 4.25)
    {
        freqN3 = 0;
    }
    
    if (time > 1.5)
    {
        
        if (time < 2.5)
        {
            seaweedforce = MAXSEAWEEDFORCE * ((a - .005)/.003)*(time - 1.5);
        }
        
        else
        {
            seaweedforce =  MAXSEAWEEDFORCE*((a - .005)/.003);
        }
        
        
    }
}

bool openAndRead(string behaviorType)
{
    if(behaviorType == "Dynamic"){
        if(numberColumns == 15){
            inputFile.open("SlugInput2.txt");
            //If the File isn't openable/found then it will fail
            if(inputFile.fail()){
                //cout << "File Not Found  ";
            } //Talk to jeff... Not sure why but the file fails to open every time, but still is read and works as it should..?
            int count = 0; //while loop counter
            while ((inputFile >> times >> positions >> radiuss >> angles >> hingeFs >> fitnesss >> freqI2s >> freqI1I3s >> freqN3s >> freqHinges >>actI2s >> actI1I3s >> acthinges >> fitnesss >> seaweedforces))
            {
                if(count > 0){ //first (0th) line of file is text
                    timearray[count-1] = std::stod(times);
                    positionarray[count-1] = std::stod(positions);
                    radiusarray[count-1] = std::stod(radiuss);
                    anglearray[count-1] = std::stod(angles);
                    hingeFarray[count-1] = std::stod(hingeFs);
                    fitnessarray[count-1] = std::stod(fitnesss);
                    freqI2array[count-1] = std::stod(freqI2s);
                    freqI1I3array[count-1] = std::stod(freqI1I3s);
                    freqN3array[count-1] = std::stod(freqN3s);
                    freqHingearray[count-1] = std::stod(freqHinges);
                    actI2array[count-1] = std::stod(actI2s);
                    actI1I3array[count-1] = std::stod(actI1I3s);
                    acthingearray[count-1] = std::stod(acthinges);
                    seaweedforcearray[count-1] = std::stod(seaweedforces);
                }
                count++;
            }
            return true;
        }
        else if (numberColumns == 6){
            inputFile.open("SlugInput2.txt");
            //If the File isn't openable/found then it will fail
            if(inputFile.fail()){
                //cout << "File Not Found  ";
            } //Talk to jeff... Not sure why but the file fails to open every time, but still is read and works as it should..?
            int count = 0; //while loop counter
            while ((inputFile >> times >> freqI2s >> freqI1I3s >> freqN3s >> freqHinges >> seaweedforces))
            {
                if(count > 0){ //first (0th) line of file is text
                    timearray[count-1] = std::stod(times);
                    freqI2array[count-1] = std::stod(freqI2s);
                    freqI1I3array[count-1] = std::stod(freqI1I3s);
                    freqN3array[count-1] = std::stod(freqN3s);
                    freqHingearray[count-1] = std::stod(freqHinges);
                    seaweedforcearray[count-1] = std::stod(seaweedforces);
                }
                count++;
            }
            return true;
        }
        else if (numberColumns == 5){
            inputFile.open("SlugInput2.txt");
            //If the File isn't openable/found then it will fail
            if(inputFile.fail()){
                //cout << "File Not Found  ";
            } //Talk to jeff... Not sure why but the file fails to open every time, but still is read and works as it should..?
            int count = 0; //while loop counter
            while ((inputFile >> times >> freqI2s >> freqI1I3s >> freqN3s >> freqHinges))
            {
                if(count > 0){ //first (0th) line of file is text
                    timearray[count-1] = std::stod(times);
                    freqI2array[count-1] = std::stod(freqI2s);
                    freqI1I3array[count-1] = std::stod(freqI1I3s);
                    freqN3array[count-1] = std::stod(freqN3s);
                    freqHingearray[count-1] = std::stod(freqHinges);
                }
                count++;
            }
            return true;
        }
        else{
            return false; //The entire program will not run if the behavior type is "dynamic" and there is no correct input
        }
    }
    else {
        return true;
    }
}

/* timeAdjuster takes an double array "timeArray" and double variable "timeStamp" as inputs. It's purpose is to determine
 what value within an array is closest to the input timeStamp, and output the array index at which that value is located
 within the array */
int timeAdjuster(double timeArray[], double timeStamp)
{
    int count = 0;
    if(numberColumns == 5){
        while(timeStamp > timeArray[count])
        {
            count++;
        }
    }
    if(numberColumns == 6){
        while(timeStamp >= timeArray[count])
        {
            count++;
        }
    }
    return count - 1; //Fixes cases where input data does not start at zero seconds
}

/*
 Updates the inputs of the model to equal the values from the input file at a time. For seaweed force dependent models, input files must include two values under the seaweed force
 column. These values represent the times at which forces begin affecting the model, and when the force accelerates. These times for "SwallowPerturbed" are 1.5 and 2.5, respectively
 */
void updateDynamicInputs(double time, double & freqI2, double & freqHinge, double & freqI1I3, double & freqN3, double & seaweedforce, double a, double frequencyiterationtime, double frequencyiterationtime2)
{
    freqI2 = freqI2array[timeAdjuster(timearray, time)];
    freqHinge = freqHingearray[timeAdjuster(timearray, time)];
    freqI1I3 = freqI1I3array[timeAdjuster(timearray, time)];
    freqN3 = freqN3array[timeAdjuster(timearray, time)];
    if(numberColumns == 6){
        double timeForceBegins, timeForceAccelerates;
        timeForceBegins = seaweedforcearray[0];
        timeForceAccelerates = seaweedforcearray[1];
        if (time > timeForceBegins){
            if (time < timeForceAccelerates)
            {
                seaweedforce = MAXSEAWEEDFORCE * ((a - .005)/.003)*(time - 1.5);
            }
            else
            {
                seaweedforce =  MAXSEAWEEDFORCE*((a - .005)/.003);
            }
        }
    }
}

string returnFirstLine(string s){
    string str;
    for (int i = 0; i < s.size(); i++)
    {
        if(s.at(i) != 10){
            str += (s.at(i));
        }
        else{
            break;
        }
    }
    return str;
}

int counter(string str)
{
    int count = 0;
    for (int i = 0; i < str.size(); i++){
        if(str.at(i) == 9){
            count++;
        }
    }
    return count+1;
}

int columnChecker()
{
    ifstream inFile;
    inFile.open("SlugInput2.txt");
    stringstream strStream;
    strStream << inFile.rdbuf();
    string str = strStream.str();
    string str2 = returnFirstLine(str); //only the first line
    return counter(str2);
}
