/*
  Copyright (C) 2021  Jhelam N. Deshpande
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
*/

//============================================================================

#include <cstdlib>	//standard C library
#include <iostream>	//standard input output library
#include <fstream>	//file stream library
#include <numeric>	
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>   //random number generator gsl       
#include <gsl/gsl_randist.h>	//gsl random distribution
#include <algorithm>
#include <math.h>	//standard math library
using namespace std;

//________________________________________________________________________________________
//--------------------------------------------------------Size of the local adaptation GRN
const int INPUT_SIZE=1;	//size of input layer
const int REGULATORY_SIZE=3;	//number of regulatory genes
const int GRN_TIME=20;	//number of iterations the GRN expression states go through

//________________________________________________________________________________________
//----------------------------------------------------------Class defining the individuals
class TInd{       
public:
	TInd();
	//GRN encodes niche optimum for local adaptation 
	double input_weights_optimum[2][INPUT_SIZE][REGULATORY_SIZE];          //Weights input to regulatory layer niche optimum
	double regulatory_weights_optimum[2][REGULATORY_SIZE][REGULATORY_SIZE];  //Weights regulatory layer within itself niche optimum   
	double output_weights_optimum[2][REGULATORY_SIZE];      //Weights regulatory layer to output layer niche optimum
	double regulatory_threshold_optimum[2][REGULATORY_SIZE];    //thresholds of regulatory layer niche optimum
	double regulatory_slope_optimum[2][REGULATORY_SIZE];    //slope of regulatory layer niche optimum
	//single quantitative locus encodes dispersal probability
	double disp_prob[2]; //single locus encodes dispersal
};

//----------------------------------------------------------Constructor for the class TInd
TInd::TInd(){     
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<INPUT_SIZE;j++)
		{
			for(int k=0; k<REGULATORY_SIZE;k++)
			{
				input_weights_optimum[i][j][k]=0;
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			for(int k=0; k<REGULATORY_SIZE;k++)
			{
				regulatory_weights_optimum[i][j][k]=0;
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			output_weights_optimum[i][j]=0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			regulatory_threshold_optimum[i][j]=0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			regulatory_slope_optimum[i][j]=0;
		}
	}
	disp_prob[0]=0;
	disp_prob[1]=0;

}

//________________________________________________________________________________________
//--------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
	TPatch();
	vector <TInd> females;    //females in the simulation
	vector <TInd> newfemales;   //vector to store new dispersers or newborn females
	vector <TInd> males;        
	vector <TInd> newmales;
	double measured_dispersal;    //fraction of individuals dispersing
	double measured_optimum;    //measured level of adaptation to a patch
	double grn_mortality;	//tracks number dying from non-viable dispersal GRN
	double environment;	//stores value of external environment
};

//--------------------------------------------------------Constructor for the class TPatch
TPatch::TPatch(){    
	females.clear();
	males.clear();
	newfemales.clear();
	newmales.clear();
	measured_dispersal=0;
	measured_optimum=0;
	grn_mortality=0;
	environment=0;
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Model Parameters 

int No;         //Initial number of individuals
double LAMBDA;      //Growth rate
double ALPHA;     //Intra-specific competition coefficient
double DISPERSAL_PROB;  // Dispersal probability
int BURN_IN_TIME;       //Number of time steps until range expansions begin---duration of burn-in period
int REPLICATES;         //Number of replicates
double DISP_MORT; //dispersal mortality
double EXTINCTION_PROB;//probability of local patch extinction
double VARIANCE_MIN; //variance of mutation
double MUT_RATE; //mutation rate at the beginning
double MUT_RATE_MIN;	//mutation rate after 5000 time steps
double GENETIC_VARIATION; //standing genetic variation
double VARIANCE_MAX;	//mutation effect standard deviation at the beginning of the simulation
double MUT_SLOW;	//slow mutation rate set to 0
double VAR_SLOW;	//standard deviation of mutation effects slow mutations, set to 0
double WIDTH;	//niche width
double SLOPE;	//slope of external environmental gradient
const int world_size_x = 500;//size of landscape in x direction
const int world_size_y = 5;//size of landscape in y direction
const int burn_in_x=10;//size of landscape before range expansion
TPatch world[world_size_x][world_size_y];//Creating patches

const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}


const int RS = 100;                 // random seed


// random number generator
		//specify_rng(time(NULL));
		//specify_rng(RS);


//________________________________________________________________________________________
//--------------------------------------------------------Miscellaneous required functions


//----------------------------------------------------function returns mean of two numbers
double mean(double a, double b) 
{
	return ((a+b)/double(2));
}



//--------------------------------------------------function returns mutated niche optimum
double mutate(double d,int t) 
{
	double new_trait;	//stores value of an allele with mutation added
	double mut;	//stores the mutation rate
	double var;	//stores the standard deviation of mutation effect
	if(t<5000)	//for first 5000 time steps
		mut=-(t*(MUT_RATE-MUT_RATE_MIN)/double(5000))+MUT_RATE;	//mutation rate decreases linearly with time
	else mut=MUT_RATE_MIN;	//mutation rate is constant after
	//VARIANCE_MIN=VARIANCE_MAX for all simulations anyway
	if(t<5000)	//same as for mutation rate
		var=-(t*(VARIANCE_MAX-VARIANCE_MIN)/double(5000))+VARIANCE_MAX;
	else var=VARIANCE_MIN;
	if(ran()<mut)	//mutation in the allele is introduced with the calculated probability
		new_trait= d+ gauss(var);  	//mutation size is drawn for a Gaussian with specified standard deviation
	else new_trait=d;	//if there is no mutation allele is unchanged
	if(ran()<MUT_SLOW)	
		new_trait= new_trait+ gauss(VAR_SLOW);  
	return new_trait;	//return allele with mutation added

}

//----------------------------------function returns mutated dispersal probability allele
double mutate_disp(double d) 
{
  double new_trait;	//stores the mutated allelic value
  if(ran()<MUT_RATE_MIN)	//introduce mutation with a specified probability
    new_trait= d+ gauss(VARIANCE_MIN);  //mutation drawn from Gaussian with specified standard deviation
  else new_trait=d;	//if there is no mutation allele is unchanged
  return new_trait;	//return mutated allele

}


//________________________________________________________________________________________
//----------------------------------------------------------Reading and setting parameters
void set_parameters()  
{
	string para[17];															//stores model parameter values
	string line;
	ifstream myfile ("input.txt");								//read the parameter input file
	int count=0;
	if (myfile.is_open())
	{
		while ( getline (myfile,line))
		{
			if(count%2==1)
			{
				para[count/2]=line;											//store only numeric values from the input file
			}
			count++;

		}
		myfile.close();
	}
	else cout << "Unable to open file";
	No = (int) std::atof(para[0].c_str());							//sets initial population size per patch
	LAMBDA=std::atof(para[1].c_str());									//sets mean fecundity of the female
	ALPHA=std::atof(para[2].c_str());										//sets intra-specific competition coefficient of the Beverton-Holt model
	BURN_IN_TIME= (int) std::atof(para[3].c_str());			//sets the number of time steps before the beginning of range expansion
	REPLICATES=(int) std::atof(para[4].c_str());				//sets number of replicate simulations that are run
	DISPERSAL_PROB=std::atof(para[5].c_str());					//sets dispersal probability if it is not genetically encoded
	DISP_MORT=std::atof(para[6].c_str());								//sets dispersal costs
	EXTINCTION_PROB=std::atof(para[7].c_str());					//sets probability of random patch extinction per time step per patch
	VARIANCE_MIN=std::atof(para[8].c_str());						//sets mutation effects after 5000 time steps
	MUT_RATE=std::atof(para[9].c_str());								//sets maximum mutation probability before 5000 time steps
	MUT_RATE_MIN=std::atof(para[10].c_str());						//sets mutation probability after 5000 time steps
	GENETIC_VARIATION=std::atof(para[11].c_str());			//sets standing genetic variation
	VARIANCE_MAX=std::atof(para[12].c_str());						//sets maximum mutation effects after 5000 time steps
	MUT_SLOW=std::atof(para[13].c_str());								//slow mutation rate---set to 0
	VAR_SLOW=std::atof(para[14].c_str());								//effects of slow mutations---set to 0
	WIDTH=std::atof(para[15].c_str());									//niche width
	SLOPE=std::atof(para[16].c_str());									//slope of the environmental gradient
}



//________________________________________________________________________________________
//----------------------------------------------------------------------Initialising model
void init_world()
{
	for(int x=0;x<world_size_x;x++)	//go through all patches in the landscape
	{
		for(int y=0;y<world_size_y;y++) //clear all individuals in a patch
		{

			world[x][y].females.clear(); //clear all females in the patch
			world[x][y].males.clear();	//clear all males in the patch
			world[x][y].newfemales.clear(); 	//clear newfemales vector
			world[x][y].newmales.clear();	//clear newmales vector
			if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1)	//setting the external environment in the range core to a constant
				world[x][y].environment=0.1;
			if(x<world_size_x/2-burn_in_x/2)
				world[x][y].environment=0.1+SLOPE*(world_size_x/2-burn_in_x/2-x);	//setting a linear environmental gradient towards the left of the landscape
			if(x>world_size_x/2+burn_in_x/2-1)
				world[x][y].environment=0.1-SLOPE*(world_size_x/2+burn_in_x/2-1-x);	//setting a linear environmental gradient towards the right of the landscape
			world[x][y].measured_optimum=0;	//set the measure level of adaptation of each patch to 0
			world[x][y].measured_dispersal=0;	//set the measured dispersal rate in each patch to 0
		}
	}
	//population is initialised only in the range core
	for(int x=world_size_x/2-burn_in_x/2;x<world_size_x/2+burn_in_x/2;x++)	//go through all patches in the range core
	{
		for(int y=0;y<world_size_y ;y++)
		{
		world[x][y].females.clear(); //clear all individuals in a patch
		world[x][y].males.clear();
		for(int n=0; n<No; n++)	//initially create No individuals in each patch
		{
			TInd newind;  //create a new individual
			
			 //initialising  network parameters for optimum with standing genetic variation
			for (int i = 0; i < 2; i++) 
			{
				for(int j=0;j<INPUT_SIZE;j++)
				{
					for(int k=0; k<REGULATORY_SIZE;k++)
					{
						newind.input_weights_optimum[i][j][k]=gauss(GENETIC_VARIATION);
					}
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					for(int k=0; k<REGULATORY_SIZE;k++)
					{
						newind.regulatory_weights_optimum[i][j][k]=gauss(GENETIC_VARIATION);
					}
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.output_weights_optimum[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.regulatory_threshold_optimum[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.regulatory_slope_optimum[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			newind.disp_prob[0]=ran(); //initialising single locus dispersal with standing genetic variation
			newind.disp_prob[1]=ran(); //initialising single locus dispersal with standing genetic variation
			
			if(ran()<0.5)	//individual has equal chance of being male or female
			world[x][y].females.push_back(newind); //adding new females
			else world[x][y].males.push_back(newind);//adding new males
		}

	}
}
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Output genotypes 
void output_genotype(ofstream& op5, int x_min,int x_max, int r, int t)
{
	for(int x=x_min; x<x_max;x++)
	{
		for(int y=0; y<world_size_y;y++)
		{
			//output 50 percent of females from a given region
			for(int f=0;f<world[x][y].females.size();f++)
			{
				//output 50 percent of females from a given region
				if(ran()<0.5){
					op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" ";


					for(int i=0;i<REGULATORY_SIZE;i++) //output GRN for optimum
					{

						for(int j=0; j<INPUT_SIZE;j++)
						{

							op5<<mean(world[x][y].females.at(f).input_weights_optimum[0][j][i],
								world[x][y].females.at(f).input_weights_optimum[1][j][i])<<" ";
						}

						for(int j=0; j<REGULATORY_SIZE;j++)
						{

							op5<<mean(world[x][y].females.at(f).regulatory_weights_optimum[0][i][j],
								world[x][y].females.at(f).regulatory_weights_optimum[1][i][j])<<" ";
						}


						op5<<mean(world[x][y].females.at(f).regulatory_threshold_optimum[0][i],
							world[x][y].females.at(f).regulatory_threshold_optimum[1][i])<<" ";
						op5<<mean(world[x][y].females.at(f).output_weights_optimum[0][i],
							world[x][y].females.at(f).output_weights_optimum[1][i])<<" ";
						op5<<mean(world[x][y].females.at(f).regulatory_slope_optimum[0][i],
							world[x][y].females.at(f).regulatory_slope_optimum[1][i])<<" ";
					}
					op5<<mean(world[x][y].females.at(f).disp_prob[0],
						world[x][y].females.at(f).disp_prob[1])<<" ";
					
					op5<<endl;
				}
			}
			//output 50 percent of males from a given region
			for(int m=0;m<world[x][y].males.size();m++)
			{
				if(ran()<0.5){
					op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" ";
					
					for(int i=0;i<REGULATORY_SIZE;i++) //output GRN for  optimum
					{

						for(int j=0; j<INPUT_SIZE;j++)
						{

							op5<<mean(world[x][y].males.at(m).input_weights_optimum[0][j][i],
								world[x][y].males.at(m).input_weights_optimum[1][j][i])<<" ";
						}

						for(int j=0; j<REGULATORY_SIZE;j++)
						{

							op5<<mean(world[x][y].males.at(m).regulatory_weights_optimum[0][i][j],
								world[x][y].males.at(m).regulatory_weights_optimum[1][i][j])<<" ";
						}


						op5<<mean(world[x][y].males.at(m).regulatory_threshold_optimum[0][i],
							world[x][y].males.at(m).regulatory_threshold_optimum[1][i])<<" ";
						op5<<mean(world[x][y].males.at(m).output_weights_optimum[0][i],
							world[x][y].males.at(m).output_weights_optimum[1][i])<<" ";
						op5<<mean(world[x][y].males.at(m).regulatory_slope_optimum[0][i],
							world[x][y].males.at(m).regulatory_slope_optimum[1][i])<<" ";
					}
					op5<<mean(world[x][y].males.at(m).disp_prob[0],
						world[x][y].males.at(m).disp_prob[1])<<" ";
					

					op5<<endl;
				}
			}
		}  

	}
}

//________________________________________________________________________________________
//--------------------Output population density,range front position.etc for a given patch

void output_metapopulation(ofstream& op, int x,int  y, int r, int t, int margin_x_left, int margin_x_right)
{

	op <<r<<" "<<t<<" "<<x<<" "<<y<<" "<<world[x][y].females.size()+world[x][y].males.size()<<" "<<
	double(world[x][y].females.size())/double(world[x][y].females.size()+world[x][y].males.size())<<" "<<
	world[x][y].measured_dispersal<<" "<<world[x][y].grn_mortality<<" "<<world[x][y].measured_optimum<<" ";
	op<<margin_x_left<<" "<<margin_x_right;
	op<<endl;


}

//---------------------------------------------return a subset of individuals given patches
vector <TInd> subset_metapopulation(int x_min, int x_max)
{
  vector <TInd> all_individuals;//stores all individuals in patches between x_min to x_max 
  for(int x=x_min;x<x_max;x++)
  {
  	for(int y=0; y<world_size_y; y++)
  	{
  		for(int f=0;f<world[x][y].females.size();f++)
  		{
  			all_individuals.push_back(world[x][y].females.at(f));
  		}
  		for(int m=0;m<world[x][y].males.size();m++)
  		{
  			all_individuals.push_back(world[x][y].males.at(m));
  		}
  	}
  }
  int no_of_individuals=100;
  if(no_of_individuals>all_individuals.size())
  	no_of_individuals=all_individuals.size();
  vector<TInd> random_individuals;
  for(int ind=0;ind<no_of_individuals;ind++)//draw no_of_individuals individuals
  { 
    int position_random=floor(ran()*all_individuals.size());//individuals drawn at random
    random_individuals.push_back(all_individuals.at(position_random));
}
return random_individuals;	//return 100 randomly chosen individuals
}



//________________________________________________________________________________________
//--------------------------------------Deciding coordinates of new patch while dispersing
vector<int> decide_patch(int x, int y,int t) 
{ 
	//nearest neighbor 8 dipsersal
	vector<int> coordinates;
	int newx=0;
	int newy=0;
	int decider=floor(ran()*double(8));
	switch(decider){
		case 0: newx++; 
		break;
		case 1: newy++;
		break;
		case 2: newx++;
		newy++;
		break;
		case 3:newx--;
		break;
		case 4:newy--;
		break;
		case 5:newx--;
		newy--;
		break;
		case 6: newx++;
		newy--;
		break;
		case 7:newx--;
		newy++;
		break;
		default: cout<<"error in nn8"<<endl;
	}
	newx=newx+x;
	newy=newy+y;
	//before range expansion begins reflecting boundary condition in x direction
	if(t<BURN_IN_TIME)
	{
		if (newx<world_size_x/2-burn_in_x/2)	//if the patch chosen is beyond the range core in the left direction
			newx=world_size_x/2-burn_in_x/2+1;	//send the individual one patch to the right
		if(newx>world_size_x/2+burn_in_x/2-1)
			newx=world_size_x/2+burn_in_x/2-2;
	}
	if(t>=BURN_IN_TIME)	//when range expansions begin
	{
		if (newx<0)	//simulations anyway stop when the margins are reached
			newx=0;
		if(newx==world_size_x)
			newx=world_size_x-1;
	}
	//torus in y direction
	if(newy<0)	//if the individual is at the bottom of the landscape
		newy=world_size_y - 1;	//send it to the top
	if(newy==world_size_y)
		newy=0;
	coordinates.push_back(newx);	
	coordinates.push_back(newy);
	return coordinates;	//return new coordinates of the individual
}
//________________________________________________________________________________________
//------------------------------------------function returns the trait calculated from GRN


double  grn_output(double inp[INPUT_SIZE],
	double iw[2][INPUT_SIZE][REGULATORY_SIZE], 
	double rw[2][REGULATORY_SIZE][REGULATORY_SIZE],
	double ow[2][REGULATORY_SIZE], 
	double rt[2][REGULATORY_SIZE], double rs[2][REGULATORY_SIZE])
{
	double input_weights[INPUT_SIZE][REGULATORY_SIZE];						//stores mid allelic value of input matrix weights
	double regulatory_weights[REGULATORY_SIZE][REGULATORY_SIZE];	//stores mid allelic value of regulatory matrix weights
	double output_weights[REGULATORY_SIZE];												//stores mid allelic value of output matrix weights
	double regulatory_threshold[REGULATORY_SIZE];									//stores mid allelic value of regulatory thresholds
	double regulatory_slope[REGULATORY_SIZE];											//stores mid allelic value of slopes of regulatory genes

		double genes[REGULATORY_SIZE][GRN_TIME]; //stores values of genes at each iteration
		double output[GRN_TIME];//stores values of output
		output[0]=0;
		double epsilon=0.0001;//tolerance
		double var=0; //stroes the variation in the equilibrium phenotype
		for (int i = 0; i < REGULATORY_SIZE; i++) // averaging between the two alleles 
		{
			for(int j=0; j<REGULATORY_SIZE; j++)
			{
				regulatory_weights[i][j]=mean(rw[0][i][j],rw[1][i][j]);
			}
			for(int j=0; j<INPUT_SIZE; j++)
			{
				input_weights[j][i]=mean(iw[0][j][i],iw[1][j][i]);
			}
			regulatory_threshold[i]=mean(rt[0][i],rt[1][i]);
			regulatory_slope[i]=mean(rs[0][i],rs[1][i]);
			output_weights[i]=mean(ow[0][i],ow[1][i]);
			genes[i][0]=2*ran()-1; //initialising the genes randomly with a number between -1 and 1
		}
		for(int time=1;time<GRN_TIME;time++)	//iterate through the gene-regulatory network
		{ 
			double output_store=0.0; //stores grn output at each time step
			for(int i=0;i<REGULATORY_SIZE;i++)	//go through all the genes in the regulatory layer
			{
				double regulatory_output=0.0;	//stores the input to a gene
				for(int j=0;j<INPUT_SIZE;j++)
				{
					regulatory_output=regulatory_output+input_weights[j][i]*inp[j]; //sums up all the inputs from input layer
				}
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					regulatory_output=regulatory_output+regulatory_weights[j][i]*genes[j][time-1]; //sums up the input from other regulatory genes
				}
				genes[i][time]=(2.0/(1.0+exp(-(regulatory_slope[i])*(regulatory_output-regulatory_threshold[i]))))-1; //gene expression state for a given iteration
				output_store=output_store+(genes[i][time]*output_weights[i]);  //store the output
			}
			output[time]=output_store;
			if(time >10)
			{
				var=var+(output[time]-output[time-1])*(output[time]-output[time-1]);
			}
		}
		var=sqrt(var)/(GRN_TIME-10);
		if(var<epsilon)
			return output[GRN_TIME-1] ;
		else
		{ 
			return 42.0;   
		}  
	}

//________________________________________________________________________________________
//-------------------------------------------------------------Density regulation function

	double densReg(double a) 
	{
		return(1 /(1+(a)));
	}



//calculate perturbed phenotype for a given individual
	double perturbed_trait(
		double iw[2][INPUT_SIZE][REGULATORY_SIZE], 
		double rw[2][REGULATORY_SIZE][REGULATORY_SIZE],
		double ow[2][REGULATORY_SIZE], 
		double rt[2][REGULATORY_SIZE], double rs[2][REGULATORY_SIZE])
	{
	//randomly draw the position of the perturbation
		double perturbed=0;
		double inp[1];
		inp[0]=0.5;
		int pos_i=floor(ran()*2);  //choose an allele to perturb
		int pos_j=floor(ran()*(REGULATORY_SIZE+INPUT_SIZE+3));	//choose the position and what part of the GRN to perturb
		int pos_k=floor(ran()*REGULATORY_SIZE);	//choose the position in the matrix
      if(pos_j<REGULATORY_SIZE) //perturb regulatory weights if the number drawn is 0,..,REGULATORY_SIZE-1
      {
      	double rw_new[2][REGULATORY_SIZE][REGULATORY_SIZE];//copy regulatory weights into a new matrix
      	for(int i=0;i<2; i++)
      	{
      		for(int j=0;j<REGULATORY_SIZE;j++)
      		{
      			for(int k=0;k<REGULATORY_SIZE;k++)
      			{
      				rw_new[i][j][k]=rw[i][j][k];
      			}
      		}
      	}
      	rw_new[pos_i][pos_j][pos_k]=rw[pos_i][pos_j][pos_k]+gauss(VARIANCE_MIN); //perturb regulatory weights at the chosen position
      	perturbed= grn_output(inp,iw,rw_new,ow,rt,rs); //calculate perturbed phenotype
      }
      else if(pos_j<REGULATORY_SIZE+INPUT_SIZE)	//if the random number drawn is between REGULATORY_SIZE and REGULATORY_SIZE+INPUT_SIZE-1 then perturb the input weights
      {
      	pos_j=pos_j-REGULATORY_SIZE;	//get the position within the input matrix to perturb
      	double iw_new[2][INPUT_SIZE][REGULATORY_SIZE];//copy the input weights into a new matrix
      	for(int i=0;i<2; i++)
      	{
      		for(int j=0;j<INPUT_SIZE;j++)
      		{
      			for(int k=0;k<REGULATORY_SIZE;k++)
      			{
      				iw_new[i][j][k]=iw[i][j][k];
      			}
      		}
      	}
      	iw_new[pos_i][pos_j][pos_k]=iw[pos_i][pos_j][pos_k]+gauss(VARIANCE_MIN);	//add a perturbation to the chosen input weight
      	perturbed= grn_output(inp,iw_new ,rw,ow,rt,rs);
      }
      else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+1)
      {
        double ow_new[2][REGULATORY_SIZE];//stores perturbed output weights  
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		ow_new[i][j]=ow[i][j];
        	}
        }
        ow_new[pos_i][pos_k]=ow[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow_new,rt,rs);
    }
    else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+2)
    {
        double rs_new[2][REGULATORY_SIZE];//stores perturbed slope of the gene        
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		rs_new[i][j]=rs[i][j];
        	}
        }
        rs_new[pos_i][pos_k]=rs[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow,rt,rs_new);
    }
    else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+3)
    {
        double rt_new[2][REGULATORY_SIZE];//stores perturbed  threshold of the gene
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		rt_new[i][j]=rt[i][j];
        	}
        }
        rt_new[pos_i][pos_k]=rt[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow,rt_new,rs);
    }
    return perturbed;
}


vector <double> sensitivity_to_mutation(int x_min,int x_max)
{
	vector<TInd> individuals;	//stores subset of individuals
	vector <double> sensitivity;	//stores sensitivity to mutation of local adaptation trait
	sensitivity.push_back(0);
	sensitivity.push_back(0);
	individuals=subset_metapopulation(x_min,x_max);	//subset individuals in a given patch cross section
	double inp[1];	//input to GRN
	inp[0]=0.5;
	int no_perturbations=100;	//number of perturbations to each GRN
	int count_opt=0;	//count number of viable GRNs
	if(individuals.size()>0){	//if there are individuals
		for(int p=0;p<individuals.size();p++)	//go through the individuals
		{
			//calculate unperturbed phenotype
			double unperturbed=grn_output(inp,individuals.at(p).input_weights_optimum ,individuals.at(p).regulatory_weights_optimum,
				individuals.at(p).output_weights_optimum,individuals.at(p).regulatory_threshold_optimum,individuals.at(p).regulatory_slope_optimum);
			//perturbed the GRN for the specified number of times
			for(int q=0;q<no_perturbations;q++){
				double perturbed=perturbed_trait(individuals.at(p).input_weights_optimum ,individuals.at(p).regulatory_weights_optimum,
					individuals.at(p).output_weights_optimum,individuals.at(p).regulatory_threshold_optimum,individuals.at(p).regulatory_slope_optimum);
				//for viable GRNs
				if(perturbed !=42 && unperturbed !=42){
					//calculate sum of squares of difference between perturbed and unperturbed trait
					sensitivity.at(0)=sensitivity.at(0)+(perturbed-unperturbed)*(perturbed-unperturbed);
					count_opt++;
				}

			}
		}
		individuals=subset_metapopulation(x_min,x_max);

		sensitivity.at(0)=sqrt(sensitivity.at(0)/double(count_opt));	//sensitivity of local adaptation trait
	}

	return sensitivity;	//returns sensitivity of local adaptation trait and 0 for dispersal trait sensitivity (since it is not defined for single locus)
}







//________________________________________________________________________________________
//-------------------------------------------------------------------Life cycle procedures
//---------------------------------------------------------------------Dispersal procedure
void disperse(int t)
{
	for(int x=0; x<world_size_x;x++)    //clears newmales and newfemales in all patches
	{
		for(int y=0; y<world_size_y;y++)
		{
			world[x][y].newfemales.clear();
			world[x][y].newmales.clear();
		}
	}
	for(int x=0;x<world_size_x;x++)   //dispersal loop, go through all patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			int count_dispersers=0;	//count number of individuals leaving a patch
			//inputs of the GRN stored in inp
			double population_density=double(world[x][y].males.size()+world[x][y].females.size());
			double normalised_density=population_density*ALPHA/(LAMBDA-1);

			for(int f=0;f<world[x][y].females.size();f++)	//go through all the females
			{

				//dispersal

				if(ran()< mean(world[x][y].females.at(f).disp_prob[0],world[x][y].females.at(f).disp_prob[1]))	//dispersal probability is the mid allelic value
				{
					std::vector<int> coor=decide_patch(x,y,t);	//decide where to disperse
						//dispersal mortality
					if(ran()>DISP_MORT)	//if the individual does not die while dispersing
						world[coor.at(0)][coor.at(1)].newfemales.push_back(world[x][y].females.at(f));	//add it to the nefemales vector in the target patch
					world[x][y].females.erase(world[x][y].females.begin()+f);	//remove the female from its natal patch
					f--;	//next female now occupies its position
					count_dispersers++;	//count dispersing individuals
				}
				

				
			}
			for(int m=0;m<world[x][y].males.size();m++)	//same as above for males
			{


				if(ran()< mean(world[x][y].males.at(m).disp_prob[0],world[x][y].males.at(m).disp_prob[1]))
				{
					std::vector<int> coor=decide_patch(x,y,t);
					if(ran()>DISP_MORT)
						world[coor.at(0)][coor.at(1)].newmales.push_back(world[x][y].males.at(m));
					world[x][y].males.erase(world[x][y].males.begin()+m);
					m--;
					count_dispersers++;
				}
				

				
			}
			
			world[x][y].measured_dispersal=double(count_dispersers)/population_density;	//store the measured dispersal rate as a patch property



		}
		
	}
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			if(world[x][y].newfemales.size()>0)	//if there are dispersers arriving in a patch
			{
				for(int f=0;f<world[x][y].newfemales.size();f++)
				{
					world[x][y].females.push_back(world[x][y].newfemales.at(f));	//add them to the vector of females in that patch
				}
			}
			if(world[x][y].newmales.size()>0)	//same for males as females
			{
				for(int m=0;m<world[x][y].newmales.size();m++)
				{
					world[x][y].males.push_back(world[x][y].newmales.at(m));
				}
			}
			world[x][y].newfemales.clear();	//clear newmales and newfemales vectors
			world[x][y].newmales.clear();

		}
	}

}


//---------------------------------------------------------------------reproduction procedure
void reproduce(int t) //reproduction loop
{ 
	for(int x=0; x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0; y<world_size_y;y++)
		{
			world[x][y].newfemales.clear();	//clear newfemales vector just in case
			world[x][y].newmales.clear();	//clear nemales vector
			world[x][y].measured_optimum=0;	//stores the measured adaptation to a patch
			int count_females=0;	//count the number of females that reproduce
			if(world[x][y].males.size() !=0 && world[x][y].females.size() != 0)	//if the patch is not empty
			{
				int Nf=world[x][y].females.size();	//count number of females
				int Nm=world[x][y].males.size();	//count number of males
				double alpha_net=0;	//stores net competition
				for(int f=0; f<Nf;f++)      //calculating net alpha for males and females
				{
					alpha_net=alpha_net+ALPHA;
				}
				for(int m=0; m<Nm;m++)
				{
					alpha_net=alpha_net+ALPHA;
				}
				for(int f=0;f<Nf;f++)	//go through all the females in the patch
				{
					int mate_position=floor(ran()*world[x][y].males.size());	//female mates with randomly chosen male
						double mean_lambda= LAMBDA;//set the mean lambda
						double inp[1];	//input to GRN
						inp[0]=0.5;
						double female_optimum=grn_output(inp,
							world[x][y].females.at(f).input_weights_optimum,
							world[x][y].females.at(f).regulatory_weights_optimum,world[x][y].females.at(f).output_weights_optimum,
							world[x][y].females.at(f).regulatory_threshold_optimum,world[x][y].females.at(f).regulatory_slope_optimum); //calculate niche optimum for the female
						if(female_optimum !=42.0){	//if the GRN is viable
							count_females++;	//count number of reproducing females
						double mean_offspring = 2*mean_lambda * densReg(alpha_net)*exp(-(world[x][y].environment-female_optimum)*(world[x][y].environment-female_optimum)/WIDTH);//beverton holt model with density independent mortality of offspring
						world[x][y].measured_optimum=world[x][y].measured_optimum+exp(-(world[x][y].environment-female_optimum)*(world[x][y].environment-female_optimum)/WIDTH);
						int no_of_babies= poisson(mean_offspring); //number of offspring Poisson distributed
						for(int b=0;b<no_of_babies;b++)         //setting parental characteristics
						{
							TInd newind;
							for(int i=0;i<REGULATORY_SIZE;i++) //niche optimum GRN inheritance
							{
								for(int j=0;j<INPUT_SIZE;j++)
								{
									if(ran()<0.5)
										newind.input_weights_optimum[0][j][i]=mutate(world[x][y].females.at(f).input_weights_optimum[0][j][i],t);
									else newind.input_weights_optimum[0][j][i]=mutate(world[x][y].females.at(f).input_weights_optimum[1][j][i],t);
									if(ran()<0.5)
										newind.input_weights_optimum[1][j][i]=mutate(world[x][y].males.at(mate_position).input_weights_optimum[0][j][i],t);
									else newind.input_weights_optimum[1][j][i]=mutate(world[x][y].males.at(mate_position).input_weights_optimum[1][j][i],t); 
								}
								for(int j=0;j<REGULATORY_SIZE;j++)
								{
									if(ran()<0.5)
										newind.regulatory_weights_optimum[0][i][j]=mutate(world[x][y].females.at(f).regulatory_weights_optimum[0][i][j],t);
									else newind.regulatory_weights_optimum[0][i][j]=mutate(world[x][y].females.at(f).regulatory_weights_optimum[1][i][j],t);
									if(ran()<0.5)
										newind.regulatory_weights_optimum[1][i][j]=mutate(world[x][y].males.at(mate_position).regulatory_weights_optimum[0][i][j],t);
									else newind.regulatory_weights_optimum[1][i][j]=mutate(world[x][y].males.at(mate_position).regulatory_weights_optimum[1][i][j],t);
								}
								if(ran()<0.5)
									newind.output_weights_optimum[0][i]=mutate(world[x][y].females.at(f).output_weights_optimum[0][i],t);
								else newind.output_weights_optimum[0][i]=mutate(world[x][y].females.at(f).output_weights_optimum[1][i],t);
								if(ran()<0.5)
									newind.output_weights_optimum[1][i]=mutate(world[x][y].males.at(mate_position).output_weights_optimum[0][i],t);
								else newind.output_weights_optimum[1][i]=mutate(world[x][y].males.at(mate_position).output_weights_optimum[1][i],t); 
								if(ran()<0.5)
									newind.regulatory_threshold_optimum[0][i]=mutate(world[x][y].females.at(f).regulatory_threshold_optimum[0][i],t);
								else newind.regulatory_threshold_optimum[0][i]=mutate(world[x][y].females.at(f).regulatory_threshold_optimum[1][i],t);
								if(ran()<0.5)
									newind.regulatory_threshold_optimum[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_threshold_optimum[0][i],t);
								else newind.regulatory_threshold_optimum[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_threshold_optimum[1][i],t); 
								if(ran()<0.5)
									newind.regulatory_slope_optimum[0][i]=mutate(world[x][y].females.at(f).regulatory_slope_optimum[0][i],t);
								else newind.regulatory_slope_optimum[0][i]=mutate(world[x][y].females.at(f).regulatory_slope_optimum[1][i],t);
								if(ran()<0.5)
									newind.regulatory_slope_optimum[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_slope_optimum[0][i],t);
								else newind.regulatory_slope_optimum[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_slope_optimum[1][i],t); 

							}
							if(ran()<0.5)	//dispersal probability inheritance
								newind.disp_prob[0]=mutate_disp(world[x][y].females.at(f).disp_prob[0]);
							else newind.disp_prob[0]=mutate_disp(world[x][y].females.at(f).disp_prob[1]);

							if(ran()<0.5)
								newind.disp_prob[1]=mutate_disp(world[x][y].males.at(mate_position).disp_prob[0]);
							else newind.disp_prob[1]=mutate_disp(world[x][y].males.at(mate_position).disp_prob[1]); 


							if(ran()<0.5)
								world[x][y].newfemales.push_back(newind);
							else world[x][y].newmales.push_back(newind);

						}
					}
				}
				world[x][y].measured_optimum=world[x][y].measured_optimum/count_females;

			}
			else 
			{
				world[x][y].females.clear();
				world[x][y].males.clear();
			}

		}
	}

}




void death()	//death procedure
{
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			world[x][y].females.clear();	//clear the previous generation males and females
			world[x][y].males.clear();
			world[x][y].females=world[x][y].newfemales;	//put the offspring in the males and females vector
			world[x][y].males=world[x][y].newmales;
			world[x][y].newfemales.clear();	//clear the newfemales vector
			world[x][y].newmales.clear();
		}
	}
}

void patch_extinction()	//random patch extinction
{
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0; y<world_size_y;y++)
		{
			if(ran()<EXTINCTION_PROB)	//the patch is cleared with a probability EXTINCTION_PROB
			{
				world[x][y].females.clear();
				world[x][y].males.clear();
			}
		}
	}
}

int main()	//main function
{
	//output the population size, sex ratio, measured dispersal at each patch
	ofstream op;
	op.open("output.txt");
	ofstream op1;
	op1.open("genotype_properties.txt");	//output the sensitivity to mutation of the local adaptation trait
	ofstream op5;	
	op5.open("output1.txt");	//output the genotypes
	specify_rng(RS);	//set seed for the random number generator
	set_parameters();	//set the model parameters
	op5 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y ";	//headers of output files

	for(int i=0;i<REGULATORY_SIZE;i++)
	{
		for(int j=0;j<INPUT_SIZE;j++)
		{
			op5<<"iw_o_"<<j+1<<i+1<<" ";
		}
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			op5<<"rw_o_"<<i+1<<j+1<<" ";
		}
		op5<<"rt_o_"<<i+1<<" ";
		op5<<"ow_o_"<<i+1<<" ";
		op5<<"rs_o_"<<i+1<<" ";
	}
	op5<<"disp_prob";

	op5<<endl;
	op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" "<<"N"<<" "<<"sex_ratio"<<" "<<"disp_rate grn_mortality measured ";
	op<<"margin_x_left margin_x_right";
	op<<endl;
	op1<<"rep t sensitivity_core sensitivity_core_disp sensitivity_front_left sensitivity_front_left_disp sensitivity_front_right sensitivity_front_right_disp"<<endl;
	for(int r=0; r<REPLICATES; r++)     //replicates
	{
		init_world();	//initialise the landscape
		int t=0;	//set timer to 0
		int margin_x_left=world_size_x/2-burn_in_x/2;	//store the left boundary of range from
		int margin_x_right=world_size_x/2+burn_in_x/2-1;	//stores the right boundary of the range front
		do{	//loop through time
			//output
			for(int x=0; x<world_size_x;x++)	//output patch properties
			{
				for(int y=0; y<world_size_y ;y++)
				{
					if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1 && t>BURN_IN_TIME-50 && t%50==0 )
						output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);	//output every 50 time steps since range expansions begin
					if(x<margin_x_left && world[x][y].females.size()+world[x][y].males.size()>0)
						margin_x_left=x;	//update the range front position based on farthest occupied patch
					if(x>margin_x_right && world[x][y].females.size()+world[x][y].males.size()>0)
						margin_x_right=x;
				} 
			}

			if(t>=BURN_IN_TIME && t%50==0)	//output range front properties
			{
				for(int x=margin_x_left; x<margin_x_left+5;x++)
				{
					for(int y=0; y<world_size_y ;y++)
					{

						output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);

					}
				} 

				for(int x=margin_x_right-4; x<margin_x_right+1;x++)
				{
					for(int y=0; y<world_size_y ;y++)
					{

						output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);

					}
				}
			}

			if(t<=BURN_IN_TIME)	//output sensitivity to mutation
			{
				if(t%2500==0)
				{
					vector <double> sensitivity;
					sensitivity=sensitivity_to_mutation(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2);

					op1<<r<<" "<<t<<" "<<sensitivity.at(0)<<" "<<sensitivity.at(1)<<" ";
					op1<<"NA NA NA NA"<<endl;
				}
			}
			if(t>=BURN_IN_TIME && t%50==0)
			{
				vector <double> sensitivity;
				sensitivity=sensitivity_to_mutation(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2);
				op1<<r<<" "<<t<<" "<<sensitivity.at(0)<<" "<<sensitivity.at(1)<<" ";
				sensitivity=sensitivity_to_mutation(margin_x_left,margin_x_left+5);
				op1<<sensitivity.at(0)<<" "<<sensitivity.at(1)<<" ";
				sensitivity=sensitivity_to_mutation(margin_x_right-4,margin_x_right+1);
				op1<<sensitivity.at(0)<<" "<<sensitivity.at(1)<<" "<<endl;
				
			} 

			if(t==BURN_IN_TIME-1)	//output genotypes at the end of the burn in period and at the end of the simulation
				output_genotype(op5,world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2,r,t);
			
//output the genotypes at the range front at the end of range expansion
			if(margin_x_right==world_size_x-1 && margin_x_left==0){
				output_genotype(op5,margin_x_right-4,margin_x_right+1,r,t);
				output_genotype(op5,margin_x_left,margin_x_left+5,r,t);
			}




			//life cycle
			disperse(t);	//dispersal is natal
			reproduce(t);	//individuals reproduce in target patches
			death();	//individuals die and are replaced by their offspring in the landscape
			patch_extinction();	//entire patches may go extinct with a specified probability
			t++;	//increase timer
			if(t>BURN_IN_TIME+40000)	//leave the loop if simulations take too long
				break;

		}
		while(margin_x_right!=world_size_x-1 || margin_x_left!=0);	//run simulations till the entire landscape is occupied
	}
	op.close();	//close output files
	op1.close();
	op5.close();
	return 0;
}







