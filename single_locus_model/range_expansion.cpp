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


#include <cstdlib>    //standard C library
#include <iostream>   //standard input output function library
#include <fstream>    //file stream library
#include <string>
#include <numeric>
#include <vector>
#include <gsl/gsl_rng.h>    //random number generator gsl library      
#include <gsl/gsl_randist.h>  //gsl random distribution library
#include <math.h>   //standard math library
#include <algorithm>
using namespace std;

const int NO_OF_LOCI=1;   //number of loci per trait

//________________________________________________________________________________________
//----------------------------------------------------------Class defining the individuals
class TInd{       
public:
  TInd();
  double lambda[2]; //intrinsic growth rate---Beverton-Holt model
  double alpha[2];  //intra-specific competition coefficient---Beverton-Holt model
  double disp_prob[2][NO_OF_LOCI];  //2 alleles and NO_OF_LOCI loci for dispersal trait
  double optimum[2][NO_OF_LOCI];    //2 alleles and NO_OF_LOCI loci for niche optimum trait
};

//----------------------------------------------------------Constructor for the class TInd
TInd::TInd(){     
  lambda[0]=0;      
  lambda[1]=0;        
  alpha[0]=0;
  alpha[1]=0;
  for(int i=0;i<NO_OF_LOCI;i++)
  {
    disp_prob[0][i]=0;
    disp_prob[1][i]=0;
    optimum[0][i]=0;
    optimum[1][i]=0;

  }

  //dispersal_probability=0;
}

//________________________________________________________________________________________
//--------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
  TPatch();
  vector <TInd> females;    //females in the simulation
  vector <TInd> newfemales;   //vector to store new dispersers or new born females
  vector <TInd> males;        
  vector <TInd> newmales;
  double measured_dispersal;    //fraction of individuals dispersing from a patch
  double environment;   //external environmental gradient value in a patch
  double measured;      //stores the measured level of adaptation in a patch
};

//--------------------------------------------------------Constructor for the class TPatch
TPatch::TPatch(){    
  females.clear();
  males.clear();
  newfemales.clear();
  newmales.clear();
  measured_dispersal=0;
  environment=0;
  measured=0;
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Model Parameters 

int No;         //Initial number of individuals
double LAMBDA;      //Intrinsic growth rate---Beverton-Holt model
double ALPHA;    //Intra-specific competition coefficient---Beverton-Holt model
double DISPERSAL_PROB;  // Dispersal probability
int BURN_IN_TIME;       //Number of time steps before beginning range expansion
int REPLICATES;         //Number of replicates
double DISP_MORT; //dispersal costs
double EXTINCTION_PROB;//probability of local patrch extinction
double VARIANCE; //standard deviation of mutation effects for both traits
double MUT_RATE; //mutation rate for both traits
double WIDTH; //niche width
double SLOPE; //slope of environmental gradient
const int world_size_x = 500;//size of landscape in x direction
const int world_size_y = 5;//size of landscape in y direction
const int burn_in_x=10;//size of landscape before range expansion---range core size
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


//--------------------------------------------------function returns mutated allelic value
double mutate(double d) 
{
  double new_trait; 
  if(ran()<MUT_RATE)  //mutation added to an allele with a probability MUT_RATE
    new_trait= d+ gauss(VARIANCE);  //mutation drawn from Gaussian
  else new_trait=d; //if there is no mutation no change to allele
  return new_trait;   //return mutated allele

}


//________________________________________________________________________________________
//--------------------------------------------------------Inputting and setting parameters
void set_parameters()  
{
  string para[12];  //stores model parameter values
  string line;      
  ifstream myfile ("input.txt");    //read the parameter input file
  int count=0;
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      if(count%2==1)
      {
        para[count/2]=line;   //store only numeric values from the file
      }
      count++;

    }
    myfile.close();
  }
  else cout << "Unable to open file";
  No = (int) std::atof(para[0].c_str());    //set initial number of individuals in landscape
  LAMBDA=std::atof(para[1].c_str());  //set intrinsic growth rate
  ALPHA=std::atof(para[2].c_str());   //set intra-specific competition coefficient of Beverton-Holt model
  BURN_IN_TIME= (int) std::atof(para[3].c_str());   //set duration of burn-in period
  REPLICATES=(int) std::atof(para[4].c_str());    //set number of replicates
  DISPERSAL_PROB=std::atof(para[5].c_str());    //set dispersal probability
  DISP_MORT=std::atof(para[6].c_str());       //set dispersal costs
  EXTINCTION_PROB=std::atof(para[7].c_str()); //set probability of random patch extinctions
  VARIANCE=std::atof(para[8].c_str());    //set sd of mutation effects
  MUT_RATE=std::atof(para[9].c_str());    //set mutation rate
  WIDTH=std::atof(para[10].c_str());  //set niche width
  SLOPE=std::atof(para[11].c_str());  //set slope of environmental gradient
}


//________________________________________________________________________________________
//----------------------------------------------------------------------Initialising model
void init_world()
{
  for(int x=0;x<world_size_x;x++)
  {
    for(int y=0;y<world_size_y ;y++)
    {
    //clear all individuals in a patch
      world[x][y].females.clear(); 
      world[x][y].males.clear();
      world[x][y].newfemales.clear(); 
      world[x][y].newmales.clear();
      if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1)   //set environmental gradient to constant in core
        world[x][y].environment=0.1;
      if(x<world_size_x/2-burn_in_x/2)
        world[x][y].environment=0.1+SLOPE*(world_size_x/2-burn_in_x/2-x); //linear increase in environment to the left of range core
      if(x>world_size_x/2+burn_in_x/2-1)
        world[x][y].environment=0.1-SLOPE*(world_size_x/2+burn_in_x/2-1-x); //linear increase in environment to the right of range core
      world[x][y].measured=0;
      world[x][y].measured_dispersal=0;


    }
  }
 //seed population only in the central 10 patches in x direction
  for(int x=world_size_x/2-burn_in_x/2;x<world_size_x/2+burn_in_x/2;x++)
  {
    for(int y=0;y<world_size_y ;y++)
    {
    //clear all individuals in a patch
      for(int n=0; n<No; n++)
      {
          TInd newind;      //stores new individuals
          newind.lambda[0]=LAMBDA;
          newind.lambda[1]=LAMBDA;
          newind.alpha[0]=ALPHA;
          newind.alpha[1]=ALPHA;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            newind.disp_prob[0][i]=ran();   // initialise dispersal probability with a uniform distribution
            newind.disp_prob[1][i]=ran();
            newind.optimum[0][i]=ran();
            newind.optimum[1][i]=ran(); // initialise measured optimum with a uniform distribution
          }
          if(ran()<0.5) //individual has equal chance of being male or female
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
    for(int x=x_min; x<x_max;x++) //go through all patches
    {
      for(int y=0; y<world_size_y ;y++)
      {
       for(int f=0;f<world[x][y].females.size();f++)  //output genotypes of half the individuals
       {
        if(ran()<0.5){
          op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" " ;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].females.at(f).disp_prob[0][i],world[x][y].females.at(f).disp_prob[1][i])<<" ";
          }
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].females.at(f).optimum[0][i],world[x][y].females.at(f).optimum[1][i])<<" ";
          }
          op5<<endl;
        }
      }
      for(int m=0;m<world[x][y].males.size();m++)
      {
        if(ran()<0.5){
          op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" " ;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].males.at(m).disp_prob[0][i],world[x][y].males.at(m).disp_prob[1][i])<<" ";
          }
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].males.at(m).optimum[0][i],world[x][y].males.at(m).optimum[1][i])<<" ";
          }

          op5<<endl;
        }
      }

    }
  }



}

//________________________________________________________________________________________
//----------------------------------------------------------Output population density.etc

void output_metapopulation(ofstream& op, int x,int  y, int r, int t, int margin_x_left, int margin_x_right)
{
  op <<r<<" "<<t<<" "<<x<<" "<<y<<" "<<world[x][y].females.size()+world[x][y].males.size()<<" "<<
  double(world[x][y].females.size())/double(world[x][y].females.size()+world[x][y].males.size())<<" "<<
  world[x][y].measured_dispersal<<" "<<world[x][y].measured<<" "<<margin_x_left<<" "<<margin_x_right<<" ";
  op<<endl;


}




//________________________________________________________________________________________
//--------------------------------------Deciding coordinates of new patch while dispersing
vector<int> decide_patch(int x, int y,int t) 
{ 
  //nearest neighbor 8 dispersal
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
    if (newx<world_size_x/2-burn_in_x/2)  //if the target patch is to the left of the range core, send the individual to the right
      newx=world_size_x/2-burn_in_x/2+1;
    if(newx>world_size_x/2+burn_in_x/2-1)//same as above for right side
      newx=world_size_x/2+burn_in_x/2-2;
  }
  if(t>=BURN_IN_TIME) //after range expansions begin, does not matter
  {
    if (newx<0)
      newx=0;
    if(newx==world_size_x)
      newx=world_size_x-1;
  }
  //torus in y direction
  if(newy<0)
    newy=world_size_y - 1;
  if(newy==world_size_y)
    newy=0;
  coordinates.push_back(newx);
  coordinates.push_back(newy);
  return coordinates;

}


//________________________________________________________________________________________
//----------------------------function returns the dispersal rate calculated from genotype

double  dispersal_probability_calc(double c[2][NO_OF_LOCI], double density)
{
  double d_mean=0;
  for(int i=0;i<NO_OF_LOCI;i++)
  {
    d_mean=d_mean+mean(c[0][i],c[1][i]);  //mid allelic value for each locus
  }
  return (d_mean/double(NO_OF_LOCI)); //average trait value over all loci
}


//________________________________________________________________________________________
//-----------------------------function returns the niche optimum calculated from genotype


double  optimum_calc(double c[2][NO_OF_LOCI])
{
  double d_mean=0;
  for(int i=0;i<NO_OF_LOCI;i++)
  {
    d_mean=d_mean+mean(c[0][i],c[1][i]);
  }
  return (d_mean/double(NO_OF_LOCI));
}



//________________________________________________________________________________________
//-------------------------------------------------------------Density regulation function

double densReg(double a) 
{
  return(1 /(1+(a)));
}



//________________________________________________________________________________________
//-------------------------------------------------------------------Life cycle procedures


//---------------------------------------------------------------------Dispersal procedure
void disperse(int t)
{
  for(int x=0; x<world_size_x;x++)    //clears newmales and newfemales vectors
  {
    for(int y=0; y<world_size_y ;y++)
    {
      world[x][y].newfemales.clear();
      world[x][y].newmales.clear();
    }
  }
  for(int x=0;x<world_size_x;x++)   //go through all patches
  {
    for(int y=0;y<world_size_y ;y++)
    {
      int count_dispersers=0; //stores the number of individuals dispersing from a patch
      double population_density=double(world[x][y].males.size()+world[x][y].females.size());
      for(int f=0;f<world[x][y].females.size();f++)
      {
        //dispersal probability calculation
        double dispersal_probability; //stores dispersal probability of each individual    
        //calculate dispersal probability
        dispersal_probability=dispersal_probability_calc(world[x][y].females.at(f).disp_prob,population_density*ALPHA/(LAMBDA-1));
        //dispersal
        if(ran()< dispersal_probability)  //individuals disperse with a probability
        {
          std::vector<int> coor=decide_patch(x,y,t);  //stores where individuals disperse
            //dispersal mortality
          if(ran()>DISP_MORT) //if individuals do not die during dispersal
            world[coor.at(0)][coor.at(1)].newfemales.push_back(world[x][y].females.at(f));  //send them to target patch
          world[x][y].females.erase(world[x][y].females.begin()+f); //remove this individual from present patch in any case
          f--;  //next female now is in the same position as this female
          count_dispersers++; //count number of dispersing females
        }
        
      }
      //same for males
      for(int m=0;m<world[x][y].males.size();m++)
      {
        double dispersal_probability; //stores dispersal probability of each individual
        dispersal_probability=dispersal_probability_calc(world[x][y].males.at(m).disp_prob,population_density*ALPHA/(LAMBDA-1));
        if(ran()< dispersal_probability)
        {
          std::vector<int> coor=decide_patch(x,y,t);
          if(ran()>DISP_MORT)
            world[coor.at(0)][coor.at(1)].newmales.push_back(world[x][y].males.at(m));
          world[x][y].males.erase(world[x][y].males.begin()+m);
          m--;
          count_dispersers++;
        }
      }

      //store the measured dispersal rate
      world[x][y].measured_dispersal=double(count_dispersers)/population_density;

    }
    
  }
  //add disperser to their target patch
  for(int x=0;x<world_size_x;x++)
  {
    for(int y=0;y<world_size_y ;y++)
    {
      if(world[x][y].newfemales.size()>0)
      {
        for(int f=0;f<world[x][y].newfemales.size();f++)
        {
          world[x][y].females.push_back(world[x][y].newfemales.at(f));
        }
      }
      if(world[x][y].newmales.size()>0)
      {
        for(int m=0;m<world[x][y].newmales.size();m++)
        {
          world[x][y].males.push_back(world[x][y].newmales.at(m));
        }
      }
      world[x][y].newfemales.clear();
      world[x][y].newmales.clear();

    }
  }

}

//---------------------------------------------------------------------reproduction procedure
void reproduce() //reproduction loop
{ 
  for(int x=0; x<world_size_x;x++)
  {
    for(int y=0; y<world_size_y ;y++)
    {
      double measured=0;  //stores measured level of adaptation in a patch
      world[x][y].newfemales.clear(); //clear newmales and newfemales vectors
      world[x][y].newmales.clear();
      if(world[x][y].males.size() !=0 && world[x][y].females.size() != 0) //if there are individuals in a patch
      {
        int Nf=world[x][y].females.size();  //count number of females
        int Nm=world[x][y].males.size();  //count number of males
        double alpha_net=0;
        for(int f=0; f<Nf;f++)      //calculating net alpha for males and females, here equal to ALPHA*N because all have same alpha value
        {
          alpha_net=alpha_net+mean(world[x][y].females.at(f).alpha[0],world[x][y].females.at(f).alpha[1]);
        }
        for(int m=0; m<Nm;m++)
        {
          alpha_net=alpha_net+mean(world[x][y].males.at(m).alpha[0],world[x][y].males.at(m).alpha[1]);
        }
        for(int f=0;f<Nf;f++) //go through all the females
        {
          int mate_position=floor(ran()*world[x][y].males.size());  //randomly choose a mate
            double mean_lambda= mean(world[x][y].females.at(f).lambda[0],world[x][y].females.at(f).lambda[1]);// setting mean as intrinsic growth rate
            double mean_optimum=optimum_calc(world[x][y].females.at(f).optimum);  //calculate niche optimum of the female
            //mean number of offspring after density regulation and density independent mortality due to maladaptation
            double mean_offspring = 2*mean_lambda * densReg(alpha_net)*exp(-(world[x][y].environment-mean_optimum)*(world[x][y].environment-mean_optimum)/WIDTH);//beverton holt model
            measured=measured+exp(-(world[x][y].environment-mean_optimum)*(world[x][y].environment-mean_optimum)/WIDTH);  //calculate mean level of adaptation to patch
            int no_of_babies= poisson(mean_offspring); //number of offspring Poisson distributed
            for(int b=0;b<no_of_babies;b++)         //setting parental characteristics
            {
              TInd newind;  //create the offspring
              if(ran()<0.5) //lambda and alpha traits cannot mutate, are just fixed parameters
                newind.lambda[0]=world[x][y].females.at(f).lambda[0];
              else newind.lambda[0]=world[x][y].females.at(f).lambda[1];
              if(ran()<0.5)
                newind.lambda[1]=world[x][y].males.at(mate_position).lambda[0];
              else newind.lambda[1]=world[x][y].males.at(mate_position).lambda[1];
              if(ran()<0.5)
                newind.alpha[0]=world[x][y].females.at(f).alpha[0];
              else newind.alpha[0]=world[x][y].females.at(f).alpha[1];
              if(ran()<0.5)
                newind.alpha[1]=world[x][y].males.at(mate_position).alpha[0];
              else newind.alpha[1]=world[x][y].males.at(mate_position).alpha[1];
              //inheritance of dispersal
              for(int i=0;i<NO_OF_LOCI;i++)
              {
                if(ran()<0.5)
                  newind.disp_prob[0][i]=mutate(world[x][y].females.at(f).disp_prob[0][i]);
                else newind.disp_prob[0][i]=mutate(world[x][y].females.at(f).disp_prob[1][i]);
                if(ran()<0.5)
                  newind.disp_prob[1][i]=mutate(world[x][y].males.at(mate_position).disp_prob[0][i]);
                else newind.disp_prob[1][i]=mutate(world[x][y].males.at(mate_position).disp_prob[1][i]);
              }
              //inheritance of local adaptation
              for(int i=0;i<NO_OF_LOCI;i++)
              {
                if(ran()<0.5)
                  newind.optimum[0][i]=mutate(world[x][y].females.at(f).optimum[0][i]);
                else newind.optimum[0][i]=mutate(world[x][y].females.at(f).optimum[1][i]);
                if(ran()<0.5)
                  newind.optimum[1][i]=mutate(world[x][y].males.at(mate_position).optimum[0][i]);
                else newind.optimum[1][i]=mutate(world[x][y].males.at(mate_position).optimum[1][i]);
              }
              if(ran()<0.5)
                world[x][y].newfemales.push_back(newind);
              else world[x][y].newmales.push_back(newind);
            }
          }
          world[x][y].measured=measured/double(world[x][y].females.size());   //store measured level of adaptation to the patch
        }
        else 
        {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }


      }

    }

  }

//---------------------------------------------------------------------death procedure

  void death()
  {
    for(int x=0;x<world_size_x;x++) //go through all patches
    {
      for(int y=0;y<world_size_y ;y++)
      {
        world[x][y].females.clear();  //old generation dies
        world[x][y].males.clear();
        world[x][y].females=world[x][y].newfemales; //old generation is replaced by offspring
        world[x][y].males=world[x][y].newmales;
        world[x][y].newfemales.clear(); //clear newfemales and newmales vector
        world[x][y].newmales.clear();
      }
    }
  }

  void patch_extinction() //go through all patches
  {
    for(int x=0;x<world_size_x;x++)
    {
      for(int y=0; y<world_size_y ;y++)
      {
        if(ran()<EXTINCTION_PROB) //entire patches go extinct with this probability
        {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }
      }
    }
  }


  int main()
  {
  //output the population size, sex ratio, measured dispersal at each patch
    ofstream op;
    op.open("output.txt");
    ofstream op5;
    op5.open("output1.txt");

    //column names of output files
    op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" "<<"N"<<" "<<"sex_ratio"<<" "<<"disp_rate measured margin_x_left margin_x_right ";
    /*for(int i=0;i<NO_OF_LOCI;i++)
    {
      op<<"d_"<<i<<" ";
    }*/
    op<<endl;

    op5 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" ";
    for(int i=0;i<NO_OF_LOCI;i++)
    {
      op5<<"d_"<<i<<" ";
    }
    for(int i=0;i<NO_OF_LOCI;i++)
    {
      op5<<"o_"<<i<<" ";
    }
    op5<<endl;
    //seed 

    specify_rng(RS);
    //read parameters from input file
    set_parameters();
  for(int r=0; r<REPLICATES; r++)     //run replicates
  {
    init_world(); //initialise the landscape
    int t=0;  //set time count to 0
    int margin_x_left=world_size_x/2-burn_in_x/2; //track range front position
    int margin_x_right=world_size_x/2+burn_in_x/2-1;;
    do{ //start the time loop
      //output patch properties
      for(int x=0; x<world_size_x;x++)
      {
        for(int y=0; y<world_size_y ;y++)
        {

         if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1 && t>BURN_IN_TIME-50 && t%50==0 )//output the central 10x5 patches
          output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);
        if(x<margin_x_left && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_left=x;
        if(x>margin_x_right && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_right=x;
      } 
    }

  if(t>=BURN_IN_TIME && t%50==0)
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
  //output genotypes at the range core at the beginning of range expansion
  if(t==BURN_IN_TIME-1)
    output_genotype(op5,world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2,r,t);

//output the genotypes at the range front at the end of range expansion
  if(margin_x_right==world_size_x-1 && margin_x_left==0){
    output_genotype(op5,margin_x_right-4,margin_x_right+1,r,t);
    output_genotype(op5,margin_x_left,margin_x_left+5,r,t);
  }


  //life cycle
  disperse(t);  //dispersal is natal
  reproduce();  //individuals reproduce in target patches
  death();  //old generation dies and is replaced by offspring
  patch_extinction(); //random patch extinction
  t++;  //increase time count
  if(t>BURN_IN_TIME+40000)  //if range expansions take too long
    break;  //get out of the time loop

}
while(margin_x_right!=world_size_x-1 || margin_x_left!=0); //run till the farthest patches are occupied

}
//close output files
op.close();
op5.close();
return 0;

}







