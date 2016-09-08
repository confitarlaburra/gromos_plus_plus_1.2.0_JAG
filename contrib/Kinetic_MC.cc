#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include "../src/args/Arguments.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/fit/AtomDistances.h"
#include "../src/gmath/Stat.h"
#include <gsl/gsl_rng.h>


using namespace std;
using namespace gmath;
using namespace args;


// main

enum agents_enum {FULL=0,LEFT=1,RIGHT=2,SOLUTE=3};
enum IN_OUT{OUT=0,IN=1};
// Defines each reaction effect
void reactions(int i,int j,  vector<int> & agents);
double permeation(int i,int j,  int & lr_c, int & rl_c, list<double> &lt,double t);


int main(int argc, char **argv)  {
  
  Argument_List knowns;
  knowns <<"cycles"<<"kOn"<<"kOff"<<"kHop"<<"bins"
	 <<"PoreLoad"<<"Volume"<<"Nmol";  
         
  string usage = "# " + string(argv[0]);
  usage += "\n\t@cycles number of cycles \n";
  try {
    Arguments args(argc, argv, knowns, usage);

    // get the @cycles argument
    int cycles=0;
    if ( args.count("cycles") > 0)
      cycles=args.getValue<int>("cycles");
      
    cout<<"Number of Cycles: "<<cycles<<endl;
    
    int PoreLoad=args.getValue<int>("PoreLoad");
    cout<<"Pore Load       : "<<PoreLoad<<endl;
    
    int Nmol=args.getValue<int>("Nmol");
    cout<<"N mol           : "<<Nmol<<endl;
    

    int bins=args.getValue<int>("bins");
    cout<<"bins            : "<<bins<<endl;


    double Volume=args.getValue<double>("Volume");
    cout<<"Volume          : "<<Volume<<endl;
    
    double kOn=args.getValue<double>("kOn");
    cout<<"kOn             : "<<kOn<<endl;
    

    double kOff=args.getValue<double>("kOff");
    cout<<"kOff            : "<<kOff<<endl;
    

    double kHop=args.getValue<double>("kHop");
    cout<<"kHop            : "<<kHop<<endl;
    
        
    vector<Stat<double> >PermeationRates(2);
    Stat<double>  Survival_Times;
 

        
    // Tube state  vector
    vector<string> states(3); // FULL, LEFT, RIGHT
    states[0]="FULL";
    states[1]="LEFT";
    states[2]="RIGHT";
    
    // Agents Vector
    vector<int> agents(4,0); // tube and molecules to be loaded
    agents[FULL]=1; // 0 agent
    agents[LEFT]=0; // 1 agent
    agents[RIGHT]=0; // 2 agent 
    agents[SOLUTE]=Nmol; // [nm^3]^-1 or 0.33 mol/L
    
  
    
    gsl_matrix * rates=  gsl_matrix_calloc(states.size(), states.size());
    // set initial propensities matrix elements
    gsl_matrix * propensities=  gsl_matrix_calloc(states.size(), states.size());



    // set rates matrix elements
    gsl_matrix_set (rates, FULL, LEFT, kOn);
    gsl_matrix_set (rates, FULL, RIGHT, kOn);
    
    gsl_matrix_set (rates, LEFT, FULL, kOff);
    gsl_matrix_set (rates, LEFT, RIGHT, kHop);

    gsl_matrix_set (rates, RIGHT, FULL, kOff);
    gsl_matrix_set (rates, RIGHT, LEFT, kHop);
    
   
    // Permeation variabless
    int permeation_counter_RL=PoreLoad/2;
    int permeation_counter_LR=PoreLoad/2 + PoreLoad%2;
    int Total_permeation_LR=0;
    int Total_permeation_RL=0;
    list<double> hop_times(PoreLoad,0.0);

    //gillespie algorithm variables
    double Time=0.0;
    double a0 =0;
    double partial_a0=0;
    double tau = 0.0;
    
    int counter=0;
    // Random seed environment
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    srand((time(0))); // srand & time are built-in
    int seed = random();

    gsl_rng_set(r,seed); // and feed with a random seed

    double  r1 = gsl_rng_uniform (r);
    double  r2 = gsl_rng_uniform (r);
    

    for (int i = 0; i < cycles; i++) {
      a0=0;
      partial_a0=0;
      r1 = gsl_rng_uniform (r);
      r2 = gsl_rng_uniform (r);
      // compute propensities (alpha 0,0....alpha n,n)
      
      gsl_matrix_set (propensities, FULL, LEFT, gsl_matrix_get(rates,FULL,LEFT)*agents[FULL]*(agents[SOLUTE]/Volume));
      gsl_matrix_set (propensities, FULL, RIGHT, gsl_matrix_get(rates,FULL,RIGHT)*agents[FULL]*(agents[SOLUTE]/Volume));
      
      gsl_matrix_set (propensities, LEFT, FULL, gsl_matrix_get(rates,LEFT,FULL)*agents[LEFT]);
      gsl_matrix_set (propensities, LEFT, RIGHT, gsl_matrix_get(rates,LEFT,RIGHT)*agents[LEFT]);
      
      gsl_matrix_set (propensities, RIGHT, FULL, gsl_matrix_get(rates,RIGHT,FULL)*agents[RIGHT]);
      gsl_matrix_set (propensities, RIGHT, LEFT, gsl_matrix_get(rates,RIGHT,LEFT)*agents[RIGHT]);
      
      //Get alpha 0 
      for (int i = 0; i < states.size(); i++) 
	for (int j = 0; j < states.size(); j++)
	  a0+=gsl_matrix_get (propensities, i, j);
      
      tau=(1/a0)*log(1/r1);
      Time+=tau;
      
      // Select reaction to be performed 
      for (int i = 0; i < states.size(); i++) {
	for (int j = 0; j < states.size(); j++) {
	  counter=0;
	  if ((1/a0)*partial_a0 <= r2)
	    counter++;
	  partial_a0+=gsl_matrix_get (propensities, i, j);
	  if ((1/a0)*partial_a0 > r2)
	    counter++;
	  if (counter==2) {
	    reactions(i,j,agents);
	    double surv_time = permeation(i,j,permeation_counter_LR,permeation_counter_RL,hop_times,Time);
	    //permeation_event
	    if ( permeation_counter_LR == PoreLoad+1 ) {
	      PermeationRates[0].addval(Time - surv_time);
	      permeation_counter_LR=PoreLoad;
              permeation_counter_RL=0;
	      Total_permeation_LR++;
	    }
	    else if ( permeation_counter_RL == (PoreLoad+1) ) {
	      PermeationRates[1].addval(Time - surv_time);
	      permeation_counter_RL=PoreLoad;
              permeation_counter_LR=0;
	      Total_permeation_RL++;
	    }
            //survival event
            if(surv_time != 0.0)
              Survival_Times.addval(Time-surv_time);
	    break;
	  }
	}
	if (counter==2)
	  break;
      }
    }

    cout<<"Total Permeations/unit Time= "<<(Total_permeation_LR+Total_permeation_RL)/Time << endl;
    cout<< "Permeations Average LR: " << PermeationRates[0].ave()<<" +- "<<PermeationRates[0].ee()<< endl;
    cout<< "Permeations Average RL: " <<PermeationRates[1].ave()<<" +- "<<PermeationRates[1].ee()<<endl;
    cout<< "Survival Time Average: " << Survival_Times.ave() << " +- " << Survival_Times.ee() << endl;
    
    
    Survival_Times.dist_init(Survival_Times.min(),
				 Survival_Times.max(),
				 bins);  

    ofstream output;
    output.open("Survival.out");
    output<<"#Time    Survival time"<<endl;
    Survival_Times.distribution().write_normalized(output);
    output.close();

    PermeationRates[0].dist_init(PermeationRates[0].min(),
				 PermeationRates[0].max(),
				 bins);   

    output.open("LR.out");
    output<<"#Time    Perm Rate Prob"<<endl;
    PermeationRates[0].distribution().write_normalized(output);
    output.close();

    output.open("RL.out");
    output<<"#Time    Perm Rate Prob"<<endl;
    PermeationRates[1].dist_init(PermeationRates[1].min(),
				 PermeationRates[1].max(),
				 bins);   
    
    PermeationRates[1].distribution().write_normalized(output);
    output.close();
   
    
  } catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

//functions and objects implementations

void reactions(int i,int j , vector<int> & agents) {
  if (i==FULL) {
    if (j==RIGHT) {
      agents[FULL]--;
      agents[RIGHT]++;
      agents[SOLUTE]--;
    }
    if (j==LEFT) {
      agents[FULL]--;
      agents[LEFT]++;
      agents[SOLUTE]--;
    }
  }

  if (i==RIGHT) {
    if (j==LEFT) {
      agents[RIGHT]--;
      agents[LEFT]++;
    }
    if (j==FULL) {
      agents[FULL]++;
      agents[RIGHT]--;
      agents[SOLUTE]++;
    }
  }

  if (i==LEFT) {
    if (j==RIGHT) {
      agents[RIGHT]++;
      agents[LEFT]--;
    }
    if (j==FULL) {
      agents[FULL]++;
      agents[LEFT]--;
      agents[SOLUTE]++;
    }
  }


}



double permeation(int i,int j , int & permeation_counterLR, int & permeation_counterRL, list<double> & hop_t, double t) {
  double val = 0.0;
  if (i==RIGHT && j==LEFT) {
    permeation_counterRL++;
    if(permeation_counterLR > 0)
      permeation_counterLR--;
    hop_t.push_back(t);
    val = hop_t.front();
    hop_t.pop_front();
  }
  else if (i==LEFT && j==RIGHT) {
    permeation_counterLR++;
    if(permeation_counterRL > 0)
      permeation_counterRL--;
    hop_t.push_front(t);
    val = hop_t.back();
    hop_t.pop_back();
  }
  return val;
}


