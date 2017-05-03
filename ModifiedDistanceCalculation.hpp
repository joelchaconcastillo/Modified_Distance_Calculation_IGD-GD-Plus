/*
   Autor: Joel Chacón Castillo
   Date: 02/may/2017

   This is based in the paper:
   "Modified Distance Calculation in Generational Distance and Inverted Generational Distance"
       by Hisao Ishibuchi, Hiroyuki Masuda, Yuki Tanigaki, and Yusuke Nojima.

*/
#ifndef MODIFIED_DISTANCES_HPP
#define MODIFIED_DISTANCES_HPP
#define MAXIMIZATION 1
#define MINIMIZATION 2
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

using namespace std;
class DistanceCalculation
{
   public:
     typedef vector<vector<double> > vvd;
      DistanceCalculation();
      double do_IGD_Plus(string File_ParetoFront, string File_Approximation);
      double do_GD_Plus(string File_ParetoFront, string File_Approximation);
      double InvertedGenerationalDistance(vvd &ParetoFront, vvd &Approximation);
      double GenerationalDistance(vvd &ParetoFront, vvd &Approximation);
      void NormalizeDataSets(vvd &ParetoFront, vvd &Approximation, int M);
      vector<double> ComputeMinDistances(vvd &SetA, vvd &SetB);
      double ModifiedDistance(vector<double> &x, vector<double> &y);
      void ReadData(string FileName, vvd &Set, int &cols);
      inline void setTypeProblem(int type){this->TypeProblem = type;}
      inline void setP(int p){ this->p = p;}
      inline void setNormalize(bool Normalize) { this->Normalize = Normalize;}

   private:
      double p=1;
      int TypeProblem = MINIMIZATION; 
      bool Normalize = false;
      
};
#endif

DistanceCalculation::DistanceCalculation()
{
}
void DistanceCalculation::ReadData(string FileName, vvd &Set, int &cols)
{
    FILE *instream;
    int retval = 1;
    char newline[2];
    double number;
    int error = 0;

    if (FileName.empty()) {
        instream = stdin;
        FileName = "<stdin>"; /* used to diagnose errors.  */
    }
    else if (NULL == (instream = fopen (FileName.c_str(),"r"))) {
   //     errprintf ("%s: %s\n", FileName, strerror (errno));
        exit (EXIT_FAILURE);
    }

   do
   {
	cols = 0;
	vector<double> row;

	/* Beginning of a row */
	do
	{
	   /* New column */
	   cols++;

		if (fscanf (instream, "%lf", &number) != 1) {
                    char buffer[64];
                    fscanf (instream, "%60[^ \t\r\n]", buffer);
                    exit (EXIT_FAILURE);
                }
		row.push_back(number);

	   fscanf(instream, "%*[ \t]");
	   retval = fscanf (instream, "%1[\r\n]", newline);
	
	}while(!retval);


	Set.push_back(row);
	/* skip over successive empty lines */
		do { 
		    if (!fscanf (instream, "%1[#]%*[^\n\r]", newline))
			fscanf (instream, "%*[ \t]");
		    retval = fscanf (instream, "%1[\r\n]", newline);
		} while (retval == 1);
   }while(retval != EOF);
}
double DistanceCalculation::do_GD_Plus(string File_ParetoFront, string File_Approximation)
{
   vvd ParetoFront, Approximation;
   int M;
   ReadData(File_ParetoFront, ParetoFront, M);
   ReadData(File_Approximation, Approximation, M);
   if(Normalize)
   NormalizeDataSets(ParetoFront, Approximation, M);
   return GenerationalDistance(ParetoFront, Approximation);
}
double DistanceCalculation::do_IGD_Plus(string File_ParetoFront, string File_Approximation)
{
   vvd ParetoFront, Approximation;
   int M;
   ReadData(File_ParetoFront, ParetoFront, M);
   ReadData(File_Approximation, Approximation, M);
   if(Normalize)
   NormalizeDataSets(ParetoFront, Approximation, M);
   return InvertedGenerationalDistance(ParetoFront, Approximation);
}
double DistanceCalculation::InvertedGenerationalDistance(vvd &ParetoFront, vvd &Approximation)
{
   //Compute minimal distances...
   vector<double> MinDistances(ParetoFront.size(), INFINITY);
   for(int i = 0; i < ParetoFront.size(); i++)
   {
	for(int j = 0; j < Approximation.size(); j++)
	{
		double distance = 0;

		for(int m = 0 ; m < ParetoFront[i].size(); m++)
		    if(TypeProblem == MAXIMIZATION)
		       distance += max( ParetoFront[i][m]-Approximation[j][m],0.0)*max( ParetoFront[i][m]-Approximation[j][m] , 0.0); // maximization max{z_i - a_i, 0}
		    else if(TypeProblem == MINIMIZATION)
		       distance += max( Approximation[j][m]-ParetoFront[i][m], 0.0)*max( Approximation[j][m]-ParetoFront[i][m] , 0.0); //minimization max{a_i - z_i, 0}
		   distance = sqrt(distance);
		MinDistances[i] = min(distance, MinDistances[i] );
	}
   }
   int N = ParetoFront.size();
   double Sum = 0;
   for(int i =0; i < N; i++)
	Sum += pow(MinDistances[i], this->p);
   return pow( Sum/N , 1.0/this->p  ); 
}
double DistanceCalculation::GenerationalDistance(vvd &ParetoFront, vvd &Approximation)
{

 //Compute minimal distances...
   vector<double> MinDistances(Approximation.size(), INFINITY);
   for(int i = 0; i < Approximation.size(); i++)
   {
	for(int j = 0; j < ParetoFront.size(); j++)
	{
		double distance = 0;

		for(int m = 0 ; m < Approximation[i].size(); m++)
		    if(TypeProblem == MAXIMIZATION)
		       distance += max( ParetoFront[j][m]-Approximation[i][m],0.0)*max( ParetoFront[j][m]-Approximation[i][m] , 0.0); // maximization max{z_i - a_i, 0}
		    else if(TypeProblem == MINIMIZATION)
		       distance += max( Approximation[i][m]-ParetoFront[j][m], 0.0)*max( Approximation[i][m]-ParetoFront[j][m] , 0.0); //minimization max{a_i - z_i, 0}
		   distance = sqrt(distance);
		MinDistances[i] = min(distance, MinDistances[i] );
	}
   }
   double Sum = 0;
   int N = Approximation.size();
   for(int i =0; i < N ; i++)
	Sum += pow(MinDistances[i], this->p);
   return pow( Sum/N , 1.0/this->p  );
}
void DistanceCalculation::NormalizeDataSets(vvd &ParetoFront, vvd &Approximation, int M)
{
  vector<double> Min(M, INFINITY), Max(M, -INFINITY);

   for(int  i = 0; i < ParetoFront.size(); i++)
   {
       for(int j = 0; j < M; j++)
	{
	   Min[j] = min(Min[j], ParetoFront[i][j]);
	   Max[j] = max(Max[j], ParetoFront[i][j]);
	}
   }
   for(int i = 0 ; i < ParetoFront.size(); i++)
   {
      for(int j = 0; j < M; j++)
	{
	   ParetoFront[i][j] = (ParetoFront[i][j] - Min[j]) / (Max[j] - Min[j]);
	}
   }
   for(int i = 0; i < Approximation.size(); i++)
   {
      for(int j = 0; j < M; j++)
	{
	  Approximation[i][j] = (Approximation[i][j]) / (Max[j] - Min[j]);
	}
   }
}
vector<double> DistanceCalculation::ComputeMinDistances(vvd &SetA, vvd &SetB)
{
   vector<double> MinDistances(SetA.size(), INFINITY);
   for(int i = 0; i < SetA.size(); i++)
   {
	for(int j = 0; j < SetB.size(); j++)
	{
		double distance = ModifiedDistance(SetA[i], SetB[j]);
		MinDistances[i] = min(distance, MinDistances[i] );
	}
   }
   return MinDistances;
}
/**
  The modified distance is considered as weakly pareto compliance...
**/
double DistanceCalculation::ModifiedDistance(vector<double> &x, vector<double> &y)
{
   double Sum = 0;
   for(int i = 0 ; i < x.size(); i++)
    if(TypeProblem == MAXIMIZATION)
       Sum += max( x[i]-y[i],0.0)*max( x[i]-y[i] , 0.0); //This is the modified distance 
    else if(TypeProblem == MINIMIZATION)
       Sum += max( y[i]-x[i],0.0)*max( y[i]-x[i] , 0.0); 
   return sqrt(Sum);
}

