/**
   Autor: Joel Chacón Castillo
   Date: 03/may/2017
   
   The IGD and GD plus versions are presented.

   This is based in the paper:
   "Modified Distance Calculation in Generational Distance and Inverted Generational Distance"
   by Hisao Ishibuchi, Hiroyuki Masuda, Yuki Tanigaki, and Yusuke Nojima.
**/
#include <iostream>
#include <string>
#include "ModifiedDistanceCalculation.hpp"

using namespace std;
 
string verbose()
{
   string message = "";
   message += "--h --help print this summary and exit\n";
   message += "--t --Type Max or Min (All objectives are setted default Min)\n";
   message += "--IGD compute the Inverted Generational Distance Plus\n";
   message += "--GD compute the Generational Distance Plus\n";
   message += "--P power of the formula (default 1 as the paper)\n";
   message += "Only can be used --IGD or --GD no both\n";
   message += "--r the filename of reference points\n";
   message += "--n Normalize the sets by the reference points provided\n";
   message += "--s the filename of data approximations, alternativelly this can be calculated as \" cat file | ./MPlus --GD --r reference.txt --t Max  \"  \n";
   message += "Example:\n  echo -e \"2 2\\n\" | ./MPlus --r reference.txt --ID";
   return message;
}
int main(int argc, char *argv[])
{

  if(argc<2)
         {
	    
	    cout << "Unknown Argument.."<<endl;
	    cout << verbose();
	    exit(0);
	 }

   int Type=MINIMIZATION; //max or min..
   double p = 1;
   bool Normalize = false;
   string Metric ="IGD", filenamePareto="", filenameSet="";
   for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--h")
			cout << verbose<<endl;
		else if( Terminal == "--t")
		{
			if(string(argv[++i]) == "Max")
			 Type = MAXIMIZATION;
			if(string(argv[++i]) == "Min")
			 Type = MINIMIZATION;
		}
		else if(Terminal == "--IGD")
			Metric = "IGD";
		else if(Terminal == "--GD")
			Metric = "GD";
		else if(Terminal == "--r" )
			filenamePareto = string(argv[++i]);
		else if(Terminal == "--s" )
			filenameSet = string(argv[++i]);
		else if(Terminal == "--P" )
			p = atof(argv[++i]);
		else if(Terminal == "--n" )
			Normalize = true;
		else
		{
		   cout << "Unknown Argument...";
		   exit(0);
		}

	}
  
   if(filenamePareto.empty()) {cout << "Reference points are needed"; exit(0);}

   DistanceCalculation ObjDistance; 
   ObjDistance.setP(p);
   ObjDistance.setTypeProblem(Type);
   ObjDistance.setNormalize(Normalize);
   if(Metric == "IGD")
      cout << ObjDistance.do_IGD_Plus(filenamePareto, filenameSet) << endl;
   else if(Metric =="GD")
      cout << ObjDistance.do_GD_Plus(filenamePareto, filenameSet) << endl;
   
   return 0;
}
