#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Code/PersCoH.hh>
//For parallelization
#include <Code/par_for.hh>


int main()
{
    //load dataset, e.g. from .csv file, to format:
    std::vector<std::vector<double> > data;

    //specify max dimension for homology
    int n;
    //specify local radius for local homology computation.
    long double rad;
    //Note that rad/2 will determine the maximal scale of the Vietroris-Rips complexes involved

    pl::async_par_for(0, data.size(), [&](unsigned i)
    {
      for(int i = 0; i < data.size(); i++)
      {
          vector<long double> z = data[i];
  
          vector < vector< pair<long double, long double > > > LCP
          = BndLocGCopairings(data,n,z,rad,false);

          //or alternatively:
          vector < vector< pair<long double, long double > > > LCP
          = CLocGCopairings(data,n,z,rad,false);

          //one may proceed as follows to safe the results as a csv file:
          ofstream myfile;
          string filename;
          myfile.open(path_out + filename);
  
          for (int j = 0; j < LCP.size(); j++)
          {   
              if(LCP[j].size())
              {
                  std::cout << "Total amount of local "<< j << "-pairs: " << LCP[j].size() << std::endl;
                  for (int l = 0; l < LCP[j].size(); l++)
                  {
                      myfile << LCP[j][l].first << "," << LCP[j][l].second << "," << j << std::endl;
                  }
              }
          }
          myfile.close();
      }
    });
    //change ending to },false); for an ordinary for-loop

    //Now, we can use the results (safed as .csv files to not overload the RAM) to determine the Phi_PL values.
    //See LocCoH_plot.py for a visualization script


}
