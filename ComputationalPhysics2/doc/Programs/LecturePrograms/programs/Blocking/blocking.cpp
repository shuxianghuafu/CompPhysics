// Blocking program
// Run as ./mcint_blocking.x <min_bloc_size> <max_block_size> <n_block_samples>

#include <cmath>
#include <iostream> // std io
#include <iomanip>  // std io
#include <fstream>  // file io
#include <sstream>  // string io
#include <sys/stat.h> // to get size of file
using namespace std;

// find mean of values in vals
double mean(double *vals, int n_vals){
  
  double m=0;
  for(int i=0; i<n_vals; i++){
    m+=vals[i];
  }

  return m/double(n_vals);
}

// calculate mean and variance of vals, results stored in res
void meanvar(double *vals, int n_vals, double *res){
  double m2=0, m=0, val;
  for(int i=0; i<n_vals; i++){
    val=vals[i];
    m+=val;
    m2+=val*val;
  }

  m  /= double(n_vals);
  m2 /= double(n_vals);

  res[0] = m;
  res[1] = m2-(m*m);

}

// find mean and variance of blocks of size block_size.
// mean and variance are stored in res
void blocking(double *vals, int n_vals, int block_size, double *res){

  // note: integer division will waste some values
  int n_blocks = n_vals/block_size;
  double* block_vals = new double[n_blocks];

  for(int i=0; i<n_blocks; i++){
    block_vals[i] = mean(vals+i*block_size, block_size);
  }
  meanvar(block_vals, n_blocks, res);
  delete block_vals;
}

int main (int nargs, char* args[])
{

  int min_block_size, max_block_size, n_block_samples;
  //   Read from screen a possible new vaue of n
  if (nargs > 3) {
    min_block_size  = atoi(args[1]);
    max_block_size  = atoi(args[2]);
    n_block_samples = atoi(args[3]); 
  }
  else{
    cerr << "usage: ./mcint_blocking.x <min_bloc_size> " 
	 << "<max_block_size> <n_block_samples>" << endl;
    exit(1);
  }

  // get file size using stat
  struct stat result;

  int n;
  if(stat("result.dat", &result) == 0){
    n = result.st_size;
  }
  else{
    cerr << "error in getting file size" << endl;
    exit(1);
  }

  // get all mc results from files
  double* mc_results = new double[n];
  ostringstream ost;
  ifstream infile;
  infile.open(ost.str().c_str(), ios::in | ios::binary);
  infile.read((char*)&(mc_results[n]),result.st_size);
  infile.close();
   // and compute total mean and variance
  double mean, sigma;
  double res[2];
  meanvar(mc_results, n, res);
  mean = res[0]; sigma= res[1];
  cout << "Value of mean = " <<  mean <<  endl;
  cout << "Value of variance = " << sigma << endl;
  cout << "Standard deviation = " <<  sqrt(sigma/(n-1.0)) << endl;

  // Open file for writing, writing results in formated output for plotting:
  ofstream outfile;
  outfile.open("blocking.dat", ios::out);
  outfile << setprecision(10);
  double* block_results = new double[n_block_samples];
  int block_size, block_step_length;

  block_step_length = (max_block_size-min_block_size)/n_block_samples;

  // loop over block sizes
  for(int i=0; i<n_block_samples; i++){
    block_size = min_block_size+i*block_step_length;
    blocking(mc_results, n, block_size, res);
    mean  = res[0];
    sigma = res[1];
    // formated output
    outfile << block_size << "\t" << sqrt(sigma/((n/block_size)-1.0))<< endl;
  }
  outfile.close();
  return 0;
}
