
// redefining random generatorm, should be fixed for the sace of aesthetics

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>
#include <sstream>
#include <string>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace boost::math;

// constants
double pi = 4*atan(1);
double hbar = 1.0545718e-34;

// derived constants
double reshbar = 1/hbar;

class particle
{
	// data members
	double pos;
	double energy;
	double prob;
	double alpha;

	// derived data members
	//double omega;

	// some constants
	//double c1 = pow(mass*omega/(pi*hbar),0.25);
	//double c2 = sqrt(mass*omega/hbar);

	// function members

	public:
	
	// public functions
	void dice();
	int tokyo(double);
	void redbull();
	
	// accessors 
	void setpos (double x) {pos = x;};
	double getpos() {return pos;};
	double geten() {return energy;};

	// constructor
	particle(double ipos, double ialpha)
		: pos(ipos), alpha(ialpha)
	{
	//	omega = sqrt(spring/mass);
	}
};

void particle::dice()
{
	//double psi = 1/sqrt(pow(2,enlevel)*factorial<double>(enlevel))*c1*exp(-0.5*mass*omega*pow(pos,2)*reshbar)*hermite(enlevel,c2*pos);
	double psi = exp(-alpha*pow(pos,2));
	prob = pow(abs(psi),2);

	//cout << prob << endl;
}

void particle::redbull()
{
	energy = alpha + pow(pos,2)*(0.5-2*pow(alpha,2));
}

int particle::tokyo(double maxdisp)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-maxdisp,maxdisp);

	double olx = pos;
	double olp = prob;
	
	double dx = dis(gen);
	pos += dx;

	dice();
	prob /= olp;

	//cout << olp << " "  <<  prob << endl;

	bernoulli_distribution accprob(prob);

	if (prob >= 1)
		return 1;
    else if(accprob(gen))
		return 1;
    else
		{pos = olx; return 0;} 
}

double shiny(double x1, double x2, double x3, double tol, double* suchfun(double))
{
	double y1 = *suchfun(x1);
	double y2 = *suchfun(x2);
	double y3 = *suchfun(x3);

	if (~(y2 < y1 && y2 < y3))
		{cout << "bad initial conditions, no minimum in range" << endl; return 0;}
	
	double x4 = x1 + (x3 - x2);
	double y4 = *suchfun(x4);

	if (abs(x3-x1) < tol*(fabs(x2)+fabs(x4)))
		return min(y2,y4);
	else if (y2 > y4)
		return shiny(x2,x4,x3,tol,suchfun);
	else if (y4 > y2)
		return shiny(x1,x2,x4,tol,suchfun);
	else 
		{cout << "bad comparison, unable to determine next move" << endl; return 0;}
}

double integrate (double x1, double x2, int n, function<double(double)> f)
{
	vector<double> x;  
	vector<double> y;  
	double h=(x2-x1)/n;    

	for (int i=0;i<=n;i++)		
	{    
		x.push_back(x1+i*h);
		y.push_back(f(x[i]));
	}   

	double sum = 0;
	for (auto& el : y)    
	{   
		sum += h*el;
	}   

   return h/2.0*(y[0]+y[n])+sum;    
}

int main(int argc, char** argv)
{
	int runs = 100;
	int duration = 1000;
	double maxpos = 2;
	double maxdisp = 0.01;
	double en, anen;
	int steps = 1e5;
	double scale = 100;
	double alpha = atof(argv[1]) / scale;
	
	argc == 2 ? cout << "alpha used: " << alpha << endl : cout << "too many or too few input" << endl;

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> ranpos(-maxpos,maxpos);
	
	ostringstream namestream;
	namestream << "1denergies" << fixed << setprecision(2) << alpha << ".txt";
	string filename = namestream.str();
	
	ofstream fsen;
	fsen.open(filename);

	double account = 0;
	for (int l = 0; l < runs; l++)
	{	
		double pos = ranpos(gen);
		particle jojo(pos, alpha);
		jojo.dice();
		
		double count = 0;
	
		for (int k = 0; k < duration; k++)
		{
			int move = jojo.tokyo(maxdisp);
			if (move == 1)
				count++;
	
			jojo.redbull();
			en = jojo.geten();
	
			fsen << en << endl;
		}

		account += count;
	}
	
	fsen.close();
	
	cout << "over all acceptance probability: " << account / (duration*runs) << endl;
		
	function<double(double)> f = [alpha](double z){return exp(-alpha*pow(z,2)) * 0.5*exp(-alpha*pow(z,2))*(2*alpha+pow(z,2)-4*pow(alpha*z,2)) / pow(exp(-alpha*pow(z,2)),2);}; //This si wrong, gives incorrect values, also there should be two integrals, not one...!

	anen = integrate(-maxdisp,maxdisp,steps,f);
	cout << "analytic mean: " << anen << endl;

	//fsen << en << " " << anen << endl;
		
	
	// add loop over alpha
}









