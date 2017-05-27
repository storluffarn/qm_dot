
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace boost::math;

// constants
double pi = 4*atan(1);
double hbar = 1.0545718e-34;

// derived constants
double reshbar = 1/hbar;

class oscillator
{
	class particle
	{
		double pox;
		double poy;

		friend oscillator;
	};

	//datamembers
	
	//dynamic variables
	double psi;
	double prob;
	double r;
	double energy;

	// constants
	double alpha;
	double lambda;
	
	particle p1;
	particle p2;
	vector <particle*> partvec = {&p1,&p2};

	// function members
	double diffpart(bool, bool);

	public:
	
	void dist();
	void calcpsi();
	void calcprob();
	void localen();
	int mc(double);

	// accessors
	void setpos(vector<double>);
	double geten(){return energy;}
	double getprob(){return prob;}

	// constructor
	oscillator(double ialpha, double ilambda)
		: alpha(ialpha), lambda(ilambda)
	{}
};
	
void oscillator::setpos(vector <double> pos)
{
	p1.pox = pos[0];
	p2.pox = pos[1];
	p1.poy = pos[2];
	p2.poy = pos[3];		
}

void oscillator::dist()
{ 
	r = sqrt(pow(p1.pox-p2.pox,2) + pow(p1.poy-p2.poy,2));

	//cout << r << endl;
}

void oscillator::calcpsi()
{
	//psi = exp(-pow(p1.pox,2)+pow(p2.pox,2)+pow(p2.pox,2)+pow(p2.poy,2))*lambda*sqrt(pow(p1.pox-p2.pox,2)+pow(p1.poy-p2.poy,2))/(1+alpha*sqrt(pow(p1.pox-p2.pox,2)+pow(p1.poy-p2.poy,2)));
	psi = exp(-0.5*(pow(p1.pox,2)+pow(p2.pox,2)+pow(p1.poy,2)+pow(p2.poy,2)) + lambda*r/(1+alpha*r));
}

void oscillator::calcprob()
{
	//psi = exp(-0.5*(pow(p1.pox,2)+pow(p2.pox,2)+pow(p1.poy,2)+pow(p2.poy,2)) + lambda*r/(1+alpha*r));
	prob = pow(psi,2);

	//cout << psi << endl << prob << endl;
}	
double oscillator::diffpart(bool state1, bool state2)
{
	double x1,x2,y1,y2;
	x1 = p1.pox; x2 = p2.pox; y1 = p1.poy; y2 = p2.poy;

	if (state1)
		{swap(x1,y1);swap(x2,y2);}
	if (state2)
		swap(x1,x2);

	//cout << x1 << " " << y1 << endl << x2 << " " << y2 << endl;

	double diff = exp(lambda*r/(1+alpha*r) - 0.5*(pow(x1,2)+pow(x2,2)+pow(y1,2)+pow(y2,2))) * (-1 + lambda/(pow(1+alpha*r,2)*r) - lambda*pow(x1-x2,2)/pow((1+alpha*r)*r,3)*(1+3*alpha*r) + pow(-x1+lambda*(x1-x2)/(pow(1 + alpha*r,2)*r),2));
	//double diff = exp(-0.5*lambda*r/(1+alpha*r)*(pow(p1.pox,2)+pow(p2.pox,2)+pow(p1.poy,2)+pow(p2.poy,2))) * (-1 + lambda/(pow(1+alpha*r,2)*r) - lambda*pow(p1.pox-p2.pox,2)/pow((1+alpha*r))*r,3)*(1+3*alpha*r) + pow(-x1+lambda*(p1.pox-p2.pox)/(part*pow(1+alpha*r,2)),2));

	return diff;
}

void oscillator::localen()
{
	double diff1 = diffpart(false,false);
	double diff2 = diffpart(false,true);
	double diff3 = diffpart(true,false);
	double diff4 = diffpart(true,true);
	
	//cout << "diffs: "  << diff1 << " " << diff2 << " " << diff3 << " " << diff4 << " " << diff1+diff2+diff3+diff4 << endl;
	//cout << "sqrs: " << pow(p1.pox,2) << " " << pow(p2.pox,2) << " " <<  pow(p1.poy,2) << " " << pow(p2.poy,2) << " " << pow(p1.pox,2) + pow(p2.pox,2) + pow(p1.poy,2) + pow(p1.poy,2) << endl;
	//cout<< "psi: " << psi << endl;

	energy = 0.5*(-(diff1 + diff2 + diff3 + diff4)/psi + pow(p1.pox,2) + pow(p1.poy,2) + pow(p2.pox,2) + pow(p2.poy,2) + 2*lambda/r);

	//cout << "energy: " << energy << endl;
}


int oscillator::mc(double maxdisp)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-maxdisp,maxdisp);

	uniform_int_distribution<> randpart(0,1);

	int k = randpart(gen);

	double oldpox = (*partvec[k]).pox;
	double oldpoy = (*partvec[k]).poy;
	double olp = prob;
	
	double dx = dis(gen);
	double dy = dis(gen);
	(*partvec[k]).pox += dx;
	(*partvec[k]).poy += dy;

	dist();			
	calcpsi();
	calcprob();
	
	//cout << "old prob: " << olp << endl << "new prob: " << prob << endl;
	
	prob /= olp;

	bernoulli_distribution accprob(prob);


	if (prob >= 1)
		return 1;
    else if(accprob(gen))
		return 1;
    else
		{(*partvec[k]).pox = oldpox; (*partvec[k]).poy = oldpoy; return 0;}
}

int main(int argc, char** argv)
{
	int runs = 100;
	int duration = 1000;
	double maxpos = 2;
	double en;
	double maxdisp = 0.01;
	double ascale = 10;
	double lscale = 1;
	double alpha = atof(argv[1]) / ascale;
	double lambda = atof(argv[2]) / lscale;

	argc == 3 ? cout << "alpha used: " << alpha << endl : cout << "too many or too few input" << endl;

	random_device rd; 
	mt19937 gen(rd());
	uniform_real_distribution<> ranpos(0,maxpos);

	oscillator harmonic(alpha, lambda);

	ostringstream namestream;
	namestream << "2denergies" << fixed << setprecision(2) << alpha << ".txt";
	string filename = namestream.str();

	ofstream fsen;
	fsen.open(filename);
	
	double account = 0;
	for (int l = 0; l < runs; l++) 
	{
		oscillator harmonic(alpha, lambda);

		vector <double> initpos = {-ranpos(gen),ranpos(gen),-ranpos(gen),ranpos(gen)};
		
		//vector <double> initpos (4);
		
		//for(auto& el : initpos)
		//	el = ranpos(gen);
		
		harmonic.setpos(initpos);
		harmonic.dist();
		harmonic.calcprob();

		//cout << "initprob: " <<  harmonic.getprob() << endl;

		double count = 0;

		for (int k = 0; k < duration; k++)
		{	
			int move = harmonic.mc(maxdisp);
			if (move == 1)
				count++;

			harmonic.dist();		// a little bit inefficientm should probably be handled in mc, don't need to caclulate distance if move was not accepted
			harmonic.calcpsi();
			harmonic.localen();
			
			en = harmonic.geten();
			
			//cout << en << endl;
			fsen << en << endl;
		}

		account += count;
	}	

	fsen.close();
	
	cout << "over all acceptance probability: " << account/(duration*runs) << endl;

}









