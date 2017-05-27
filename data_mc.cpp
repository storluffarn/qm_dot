
// compile with -lboost_filesystem -lboost_system

#define BOOST_FILESYSTEM_VERSION 3
#define BOOST_FILESYSTEM_NO_DEPRECATED 
#include <boost/filesystem.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iterator>

using namespace std;
namespace fs = ::boost::filesystem;

class dataclass                                       // data to be analysed
{
	// private section
	
	// data members
	
	string filename;                                  // name of file to read
	vector <double> data;                             // data read
	int size;                                         // numbers written/read

	// function members
	vector<double> readdata(string);                 // function for reading data
	void writedata(vector<double>*,string);          // function for writing data
	double average;
		
	// public section

	public:

    //function members
	
	vector <double> calcerror ();
	
	// inline functions
	
	void calcavg()			            // calculates mean of data
	{
		average = accumulate(data.begin(),data.end(),0.0)/data.size();
	}	

	// accessors
	
	vector <double> getdata(){return data;}
	double getavg(){return average;}
	
	// constructor

	dataclass(string s)
        : filename(s)                   // setting initial values
    {
        data = readdata(filename);		// read data from file
        size = data.size();				// save size of data
		calcavg();
    }
};


vector <double> dataclass::readdata(string filename)           // reads data from file
{
	ifstream readdata(filename);	      // read from start to end to vector
	istream_iterator<double> start(readdata), end;
	vector<double> data(start, end);

	cout << "Read " << filename << " " << data.size() << " numbers" << endl;    // inform user of read date

	return data;
}

void dataclass::writedata(vector <double>* writedata, string filename)  // wrotes data to file
{
	ofstream filestream;
	filestream.open(filename);

	for(unsigned int k = 0; k < writedata->size(); k++)         // write all data to file
	{
		filestream << (*writedata)[k] << endl;
	}

	filestream.close();

	cout << "Wrote " << writedata->size() << " numbers to " << filename << endl;    // inform user of what happens
}
			   
vector <double> dataclass::calcerror()
{
	vector <double> shifted (data.size());

	for (auto& el: shifted)
	{   
		auto k = &el - &shifted[0];
		el = pow(data[k] - average,2);
	}
	
	double variance = accumulate(shifted.begin(),shifted.end(),0.0)/shifted.size();
	double stdev = sqrt(variance);
	double error = stdev/sqrt(data.size()-1);

	vector <double> retvals = {variance,stdev,error};

	return retvals;
}

//
// usefull stuff, searching directory and sub directories for files of a specific type
//
//vector<fs::path> get_all(const fs::path& root, const string& ext)	// searches directory and returns all filenames of given extention
//{
//	vector<fs::path> ret;
//	if(!fs::exists(root) || !fs::is_directory(root)) 
//		{cout << "no such file or directory" << endl; return ret;}
//
//	fs::directory_iterator it(root);
//	fs::directory_iterator endit;
//	
//	while(it != endit)
//	{ 
//		if(fs::is_regular_file(*it) && it->path().extension() == ext) 
//			ret.push_back(it->path().filename());
//		
//		++it;
//		cout << it->path().string() << endl;
//	} 
//
//	return ret;
//}

vector<fs::path> get_all(const fs::path& root, const string& ext)	// searches directory and returns all filenames of given extention
{
	vector<fs::path> ret;
	if(!fs::exists(root) || !fs::is_directory(root)) 
		{cout << "no such file or directory" << endl; return ret;}

	fs::directory_iterator it(root);
	
	for	(fs::directory_entry& file : it)
	{ 
		char test = '2';
		if(fs::is_regular_file(*it) && it->path().extension() == ext && file.path().filename().string().at(0) == test) 
			ret.push_back(file.path().filename());
	
		//cout << file.path().filename()<< endl;
	} 

	return ret;
}


int main()
{
	vector <fs::path> txtfiles = get_all("/home/storluffarn/phd-stuff/courses/comp_phys", ".txt");
	sort(txtfiles.begin(),txtfiles.end());

	vector <dataclass> alldata;
	
	for (auto& el : txtfiles)
		{alldata.push_back(el.string()); cout << el.string() << endl;}

	ofstream output;
	output.open("analysis.csv");
	
	for (auto& el : alldata)
		output << el.getavg() << "," << el.calcerror()[0] << endl;

	output.close();	
}























