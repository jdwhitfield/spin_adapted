#include <random>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	int d=5; // in future version read from commandline
	if(argc>1)
	{
		d = atoi(argv[1]);
	}

	std::fstream file("matrix.dat", std::ofstream::out | std::ofstream::trunc);

	std::mt19937 rng;
	rng.seed(std::random_device()());
	std::normal_distribution<double> normal_dist(0, 1); 

	file << d << std::endl;

	for(int i=0; i<d*(d+1)/2; i++)
		file << normal_dist(rng) <<  "  ";

}
