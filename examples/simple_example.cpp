// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
using namespace std;

#include <easy-ip.h>

int main()
{
	IP ip;
	auto x = ip.add_vector(5, IP::Binary);

	ip.add_objective( -2*x[0] - 5*x[1] - x[2] - 3*x[3] - x[4] );

	ip.add_constraint( x[0] + x[1] + x[2] + x[3] + x[4] <= 3 );

	ip.solve();

	for (auto xi : x) {
		cout << xi << endl;
	}
}

