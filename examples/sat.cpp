// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <easy-ip.h>

int main()
{
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	auto u = ip.add_boolean();

	ip.add_constraint( ( x ||  y ||  z) &&
	                   ( z ||  w ||  u) &&
	                   (!x || !y || !w) &&
	                   (!x || !y || !u) &&
	                   (!x ||  z)         );

	ip.solve();

	cout << x << " " << y << " " << z << " "
	     << w << " " << u << endl;
}

