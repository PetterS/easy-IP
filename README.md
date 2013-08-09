[![Build Status](https://travis-ci.org/PetterS/easy-IP.png)](https://travis-ci.org/PetterS/easy-IP)

This is an experimental modelling library for integer (and linear) programming. Example programs look like this:

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

and

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

Requirements
------------
 * Cbc, see https://projects.coin-or.org/Cbc
 * C++11
 
Ubuntu
------
To install Cbc and its dependencies on Ubuntu, issue the following command:

	sudo apt-get install coinor-libcbc-dev coinor-libclp-dev coinor-libosi-dev coinor-libcgl-dev

CMake should then find and build everything.
