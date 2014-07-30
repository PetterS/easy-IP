// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>

#include <easy-ip.h>

void main_program()
{
	using namespace std;

	IP ip;
	int n = 3;
	auto P = ip.add_boolean_cube(n*n, n*n, n*n);

	// Exactly one indicator equal to 1.
	cerr << "Indicator variables\n";
	for (int i = 0; i < n*n; ++i) {
		for (int j = 0; j < n*n; ++j) {
			Sum k_sum;
			for (int k = 0; k < n*n; ++k) {
				k_sum += P[i][j][k];
			}
			ip.add_constraint(k_sum == 1);

			// Advanced tip: One can use std::move to avoid a copy
			// here:
			//        ip.add_constraint(move(k_sum) == 1);
		}
	}

	// All rows have every number.
	cerr << "Row constraints\n";
	for (int i = 0; i < n*n; ++i) {
		for (int k = 0; k < n*n; ++k) {
			Sum row_k_sum;
			for (int j = 0; j < n*n; ++j) {
				row_k_sum += P[i][j][k];
			}
			ip.add_constraint(row_k_sum == 1);
		}
	}

	// All columns have every number.
	cerr << "Columns constraints\n";
	for (int j = 0; j < n*n; ++j) {
		for (int k = 0; k < n*n; ++k) {
			Sum col_k_sum;
			for (int i = 0; i < n*n; ++i) {
				col_k_sum += P[i][j][k];
			}
			ip.add_constraint(col_k_sum == 1);
		}
	}

	// The n*n subsquares have every number.
	cerr << "Square constraints\n";
	for (int i1 = 0; i1 < n; ++i1) {
		for (int j1 = 0; j1 < n; ++j1) {
			for (int k = 0; k < n*n; ++k) {
				Sum square_k_sum;
				for (int i2 = 0; i2 < n; ++i2) {
					for (int j2 = 0; j2 < n; ++j2) {
						square_k_sum += P[n*i1 + i2][n*j1 + j2][k];
					}
				}
				ip.add_constraint(square_k_sum == 1);
			}
		}
	}

	// A very difficult sudoku…
	const char* given[] = {"1**" "***" "7**",
	                       "**7" "1*9" "***",
	                       "68*" "*7*" "***",

	                       "**1" "*9*" "6**",
	                       "***" "3**" "*2*",
	                       "*4*" "***" "**3",

	                       "**8" "*6*" "1**",
	                       "5**" "***" "*4*",
	                       "***" "**2" "**5"};

	// http://www.mirror.co.uk/news/weird-news/worlds-hardest-sudoku-puzzle-ever-942299
	//const char* given[] = {"8**" "***" "***",
	//                       "**3" "6**" "***",
	//                       "*7*" "*9*" "2**",

	//                       "*5*" "**7" "***",
	//                       "***" "*45" "7**",
	//                       "***" "1**" "*3*",

	//                       "**1" "***" "*68",
	//                       "**8" "5**" "*1*",
	//                       "*9*" "***" "4**"};

	//const char* given[] = {"**" "**",
	//                       "**" "**",

	//                       "**" "12",
	//                       "**" "34"};

	cerr << "Preassignments\n";
	for (int i = 0; i < n*n; ++i) {
		for (int j = 0; j < n*n; ++j) {
			if (given[i][j] != '*') {
				int k = given[i][j] - '1';
				attest(0 <= k && k < n*n);
				ip.add_constraint( P[i][j][k] );
			}
		}
	}

	auto print_solution = [&P, n]()
	{
		cout << endl;
		for (int i = 0; i < n*n; ++i) {
			for (int j = 0; j < n*n; ++j) {
				for (int k = 0; k < n*n; ++k) {
					if (P[i][j][k].value()) {
						cout << k+1;
					}
				}
				if (j%n == n-1) {
					cout << ' ';
				}
			}
			cout << endl;
			if (i%n == n-1) {
				cout << endl;
			}
		}
	};

	#ifndef HAS_MINISAT
	#define HAS_MINISAT false
	#endif
	if (HAS_MINISAT) {
		cerr << "Solving with Minisat..." << endl;		
		ip.set_external_solver(IP::Minisat);
		attest(ip.solve());
		print_solution();
	}

	ip.set_external_solver(IP::Default);
	//ip.set_external_solver(IP::MOSEK);
	ip.solve();

	do {
		print_solution();
	} while (false /*ip.next_solution()*/);
}

int main()
{
	try {
		main_program();
	}
	catch (std::exception& err) {
		std::cerr << "EXCEPTION: " << err.what() << std::endl;
	}
}
