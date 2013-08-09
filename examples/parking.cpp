// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
#include <stdexcept>
#include <vector>
using std::vector;
using std::size_t;

#include <easy-ip.h>

class MyIP : 
	public IP
{
public:
	void add_connected_constraint(const vector<vector<BooleanVariable>>& grid)
	{
		auto m = grid.size();
		if (m == 0) return;
		auto n = grid[0].size();

		// Pad grid with zeroes to be able to add constraints everywhere.
		vector<vector<BooleanVariable>> padded_grid;
		for (int i = 0; i < m + 2; ++i) {
			padded_grid.push_back(vector<BooleanVariable>());
			for (int j = 0; j < n + 2; ++j) {
				if (i == 0 || i == m + 1 || j == 0 || j == n + 1) {
					auto zero_var = add_boolean();
					set_bounds(0, zero_var, 0);
					padded_grid.back().push_back(zero_var);
				}
				else {
					padded_grid.back().push_back(grid[i-1][j-1]);
				}
			}
		}

		Sum turn_sum;

		auto add_triple_constraints = 
			[&]
			(const BooleanVariable& x1, const BooleanVariable& x2, const BooleanVariable& y)
			{
				auto v1 = add_boolean();
				// v1  <=>  y && !x1 && !x2.
				// <=
				add_constraint( y + (1-x1) + (1-x2) <= 2 + v1 ); 
				// =>
				add_constraint( implication(v1, !x1) );
				add_constraint( implication(v1, !x2) );
				add_constraint( implication(v1, y) );
				turn_sum += v1;

				auto v2 = add_boolean();
				// v2  <=>  !y && x1 && x2.
				// <=
				add_constraint( (1-y) + x1 + x2 <= 2 + v2 );
				// =>
				add_constraint( implication(v2, x1) );
				add_constraint( implication(v2, x2) );
				add_constraint( implication(v2, !y) );
				turn_sum -= v2;
			};

		for (int i = 0; i < m + 1; ++i) {
			for (int j = 0; j < n + 1; ++j) {

				//  X
				//  YX
				{
					auto& x1 = padded_grid[i][j];
					auto& x2 = padded_grid[i+1][j+1];
					auto& y  = padded_grid[i][j+1];
					add_triple_constraints(x1, x2, y);
				}

				// XY
				//  X
				{
					auto& x1 = padded_grid[i][j];
					auto& x2 = padded_grid[i+1][j+1];
					auto& y  = padded_grid[i+1][j];
					add_triple_constraints(x1, x2, y);
				}

				// YX
				// X
				{
					auto& x1 = padded_grid[i+1][j];
					auto& x2 = padded_grid[i][j+1];
					auto& y  = padded_grid[i][j];
					add_triple_constraints(x1, x2, y);
				}

				//  X
				// XY
				{
					auto& x1 = padded_grid[i+1][j];
					auto& x2 = padded_grid[i][j+1];
					auto& y  = padded_grid[i+1][j+1];
					add_triple_constraints(x1, x2, y);
				}
			}
		}

		add_constraint(4, turn_sum, 4);
	}
};

int main_program()
{
	using namespace std;

	int n = 6;

	MyIP ip;

	// Parking spots. We want to maximize the number of these.
	auto P = ip.add_boolean_grid(n, n, -1.0);
	// Transporation to parking spots.
	auto T = ip.add_boolean_grid(n, n);

	// Entrance to the parking lot.
	ip.add_constraint( T[0][1] );

	// T[0][0] == 1
	//. P P P P P
	//. . . . . .
	//. P P P P P
	//. P P P P P
	//. . . . . .
	//P P P P P P
	// 21

	// T[0][1] == 1
	//P . P P P P
	//P . . . . .
	//P . P P P P
	//P . P P P P
	//P . . . . .
	//  P P P P P
	// 22

	// T[0][2] == 1
	//P . . . . P
	//P . P P . P
	//P . P P . P
	//P . P P . P
	//P . P P . P
	//P . P P . P
	// 22


	// We can not have both parking and transportation.
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			ip.add_constraint(P[i][j] + T[i][j] <= 1);
		}
	}

	// Every parking lot needs to be adjacent to a transport square.
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Sum neighborhood;
			if (i > 0) {
				neighborhood += T[i-1][j];
			}
			if (i < n - 1) {
				neighborhood += T[i+1][j];
			}
			if (j > 0) {
				neighborhood += T[i][j-1];
			}
			if (j < n - 1) {
				neighborhood += T[i][j+1];
			}

			ip.add_constraint( P[i][j] <= neighborhood );
			ip.add_constraint( T[i][j] <= neighborhood + P[i][j] );
		}
	}

	// The hardest part: the transport squares need to be
	// connected to each other.
	ip.add_connected_constraint(T);

	auto print_solution = [&] () 
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				auto p = ip.get_solution(P[i][j]);
				auto t = ip.get_solution(T[i][j]);
				if (p && t) {
					// Infeasible.
					cout << '#';
				}
				else if (!p && !t) {
					cout << ' ';
				}
				else if (p) {
					cout << 'P';
				}
				else {
					cout << '.';
				}
				cout << ' ';
			}
			cout << endl;
		}
		cout << endl;
	};

	//ip.set_external_solver(IP::MOSEK);
	attest(ip.solve());

	do {
		print_solution();
	} while (ip.next_solution());


	return 0;
}

int main()
{
	try {
		return main_program();
	}
	catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
	}
}
