// Petter Strandmark 2013–2014.

#include <catch.hpp>

#include <easy-ip.h>

using namespace std;

#define REQUIRE_NUM_SOLUTIONS(ip, num) \
	REQUIRE(ip.solve()); \
	int num_solutions = 0; \
	do { \
		num_solutions++; \
	} while (ip.next_solution()); \
	REQUIRE(num_solutions == num);

auto create_soduku_IP(IP& ip, int n = 3)
	-> decltype(ip.add_boolean_cube(9, 9, 9))
{
	auto P = ip.add_boolean_cube(n*n, n*n, n*n);

	// Exactly one indicator equal to 1.
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
	return P;
}

TEST_CASE("minisat")
{
	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y == 1);
		REQUIRE(ip.solve());
		CHECK((x.value() + y.value() == 1));

		auto z = ip.add_boolean();
		ip.add_constraint(x + z == 1);
		ip.add_constraint(y + z == 1);
		CHECK( ! ip.solve());
		ip.add_constraint(x + y + z == 1);
		CHECK( ! ip.solve());
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		ip.add_objective(-x);
		CHECK(ip.solve());
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_objective(2*x + y);
		ip.add_constraint(x + y <= 1);
		ip.allow_ignoring_cost_function();
		CHECK(ip.solve());	
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y <= 1);
		CHECK(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 3);
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		auto u = ip.add_boolean();
		auto v = ip.add_boolean();
		ip.add_constraint(x + y + z + w + u + v <= 5);
		ip.add_constraint(x + y + z + w + u + v >= 2);
		CHECK(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 56);
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y == 2);
		CHECK(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 1);
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + 2*y == 1);
		CHECK(ip.solve());
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x - y == 1);
		CHECK(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 1);
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		ip.add_constraint(x + y - z - w == 0);
		CHECK(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 6);
	}

	{
		IP ip;
		ip.set_external_solver(IP::Minisat);
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		ip.add_constraint(x + y + z == 1);
		REQUIRE(ip.solve());
		int num_solutions = 0;
		do {
			num_solutions++;
		} while (ip.next_solution());
		CHECK(num_solutions == 3);
	}

	{
		IP ip;
		int n = 3;
		auto P = create_soduku_IP(ip, n);

		ip.set_external_solver(IP::Minisat);
		REQUIRE(ip.solve());
		 
		vector<vector<int>> solution(n*n);

		cout << endl;
		for (int i = 0; i < n*n; ++i) {
			for (int j = 0; j < n*n; ++j) {
				solution[i].emplace_back();

				for (int k = 0; k < n*n; ++k) {
					if (P[i][j][k].bool_value()) {
						solution[i][j] = k + 1;
					}
				}
			}
		}

		for (int i = 0; i < n*n; ++i) {
			int row_sum = 0;
			int col_sum = 0;
			for (int j = 0; j < n*n; ++j) {
				row_sum += solution[i][j];
				col_sum += solution[j][i];
			}
			CHECK(row_sum == 45);
			CHECK(col_sum == 45);
		}
	}
}

TEST_CASE("minisat-large")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	ip.allow_ignoring_cost_function();
	Sum sum = 0;
	for (int i = 0; i < 100; ++i) {
		auto x = ip.add_boolean(1);
		sum += x;
	}
	ip.add_constraint(sum == 50);
	CHECK(ip.solve());
	CHECK(sum.value() == 50);
}

TEST_CASE("sat-objective")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	ip.add_objective(4*x + y + 2*z + w);
	ip.add_constraint(x + y + z + w == 2);
	CHECK(ip.solve());
	CHECK(!x.bool_value());
	CHECK( y.bool_value());
	CHECK(!z.bool_value());
	CHECK( w.bool_value());
}

TEST_CASE("sat-constraints_coefs_larger_than_1")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(4*x - y - z >= 0);
	REQUIRE(ip.solve());

	REQUIRE_NUM_SOLUTIONS(ip, 5);
}

TEST_CASE("sat-constraints_coefs_zero")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	ip.add_constraint(0*x >= 0);

	REQUIRE_NUM_SOLUTIONS(ip, 2);
}

TEST_CASE("sat-constraints_fractional_coeff")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	ip.add_constraint(0.5*x >= 0);

	CHECK_THROWS(ip.solve());
}

TEST_CASE("sat-constraints_coefs_less_than_minus_1")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(8*x - 2*y - 2*z >= 0);

	REQUIRE_NUM_SOLUTIONS(ip, 5);
}

TEST_CASE("sat-constraints_coefs_negative_multiple_times")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(8*x - y - y - z - z >= 0);

	REQUIRE_NUM_SOLUTIONS(ip, 5);
}

TEST_CASE("sat-negative-objective")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	ip.add_objective(-4 * x - y - 2 * z - w);
	ip.add_constraint(x + y + z + w == 2);
	CHECK(ip.solve());
	CHECK(x.bool_value());
	CHECK(!y.bool_value());
	CHECK(z.bool_value());
	CHECK(!w.bool_value());
}

TEST_CASE("sat-objective-next_solution")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	auto u = ip.add_boolean();
	auto v = ip.add_boolean();
	auto obj = 4*x + y + 2*z + w + u + v;
	ip.add_objective(obj);
	ip.add_constraint(x + y + z + w + u + v == 2);
	
	REQUIRE_NUM_SOLUTIONS(ip, 6);
}

TEST_CASE("add_min_consequtive_constraints-2-true")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 3; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 2);
	ip.add_min_consequtive_constraints(2, x, true);

	// 1 1 0
	// 0 1 1
	// 1 0 1
	REQUIRE_NUM_SOLUTIONS(ip, 3);
}

TEST_CASE("add_min_max_consequtive_constraints-4")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_min_consequtive_constraints(4, x, false);
	ip.add_max_consequtive_constraints(4, x);

	int num_solutions = 0;
	REQUIRE(ip.solve(nullptr, true));
	do {
		CHECK(x_sum.value() <= 8);
		int consequtive = 0;
		for (auto& xx : x) {
			if (xx.value() > 0.5) {
				consequtive++;
			}
			else {
				consequtive = 0;
			}
			CHECK(consequtive <= 4);
		}
		num_solutions++;
	} while (ip.next_solution());
	CHECK(num_solutions == 20);
}

TEST_CASE("add_min_max_consequtive_constraints-5")
{
	IP ip;
	ip.set_external_solver(IP::Minisat);
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_min_consequtive_constraints(5, x, false);
	ip.add_max_consequtive_constraints(5, x);
	int num_solutions = 0;
	REQUIRE(ip.solve(nullptr, true));
	do {
		int consequtive = 0;
		for (auto& xx : x) {
			if (xx.value() > 0.5) {
			}
			else {
				consequtive = 0;
			}
			CHECK(consequtive <= 5);
		}
		num_solutions++;
	} while (ip.next_solution());
	CHECK(num_solutions == 12);
}
