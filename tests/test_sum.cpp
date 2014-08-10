// Petter Strandmark 2013–2014.

#include <catch.hpp>

#define DISALLOW_COPY_CONSTRUCTION
#include <easy-ip.h>

using namespace std;

class IPTestFixture
{
public:
	IP ip;
	Variable x0; // == 0
	Variable x1; // == 1
	Variable x2;
	Variable x3;
	Variable x4;
	Variable x5; // == 5

	IPTestFixture()
		: ip(),
		  x0(ip.add_variable(IP::Integer)),
		  x1(ip.add_variable(IP::Integer)),
		  x2(ip.add_variable(IP::Integer)),
		  x3(ip.add_variable(IP::Integer)),
		  x4(ip.add_variable(IP::Integer)),
		  x5(ip.add_variable(IP::Integer))
	{
		ip.set_bounds(0, x0, 0);
		ip.set_bounds(1, x1, 1);
		ip.set_bounds(2, x2, 2);
		ip.set_bounds(3, x3, 3);
		ip.set_bounds(4, x4, 4);
		ip.set_bounds(5, x5, 5);
		ip.solve();
	}
};

TEST_CASE_METHOD(IPTestFixture, "test_variables")
{
	CHECK(ip.get_solution(x0) == 0);
	CHECK(ip.get_solution(x1) == 1);
	CHECK(ip.get_solution(x2) == 2);
	CHECK(ip.get_solution(x3) == 3);
	CHECK(ip.get_solution(x4) == 4);
	CHECK(ip.get_solution(x5) == 5);
}

TEST_CASE_METHOD(IPTestFixture, "no_coefficients")
{
	CHECK(ip.get_solution(x0 + x0 + x0) == 0);
	CHECK(ip.get_solution(x1 + x0 + x1) == 2);
	CHECK(ip.get_solution(x1 + x2 + x3) == 6);
	CHECK(ip.get_solution(x0 + x0 + x5 + x0 + x0 + x1 + x1) == 7);

	Sum s = x1 + x1;
	CHECK(ip.get_solution((x1+x0) + (x1+x0)) == 2);
	CHECK(ip.get_solution(s + (x1+x0)) == 3);
	CHECK(ip.get_solution((x1+x0) + s) == 3);
	CHECK(ip.get_solution(s + s) == 4);
}

TEST_CASE_METHOD(IPTestFixture, "coefficients")
{
	CHECK(ip.get_solution(1.0*x0) == 0);
	CHECK(ip.get_solution(0.0*x0) == 0);
	CHECK(ip.get_solution(1.0*x2) == 2);
	CHECK(ip.get_solution(5.0*x4) == 20);
	CHECK(ip.get_solution(-2.0*x5) == -10);
}

TEST_CASE_METHOD(IPTestFixture, "coefficients_with_constants")
{
	CHECK(ip.get_solution(1.0*x0 + 1) == 1);
	CHECK(ip.get_solution(1 + 0.0*x0 + 1) == 2);
	CHECK(ip.get_solution(5.0*x4 - 1) == 19);
	CHECK(ip.get_solution(-2.0*x5 + 10) == 0);
}

TEST_CASE_METHOD(IPTestFixture, "linear_combinations")
{
	CHECK(ip.get_solution(1.0*x0 + 3.0*x2) == 6);
	CHECK(ip.get_solution(x0 + x0 + 4.0*x5 + x0 + x0 + 3.0*x1 + x1) == 24);
}

TEST_CASE_METHOD(IPTestFixture, "unary_minus")
{
	CHECK(ip.get_solution(4.0*x5 + (-(3.0*x1)) ) == 17);
}

TEST_CASE_METHOD(IPTestFixture, "unary_minus_no_coefficient")
{
	CHECK(ip.get_solution(4.0*x5 + (-x1)) == 19);
	CHECK(ip.get_solution(4.0*x5 + (-(x1 + x0))) == 19);
}

TEST_CASE_METHOD(IPTestFixture, "binary_minus")
{
	CHECK(ip.get_solution(4.0*x5 - 3.0*x1) == 17);
	CHECK(ip.get_solution(4.0*x5 - x1) == 19);

	Sum s = x2 + x0;
	CHECK(ip.get_solution(s - (x1+x0)) == 1);
	CHECK(ip.get_solution((x1+x0) - s) == -1);
	CHECK(ip.get_solution(s - s) == 0);
}

TEST_CASE_METHOD(IPTestFixture, "scalar_multiplication")
{
	CHECK(ip.get_solution( 0 * (x1 + x2 + x3) ) == 0);
	CHECK(ip.get_solution( (x1 + x2 + x3) * 0 ) == 0);
}
