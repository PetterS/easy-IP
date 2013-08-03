// Petter Strandmark 2013.

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <easy-ip.h>

using namespace std;

class IPTestFixture
{
public:
	IP ip;
	BooleanVariable x0; // == 0
	BooleanVariable x1; // == 1
	BooleanVariable y0; // == 0
	BooleanVariable y1; // == 1

	IPTestFixture()
		: x0(ip.add_boolean()),
		  x1(ip.add_boolean()),
		  y0(ip.add_boolean()),
		  y1(ip.add_boolean())
	{
		ip.set_bounds(0, x0, 0);
		ip.set_bounds(1, x1, 1);
		ip.set_bounds(0, y0, 0);
		ip.set_bounds(1, y1, 1);
		ip.solve();
	}
};

TEST_CASE_METHOD(IPTestFixture, "test_variables")
{
	CHECK(ip.get_solution(x0) == false);
	CHECK(ip.get_solution(x1) == true);
	CHECK(ip.get_solution(y0) == false);
	CHECK(ip.get_solution(y1) == true);
}

TEST_CASE_METHOD(IPTestFixture, "test_not")
{
	CHECK(ip.get_solution(!x0) == true);
	CHECK(ip.get_solution(!x1) == false);
}

TEST_CASE_METHOD(IPTestFixture, "test_or")
{
	CHECK(ip.get_solution(x1 || y1) == true);
	CHECK(ip.get_solution(x1 || x1) == true);
	CHECK(ip.get_solution(x1 || x1 || y1) == true);
	CHECK(ip.get_solution(x1 || y1 || x1) == true);
	CHECK(ip.get_solution(y1 || x1 || x1) == true);
	CHECK(ip.get_solution(x0 || x1 || y0) == true);
	CHECK(ip.get_solution(x1 || x0 || y0) == true);
	CHECK(ip.get_solution(y0 || x0 || y1) == true);

	CHECK(ip.get_solution(x0 || y0) == false);
	CHECK(ip.get_solution(y0 || x0 || y0) == false);
	CHECK(ip.get_solution(x0 || x0 || y0) == false);
}
