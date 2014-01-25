// Petter Strandmark 2014.

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <first_order_solver.h>

TEST_CASE("Strange_convergence")
{
	using namespace Eigen;
	using namespace std;

	int m = 4;
	int n = 8;
	MatrixXd Adense(m, n);
	Adense.row(0) << 1, 5,  6,  1, -5, 1, 6, 0;
	Adense.row(1) << 2, 5,  5, -1,  5, 1, 6, 0;
	Adense.row(2) << 3, 5, -4, -1,  5, 2, 0, 9;
	Adense.row(3) << 4, 5,  1,  1, -5, 3, 0, 9;
	SparseMatrix<double> A;
	A = Adense.sparseView();

	VectorXd x(n);
	x << 1, 1, 0, 0, 1, 1, 0, 0;

	VectorXd b(m);
	b = A*x;
	x.setZero();

	VectorXd c(n);
	c << 1, -1, 2, -3, -5, 7, 9, -6;

	VectorXd y(m);
	y.setZero();

	VectorXd lb(n);
	lb.setZero();

	VectorXd ub(n);
	ub.setConstant(1.0);

	FirstOrderOptions options;

	// From my first implementation, which obtains similar
	// convergence plots as the paper.
	vector<double> ground_truth = 
		{-5.550000,
		-11.359921,
		-14.174585,
		-14.017215,
		-13.926363,
		-13.702674,
		-13.234735,
		-12.690176,
		-12.339068,
		-12.289323,
		-12.469620,
		-12.722200,
		-12.900454,
		-12.898698,
		-12.640000,
		-12.260535,
		-11.814795,
		-11.328464,
		-10.804136,
		-10.233576};

	for (size_t i = 0; i < ground_truth.size(); ++i) {
		options.maximum_iterations = i + 1;
		x.setZero();
		y.setZero();
		first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, options);
		CAPTURE(i);
		CAPTURE(c.dot(x));
		CAPTURE(ground_truth[i]);
		REQUIRE(abs(c.dot(x) - ground_truth[i]) < 1e-6);
	}

	options.maximum_iterations = 20000;
	options.print_interval = 1000;
	options.log_function = [](const string& str) { cerr << str << endl; };
	x.setZero();
	y.setZero();
	first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, options);
	
	cerr << "x     = " << x.transpose() << endl;
	cerr << "<c|x> = " << c.dot(x) << endl;
	
	CHECK( abs(c.dot(x) - 2.0) <= 1e-6);
}

TEST_CASE("FirstOrderProblem/tiny")
{
	using namespace std;

	FirstOrderProblem problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);
	
	problem.set_bounds(0.0, x, 1.0);
	problem.set_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y == 1.0);
	problem.add_objective(-2.0*x - 1.0*y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 1;
	options.log_function = [](const string& str) { cerr << str << endl; };
	REQUIRE(problem.solve_first_order(options));
	CHECK(x.value() >= 0.999);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/tiny/inequality1")
{
	using namespace std;

	FirstOrderProblem problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.set_bounds(0.0, x, 1.0);
	problem.set_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y <= 1.0);
	problem.add_objective(-2.0*x - 1.0*y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 1;
	options.log_function = [](const string& str) { cerr << str << endl; };
	REQUIRE(problem.solve_first_order(options));
	CHECK(x.value() >= 0.999);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/tiny/inequality2")
{
	using namespace std;

	FirstOrderProblem problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.set_bounds(0.0, x, 1.0);
	problem.set_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y <= 1.0);
	problem.add_objective(2.0*x + 1.0*y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 10;
	options.log_function = [](const string& str) { cerr << str << endl; };
	REQUIRE(problem.solve_first_order(options));
	CHECK(x.value() <= 0.001);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/no_integer_variables")
{
	FirstOrderProblem problem;
	auto x1 = problem.add_variable(IP::Integer);
	FirstOrderOptions options;
	CHECK_THROWS(problem.solve_first_order(options));
}

TEST_CASE("FirstOrderProblem/inequalities")
{
	using namespace std;

	FirstOrderProblem problem;
	auto x1 = problem.add_variable(IP::Real);
	auto x2 = problem.add_variable(IP::Real);
	auto x3 = problem.add_variable(IP::Real);
	auto x4 = problem.add_variable(IP::Real);
	auto x5 = problem.add_variable(IP::Real);

	problem.set_bounds(-10.0, x1, 10.0);
	problem.set_bounds(-10.0, x2, 10.0);
	problem.set_bounds(-10.0, x3, 10.0);
	problem.set_bounds(-10.0, x4, 10.0);
	problem.set_bounds(-10.0, x5, 10.0);

	problem.add_constraint(x1   + 2*x5 <= 5.0);
	problem.add_constraint(x2   + x4   <= 7.0);
	problem.add_constraint(2*x3 + x4   <= 4.0);
	problem.add_constraint(x1   + x4   <= 5.0);
	problem.add_constraint(x3   + 2*x5 <= 4.0);
	problem.add_constraint(x1   + x2   <= 6.0);

	auto obj = -2.0*x1 - 1.0*x2 - 3.0*x3 - 1.0*x4 - 2.0*x5;
	problem.add_objective(obj);

	problem.solve();
	auto optimal_objective = obj.value();
	cerr << "Optimal value: " << optimal_objective << endl;

	FirstOrderOptions options;
	options.maximum_iterations = 20000;
	options.print_interval = 10;
	options.log_function = [](const string& str) { cerr << str << endl; };

	REQUIRE(problem.solve_first_order(options));
	CHECK(abs(obj.value() - optimal_objective) <= 1e-4);
}

TEST_CASE("FirstOrderProblem/infeasible")
{
	using namespace std;

	FirstOrderProblem problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.set_bounds(0.0, x, 1.0);
	problem.set_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y == 1.0);
	problem.add_constraint(x + y == 2.0);
	problem.add_objective(-2.0*x - 1.0*y);

	FirstOrderOptions options;
	options.maximum_iterations = 2000;
	options.print_interval = 100;
	options.log_function = [](const string& str) { cerr << str << endl; };
	REQUIRE_FALSE(problem.solve_first_order(options));
}
