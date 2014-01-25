// Petter Strandmark 2014
// petter.strandmark@gmail.com
//
// [1] Chambolle, A., & Pock, T. (2011). A first-order primal-dual algorithm for convex
//     problems with applications to imaging. Journal of Mathematical Imaging and Vision,
//     40(1), 120-145.
//
// [2] Pock, T., & Chambolle, A. (2011, November). Diagonal preconditioning for first
//     order primal-dual algorithms in convex optimization. In Computer Vision (ICCV),
//     2011 IEEE International Conference on (pp. 1762-1769). IEEE.
//

#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <first_order_solver.h>

double get_feasibility_error(const Eigen::VectorXd& x,
                             Eigen::VectorXd* temp_storage,
                             const Eigen::SparseMatrix<double>& A,
                             const Eigen::VectorXd& b)
{
	using namespace std;
	const auto m = A.rows();

	double feasibility_error = 0;
	*temp_storage = A*x - b;
	for (ptrdiff_t i = 0; i < m; ++i) {
		feasibility_error = max(feasibility_error, abs((*temp_storage)(i)));
	}
	return feasibility_error;
}

// Will use x_prev as a temporary storage after
// examining it.
bool check_convergence_and_log(std::size_t iteration,
                               const Eigen::VectorXd& x,
                               Eigen::VectorXd* x_prev,
                               const Eigen::VectorXd& y,
                               const Eigen::VectorXd& y_prev,
                               const Eigen::VectorXd& c,
                               const Eigen::SparseMatrix<double>& A,
                               const Eigen::VectorXd& b,
                               const FirstOrderOptions& options)
{
	using namespace Eigen;
	using namespace std;

	const auto n = x.size();
	const auto m = y.size();

	bool is_converged = false;

	//cerr << "x      = " << (*x).transpose() << endl;
	//cerr << "Ax     = " << (A*(*x)).transpose() << endl;
	//cerr << "Ax - b = " << (A*(*x) - b).transpose() << endl;

	auto get_relative_change = [iteration](const Eigen::VectorXd& x, const Eigen::VectorXd& x_prev) -> double
	{
		double relative_change = (x - x_prev).norm() / (x.norm() + x_prev.norm());
		if (relative_change != relative_change) {
			// Both x and x_prev were null vectors.
			relative_change = 0;
		}
		if (iteration == 0) {
			relative_change = std::numeric_limits<double>::quiet_NaN();
		}
		return relative_change;
	};

	double relative_change_x = get_relative_change(x, *x_prev);
	double relative_change_y = get_relative_change(y, y_prev);
	double feasibility_error = get_feasibility_error(x, x_prev, A, b);

	if (relative_change_x < options.tolerance && relative_change_y < options.tolerance) {
		is_converged = true;
	}

	if (options.log_function) {
		ostringstream message;
		message << setw(9);
		if (iteration == 0) {
			message << "end";
		}
		else {
			message << iteration;
		}
		message << "   "
				<< setw(15) << setprecision(6) << scientific << c.dot(x) << " "
				<< setw(12) << setprecision(3) << scientific << relative_change_x << " "
				<< setw(12) << setprecision(3) << scientific << relative_change_y << " "
				<< setw(15) << setprecision(6) << scientific << feasibility_error;
		options.log_function(message.str());
	}

	return is_converged;
}

/// Solves the linear program
///
///   minimize c·x
///   such that Ax = b,
///             l ≤ x ≤ u.
bool first_order_primal_dual_solve(Eigen::VectorXd* x_ptr,    /// Primal variables (in/out).
                                   Eigen::VectorXd* y_ptr,    /// Dual variables (in/out).
                                   const Eigen::VectorXd& c,  /// Objective function.
                                   const Eigen::VectorXd& lb, /// Lower bound on x.
                                   const Eigen::VectorXd& ub, /// Upper bound on x.
                                   const Eigen::SparseMatrix<double>& A,   /// Equality constraint matrix.
                                   const Eigen::VectorXd& b,  /// Right-hand side of constraints.
                                   const FirstOrderOptions& options)
{
	using namespace Eigen;
	using namespace std;

	VectorXd& x = *x_ptr;
	VectorXd& y = *y_ptr;

	const auto n = x.size();
	const auto m = y.size();
	attest(c.size() == n);
	attest(A.rows() == m);
	attest(A.cols() == n);

	VectorXd x_prev(n);
	VectorXd y_prev(m);

	//cerr << "m=" << m << ", n=" << n << endl;
	//cerr << "A = \n" << A << endl;
	//cerr << "b = \n" << b  << endl;
	//cerr << "lb = \n" << lb.transpose() << endl;
	//cerr << "ub = \n" << ub.transpose() << endl;

	const SparseMatrix<double> AT = A.transpose();

	if (options.log_function) {
		options.log_function("   Iter         Objective     Rel. ch. x   Rel. ch. y   ||Ax - b||_inf");
		options.log_function("----------------------------------------------------------------------");
	}

	// Compute preconditioners as in eq. (10) from [2], with alpha = 1.
	VectorXd Tvec(n);
	VectorXd Svec(m);
	Tvec.setZero();
	Svec.setZero();
	for (int k = 0; k < A.outerSize(); ++k) {
		for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
			auto i = it.row();
			auto j = it.col();
			auto value = abs(it.value());
			Tvec(j) += value;
			Svec(i) += value;
		}
	}

	for (ptrdiff_t j = 0; j < n; ++j) {
		Tvec(j) = 1.0 / Tvec(j);
	}

	for (ptrdiff_t i = 0; i < m; ++i) {
		Svec(i) = 1.0 / Svec(i);
	}


	size_t iteration;
	for (iteration = 1; iteration <= options.maximum_iterations; ++iteration) {
		bool should_check_convergence = options.print_interval <= 1 || iteration % options.print_interval == 1;

		x_prev = x;
		if (should_check_convergence) {
			y_prev = y;
		}

		// See eq. (18) from [2].

		x = x - Tvec.asDiagonal() * ( AT*y + c);

		for (ptrdiff_t i = 0; i < n; ++i) {
			x(i) = max(lb(i), min(ub(i), x(i)));
		}
		
		y = y + Svec.asDiagonal() * (A*(2 * x - x_prev) - b);

		if (should_check_convergence) {
			if (check_convergence_and_log(iteration, x, &x_prev, y, y_prev, c, A, b, options)) {
				break;
			}
		}
	}

	check_convergence_and_log(0, x, &x_prev, y, y_prev, c, A, b, options);

	double feasibility_error = get_feasibility_error(x, &x_prev, A, b);
	return feasibility_error < 100*options.tolerance;
}


void FirstOrderProblem::check_invariants()
{
	auto n = get_cost().size();
	attest(get_var_lb().size() == n);
	attest(get_var_ub().size() == n);

	attest(get_rows().size() == get_values().size());
	attest(get_cols().size() == get_values().size());

	auto m = get_rhs_lower().size();
	attest(get_rhs_upper().size() == m);

	attest(get_integer_variables().empty());
}


void FirstOrderProblem::convert_into_equality_constrained_problem()
{
	auto m = get_rhs_lower().size();

	size_t constraints_added = 0;
	for (size_t i = 0; i < m; ++i) {
		auto& lb = get_rhs_lower()[i];
		auto& ub = get_rhs_upper()[i];
		if (lb != ub) {
			auto slack = add_variable(IP::Real);
			set_bounds(0.0, slack, std::numeric_limits<double>::infinity());
			get_rows().push_back(int(i));
			get_cols().push_back(int(get_variable_index(slack)));

			if (ub < 1e100) {
				attest(lb <= -1e100); // Cannot convert if both lb and ub exists.
				lb = ub;
				get_values().push_back(1.0);
			}
			else {
				attest(lb > -1e100);
				ub = lb;
				get_values().push_back(-1.0);
			}
			constraints_added++;
		}
	}

	check_invariants();
}


bool FirstOrderProblem::solve_first_order(const FirstOrderOptions& options)
{
	using namespace Eigen;
	using namespace std;

	auto n = get_cost().size();
	auto m = get_rhs_lower().size();
	check_invariants();

	// First, convert all inequality constraints to
	// equality constraints.
	convert_into_equality_constrained_problem();
	size_t constraints_added = get_cost().size() - n;
	n = get_cost().size();
	if (options.log_function && constraints_added > 0) {
		ostringstream sout;
		sout << constraints_added << " inequality constraints converted.";
		options.log_function(sout.str());
	}

	typedef unsigned int index;
	vector<Triplet<double, index>> sparse_indices;
	for (size_t ind = 0; ind < get_rows().size(); ++ind) {
		size_t i = get_rows()[ind];
		size_t j = get_cols()[ind];
		auto value = get_values()[ind];
		attest(i < m);
		attest(j < n);
		sparse_indices.emplace_back(index(i), index(j), value);
	}

	// Check that we only have equality constraints left.
	for (size_t i = 0; i < m; ++i) {
		auto lb = get_rhs_lower()[i];
		auto ub = get_rhs_upper()[i];
		attest(lb == ub);
	}

	SparseMatrix<double> A(static_cast<int>(m), static_cast<int>(n));
	A.setFromTriplets(sparse_indices.begin(), sparse_indices.end());

	Map<const VectorXd> lb(get_var_lb().data(), n);
	Map<const VectorXd> ub(get_var_ub().data(), n);
	Map<const VectorXd> b(get_rhs_upper().data(), m);
	Map<const VectorXd> c(get_cost().data(), n);
	VectorXd x(n); x.setZero();
	VectorXd y(m); y.setZero();

	bool feasible = first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, options);

	get_solution().resize(n);
	for (size_t j = 0; j < n; ++j) {
		get_solution()[j] = x(j);
	}

	return feasible;
}
