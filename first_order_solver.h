// Petter Strandmark 2014
// petter.strandmark@gmail.com
//
// Implements the first-order solver described in [1], with the preconditioner of [2].
//
// [1] Chambolle, A., & Pock, T. (2011). A first-order primal-dual algorithm for convex
//     problems with applications to imaging. Journal of Mathematical Imaging and Vision,
//     40(1), 120-145.
//
// [2] Pock, T., & Chambolle, A. (2011, November). Diagonal preconditioning for first
//     order primal-dual algorithms in convex optimization. In Computer Vision (ICCV),
//     2011 IEEE International Conference on (pp. 1762-1769). IEEE.
//

#ifndef FIRST_ORDER_SOLVER_HEADER
#define FIRST_ORDER_SOLVER_HEADER

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <easy-ip.h>

struct FirstOrderOptions
{
	std::size_t maximum_iterations = 5000;
	// Print interval determines how often information should be logged
	// and how often convergence should be checked.
	std::size_t print_interval = 100;
	std::function<void(const std::string&)> log_function = nullptr;

	// Stops the solver if the relative changes in the primal and dual
	// variables are both lower than this value.
	// Also the solution is considered feasible if the maximum error
	// is less than 100 times this value.
	double tolerance = 1e-9;
};

/// Solves the linear program
///
///   minimize c·x
///   such that Ax = b,
///             l ≤ x ≤ u.
bool EASY_IP_API first_order_primal_dual_solve(Eigen::VectorXd* x,        /// Primal variables (in/out).
                                               Eigen::VectorXd* y,        /// Dual variables (in/out).
                                               const Eigen::VectorXd& c,  /// Objective function.
                                               const Eigen::VectorXd& lb, /// Lower bound on x.
                                               const Eigen::VectorXd& ub, /// Upper bound on x.
                                               const Eigen::SparseMatrix<double>& A,   /// Equality constraint matrix.
                                               const Eigen::VectorXd& b,  /// Right-hand side of constraints.
                                               const FirstOrderOptions& options);

bool EASY_IP_API first_order_admm_solve(Eigen::VectorXd* x,        /// Primal variables (in/out).
                                        const Eigen::VectorXd& c,  /// Objective function.
                                        const Eigen::VectorXd& lb, /// Lower bound on x.
                                        const Eigen::VectorXd& ub, /// Upper bound on x.
                                        const Eigen::SparseMatrix<double>& A,   /// Equality constraint matrix.
                                        const Eigen::VectorXd& b,  /// Right-hand side of constraints.
                                        const FirstOrderOptions& options);

class EASY_IP_API FirstOrderProblem
	: public IP
{
public:
	bool solve_first_order(const FirstOrderOptions& options);
protected:
	void convert_into_equality_constrained_problem();
private:
	void check_invariants();
};

#endif
