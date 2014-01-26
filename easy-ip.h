// Petter Strandmark 2013
// petter.strandmark@gmail.com
//
/// This header provides a simple modelling language in C++.
///
/// Example:
///
///		IP ip;
///		auto x = ip.add_variable();
///		auto y = ip.add_variable();
///
///		ip.add_objective(3.0*x + y);
///		ip.add_constraint(x + y <= 1);
///
///		ip.solve();
///		cout << x << " " << y << endl;
///

#ifndef EASY_IP_HEADER
#define EASY_IP_HEADER

#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>
using std::vector;
using std::size_t;

class OsiSolverInterface;

#ifdef _WIN32
#	ifdef easyip_EXPORTS
#		define EASY_IP_API __declspec(dllexport)
#		define EASY_IP_API_EXTERN_TEMPLATE
#	else
#		define EASY_IP_API __declspec(dllimport)
#		define EASY_IP_API_EXTERN_TEMPLATE extern
#	endif
#else
#	define EASY_IP_API
#	define EASY_IP_API_EXTERN_TEMPLATE
#endif // WIN32

void EASY_IP_API check(bool expr, const char* message);
void EASY_IP_API assertion_failed(const char* expr, const char* file, int line);

#define attest(expr) if (!(expr)) { assertion_failed(#expr, __FILE__, __LINE__); }
#ifndef NDEBUG
	#define dattest(expr) if (!(expr)) { assertion_failed(#expr, __FILE__, __LINE__); }
#else
	#define dattest(expr) ((void)0)
#endif

class IP;

/// Represents a variable in the optimization problem. It is
/// only created by the IP class.
///
class EASY_IP_API Variable
{
friend class IP;
friend class Sum;
friend class LogicalExpression;
friend class BooleanVariable;
public:
	double value() const;
private:
	Variable(size_t index_, const IP* creator_)
		: index(index_), creator(creator_)
	{ }

	size_t index;
	const IP* creator;

	// Everything below ONLY exists so that the variable
	// may be stored in an std::map.
	Variable() 
		: index(-1), creator(nullptr)
	{ }
	template <typename A, typename B, typename C, typename D>
	friend class std::map;
	template <typename A, typename B>
	friend struct std::pair;
};

/// Represents a boolean variable in order to enforce
/// compile-time restrictions to boolean expressions.
class EASY_IP_API BooleanVariable
	: public Variable
{
friend class IP;
public:
	bool value() const;
private:
	BooleanVariable(const Variable& variable)
		: Variable(variable)
	{ }
};

/// Prints the value of the variable after a solution
/// has been obtained.
EASY_IP_API std::ostream& operator << (std::ostream& out, const Variable& variable);
EASY_IP_API std::ostream& operator << (std::ostream& out, const BooleanVariable& variable);

/// A weighted sum of several variables (plus a constant).
///
/// It is created by the user in several ways:
///
///		Sum  s1 = 0;
///		Sum  s2 = x;
///		auto s3 = x + y;
///		auto s4 = 2.0*x - y + 3;
///
class EASY_IP_API Sum
{
friend class Constraint;
friend class IP;
friend class LogicalExpression;
public:

	Sum();
	Sum(const Sum& sum);
	Sum& operator = (const Sum& sum);
	Sum(Sum&& sum);
	Sum& operator = (Sum&& sum);
	Sum(double constant_);
	Sum(const Variable& variable);
	~Sum();

	void add_term(double coeff, const Variable& variable);

	Sum& operator += (const Sum& rhs);
	Sum& operator -= (const Sum& rhs);
	Sum& operator *= (double coeff);
	Sum& operator /= (double coeff);
	void negate();

	double value() const;

protected:
	class Implementation;
	Implementation* impl;

	void match_solvers(const Sum& sum);
};

EASY_IP_API Sum operator * (double coeff, const Variable& variable);
EASY_IP_API Sum operator * (double coeff, Sum sum);
EASY_IP_API Sum operator * (Sum sum, double coeff);

EASY_IP_API Sum operator / (Sum sum, double coeff);

EASY_IP_API Sum operator + (const Sum& lhs, const Sum& rhs);
EASY_IP_API Sum operator + (const Sum& lhs, Sum&& rhs);
EASY_IP_API Sum operator + (Sum&& lhs, const Sum& rhs);
EASY_IP_API Sum operator + (Sum&& lhs, Sum&& rhs);

EASY_IP_API Sum operator - (const Variable& variable);
EASY_IP_API Sum operator - (Sum rhs);

EASY_IP_API Sum operator - (const Sum& lhs, const Sum& rhs);
EASY_IP_API Sum operator - (const Sum& lhs, Sum&& rhs);
EASY_IP_API Sum operator - (Sum&& lhs, const Sum& rhs);
EASY_IP_API Sum operator - (Sum&& lhs, Sum&& rhs);

namespace
{
	std::ostream& operator << (std::ostream& out, const Sum& sum)
	{
		out << sum.value();
		return out;
	}
}

/// A weighted sum of several variables (plus a constant).
///
/// It is created by the user in two ways:
///
///		auto  l1 = !x;
///		auto  l2 = x || y;
///		auto  l3 = x || !y || !z || w;
///
class EASY_IP_API LogicalExpression
{
friend class Constraint;
friend class IP;
public:
	LogicalExpression();
	LogicalExpression(const LogicalExpression& expr);
	LogicalExpression& operator = (const LogicalExpression& expr);
	LogicalExpression(LogicalExpression&& expr);
	LogicalExpression& operator = (LogicalExpression&& expr);
	LogicalExpression(const BooleanVariable& variable, bool is_negated = false);

	LogicalExpression& operator |= (const LogicalExpression& lhs);

private:
	Sum sum;
};

EASY_IP_API LogicalExpression operator || (const LogicalExpression& lhs, const LogicalExpression& rhs);
EASY_IP_API LogicalExpression operator || (LogicalExpression&& lhs, const LogicalExpression& rhs);
EASY_IP_API LogicalExpression operator || (const LogicalExpression& lhs, LogicalExpression&& rhs);
EASY_IP_API LogicalExpression operator || (LogicalExpression&& lhs, LogicalExpression&& rhs);

EASY_IP_API LogicalExpression operator ! (const BooleanVariable& variable);

EASY_IP_API LogicalExpression implication(const BooleanVariable& antecedent, const LogicalExpression& consequent);
EASY_IP_API LogicalExpression implication(const BooleanVariable& antecedent, LogicalExpression&& consequent);

/// Represents a linear constraint.
///
/// It is created from a sum and one of the
/// operators <=, >= or ==.
///
///		auto c1 = x >= 0;
///		auto c2 = 3*x + 4*y == 4; 
///
class EASY_IP_API Constraint
{
friend class IP;
friend EASY_IP_API Constraint operator <= (Sum lhs, const Sum& rhs);
friend EASY_IP_API Constraint operator >= (Sum lhs, const Sum& rhs);
friend EASY_IP_API Constraint operator == (Sum lhs, const Sum& rhs);
public:
	Constraint();  // For std::vector.
	Constraint(const LogicalExpression& expression);
	Constraint(LogicalExpression&& expression);
private:
	Constraint(double lower_bound_, const Sum& sum_, double upper_bound_);
	Constraint(double lower_bound_, Sum&& sum_, double upper_bound_);
	double lower_bound, upper_bound;
	Sum sum;
};

// TODO: More overloads here will reduce temporaries.
EASY_IP_API Constraint operator <= (Sum lhs, const Sum& rhs);
EASY_IP_API Constraint operator >= (Sum lhs, const Sum& rhs);
EASY_IP_API Constraint operator == (Sum lhs, const Sum& rhs);

class EASY_IP_API ConstraintList
{
friend class IP;
public:
	ConstraintList(Constraint&& constraint);
	ConstraintList(ConstraintList&&);
	~ConstraintList();

	ConstraintList& operator &= (Constraint&&);
private:
	ConstraintList();
	ConstraintList(const ConstraintList&);
	ConstraintList& operator = (const ConstraintList&);

	class Implementation;
	Implementation* impl;
};

EASY_IP_API ConstraintList operator && (Constraint&&, Constraint&&);
EASY_IP_API ConstraintList operator && (ConstraintList&&, Constraint&&);

/// Represents an integer (linear) program.
///
///
class EASY_IP_API IP
{
public:
	enum VariableType {Boolean, Binary = Boolean, Integer, Real};

	IP();
	IP(const IP&) = delete;
	IP(IP&&);
	~IP();

	/// Adds a variable to the optimization problems. Variables must not
	/// have their values queried after the creating IP class has been
	/// destroyed.
	Variable add_variable(VariableType type = Boolean, double this_cost = 0.0);
	/// Adds a boolean variable to the optimization problems.
	BooleanVariable add_boolean(double this_cost = 0.0);

	/// Creates a vector of variables.
	vector<Variable> add_vector(int n, VariableType type = Boolean, double this_cost = 0.0);
	/// Creates a vector of logical variables.
	vector<BooleanVariable> add_boolean_vector(int n, double this_cost = 0.0);

	/// Creates a grid of variables.
	vector<vector<Variable>> add_grid(int m, int n, VariableType type = Boolean, double this_cost = 0.0);
	/// Creates a grid of logical variables.
	vector<vector<BooleanVariable>> add_boolean_grid(int m, int n, double this_cost = 0.0);

	/// Creates a 3D grid of variables.
	vector<vector<vector<Variable>>> add_cube(int m, int n, int o, VariableType type = Boolean, double this_cost = 0.0);
	/// Creates a 3D grid of logical variables.
	vector<vector<vector<BooleanVariable>>> add_boolean_cube(int m, int n, int o, double this_cost = 0.0);

	/// Adds the constraint
	///    L <= constraint <= U
	void add_constraint(double L, const Sum& sum, double U);

	void add_constraint(const Constraint& constraint);
	void add_constraint(const ConstraintList& list);

	void add_constraint(const BooleanVariable& variable);

	//void add_constraint(const LogicalExpression& expression);

	/// Adds the constraint
	///    L <= variable <= U
	void set_bounds(double L, const Variable& variable, double U);

	/// Adds a linear function (plus a constant) to the
	/// objective function.
	void add_objective(const Sum& sum);

	enum Solver {Default, CPLEX, MOSEK};
	/// Switches to an external solver (if available).
	void set_external_solver(Solver solver);

	typedef std::function<void()> CallBack;
	/// Solves the integer program. Returns false if the program
	/// is infeasible or unbounded.
	///
	/// If a callback is used, preprocessing will be turned off.
	bool solve(const CallBack& callback_function = nullptr);

	bool next_solution();

	/// Retrieves the value of a variable in the solution.
	double get_solution(const Variable& variable) const;
	/// Retrieves the value of a variable in the solution.
	bool get_solution(const BooleanVariable& variable) const;

	/// Evaluates a sum of variables in the solution.
	double get_solution(const Sum& sum) const;
	/// Evaluates an expression with the solution.
	bool get_solution(const LogicalExpression& expression) const;

	size_t get_number_of_variables() const;

	// Resets the optimization problem and starts over.
	void clear();


	//
	// More advanced functionality
	//

	/// Creates and returns a pointer to a solver. 
	///
	/// If existing_solver is not empty, it returns it after loading the problem. 
	/// This allows any interface to be created, e.g. OsiGrbSolverInterface,
	/// which this library might not know about.
	std::unique_ptr<OsiSolverInterface> get_problem(std::unique_ptr<OsiSolverInterface> existing_solver = nullptr);


protected:
	vector<double>& get_rhs_lower();
	const vector<double>& get_rhs_lower() const;
	vector<double>& get_rhs_upper();
	const vector<double>& get_rhs_upper() const;
	vector<int>& get_rows();
	const vector<int>& get_rows() const;
	vector<int>& get_cols();
	const vector<int>& get_cols() const;
	vector<double>& get_values();
	const vector<double>& get_values() const;

	const vector<double>& get_var_lb() const;
	const vector<double>& get_var_ub() const;
	const vector<double>& get_cost() const;
	const vector<std::size_t>& get_integer_variables() const;

	size_t get_variable_index(const Variable& x) const { return x.index; }

	vector<double>& get_solution();

private:
	class Implementation;
	Implementation* impl;
};

#endif
