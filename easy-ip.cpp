// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

#include <easy-ip.h>
#include <easy-ip-internal.h>


void check(bool expr, const char* message)
{
	if (!expr) {
		throw std::invalid_argument(message);
	}
}

void assertion_failed(const char* expr, const char* file_cstr, int line)
{
	using namespace std;

	// Extract the file name only.
	string file(file_cstr);
	auto pos = file.find_last_of("/\\");
	if (pos == string::npos) {
		pos = 0;
	}
	file = file.substr(pos + 1);  // Returns empty string if pos + 1 == length.

	stringstream sout;
	sout << "Assertion failed: " << expr << " in " << file << ":" << line << ".";
	throw runtime_error(sout.str().c_str());
}

double EASY_IP_API wall_time()
{
	#ifdef HAS_OPENMP
		return ::omp_get_wtime();
	#else
		return 0;
	#endif
}

double Variable::value() const
{
	return creator->get_solution(*this);
}

bool BooleanVariable::bool_value() const
{
	return creator->get_solution(*this);
}

std::ostream& operator << (std::ostream& out, const Variable& variable)
{
	out << variable.value();
	return out;
}

std::ostream& operator << (std::ostream& out, const BooleanVariable& variable)
{
	out << variable.bool_value();
	return out;
}

Sum::Sum()
	: impl(new Implementation)
{ }

Sum::Sum(const Sum& sum)
	: impl(new Implementation(*sum.impl))
{ }

Sum& Sum::operator = (const Sum& sum)
{
	*impl = *sum.impl;
	return *this;
}

Sum::Sum(Sum&& sum)
	: impl(nullptr)
{
	impl = sum.impl;
	sum.impl = nullptr;
}

Sum& Sum::operator = (Sum&& sum)
{
	if (impl) {
		delete impl;
	}

	impl = sum.impl;
	sum.impl = nullptr;
	return *this;
}

Sum::Sum(double constant_)
	: impl(new Implementation(constant_))
{ }

Sum::Sum(const Variable& variable)
	: impl(new Implementation(variable))
{ }

Sum::~Sum()
{
	if (impl) {
		delete impl;
	}
}

void Sum::add_term(double coeff, const Variable& variable)
{
	impl->cols.push_back(static_cast<int>(variable.index));
	impl->values.push_back(coeff);
}

double Sum::value() const
{
	if (!impl->creator) {
		// This happens if Sum is constant. No variables
		// have been added.
		attest(impl->cols.empty());
		return impl->constant;
	}

	return impl->creator->get_solution(*this);
}

Sum& Sum::operator += (const Sum& rhs)
{
	match_solvers(rhs);

	impl->constant += rhs.impl->constant;
	for (size_t i = 0; i < rhs.impl->cols.size(); ++i) {
		impl->cols.push_back(rhs.impl->cols[i]);
		impl->values.push_back(rhs.impl->values[i]);
	}
	return *this;
}

Sum& Sum::operator -= (const Sum& rhs)
{
	match_solvers(rhs);

	impl->constant -= rhs.impl->constant;
	for (size_t i = 0; i < rhs.impl->cols.size(); ++i) {
		impl->cols.push_back(rhs.impl->cols[i]);
		impl->values.push_back(-rhs.impl->values[i]);
	}
	return *this;
}

Sum& Sum::operator *= (double coeff)
{
	if (coeff == 0.0) {
		impl->values.clear();
		impl->cols.clear();
		impl->constant = 0.0;
		return *this;
	}

	for (auto& value : impl->values) {
		value *= coeff;
	}
	impl->constant *= coeff;
	return *this;
}

Sum& Sum::operator /= (double coeff)
{
	if (coeff == 0.0) {
		throw std::runtime_error("Sum: Division by zero.");
	}

	*this *= (1.0 / coeff);
	return *this;
}

Sum operator * (double coeff, const Variable& variable)
{
	Sum sum(variable);
	sum *= coeff;
	return sum;
}

Sum operator *  (double coeff, Sum sum)
{
	sum *= coeff;
	return sum;
}

Sum operator *  (Sum sum, double coeff)
{
	sum *= coeff;
	return sum;
}

Sum operator /  (Sum sum, double coeff)
{
	sum /= coeff;
	return sum;
}

Sum operator + (const Sum& lhs, const Sum& rhs)
{
	// Copy needed.
	Sum sum(rhs);
	sum += lhs;
	return sum;
}

Sum operator + (const Sum& lhs, Sum&& rhs)
{
	Sum sum(std::move(rhs));
	sum += lhs;
	return sum;
}

Sum operator + (Sum&& lhs, const Sum& rhs)
{
	Sum sum(std::move(lhs));
	sum += rhs;
	return sum;
}

Sum operator + (Sum&& lhs, Sum&& rhs)
{
	Sum sum(std::move(lhs));
	sum += rhs;
	return sum;
}

Sum operator - (Sum sum)
{
	sum.negate();
	return sum;
}

Sum operator - (const Variable& variable)
{
	Sum sum(variable);
	return -sum;
}

void Sum::negate()
{
	impl->constant = -impl->constant;
	for (auto& value : impl->values) {
		value = -value;
	}
}

Sum operator - (const Sum& lhs, const Sum& rhs)
{
	// Copy needed.
	Sum sum(lhs);
	sum -= rhs;
	return sum;
}

Sum operator - (const Sum& lhs, Sum&& rhs)
{
	Sum sum(std::move(rhs));
	sum.negate();
	sum += lhs;
	return sum;
}

Sum operator - (Sum&& lhs, const Sum& rhs)
{
	Sum sum(std::move(lhs));
	sum -= rhs;
	return sum;
}

Sum operator - (Sum&& lhs, Sum&& rhs)
{
	Sum sum(std::move(rhs));
	sum.negate();
	sum += lhs;
	return sum;
}

void Sum::match_solvers(const Sum& sum)
{
	check(impl->creator == nullptr || sum.impl->creator == nullptr || impl->creator == sum.impl->creator,
	      "Variables from different solver can not be mixed.");
	if (impl->creator == nullptr) {
		impl->creator = sum.impl->creator;
	}
}

LogicalExpression operator ! (const BooleanVariable& variable)
{
	// TODO: check that variable is integer through creator.
	return LogicalExpression(variable, true);
}

LogicalExpression::LogicalExpression(const BooleanVariable& variable, bool is_negated)
{
	if (is_negated) {
		sum += 1;
		sum -= variable;
	}
	else {
		sum += variable;
	}
}

LogicalExpression::LogicalExpression()
{ }

LogicalExpression::LogicalExpression(const LogicalExpression& expr)
{
	*this = expr;
}

LogicalExpression& LogicalExpression::operator = (const LogicalExpression& expr)
{
	sum = expr.sum;
	return *this;
}


LogicalExpression::LogicalExpression(LogicalExpression&& expr)
{
	*this = std::move(expr);
}

LogicalExpression& LogicalExpression::operator = (LogicalExpression&& expr)
{
	sum = std::move(expr.sum);
	return *this;
}

LogicalExpression& LogicalExpression::operator |= (const LogicalExpression& lhs)
{
	sum.match_solvers(lhs.sum);

	sum += lhs.sum;
	return *this;
}

LogicalExpression operator || (const LogicalExpression& lhs, const LogicalExpression& rhs)
{
	// Copying necessary.
	LogicalExpression result(lhs);
	result |= rhs;
	return result;
}

EASY_IP_API LogicalExpression operator || (LogicalExpression&& lhs, const LogicalExpression& rhs)
{
	LogicalExpression result(std::move(lhs));
	result |= rhs;
	return result;
}

EASY_IP_API LogicalExpression operator || (const LogicalExpression& lhs, LogicalExpression&& rhs)
{
	LogicalExpression result(std::move(rhs));
	result |= lhs;
	return result;
}

EASY_IP_API LogicalExpression operator || (LogicalExpression&& lhs, LogicalExpression&& rhs)
{
	LogicalExpression result(std::move(lhs));
	result |= rhs;
	return result;
}

LogicalExpression implication(const BooleanVariable& antecedent, const LogicalExpression& consequent)
{
	return !antecedent || consequent;
}

LogicalExpression implication(const BooleanVariable& antecedent, LogicalExpression&& consequent)
{
	return !antecedent || std::move(consequent);
}

Constraint::Constraint(const LogicalExpression& expression)
	: lower_bound(1), upper_bound(static_cast<double>(expression.sum.impl->values.size())), sum(expression.sum)
{ }

Constraint::Constraint(LogicalExpression&& expression)
	: lower_bound(1), upper_bound(static_cast<double>(expression.sum.impl->values.size())), sum(std::move(expression.sum))
{ }

Constraint::Constraint(double lower_bound_, const Sum& sum_, double upper_bound_)
		: lower_bound(lower_bound_), upper_bound(upper_bound_), sum(sum_)
{ }

Constraint::Constraint(double lower_bound_, Sum&& sum_, double upper_bound_)
		: lower_bound(lower_bound_), upper_bound(upper_bound_), sum(std::move(sum_))
{ }

Constraint operator <= (Sum lhs, const Sum& rhs)
{
	lhs -= rhs;
	return Constraint(-1e100, std::move(lhs), 0.0);
}

Constraint operator >= (Sum lhs, const Sum& rhs)
{
	lhs -= rhs;
	return Constraint(0.0, std::move(lhs), 1e100);
}

Constraint operator == (Sum lhs, const Sum& rhs)
{
	lhs -= rhs;
	return Constraint(0.0, std::move(lhs), 0.0);
}

class ConstraintList::Implementation
{
public:
	Implementation() { }

	template<typename T>
	Implementation(T&& t) : constraints(std::forward<T>(t)) { }

	vector<Constraint> constraints;
};

ConstraintList::ConstraintList()
	: impl(new Implementation)
{ }

ConstraintList::ConstraintList(ConstraintList&& rhs)
	: impl(rhs.impl)
{
	rhs.impl = nullptr;
}

ConstraintList::ConstraintList(Constraint&& constraint)
	: impl(new Implementation)
{
	impl->constraints.push_back(std::move(constraint));
}

ConstraintList& ConstraintList::operator &= (Constraint&& rhs)
{
	impl->constraints.push_back(std::move(rhs));
	return *this;
}

ConstraintList::~ConstraintList()
{
	if (impl) {
		delete impl;
	}
}

ConstraintList operator && (Constraint&& lhs, Constraint&& rhs)
{
	ConstraintList list(std::move(lhs));
	list &= std::move(rhs);
	return list;
}

ConstraintList operator && (ConstraintList&& lhs, Constraint&& rhs)
{
	ConstraintList list(std::move(lhs));
	list &= std::move(rhs);
	return list;
}

IP::IP() 
	: impl(new Implementation(this))
{
};

IP::IP(IP&& lhs) 
	: impl(lhs.impl)
{
	impl->creator = this;
	lhs.impl = nullptr;
};

IP::~IP()
{
	if (impl) {
		delete impl;
		impl = nullptr;
	}
}

Variable IP::add_variable(VariableType type, double this_cost)
{
	if (type == Boolean) {
		impl->var_lb.push_back(0.0);
		impl->var_ub.push_back(1.0);
	}
	else {
		impl->var_lb.push_back(-1e100);
		impl->var_ub.push_back(1e100);
	}

	impl->cost.push_back(this_cost);
	auto index = impl->cost.size() - 1;

	if (type == Boolean || type == Integer) {
		impl->integer_variables.push_back(index);
	}

	attest(impl->var_lb.size() == impl->cost.size());
	attest(impl->var_ub.size() == impl->cost.size());
	return Variable(index, this);
}

BooleanVariable IP::add_boolean(double this_cost)
{
	auto variable = add_variable(Boolean, this_cost);
	return BooleanVariable(variable);
}

Sum IP::add_variable_as_booleans(int lower_bound, int upper_bound)
{
	Sum sum = 0;
	Sum constraint = 0;
	for (int i = lower_bound; i <= upper_bound; ++i) {
		auto var = add_boolean();
		sum += i * var;
		constraint += var;
	}
	add_constraint(constraint == 1);
	return sum;
}

Sum IP::add_variable_as_booleans(const std::initializer_list<int>& values)
{
	Sum sum = 0;
	Sum constraint = 0;
	for (int i: values) {
		auto var = add_boolean();
		sum += i * var;
		constraint += var;
	}
	add_constraint(constraint == 1);
	return sum;
}

vector<Variable> IP::add_vector(int n, VariableType type, double this_cost)
{
	vector<Variable> v;
	for (int i = 0; i < n; ++i) {
		v.push_back(add_variable(type, this_cost));
	}
	return v;
}

vector<BooleanVariable> IP::add_boolean_vector(int n, double this_cost)
{
	vector<BooleanVariable> v;
	for (int i = 0; i < n; ++i) {
		v.push_back(add_boolean(this_cost));
	}
	return v;
}

template<typename Var>
vector<vector<Var>> grid_creator(int m, int n, std::function<Var()> var_creator)
{
	attest(m >= 0 && n >= 0);
	vector<vector<Var>> grid;

	for (int i = 0; i < m; ++i) {
		grid.push_back(vector<Var>());
		for (int j = 0; j < n; ++j) {
			grid.back().push_back(var_creator());
		}
	}

	return grid;
}

vector<vector<Variable>> IP::add_grid(int m, int n, VariableType type, double this_cost)
{
	auto var_creator = [&]() { return add_variable(type, this_cost); };
	return grid_creator<Variable>(m, n, var_creator);
}

vector<vector<BooleanVariable>> IP::add_boolean_grid(int m, int n, double this_cost)
{
	auto var_creator = [&]() { return add_boolean(this_cost); };
	return grid_creator<BooleanVariable>(m, n, var_creator);
}

template<typename Var>
vector<vector<vector<Var>>> cube_creator(int m, int n, int o, std::function<Var()> var_creator)
{
	attest(m >= 0 && n >= 0 && o >= 0);
	vector<vector<vector<Var>>> grid;

	for (int i = 0; i < m; ++i) {
		grid.push_back(vector<vector<Var>>());
		for (int j = 0; j < n; ++j) {
			grid.back().push_back(vector<Var>());
			for (int k = 0; k < o; ++k) {
				grid.back().back().push_back(var_creator());
			}
		}
	}

	return grid;
}

vector<vector<vector<Variable>>> IP::add_cube(int m, int n, int o, VariableType type, double this_cost)
{
	auto var_creator = [&]() { return add_variable(type, this_cost); };
	return cube_creator<Variable>(m, n, o, var_creator);
}

vector<vector<vector<BooleanVariable>>> IP::add_boolean_cube(int m, int n, int o, double this_cost)
{
	auto var_creator = [&]() { return add_boolean(this_cost); };
	return cube_creator<BooleanVariable>(m, n, o, var_creator);
}

// Adds the constraint
//    L <= constraint <= U
int IP::add_constraint(double L, const Sum& sum, double U)
{
	impl->check_creator(sum);

	//std::cerr << L-sum.constant << " <= ";
	//for (int i = 0; i < sum.cols.size(); ++i) {
	//	std::cerr << sum.values[i] << "*x" << sum.cols[i] << " ";
	//}
	//std::cerr << " <= " << U-sum.constant << std::endl;

	if (sum.impl->cols.empty()) {
		check(L <= sum.impl->constant && sum.impl->constant <= U,
			"A constraint that is always false may not be added.");
		return 0;
	}

	impl->rhs_lower.push_back(L - sum.impl->constant);
	impl->rhs_upper.push_back(U - sum.impl->constant);
	auto row_index = impl->rhs_upper.size() - 1;
	for (size_t i = 0; i < sum.impl->cols.size(); ++i) {
		impl->rows.push_back(static_cast<int>(row_index));
		impl->cols.push_back(sum.impl->cols[i]);
		impl->values.push_back(sum.impl->values[i]);
	}

	attest(impl->rows.size() == impl->values.size());
	attest(impl->cols.size() == impl->values.size());
	attest(impl->rhs_lower.size() == impl->rhs_upper.size());
	return 1;
}

int IP::add_constraint(const Constraint& constraint)
{
	return add_constraint(constraint.lower_bound, constraint.sum, constraint.upper_bound);
}

int IP::add_constraint(const ConstraintList& list)
{
	int added = 0;
	for (auto& constraint: list.impl->constraints) {
		added += add_constraint(constraint);
	}
	return added;
}

void IP::add_constraint(const BooleanVariable& variable)
{
	set_bounds(1, variable, 1);
}

// Adds the constraint
//    L <= variable <= U
void IP::set_bounds(double L, const Variable& variable, double U)
{
	impl->check_creator(variable);

	impl->var_lb[variable.index] = L;
	impl->var_ub[variable.index] = U;
}

void IP::add_objective(const Sum& sum)
{
	impl->check_creator(sum);

	for (size_t i = 0; i < sum.impl->cols.size(); ++i) {
		impl->cost.at(sum.impl->cols[i]) += sum.impl->values[i];
	}
	// TODO: handle constant.
}

void IP::Implementation::check_creator(const Variable& t) const
{
	check(t.creator == nullptr || t.creator == creator,
	      "Variable comes from a different solver.");
}

void IP::Implementation::check_creator(const Sum& t) const
{
	check(t.impl->creator == nullptr || t.impl->creator == creator,
		"Sum comes from a different solver.");
}

double IP::get_solution(const Variable& variable) const
{
	impl->check_creator(variable);

	return impl->solution.at(variable.index);
}

bool IP::get_solution(const BooleanVariable& variable) const
{
	impl->check_creator(variable);

	return impl->solution.at(variable.index) > 0.5;
}

double IP::get_solution(const Sum& sum) const
{
	impl->check_creator(sum);

	double value = sum.impl->constant;
	for (size_t i = 0; i < sum.impl->cols.size(); ++i) {
		value += sum.impl->values[i] * impl->solution.at(sum.impl->cols[i]);
	}
	return value;
}

bool IP::get_solution(const LogicalExpression& expression) const
{
	return get_solution(expression.sum) > 0.5;
}

void IP::set_external_solver(Solver solver)
{
	impl->external_solver = solver;

	impl->minisat_solver.release();
	impl->literals.clear();
	impl->objective_function_literals.clear();
	impl->objective_function_slack_literals.clear();

	if (solver == CPLEX) {
		#ifndef HAS_CPLEX
			throw std::runtime_error("IP::set_external_solver: CPLEX not installed.");
		#endif
	}
	else if (solver == MOSEK) {
		#ifndef HAS_MOSEK
			throw std::runtime_error("IP::set_external_solver: MOSEK not installed.");
		#endif
	}
}

void IP::set_time_limit(double seconds)
{
	impl->time_limit_in_seconds = seconds;
}

void IP::allow_ignoring_cost_function()
{
	impl->allow_ignoring_cost_function = true;
}

size_t IP::get_number_of_variables() const
{
	return impl->cost.size();
}

void IP::clear()
{
	impl->rhs_lower.clear();
	impl->rhs_upper.clear();
	impl->rows.clear();
	impl->cols.clear();
	impl->values.clear();

	impl->var_lb.clear();
	impl->var_ub.clear();
	impl->cost.clear();

	impl->solution.clear();

	impl->integer_variables.clear();
}

int IP::add_max_consequtive_constraints(int N, const std::vector<Sum>& variables)
{
	if (N >= variables.size()) {
		return 0;
	}

	int constraints_added = 0;

	for (int d = 0; d < variables.size() - N; ++d) {
		Sum active_in_window = 0;
		for (int d2 = d; d2 < d + N + 1; ++d2) {
			active_in_window += variables.at(d2);
		}
		constraints_added += add_constraint(active_in_window <= N);
	}
	return constraints_added;
}

int IP::add_min_consequtive_constraints(int N, const std::vector<Sum>& variables, bool OK_at_the_border)
{
	if (N <= 1) {
		return 0;
	}
	attest(N <= variables.size());
	if (N == variables.size()) {
		for (auto& var: variables) {
			add_constraint(var == 1);
		}
		return 0;
	}

	int constraints_added = 0;

	for (int window_size = 1; window_size <= N - 1; ++window_size) {
		// Look for windows of size minimum - 1 along with the
		// surrounding slots.
		//  
		// […] [x1] [y1] [x2] […]
		//
		// x1 = 0 ∧ x2 = 0 ⇒ y1 = 0
		// ⇔
		// x1 + x2 - y1 ≥ 0
		//
		// Then add windows with more y variables. E.g. 
		//
		// x1 + x2 - y1 - y2 - y3 ≥ -2.

		for (int window_start = 0; window_start < variables.size() - window_size + 1; ++window_start) {

			Sum constraint = 0;

			if (window_start - 1 >= 0) {
				constraint += variables.at(window_start - 1);
			}
			else if (OK_at_the_border) {
				continue;
			}

			for (int i = window_start; i < window_start + window_size; ++i) {
				constraint -= variables.at(i);
			}

			if (window_start + window_size < variables.size()) {
				constraint += variables.at(window_start + window_size);
			}
			else if (OK_at_the_border) {
				continue;
			}

			constraints_added += add_constraint(constraint >= -window_size + 1);
		}
	}
	return constraints_added;
}

void IP::save_MPS(const std::string& file_name)
{
	OsiSolverInterface* solver = nullptr;
	std::unique_ptr<OsiSolverInterface> new_problem;

	if (impl->use_osi() && impl->problem) {
		solver = impl->problem.get();
	}
	else if (impl->model) {
		attest(impl->model);
		solver = impl->model->solver();
	}
	else {
		get_problem(new_problem);
		solver = new_problem.get();
	}
	solver->writeMps(file_name.c_str());
}

void IP::save_CNF(const std::string& file_name)
{
	impl->convert_to_minisat();
	auto f = std::fopen(file_name.c_str(), "w");
	impl->minisat_solver->toDimacs(f, {});
	std::fclose(f);
}

vector<double>& IP::get_rhs_lower() { return impl->rhs_lower; }
vector<double>& IP::get_rhs_upper() { return impl->rhs_upper; }
vector<int>& IP::get_rows() { return impl->rows; }
vector<int>& IP::get_cols() { return impl->cols; }
vector<double>& IP::get_values() { return impl->values; }

const vector<double>& IP::get_var_lb() const { return impl->var_lb; }
const vector<double>& IP::get_var_ub() const { return impl->var_ub; }
const vector<double>& IP::get_cost() const { return impl->cost; }
const vector<std::size_t>& IP::get_integer_variables() const { return impl->integer_variables; }

vector<double>& IP::get_solution() { return impl->solution; }
