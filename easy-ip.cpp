// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
//#include <typeid>
#include <vector>

// Using Clp as the solver
#include <coin/OsiClpSolverInterface.hpp>

#ifdef HAS_CPLEX
	#include <coin/OsiCpxSolverInterface.hpp>
#endif

#ifdef HAS_MOSEK
	#include <coin/OsiMskSolverInterface.hpp>
#endif


#include "minisat/core/Solver.h"


#include <coin/CbcModel.hpp>
#include <coin/CbcEventHandler.hpp>
#include <coin/CbcCutGenerator.hpp>
#include <coin/CbcStrategy.hpp>
#include <coin/CbcHeuristicLocal.hpp>

#include <coin/CglGomory.hpp>
#include <coin/CglProbing.hpp>
#include <coin/CglKnapsackCover.hpp>
#include <coin/CglRedSplit.hpp>
#include <coin/CglClique.hpp>
#include <coin/CglOddHole.hpp>
#include <coin/CglFlowCover.hpp>
#include <coin/CglMixedIntegerRounding2.hpp>
#include <coin/CglPreProcess.hpp>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

#include <easy-ip.h>

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

bool BooleanVariable::value() const
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
	out << variable.value();
	return out;
}

class Sum::Implementation
{
public:
	Implementation()
		: constant(0.0),
		  creator(nullptr)
	{ }

	Implementation(const Implementation& rhs)
	{
		*this = rhs;
	}

	Implementation& operator = (const Implementation& rhs)
	{
		cols = rhs.cols;
		values = rhs.values;
		constant = rhs.constant;
		creator = rhs.creator;
		//std::cerr << "Copy sum.\n";
		//throw std::runtime_error("Copying took place.");
		return *this;
	}

	Implementation(double constant_)
		: constant(constant_),
		  creator(nullptr)
	{ }

	Implementation(const Variable& variable)
		: constant(0.0),
		  creator(variable.creator)
	{
		cols.push_back(int(variable.index));
		values.push_back(1.0);
	}

	double constant;
	vector<int> cols;
	vector<double> values;

	const IP* creator;
};

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

class IP::Implementation
{
public:
	Implementation(const IP* creator_)
		: external_solver(IP::Default), creator(creator_)
	{ }

	bool parse_solution();
	
	bool use_osi() const;

	bool solve_minisat();
	bool next_minisat();

	vector<double> rhs_lower;
	vector<double> rhs_upper;
	vector<int> rows;
	vector<int> cols;
	vector<double> values;

	vector<double> var_lb;
	vector<double> var_ub;
	vector<double> cost;

	vector<double> solution;

	vector<size_t> integer_variables;

	Solver external_solver;

	double time_limit_in_seconds = -1;

	bool preprocess;
	std::unique_ptr<OsiSolverInterface> problem;
	std::unique_ptr<CbcModel> model;
	
	std::unique_ptr<Minisat::Solver> minisat_solver;
	vector<Minisat::Lit> literals;
	vector<Minisat::Lit> objective_function_literals;
	vector<Minisat::Lit> objective_function_slack_literals;
	int sat_objective_offset = 0;

	std::vector<std::unique_ptr<CglCutGenerator>> generators;

	bool allow_ignoring_cost_function = false;

	void check_creator(const Variable& t) const;
	void check_creator(const Sum& t) const;

	const IP* creator;
};

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
void IP::add_constraint(double L, const Sum& sum, double U)
{
	impl->check_creator(sum);

	//std::cerr << L-sum.constant << " <= ";
	//for (int i = 0; i < sum.cols.size(); ++i) {
	//	std::cerr << sum.values[i] << "*x" << sum.cols[i] << " ";
	//}
	//std::cerr << " <= " << U-sum.constant << std::endl;

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
}

void IP::add_constraint(const Constraint& constraint)
{
	add_constraint(constraint.lower_bound, constraint.sum, constraint.upper_bound);
}

void IP::add_constraint(const ConstraintList& list)
{
	for (auto& constraint: list.impl->constraints) {
		add_constraint(constraint);
	}
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


class MyEventHandler
	: public CbcEventHandler 
{

public:

	virtual CbcAction event(CbcEvent whichEvent)
	{
		// If in sub tree carry on
		if (!model_->parentModel()) {
			if (whichEvent == CbcEventHandler::solution ||
			    whichEvent == CbcEventHandler::heuristicSolution) {

				if (callback_function) {

					int org_n; 
					auto n = model->getNumCols();
					auto best_solution = model->bestSolution();
					const int* org_columns;

					org_n = n;
					org_columns = model->originalColumns();

					solution->clear();
					solution->resize(org_n, 0.0);
					if (org_columns) {
						for (int i = 0; i < n; ++i) {
							solution->at(org_columns[i]) = best_solution[i];
						}
					}
					else {
						std::copy(best_solution, best_solution + n, solution->begin());
					}
				
					callback_function();
				}

				return noAction; // carry on
			} else {
				return noAction; // carry on
			}
		} else {
			return noAction; // carry on
		}
	}

	MyEventHandler(const IP::CallBack& callback_function_,
	               std::vector<double>* solution_,
	               const CbcModel* model_)
		: callback_function(callback_function_),
		  solution(solution_),
		  model(model_)
	{ }

	MyEventHandler(const MyEventHandler& rhs)
		: CbcEventHandler(rhs),
		  callback_function(rhs.callback_function),
		  solution(rhs.solution),
		  model(rhs.model)
	{ }

	virtual ~MyEventHandler()
	{ }

	MyEventHandler& operator=(const MyEventHandler & rhs)
	{
		if (this != &rhs) {
			callback_function = rhs.callback_function;
			solution = rhs.solution;
			model = rhs.model;
		}
		return *this;
	}


	virtual CbcEventHandler* clone() const
	{
		return new MyEventHandler(*this);
	}

protected:
	IP::CallBack callback_function;
	std::vector<double>* solution;
	const CbcModel* model;
};

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

void IP::get_problem(std::unique_ptr<OsiSolverInterface>& problem)
{
	if (!problem) {
		auto clp_problem = new OsiClpSolverInterface;
		problem.reset(clp_problem);
		// Turn off information from the LP solver if we are
		// using branch and cut/bound since this means a lot
		// of LPs.
		if ( ! impl->integer_variables.empty()) {
			clp_problem->setLogLevel(0);
		}
	}

	attest(impl->var_lb.size() == impl->cost.size());
	attest(impl->var_ub.size() == impl->cost.size());

	attest(impl->rows.size() == impl->values.size());
	attest(impl->cols.size() == impl->values.size());
	attest(impl->rhs_lower.size() == impl->rhs_upper.size());

	// Check if last_index is present.
	auto last_index = impl->cost.size() - 1;
	bool last_index_present = false;
	for (auto var : impl->cols) {
		if (var == last_index) {
			last_index_present = true;
			break;
		}
	}
	// If not, add a dummy constraint to satisfy Clp.
	if (!last_index_present) {
		add_constraint(-1e100, 100*Variable(last_index, this), 1e100);
	}

	CoinPackedMatrix coinMatrix(false,
	                            &impl->rows[0],
	                            &impl->cols[0],
	                            &impl->values[0],
	                            CoinBigIndex(impl->values.size()) );

	problem->loadProblem(coinMatrix,
	                     &impl->var_lb[0],
	                     &impl->var_ub[0],
	                     &impl->cost[0],
	                     &impl->rhs_lower[0],
	                     &impl->rhs_upper[0]);

	for (auto index: impl->integer_variables) {
		problem->setInteger(static_cast<int>(index));
	}

	//problem->writeMps("problem");
}

bool IP::Implementation::use_osi() const
{
	return external_solver != IP::Default || integer_variables.empty();
}

bool IP::Implementation::parse_solution()
{
	OsiSolverInterface* solved_problem = nullptr;

	if (use_osi()) {
		solved_problem = problem.get();

		if (solved_problem->isAbandoned()) {
			std::cerr << "-- Abandoned." << std::endl;
			return false;
		}
		else if (solved_problem->isProvenPrimalInfeasible()) {
			std::cerr << "-- Infeasible." << std::endl;
			return false;
		}
		else if (solved_problem->isProvenDualInfeasible()) {
			std::cerr << "-- Unbounded." << std::endl;
			return false;
		}
		else if (solved_problem->isPrimalObjectiveLimitReached()) {
			std::cerr << "-- Primal objective limit." << std::endl;
			return false;
		}
		else if (solved_problem->isDualObjectiveLimitReached()) {
			std::cerr << "-- Dual objective limit." << std::endl;
			return false;
		}
		else if (solved_problem->isIterationLimitReached()) {
			std::cerr << "-- Iteration limit." << std::endl;
			return false;
		}
		else if (!integer_variables.empty() && !solved_problem->isProvenOptimal()) {
			std::cerr << "-- Not optimal." << std::endl;
			return false;
		}
	}
	else {
		if (model->isProvenInfeasible()) {
			//throw std::runtime_error("Problem infeasible.");
			return false;
		}
	
		if (model->isProvenDualInfeasible()) {
			//throw std::runtime_error("Problem unbounded.");
			return false;
		}

		if (!model->isProvenOptimal()) {
			//throw std::runtime_error("Time limit reached.");
			return false;
		}

		solved_problem = model->solver();
	}

	int numberColumns = solved_problem->getNumCols();
	const double* raw_solution = solved_problem->getColSolution();

	solution.clear();
	for (size_t i = 0; i < numberColumns; ++i) {
		solution.push_back(raw_solution[i]);
	}

	return true;
}

void IP::set_time_limit(double seconds)
{
	impl->time_limit_in_seconds = seconds;
}

bool IP::solve(const CallBack& callback_function, bool silent_mode)
{
	if (impl->external_solver == IP::Minisat) {
		return impl->solve_minisat();
	}
	if(impl->external_solver == IP::CPLEX) {
		#ifdef HAS_CPLEX
			auto cplex_solver = std::unique_ptr<OsiSolverInterface>(new OsiCpxSolverInterface);
			get_problem(cplex_solver);
			impl->problem = std::move(cplex_solver);
		#endif
	}
	else if(impl->external_solver == IP::MOSEK) {
		#ifdef HAS_MOSEK
			auto mosek_solver = std::unique_ptr<OsiSolverInterface>(new OsiMskSolverInterface);
			get_problem(mosek_solver);
			impl->problem = std::move(mosek_solver)
		#endif
	}
	else {
		get_problem(impl->problem);
	}

	if (impl->external_solver != IP::Default || impl->integer_variables.empty()) {
		
		if (impl->integer_variables.empty()) {
			impl->problem->initialSolve();
		}
		else {
			impl->problem->branchAndBound();
		}
	}
	else {
		// Pass the solver with the problem to be solved to CbcModel 
		impl->model.reset(new CbcModel(*impl->problem.get()));

		// Only the most important log messages.
		impl->model->setLogLevel(1);
		if (silent_mode) {
			impl->model->setLogLevel(0);
		}

		if (impl->time_limit_in_seconds > 0) {
			impl->model->setDblParam(CbcModel::CbcMaximumSeconds, impl->time_limit_in_seconds);
		}

		if (callback_function) {
			MyEventHandler my_event_handler(callback_function,
			                                &impl->solution,
			                                impl->model.get());
			impl->model->passInEventHandler(&my_event_handler);
		}

		impl->generators.clear();

		// Add in generators
		// Experiment with -1 and -99 etc
		int how_often = -1;

		auto generator1 = new CglProbing;
		impl->generators.emplace_back(generator1);
		generator1->setUsingObjective(true);
		generator1->setMaxPass(1);
		generator1->setMaxPassRoot(5);
		// Number of unsatisfied variables to look at
		generator1->setMaxProbe(10);
		generator1->setMaxProbeRoot(1000);
		// How far to follow the consequences
		generator1->setMaxLook(50);
		generator1->setMaxLookRoot(500);
		// Only look at rows with fewer than this number of elements
		generator1->setMaxElements(200);
		generator1->setRowCuts(3);
		
		impl->model->addCutGenerator(generator1, how_often, "Probing");

		auto generator2 = new CglGomory;
		impl->generators.emplace_back(generator2);
		// try larger limit
		generator2->setLimit(300);
		impl->model->addCutGenerator(generator2, how_often, "Gomory");

		auto generator3 = new CglKnapsackCover;
		impl->generators.emplace_back(generator3);
		impl->model->addCutGenerator(generator3, how_often, "Knapsack");

		auto generator4 = new CglRedSplit;
		impl->generators.emplace_back(generator4);
		// try larger limit
		generator4->setLimit(200);
		//impl->model->addCutGenerator(generator4, how_often, "RedSplit");

		auto generator5 = new CglClique;
		impl->generators.emplace_back(generator5);
		generator5->setStarCliqueReport(false);
		generator5->setRowCliqueReport(false);
		impl->model->addCutGenerator(generator5, how_often, "Clique");

		auto generator6 = new CglOddHole;
		impl->generators.emplace_back(generator6);
		generator6->setMinimumViolation(0.005);
		generator6->setMinimumViolationPer(0.00002);
		// try larger limit
		generator6->setMaximumEntries(200);
		//impl->model->addCutGenerator(generator6, how_often, "OddHole");

		auto mixedGen = new CglMixedIntegerRounding2;
		impl->generators.emplace_back(mixedGen);
		impl->model->addCutGenerator(mixedGen, how_often, "MixedIntegerRounding");
		
		auto flowGen = new CglFlowCover;
		impl->generators.emplace_back(flowGen);
		impl->model->addCutGenerator(flowGen, how_often, "FlowCover");

		int numberGenerators = impl->model->numberCutGenerators();
		for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
			CbcCutGenerator * generator = impl->model->cutGenerator(iGenerator);
			generator->setTiming(true);
		}

		/*
		CbcRounding heuristic1(*impl->model);
		impl->model->addHeuristic(&heuristic1);
		// And local search when new solution found
		CbcHeuristicLocal heuristic2(*impl->model);
		impl->model->addHeuristic(&heuristic2);
		
		// Do initial solve to continuous
		impl->model->initialSolve();

		// Could tune more
		double objValue = impl->model->solver()->getObjSense() * impl->model->solver()->getObjValue();
		double minimumDropA=CoinMin(1.0,fabs(objValue)*1.0e-3+1.0e-4);
		double minimumDrop= fabs(objValue)*1.0e-4+1.0e-4;
		printf("min drop %g (A %g)\n",minimumDrop,minimumDropA);
		impl->model->setMinimumDrop(minimumDrop);
		*/

		/*
		// Default strategy will leave cut generators as they exist already
		// so cutsOnlyAtRoot (1) ignored
		// numberStrong (2) is 5 (default)
		// numberBeforeTrust (3) is 5 (default is 0)
		// printLevel (4) defaults (0)
		CbcStrategyDefault strategy(true,5,5);
		// Set up pre-processing to find sos if wanted
		strategy.setupPreProcessing(2);
		impl->model->setStrategy(strategy);
		*/

		// Do complete search
		impl->model->branchAndBound();
	}

	return impl->parse_solution();
}

void IP::allow_ignoring_cost_function()
{
	impl->allow_ignoring_cost_function = true;
}

void add_at_most_k_constraint_binomial(Minisat::Solver* solver,
                                       const vector<Minisat::Lit>& literals,
                                       int k)
{
	vector<int> index_set;
	for (size_t ix = 0; ix < literals.size(); ++ix) {
		index_set.emplace_back(ix);
	}

	vector<vector<int>> subsets;
	generate_subsets(index_set, k + 1, &subsets);

	Minisat::vec<Minisat::Lit> clause(k + 1);
	for (auto& subset: subsets) {
		int clause_index = 0;
		for (int ix: subset) {
			clause[clause_index++] = ~ literals[ix];
		}
		solver->addClause(clause);
	}
}

void add_at_most_k_constraint(Minisat::Solver* solver,
                              const vector<Minisat::Lit>& X,
                              int k)
{
	attest(k >= 0);
	auto n = X.size();

	// Perhaps the binomial is better when e.g. n <= 7?
	if (n <= 1 || k == 0) {
		add_at_most_k_constraint_binomial(solver, X, k);
		return;
	}

	// This implementation follows
	// Carsten Sinz,
	// “Towards an Optimal CNF Encoding of Boolean Cardinality Constraints,”
	// Principles and Practice of Constraint Programming, 2005.
	// §2, page 2.

	vector<vector<Minisat::Lit>> s;
	for (size_t i = 0; i < n; ++i) {
		s.emplace_back();
		for (size_t j = 0; j < k; ++j) {
			s.back().emplace_back(Minisat::mkLit(solver->newVar()));
		}
	}

	// (¬x1 ∨ s1,1)
	solver->addClause(~X[0], s[0][0]);

	// for 1 < j ≤ k
	for (size_t j = 1; j < k; ++j) {
		// (¬s1, j) 
		solver->addClause(~s[0][j]);
	}

	// for 1 < i < n
	for (size_t i = 1; i < n - 1; ++i) {
		// (¬xi ∨ si,1)
		solver->addClause(~X[i], s[i][0]);
		// (¬si−1,1 ∨ si,1)
		solver->addClause(~s[i-1][0], s[i][0]);

		// for 1 < j ≤ k
		for (size_t j = 1; j < k; ++j) {
			// (¬xi ∨ ¬si−1,j−1 ∨ si,j)
			solver->addClause(~X[i], ~s[i-1][j-1], s[i][j]);
			// (¬si−1,j ∨ si,j)
			solver->addClause(~s[i-1][j], s[i][j]);
		}
		// (¬xi ∨ ¬si−1,k)
		solver->addClause(~X[i], ~s[i-1][k-1]);
	}
	// (¬xn ∨ ¬sn−1,k)
	solver->addClause(~X[n-1], ~s[n-2][k-1]);
}

bool IP::Implementation::solve_minisat()
{
	using namespace std;

	minisat_solver.reset(new Minisat::Solver);
	literals.clear();
	for (size_t j = 0; j < cost.size(); ++j) {
		auto lb = var_lb.at(j);
		auto ub = var_ub.at(j);
		check( (lb == 0 || lb == 1) && (ub == 0 || ub == 1),
		       "SAT solver requires boolean variables.");
		attest(lb <= ub);

		literals.push_back(Minisat::mkLit(minisat_solver->newVar()));
		if (lb == 1) {
			minisat_solver->addClause(literals.back());
		}
		if (ub == 0) {
			minisat_solver->addClause( ~literals.back());
		}

		if (!allow_ignoring_cost_function) {
			if (cost.at(j) != 0) {
				int icost = cost.at(j);
				check(icost == cost.at(j), "SAT requires integer costs.");

				// Add new literals equivalent to the variable and add them to
				// the vector of cost literals.
				for (int count = 1; count <= std::abs(icost); ++count) {
					auto lit = Minisat::mkLit(minisat_solver->newVar());
					minisat_solver->addClause(literals.back(), ~lit);
					minisat_solver->addClause(~literals.back(), lit);
					if (icost < 0) {
						lit = ~lit;
						sat_objective_offset -= 1;
					}
					objective_function_literals.emplace_back(lit);
				}
			}
		}
	}

	if (objective_function_literals.size() > 0) {
		// An objective function of
		//   3x + y
		// is modelled as 
		//   x1 + x2 + x3 + y1
		// where
		//   x1 ⇔ x
		//   x2 ⇔ x
		//   x3 ⇔ x
		//   y1 ⇔ y.
		// Then slack variables are added with an upper bound:
		//   x1 + x2 + x3 + y1 + s1 + s2 + s3 + s4 ≥ 4.
		// By assuming a different number of slack variables = 1, different
		// objective functions value can be tested for satisfiability.

		// First, we need to add the slack literals to the objective function.
		for (size_t i = 0; i < objective_function_literals.size(); ++i) {
			auto lit = Minisat::mkLit(minisat_solver->newVar());
			objective_function_slack_literals.emplace_back(lit);
		}

		std::vector<Minisat::Lit> objective_clause;
		for (auto& lit: objective_function_literals) {
			objective_clause.emplace_back(lit);
		}
		for (auto& lit: objective_function_slack_literals) {
			objective_clause.emplace_back(lit);
		}
		add_at_most_k_constraint(minisat_solver.get(), objective_clause, objective_function_literals.size());
	}

	auto num_constraints = rhs_lower.size();
	std::vector<int> lower(num_constraints);
	std::vector<int> upper(num_constraints);
	for (size_t i = 0; i < num_constraints; ++i) {
		auto to_int = [](double rhs)
		{
			const int limit = 1000 * 1000 * 1000;
			if (rhs > limit) {
				return limit;
			}
			else if (rhs < -limit) {
				return -limit;
			}
			else {
				int irhs = rhs;
				attest(rhs == irhs);
				return irhs;
			}
		};

		lower[i] = to_int(rhs_lower.at(i));
		upper[i] = to_int(rhs_upper.at(i));
	}

	vector<vector<Minisat::Lit>> lit_rows(num_constraints);
	for (size_t ind = 0; ind < rows.size(); ++ind) {
		auto var = literals.at(cols.at(ind));
		auto coeff = values.at(ind);
		check(coeff == 1 || coeff == -1, "SAT solver requires constraint coefficients of +-1.");

		if (coeff == 1) {
			lit_rows.at(rows.at(ind)).emplace_back(var);
		}
		else {
			lower[rows.at(ind)] += 1;
			upper[rows.at(ind)] += 1;
			lit_rows.at(rows.at(ind)).emplace_back( ~ var);
		}
	}

	vector<int> index_set;
	vector<vector<int>> subsets;
	for (size_t i = 0; i < num_constraints; ++i) {

		auto num_literals = lit_rows.at(i).size();

		index_set.clear();
		for (size_t ix = 0; ix < num_literals; ++ix) {
			index_set.emplace_back(ix);
		}

		if (lower[i] > 0) {
			auto neg_lit_row = lit_rows[i];
			for (auto& lit: neg_lit_row) {
				lit = ~ lit;
			}
			add_at_most_k_constraint(minisat_solver.get(), neg_lit_row, neg_lit_row.size() - lower[i]);
		}

		if (upper[i] < num_literals) {
			add_at_most_k_constraint(minisat_solver.get(), lit_rows[i], upper[i]);
		}
	}

	solution.clear();
	return next_minisat();
}

bool IP::Implementation::next_minisat()
{
	attest(minisat_solver);

	if (!solution.empty()) {
		// Forbid previous solution.
		attest(solution.size() == literals.size());
		Minisat::vec<Minisat::Lit> negated_solution;
		for (size_t j = 0; j < solution.size(); ++j) {
			if (solution[j] == 1) {
				negated_solution.push( ~literals[j]);
			}
			else {
				negated_solution.push(literals[j]);	
			}
		}
		minisat_solver->addClause(negated_solution);
	}

	bool result = minisat_solver->solve();
	if (!result) {
		solution.clear();
		return false;
	}

	if (solution.empty() && objective_function_literals.size() > 0) {
		int upper = objective_function_slack_literals.size();
		int lower = 0;
		bool ok;
		
		do {
			int current = (lower + upper) / 2;
			std::clog << "Objective value in [" << lower + sat_objective_offset << ", " << upper + sat_objective_offset << "]." << std::endl;
			if (lower >= upper) {
				break;
			}
			std::clog << "-- Trying " << current + sat_objective_offset << "... ";

			Minisat::vec<Minisat::Lit> assumptions;
			for (int i = 1; i <= objective_function_literals.size() - current; ++i) {
				assumptions.push(objective_function_slack_literals[i]);
			}
			ok = minisat_solver->solve(assumptions);

			if (ok) {
				upper = current;
				std::clog << "SAT." << std::endl;
			}
			else {
				lower = current + 1;
				std::clog << "UNSAT." << std::endl;
			}
		} while (true);

		// Add this objective as an assumption when resolving.
		std::vector<Minisat::Lit> neg_objective_function_literals = objective_function_literals;
		for (auto& lit: neg_objective_function_literals) {
			lit = ~lit;
		}
		add_at_most_k_constraint(minisat_solver.get(),     objective_function_literals, upper);
		add_at_most_k_constraint(minisat_solver.get(), neg_objective_function_literals, objective_function_literals.size() - upper);
	}

	//minisat_solver->printStats();

	solution.clear();
	for (size_t j = 0; j < cost.size(); ++j) {
		auto value = minisat_solver->modelValue(literals.at(j));
		attest(value == Minisat::l_True || value == Minisat::l_False);
		solution.push_back(value == Minisat::l_True ? 1 : 0);
	}

	// Check feasibility just to make sure everything is alright.
	auto num_constraints = rhs_lower.size();
	vector<double> row_sums(num_constraints, 0);
	for (size_t ind = 0; ind < rows.size(); ++ind) {
		auto var   = solution.at(cols.at(ind));
		auto coeff = values.at(ind);
		row_sums.at(rows.at(ind)) += coeff * var;
	}
	for (size_t i = 0; i < num_constraints; ++i) {
		auto lower = rhs_lower.at(i);
		auto upper = rhs_upper.at(i);
		attest(lower - 1e-9   <= row_sums.at(i));
		attest(row_sums.at(i) <= upper + 1e-9);
	}

	return true;
}

bool IP::next_solution()
{
	check( ! impl->integer_variables.empty(), "next_solution(): Need integer variables.");

	if (impl->external_solver == Minisat) {
		return impl->next_minisat();
	}

	OsiSolverInterface * refSolver = nullptr;
	OsiSolverInterface* solver = nullptr;
	const double * objective = nullptr;

	if (impl->use_osi()) {
		refSolver = impl->problem.get();
		solver = impl->problem.get();

		objective = solver->getObjCoefficients();
	}
	else {
		attest(impl->model);
		refSolver = impl->model->referenceSolver();
		solver = impl->model->solver();

		objective = refSolver->getObjCoefficients();	
	}

	//
	// We add two new rows to the problem in order to get the
	// next solution. If the current solution is x = (1, 0, 1),
	//
	//    (a)  (1 - x1) + x2 + (1 - x3) >= 1
	//    (b)  objective(x) == *optimal*.
	//
	CoinPackedVector solution_cut, objective_cut;
	double best_objective = 0;
	double solution_rhs = 1.0;
	for (int iColumn = 0; iColumn < impl->solution.size(); iColumn++) {
		double value = impl->solution[iColumn];
		if (solver->isInteger(iColumn)) {
			// only works for 0-1 variables
			attest(impl->var_lb[iColumn] == 0.0 || impl->var_lb[iColumn] == 1.0);
			attest(impl->var_ub[iColumn] == 0.0 || impl->var_ub[iColumn] == 1.0);
			// double check integer
			attest (fabs(floor(value+0.5)-value)<1.0e-5);
			if (value>0.5) {
			// at 1.0
			solution_cut.insert(iColumn,-1.0);
			solution_rhs -= 1.0;
			} else {
				// at 0.0
				solution_cut.insert(iColumn,1.0);
			}
		}

		best_objective += value * objective[iColumn];
		objective_cut.insert(iColumn, objective[iColumn]);
		refSolver->setObjCoeff(iColumn, 0.0);
	}

    // now add cut
	refSolver->addRow(solution_cut, solution_rhs, COIN_DBL_MAX);
	refSolver->addRow(objective_cut, best_objective, best_objective);

	if (impl->use_osi()) {
		refSolver->branchAndBound();
	}
	else {
		impl->model->resetToReferenceSolver();
		//impl->model->setHotstartSolution(bestSolution, nullptr);

		// Do complete search
		impl->model->branchAndBound();
	}

	return impl->parse_solution();
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
		add_constraint(active_in_window <= N);
		constraints_added++;
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

			add_constraint(constraint >= -window_size + 1);
			constraints_added++;
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

void internal_subset(const std::vector<int>& set, int left, int index, std::vector<int>* scratch_space, std::vector<std::vector<int>>* all_subsets){
	if (left == 0){
		all_subsets->push_back(*scratch_space);
		return;
	}
	if (left > set.size() - index) {
		// We don’t have enough elements left to create a subset.
		return;
	}
	for (std::size_t i = index; i < set.size(); i++){
		scratch_space->push_back(set[i]);
		internal_subset(set, left - 1, i + 1, scratch_space, all_subsets);
		scratch_space->pop_back();
	}
}

size_t choose(size_t n, size_t k)
{
	if (k == 0) {
		return 1;
	}
	return  (n * choose(n - 1, k - 1)) / k;
}

void generate_subsets(const std::vector<int>& set, int subset_size, std::vector<std::vector<int>>* output)
{
	size_t num_subsets = choose(set.size(), subset_size);
	if (num_subsets > 50000000) {
		// Maybe change this limit in the future.
		throw std::runtime_error("Too many subsets. Choose a better algorithm.");
	}

	output->clear();
	output->reserve(num_subsets);
	std::vector<int> scratch_space;
	scratch_space.reserve(subset_size);
	internal_subset(set, subset_size, 0, &scratch_space, output);
}
