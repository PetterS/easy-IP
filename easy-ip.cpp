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
	// Copying neccessary.
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
		: creator(creator_), external_solver(IP::Default)
	{ }

	bool parse_solution();
	
	bool use_osi() const;

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

	bool preprocess;
	std::unique_ptr<OsiSolverInterface> problem;
	std::unique_ptr<CbcModel> model;
	std::vector<std::unique_ptr<CglCutGenerator>> generators;

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

std::unique_ptr<OsiSolverInterface> IP::get_problem(std::unique_ptr<OsiSolverInterface> existing_solver)
{
	std::unique_ptr<OsiSolverInterface> problem;

	if (existing_solver) {
		// Take ownership.
		problem = std::move(existing_solver);
	}
	else {
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

	return std::move(problem);
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


bool IP::solve(const CallBack& callback_function)
{
	if(impl->external_solver == IP::CPLEX) {
		#ifdef HAS_CPLEX
			auto cplex_solver = std::unique_ptr<OsiSolverInterface>(new OsiCpxSolverInterface);
			impl->problem = get_problem(std::move(cplex_solver));
		#endif
	}
	else if(impl->external_solver == IP::MOSEK) {
		#ifdef HAS_MOSEK
			auto mosek_solver = std::unique_ptr<OsiSolverInterface>(new OsiMskSolverInterface);
			impl->problem = get_problem(std::move(mosek_solver));
		#endif
	}
	else {
		impl->problem = get_problem();
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

bool IP::next_solution()
{
	check( ! impl->integer_variables.empty(), "next_solution(): Need integer variables.");

	OsiSolverInterface * refSolver = nullptr;
	OsiSolverInterface* solver = nullptr;
	const double * objective = nullptr;

	if (impl->use_osi()) {
		refSolver = impl->problem.get();
		solver = impl->problem.get();

		objective = solver->getObjCoefficients();
	}
	else {
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
