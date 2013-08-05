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

Sum::Sum()
	: constant(0.0),
	  creator(nullptr)
{ }

Sum::Sum(const Sum& sum) 
{
	*this = sum;
}

Sum& Sum::operator = (const Sum& sum)
{
	cols = sum.cols;
	values = sum.values;
	constant = sum.constant;
	creator = sum.creator;
	//std::cerr << "Copy sum.\n";
	//throw std::runtime_error("Copying took place.");
	return *this;
}

Sum::Sum(Sum&& sum)
{
	*this = std::move(sum);
}

Sum& Sum::operator = (Sum&& sum)
{
	cols = std::move(sum.cols);
	values = std::move(sum.values);
	constant = sum.constant;
	creator = sum.creator;
	//std::cerr << "Move sum.\n";
	return *this;
}

Sum::Sum(double constant_)
	: constant(constant_),
	  creator(nullptr)
{ }

Sum::Sum(const Variable& variable)
	: constant(0.0), 
	  creator(variable.creator)
{
	cols.push_back(int(variable.index));
	values.push_back(1.0);
}

Sum::~Sum()
{ }

void Sum::add_term(double coeff, const Variable& variable)
{
	cols.push_back(static_cast<int>(variable.index));
	values.push_back(coeff);
}

double Sum::value() const
{
	return creator->get_solution(*this);
}

Sum& Sum::operator += (const Sum& rhs)
{
	match_solvers(rhs);

	constant += rhs.constant;
	for (size_t i = 0; i < rhs.cols.size(); ++i) {
		cols.push_back(rhs.cols[i]);
		values.push_back(rhs.values[i]);
	}
	return *this;
}

Sum& Sum::operator -= (const Sum& rhs)
{
	match_solvers(rhs);

	constant -= rhs.constant;
	for (size_t i = 0; i < rhs.cols.size(); ++i) {
		cols.push_back(rhs.cols[i]);
		values.push_back(-rhs.values[i]);
	}
	return *this;
}

Sum& Sum::operator *= (double coeff)
{
	if (coeff == 0.0) {
		values.clear();
		cols.clear();
		constant = 0.0;
		return *this;
	}

	for (auto& value : values) {
		value *= coeff;
	}
	constant *= coeff;
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
	constant = -constant;
	for (auto& value: values) {
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
	check(creator == nullptr || sum.creator == nullptr || creator == sum.creator,
	      "Variables from different solver can not be mixed.");
	if (creator == nullptr) {
		creator = sum.creator;
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
	: lower_bound(1), upper_bound(static_cast<double>(expression.sum.values.size())), sum(expression.sum)
{ }

Constraint::Constraint(LogicalExpression&& expression)
	: lower_bound(1), upper_bound(static_cast<double>(expression.sum.values.size())), sum(std::move(expression.sum))
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
	Implementation()
		: external_solver(IP::Default)
	{ }

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

	template<typename T>
	void check_creator(const T& t) const;
};

IP::IP() 
	: impl(new Implementation)
{
};

IP::IP(IP&& lhs) 
	: impl(lhs.impl)
{
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

	impl->rhs_lower.push_back(L - sum.constant);
	impl->rhs_upper.push_back(U - sum.constant);
	auto row_index = impl->rhs_upper.size() - 1;
	for (size_t i = 0; i < sum.cols.size(); ++i) {
		impl->rows.push_back(static_cast<int>(row_index));
		impl->cols.push_back(sum.cols[i]);
		impl->values.push_back(sum.values[i]);
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

	for (size_t i = 0; i < sum.cols.size(); ++i) {
		impl->cost.at(sum.cols[i]) += sum.values[i];
	}
	// TODO: handle constant.
}

template<typename T>
void IP::Implementation::check_creator(const T& t) const
{
	check(t.creator == nullptr || t.creator->impl == this,
	      "Variable comes from a different solver.");
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

	double value = sum.constant;
	for (size_t i = 0; i < sum.cols.size(); ++i) {
		value += sum.values[i] * impl->solution.at(sum.cols[i]);
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

					if (process) {
						org_n = process->originalModel()->getNumCols();
						org_columns = process->originalColumns();
					}
					else {
						org_n = n;
						org_columns = model->originalColumns();
					}

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
	               const CbcModel* model_,
	               CglPreProcess* process_)
		: callback_function(callback_function_),
		  solution(solution_),
		  model(model_),
		  process(process_)
	{ }

	MyEventHandler(const MyEventHandler& rhs)
		: CbcEventHandler(rhs),
		  callback_function(rhs.callback_function),
		  solution(rhs.solution),
		  model(rhs.model),
		  process(rhs.process)
	{ }

	virtual ~MyEventHandler()
	{ }

	MyEventHandler& operator=(const MyEventHandler & rhs)
	{
		if (this != &rhs) {
			callback_function = rhs.callback_function;
			solution = rhs.solution;
			model = rhs.model;
			process = rhs.process;
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
	CglPreProcess* process;
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

bool IP::solve(const CallBack& callback_function)
{
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

	std::unique_ptr<OsiSolverInterface> problem;
	if(impl->external_solver == IP::CPLEX) {
		#ifdef HAS_CPLEX
			problem.reset(new OsiCpxSolverInterface);
		#endif
	}
	else if(impl->external_solver == IP::MOSEK) {
		#ifdef HAS_MOSEK
			problem.reset(new OsiMskSolverInterface);
		#endif
	}
	else {
		auto clp_problem = new OsiClpSolverInterface;
		problem.reset(clp_problem);
		// Turn off information from the LP solver.
		clp_problem->setLogLevel(0);
	}

	problem->loadProblem(coinMatrix,
	                     &impl->var_lb[0],
	                     &impl->var_ub[0],
	                     &impl->cost[0],
	                     &impl->rhs_lower[0],
	                     &impl->rhs_upper[0]);

	for (auto index: impl->integer_variables) {
		problem->setInteger(static_cast<int>(index));
	}

	std::unique_ptr<CglPreProcess> process(nullptr);
	OsiSolverInterface* preprocessed_problem;
	OsiSolverInterface* solved_problem;

	if (impl->external_solver != IP::Default) {
		problem->branchAndBound();
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
		else if (!impl->integer_variables.empty() && !solved_problem->isProvenOptimal()) {
			std::cerr << "-- Not optimal." << std::endl;
			return false;
		}
	}
	else {
		bool preprocess = true;
		if (callback_function) {
			// Intermediate solutions are not correct
			// when using preprocessing.
			preprocess = false;
		}

		if (preprocess) {
			process.reset(new CglPreProcess);
			problem->initialSolve();
			preprocessed_problem = process->preProcess(*problem, false, 20);
			if (!preprocessed_problem) {
				//throw std::runtime_error("Problem infeasible.");
				return false;
			}
			preprocessed_problem->resolve();
		}
		else {
			preprocessed_problem = problem.get();
		}

		// Pass the solver with the problem to be solved to CbcModel 
		CbcModel model(*preprocessed_problem);

		// Only the most important log messages.
		model.setLogLevel(1);

		if (callback_function) {
			MyEventHandler my_event_handler(callback_function,
			                                &impl->solution,
			                                &model,
			                                process.get());
			model.passInEventHandler(&my_event_handler);
		}

		// Do complete search
		model.branchAndBound();

		if (model.isProvenInfeasible()) {
			//throw std::runtime_error("Problem infeasible.");
			return false;
		}
	
		if (model.isProvenDualInfeasible()) {
			//throw std::runtime_error("Problem unbounded.");
			return false;
		}

		if (!model.isProvenOptimal()) {
			throw std::runtime_error("Could not solve IP.");
		}

		if (preprocess) {
			process->postProcess(*model.solver());
			solved_problem = problem.get();
		}
		else {
			solved_problem = model.solver();
		}
	}

	int numberColumns = solved_problem->getNumCols();
	const double* raw_solution = solved_problem->getColSolution();

	impl->solution.clear();
	for (size_t i = 0; i < numberColumns; ++i) {
		impl->solution.push_back(raw_solution[i]);
	}

	return true;
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
