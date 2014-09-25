// Petter Strandmark 2014
// petter.strandmark@gmail.com

#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
//#include <typeid>
#include <vector>

// Using Clp as the solver
#include <coin/OsiClpSolverInterface.hpp>

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

#ifdef HAS_CPLEX
	#include <coin/OsiCpxSolverInterface.hpp>
#endif

#ifdef HAS_MOSEK
	#include <coin/OsiMskSolverInterface.hpp>
#endif

#include <easy-ip.h>
#include <easy-ip-internal.h>

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

bool IP::solve(const CallBack& callback_function, bool silent_mode)
{
	if (impl->external_solver == IP::Minisat) {
		return impl->solve_minisat();
	}
	else if (impl->external_solver == IP::Gecode) {
		return impl->solve_gecode();
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
}

bool IP::solve_relaxation()
{
	attest(impl->external_solver != IP::Minisat);
	std::vector<size_t> integer_variables_copy;
	integer_variables_copy.swap(impl->integer_variables);

	easyip_at_scope_exit(impl->integer_variables.swap(integer_variables_copy));
	easyip_at_scope_exit(impl->problem.release());

	return solve();
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

bool IP::Implementation::use_osi() const
{
	return external_solver != IP::Default || integer_variables.empty();
}
