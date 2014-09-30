// Petter Strandmark 2014
// petter.strandmark@gmail.com
//
// Converts an integer program to a SAT problem.
//

#include <easy-ip.h>
#include <easy-ip-internal.h>

// The binomial approach is only effective for small k.
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
	convert_to_minisat();
	solution.clear();
	return next_minisat();
}

void IP::Implementation::convert_to_minisat()
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
		check(coeff >= 1 || coeff == -1, "SAT solver requires constraint coefficients of +-1.");

		if (coeff == 1) {
			lit_rows.at(rows.at(ind)).emplace_back(var);
		}
		else if (coeff > 1) {
			// Add new literals equivalent to the variable and add them to
			// the vector of cost literals.
			for (int count = 1; count <= coeff; ++count) {
				auto lit = Minisat::mkLit(minisat_solver->newVar());
				minisat_solver->addClause(var, ~lit);
				minisat_solver->addClause(~var, lit);
				lit_rows.at(rows.at(ind)).emplace_back(lit);
			}
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
