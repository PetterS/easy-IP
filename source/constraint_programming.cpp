// Petter Strandmark 2014
// petter.strandmark@gmail.com

#include <gecode/int.hh>
#include <gecode/search.hh>

#include <easy-ip.h>
#include <easy-ip-internal.h>

namespace
{
	int gecode_round(double value)
	{
		int intval = [&]()
			{
				if (value >= 0) {
					return int(value + 0.5);
				}
				else {
					return int(value - 0.5);
				}
			}();
		check(intval == value, "Constraint programming requires integers.");
		return intval;
	}
}

class GecodeIPModel:
	public Gecode::Space
{
protected:
	Gecode::IntVarArray vars;
public:
	GecodeIPModel(
		const vector<double>& rhs_lower,
		const vector<double>& rhs_upper,
		const vector<int>& rows,
		const vector<int>& cols,
		const vector<double>& values,
		const vector<double>& var_lb,
		const vector<double>& var_ub,
		const vector<double>& cost)	
	: vars(*this, int(var_lb.size()))
	{
		using namespace Gecode;
		using namespace std;

		auto n = var_lb.size();
		for (int j = 0; j < n; ++j) {
			vars[j] = IntVar(*this,
			                 gecode_round(var_lb[j]),
			                 gecode_round(var_ub[j]));
		}

		auto m = rhs_lower.size();
		vector<vector<int>> row_coeffs(m);
		vector<vector<int>> row_vars(m);

		for (size_t ind = 0; ind < rows.size(); ++ind) {
			auto row = rows.at(ind);
			auto col = cols.at(ind);
			auto val = gecode_round(values.at(ind));
			row_coeffs[row].emplace_back(val);
			row_vars[row].emplace_back(col);
		}

		for (size_t i = 0; i < m; ++i) {
			IntArgs c(int(row_coeffs[i].size()));
			IntVarArgs x(int(row_coeffs[i].size()));
			for (int ind = 0; ind < int(row_coeffs[i].size()); ++ind) {
				c[ind] = row_coeffs[i][ind];
				x[ind] = vars[row_vars[i][ind]];
			}

			if (rhs_lower[i] >= rhs_upper[i]) {
				linear(*this, c, x, IRT_EQ, gecode_round(rhs_lower[i]));
			}
			else {
				if (rhs_lower[i] > -1e9) {
					linear(*this, c, x, IRT_GQ, gecode_round(rhs_lower[i]));
				}
				if (rhs_upper[i] < 1e9) {
					linear(*this, c, x, IRT_LQ, gecode_round(rhs_upper[i]));
				}
			}
		}

		branch(*this, vars, INT_VAR_NONE(), INT_VAL_MIN());
	}

	// Constructor for cloning \a s
	GecodeIPModel(bool share, GecodeIPModel& s)
		: Space(share, s)
	{
		vars.update(*this, share, s.vars);
	}

	// Copy during cloning
	virtual Space* copy(bool share)
	{
		return new GecodeIPModel(share, *this);
	}

	void print() const
	{
		std::cout << vars << std::endl;
	}

	void get_solution(std::vector<double>* x) const
	{
		auto n = vars.size();
		x->resize(n);
		for (int j = 0; j < n; ++j) {
			(*x)[j] = vars[j].val();
		}
	}
};

bool IP::Implementation::solve_gecode()
{
	using namespace Gecode;

	auto model = new GecodeIPModel(rhs_lower, rhs_upper, rows, cols, values, var_lb, var_ub, cost);
	DFS<GecodeIPModel> dfs(model);
	delete model;

	if (GecodeIPModel* m = dfs.next()) {
		m->get_solution(&solution);
		delete m;
		return true;
	}

	return false;
}
