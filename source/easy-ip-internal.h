#ifndef EASY_IP_INTERNAL_HEADER
#define EASY_IP_INTERNAL_HEADER

#include <coin/CbcModel.hpp>
#include <coin/CbcCutGenerator.hpp>

#include "minisat/core/Solver.h"

#include <easy-ip.h>

void check(bool expr, const char* message);
void assertion_failed(const char* expr, const char* file_cstr, int line);

//
// easyip_at_scope_exit( statement; ) executes statement at the end 
// of the current scope.
//
template <typename F>
class ScopeGuard
{
public:
    ScopeGuard(F&& f)
		: f(std::forward<F>(f))
	{}

	ScopeGuard(ScopeGuard&& guard)
		: f(std::move(guard.f)), active(guard.active)
	{
		guard.dismiss();
	}

    ~ScopeGuard()
	{
		if (active) { 
			f();
		}
	}

	ScopeGuard(const ScopeGuard&) = delete;
	ScopeGuard& operator = (const ScopeGuard&) = delete;

	void dismiss()
	{
		active = false;
	}

private:
    F f;
	bool active = true;
};

template <typename F>
ScopeGuard<F> make_scope_guard(F&& f) {
    return std::move(ScopeGuard<F>(std::forward<F>(f)));
};

#define EASYIP_JOIN_PP_SYMBOLS_HELPER(arg1, arg2) arg1 ## arg2
#define EASYIP_JOIN_PP_SYMBOLS(arg1, arg2) EASYIP_JOIN_PP_SYMBOLS_HELPER(arg1, arg2)
#define easyip_at_scope_exit(code) \
    auto EASYIP_JOIN_PP_SYMBOLS(easyip_scope_exit_guard_, __LINE__) = ::make_scope_guard([&](){code;})


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
		check(variable.creator != nullptr, "Variables used in sums must be created by an IP object.");
		cols.push_back(int(variable.index));
		values.push_back(1.0);
	}

	double constant;
	vector<int> cols;
	vector<double> values;

	const IP* creator;
};

class GecodeIPModel;
class GecodeContainer
{
public:
	GecodeContainer();
	GecodeContainer(GecodeIPModel* model);
	GecodeContainer(const GecodeContainer&)  = delete;
	void operator = (const GecodeContainer&) = delete;
	GecodeContainer(GecodeContainer&&)       = delete;
	void operator = (GecodeContainer&&);
	~GecodeContainer();
	bool next_solution(std::vector<double>* x);
private:
	class Implementation;
	Implementation* impl;
};

class IP::Implementation
{
public:
	Implementation(const IP* creator_)
		: external_solver(IP::Default), creator(creator_)
	{ }

	bool parse_solution();
	
	bool use_osi() const;

	void convert_to_minisat();
	bool solve_minisat();
	bool next_minisat();

	bool solve_gecode();

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

	GecodeContainer gecode_container;

	std::vector<std::unique_ptr<CglCutGenerator>> generators;

	bool allow_ignoring_cost_function = false;

	void check_creator(const Variable& t) const;
	void check_creator(const Sum& t) const;

	const IP* creator;
};

#endif
