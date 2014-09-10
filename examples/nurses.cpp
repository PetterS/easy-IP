// Petter Strandmark 2014.
//
// Reads and solves problem files from NSPLib
// http://www.projectmanagement.ugent.be/?q=research/personnel_scheduling/nsp
//
// First argument to the program is the problem .nsp file. The second
// argument is the .gen-file with the constraints to use.
// 
#include <fstream>
#include <string>
#include <vector>

#ifdef USE_OPENMP
	#include <omp.h>
#else
	namespace { double omp_get_wtime() { return 0; } }
#endif

#include <easy-ip.h>

class NurseProblem
{
public:
	NurseProblem(std::string problem_file, std::string cases_file)
	{
		std::ifstream problem(problem_file);
		attest(problem);

		problem >> num_nurses >> num_days >> num_shifts;
		attest(problem);

		assignments.resize(num_nurses, std::vector<int>(num_days, 0));

		for (int d = 0; d < num_days; ++d) {
			day_coverage.emplace_back();
			for (int s = 0; s < num_shifts; ++s) {
				int coverage = -1;
				problem >> coverage;
				attest(problem);
				day_coverage.back().emplace_back(coverage);
			}
		}

		for (int n = 0; n < num_nurses; ++n) {
			shift_preferences.emplace_back();
			for (int d = 0; d < num_days; ++d) {
				shift_preferences.back().emplace_back();
				for (int s = 0; s < num_shifts; ++s) {
					int preference = -1;
					problem >> preference;
					attest(problem);
					shift_preferences.back().back().emplace_back(preference);
				}
			}
		}

		// Check that there isn’t more information.
		int dummy = 0;
		problem >> dummy;
		attest(!problem);

		std::ifstream cases(cases_file);
		int num_days_tmp = -1;
		int num_shifts_tmp = -1;
		cases >> num_days_tmp >> num_shifts_tmp;
		attest(cases);
		attest(num_days_tmp == num_days);
		attest(num_shifts_tmp == num_shifts);

		cases >> min_assignments >> max_assignments;
		attest(cases);

		cases >> min_consequtive >> max_consequtive;
		attest(cases);

		for (int s = 0; s < num_shifts; ++s) {
			min_consequtive_shifts.emplace_back();
			cases >> min_consequtive_shifts.back();

			max_consequtive_shifts.emplace_back();
			cases >> max_consequtive_shifts.back();

			min_assignments_per_shift.emplace_back();
			cases >> min_assignments_per_shift.back();

			max_assignments_per_shift.emplace_back();
			cases >> max_assignments_per_shift.back();
			attest(cases);
		}

		// Check that there isn’t more information.
		cases >> dummy;
		attest(!cases);
	}

	int get_num_nurses() const
	{
		return num_nurses;
	}

	int get_num_days() const
	{
		return num_days;
	}

	int get_num_shifts() const
	{
		return num_shifts;
	}

	int get_coverage(int day, int shift) const
	{
		return day_coverage.at(day).at(shift);
	}

	int get_preference(int nurse, int day, int shift) const
	{
		return shift_preferences.at(nurse).at(day).at(shift);
	}

	int get_min_assignments() const { return min_assignments; }
	int get_max_assignments() const { return max_assignments; }

	int get_min_consequtive() const { return min_consequtive; }
	int get_max_consequtive() const { return max_consequtive; }
	int get_min_consequtive(int shift) const { return min_consequtive_shifts.at(shift); }
	int get_max_consequtive(int shift) const { return max_consequtive_shifts.at(shift); }

	int get_min_assignments_per_shift(int shift) const { return min_assignments_per_shift.at(shift); }
	int get_max_assignments_per_shift(int shift) const { return max_assignments_per_shift.at(shift); }

	void assign(int nurse, int day, int shift)
	{
		assignments.at(nurse).at(day) = shift;
	}

	int objective() const
	{
		int obj = 0;
		for (int n = 0; n < get_num_nurses(); ++n) {
			for (int d = 0; d < get_num_days(); ++d) {
				obj += get_preference(n, d, assignments.at(n).at(d));
			}
		}
		return obj;
	}

	void write_to_stream(std::ofstream& solution) const
	{
		solution << get_num_nurses() << " "
		         << get_num_days() << " "
		         << get_num_shifts() << std::endl;
		for (int n = 0; n < get_num_nurses(); ++n) {
			for (int d = 0; d < get_num_days(); ++d) {
				solution << assignments.at(n).at(d) << " ";
			}
			solution << std::endl;
		}
	}

protected:
	int num_nurses;
	int num_days;
	int num_shifts;

	int min_assignments = 0;
	int max_assignments = 1000000;

	int min_consequtive = 0;
	int max_consequtive = 1000000;

	std::vector<std::vector<int>> day_coverage;
	std::vector<std::vector<std::vector<int>>> shift_preferences;

	std::vector<int> min_assignments_per_shift;
	std::vector<int> max_assignments_per_shift;
	std::vector<int> min_consequtive_shifts;
	std::vector<int> max_consequtive_shifts;

	std::vector<std::vector<int>> assignments;
};

int main_program(int num_args, char* args[])
{
	using namespace std;

	attest(num_args >= 3);
	NurseProblem problem(args[1], args[2]);
	string solution_file_name = "";
	if (num_args >= 4) {
		solution_file_name = args[3];
	}

	string problem_number = args[1];
	auto p1 = problem_number.find_last_of("/\\");
	if (p1 != string::npos) {
		problem_number = problem_number.substr(p1 + 1);
	}
	auto p2 = problem_number.find(".");
	if (p2 != string::npos) {
		problem_number = problem_number.substr(0, p2);
	}

	clog << "Number of nurses: " << problem.get_num_nurses() << endl;
	clog << "Number of days  : " << problem.get_num_days() << endl;
	clog << "Number of shifts: " << problem.get_num_shifts() << endl;

	double start_time = omp_get_wtime();

	IP ip;

	// The assignment variables.
	//   x[n][d][s] = 1 ⇔ Nurse n is working day d, shift s.
	auto x = ip.add_boolean_cube(problem.get_num_nurses(), problem.get_num_days(), problem.get_num_shifts());

	// Returns a variable (Sum) that is 1 iff the nurse is working
	// on a particular day.
	auto working_on_day = [&x, &problem](int nurse, int day)
	{
		Sum working_on_day = 0;
		// Do not count the empty shift here, hence - 1.
		for (int s = 0; s < problem.get_num_shifts() - 1; ++s) {
			working_on_day += x[nurse][day][s];
		}
		return working_on_day;
	};

	// Exactly one shift per day per nurse (last shift is the empty shift).
	for (int n = 0; n < problem.get_num_nurses(); ++n) {
		for (int d = 0; d < problem.get_num_days(); ++d) {
			Sum num_day_shifts = 0;
			for (int s = 0; s < problem.get_num_shifts(); ++s) {
				num_day_shifts += x[n][d][s];
			}
			ip.add_constraint(num_day_shifts == 1);
		}
	}

	// Shift coverage.
	for (int d = 0; d < problem.get_num_days(); ++d) {
		for (int s = 0; s < problem.get_num_shifts(); ++s) {
			Sum shift_workers = 0;
			for (int n = 0; n < problem.get_num_nurses(); ++n) {
				shift_workers += x[n][d][s];
			}
			ip.add_constraint(shift_workers >= problem.get_coverage(d, s));
		}
	}

	// Assignments per nurse.
	for (int n = 0; n < problem.get_num_nurses(); ++n) {
		Sum num_assignments = 0;
		for (int d = 0; d < problem.get_num_days(); ++d) {
			num_assignments += working_on_day(n, d);
		}
		ip.add_constraint(problem.get_min_assignments(), num_assignments, problem.get_max_assignments());
	}

	// Min/max consequtive days working.
	{
		int constraints_added = 0;
		for (int n = 0; n < problem.get_num_nurses(); ++n) {

			vector<Sum> working_on_days;
			for (int d = 0; d < problem.get_num_days(); ++d) {
				working_on_days.emplace_back(working_on_day(n, d));
			}

			// I think this constraint should have ’false,’ i.e. singles at borders can be disallowed.
			constraints_added +=
				ip.add_min_consequtive_constraints(problem.get_min_consequtive(), working_on_days, false);

			constraints_added +=
				ip.add_max_consequtive_constraints(problem.get_max_consequtive(), working_on_days);
		}
		clog << "Added " << constraints_added << " constraints for minimum/maximum consequtive days." << endl;
	}

	// Min/max consequtive days working a specific shift.
	{
		int constraints_added = 0;
		for (int s = 0; s < problem.get_num_shifts(); ++s) {
			for (int n = 0; n < problem.get_num_nurses(); ++n) {

				vector<Sum> working_on_days;
				for (int d = 0; d < problem.get_num_days(); ++d) {
					working_on_days.emplace_back(x[n][d][s]);
				}

				// I think this constraint should have ’false,’ i.e. singles at borders can be disallowed.
				constraints_added +=
					ip.add_min_consequtive_constraints(problem.get_min_consequtive(s), working_on_days, false);

				constraints_added +=
					ip.add_max_consequtive_constraints(problem.get_max_consequtive(s), working_on_days);
			}
		}
		if (constraints_added > 0) {
			clog << "Added " << constraints_added << " constraints for minimum/maximum consequtive days of specific shifts." << endl;
		}
	}

	// Assignments per shift.
	for (int n = 0; n < problem.get_num_nurses(); ++n) {
		for (int s = 0; s < problem.get_num_shifts() - 1; ++s) {
			Sum num_assignments = 0;
			for (int d = 0; d < problem.get_num_days(); ++d) {
				num_assignments += x[n][d][s];
			}
			ip.add_constraint(problem.get_min_assignments_per_shift(s), num_assignments, problem.get_max_assignments_per_shift(s));
		}
	}

	// Preferences.
	Sum preferences = 0;
	for (int n = 0; n < problem.get_num_nurses(); ++n) {
		for (int d = 0; d < problem.get_num_days(); ++d) {
			for (int s = 0; s < problem.get_num_shifts(); ++s) {
				preferences += problem.get_preference(n, d, s) * x[n][d][s];
			}
		}
	}
	ip.add_objective(preferences);

	double startup_time = omp_get_wtime() - start_time;

	clog << "IP created." << endl;

	auto assign_solution = [&x, &problem]()
	{
		for (int n = 0; n < problem.get_num_nurses(); ++n) {
			for (int d = 0; d < problem.get_num_days(); ++d) {
				for (int s = 0; s < problem.get_num_shifts(); ++s) {
					if (x[n][d][s].value()) {
						problem.assign(n, d, s);
					}
				}
			}
		}
	};

	start_time = omp_get_wtime();

	// Try SAT solver first.
	bool try_SAT_solver = true;
	if (try_SAT_solver) {
		start_time = omp_get_wtime();
		ip.set_external_solver(IP::Minisat);
		ip.allow_ignoring_cost_function();
		ip.solve();
		int num_solutions = 0;
		do {
			double current_time = omp_get_wtime() - start_time;
			num_solutions++;
			assign_solution();
			clog << "SAT solver found solution of " << problem.objective() << " in " << current_time << " seconds." << endl;
			if (num_solutions >= 1) {  // Increase for more solutions.
				break;
			}
		} while (ip.next_solution());
		ip.set_external_solver(IP::Default);
	}

	// Solve in completely silent mode to keep stdout clean.
	ip.set_time_limit(5.0);

	if (ip.solve(nullptr, true)) {
		clog << "IP solved." << endl;
		assign_solution();
		attest(preferences.value() == problem.objective());
	}
	else if (!try_SAT_solver) {
		clog << "IP was not solved." << endl;
		cout << problem_number << "\t"
		     << "infeasible/timelimit" << "\t"
		     << 0 << endl;
		return 2;
	}

	double solver_time = omp_get_wtime() - start_time;

	if (solution_file_name != "") {
		ofstream solution(solution_file_name);
		problem.write_to_stream(solution);
	}

	clog << "Obtained objective function " << problem.objective() << " in " << solver_time << " seconds." << endl;

	cout << problem_number << "\t"
	     << problem.objective() << "\t"
	     << startup_time + solver_time << endl;
	
	clog << endl;
	return 0;
}

int main(int num_args, char* args[])
{
	try {
		return main_program(num_args, args);
	}
	catch (std::exception& err) {
		std::cerr << "Error: " << err.what() << std::endl;
		return 1;
	}
}
