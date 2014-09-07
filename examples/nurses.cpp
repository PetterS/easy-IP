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

#include <omp.h>

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
};


// Adds constraints forbidding less than N consequtive variables in a row. 
int add_min_consequtive_constraints(IP* ip, int N, const std::vector<Sum>& variables)
{
	int constraints_added = 0;

	for (int window_size = 1; window_size <= N - 1; ++window_size) {
		// Look for windows of size minimum - 1 along with the
		// surrounding slots.
		//  
		// [*] [x1] [y1] [x2] [*]
		//
		// x1 = 0 ∧ x2 = 0 ⇒ y1 = 0
		// ⇔
		// x1 + x2 - y1 ≥ 0
		//
		// Then add windows with more y variables.

		for (int window_start = 0; window_start < variables.size() - window_size; ++window_start) {

			Sum constraint = 0;

			if (window_start - 1 >= 0) {
				constraint += variables.at(window_start - 1);
			}

			for (int i = window_start; i < window_start + window_size; ++i) {
				constraint -= variables.at(i);
			}

			if (window_start + window_size < variables.size()) {
				constraint += variables.at(window_start + window_size);
			}

			ip->add_constraint(constraint >= 0);
			constraints_added++;
		}
	}
	return constraints_added;
}

int main_program(int num_args, char* args[])
{
	using namespace std;

	attest(num_args == 3);
	NurseProblem problem(args[1], args[2]);

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

	// x[n][d][s] = 1 ⇔ Nurse n is working day d, shift s.
	auto x = ip.add_boolean_cube(problem.get_num_nurses(), problem.get_num_days(), problem.get_num_shifts());

	// Returns a variable (sum) that is 1 iff the nurse is working
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

	// Min consequtive days working.
	if (problem.get_min_consequtive() > 1) {
		int constraints_added = 0;
		for (int n = 0; n < problem.get_num_nurses(); ++n) {

			vector<Sum> working_on_days;
			for (int d = 0; d < problem.get_num_days(); ++d) {
				working_on_days.emplace_back(working_on_day(n, d));
			}
			constraints_added +=
				add_min_consequtive_constraints(&ip, problem.get_min_consequtive(), working_on_days);
		}
		clog << "Added " << constraints_added << " constraints for minimum consequtive days." << endl;
	}

	// Min consequtive days working a specific shift.
	{
		int constraints_added = 0;
		for (int s = 0; s < problem.get_num_shifts(); ++s) {
			if (problem.get_min_consequtive(s) > 1) {

				for (int n = 0; n < problem.get_num_nurses(); ++n) {

					vector<Sum> working_on_days;
					for (int d = 0; d < problem.get_num_days(); ++d) {
						working_on_days.emplace_back(x[n][d][s]);
					}
					constraints_added +=
						add_min_consequtive_constraints(&ip, problem.get_min_consequtive(s), working_on_days);
				}
			}
		}
		if (constraints_added > 0) {
			clog << "Added " << constraints_added << " constraints for minimum consequtive days of specific shifts." << endl;
		}
	}

	// Max consequtive days working.
	if (problem.get_max_consequtive() < problem.get_max_assignments()) {
		int constraints_added = 0;
		for (int n = 0; n < problem.get_num_nurses(); ++n) {
			for (int d = 0; d < problem.get_num_days() - problem.get_max_consequtive(); ++d) {
				Sum working_in_window = 0;
				for (int d2 = d; d2 < d + problem.get_max_consequtive(); ++d2) {
					working_in_window += working_on_day(n, d2);
				}
				ip.add_constraint(working_in_window <= problem.get_max_consequtive());
				constraints_added++;
			}
		}
		clog << "Added " << constraints_added << " constraints for maximum consequtive days." << endl;
	}

	// Max consequtive days working a specific shift.
	{
		int constraints_added = 0;
		for (int s = 0; s < problem.get_num_shifts(); ++s) {
			if (problem.get_max_consequtive(s) < problem.get_max_assignments()) {
				for (int n = 0; n < problem.get_num_nurses(); ++n) {
					for (int d = 0; d < problem.get_num_days() - problem.get_max_consequtive(s); ++d) {
						Sum working_in_window = 0;
						for (int d2 = d; d2 < d + problem.get_max_consequtive(s); ++d2) {
							working_in_window += x[n][d2][s];
						}
						ip.add_constraint(working_in_window <= problem.get_max_consequtive(s));
						constraints_added++;
					}
				}
			}
		}
		if (constraints_added > 0) {
			clog << "Added " << constraints_added << " constraints for maximum consequtive days of specific shifts." << endl;
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
	clog << "IP created." << endl;

	// Solve in completely silent mode to keep stdout clean.
	if (!ip.solve(nullptr, true)) {
		cout << problem_number << "\t"
		     << "infeasible" << "\t"
		     << 0 << endl;
	}

	double elapsed_time = omp_get_wtime() - start_time;

	clog << endl << endl << "IP solved." << endl;
	//for (int n = 0; n < problem.get_num_nurses(); ++n) {
	//	clog << "Nurse " << n << endl;
	//	for (int d = 0; d < problem.get_num_days(); ++d) {
	//		for (int s = 0; s < problem.get_num_shifts() - 1; ++s) {
	//			if (x[n][d][s].value() > 0.5) {
	//				clog << "-- Works day " << d << ", shift " << s << endl;
	//			}
	//		}
	//	}
	//}

	clog << endl << "Preferences are " << preferences.value() << endl;

	cout << problem_number << "\t"
	     << preferences.value() << "\t"
	     << elapsed_time << endl;

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
