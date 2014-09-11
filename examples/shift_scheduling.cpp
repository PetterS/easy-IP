// Petter Strandmark 2014.
//
// Reads and solves problem files from “Employee Shift Scheduling Benchmark
// Data Sets”
// http://www.cs.nott.ac.uk/~tec/NRP/
//
// First argument to the program is the problem .txt file.
// 
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#ifdef USE_OPENMP
	#include <omp.h>
#else
	namespace { double omp_get_wtime() { return 0; } }
#endif

#include <easy-ip.h>

#include "string_utils.h"
using namespace spii;

string read_line(istream& file)
{
	string line;
	getline(file, line);
	if (!file) {
		return "";
	}
	int i = 0;
	
	// Check for comment.
	for (auto ch: line) {
		if (ch == '#') {
			return read_line(file);
		}
		else if (ch != ' ') {
			break;
		}
	}

	return line;
}

class ShiftSchedulingProblem
{
public:
	ShiftSchedulingProblem(string txt_file_name);

	struct Shift
	{
		string ID;
		int duration;
		vector<int> forbidden_after;
	};

	struct Staff
	{
		string ID = "";
		vector<int> max_shifts;
		int max_total_minutes = 0;
		int min_total_minutes = 0;
		int max_consequtive_shifts   = 0;
		int min_consequtive_shifts   = 0;
		int min_consequtive_days_off = 0;
		int max_weekends = 0;
		vector<int> days_off;
	};

	struct ShiftRequest
	{
		int staff;
		int day;
		int shift;
		int weight;
	};

	struct CoverConstraint
	{
		int day;
		int shift;
		int requirement;
		int under_weight;
		int over_weight;
	};

	int get_num_staff()  const { return all_staff.size(); }
	int get_num_days()   const { return num_days; }
	int get_num_shifts() const { return shifts.size(); }

	const Staff& get_staff(int p) const { return all_staff.at(p); }
	const Shift& get_shift(int s) const { return shifts.at(s); }

	const vector<ShiftRequest>& get_shift_on_requests()  const { return shift_on_requests; }
	const vector<ShiftRequest>& get_shift_off_requests() const { return shift_off_requests; }


	const vector<CoverConstraint>& get_cover_constraints() const { return cover_constraints; }

	Sum create_ip(IP& ip, const vector<vector<vector<BooleanVariable>>>& x) const;

protected:
	vector<Shift> shifts;
	map<string, size_t> id_to_shift;

	vector<Staff> all_staff;
	map<string, size_t> id_to_staff;
	
	vector<ShiftRequest> shift_on_requests;
	vector<ShiftRequest> shift_off_requests;

	vector<CoverConstraint> cover_constraints;

	int num_days   = 0;
};

ShiftSchedulingProblem::ShiftSchedulingProblem(string txt_file_name)
{
	ifstream file(txt_file_name);

	attest(read_line(file) == "SECTION_HORIZON");
	num_days = from_string<int>(read_line(file));
	clog << "Number of days: " << num_days << endl;
	attest(read_line(file) == "");

	attest(read_line(file) == "SECTION_SHIFTS");
	auto shift_string = read_line(file);
	vector<pair<int, string>> forbidden_after_pairs;
	while (shift_string != "") {
		auto parts = split(shift_string, ',');
		Shift shift;

		int e = 0;
		shift.ID = parts.at(e++);
		shift.duration = from_string<int>(parts.at(e++));
		if (e < parts.size()) {
			forbidden_after_pairs.emplace_back(shifts.size(), parts.at(e++));
		}
		attest(e == parts.size());

		clog << "-- Shift " << shift.ID << ", "
			    << shift.duration << " min. "
			    << "Forbidden after: " << to_string(shift.forbidden_after) << endl;
		shifts.emplace_back(shift);
		id_to_shift[shift.ID] = shifts.size() - 1;

		shift_string = read_line(file);
	}
	clog << "Number of shifts: " << shifts.size() << endl;

	// We need to parse the forbidden combinations afterwards,
	// when we have seen all shift IDs.
	for (auto& ix_string: forbidden_after_pairs) {
		auto forbidden_after = split(ix_string.second, '|');
		for (auto& forbidden: forbidden_after) {
			auto itr = id_to_shift.find(forbidden);
			attest(itr != id_to_shift.end());
			shifts[ix_string.first].forbidden_after.emplace_back(itr->second);
		}
	}

	attest(read_line(file) == "SECTION_STAFF");
	auto staff_string = read_line(file);
	while (staff_string != "") {
		auto parts = split(staff_string, ',');
		Staff staff;
		int e = 0;
		staff.ID = parts.at(e++);

		auto max_shifts = split(parts.at(e++), '|');
		staff.max_shifts.resize(shifts.size(), -1);
		for (auto& max_shift: max_shifts) {
			auto id_value = split(max_shift, '=');
			staff.max_shifts.at(id_to_shift[id_value.at(0)]) = from_string<int>(id_value.at(1));
		}

		staff.max_total_minutes        = from_string<int>(parts.at(e++));
		staff.min_total_minutes        = from_string<int>(parts.at(e++));
		staff.max_consequtive_shifts   = from_string<int>(parts.at(e++));
		staff.min_consequtive_shifts   = from_string<int>(parts.at(e++));
		staff.min_consequtive_days_off = from_string<int>(parts.at(e++));
		staff.max_weekends             = from_string<int>(parts.at(e++));
		attest(e == parts.size());

		clog << "-- Staff " << staff.ID << ", " << to_string(staff.max_shifts) << " max shifts." << endl;
		all_staff.emplace_back(staff);
		id_to_staff[staff.ID] = all_staff.size() - 1;

		staff_string = read_line(file);
	}
	clog << "Number of staff: " << all_staff.size() << endl;

	attest(read_line(file) == "SECTION_DAYS_OFF");
	auto days_off_string = read_line(file);
	while (days_off_string != "") {
		auto days_off = split(days_off_string, ',');
		auto s = id_to_staff[days_off.at(0)];
		for (auto itr = days_off.begin() + 1; itr != days_off.end(); ++itr) {
			all_staff.at(s).days_off.emplace_back(from_string<int>(*itr));
		}
		days_off_string = read_line(file);
	}

	for (const auto& staff: all_staff) {
		clog << "-- Staff " << staff.ID << " has " << staff.days_off.size() << " day(s) off." << endl;
	}

	attest(read_line(file) == "SECTION_SHIFT_ON_REQUESTS");
	auto shift_on_request_string = read_line(file);
	while (shift_on_request_string != "") {
		auto shift_on_request = split(shift_on_request_string, ',');
		attest(shift_on_request.size() == 4);
		shift_on_requests.emplace_back();
		shift_on_requests.back().staff  = id_to_staff[shift_on_request.at(0)];
		shift_on_requests.back().day    = from_string<int>(shift_on_request.at(1));
		shift_on_requests.back().shift  = id_to_shift[shift_on_request.at(2)];
		shift_on_requests.back().weight = from_string<int>(shift_on_request.at(3));

		shift_on_request_string = read_line(file);
	}
	clog << "Read " << shift_on_requests.size() << " shift on requests." << endl;

	attest(read_line(file) == "SECTION_SHIFT_OFF_REQUESTS");
	auto shift_off_request_string = read_line(file);
	while (shift_off_request_string != "") {
		auto shift_off_request = split(shift_off_request_string, ',');
		attest(shift_off_request.size() == 4);
		shift_off_requests.emplace_back();
		shift_off_requests.back().staff  = id_to_staff[shift_off_request.at(0)];
		shift_off_requests.back().day    = from_string<int>(shift_off_request.at(1));
		shift_off_requests.back().shift  = id_to_shift[shift_off_request.at(2)];
		shift_off_requests.back().weight = from_string<int>(shift_off_request.at(3));

		shift_off_request_string = read_line(file);
	}
	clog << "Read " << shift_off_requests.size() << " shift off requests." << endl;

	attest(read_line(file) == "SECTION_COVER");
	auto cover_string = read_line(file);
	while (cover_string != "") {
		auto cover = split(cover_string, ',');
		cover_constraints.emplace_back();
		int e = 0;
		cover_constraints.back().day   = from_string<int>(cover.at(e++));
		cover_constraints.back().shift = id_to_shift[cover.at(e++)];
		cover_constraints.back().requirement  = from_string<int>(cover.at(e++));
		cover_constraints.back().under_weight = from_string<int>(cover.at(e++));
		cover_constraints.back().over_weight  = from_string<int>(cover.at(e++));
		attest(e == cover.size());

		cover_string = read_line(file);
	}
	clog << "Read " << cover_constraints.size() << " cover constraints." << endl;

	attest(id_to_staff.size() == all_staff.size());
	attest(id_to_shift.size() == shifts.size());
}

Sum ShiftSchedulingProblem::create_ip(IP& ip, const vector<vector<vector<BooleanVariable>>>& x) const
{
	// Returns a variable (Sum) that is 1 iff a person is working
	// on a particular day.
	auto working_on_day = [&](int person, int day)
	{
		Sum working_on_day = 0;
		for (int s = 0; s < get_num_shifts(); ++s) {
			working_on_day += x[person][day][s];
		}
		return working_on_day;
	};

	// At most one shift per day per person.
	for (int p = 0; p < get_num_staff(); ++p) {
		for (int d = 0; d < get_num_days(); ++d) {
			Sum num_day_shifts = 0;
			for (int s = 0; s < get_num_shifts(); ++s) {
				num_day_shifts += x[p][d][s];
			}
			ip.add_constraint(num_day_shifts <= 1);
		}
	}

	// Maximum shifts.
	for (int s = 0; s < get_num_shifts(); ++s) {
		for (int p = 0; p < get_num_staff(); ++p) {
			Sum working_shifts = 0;
			for (int d = 0; d < get_num_days(); ++d) {
				working_shifts += x[p][d][s];
			}
			ip.add_constraint(working_shifts <= get_staff(p).max_shifts.at(s));
		}
	}

	// Max/min total minutes.
	for (int p = 0; p < get_num_staff(); ++p) {
		Sum working_minutes = 0;
		for (int d = 0; d < get_num_days(); ++d) {
			for (int s = 0; s < get_num_shifts(); ++s) {
				working_minutes += get_shift(s).duration * x[p][d][s];
			}
		}
		ip.add_constraint(get_staff(p).min_total_minutes, working_minutes, get_staff(p).max_total_minutes);
	}

	// Shifts that cannot follow other shifts.
	{
		int constraints_added = 0;
		for (int s1 = 0; s1 < get_num_shifts(); ++s1) {
			for (int s2: get_shift(s1).forbidden_after) {
				for (int p = 0; p < get_num_staff(); ++p) {
					for (int d = 0; d < get_num_days() - 1; ++d) {
						ip.add_constraint(x[p][d][s1] + x[p][d+1][s2] <= 1);
						constraints_added++;
					}
				}
			}
		}
		clog << "Added " << constraints_added << " constraints for impossible shift combinations." << endl;
	}

	// Min consequtive days working.
	{
		int constraints_added = 0;
		for (int p = 0; p < get_num_staff(); ++p) {
			int min_consequtive = get_staff(p).min_consequtive_shifts;
			int max_consequtive = get_staff(p).max_consequtive_shifts;

			vector<Sum> working_on_days;
			for (int d = 0; d < get_num_days(); ++d) {
				working_on_days.emplace_back(working_on_day(p, d));
			}
			constraints_added +=
				ip.add_min_consequtive_constraints(min_consequtive, working_on_days, true);

			constraints_added +=
				ip.add_max_consequtive_constraints(max_consequtive, working_on_days);
		}
		clog << "Added " << constraints_added << " constraints for minimum/maximum consequtive shifts." << endl;
	}

	// Min consequtive days off.
	{
		int constraints_added = 0;
		for (int p = 0; p < get_num_staff(); ++p) {
			int min_consequtive = get_staff(p).min_consequtive_days_off;
			if (min_consequtive > 1) {
				vector<Sum> day_off;
				for (int d = 0; d < get_num_days(); ++d) {
					day_off.emplace_back(1 - working_on_day(p, d));
				}
				constraints_added +=
					ip.add_min_consequtive_constraints(min_consequtive, day_off, true);
			}
		}
		clog << "Added " << constraints_added << " constraints for minimum consequtive days off." << endl;
	}

	// Max working weekends.
	for (int p = 0; p < get_num_staff(); ++p) {
		Sum working_weekends = 0;
		for (int d = 5; d < get_num_days(); d += 7) {
			auto saturday = working_on_day(p, d);
			attest(d + 1 < get_num_days());
			auto sunday   = working_on_day(p, d + 1);

			auto weekend = ip.add_boolean();
			ip.add_constraint(saturday <= weekend);
			ip.add_constraint(sunday   <= weekend);
			working_weekends += weekend;
		}
		ip.add_constraint(working_weekends <= get_staff(p).max_weekends);
	}

	// Days off.
	for (int p = 0; p < get_num_staff(); ++p) {
		for (int day_off: get_staff(p).days_off) {
			ip.add_constraint(working_on_day(p, day_off) == 0);
		}
	}

	Sum objective = 0;

	// Shift on requests.
	for (auto request: get_shift_on_requests()) {
		objective += request.weight * (1 - x[request.staff][request.day][request.shift]);
	}

	// Shift off requests.
	for (auto request: get_shift_off_requests()) {
		objective += request.weight * x[request.staff][request.day][request.shift];
	}

	// Cover
	for (auto constraint: get_cover_constraints()) {
		Sum num_working = 0;
		for (int p = 0; p < get_num_staff(); ++p) {
			num_working += x[p][constraint.day][constraint.shift];
		}
		auto under_slack = ip.add_variable(IP::Real);
		objective += constraint.under_weight * under_slack;
		auto over_slack = ip.add_variable(IP::Real);
		objective += constraint.over_weight * over_slack;
		ip.add_constraint(over_slack >= 0);
		ip.add_constraint(under_slack >= 0);
		ip.add_constraint(num_working + under_slack - over_slack == constraint.requirement);
	}

	ip.add_objective(objective);

	return objective;
}

int main_program(int num_args, char* args[])
{
	using namespace std;

	attest(num_args >= 2);
	ShiftSchedulingProblem problem(args[1]);

	double start_time = omp_get_wtime();

	IP ip;

	// The assignment variables.
	//   x[p][d][s] = 1 ⇔ Person p is working day d, shift s.
	auto x = ip.add_boolean_cube(problem.get_num_staff(), problem.get_num_days(), problem.get_num_shifts());
	auto objective = problem.create_ip(ip, x);

	bool first_solution = true;
	auto print_solution = [&] () 
	{
		double elapsed_time = omp_get_wtime() - start_time;
		clog << "Integer solution with objective " << objective.value() << " in " << elapsed_time << " seconds." << endl;

		if (!first_solution) {
			cout << "======" << endl;
		}
		first_solution = false;

		for (int p = 0; p < problem.get_num_staff(); ++p) {
			for (int d = 0; d < problem.get_num_days(); ++d) {
				bool printed = false;
				for (int s = 0; s < problem.get_num_shifts(); ++s) {
					if (x[p][d][s].value()) {
						cout << problem.get_shift(s).ID;
						printed = true;
					}
				}
				cout << '\t';
			}
			cout << endl;
		}
	};

	//ip.save_MPS("scheduling");
	clog << "Starting search. Will print integer solutions to stdout." << endl;

	bool silent_mode = false;
	if (ip.solve(print_solution, silent_mode)) {
		double elapsed_time = omp_get_wtime() - start_time;
		clog << "Objective " << objective.value() << " proven optimal in " << elapsed_time << " seconds." << endl;
	}
	else {
		clog << "IP was not solved." << endl;
		return 2;
	}

	print_solution();

	return 0;
}

int main(int num_args, char* args[])
{
	try {
		return main_program(num_args, args);
	}
	catch (std::exception& err) {
		cerr << "Error: " << err.what() << endl;
		return 1;
	}
}
