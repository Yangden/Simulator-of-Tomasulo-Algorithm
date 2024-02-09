#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <limits>
#include <string>

using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::vector;

string inputtracename = "trace.txt";
// remove the ".txt" and add ".out.txt" to the end as output name
string outputtracename = inputtracename.substr(0, inputtracename.length() - 4) + ".out.txt";
string hardwareconfigname = "config.txt";

enum Operation
{
	ADD,
	SUB,
	MULT,
	DIV,
	LOAD,
	STORE
};
// The execute cycle of each operation: ADD, SUB, MULT, DIV, LOAD, STORE
const int OperationCycle[6] = {2, 2, 10, 40, 2, 2};

struct HardwareConfig
{
	int LoadRSsize;	 // number of load reservation stations
	int StoreRSsize; // number of store reservation stations
	int AddRSsize;	 // number of add reservation stations
	int MultRSsize;	 // number of multiply reservation stations
	int FRegSize;	 // number of fp registers
};

// We use the following structure to record the time of each instruction
struct InstructionStatus
{
	int cycleIssued;
	int cycleExecuted; // execution completed
	int cycleWriteResult;
        InstructionStatus() : cycleIssued(0), cycleExecuted(0), cycleWriteResult(0) {}
};

// Register Result Status structure
struct RegisterResultStatus
{
	string ReservationStationName;
	bool dataReady;
    // initialize as "all not ready"
	RegisterResultStatus() : dataReady(0) {}
};

/*********************************** ↓↓↓ Todo: Implement by you ↓↓↓ ******************************************/
struct Instruction
{
	string op, dest;
	string source1, source2; 

	Instruction() {}
	Instruction(string operation, string destination, 
	            string s1, string s2 )
	{
		op = operation;
		dest = destination;
		source1 = s1;
		source2 = s2;
	}
};

class RegisterResultStatuses
{
public:
    //constructor 
	RegisterResultStatuses () {}

	RegisterResultStatuses(HardwareConfig config) {
        _registers.resize(config.FRegSize, RegisterResultStatus());
    }
    

	// write to the status
	void write(Instruction instr, string RSname) {
       int pos = stoi(instr.dest.substr(1)); // get the position of reg
	   _registers.at(pos).ReservationStationName = RSname;
	   _registers.at(pos).dataReady = 0;
	}
    
	// issue stage : check if the data is ready
	bool if_ready(string r) {
        int pos = stoi(r.substr(1));
		bool ready = _registers.at(pos).dataReady;
		return ready;
	}

	string get_RSname(string r) {
		int pos = stoi(r.substr(1));
		return _registers.at(pos).ReservationStationName;
	}

	// broadcast stage : switch the status to ready
	void switch_status(string RSname) {
        for (RegisterResultStatus& reg : _registers) {
            if (reg.ReservationStationName == RSname && reg.dataReady == 0) {
				reg.dataReady = 1;
				break;
			}
		}
	}

	 

/*********************************** ↑↑↑ Todo: Implement by you ↑↑↑ ******************************************/
	/*
	Print all register result status. It is called by PrintRegisterResultStatus4Grade() for grading.
	If you don't want to write such a class, then make sure you call the same function and print the register
	result status in the same format.
	*/
	string _printRegisterResultStatus() const
	{
		std::ostringstream result;
		for (int idx = 0; idx < _registers.size(); idx++)
		{
			result << "F" + std::to_string(idx) << ": ";
			result << _registers[idx].ReservationStationName << ", ";
			result << "dataRdy: " << (_registers[idx].dataReady ? "Y" : "N") << ", ";
			result << "\n";
		}
		return result.str();
	}
/*********************************** ↓↓↓ Todo: Implement by you ↓↓↓ ******************************************/
private:
	vector<RegisterResultStatus> _registers;
	
};

// Define your Reservation Station structure
struct ReservationStation
{   
	// the name of the RS 
	string name;
	// when release : switch val of buzy to 0
	bool buzy;
	// pointer to other RS
	string Qj,Qk;
	// mark if source oprands are ready
	bool Vj, Vk;

	// index of the instruction
	int index;

    /* record total cycles for the instruction to complete.
	   initialize as op cycle of the instruction  */ 
	int countdown; 
    
	/*record remaining cycles for instructions to be completed
	 Initialize as -1, if source operand ready, be assigned with
	value of countdown. Decreasing to zero: mark complete.
	when release the reservation station, reinitialize to -1*/
	int timer; 

    ReservationStation() : buzy(0), Vj(0), Vk(0), timer(-1) {}
};

class CommonDataBus
{
public:
	// record the RS that broadcast
	CommonDataBus() { ready_RS = "";}
    
	// check if data available in CDB 
	bool check_available () {
        if (!ready_RS.empty()) {
            return true;
		} else {
			return false;
		}
	}

	// get the data in CDB
	string get_data() {
		return ready_RS;
	}

	// write to CDB
	void write(string RSname) {
		ready_RS = RSname;
	}

	// After forwarding data, clearing CDB
	void clear() {
		ready_RS.clear();
	}

private:
    string ready_RS;

};

/* */
class ReservationStations
{
public:
	// constructor 
	ReservationStations() {}
	ReservationStations(HardwareConfig config) {
        int totalsize = config.AddRSsize + config.LoadRSsize \
		+ config.StoreRSsize + config.MultRSsize;

		load_start = 0;
		store_start = load_start + config.LoadRSsize;
		add_start = store_start + config.StoreRSsize;
		mult_start = add_start + config.AddRSsize;
        
		_stations.resize(totalsize, ReservationStation());

        // tag names of each RS
		tag_name(load_start, store_start, "Load");
		tag_name(store_start, add_start, "Store");
		tag_name(add_start, mult_start, "Add");
		tag_name(mult_start, totalsize, "Mult");
    
	}
    
	// subroutine of 'occupy' 
	bool update_reg(Instruction instr, int start, int end, RegisterResultStatuses& reg,
	              int count, int issue_index, string current_result) 
	{
        string RSname; 
		bool if_issue = 0; // determine whether to issue the instruction 
                
		for (int i = start; i < end; i++) {
			if (_stations.at(i).buzy == 0 && current_result != _stations.at(i).name) {  // locate the available slot 
			    if_issue = 1;
			    _stations.at(i).buzy = 1;
                _stations.at(i).index = issue_index; // record order of issuing 
				_stations.at(i).countdown = count;  // record the execution cycle
                _stations.at(i).timer = -1;
				// special case : store instruction
				if (instr.op == "STORE") {
					_stations.at(i).Vk = 1;
					if (reg.if_ready(instr.dest)) {
						_stations.at(i).Vj = 1;
						_stations.at(i).timer = _stations.at(i).countdown;
					} else {
						_stations.at(i).Vj = 0;
						_stations.at(i).Qj = reg.get_RSname(instr.dest);
					}
                        break;
				}

				// write RS name
			   // RSname = _stations.at(i).name;
				
				// update register 
               // reg.write(instr,RSname);

				// special case : load instruction 
				if (instr.op == "LOAD") {
					_stations.at(i).Vj = _stations.at(i).Vk = 1;
                                        RSname = _stations.at(i).name;
                                        reg.write(instr, RSname);
					_stations.at(i).timer = _stations.at(i).countdown;
					break;
				}

				// if source oprand ready 
                if (reg.if_ready(instr.source1)) {
					_stations.at(i).Vj = 1;
				} else {  // not ready
					_stations.at(i).Qj = reg.get_RSname(instr.source1);
					_stations.at(i).Vj = 0;
				}
				if (reg.if_ready(instr.source2)) { // for operand 'k' 
					_stations.at(i).Vk = 1;
				} else {
					_stations.at(i).Qk = reg.get_RSname(instr.source2);
					_stations.at(i).Vk = 0; 
				}
				if (_stations.at(i).Vj && _stations.at(i).Vk) {
					_stations.at(i).timer = _stations.at(i).countdown;
				}
                                RSname = _stations.at(i).name;
                                reg.write(instr,RSname);

				break;
			}
		}
		return if_issue;
	}



    // issue stage : occupy RS resource 
	bool occupy(Instruction instr, RegisterResultStatuses& reg, int issue_index, 
	            string current_result) {
		int start, end, count;
		Operation opcode;
		bool if_issue;

		if (instr.op == "LOAD") {
            start = load_start; 
			end = store_start;
			opcode = LOAD;
            count = OperationCycle[opcode];
		}
		if (instr.op == "STORE") {
			start = store_start;
			end = add_start;
			opcode = STORE;
			count = OperationCycle[opcode];
		}
		if (instr.op == "ADD" || instr.op == "SUB") {
			start = add_start;
			end = mult_start;
			opcode = ADD;
			count = OperationCycle[opcode];
		}
		if (instr.op == "MULT" || instr.op == "DIV") {
			start = mult_start;
			end = _stations.size();
			if (instr.op == "MULT") {
				opcode = MULT;
			} else {
				opcode = DIV;
			}
			count = OperationCycle[opcode];
		}
        
		if_issue = update_reg(instr, start, end, reg, count, issue_index, current_result);
        return if_issue;
	}	
    
	// update status based on CDB 
	void update_status(string RSname) {
		for (ReservationStation &r : _stations) {
			if (r.buzy) {
				if (!r.Vj) {
					if (r.Qj == RSname) {
						r.Vj = 1;
					}
				}

				if (!r.Vk) {
					if(r.Qk == RSname) {
						r.Vk = 1;
					}
				}
			}
		}
	}
    
	/* After write : release the RS resources
	   reinitialize timer as -1
	   record the write cycle 
	*/
    void release(int pos, int cycle, vector<InstructionStatus>& track_instr) {
        track_instr[_stations.at(pos).index].cycleWriteResult = cycle;
		_stations.at(pos).timer = -1;
		_stations.at(pos).Vj = _stations.at(pos).Vk = 0;
		_stations.at(pos).buzy = 0;
	}

	// check if operands ready for execution
	bool check_ready(ReservationStation& r) {
		if (r.Vj == 1 && r.Vk == 1 && r.buzy == 1) {
			return true;
		} else {
			return false;
		}
	}
    
	/* decrement timer if it is ready. if timer == 0: mark complete
	for all RS with timer == 0: select the brocasted one. if none to broacast,
	return -1. otherwise, the position of broadcasted RS will return
	*/
	int execution_track(CommonDataBus& CDB, int cycle, vector<InstructionStatus>& track_instr) {
        string broadcast_result = "";
        int min_index = cycle; // initialize the min index
        int pos = -1;
		int i;

        // execution        
		for ( i = 0; i < _stations.size(); i++) {
			if (check_ready(_stations.at(i)) && _stations.at(i).timer > 0) { 
				_stations.at(i).timer --;
            }

			if (check_ready(_stations.at(i)) && _stations.at(i).timer < 0) {
				_stations.at(i).timer = _stations.at(i).countdown;
			} 

			if (_stations.at(i).timer == 0) {
				if (track_instr[_stations.at(i).index].cycleExecuted == 0) {
                    track_instr[_stations.at(i).index].cycleExecuted = cycle;
				}
				
			}

		} 
        
		// check what to broadcast
		for (i = 0; i < _stations.size(); i++) {
			if (check_ready(_stations.at(i)) && _stations.at(i).timer == 0) {
				if (min_index > _stations.at(i).index) {
                    min_index = _stations.at(i).index;
				    pos = i;
				}
			}
		}

                
		// broadcast to CDB 
		if (pos > -1) {
			broadcast_result = _stations.at(pos).name;
            CDB.write(broadcast_result);
		}
              
	    return pos;
	}


private:
	vector<ReservationStation> _stations;
	int load_start, store_start, add_start, mult_start;
	/* helper function :
	  tag the name for each RS*/
	void tag_name(int start, int end, string n) {
		for (int i = start; i < end; i++) {
			_stations.at(i).name = n + to_string(i - start);
		}
	}
};

// declaration of function 
void PrintRegisterResultStatus4Grade(const string &filename,
									 const RegisterResultStatuses &registerResultStatus,
									 const int thiscycle);

void simulateTomasulo(RegisterResultStatuses& Reg, ReservationStations& RS, 
                      vector<Instruction>& instructions, vector<InstructionStatus>& track_instr,
					  CommonDataBus& CDB);

void PrintResult4Grade(const string &filename, const vector<InstructionStatus> &instructionStatus);

// determine if all instructions finished
bool if_finished(vector<InstructionStatus>& track_instr) {
    bool flag = 1;
	for (InstructionStatus i : track_instr) {
		if (i.cycleWriteResult == 0) {
			flag = 0;
		}
	}

	return flag;
}

// Function to simulate the Tomasulo algorithm
void simulateTomasulo(RegisterResultStatuses& Reg, ReservationStations& RS, 
                      vector<Instruction>& instructions, vector<InstructionStatus>& track_instr,
					  CommonDataBus& CDB)

{

	int thiscycle = 1; // start cycle: 1
	int issue_index = 0; // record the order of issue instruction 
	Instruction instruction;
	bool if_issue = 0; // Determine whether issue a new instruction -> increase issue_index 
	string ready_RS;
	int pos = -1; 
	string current_result = "";

	while (thiscycle < 100000000)
	{    
        // section 1 : data bus forwards data
        current_result = CDB.get_data();
		if (pos != -1) { // There is data in the data bus : forward
            ready_RS = CDB.get_data(); // forward data
		    RS.update_status(ready_RS);
			Reg.switch_status(ready_RS);
			RS.release(pos, thiscycle, track_instr);
			CDB.clear();            // After forwarding, clear the data buse
		} 

		// section 2 : Reservation Stations should be updated every cycle, and broadcast to Common Data Bus
		pos = RS.execution_track(CDB, thiscycle, track_instr); 

		// Issue new instruction in each cycle
		if (issue_index < instructions.size()) {
            instruction = instructions.at(issue_index);
		    // count down cycle : solve "determine whether issue "
		    if_issue = RS.occupy(instruction, Reg, issue_index, current_result);
		    if (if_issue) {
			    track_instr[issue_index].cycleIssued = thiscycle;
			    issue_index++;
		    }
		}
		

		// At the end of this cycle, we need this function to print all registers status for grading
		PrintRegisterResultStatus4Grade(outputtracename, Reg, thiscycle);

		++thiscycle;

		// The simulator should stop when all instructions are finished.
		if (if_finished(track_instr)) {
			break;
		}
	}
};

/*********************************** ↑↑↑ Todo: Implement by you ↑↑↑ ******************************************/

/*
print the instruction status, the reservation stations and the register result status
@param filename: output file name
@param instructionStatus: instruction status
*/
void PrintResult4Grade(const string &filename, const vector<InstructionStatus> &instructionStatus)
{
	std::ofstream outfile(filename, std::ios_base::app); // append result to the end of file
	outfile << "Instruction Status:\n";
	for (int idx = 0; idx < instructionStatus.size(); idx++)
	{
		outfile << "Instr" << idx << ": ";
		outfile << "Issued: " << instructionStatus[idx].cycleIssued << ", ";
		outfile << "Completed: " << instructionStatus[idx].cycleExecuted << ", ";
		outfile << "Write Result: " << instructionStatus[idx].cycleWriteResult << ", ";
		outfile << "\n";
	}
	outfile.close();
}

/*
print the register result status each 5 cycles
@param filename: output file name
@param registerResultStatus: register result status
@param thiscycle: current cycle
*/
void PrintRegisterResultStatus4Grade(const string &filename,
									 const RegisterResultStatuses &registerResultStatus,
									 const int thiscycle)
{
	if (thiscycle % 5 != 0)
		return;
	std::ofstream outfile(filename, std::ios_base::app); // append result to the end of file
	outfile << "Cycle " << thiscycle << ":\n";
	outfile << registerResultStatus._printRegisterResultStatus() << "\n";
	outfile.close();
}

int main(int argc, char **argv)
{
	if (argc > 1)
	{
		hardwareconfigname = argv[1];
		inputtracename = argv[2];
	}

	HardwareConfig hardwareConfig;
	std::ifstream config;
	config.open(hardwareconfigname);
	config >> hardwareConfig.LoadRSsize;  // number of load reservation stations
	config >> hardwareConfig.StoreRSsize; // number of store reservation stations
	config >> hardwareConfig.AddRSsize;	  // number of add reservation stations
	config >> hardwareConfig.MultRSsize;  // number of multiply reservation stations
	config >> hardwareConfig.FRegSize;	  // number of fp registers
	config.close();

/*********************************** ↓↓↓ Todo: Implement by you ↓↓↓ ******************************************/

	// Read instructions from a file (replace 'instructions.txt' with your file name)
	ifstream traces;
	string line, d, s1, s2, op;
	traces.open(inputtracename); // open the file
	vector<Instruction> instructions; // store all instructions

	while (getline(traces, line)) { // read line by line
        istringstream iss(line);
        iss >> op >> d >> s1 >> s2;
		Instruction instr(op, d, s1, s2);
        instructions.push_back(instr);
        if (traces.eof()) {
            break;

	}
       }
	 
    // Initialize the register result status and RS
	RegisterResultStatuses Reg(hardwareConfig);
    ReservationStations RS(hardwareConfig);
      
	// Initialize common data bus 
	CommonDataBus CDB = CommonDataBus();

	// Initialize the instruction status table
	vector<InstructionStatus> track_instr;
	track_instr.resize(instructions.size(), InstructionStatus());
	// ..
     
	// Simulate Tomasulo:
	// simulateTomasulo(registerResultStatus, instructionStatus, ...);
	simulateTomasulo(Reg, RS, instructions, track_instr, CDB);

/*********************************** ↑↑↑ Todo: Implement by you ↑↑↑ ******************************************/

	// At the end of the program, print Instruction Status Table for grading
	PrintResult4Grade(outputtracename, track_instr);

	return 0;
}
