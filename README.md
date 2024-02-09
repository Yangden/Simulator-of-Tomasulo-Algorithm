# test cases for the simulator
Sample test cases are contained in the file "trace.txt". The format of a test case is :
```
<operation> <destination floating point register> <source floating point register 1> <source floating point register 2/immediate value>
```
For example, for the Add operation:
```
Add F6 F8 F2
```
The simulator supports four types of basic arithmetic operations on floating point numbers : 
1. ADD
2. SUB
3. MULT
4. DIV

Operations concerning memory: 
1. LOAD
2. STORE


# the configuration file 
The configuration file 'config.txt' is used to configure the layout of the **Reservation Stations**. 
The format of the file is written as : \
< Reservation Station size for LOAD > \
< Reservation Station size for STORE > \
< Reservation Station size for ADD > \
< Reservation Station size for MULT > \
< floating point register file size >

# output of the simulator 
The output will be writtien to 'trace.out.txt', which includes the following aspects : \
1. Insruction Status : issue cycle, complete cycle and cycle when broadcasting the results of arithmetic operations
2. Register Status : The statuses of reigsters after each cycle of executions.

# Hazards handling
This simulator focuses on data harzards including **RAW**, **WAW**, **WAR** and **WAR**. Memory Hazards are taken into account. 
