FUNCTION_BLOCK constgen_fb (*constgen_fb function block*)
VAR_INPUT
timestep : LREAL := 0.001;
timesec : LREAL := 0;
init : BOOL;
stop : BOOL;
END_VAR
VAR_OUTPUT
reg_close :LREAL := 1;
reg_open :LREAL := 1;
reg_v :LREAL := 4;
reg_level :LREAL := 2;
END_VAR
VAR
constgenv6_out_0 :LREAL := 2; (*Language out*)
constgenv6_out_1 :LREAL := 4; (*Language out*)
constgenv6vsection : ARRAY [0..3] OF LREAL := [4,0,0,0];
constgenv6vtotal :LREAL := 17;
constgenv6vsum : ARRAY [0..3] OF LREAL := [17,0,0,0];
constgenv6g :LREAL := 9.81;
constgenv6i :LREAL := 0;
constgenv6dh :LREAL := 0;
constgenv6outrate :LREAL := 0;
constgenv6d_rate :LREAL := 1;
constgenv6constgenv6_sinit___const_1 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6constgenv6_sinit___const_2 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6constgenv6_sinit___const_3 : ARRAY [0..4] OF LREAL := [1,3,3,1,1];
constgenv6constgenv6_sinit___const_4 : ARRAY [0..4] OF LREAL := [1,3,3,1,1];
constgenv6func1volume_onh :LREAL := 17;
constgenv6func1h_m :LREAL := 2;
constgenv6func1l :LREAL := 2;
constgenv6func1v0 :LREAL := 17;
constgenv6func2section_onh :LREAL := 2;
constgenv6func2h_level :LREAL := 2;
constgenv6func2constgenv6_sinit_func2__const_5 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6func1constgenv6_sinit_func1__const_6 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6func1constgenv6_sinit_func1__const_7 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6func1constgenv6_sinit_func1__const_8 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6func3hight_onv :LREAL := 2;
constgenv6func3volume :LREAL := 17;
constgenv6func3l :LREAL := 2;
constgenv6func3v0 :LREAL := 17;
constgenv6func3h0 :LREAL := 2;
constgenv6func4section_onv :LREAL := 2;
constgenv6func4volume :LREAL := 17;
constgenv6func4vtemp :LREAL := 13;
constgenv6func3constgenv6__func3__const_9 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv6func3constgenv6__func3__const_10 : ARRAY [0..4] OF LREAL := [0,2,4,5,6];
constgenv4_out_0 :LREAL; (*Multiplication output*)
constgenv5_out_0 :LREAL; (*Multiplication output*)
constgenv6v :LREAL := 17;
constgenv6v_deri :LREAL;
END_VAR
VAR CONSTANT
constgenv0_a :LREAL := 1; (*Constant value*)
constgenv1_a :LREAL := 1; (*Constant value*)
constgenv2_a :LREAL := 1; (*Constant value*)
constgenv3_a :LREAL := 1; (*Constant value*)
END_VAR
END_FUNCTION_BLOCK
