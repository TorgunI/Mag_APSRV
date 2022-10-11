/*  ------------------------------------------------------
     Routine name:  my_diagram
     Description:   
     Project file:  ciclo1.prt

------------------------------------------------------  */

/*       Local stack variables                */
int i;
int j;
int c;
int itmp1;
int itmp2;
int itmp3;
double tmp1;
double tmp2;
double tmp3;
double tmp4;
double tmp5;
double tmp6;
double tmp7;
char f;
char tmp_f_1;
char u_s;
char u_r;
ret = 0;

switch (action){
case f_Stop:{
/* Index=1
   UID=1
   GeneratorClassName=TLanguage
   Name=LangBlock22
   Type=Язык программирования */

my_diagramv1_sfinal__0:
;
};break;
case f_GetDeri:{
};break;
case f_GetAlgFun:{
};break;
default:{

/* Index=0
   UID=0
   GeneratorClassName=TInputPin
   Name=Const_source18
   Type=Входной контакт */

state_vars->my_diagramv0_out_0_ = ((input_0));

/* Index=1
   UID=1
   GeneratorClassName=TLanguage
   Name=LangBlock22
   Type=Язык программирования */

symbol = (enum WeldingInputSymbol)state_vars->my_diagramv0_out_0_;
welding_engine(&engine, symbol);

output1 = (enum GazConsume)engine.gaz;
state_vars->my_diagramv1_out_0_ = output1;

output2 = (enum WeldingCurrent)engine.current;
state_vars->my_diagramv1_out_1_ = output1;

output3 = (enum ArcVoltage)engine.arcv;
state_vars->my_diagramv1_out_2_ = output1;

output4 = (enum WeldingSpeed)engine.speed;
state_vars->my_diagramv1_out_3_ = output1;

output5 = (enum OscVoltage)engine.oscv;
state_vars->my_diagramv1_out_4_ = output1;
/* Index=1
   UID=1
   GeneratorClassName=TLanguage
   Name=LangBlock22
   Type=Язык программирования */

my_diagramv1___0:
;

/* Index=2
   UID=2
   GeneratorClassName=TOutPin
   Name=OutPin15
   Type=Выходной контакт */

if(isfinite(state_vars->my_diagramv1_out_0_)){
out_0 = state_vars->my_diagramv1_out_0_;
};

/* Index=3
   UID=3
   GeneratorClassName=TOutPin
   Name=OutPin16
   Type=Выходной контакт */

if(isfinite(state_vars->my_diagramv1_out_1_)){
out_1 = state_vars->my_diagramv1_out_1_;
};

/* Index=4
   UID=4
   GeneratorClassName=TOutPin
   Name=OutPin17
   Type=Выходной контакт */

if(isfinite(state_vars->my_diagramv1_out_2_)){
out_2 = state_vars->my_diagramv1_out_2_;
};

/* Index=5
   UID=5
   GeneratorClassName=TOutPin
   Name=OutPin18
   Type=Выходной контакт */

if(isfinite(state_vars->my_diagramv1_out_3_)){
out_3 = state_vars->my_diagramv1_out_3_;
};

/* Index=6
   UID=6
   GeneratorClassName=TOutPin
   Name=OutPin19
   Type=Выходной контакт */

if(isfinite(state_vars->my_diagramv1_out_4_)){
out_4 = state_vars->my_diagramv1_out_4_;
};
};break;
};
