/*  ------------------------------------------------------
     Routine name:  constgen
     Description:   
     Project file:  5.1.prt

------------------------------------------------------  */

#include <bur/plctypes.h>
#ifdef __cplusplus
	extern "C"
	{
#endif

	#include "constgen.h"

#ifdef __cplusplus
	};
#endif

/*  --- Source model preferences --- */
/* Minimum integration step */
#define INTEGRATION_MIN_STEP 0.001
/* Maximum integration step */
#define INTEGRATION_MAX_STEP 0.001
/* Integration synchronization step */
#define INTEGRATION_SYNC_STEP 0
/* Model integration method */
#define INTEGRATION_METHOD 0
/* Model relative error */
#define INTEGRATION_RELATIVE_ERROR 0.0001
/* Model absolute error */
#define INTEGRATION_ABSOLUTE_ERROR 1E-6
/* Model end time */
#define INTEGRATION_END_TIME 1000
/* Model maximum iteration count */
#define INTEGRATION_MAX_LOOP_ITER_COUNT 10
/* Real time synchronization flag */
#define MODEL_REAL_TIME_SYNC_FLAG 1
/* Real time synchronization gain */
#define MODEL_REAL_TIME_SYNC_GAIN 1


/* Summary dimension for integration method */
#define INTEGRATION_FULL_DIMENSION 1

/* Standart algorithms defines */
#include "c_types.h"
/*  --- Base generator data types --- */
/* Real data type */
typedef LREAL g_real_type;
/* Integer data type */
typedef DINT g_int_type;
/* Boolean data type */
typedef BOOL g_boolean_type;
/* Complex data type */
typedef complex_64 g_complex_type;

#include "integr_method.h"
/* Model 1 function block */
void constgen_fb(struct constgen_fb* inst)
{

/*       Local stack variables                */
DINT i;
DINT j;
DINT c;
DINT itmp1;
DINT itmp2;
DINT itmp3;
LREAL tmp1;
LREAL tmp2;
LREAL tmp3;
LREAL tmp4;
LREAL tmp5;
LREAL tmp6;
LREAL tmp7;
BOOL f;
BOOL tmp_f_1;
BOOL u_s;
BOOL u_r;
DINT ret;
LREAL step;
LREAL modeltime;
DINT action;
ret = 0;

step = inst->timestep;
modeltime = inst->timesec;
action = f_GoodStep;
if(inst->init){
action = f_InitState;
inst->constgenv6g=(9.81);
inst->constgenv6i=(1);
inst->constgenv6i=(1);
constgenv6_sinit__5:
if(inst->constgenv6i > ((5)-(1)))
goto constgenv6_sinit__10;
inst->constgenv6vsection[(int)inst->constgenv6i - 1]=(((inst->constgenv6constgenv6_sinit___const_1[(int)(inst->constgenv6i+(g_real_type)((1))) - 1]-inst->constgenv6constgenv6_sinit___const_2[(int)inst->constgenv6i - 1])*(inst->constgenv6constgenv6_sinit___const_3[(int)inst->constgenv6i - 1]+inst->constgenv6constgenv6_sinit___const_4[(int)(inst->constgenv6i+(g_real_type)((1))) - 1]))/(2));
inst->constgenv6vtotal=(inst->constgenv6vtotal+inst->constgenv6vsection[(int)inst->constgenv6i - 1]);
inst->constgenv6vsum[(int)inst->constgenv6i - 1]=inst->constgenv6vtotal;
inst->constgenv6i = inst->constgenv6i + 1;
goto constgenv6_sinit__5;
constgenv6_sinit__10:
inst->constgenv6func1h_m=(2);
inst->constgenv6func2h_level=inst->constgenv6func1h_m;
inst->constgenv6i=((5)-(1));
inst->constgenv6i=((5)-(1));
constgenv6_sinit_func2_2:
if(((-1) >= 0)&&(inst->constgenv6i > (1))) {
goto constgenv6_sinit_func2_6;
} else {
if(((-1) < 0)&&(inst->constgenv6i < (1)))
goto constgenv6_sinit_func2_6;
};
if((inst->constgenv6func2h_level<inst->constgenv6func2constgenv6_sinit_func2__const_5[(int)(inst->constgenv6i+(g_real_type)((1))) - 1])){goto constgenv6_sinit_func2_4;} else {goto constgenv6_sinit_func2_5;};
constgenv6_sinit_func2_4:
inst->constgenv6func2section_onh=inst->constgenv6i;
constgenv6_sinit_func2_5:
if((-1) != 0){
inst->constgenv6i = inst->constgenv6i + (-1);
goto constgenv6_sinit_func2_2;
};
constgenv6_sinit_func2_6:
;
inst->constgenv6func1l=inst->constgenv6func2section_onh;
if((inst->constgenv6func1l==(g_real_type)((1)))){goto constgenv6_sinit_func1_2;} else {goto constgenv6_sinit_func1_4;};
constgenv6_sinit_func1_2:
inst->constgenv6func1v0=(0);
goto constgenv6_sinit_func1_5;
constgenv6_sinit_func1_4:
inst->constgenv6func1v0=inst->constgenv6vsum[(int)(inst->constgenv6func1l-(g_real_type)((1))) - 1];
constgenv6_sinit_func1_5:
inst->constgenv6func1volume_onh=(inst->constgenv6func1v0+((inst->constgenv6vsection[(int)inst->constgenv6func1l - 1]*(inst->constgenv6func1h_m-inst->constgenv6func1constgenv6_sinit_func1__const_6[(int)inst->constgenv6func1l - 1]))/(inst->constgenv6func1constgenv6_sinit_func1__const_7[(int)(inst->constgenv6func1l+(g_real_type)((1))) - 1]-inst->constgenv6func1constgenv6_sinit_func1__const_8[(int)inst->constgenv6func1l - 1])));
constgenv6_sinit_func1_6:
;
inst->constgenv6_out_1=inst->constgenv6func1volume_onh;
constgenv6_sinit__13:
;
integr_method_init(&inst->constgenv6v, &inst->constgenv6v_deri, step, 1, 0, INTEGRATION_METHOD);
};
if(inst->stop){
action = f_Stop;
/* Index=6
   UID=6
   GeneratorClassName=TLanguage
   Name=Tank
   Type=Язык программирования */

constgenv6_sfinal__0:
;
} else {

/* Index=4
   UID=4
   GeneratorClassName=TMulDbl
   Name=Mul_oper8
   Type=Перемножитель */

inst->constgenv4_out_0 = inst->constgenv2_a*inst->constgenv0_a;

/* Index=5
   UID=5
   GeneratorClassName=TMulDbl
   Name=Mul_oper9
   Type=Перемножитель */

inst->constgenv5_out_0 = inst->constgenv3_a*inst->constgenv1_a;
/* Index=6
   UID=6
   GeneratorClassName=TLanguage
   Name=Tank
   Type=Язык программирования */

inst->constgenv6func3volume=inst->constgenv6_out_1;
inst->constgenv6func4volume=inst->constgenv6func3volume;
inst->constgenv6func4vtemp=inst->constgenv6vtotal;
inst->constgenv6i=((5)-(1));
inst->constgenv6i=((5)-(1));
constgenv6__func4_3:
if(((-1) >= 0)&&(inst->constgenv6i > (1))) {
goto constgenv6__func4_8;
} else {
if(((-1) < 0)&&(inst->constgenv6i < (1)))
goto constgenv6__func4_8;
};
if((inst->constgenv6func4volume<inst->constgenv6func4vtemp)){goto constgenv6__func4_5;} else {goto constgenv6__func4_6;};
constgenv6__func4_5:
inst->constgenv6func4section_onv=inst->constgenv6i;
constgenv6__func4_6:
inst->constgenv6func4vtemp=(inst->constgenv6func4vtemp-inst->constgenv6vsection[(int)inst->constgenv6i - 1]);
if((-1) != 0){
inst->constgenv6i = inst->constgenv6i + (-1);
goto constgenv6__func4_3;
};
constgenv6__func4_8:
;
inst->constgenv6func3l=inst->constgenv6func4section_onv;
if((inst->constgenv6func3l==(g_real_type)((1)))){goto constgenv6__func3_2;} else {goto constgenv6__func3_4;};
constgenv6__func3_2:
inst->constgenv6func3v0=(0);
goto constgenv6__func3_5;
constgenv6__func3_4:
inst->constgenv6func3v0=inst->constgenv6vsum[(int)(inst->constgenv6func3l-(g_real_type)((1))) - 1];
constgenv6__func3_5:
inst->constgenv6func3h0=inst->constgenv6func3constgenv6__func3__const_9[(int)inst->constgenv6func3l - 1];
inst->constgenv6func3hight_onv=(inst->constgenv6func3h0+(((inst->constgenv6func3volume-inst->constgenv6func3v0)*(inst->constgenv6func3constgenv6__func3__const_10[(int)(inst->constgenv6func3l+(g_real_type)((1))) - 1]-inst->constgenv6func3h0))/inst->constgenv6vsection[(int)inst->constgenv6func3l - 1]));
constgenv6__func3_7:
;
inst->constgenv6_out_0=inst->constgenv6func3hight_onv;
inst->constgenv6dh=(inst->constgenv6_out_0-(1));
if((inst->constgenv6dh<(g_real_type)((0)))){goto constgenv6___3;} else {goto constgenv6___4;};
constgenv6___3:
inst->constgenv6dh=(0);
constgenv6___4:
inst->constgenv6outrate=((sqrt((((g_real_type)((2))*inst->constgenv6g)*inst->constgenv6dh))*(0.05))*inst->constgenv5_out_0);
inst->constgenv6d_rate=(inst->constgenv4_out_0-inst->constgenv6outrate);
inst->constgenv6v_deri=inst->constgenv6d_rate;
constgenv6___7:
;

/* Index=8
   UID=8
   GeneratorClassName=TOutPin
   Name=OutPin15
   Type=Выходной контакт */

if(isfinite(inst->constgenv5_out_0)){
inst->reg_close = inst->constgenv5_out_0;
};

/* Index=9
   UID=9
   GeneratorClassName=TOutPin
   Name=OutPin16
   Type=Выходной контакт */

if(isfinite(inst->constgenv4_out_0)){
inst->reg_open = inst->constgenv4_out_0;
};

/* Index=10
   UID=10
   GeneratorClassName=TOutPin
   Name=OutPin17
   Type=Выходной контакт */

if(isfinite(inst->constgenv6_out_1)){
inst->reg_v = inst->constgenv6_out_1;
};

/* Index=11
   UID=11
   GeneratorClassName=TOutPin
   Name=OutPin18
   Type=Выходной контакт */

if(isfinite(inst->constgenv6_out_0)){
inst->reg_level = inst->constgenv6_out_0;
};
integr_method_call(&inst->constgenv6v, &inst->constgenv6v_deri, step, 1, 0, INTEGRATION_METHOD);
  };
inst->timesec = inst->timesec + inst->timestep;
};
