
/*  ------------------------------------------------------
     Routine name:  my_diagram
     Description:   
     Project file:  code.prt

------------------------------------------------------  */

/*  --- Base generator data types --- */
/* Real data type */
typedef double g_real_type;
/* Integer data type */
typedef int g_int_type;
/* Boolean data type */
typedef char g_boolean_type;
/* Complex data type */
typedef complex_64 g_complex_type;


/* Default initialization values */
const double input_0_default=0;
const double out_0_default=0;
const double my_diagramv0_out_0_default=0;
const double my_diagramv1_out_0_default=0;
const double my_diagramv2_out_0_default=0;
const double my_diagramv3_out_0_default=0;
const double my_diagramv4_out_0_default=0;
const double my_diagramv4_out_1_default=0;
const double my_diagramv4_out_2_default=0;

/* Project signal database hash */
const unsigned int sp_database_hash_32=0;
/* Project sheme structure hash */
const unsigned int sp_sheme_hash_32=2293132410;


/* Structures defines */
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
#define INTEGRATION_END_TIME 10
/* Model maximum iteration count */
#define INTEGRATION_MAX_LOOP_ITER_COUNT 10

/* External variables */
#define ext_vars_count 2

#define input_0 (*(double*)(ext_vars_addr[0]))
#define out_0 (*(double*)(ext_vars_addr[1]))

const ext_var_info_record ext_vars_names[ext_vars_count] = {
  {"input:0", "", (void*)&input_0_default, vt_double,   {1}, 0, dir_input, sizeof(double)} ,
  {"out:0", "", (void*)&out_0_default, vt_double,   {1}, 1, dir_out, sizeof(double)} 
};
/* State variables */
#define state_vars_count 7


typedef struct {
  double my_diagramv0_out_0_;
  double my_diagramv1_out_0_;
  double my_diagramv2_out_0_;
  double my_diagramv3_out_0_;
  double my_diagramv4_out_0_;
  double my_diagramv4_out_1_;
  double my_diagramv4_out_2_;
} t_state_vars;



const ext_var_info_record state_vars_names[state_vars_count] = {
  {"my_diagramv0_out_0", "Input pin state variable", (void*)&my_diagramv0_out_0_default, vt_double,   {1}, 0, dir_input, sizeof(double)} ,
  {"my_diagramv1_out_0", "Input pin state variable", (void*)&my_diagramv1_out_0_default, vt_double,   {1}, 8, dir_input, sizeof(double)} ,
  {"my_diagramv2_out_0", "Input pin state variable", (void*)&my_diagramv2_out_0_default, vt_double,   {1}, 16, dir_input, sizeof(double)} ,
  {"my_diagramv3_out_0", "Input pin state variable", (void*)&my_diagramv3_out_0_default, vt_double,   {1}, 24, dir_input, sizeof(double)} ,
  {"my_diagramv4_out_0", "Language out", (void*)&my_diagramv4_out_0_default, vt_double,   {1}, 32, dir_input, sizeof(double)} ,
  {"my_diagramv4_out_1", "Language out", (void*)&my_diagramv4_out_1_default, vt_double,   {1}, 40, dir_input, sizeof(double)} ,
  {"my_diagramv4_out_2", "Language out", (void*)&my_diagramv4_out_2_default, vt_double,   {1}, 48, dir_input, sizeof(double)} 
};
typedef char t_consts;
typedef char t_local;
#include "C:\stuff\Hleb\9\pneumo_ctrl.h"
#include "C:\stuff\Hleb\9\pneumo_ctrl.c"
struct PneumoEngine engine;
