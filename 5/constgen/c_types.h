
  /* ������������� ����� ������ */
#define    vt_double         0      /* ������������ */
#define    vt_bool           1      /* ��������� */
#define    vt_int            2      /* ����� */
#define    vt_pointer        3      /* ��� ��� - ����� ����� ��������� �� ����������� ������� (��� ��������� ����������-���������� ��������) */

  /* ���������������
  ������������ */
#define    vt_float          4      /* 4 �������� � �� */
#define    vt_fix32          5      /* 4 �������� fixpoint(16.16) */
  /* ����� */
#define    vt_int64          6      /* 64 ������ ����� */
#define    vt_int16          7      /* 16 ������ ����� */

#define    vt_complex        8      /* ����������� �� 2 64 ������ ������������ */
#define    vt_complex_32     9      /* ����������� �� 2 32 ������ ������������ */
#define    vt_complex_fix16 10      /* ����������� �� 2 32 ������ ������������ � ������������� ������ */

  /* �������������� ��� - 8 ������ �������� ����� */
#define    vt_int8          11

  /* ���������� ����� ���� ������ */
#define    vt_uint8         12
#define    vt_uint16        13
#define    vt_uint32        14
#define    vt_uint64        15

  /* ����������� ���������� */
#define    dir_input         0
#define    dir_out           1
#define    dir_inout         2

  /* ����� ������ run-������� ����� */
#define    f_InitState       1      /* ������ ��������� ��������� */
#define    f_UpdateOuts      2      /* �������� ������ �� ��������������� ���� */
#define    f_GoodStep        3      /* �������� ������ �� "�������" ���� */
#define    f_GetDeri         4      /* ��������� �������� ������ ������ ���������������� ��������� */
#define    f_GetAlgFun       5      /* ��������� �������� ������ ������ �������������� ��������� */
#define    f_SetState        6      /* ��������� �������� ���������� ���������� ��������� (����� ���� ��������������) */
#define    f_UpdateProps     7      /* �������� ������ ���������� (� ������ ����� ������������) */
#define    f_GetJacobyState  8      /* ��������� �������� ���������� ���������� ��������� ��� ������� �������� */
#define    f_UpdateJacoby    9      /* �������� ������� ����� */
#define    f_RestoreOuts    10      /* �������� ������ ����� �������� (������ ���� ����� ����, �.�. ������ �� ����� ����� ������������) */
#define    f_SetAlgOut      11      /* ��������� ������ �����, ���������� �������������� ���������� */
#define    f_InitAlgState   12      /* ��������� ��������� ����������� ��� �������������� ���������� */
#define    f_Stop           13      /* ���������� ��� ��������� ������� (����� �������������) */
#define    f_InitObjects    14      /* ������������� ��������, �������� � �.�. (����� ����� ����������) (������ �������������) */

  /* ����� ������� �������������� ������� ����� */ 
#define    i_GetBlockType    1      /* �������� ��� ����� (��������, ������������ � �.� */
#define    i_GetDifCount     2      /* �������� ����� ���������������� ���������� */
#define    i_GetDisCount     3      /* ���� Result > 0 �� �������������� ���� f_SetState */
#define    i_GetAlgCount     4      /* �������� ����� �������������� ���������� */
#define    i_GetCount        5      /* �������� ����������� ������\������� */
#define    i_GetInit         6      /* �������� ���� ����������� ������� �� ������ */
#define    i_GetPropErr      7      /* �������� ������������ ������� ���������� ����� (����� �����������) */
#define    i_HaveSpetialEditor  8   /* ���� - run-������ ����� ������������������ �������� ����� */

  /* ���� ������ (��� ����������, ���������� �������, �������) */
#define    t_none              0    /* ��������� ����, � ������� �� ��������� */
#define    t_src               1    /* ����-�������� ������� */
#define    t_fun               2    /* �������������� ���� */
#define    t_dst               3    /* ����-�������� ���������� */
#define    t_del               4    /* ����� ������������ */
#define    t_ext               5    /* �����-�������������� */
#define    t_der               6    /* �����-����������� */
#define    t_imp               7    /* �����-��������� ������ */
#define    t_exp               8    /* �����-���������� ������ */

  /* ��������� ���������� ������� */
#define    r_Success        0       /* ��� ������ */
#define    r_Fail           1       /* �������� ������ */

  /* ������ - ������� ������� � DLL */
#define EXPORTED_FUNC int

  /* ����������� ������� �������� */
typedef struct { 
 double re;
 double im;
} complex_64;

  /* ����������� ��������� �������� */
typedef struct { 
 float re;
 float im;
} complex_32;

  /* ������ ���������� */
typedef struct {
  char* name;
  int   data_type;
  int   dim[3];
  int   index;
  int   direction;
  char* description;
  void* default_ptr;
  int   data_size;
} ext_var_info_record;

typedef int (*t_glob_obj_destructor_proc)(void* aGlobalObjPtr);

  /* �������� ��������� ��� ������� � ����������� ���������� � ������� �������� */
typedef struct {
  void*   LayerContext;                  /* ��������� �� �������� �������� */
    /* ��������� �������������� - ���������� ����� */
  char    IntMet;                        /* ����� �������������� */
  char    LoopMet;                       /* ����� ������� ������� ��� */
  char    IsLoop;                        /* ���� ������ �������������� (����� ��� �������) - ��� ����� - True */
  int     MaxLoopIt;                     /* ������������ ����� �������� ��� ������� ������� ��� */
  double  AbsErr;                        /* ���������� ������ */
  double  RelErr;                        /* ������������� ������ */
  double  Hmin;                          /* ����������� ��� �������������� */
  double  Hmax;                          /* ������������ ��� �������������� */
    /*���������, �������� ��� ������������� ������� */
  char    fPrecition;                    /* ���� ������ ������������� */
  char    fOneStep;                      /* ���� ���������� ������ ���� �������������� */
  char    fFirstStep;                    /* ���� ������� ���� ������� */
    /*���������� ���������� ��������� */
  double  newstep;                       /* ����� ���������� ��� �������������� */
  char    fsetstep;                      /* ���� - ���������� ����� ��� �������������� */
    /* ��� ������ ������������ ��� ���� ����� ���������������� � �������
    ������������������ ������� ����� (�������� ������������� �������� ���. ���������)
    ����� ���������� ������ �� �����    */
  void*   (*FindGlobalObject)(void* ALayerContext,char* aGlobObjectName);
    /* ���������������� ����� ���������� ������ */
  void    (*RegisterGlobalObject)(void* ALayerContext,char* aGlobObjectName,void* aNewObject,t_glob_obj_destructor_proc destructor_proc_ptr);
    /* ����������� ������ ���������� � ��������� �� �� ������� */
  unsigned int (*DoLoadNeedPlugin)(char* aPluginName);    
    /* ���� ������������� ���������� ����  */
  char    fNeedIter;  
    /* ����������� ����� ��������� ������� ����� (��� �������� ���������� �� ��������������� �����) */
  unsigned int  ShemeHash;
    /* ������ �� ���������� ������ ����������� ������������ � �����-� ���-���� */
  void*   GlobalDepList;
  void*   GlobalDepHash;
    /* ���� ����� ����������� ������������ ���������� ��� ������ ������ � ������ �������� */
  char    UseSignalExtendedSort;
    /* ���� - ������ ��� ������� �������������� ������ ������������ ������� ������ � ������ �������� */
  char    ErrorOnSignalLoop;
    /* ���� ������������ ���������� ������ ���� "������� ���������� ���������" */
  char    UseConditionsExtendedSort;
    /* ������� ������ ������������ ��� ������ �������������� ������ */
  void*   CurentDepList;
    /* ������ �� ������ ������ ������ ��� ����������� �������� ������ */
  void*   GlobalWherewithList;
  void*   GlobalWherewithHash;
    /* ������������ ������ ������ ��� ����������� ������ "������ ��������"-"������ ��������" */
  char    UseSignalsPortReconnection;
    /* ���� ������ �������� �������� */
  char    fConstantCheckMode;
    /* ���� ������ ��������� ���� - ��� ���������� ������ ��������� ������ */
  char    fCodeGenMode;
    /* ��������� ��������� ������, ����������� ��� ��������� ������� */
  void*   TaskContext;
    /* ������� ������ ��������� �� ������ �� ����� �������, ������� = ��� ������, � ��������� */
  unsigned char  (*GetDataPtr)(void* TaskContext, char* aSignalName, void** DataPtr, int* dimension); 
    /* �������� ������������� �������������� ���������  */
  char    (*StopCheck)(void* TaskContext);  
    /* ���� - ����� ������ ��������� ������ �� ���� ������������� */
  char    WriteSignalsOnSyncStep;
    /* ���� - ������������� ������� ��\� ������� �������������� ������� */
  char    common_translation_flag;
    /* ��������� �� ���������� ���������� �������� ������� */
  double* fCurentTime;
    /* ��������� �� ���������� ���������� - ������� ��� �������� */
  double* fCurentStep;
    /* ���� - ���������� ������������� ������ � ������� ����� ������ */
  char    fStartAgain;
    /* ���� - ���������� ��������� ��������� ��������� ������ */
  char    fSaveModelState;
    /* ���� - ����� �������������� �������� �������� �� ������� ��� ������ f_InitState */
  char    fWriteSignalsOnInitState;
    /* ���� - ������������ �������� �������� ��� �������������� ���������� ��� DIRK � ����� ������� */
  char    UseAlgVarsStepControl;
    /* ���������� ������� ��������� �������� (��� �������������� ������� ��������) */
  int     NLocalIter;
    /* ���� - ������ �������������� ������� ���� ��� ���������� �������� ������
       ��� ������� ������ ������� �������� � ������� ������� ���� ��������� ���� ����, �����
       ����� ����������� f_GoodStep ������� f_UpdateOuts � ��� �� ����� ��������������  */
  char    fNeedUpdateOutsBeforeGoodStep;  
    /* ���� - ������������ ��������� ���� ��� ��������� ���������� ������� */
  char    fPreciseSrcStep;   

    /* ��� ���������� ������� ����������� ����  */  
  char*   DefaultLAESolverLibraryName;
    
} solver_struct;





