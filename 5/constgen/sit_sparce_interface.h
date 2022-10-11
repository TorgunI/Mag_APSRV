  
  /* ��� ��������� �� ������-�������� ����������� ���� - �� ��������� ��� ������� ��������� */
typedef void* sp_solve_id;

#define sparce_solver_call_type 

  /* �������� ���������� �������� ����������� ������ */
typedef struct {

  //�������� ���������� �������� ����������� ������ � �������� �����������
  sp_solve_id (sparce_solver_call_type *lae_solver_create)();
  //����������� ���������� �������-�������� ����������� �������  
  int (sparce_solver_call_type *lae_solver_free)(sp_solve_id lae_solver_id);  
    
  //���������� ���������� ��������� ���������� �������� �������
  int (sparce_solver_call_type *lae_solver_getlaecount)(sp_solve_id lae_solver_id);
  //����������� ���������� ��������� ���������� �������� �������
  int (sparce_solver_call_type *lae_solver_addlaecount)(sp_solve_id lae_solver_id, int added_count);
  
  //�������� ��������� �� ����������� ����
  int (sparce_solver_call_type *lae_solver_setselfptr)(sp_solve_id lae_solver_id, void* self_ptr);
  //�������� ������������� ���������� �����
  void* (sparce_solver_call_type *lae_solver_getlastptr)(sp_solve_id lae_solver_id);
  //�������� ������������� ������� �����
  void* (sparce_solver_call_type *lae_solver_getfirstptr)(sp_solve_id lae_solver_id);
    
  //������������� �������� ������� 
  int (sparce_solver_call_type *lae_solver_begininit)(sp_solve_id lae_solver_id);
  //������������� ������� � ������� ������ ���������
  int (sparce_solver_call_type *lae_solver_initrow)(sp_solve_id lae_solver_id, int n_row, int n_col, int* col_indexes);  
  //���������� � ����� ���������� ��������� (��� ������������ �������)
  int (sparce_solver_call_type *lae_solver_setinitresult)(sp_solve_id lae_solver_id, int index, double x0);
  //������������� LU-����������
  int (sparce_solver_call_type *lae_solver_endinit)(sp_solve_id lae_solver_id);
   
  //���������� �������� ������ ����� � ������ ���������
  int (sparce_solver_call_type *lae_solver_setrightcol)(sp_solve_id lae_solver_id,int n_row,double bvalue);
  //�������� ������������ ������ � 0
  int (sparce_solver_call_type *lae_solver_resetrow)(sp_solve_id lae_solver_id,int n_row,int n_col);
  //���������� �������� ������������ ������ � �������� ���������� ������� order_index
  int (sparce_solver_call_type *lae_solver_fillitem)(sp_solve_id lae_solver_id,int n_row,int index,double coef,int sum_cols,int order_index);
  
  //LU-����������
  int (sparce_solver_call_type *lae_solver_ludecomp)(sp_solve_id lae_solver_id);
  //������� �������
  int (sparce_solver_call_type *lae_solver_solve)(sp_solve_id lae_solver_id);
  //���������� ����������� ������� �������� �������
  double (sparce_solver_call_type *lae_solver_getresult)(sp_solve_id lae_solver_id,int index);
  
  //���������� ���� �������������� ������� (� ���������� ��� ������ ����� ����� ������ �� �����������, �.�. �������� �� ���������������� ����� ����������)
  int (sparce_solver_call_type *lae_solver_setneedsort)(sp_solve_id lae_solver_id);
  //���������� ���� ������������� LU-���������� ������� ������� ���������
  int (sparce_solver_call_type *lae_solver_setneedlu)(sp_solve_id lae_solver_id);
  //���������� ���� ������������� LU-���������� ������� ������� ��� ��������� �-�� ������� ���������
  int(sparce_solver_call_type *lae_solver_setzeroschange)(sp_solve_id lae_solver_id);
    
} TSparceSolverAbstractInterface;

/* ������ �� ��������� �������� ����������� ������ */
typedef TSparceSolverAbstractInterface* p_sparce_solver_interface;

  //������� ��������� ��� �������, �������� ���������
typedef struct {

	//������ �� ��������� ��������  -  ���� ���� ������ ���� ������ ������, �.�. �� ���� �� ���������� ��� ��� �������� !!!
	p_sparce_solver_interface  solver_interface;

} TSolverInterfacedObject;

  /* ��������� ���������� �� ��������� ������� */
p_sparce_solver_interface lae_solver_getinterface(sp_solve_id lae_solver_id)
{
	void* result;

	result = ((TSolverInterfacedObject*)(lae_solver_id))->solver_interface;

	return result;
};