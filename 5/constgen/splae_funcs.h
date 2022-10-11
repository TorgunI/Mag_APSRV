
 //**************************************************************************//
 //      Функции решения разреженных СЛАУ (изолированные)                    //
 //     при помощи простого метода для разреженных матриц                    //
 //**************************************************************************//
 //             Программист:        Тимофеев К.А.                            //
 //**************************************************************************//

#include "sit_sparce_interface.h"
 
 /* Сортировка разреженной СЛАУ */
int SparceSort(int ind,int NOfE,int p,int** Wi,g_real_type** U,
               int** CU,int* U_Count,int* RU_Count,g_real_type TU)
{		 
 int ir;
 int ic;
 int r;
 int i;
 int j;
 int k;
 int m;
 int n;
 int ii;
 int marc;
 int row; 
 g_real_type amax;
 g_real_type xmax;
 int Result;
 
 Result = 0;
 ir = 0;
 ic = 0;
 ii = 0;
 xmax = 0;
 
 /* число просматриваемых строк */
 if (NOfE > 0) {
   n = NOfE - ind;
   if(n > p) n = p;   
 }else{
   n = 0;              /* Если у нас уравнений нет, то и не надо ничего там смотреть !!! */
   return Result;
 };
 
 /* отсортировать оставшиеся после k-го шага уравнения так, чтобы n первых
    уравнений имели минимальное число ненулевых элементов */
 for(j=0;j<=n-1;j++){
  r = 0x0FFFFFF;
  row = ind+j;
  for(i=ind+j;i<= NOfE - 1;i++){
   k = U_Count[(Wi[0][i]-1)];     //число Н.Н.Э. в строке row
   if (k < r) {
    r = k;
    row = i;
   };
  };
  m = Wi[0][ind+j];
  Wi[0][ind+j] = Wi[0][row];
  Wi[0][row] = m;
 };

 r = 0x0FFFFFF;
 for(j=ind;j<=ind+n-1;j++){
  row = Wi[0][j] - 1;

  m = U_Count[row] - 1;
  if (m < 0) { /* число ннэ = 0 */
    Result = 1;
    ii = -1;
    ir = j;
    break;
  };

  /* найти для каждой из n просматриваемых строк максимальный по значению элемент */
  amax = 0;
  for(i=0;i <= m;i++) amax = fmax(amax,fabs(U[row][i]));
  
  /* если ненулевой член в просматриваемой строке удовлетворяет
   условию устойчивости, то выбрать или не выбрать его в качестве главного
   по критерию Марковица */
  for(i=0; i <= m;i++) {
   k = CU[row][i] - 1;         /* номер неизвестной */

   marc = m*(RU_Count[k] - 1); /* цена Марковица */
   if (fabs(TU*U[row][i]) >= amax) { /* удовлетворяет устойчивости */
    if (marc < r) {
      /* если цена Марковица элемента меньше текущей минимальной,
      то выбрать данный элемент в качестве ведущего */
      r = marc;/* минимальная найденная цена Марковица */
      ii = i;  /* индекс ведущего элемента в строке */
      ic = k;  /* номер ведущего столбца в строке */
      ir = j;  /* номер ведущей строки*/
      xmax = fabs(U[row][i]); /* максимальный элемент среди просматриваемых строк */
    }
    else{ 
     if (marc == r) {
      /* если цена Марковица равна текущей минимальной цене, но
         элемент больше по значению, чем текущий ведущий, то
         выбрать данный элемент в качестве ведущего */
       if (fabs(U[row][i]) > xmax) {
         r = marc; /* минимальная найденная цена Марковица */
         ii = i;   /* индекс ведущего элемента в строке */
         ic = k;   /* номер ведущего столбца в строке */
         ir = j;   /* номер ведущей строки */
         xmax = fabs(U[row][i]); /* максимальный элемент среди просматриваемых строк */
       };
     };
    };
   };
  };
 };
 
 row = Wi[0][ir];
 Wi[0][ir] = Wi[0][ind];
 Wi[0][ind] = row;
 Wi[1][ind] = ic+1;
 Wi[2][ind] = ii;
 
 return Result;
}; 


int changecount_r(g_real_type** data, int new_count)
{
	
  size_t allocated_size;	
  size_t new_size;

  new_size = sizeof(g_real_type)*new_count;
  /* определяем сколько раньше выделяли */
  if (*data != 0) {
	  allocated_size = malloc_usable_size(*data);
	  /* если выделяли больше чем надо, то ничего не делаем */
	  if (new_size > allocated_size) {
		  *data = realloc(*data, new_size);
	  };
  }
  else {
	  *data = malloc(new_size);
  };
  
  return new_count;
};

int changecount_i(int** data, int new_count)
{
	
  size_t allocated_size;	
  size_t new_size;

  new_size = sizeof(int)*new_count;
  /* определяем сколько раньше выделяли */
  if (*data != 0) {
	  allocated_size = malloc_usable_size(*data);
	  /* если выделяли больше чем надо, то ничего не делаем */
	  if (new_size > allocated_size) {
		  *data = realloc(*data, new_size);
	  };
  }
  else {
	  *data = malloc(new_size);
  };
  
  return new_count;
};

int del_element_r(g_real_type* data, int count_data, int element_index)
{

  int i;

  for(i=element_index + 1;i < count_data;i++) {
    data[i - 1] = data[i];
  };	  	
  count_data = count_data - 1;	  
  return count_data;
};

int del_element_i(int* data, int count_data, int element_index)
{

  int i;

  for(i=element_index + 1;i < count_data;i++) {
    data[i - 1] = data[i];
  };	  	
  count_data = count_data - 1;	  
  return count_data;
};

 /* процедура выполнения k-го шага исключения Гаусса */
int   ElimStep(g_real_type** U, int* U_Count, g_real_type** L, int* L_Count, int** CU, int** RU, int* RU_Count, int** CL,int** Wi,int ind,g_real_type TB)
{
   g_real_type akk;
   g_real_type aik;
   g_real_type aa;
   int i;
   int i1;
   int i2;
   int i3;
   int ik;
   int ki;
   int ii;
   int j1;
   int countk;
   int counti;
   int countr;
   char iflag;
   int k;
   int m;
   int n;
   int Result;
  
   Result = 0;
   ik = 0;   /* по умолчанию 0 */
   aik = 0;

   k = Wi[0][ind] - 1;  /* номер ведущей строки */
   m = Wi[1][ind] - 1;  /* номер ведущего столбца */
   n = Wi[2][ind];      /* индекс ведущего столбца в массиве U[k] */

   if (n < 0) { 
     return Result;
   };

   akk = U[k][n];     /* главный член к-го шага */

   if (akk == 0) {
     Result = 1;
     return Result;
   };

   countk   = U_Count[k];        /* число ненулевых членов в к-ом уравнении */
   countr   = RU_Count[m];       /* число ненулевых членов в m-ом столбце матрицы */
   if ((countk == 0) || (countr == 0)) {
	  return Result;
   };

   changecount_r(&L[k],countr - 1);
   L_Count[k] = changecount_i(&CL[k],countr - 1);
   
   ki = 0;
   for (i =0; i <= countr - 1;i++){
    iflag = 0;
    ii    = RU[m][i] - 1;  /* номер i-го уравнения */
    if (ii == k) { 
	  continue;
	};
    counti = U_Count[ii];  /* число ненулевых членов в i-ом уравнении */

    /* Найти элемент в ii-ом уравнении , находящийся в столбце m */
    for(i1 = 0; i1 <= counti - 1;i1++){
      if (CU[ii][i1] == m + 1) {
        ik  = i1;
        aik = U[ii][ik];
        break;
	  };
    };
	
    aa  = -aik/akk;        /* множитель для ii-го уравнения */
    L[k][ki] = aa;         /* элемент нижней треугольной матрицы (k,ii) */
    CL[k][ki] = ii + 1;    /* индекс элемента L(k,ii) */
    ki = ki + 1;
	
    /* цикл по всем ненулевым элементам ведущей (k-ой) строки */
    for (i1 = 0; i1 <= countk - 1;i1++){
      if (i1 == n) { 
	    continue;
      };			  
      i3 = 0;
      /* найти элемент в ii-ой строке, стоящий в том же столбце,
       что и текущий элемент k-ой строки; если найден, то заполнения
       нет, изменить значение элемента ii-ой строки,
       если не найден, то есть заполнение */
      for (i2 = 0; i2 <= counti - 1;i2++) {
        if (CU[k][i1] == CU[ii][i2]) {
          U[ii][i2] = U[ii][i2] + U[k][i1]*aa;
          i3 = 0;
		  break;
        } 
		else {
		  i3 = 1;
		};  
	  };
	  

      /* есть заполнение, если это первое заполнение
      для ii-ой строки, оно помещается в исключаемый элемент ii-ой
      строки, если не первое, происходит выделение дополнительной памяти */
      if (i3 == 1) {		  
      /* если заполнение меньше барьера, заполнения нет */	 
       if (fabs(U[k][i1]*aa) >= TB) {
        if (iflag) {
         j1 = U_Count[ii];
		 		 
         changecount_r(&U[ii], j1+1 );
         U_Count[ii] = changecount_i(&CU[ii], j1+1 );
         RU_Count[CU[k][i1] - 1] = changecount_i(&RU[CU[k][i1] - 1], RU_Count[CU[k][i1] - 1] + 1 );
		 
         U[ii][j1] = U[k][i1]*aa;
         CU[ii][j1] = CU[k][i1];
         RU[ CU[k][i1] - 1 ][ RU_Count[CU[k][i1] - 1] - 1 ] = ii + 1;
        }
        else {
         iflag = 1;
		 		 
         RU_Count[CU[k][i1] - 1] = changecount_i(&RU[CU[k][i1] - 1], RU_Count[CU[k][i1] - 1] + 1  );
		 		 
         U[ii][ik] = U[k][i1]*aa;
         CU[ii][ik] = CU[k][i1];
         RU[ CU[k][i1]-1 ][ RU_Count[CU[k][i1] - 1] - 1] = ii + 1;
        };
	   };	
	  };
	  
    };

    if (!iflag) {
      del_element_r(U[ii], U_Count[ii],  ik);                //удалить элемент с индексом ik
      U_Count[ii] = del_element_i(CU[ii], U_Count[ii], ik);  //удалить элемент с индексом ik
    };
	
   };

   for (i = 0; i <= countk - 1;i++){
    ki = CU[k][i] - 1;
    for(i1 = 0;i1 <= RU_Count[ki] - 1;i1++){
      if (RU[ki][i1] == k + 1){
        RU_Count[ki] = del_element_i(RU[ki], RU_Count[ki], i1);
        break;
      };
    };
   };
   
   RU_Count[m] = changecount_i(&RU[m],0);
   
   return Result;
};
 
int  SparceLUDecomp(int NOfE, int Rows, int** Wi,
          g_real_type** U, int* U_Count, g_real_type** L, int* L_Count, int** CU, int** RU, int* RU_Count, int** CL,g_real_type TB,g_real_type TU,char first)
{
	
 int i;
 int j;
 int Result;

 Result = 0;

 if (first) {
  for (i = 0; i <= NOfE-1; i++){
    for (j = 0 ; j < 2 ; j++) {
    	Wi[j][i] = i + 1;		
	};
  }; 
 }  

 for (i = 0; i <= NOfE - 2; i++) {
 
  if(first) {
    Result = SparceSort(i,NOfE,Rows,Wi,U,CU,U_Count,RU_Count,TU);
    if (Result > 0) {
		return Result;
	};
  };

  /* выполнить i-ый шаг гауссова исключения */
  Result = ElimStep(U,U_Count,L,L_Count,CU,RU,RU_Count,CL,Wi,i,TB);
  if (Result > 0) {
	  return Result;
  };
  
 };

 if (first) {
   Result = SparceSort(NOfE-1,NOfE,Rows, Wi, U, CU, U_Count, RU_Count, TU);
   if (Result > 0) {
	   return Result;
   }
 };
 
 return Result;
		
};		  

int  SparceSolve(int NOfE,g_real_type* X,g_real_type* F,int** Wi,
          g_real_type** U, int* U_Count, g_real_type** L, int* L_Count, int** CU, int** CL)
{
	
	int m;
	int n;
	int k;
	int j;
    int i;
	int Result;
	g_real_type sum;
	
	Result = 0;

	/* определение новых правых частей при гауссовом исключении */
	for (i = 1; i <= NOfE - 1; i++) {
		k = Wi[0][i - 1] - 1;
		for (j = 0; j <= L_Count[k] - 1;j++) {
			m = CL[k][j] - 1;
			if (m >= 0) F[m] = F[m] + F[k] * L[k][j];
		};
	};

	/*выполнить обратную прогонку, найти начальное приближение вектора решений*/
	for (i = NOfE - 1; i >= 0; i--) {
		k = Wi[0][i] - 1;
		m = Wi[1][i] - 1;
		n = Wi[2][i];
		sum = F[k];

		for (j = 0; j <= U_Count[k] - 1; j++) {
			if (j != n) {
				sum = sum - U[k][j] * X[CU[k][j] - 1];
			};
		};

		if ((n < 0) || (U[k][n] == 0)) {
			X[m] = 1;
			Result = 1;
		}
		else {
			X[m] = sum / U[k][n];
		};

	};
	
	return Result;
	
};		  

int  SparceIter(int NOfE, g_real_type* X, g_real_type* Xr, g_real_type* Fr, g_real_type* F, int** Wi,
	g_real_type** U, int* U_Count, g_real_type** Us, int* Us_Count,  g_real_type** L, int* L_Count, int** CU, int** CUs, int** CL,
	g_real_type eps,
	g_real_type abserr,
	g_real_type errdxstart,
	int maxiter
	)
{
	g_real_type	errdx;
	g_real_type errdx1;
	int i;
	int j;
	int k;
	int iter;
	int Result;
	char flag;

	/* Итерационный процесс состоит из трех этапов
	 Этап 2 : d(i) = U ^ (-1)*L ^ (-1)*r(i) }
	 Этап 1 : r(i) = f(i) - A*x(i)}
	 Этап 3 : x(i + 1) = x(i) + d(i) */

	iter = 0;
	errdx1 = errdxstart;
	flag = 0;

	do {

		for (i = 0; i <= NOfE - 1; i++) {
			F[i] = Fr[i];
			for (j = 0; j <= Us_Count[i] - 1; j++) F[i] = F[i] - Us[i][j] * X[CUs[i][j] - 1];
		};

		Result = SparceSolve(NOfE, Xr, F, Wi, U, U_Count, L, L_Count, CU, CL);
		if (Result > 0) {
			return Result;
		};

		errdx = 0;
		for (k = 0; k < NOfE - 1; k++) {

			/* errx: = errx + ABS(X[k]);
			   errdx: = errdx + ABS(Xr[k]); */

			errdx = fmax(errdx, fabs(Xr[k]) / (fabs(X[k]) + abserr));
			X[k] = X[k] + Xr[k];
		};

		if (errdx <= eps) {
			flag = 1;
		};

		/* if errdx <= eps*amax1(errx, abserr) then flag : = True; */

		if (iter > 0) {
			if (errdx >= 2 * errdx1) {
				/* tells(0,'Итерационное уточнение системы ЛАУ расходится'); */
				flag = 1;
			};
		};

		/*   tells(0, 'RelErr=' + FloatToStrF(errdx / amax1(errx, abserr), ffGeneral, 6, 2) + ' d(i)=' + FloatToStrF(errdx, ffGeneral, 6, 2) +
			' d(i-1)=' + FloatToStrF(errdx1, ffGeneral, 6, 2) + ' iter=' + IntToStr(iter)); */

		errdx1 = errdx;

		iter = iter + 1;

	}  while (!(flag || (iter >= maxiter)));

	return Result;
};
 
  /* Данные для решателя разреженной СЛАУ */
typedef struct {

	//Ссылка на интерфейс решателя  -  этот член должен быть всегда первым, т.к. по нему мы определяем что нам вытащить !!!
	p_sparce_solver_interface  solver_interface;

	g_real_type* WR[4];       //Рабочие массивы для метода решения разреженных матриц
	int* WI[4];
	g_real_type** WRA[3];
	int* WRA_Count[3];
	int** WIA[4];
	int* WIA_Count[4];
							  // ---- Настройки метода решения ----
	char Iter_flag;           // Флаг итерационного уточнения решения
	int  Rows;                // число просматриваемых строк
	int  sp_maxiter;          // максималльное к-во итераций для SparceIter
	g_real_type TB;           // барьер
	g_real_type TU;           // коэффициент численной устойчивости
	g_real_type sp_eps;       // допустимая относительная ошибка итераций
	g_real_type sp_abserr;    // допустимая абсолютная ошибка
	g_real_type sp_errdxstart;// максимальная начальная ошибка

    int fEquCount;            //К-во уравнений и правых частей
    void* FirstCoefBlock;     //Первый блок присвоения коэффициентов
    void* LastCoefBlock;      //Последний блок присвоения коэффициентов
    char fNeedSort;           //Флаг пересортировки
    char fNeedLU;             //Флаг необходимости сделать заново LU-разложение
    char fLUSuccess;          //Флаг корректности LU-разложения

	g_real_type* X;           //Массив правых частей
	g_real_type* B;       

	char allocated;           //Флаг - память была выделена

} TSparceSolverData;

TSparceSolverData* lae_solver_create();

// Уничтожение объекта - решателя СЛАУ 
int lae_solver_free(TSparceSolverData* obj_ptr)
{

	if (obj_ptr->allocated) lae_solver_freememory(obj_ptr);

	free(obj_ptr);

	return 0;
};


void free_and_null(void** aptr)
{

	if (*aptr != 0) {
		free(*aptr);
		*aptr = 0;
	}

};


/* Уничтожение памяти, выделенной под решатель */
int lae_solver_freememory(TSparceSolverData* solver_context)
{
	int i,j;

	free_and_null(&solver_context->X);
	free_and_null(&solver_context->B);

	for (i = 0; i < 4; i++) free_and_null(&solver_context->WR[i]);
	for (i = 0; i < 4; i++) free_and_null(&solver_context->WI[i]);
	for (i = 0; i < 3; i++) free_and_null(&solver_context->WRA_Count[i]);
	for (i = 0; i < 4; i++) free_and_null(&solver_context->WIA_Count[i]);

	for (i = 0; i < 3; i++)
		if (solver_context->WRA[i] != 0) {
			for (j = 0; j < solver_context->fEquCount; j++)
				free_and_null(&solver_context->WRA[i][j]);
			free_and_null(&solver_context->WRA[i]);
		};

	for (i = 0; i < 4; i++)
		if (solver_context->WIA[i] != 0) {
			for (j = 0; j < solver_context->fEquCount; j++)
				free_and_null(&solver_context->WIA[i][j]);
			free_and_null(&solver_context->WIA[i]);
		};

	return 0;
};
 
  //Возвращает количество уравнений глобальной линейной системы
int lae_solver_getlaecount(sp_solve_id lae_solver_id)
{

	int Result;

    Result = 0;
	if (lae_solver_id != 0) {
		Result = ((TSparceSolverData*)(lae_solver_id))->fEquCount;
	};

	return Result;
}

 //Увеличивает количество уравнений глобальной линейной системы
int lae_solver_addlaecount(sp_solve_id lae_solver_id, int added_count)
{

	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		((TSparceSolverData*)(lae_solver_id))->fEquCount = ((TSparceSolverData*)(lae_solver_id))->fEquCount + added_count;	   
		Result = ((TSparceSolverData*)(lae_solver_id))->fEquCount;
	};

	return Result;
};

 //Добавить указатели на собственный блок
int lae_solver_setselfptr(sp_solve_id lae_solver_id, void* self_ptr)
{

	int Result;

	Result = 0;
	if (lae_solver_id != 0) {

		((TSparceSolverData*)(lae_solver_id))->LastCoefBlock = self_ptr;
		if (((TSparceSolverData*)(lae_solver_id))->FirstCoefBlock == 0) {
			((TSparceSolverData*)(lae_solver_id))->FirstCoefBlock = self_ptr;
		};

	};

	return Result;
};

 //Установить флаг пересортировки системы
int lae_solver_setneedsort(sp_solve_id lae_solver_id)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		((TSparceSolverData*)(lae_solver_id))->fNeedSort = 1;
		((TSparceSolverData*)(lae_solver_id))->fNeedLU = 1;
	};

	return Result;
};


 //Установить флаг необходимости LU-разложения матрицы системы уравнений
int lae_solver_setneedlu(sp_solve_id lae_solver_id)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		((TSparceSolverData*)(lae_solver_id))->fNeedLU = 1;
	};

	return Result;
};

//Установить флаг необходимости LU-разложения матрицы системы уравнений
int lae_solver_setzeroschange(sp_solve_id lae_solver_id)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		((TSparceSolverData*)(lae_solver_id))->fNeedSort = 1;
		((TSparceSolverData*)(lae_solver_id))->fNeedLU = 1;
	};

	return Result;
};


//Получить идентификатор последнего блока
void* lae_solver_getlastptr(sp_solve_id lae_solver_id)
{
	void* Result;

	Result = 0;
	if (lae_solver_id != 0) {
		Result = ((TSparceSolverData*)(lae_solver_id))->LastCoefBlock;
	};

	return Result;
};

//Получить идентификатор первого блока
void* lae_solver_getfirstptr(sp_solve_id lae_solver_id)
{
	void* Result;

	Result = 0;
	if (lae_solver_id != 0) {
		Result = ((TSparceSolverData*)(lae_solver_id))->FirstCoefBlock;
	};

	return Result;
};

 //Установить значение правой части в строке уравнения
int lae_solver_setrightcol(sp_solve_id lae_solver_id, int n_row, g_real_type bvalue)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		((TSparceSolverData*)(lae_solver_id))->B[n_row] = bvalue;
	};

	return Result;
};

 //Установить память под коэффициенты для указанной строки
int lae_solver_resetrow(sp_solve_id lae_solver_id, int n_row, int n_col)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {

        //Выделение памяти под матрицу
		changecount_r(&((TSparceSolverData*)(lae_solver_id))->WRA[0][n_row], n_col);
		changecount_r(&((TSparceSolverData*)(lae_solver_id))->WRA[1][n_row], n_col);
		changecount_i(&((TSparceSolverData*)(lae_solver_id))->WIA[0][n_row], n_col);
		changecount_i(&((TSparceSolverData*)(lae_solver_id))->WIA[1][n_row], n_col);

        //Обнуление счётчиков рабочих массивов
		((TSparceSolverData*)(lae_solver_id))->WRA_Count[0][n_row] = 0;
		((TSparceSolverData*)(lae_solver_id))->WRA_Count[1][n_row] = 0;
		((TSparceSolverData*)(lae_solver_id))->WIA_Count[0][n_row] = 0;
		((TSparceSolverData*)(lae_solver_id))->WIA_Count[1][n_row] = 0;

	};

	return Result;
};

//Установить память под коэффициенты для указанной строки
int lae_solver_initrow(sp_solve_id lae_solver_id, int n_row, int n_col, int* col_numbers)
{
	int Result;

	Result = 0;

	return Result;
};

 //Установить значение коэффициента
int lae_solver_fillitem(sp_solve_id lae_solver_id, int n_row, int index, g_real_type coef, int sum_cols, int order_index)
{
	int Result;
	int j;
	int num;
	int new_count;

	Result = 0;
	if (lae_solver_id != 0) {

		if ((coef != 0.0) && (index >= 0)) {
		    num = index + 1;
			if (sum_cols) {
				for (j = 0; j <= ((TSparceSolverData*)(lae_solver_id))->WIA_Count[1][n_row] - 1; j++) {
					if (((TSparceSolverData*)(lae_solver_id))->WIA[1][n_row][j] == num) {
						((TSparceSolverData*)(lae_solver_id))->WRA[1][n_row][j] = ((TSparceSolverData*)(lae_solver_id))->WRA[1][n_row][j] + coef;
						return Result;
					};
				};
			};
			// add col number
			new_count = ((TSparceSolverData*)(lae_solver_id))->WIA_Count[1][n_row];
			changecount_i(  &((TSparceSolverData*)(lae_solver_id))->WIA[1][n_row],  new_count + 1);
			((TSparceSolverData*)(lae_solver_id))->WIA_Count[1][n_row] = new_count + 1;
			((TSparceSolverData*)(lae_solver_id))->WIA[1][n_row][new_count] = num;

            // add col coefficient
			new_count = ((TSparceSolverData*)(lae_solver_id))->WRA_Count[1][n_row];
			changecount_r(&((TSparceSolverData*)(lae_solver_id))->WRA[1][n_row], new_count + 1);
			((TSparceSolverData*)(lae_solver_id))->WRA_Count[1][n_row] = new_count + 1;
			((TSparceSolverData*)(lae_solver_id))->WRA[1][n_row][new_count] = coef;

		};

	};

	return Result;
};


//Инициализация линейной системы
int lae_solver_begininit(sp_solve_id lae_solver_id)
{
	int Result;
	int i,j;

	Result = 0;
	if (lae_solver_id != 0) {
	
		((TSparceSolverData*)(lae_solver_id))->TU = 1;                         //Стандартные настроечные коэффициенты метода
		((TSparceSolverData*)(lae_solver_id))->TB = 0;                         //решения разреженных матриц
		((TSparceSolverData*)(lae_solver_id))->Rows = 3;
		((TSparceSolverData*)(lae_solver_id))->sp_eps = 0.0001;
		((TSparceSolverData*)(lae_solver_id))->sp_abserr = 1.0e-15;
		((TSparceSolverData*)(lae_solver_id))->sp_maxiter = 20;
		((TSparceSolverData*)(lae_solver_id))->sp_errdxstart = 1.0e30;
		((TSparceSolverData*)(lae_solver_id))->Iter_flag = 0;

		for (i = 0; i < 4; i++) {
			changecount_r(&((TSparceSolverData*)(lae_solver_id))->WR[i], ((TSparceSolverData*)(lae_solver_id))->fEquCount);
		};
		for (i = 0; i < 4; i++) {
			changecount_i(&((TSparceSolverData*)(lae_solver_id))->WI[i], ((TSparceSolverData*)(lae_solver_id))->fEquCount);
		};
		for (i = 0; i < 3; i++) {
			changecount_i(&((TSparceSolverData*)(lae_solver_id))->WRA_Count[i], ((TSparceSolverData*)(lae_solver_id))->fEquCount);
			((TSparceSolverData*)(lae_solver_id))->WRA[i] = malloc( sizeof(void*)*((TSparceSolverData*)(lae_solver_id))->fEquCount);
			memset(((TSparceSolverData*)(lae_solver_id))->WRA[i], 0, sizeof(void*)*((TSparceSolverData*)(lae_solver_id))->fEquCount);
			for (j = 0; j < ((TSparceSolverData*)(lae_solver_id))->fEquCount; j++) {
				((TSparceSolverData*)(lae_solver_id))->WRA_Count[i][j] = 0;
			};
		};
		for (i = 0; i < 4; i++) {
			changecount_i(&((TSparceSolverData*)(lae_solver_id))->WIA_Count[i], ((TSparceSolverData*)(lae_solver_id))->fEquCount);
			((TSparceSolverData*)(lae_solver_id))->WIA[i] = malloc(sizeof(void*)*((TSparceSolverData*)(lae_solver_id))->fEquCount);
			memset(((TSparceSolverData*)(lae_solver_id))->WIA[i], 0, sizeof(void*)*((TSparceSolverData*)(lae_solver_id))->fEquCount);
			for (j = 0; j < ((TSparceSolverData*)(lae_solver_id))->fEquCount; j++) {
				((TSparceSolverData*)(lae_solver_id))->WIA_Count[i][j] = 0;
			};
		};

		changecount_r(&((TSparceSolverData*)(lae_solver_id))->B, ((TSparceSolverData*)(lae_solver_id))->fEquCount);
		changecount_r(&((TSparceSolverData*)(lae_solver_id))->X, ((TSparceSolverData*)(lae_solver_id))->fEquCount);
		((TSparceSolverData*)(lae_solver_id))->fNeedSort = 1;
		((TSparceSolverData*)(lae_solver_id))->fNeedLU = 1;

	};

	return Result;
};

  //Завершение инициализации линейной системы
int lae_solver_endinit(sp_solve_id lae_solver_id)
{
	int Result;

	Result = 0;

	return Result;
};

 //Инициализация LU-разложения
int  lae_solver_init(sp_solve_id lae_solver_id)
{
	int Result;
	int i;
	int j;
    int k;
	int N;

	Result = 0;
	if (lae_solver_id != 0) {

		N = ((TSparceSolverData*)(lae_solver_id))->fEquCount;

		//Выделение памяти под матрицу
		for (j = 0; j <= N - 1; j++) {
			((TSparceSolverData*)(lae_solver_id))->WRA_Count[0][j] = changecount_r(&((TSparceSolverData*)(lae_solver_id))->WRA[0][j], ((TSparceSolverData*)(lae_solver_id))->WRA_Count[1][j]);
			((TSparceSolverData*)(lae_solver_id))->WIA_Count[0][j] = changecount_i(&((TSparceSolverData*)(lae_solver_id))->WIA[0][j], ((TSparceSolverData*)(lae_solver_id))->WRA_Count[1][j]);
		};
		for (j = 0; j <= N - 1; j++) {
			memcpy(((TSparceSolverData*)(lae_solver_id))->WRA[0][j], ((TSparceSolverData*)(lae_solver_id))->WRA[1][j], ((TSparceSolverData*)(lae_solver_id))->WRA_Count[0][j] * sizeof(g_real_type));
			memcpy(((TSparceSolverData*)(lae_solver_id))->WIA[0][j], ((TSparceSolverData*)(lae_solver_id))->WIA[1][j], ((TSparceSolverData*)(lae_solver_id))->WRA_Count[0][j] * sizeof(int));
		};


		for (j = 0; j <= N - 1; j++) {
			((TSparceSolverData*)(lae_solver_id))->WI[3][j] = 0;
		};
		for (j = 0; j <= N - 1; j++) {
			for (i = 0; i <= ((TSparceSolverData*)(lae_solver_id))->WIA_Count[0][j] - 1; i++) {
				k = ((TSparceSolverData*)(lae_solver_id))->WIA[0][j][i] - 1;
				if (k < 0) {
					Result = 1;
					return Result;
				};
				((TSparceSolverData*)(lae_solver_id))->WI[3][k] = ((TSparceSolverData*)(lae_solver_id))->WI[3][k] + 1;
			};
		};

		for (j = 0; j <= N - 1; j++) {
			((TSparceSolverData*)(lae_solver_id))->WIA_Count[2][j] = changecount_i(&((TSparceSolverData*)(lae_solver_id))->WIA[2][j], ((TSparceSolverData*)(lae_solver_id))->WI[3][j]);
			((TSparceSolverData*)(lae_solver_id))->WI[3][j] = 0;
		};

		for (j = 0; j <= N - 1; j++) {
			for (i = 0; i <= ((TSparceSolverData*)(lae_solver_id))->WIA_Count[0][j] - 1; i++) {
			    k = ((TSparceSolverData*)(lae_solver_id))->WIA[0][j][i] - 1;
				if (k < 0) {
					Result = 1;
					return Result;
				};
				((TSparceSolverData*)(lae_solver_id))->WIA[2][k][((TSparceSolverData*)(lae_solver_id))->WI[3][k] ] = j + 1;
				((TSparceSolverData*)(lae_solver_id))->WI[3][k] = ((TSparceSolverData*)(lae_solver_id))->WI[3][k] + 1;
			};
		};

	};

	return Result;
};


 //LU-разложение
int lae_solver_ludecomp(sp_solve_id lae_solver_id) 
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {
		
	  //Если выставлен флаг необходимости разложения, то делаем его	
	  if(((TSparceSolverData*)(lae_solver_id))->fNeedLU){	

		//Инициализация разложения
		lae_solver_init(lae_solver_id);

		//LU - разложение матрицы (возможно с пересортировкой)
		Result = SparceLUDecomp(((TSparceSolverData*)(lae_solver_id))->fEquCount, 
			                    ((TSparceSolverData*)(lae_solver_id))->Rows, 
			                    &(((TSparceSolverData*)(lae_solver_id))->WI), 
			                     ((TSparceSolverData*)(lae_solver_id))->WRA[0], 
			                     ((TSparceSolverData*)(lae_solver_id))->WRA_Count[0],
		           	             ((TSparceSolverData*)(lae_solver_id))->WRA[2], 
			                     ((TSparceSolverData*)(lae_solver_id))->WRA_Count[2],
			                     ((TSparceSolverData*)(lae_solver_id))->WIA[0], 
			                     ((TSparceSolverData*)(lae_solver_id))->WIA[2], 
			                     ((TSparceSolverData*)(lae_solver_id))->WIA_Count[2],
			                     ((TSparceSolverData*)(lae_solver_id))->WIA[3], 
			                     ((TSparceSolverData*)(lae_solver_id))->TB, 
			                     ((TSparceSolverData*)(lae_solver_id))->TU, 
			                     ((TSparceSolverData*)(lae_solver_id))->fNeedSort);

		((TSparceSolverData*)(lae_solver_id))->fLUSuccess = (Result == 0);
		
		//Сброс флагов необходимости LU-разложения матрицы, если всё вычислено успешно
		if (Result == 0){
	      ((TSparceSolverData*)(lae_solver_id))->fNeedSort = 0;
		  ((TSparceSolverData*)(lae_solver_id))->fNeedLU = 0;
		};  

	  };	
		
	};

	return Result;
};


 //Решение системы
int lae_solver_solve(sp_solve_id lae_solver_id)
{
	int Result;

	Result = 0;
	if (lae_solver_id != 0) {


		//Заполняем вектор правых частей
		memcpy(((TSparceSolverData*)(lae_solver_id))->WR[0], ((TSparceSolverData*)(lae_solver_id))->B, ((TSparceSolverData*)(lae_solver_id))->fEquCount*sizeof(g_real_type));

		//Вычисляем решение системы уравнений
		if (((TSparceSolverData*)(lae_solver_id))->fLUSuccess) {


			memcpy(((TSparceSolverData*)(lae_solver_id))->WR[1], ((TSparceSolverData*)(lae_solver_id))->WR[0], ((TSparceSolverData*)(lae_solver_id))->fEquCount*sizeof(g_real_type));

		    Result = SparceSolve(
				((TSparceSolverData*)(lae_solver_id))->fEquCount,
				((TSparceSolverData*)(lae_solver_id))->WR[2],
				((TSparceSolverData*)(lae_solver_id))->WR[1],
				&((TSparceSolverData*)(lae_solver_id))->WI,
				((TSparceSolverData*)(lae_solver_id))->WRA[0],
				((TSparceSolverData*)(lae_solver_id))->WRA_Count[0],
				((TSparceSolverData*)(lae_solver_id))->WRA[2],
				((TSparceSolverData*)(lae_solver_id))->WRA_Count[2],
				((TSparceSolverData*)(lae_solver_id))->WIA[0],
				((TSparceSolverData*)(lae_solver_id))->WIA[3]);

			if ((Result == 0) && (((TSparceSolverData*)(lae_solver_id))->Iter_flag)) {

			  Result = SparceIter(
            	((TSparceSolverData*)(lae_solver_id))->fEquCount,
				  ((TSparceSolverData*)(lae_solver_id))->WR[2],
				  ((TSparceSolverData*)(lae_solver_id))->WR[3],
				  ((TSparceSolverData*)(lae_solver_id))->WR[0],
				  ((TSparceSolverData*)(lae_solver_id))->WR[1],
				  &((TSparceSolverData*)(lae_solver_id))->WI,
				  ((TSparceSolverData*)(lae_solver_id))->WRA[0],
				  ((TSparceSolverData*)(lae_solver_id))->WRA_Count[0],
				  ((TSparceSolverData*)(lae_solver_id))->WRA[1],
				  ((TSparceSolverData*)(lae_solver_id))->WRA_Count[1],
				  ((TSparceSolverData*)(lae_solver_id))->WRA[2],
				  ((TSparceSolverData*)(lae_solver_id))->WRA_Count[2],
				  ((TSparceSolverData*)(lae_solver_id))->WIA[0],
				  ((TSparceSolverData*)(lae_solver_id))->WIA[1],
				  ((TSparceSolverData*)(lae_solver_id))->WIA[3],
				  ((TSparceSolverData*)(lae_solver_id))->sp_eps,
				  ((TSparceSolverData*)(lae_solver_id))->sp_abserr,
				  ((TSparceSolverData*)(lae_solver_id))->sp_errdxstart,
				  ((TSparceSolverData*)(lae_solver_id))->sp_maxiter
				);
			};

			if (Result == 0) {
				memcpy(((TSparceSolverData*)(lae_solver_id))->X, ((TSparceSolverData*)(lae_solver_id))->WR[2], ((TSparceSolverData*)(lae_solver_id))->fEquCount*sizeof(g_real_type));
			};

		};
	};

	return Result;
};

 //Считывание результатов решения линейной системы
g_real_type lae_solver_getresult(sp_solve_id lae_solver_id,int index)
{
	g_real_type Result;

	Result = 0;
	if (lae_solver_id != 0) {
		Result = ((TSparceSolverData*)(lae_solver_id))->X[index];
	};

	return Result;
};

//Присвоение начального значения решения системы после чтения рестарта
int lae_solver_setinitresult(sp_solve_id lae_solver_id, int index, g_real_type x0)
{
	int Result;
	
	if (lae_solver_id != 0) {
      //Записываем начальное значение решения 		
      ((TSparceSolverData*)(lae_solver_id))->X[index]=x0;
	  //Выставляем флаг необходимости обновить матрицу (на всякий случай)
	  ((TSparceSolverData*)(lae_solver_id))->fNeedLU = 1;
    }; 
	
	Result = 0;

	return Result;
};

/* Инициализация интерфейса решателя разреженных матриц */
TSparceSolverAbstractInterface SP_SIMPLE_INTERFACE = {
	&lae_solver_create,
	&lae_solver_free,
	&lae_solver_getlaecount,
	&lae_solver_addlaecount,
	&lae_solver_setselfptr,
	&lae_solver_getlastptr,
	&lae_solver_getfirstptr,
	&lae_solver_begininit,
	&lae_solver_initrow,
	&lae_solver_setinitresult,
	&lae_solver_endinit,
	&lae_solver_setrightcol,
	&lae_solver_resetrow,
	&lae_solver_fillitem,
	&lae_solver_ludecomp,
	&lae_solver_solve,
	&lae_solver_getresult,
	&lae_solver_setneedsort,
	&lae_solver_setneedlu,
	&lae_solver_setzeroschange
};

// Поиcк и регистрация глобального объекта-решателя СЛАУ
// Возвращает идентификатор глобальной линейной системы lae_solver_id, если 0 - нет такой системы
TSparceSolverData* lae_solver_create()
{

	TSparceSolverData* solv_ptr;

	solv_ptr = malloc(sizeof(TSparceSolverData));

	memset(solv_ptr, 0, sizeof(TSparceSolverData));

	//Ссылка на интерфейс решалки
	solv_ptr->solver_interface = &SP_SIMPLE_INTERFACE;

	return solv_ptr;

};

  //Возвращает идентификатор глобальной линейной системы lae_solver_id, если 0 - нет такой системы
sp_solve_id lae_solver_register(solver_struct* solver_data, char* lae_solver_name, char* library_name, p_sparce_solver_interface* sp_solver_interface)
{

	sp_solve_id solv_ptr;

	solv_ptr = lae_solver_create();
	
	*sp_solver_interface = &SP_SIMPLE_INTERFACE;
	
	return solv_ptr;

};

 //Общая инициализация модуля решателя системы разреженных ЛАУ - в данной случае эта программа пустая !
void lae_solver_initialization(solver_struct* solver_data)
{


 return;
};

 //Выгрузка биболиотеки если не нужна - в данной случае эта программа пустая !
void lae_solver_finalization(solver_struct* solver_data)
{


 return;
};
