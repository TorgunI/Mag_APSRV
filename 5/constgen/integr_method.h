//  Стандартный метод интегрирования (пока тут зашит только метод Эйлера
void integr_method_call(g_real_type* x, g_real_type* dx,g_real_type step, int dimension, int common_offset, int method)
{
  	int i;
	
	for(i=0;i<dimension;i++){
		x[i]=x[i] + dx[i]*step;		
	};
	
	return;
}	
