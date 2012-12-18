void RungeKutta(double (*f[])(double t, double *x), 
		double t0, double h, double *x, 
		int num1, // ルンゲ=クッタ法で更新する変数の数
		int num2  // ルンゲ=クッタ法で更新しない変数の数
		){

  double *k1, *k2, *k3, *k4, *temp;
  int j;
  
  if((k1 = (double*)malloc(sizeof(double)*num1)) == NULL){
    fprintf(stderr, "Allocation Error\n");
    exit(1);}
  if((k2 = (double*)malloc(sizeof(double)*num1)) == NULL){
    fprintf(stderr, "Allocation Error\n");
    exit(1);}
  if((k3 = (double*)malloc(sizeof(double)*num1)) == NULL){
    fprintf(stderr, "Allocation Error\n");
    exit(1);}
  if((k4 = (double*)malloc(sizeof(double)*num1)) == NULL){
    fprintf(stderr, "Allocation Error\n");
    exit(1);}
  if((temp = (double*)malloc(sizeof(double)*(num1+num2))) == NULL){
    fprintf(stderr, "Allocation Error\n");
    exit(1);}

  for(j=0; j<num1; j++)
    k1[j] = (*f[j])(t0, x);
  for(j=0; j<num1; j++)
    temp[j] = x[j] + h*k1[j]/2;
  for(j=num1; j<num1+num2; j++)
    temp[j] = x[j];

  for(j=0; j<num1; j++)
    k2[j] = (*f[j])(t0+h/2, temp);
  for(j=0; j<num1; j++)
    temp[j] = x[j] + h*k2[j]/2;

  for(j=0; j<num1; j++)
    k3[j] = (*f[j])(t0+h/2, temp);
  for(j=0; j<num1; j++)
    temp[j] = x[j] + h*k3[j];

  for(j=0; j<num1; j++){
    k4[j] = (*f[j])(t0+h, temp);
    x[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*h/6;
  }

  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(temp);
}
