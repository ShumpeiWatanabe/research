// Hodgkin-Huxley方程式
void   setSeed() {srand((unsigned)time(NULL));}//     
double genRand() {return ((double)rand()/(RAND_MAX));}//(0,1]       
//boxMuller                
double gaussian(double myu, double sgm)//       
{
  static int    frag = 0;
  static double save = 0.0;
  
  if(frag == 0){
    double u_0 = genRand();
    double u_1 = genRand();
    double v_0 = myu + sgm*sqrt(-2*log(1-u_0))*cos(2*M_PI*u_1);
    double v_1 = myu + sgm*sqrt(-2*log(1-u_0))*sin(2*M_PI*u_1);
    save       = v_1;
    frag       = 1;
    return v_0;
  }
  else{
    frag = 0;
    return save;
  }  
}

double f1(double t, double *x){
  return(-120.0*x[1]*x[1]*x[1]*x[2]*(x[0]-50.0)//I_Na
	 -36.0*x[3]*x[3]*x[3]*x[3]*(x[0]+77.0) //I_K
	 -0.3*(x[0]+54.39)                     //I_L
	 +x[4]+x[5]);                          //I_ext+I_syn
}
double f2(double t, double *x){
  double a=-0.1*(x[0]+40.0)/(exp(-(x[0]+40.0)/10.0)-1.0);
  double b=4.0*exp(-(x[0]+65.0)/18.0);
  return(a*(1.0-x[1])-b*x[1]);
}
double f3(double t, double *x){
  double a=0.07*exp(-(x[0]+65.0)/20.0);
  double b=1.0/(exp(-(x[0]+35.0)/10.0)+1.0);
  return(a*(1.0-x[2])-b*x[2]);
}
double f4(double t, double *x){
  double a=-0.01*(x[0]+55.0)/(exp(-(x[0]+55.0)/10.0)-1.0);
  double b=0.125*exp(-(x[0]+65.0)/80.0);
  return(a*(1.0-x[3])-b*x[3]);
}

// アルファ関数
//#define TAU 0.5
double alpha(double t){
  double tau;
  tau=gaussian(0.5,0.1);
  if(tau<0)tau=0.1;  
// double tau=0.5;
  
  if(t>=0.0)return(t/tau*exp(1.0-t/tau));
  else return(0.0);
}
