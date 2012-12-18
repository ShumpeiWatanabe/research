#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include "hh.h"
#include "rungekutta.h"


#define V_th -40.0 // 発火判定閾電位．最初は-40に設定してあった．-60
#define V_rest -65
#define TIME_DIV 0.01 // RungeKutta刻み幅
#define NUM_VAR1 4 // ニューロンの状態変数数(ルンゲクッタ更新)
#define NUM_VAR2 3 // ニューロンの状態変数数(ルンゲクッタ非更新)
#define V_EXI (V_rest+15.0) // ＋なら興奮性，ーなら抑制性シナプス結合逆転電位]
#define V_INT (V_rest-15.0)
#define TIME_DIV 0.01
#define TH 0.8//ニューロン同士が興奮性でつながる確率
#define TH_INT 0.2//ニューロン同士が抑制性でつながる確率．
#define LEFTBRAIN 100
#define RIGHTBRAIN 100
#define NEURON_NUMBER LEFTBRAIN+RIGHTBRAIN
#define NEURON_VARIABLE 7
#define TIME 1000
#define G 0.3//普通は０．３ 0.6
//ニューロン1つを表す構造体

struct neuron{
  double variable[2][NEURON_VARIABLE];//ニューロン1っこの変数[0][]が現在の数値,[1][]が１つ前の数値
  int connectionFlag[NEURON_NUMBER];//どのニューロンからつなられているか
};

void update(double t,struct neuron brain[]){
  int i,j;
  static int entered=0;
  static double (*f[4])(double, double*);
  //関数が実行される最初のみ実行
  if(entered==0){
    f[0]=f1;
    f[1]=f2;
    f[2]=f3;
    f[3]=f4;
    entered =1;
  }
  //  printf("function set\n");
  //現在の変数を1つ前の変数に退避させる
  for(i=0;i < NEURON_NUMBER;i++){
      for(j=0; j < NEURON_VARIABLE; j++){
	brain[i].variable[1][j]=brain[i].variable[0][j];
      }
    }
    //  printf("variable moved\n");
    //HH方程式の4状態変数(V,m,h,n)を更新
    for(i=0;i<NEURON_NUMBER;i++){
      RungeKutta(f,t,TIME_DIV,&brain[i].variable[0][0],NUM_VAR1,NUM_VAR2);
    } 
    //  printf("rungekutta end\n");
    //発火時間変数更新

    for(i=0;i<NEURON_NUMBER;i++){
      if(brain[i].variable[1][0]<=V_th&&brain[i].variable[0][0]>V_th){
	brain[i].variable[0][6]=t;
      }
    }
    //  printf("fied time set\n");
    //シナプス性電流の更新
    
    for(i=LEFTBRAIN;i<NEURON_NUMBER;i++){
      brain[i].variable[0][5]=0.0;
      for(j=0;j<NEURON_NUMBER;j++){
	if(brain[i].connectionFlag[j]==1){
	  brain[i].variable[0][5]+=-G*alpha(t-brain[j].variable[0][6])*(brain[i].variable[0][0]-V_EXI);
	}
	if(brain[i].connectionFlag[j]==-1){
	  brain[i].variable[0][5]+=-G*alpha(t-brain[j].variable[0][6])*(brain[i].variable[0][0]-V_INT);
	}
	//	else{
	//  brain[i].variable[0][5]+=-G*alpha(t-brain[j].variable[0][6])*(brain[i].variable[0][0]-V_INT);
	//	}
	/*if(){}の部分をbrain[i].variable[0][5]+=-0.3*alpha(t-brain[j].variable[1][6])*(brain[j].variable[1][0]-V_EXI)*flag[];としてはいけない.nazeda!*/
      }
  }
    //  printf("synapus I set\n");
}

// 
int judgeConnection(double th){
   if(((double)rand() / ((double)RAND_MAX + 1))<th)return 1;
   if(((double)rand() / ((double)RAND_MAX + 1))>=th&&((double)rand() / ((double)RAND_MAX + 1))<th+TH_INT)return -1;
  else return 0;
  // return 1;
}

void init(struct neuron brain[]){
  int i,j;
  
  for(i=0;i<NEURON_NUMBER;i++){
    brain[i].variable[0][0]=-34.063412;// HH第1変数V(膜電位
    brain[i].variable[0][1]=0.906009;   // HH第2変数m
    brain[i].variable[0][2]=0.630484;   // HH第3変数h
    brain[i].variable[0][3]=1.066412;   // HH第4変数n
    brain[i].variable[0][4]=0.0;   // 外部刺激電流I_ext
    brain[i].variable[0][5]=0.0;   // シナプス性電流変数I_syn
    brain[i].variable[0][6]=-100.0;// 発火時間変数t_fire(過去の発火の影響なし)
  }
  
  printf("variable set\n");
  //左脳のニューロンの設定．
  //外部刺激電流．いつもパルスが立つようにする．
  for(i=0; i < LEFTBRAIN; i++){
    brain[i].variable[0][4]=9.0;
    brain[i].variable[1][4]=9.0;
  }
  //左脳のニューロンのつながり設定．左脳同士は繋がられていない．
  for(i=0;i<LEFTBRAIN;i++){
    for(j=0; j<NEURON_NUMBER;j++){
      brain[i].connectionFlag[j]=0;
    }
  }
  
  //  printf("left brain set\n");
  //右脳のニューロンの設定　右脳はすべてのニューロンと一定の割合で繋がられている　
  for(i=LEFTBRAIN; i<NEURON_NUMBER; i++){
    for(j=0; j<NEURON_NUMBER;j++){
      brain[i].connectionFlag[j]=judgeConnection(TH);
    }
  }
  //右脳の最後のニューロンの設定．左脳とつながっている．
  //  printf("right brain set\n");
}

int main(void){
  struct neuron brain[NEURON_NUMBER];
  double t=0;
  double sum,sum5,suml;
  int i,j;
  init(brain);
  printf("init set\n"); 
  for(t=0; t<TIME;t+=TIME_DIV){
    update(t,brain);
    //  printf("update end\n");
    sum=0;
    sum5=0;
    suml=0;
    for(i=0;i<LEFTBRAIN;i++){
      suml+=brain[1].variable[0][0];
    }
    for(i=LEFTBRAIN;i < NEURON_NUMBER; i++){
      sum+=brain[i].variable[0][0];
      sum5+=brain[i].variable[0][5];
    }

    printf("%lf %lf %lf %lf\n",t,suml/LEFTBRAIN,sum/RIGHTBRAIN,sum5/RIGHTBRAIN);
  }
   for(i=0; i<NEURON_NUMBER;i++){
     for(j=0; j<NEURON_NUMBER;j++){
     printf("%d",brain[i].connectionFlag[j]);
     }
     printf("\n");
     }
  return(0);
}
