#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include "hh.h"
#include "rungekutta.h"


#define V_th -40.0 // ȯ��Ƚ�����Ű̡��ǽ��-40�����ꤷ�Ƥ��ä���-60
#define V_rest -65
#define TIME_DIV 0.01 // RungeKutta�����
#define NUM_VAR1 4 // �˥塼���ξ����ѿ���(��󥲥��å�����)
#define NUM_VAR2 3 // �˥塼���ξ����ѿ���(��󥲥��å��󹹿�)
#define V_EXI (V_rest+15.0) // �ܤʤ鶽ʳ�������ʤ����������ʥץ�����ž�Ű�]
#define V_INT (V_rest-15.0)
#define TIME_DIV 0.01
#define TH 0.8//�˥塼���Ʊ�Τ���ʳ���ǤĤʤ����Ψ
#define TH_INT 0.2//�˥塼���Ʊ�Τ��������ǤĤʤ����Ψ��
#define LEFTBRAIN 100
#define RIGHTBRAIN 100
#define NEURON_NUMBER LEFTBRAIN+RIGHTBRAIN
#define NEURON_VARIABLE 7
#define TIME 1000
#define G 0.3//���̤ϣ����� 0.6
//�˥塼���1�Ĥ�ɽ����¤��

struct neuron{
  double variable[2][NEURON_VARIABLE];//�˥塼���1�ä����ѿ�[0][]�����ߤο���,[1][]���������ο���
  int connectionFlag[NEURON_NUMBER];//�ɤΥ˥塼��󤫤�Ĥʤ��Ƥ��뤫
};

void update(double t,struct neuron brain[]){
  int i,j;
  static int entered=0;
  static double (*f[4])(double, double*);
  //�ؿ����¹Ԥ����ǽ�Τ߼¹�
  if(entered==0){
    f[0]=f1;
    f[1]=f2;
    f[2]=f3;
    f[3]=f4;
    entered =1;
  }
  //  printf("function set\n");
  //���ߤ��ѿ���1�������ѿ������򤵤���
  for(i=0;i < NEURON_NUMBER;i++){
      for(j=0; j < NEURON_VARIABLE; j++){
	brain[i].variable[1][j]=brain[i].variable[0][j];
      }
    }
    //  printf("variable moved\n");
    //HH��������4�����ѿ�(V,m,h,n)�򹹿�
    for(i=0;i<NEURON_NUMBER;i++){
      RungeKutta(f,t,TIME_DIV,&brain[i].variable[0][0],NUM_VAR1,NUM_VAR2);
    } 
    //  printf("rungekutta end\n");
    //ȯ�л����ѿ�����

    for(i=0;i<NEURON_NUMBER;i++){
      if(brain[i].variable[1][0]<=V_th&&brain[i].variable[0][0]>V_th){
	brain[i].variable[0][6]=t;
      }
    }
    //  printf("fied time set\n");
    //���ʥץ�����ή�ι���
    
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
	/*if(){}����ʬ��brain[i].variable[0][5]+=-0.3*alpha(t-brain[j].variable[1][6])*(brain[j].variable[1][0]-V_EXI)*flag[];�Ȥ��ƤϤ����ʤ�.nazeda!*/
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
    brain[i].variable[0][0]=-34.063412;// HH��1�ѿ�V(���Ű�
    brain[i].variable[0][1]=0.906009;   // HH��2�ѿ�m
    brain[i].variable[0][2]=0.630484;   // HH��3�ѿ�h
    brain[i].variable[0][3]=1.066412;   // HH��4�ѿ�n
    brain[i].variable[0][4]=0.0;   // �����ɷ���ήI_ext
    brain[i].variable[0][5]=0.0;   // ���ʥץ�����ή�ѿ�I_syn
    brain[i].variable[0][6]=-100.0;// ȯ�л����ѿ�t_fire(����ȯ�Фαƶ��ʤ�)
  }
  
  printf("variable set\n");
  //��Ǿ�Υ˥塼�������ꡥ
  //�����ɷ���ή�����Ĥ�ѥ륹��Ω�Ĥ褦�ˤ��롥
  for(i=0; i < LEFTBRAIN; i++){
    brain[i].variable[0][4]=9.0;
    brain[i].variable[1][4]=9.0;
  }
  //��Ǿ�Υ˥塼���ΤĤʤ������ꡥ��ǾƱ�ΤϷҤ����Ƥ��ʤ���
  for(i=0;i<LEFTBRAIN;i++){
    for(j=0; j<NEURON_NUMBER;j++){
      brain[i].connectionFlag[j]=0;
    }
  }
  
  //  printf("left brain set\n");
  //��Ǿ�Υ˥塼�������ꡡ��Ǿ�Ϥ��٤ƤΥ˥塼���Ȱ���γ��ǷҤ����Ƥ��롡
  for(i=LEFTBRAIN; i<NEURON_NUMBER; i++){
    for(j=0; j<NEURON_NUMBER;j++){
      brain[i].connectionFlag[j]=judgeConnection(TH);
    }
  }
  //��Ǿ�κǸ�Υ˥塼�������ꡥ��Ǿ�ȤĤʤ��äƤ��롥
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
