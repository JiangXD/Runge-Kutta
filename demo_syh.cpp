#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "runge-kutta.h"

//-------------------------------------------------
// Name: My C++ version implemnent of syh's model
// Author: Jiang Xiao-dong
// Date: 2007-7-22   
//     For preparation of my model, I rewrote SY Huang's computional model in the paper of NeuroReport.
// The paper was very clear for the beginner to re-build the model, however the Na/Ca2+ exchanger 
// (equation 8) was confused because of inconsistant units of mM and uM. Maybe, there was also a hidden
// parameter Cell Area that was ignored or not gvien in paper. If the value of area was 0.001 cm2,or if
// the unit of equation's result mM was not converted to uM, the result could be same as that in the paper.
//     With the aid of template gereric programing technique, this program could run much more faster  
// than any others, and was easier to debug and implement.
//     Some old C++ compilers that was not complied with ISO C++ 1998 standard, such as VC 6.0 or VC 7.0,  
// could not compile this program. I build it with gcc 3.4.6 on FreeBSD. 
// 
//-------------------------------------------------


//--------------Parameter Values--------------------
namespace model_paras
{
//Ryanodine receptor
  const double v1=0.5;          
  //const double k1=0.2;
  const double k1=0.2;
//Passive Ca2+ leak
  const double v2=0.00758;
//SERCA
  const double v3=10.8;
  const double k31=0.1;
  const double k32=250;
  const double c=0.6;  //modified                    
//Ca2+ permeable AMPA receptor
  const double v4=0.000875;
  const double t4=0.1;
  const double Ca_o=2000;
  const double fai=10.0;
//Voltage-gated Ca2+ channel
  const double v5=0.0212;
  const double k5=0.3;
//Dose-response curve
  const double k6=18;
  const double t6=0.1;
//Ca2+ pump
  const double v7=5.22;               
  const double k7=0.4;
//Na+ / Ca2+ exchanger
  const double v8=0.017;
  const double Na_i=8;
  const double Na_o=120;
  const double r=0.59;
//Ca2+ buffer
  const double k9=20;
  const double Buffer_total=5;
//RT\F
  const double R_c=8314;
  const double F_c=96487;
  const double T_c=273.15+24;
};
//--------------------------------------------------


//-------double para[6] 
//-------------0. time
//-------------1. Vm;  2. m5;  3. m4;  4. Ca_ER; 5. Ca_i;
namespace model_paras
{
   enum {i_time,i_Vm,i_m5,i_m4,i_Ca_ER,i_Ca_i};
}; 

#define BEGIN_EQUATION(funname)\
class funname\
{\
    public:\
        double operator()(const double* para)\
        {

#define END_EQUATION\
        }\
 };

class Puff_AMPA
{
  public:
    static Puff_AMPA& getInstance()
    {
      static Puff_AMPA* p_AMPA=NULL;
      if(p_AMPA==NULL)
         p_AMPA=new Puff_AMPA();
      return *p_AMPA;
    }
    
    double getConcentration()
    {
      return _c_ampa;
    }

    void setConcentration(double concentration)
    {
      _c_ampa=concentration;
    }

  private:
    double _c_ampa;
    Puff_AMPA():_c_ampa(0)
    {}
    
};



//------------ equation 8 ---------
inline double fai_ex(const double* para)   
{
  using namespace model_paras;
  using std::pow;
  using std::exp;
  return (v8*pow(Na_i,3)*Ca_o/1000*exp(r*para[i_Vm]*F_c/(R_c*T_c))-v8*pow(Na_o,3)*para[i_Ca_i]/1000*exp(0-(1-r)*para[i_Vm]*F_c/(R_c*T_c)));  //----------modified by me(not x1000 to fit syh's model) ---------------
}

//------------ equation 7 ---------
inline double fai_pump(const double* para)
{
  using namespace model_paras;
  return v7*para[i_Ca_i]/(k7+para[i_Ca_i]);
}

//------------ equation 6 ---------
BEGIN_EQUATION(fun_Vm)
   using std::pow;
   using namespace model_paras;
     return ((-56.0+52*pow(Puff_AMPA::getInstance().getConcentration(),2.2)/(pow(Puff_AMPA::getInstance().getConcentration(),2.2)+pow(k6,2.2)))-para[i_Vm])/t6; 
END_EQUATION


//------------ equation 5 ---------
inline double h5 (const double* para)
{
  using namespace model_paras;
  return k5/(k5+para[i_Ca_i]);
}

inline double Beta_m5(const double* para)
{
  using namespace model_paras;
  using std::exp;
  return 3.3*exp((-65.2-para[i_Vm])/11.25);
}

inline double Alfa_m5(const double* para)
{
  using namespace model_paras;
  using std::exp;
  return 33*(92.7-para[i_Vm])/(exp((92.7-para[i_Vm])/9.6)-1);
}

BEGIN_EQUATION(fun_m5)
  using model_paras::i_m5;
  return Alfa_m5(para)*(1-para[i_m5])-Beta_m5(para)*para[i_m5];
END_EQUATION

inline double fai_VGCC(const double* para)
{
  using namespace model_paras;
  return v5*para[i_m5]*h5(para)*(Ca_o-para[i_Ca_i]);
}


//--------------- equation 4 --------
inline double M_4(const double* para)
{
  using model_paras::fai;
  return fai*Puff_AMPA::getInstance().getConcentration();
}

BEGIN_EQUATION(fun_m4)
  using namespace model_paras;
  return (M_4(para)-para[i_m4])/t4;
END_EQUATION

inline double fai_AMPA(const double* para)
{
  using namespace model_paras;
  return v4*para[i_m4]*(Ca_o-para[i_Ca_i]);
}


//-------------- equation 3 ---------
inline double fai_fill(const double* para)
{
  using namespace model_paras;
  using std::pow;
  return v3*(para[i_Ca_i]/(k31+para[i_Ca_i]))*pow(k32,6)/(pow(k32,6)+pow(para[i_Ca_ER],6));
}

inline double fai_leak(const double* para)
{
  using namespace model_paras;
  return v2*(para[i_Ca_ER]-para[i_Ca_i]);
}

inline double m1(const double* para)
{
  using namespace model_paras;
  return para[i_Ca_i]/(para[i_Ca_i]+k1);
}

inline double fai_rel(const double* para)
{
  using namespace model_paras;
  return v1*pow(m1(para),4)*(para[i_Ca_ER]-para[i_Ca_i]);
}


inline double fai_in(const double* para)
{
  return fai_AMPA(para)+fai_VGCC(para);
}

inline double fai_out(const double* para)
{
//  return fai_pump(para)+fai_ex(para);
  return fai_pump(para)-fai_ex(para);
}


//--------------- equation 2 ----------
BEGIN_EQUATION(fun_Ca_i)
  using model_paras::c;
  return (fai_rel(para)-fai_fill(para)+fai_leak(para))*c + fai_in(para) - fai_out(para);
END_EQUATION

//--------------- equation 1 ----------
BEGIN_EQUATION(fun_Ca_ER)
  return 0-(fai_rel(para)-fai_fill(para)+fai_leak(para));
END_EQUATION

//-------------- equation list --------
typedef LOKI_TYPELIST_5(fun_Vm,fun_m5,fun_m4,fun_Ca_ER,fun_Ca_i) totalequationlist;



inline double NMDA_Ca(const double* para)
{
  using model_paras::i_Vm;
  using model_paras::i_Ca_i;
  using model_paras::Ca_o;
  using model_paras::F_c;
  using std::exp;
  return -0.00000000001*10.6*Puff_AMPA::getInstance().getConcentration()*4.0*para[i_Vm]*F_c/26.5*(0.3*para[i_Ca_i]-0.3*Ca_o*1000*exp(-2.0*para[i_Vm]/26.5))/(1-exp(-2.0*para[i_Vm]/26.5));
}



int main()
{
 using Runge_Kutta::Kutta;
 using std::printf;

 double mypara[6];
 mypara[model_paras::i_time]=0;
 mypara[model_paras::i_Vm]=-56;
 mypara[model_paras::i_m4]=-0.0000000000001;
 mypara[model_paras::i_m5]=0.000631;
 mypara[model_paras::i_Ca_ER]=250.586656;
 //mypara[model_paras::i_Ca_ER]=242.077918;
 //mypara[model_paras::i_Ca_i]=0.054571;
 mypara[model_paras::i_Ca_i]=0.054841;

 Kutta<totalequationlist> mykutta(mypara,0.001);

 for(;;)
 {
 const double* ret=mykutta.run();

 if(ret[0]>1000)
    break;

 if(ret[0]>400 && ret[0] < 700)
     Puff_AMPA::getInstance().setConcentration(20);
 else
     Puff_AMPA::getInstance().setConcentration(0);

 //printf("%f  %f   %f   %f   %f   %f  %f   %f\n",ret[0],ret[1],ret[2],ret[3],ret[4],ret[5],fai_in(ret)-fai_out(ret), fai_rel(ret)-fai_fill(ret)+fai_leak(ret));
 printf("%f  %f   %f   %f   %f   %f  %f   %f\n",ret[0],ret[1],ret[2],ret[3],ret[4],ret[5],fai_AMPA(ret),fai_ex(ret));
 } 

 return 0;
}







