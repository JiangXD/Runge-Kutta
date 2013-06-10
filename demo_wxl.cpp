#include <cmath>
#include <cstdio>
#include "runge-kutta.h"
#include "equationmacros.h"

// rebuild wang xl's model in c++ by jiang xd
// date: 2007-08-19

namespace model_paras
{
  //----Endoplasmic Reticulum---
     const double fer = 0.0025;           //: ca buffering coefficient in ER from John Rinzel 1997
     const double v_rel = 10e-3;         //(/ms): 40e-3 in Keizer and Levine 1996
     const double v_fil = 0.12e-3;    //(mM/ms) : 1e-3    //(mM/ms) 
     const double k_fil = 0.3e-3;     //(mM)   : 
     const double v_lek = 0.1e-3;     //(/ms)  : 0.5e-3 in Keizer and Levine 1996 
     const double rhover = 0.15;         // : ratio of er volume to cytoplasmic volume
     const double caer0 = 0.05;       //(mM)

     const double ka_pos = 1.5e12;    //(/mM4-ms)    :ka kb and kc from --Keizer and Levine 1996
     const double ka_neg = 0.0288;    //(/ms) 
     const double kb_pos = 1.5e9;    //(/mM3-ms) 
     const double kb_neg = 0.3859;    //(/ms) 
     const double kc_pos = 0.00175;    //(/ms) 
     const double kc_neg = 1e-2;    //(/ms)

     const double FARADAY = 96485;
     const double diam = 20;      //(um)

  //---Calcium Pump---
     const double Apump = 1.2e-10;      //(/cm2/s) : Yagi 2002 1.3e-9  //-----------是否要改成ms?
     const double Kpump = .5e-3;      //(mM) : Yagi 2002 .4e-3
   
  //---Sodium-Calcium Exchanger---
     const double K = 7.47e-7;        // (mA/cm2/mM4)    :Yagi 2002 6e-8 (mA/cm2/mM4)
     const double n = 3;
     const double r = 0.59;
     const double nai = 14.1;      // mM
     const double nao = 120;       // mM
     const double cao = 2;         // mM
     const double voltage = -60;   // mV
     const double R_c=8.314;
     const double T_c=273.15+24;
 
  //---NMDA receptor---
      const double Alpha = 0.072;   //	(/ms/mM)	: forward (binding) rate
      const double Beta	= 0.0066;   // (/ms)		: backward (unbinding) rate
      const double Erev	= 0;   //	(mV)		: reversal potential
      const double gmax	= 1e-3;   //	(S/cm2)		: maximum conductance
      const double fca	= 0.75;   
 
 //---Total Calcium---
      const double fcyt = 0.01;      //: ca buffering coefficient in cytoplasm
                                     //: from John Rinzel 1997


};

WRAPPEDVAR(NMDA,0.0)

namespace model_index
{
   enum{i_time,i_PC1,i_PO2,i_PC2,i_caer,i_Ro,i_cai};
};

//--------------------------Endoplasmic Reticulum---------------- from ryrer.mod--------------

inline double get_PO1(const double* para)
{
   using namespace model_index;
   return 1.0-para[i_PC1]-para[i_PC2]-para[i_PO2];
}

BEGIN_EQUATION(fun_PC1)
   using namespace model_paras;
   using ::std::pow;
   return 0 - ka_pos * pow(para[model_index::i_cai],4.0) * para[model_index::i_PC1] + ka_neg * get_PO1(para);
END_EQUATION

BEGIN_EQUATION(fun_PO2)
   using namespace model_paras;
   using ::std::pow;
   return kb_pos * pow(para[model_index::i_cai],3.0) * get_PO1(para) - kb_neg * para[model_index::i_PO2];
END_EQUATION

BEGIN_EQUATION(fun_PC2)
   using namespace model_paras;
   return kc_pos * get_PO1(para) - kc_neg * para[model_index::i_PC2];
END_EQUATION

inline double cer_rel(const double* para)
{
   using namespace model_paras;
   using namespace model_index;
   return v_rel * (get_PO1(para) + para[i_PO2]) * (para[i_caer] - para[i_cai]);
}

inline double cer_fil(const double* para)
{
   using namespace model_paras;
   using namespace model_index;
   using ::std::pow;
   return 0 - v_fil * pow(para[i_cai],2.0) / (pow(para[i_cai],2.0) + pow(k_fil,2.0));
}

inline double cer_lek(const double* para)
{
   using namespace model_paras;
   using namespace model_index;
   return v_lek * (para[i_caer] - para[i_cai]);   
}

BEGIN_EQUATION(fun_caer)
   using namespace model_paras;
   return 0 - (cer_rel(para) + cer_lek(para) + cer_fil(para)) / (rhover / fer);
END_EQUATION

inline double get_ica_er(const double* para)
{
   using namespace model_paras;
   return 0 - 2.0 * FARADAY * (cer_rel(para) + cer_lek(para) + cer_fil(para)) * diam / 4.0 * (1e-4) / (1 + 0.15);
}


//---------------calcium pump-------------------from cap.mod--------

inline double get_ica_pump(const double* para)
{
   using namespace model_paras;
   using namespace model_index; 
   return 2.0 * FARADAY * (1e3) * (Apump * para[i_cai] / (Kpump + para[i_cai]));
}

//--------------sodium calcium exchanger---------from cae.mod-------

inline double get_v(const double* para)
{
   return model_paras::voltage;   
}

inline double para_out(const double* para)
{
   using namespace model_paras;
   using ::std::pow;
   using ::std::exp;
   return pow(nai,3.0) * cao * exp((n-2.0) * r * FARADAY * get_v(para) * (1e-3) / R_c / T_c);
}

inline double para_in(const double* para)
{
   using namespace model_paras;
   using ::std::pow;
   using ::std::exp;
   return pow(nao,3.0) * para[model_index::i_cai] * exp(0 - (n-2.0) * (1.0 -r) * FARADAY * get_v(para) * (1e-3) / R_c / T_c); 
}

inline double get_ica_exchanger(const double* para)
{
   using namespace model_paras;
   return 0 - 2.0 * K * (para_out(para) - para_in(para));
}


//--------------NMDA receptor--------------from nmda.mod----------

inline double get_Rc(const double* para)
{
   return 1.0 - para[model_index::i_Ro];
}

BEGIN_EQUATION(fun_Ro)
   using namespace model_paras;
   return Alpha * get_NMDA() * get_Rc(para) - Beta * para[model_index::i_Ro];
END_EQUATION


inline double eca(const double* para)
{
   using  namespace model_index; 
   using ::std::log;
   return 0-log(para[i_cai]/model_paras::cao)*25.0/2.0;
   //return -120;
}

inline double get_ica_nmda(const double* para)
{
   using namespace model_paras;
   return (gmax * para[model_index::i_Ro]) * fca * (get_v(para) - eca(para));
}

inline double get_ix_nmda(const double* para)
{
   using namespace model_paras;
   return (gmax * para[model_index::i_Ro]) * (1.0 - fca) * (get_v(para) - Erev); 
}


//----------------Total Calcium current-----------

inline double get_total_ca_currents(const double* para) 
{
   return get_ica_er(para) + get_ica_pump(para) + get_ica_exchanger(para) + get_ica_nmda(para);
   //return  get_ica_pump(para) + get_ica_exchanger(para) + get_ica_nmda(para);
}

//----------------Intracellular Calcium level--------

BEGIN_EQUATION(fun_cai)
   using namespace model_paras;
   return 0 - (1e4) * get_total_ca_currents(para) / 2.0 / FARADAY / diam * 4.0 / (1.0  + 0.15) * fcyt;
END_EQUATION


//   enum{i_time,i_PC1,i_PO2,i_PC2,i_caer,i_Ro,i_cai};
typedef LOKI_TYPELIST_6(fun_PC1,fun_PO2,fun_PC2,fun_caer,fun_Ro,fun_cai) myeqlist;




int main()
{
  using Runge_Kutta::Kutta;
  using namespace model_index;

  const double myinit[7] = { 0.0,           //time
                             0.999817,      //PC1
                             0.951726e-10,  //PO2
                             0.167740e-3,   //PC2
                             model_paras::caer0,  //caer
                             //5e-5,          //empty caer
                             0.0,           //Ro
                             5e-5 };        //cai

  Kutta<myeqlist> mykutta(myinit,0.001);

  for(const double* ret=myinit; ret[i_time]<5000; ret=mykutta.run(1000))
  {
     //printf("%f   %f    %f   %f    %f    %f\n",ret[i_time], ret[i_PC1], ret[i_PC2], get_PO1(ret), ret[i_PO2], ret[i_cai]);
     if(ret[i_time]>2000 && ret[i_time]<3000)
       set_NMDA(0.03);
     else 
       set_NMDA(0);
     printf("%f   %f   %f   %f    %f    %f\n",ret[i_time],ret[i_caer], ret[i_cai], cer_lek(ret), get_ica_nmda(ret), get_ix_nmda(ret));
  }
  
  return 0;
}
