#ifndef Runge_Kutta_2006_6_26_sjtu_jxd
#define Runge_Kutta_2006_6_26_sjtu_jxd

#include <cstring>
#include <loki/Typelist.h>


namespace Runge_Kutta
{

template<class base, int weishu, int jingdu >
class Varstep
{
     public:
        Varstep():sizeofsz((weishu+1)*sizeof(double)), m_jingdu(1.0/jingdu)
        {}
     protected:
        ~Varstep()
        {}
     public:
        void ok_k_fun()
        {
            P12cache();
            static_cast<base*>(this)->K_fun();
            P12ret();
            cache2P1();
            half_h();
            static_cast<base*>(this)->K_fun();
            if(!ifallmin())
              max_fact();
            else
              min_fact();
        }

    private:
        void max_fact()
        {
           do 
           {
             cache2P1();
             static_cast<base*>(this)->K_fun(); 
             P12ret(); 
             cache2P1(); 
             half_h();
             static_cast<base*>(this)->K_fun();             
           }while(!ifallmin());
        }

        void min_fact()
        {
           do
           {
             cache2P1();
             static_cast<base*>(this)->K_fun();  
             P12ret();
             cache2P1();
             double_h();
             static_cast<base*>(this)->K_fun();   
           }while(ifallmin());
           ret2P1();
           half_h();
      
        }
        
        bool ifallmin()
        {
              bool ret=true;
              for(int i=1; i<weishu+1; i++)
              {
                  if(fabs(P1_myret[i] - static_cast<base*>(this)->P1[i] ) > m_jingdu)
                  {
                     ret=false;
                     break;
                  }
              }
             return ret;
        }
 
        void P12cache()
        {
               ::std::memcpy((void*)P1_cache, (void*)(static_cast<base*>(this)->P1), sizeofsz);
        }

        void cache2P1()
        {
              ::std::memcpy((void*)(static_cast<base*>(this)->P1), (void*)P1_cache ,sizeofsz);
        }

        void P12ret()
        {
             ::std::memcpy((void*)P1_myret, (void*)(static_cast<base*>(this)->P1), sizeofsz); 
        }

        void ret2P1()
        {
            ::std::memcpy((void*)(static_cast<base*>(this)->P1), (void*)P1_myret, sizeofsz);
        }

        void half_h()
        {
            modify_h(0.5);
        }

        void double_h()
        {
            modify_h(2.0);
        }   

    private:
       double P1_cache[weishu+1];
       double P1_myret[weishu+1];
       const int sizeofsz;
       double m_jingdu;
   private:
        void modify_h(double m)
        {
            (static_cast<base*>(this)->h) *= m;
        }
    
};


template<class base, int weishu >
class Varstep<base, weishu, 0 >
{
      public:
          Varstep()
          {}
      protected:
          ~Varstep()
          {}
      public:
          void ok_k_fun()
          {
              static_cast<base*>(this)->K_fun();
          }           
};


#define BEGIN_MYFUN(funname, argument) \
            template< class MYFUNS, int i> \
            class funname\
            {\
               public:\
               void operator()(argument) \
            { 

#define END_MYFUN(funname, argument2) \
            funname<MYFUNS, i-1>()(argument2); \
             } \
            }; 

#define BEGIN_MYFUN_ZERO(funname, argument)\
            template<class MYFUNS>\
            class funname <MYFUNS, 0>\
            {\
                public:\
                void operator()(argument)\
            {

#define END_MYFUN_ZERO \
             return;\
            }\
            };

#define GETFUN(m_i, para) typename::Loki::TL::TypeAt<MYFUNS, m_i>::Result()(para)
#define k1argument double* K1, double* P1, double* P2, double h
#define k2argument double* K2, double* P1, double* P2, double* P3, double h
#define k3argument double* K3, double* P1, double* P3, double* P4, double h
#define k4argument double* K4, double* P4, double h
#define t_argument double* K1, double* K2, double* K3, double* K4, double* P1, double h
#define k1argument2  K1,  P1,  P2,  h
#define k2argument2  K2,  P1,  P2,  P3,  h
#define k3argument2  K3,  P1,  P3,  P4,  h
#define k4argument2  K4,  P4,  h
#define t_argument2  K1,  K2,  K3,  K4,  P1,  h


  BEGIN_MYFUN(myfun_k1, k1argument)
           K1[i]=GETFUN(i,P1);
           P2[i]=P1[i]+(0.5*h*K1[i-1]);
  END_MYFUN(myfun_k1, k1argument2)
  BEGIN_MYFUN_ZERO(myfun_k1, k1argument)
           K1[0]=GETFUN(0,P1);
           P2[0]=P1[0]+(0.5*h);
  END_MYFUN_ZERO

  BEGIN_MYFUN(myfun_k2, k2argument)
          K2[i]=GETFUN(i,P2);
          P3[i]=P1[i]+(0.5*h*K2[i-1]);
  END_MYFUN(myfun_k2, k2argument2)
  BEGIN_MYFUN_ZERO(myfun_k2, k2argument) 
          K2[0]=GETFUN(0,P2);
          P3[0]=P1[0]+(0.5*h);
  END_MYFUN_ZERO

  BEGIN_MYFUN(myfun_k3, k3argument)
           K3[i]=GETFUN(i,P3);
           P4[i]=P1[i]+(h*K3[i-1]);
  END_MYFUN(myfun_k3, k3argument2)
  BEGIN_MYFUN_ZERO(myfun_k3, k3argument) 
          K3[0]=GETFUN(0,P3);
          P4[0]=P1[0]+h;
  END_MYFUN_ZERO

  BEGIN_MYFUN(myfun_k4, k4argument)
          K4[i]=GETFUN(i,P4);
  END_MYFUN(myfun_k4, k4argument2)
  BEGIN_MYFUN_ZERO(myfun_k4, k4argument) 
          K4[0]=GETFUN(0,P4);
  END_MYFUN_ZERO   
      
  BEGIN_MYFUN(myfun_t, t_argument)
          P1[i]+=((K1[i-1]+2*K2[i-1]+2*K3[i-1]+K4[i-1])*h/6.0);
  END_MYFUN(myfun_t, t_argument2)
  BEGIN_MYFUN_ZERO(myfun_t, t_argument) 
          P1[0]+=h;
  END_MYFUN_ZERO   



template<class MYFUNS, int jingdu=0 >
class Kutta:
        public Varstep<Kutta<MYFUNS, jingdu >, ::Loki::TL::Length<MYFUNS>::value, jingdu>
{
    public:
        enum {weishu=::Loki::TL::Length<MYFUNS>::value};
        Kutta(const double* init_para, double step):  h(step)
        {
            for(int i=0; i<weishu+1; i++)
                P1[i]=init_para[i];
        }
        ~Kutta()
        {}
    public:
       const double* run()
       {
            static_cast<Varstep<Kutta<MYFUNS, jingdu >, weishu, jingdu>* >(this)->ok_k_fun();
          #ifdef DEBUG
             for(int i=0; i<weishu+1; i++)
                printf("Debug:   P1[%d]=%F\t", i, P1[i]);
             printf("\n");
          #endif
          return P1;
       }
       const double* run(int n)
       {
           for(int i=0; i<n ; i++)
              run();
          return P1;
       }
        

    public:
      void K_fun()
       {
           myfun_k1<MYFUNS, weishu-1>()(K1, P1, P2, h);
           P2[weishu]= P1[weishu]+(0.5*h*K1[weishu-1]);    
 
           myfun_k2<MYFUNS, weishu-1>()(K2, P1, P2, P3, h);
            P3[weishu]=P1[weishu]+(0.5*h*K2[weishu-1]);       

           myfun_k3<MYFUNS, weishu-1>()(K3, P1, P3, P4, h);
           P4[weishu]=P1[weishu]+(h*K3[weishu-1]);

           myfun_k4<MYFUNS, weishu-1>()(K4, P4, h);
 
           myfun_t<MYFUNS, weishu>()(K1, K2, K3, K4, P1, h);

          #ifdef DEBUG
             for(int i=0; i<weishu;  i++)
                 printf("K1=%F             K2=%F              K3=%F             K4=%F\n", K1[i], K2[i], K3[i], K4[i]);
          #endif         
      }  

   public:
       double P1[weishu+1];
       double h;
   private:
      double P2[weishu+1];
      double P3[weishu+1];
      double P4[weishu+1];
      double K1[weishu];
      double K2[weishu];
      double K3[weishu];
      double K4[weishu];

};


}; //---------namespace--------------end------------



#endif //Runge_Kutta_2006_6_26_sjtu_jxd
