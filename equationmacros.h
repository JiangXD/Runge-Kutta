#ifndef EquationMacros_2007_07_25_sjtu_jxd
#define EquationMacros_2007_07_25_sjtu_jxd


#define BEGIN_EQUATION(fun_name)\
        struct fun_name\
        {\
           double operator() (const double* para)\
           {
#define END_EQUATION\
           }\
        };
     

#define WRAPPEDVAR(var_name,ini_value)\
        namespace Wrapped_Variable_Namespace_Private_\
        {\
           double var_name = ini_value;\
        };\
        inline double get_##var_name()\
        {\
           return Wrapped_Variable_Namespace_Private_::var_name;\
        }\
        inline void set_##var_name(double _v)\
        {\
           Wrapped_Variable_Namespace_Private_::var_name=_v;\
        }


#endif //EquationMacros_2007_07_25_sjtu_jxd

