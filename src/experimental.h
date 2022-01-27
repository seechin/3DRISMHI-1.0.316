#ifdef _EXPERIMENTAL_
//=============================================================================
//-----------------------------------------------------------------------------
//-------------------------------- new closures -------------------------------
//-----------------------------------------------------------------------------
//=============================================================================

#define CLOSURE_LR              61
#define CLOSURE_PLLR            62
#define CLOSURE_PYHNC           64
#define CLOSURE_D2HNC           65
#define CLOSURE_YBGOZ           66

#define CLOSURE_TRIAL           90


STKeywordTableUnit CLOSURE_alias_experimental[] = {
  // key names here
     { CLOSURE_LR               , "LR" },
     { CLOSURE_PLLR             , "PLLR" },
     { CLOSURE_PYHNC            , "PYHNC" },
     { CLOSURE_D2HNC            , "D2HNC" },
     { CLOSURE_YBGOZ            , "YBG" },
   // alias
     { CLOSURE_LR               , "HNC/2" },
     { CLOSURE_LR               , "linear-response" },
     { CLOSURE_PLLR             , "partial-linear-LR" },
     { CLOSURE_PLLR             , "LR2" },
     { CLOSURE_PYHNC            , "PH" },
     { CLOSURE_D2HNC            , "DH" },
     { CLOSURE_YBGOZ            , "YBG-OZ" },
     { CLOSURE_YBGOZ            , "YBGOZ" },
     { CLOSURE_YBGOZ            , "OZ-YBG" },
     { CLOSURE_YBGOZ            , "OZYBG" },
};

template <class DT> DT * append_memory(DT * src, int nsrc, DT * append, int nappend){
    DT * out = (DT*)memalloc(sizeof(DT)*(nsrc + nappend));
    memcpy(out, src, sizeof(DT)*nsrc);
    memcpy(&out[nsrc], append, sizeof(DT)*nappend);
    return out;
}

char software_version_experimental[sizeof(software_version)+10];

void init_experimental_contants(){
    #ifndef DISTRIBUTE_VERSION
        snprintf(software_version_experimental, sizeof(software_version_experimental), "%se", software_version);
        software_version = software_version_experimental;
    #endif

  // append experimental closures to the closure list
    /*int n_CLOSURE_alias_experimental = sizeof(CLOSURE_alias_experimental)/sizeof(CLOSURE_alias_experimental[0]);
    CLOSURE_alias = append_memory<STKeywordTableUnit>(CLOSURE_alias, n_CLOSURE_alias, CLOSURE_alias_experimental, n_CLOSURE_alias_experimental);
    n_CLOSURE_alias += n_CLOSURE_alias_experimental;*/


    int icnow = 0; for (icnow=1; icnow<sizeof(CLOSURE_alias)/sizeof(CLOSURE_alias[0]); icnow++) if (CLOSURE_alias[icnow].id<=0) break;
    for (int i=0; i<sizeof(CLOSURE_alias_experimental)/sizeof(CLOSURE_alias_experimental[0]); i++) if (i+icnow<sizeof(CLOSURE_alias)/sizeof(CLOSURE_alias[0])){
        CLOSURE_alias[i+icnow].id = CLOSURE_alias_experimental[i].id;
        CLOSURE_alias[i+icnow].name = CLOSURE_alias_experimental[i].name;
    }
 // register experimental closures to closure name list
    for (int i=0; i<sizeof(CLOSURE_name)/sizeof(CLOSURE_name[0]); i++) if (!CLOSURE_name[i]||!CLOSURE_name[i][0]){
        for (int j=0; j<sizeof(CLOSURE_alias_experimental)/sizeof(CLOSURE_alias_experimental[0]); j++) if (i==CLOSURE_alias_experimental[j].id){
            CLOSURE_name[i] = CLOSURE_alias_experimental[j].name; break;
        }
    }
}

double experimental_calculate_t_over_ch_by_chuv(int closure, double factor, double uuv, double huv, double cuv, double hlr){
    double chuv = huv - cuv; double t_over_ch = 0; double s, t;
    switch (closure) {
        case CLOSURE_LR             : t_over_ch = 0.5; break;
        case CLOSURE_PLLR           : t_over_ch = 0.5; break;
        case CLOSURE_PYHNC          : t_over_ch = chuv<=0? 1 : ln(1+chuv)/chuv; break;
        case CLOSURE_D2HNC          : t_over_ch = chuv<=0? 1 : 1 - chuv/2; break;
        default: t_over_ch = 1; break;
    }
    return t_over_ch;
}

//=============================================================================
//-----------------------------------------------------------------------------
//--------------------------------- revised HI --------------------------------
//-----------------------------------------------------------------------------
//=============================================================================

#define EOS_LIQUID          1   // Original LSE that ignores lambda settings
#define EOS_LIQUID_LAMBDA   2   // Using lambda to make LSE compatible with gas
#define EOS_GAS_OR_LIQUID   3   //

const char * EOS_equation_names[] = { "LSE/unknown",
    "LSE/liquid",
    "LSE/liquid-extend-to-gas",
    "LSE/liquid-or-gas"
};

class IET_Param_Exp {
  public:
    int lse_equation;
    double lse_lambda;
    //double lse_lnlambda_surplus[MAX_SOL]; int n_lse_lnlambda_surplus;
    double theta_expand;
  public:
    void init(){
        lse_equation = EOS_LIQUID_LAMBDA;
        lse_lambda = 1;
        theta_expand = 0;
    }
    double HIEquationSolver_f_exp(double x, double lse_a, double lse_b){
        switch (lse_equation) {
          case EOS_LIQUID:
            return ln(x+1e-15) + lse_b * lse_a * exp((1-1/x)/lse_a);
          case EOS_GAS_OR_LIQUID: {
            double pre = lse_b * lse_a * exp((1-1/x)/lse_a);
            double pre2 = ln(1/lse_lambda);
            return ln((x+1e-150) / lse_lambda) + (pre>pre2?pre:pre2);
          }
          default:
          case EOS_LIQUID_LAMBDA:
            return ln((x+1e-150) / lse_lambda) + lse_b * lse_a * exp((1-1/x)/lse_a);
        }
        // double pre = lse_b * lse_a * exp((1-1/x)/lse_a); return ln((x+1e-150) / lse_lambda) + pre;
    }
};

//=============================================================================
//-----------------------------------------------------------------------------
//-------------------------------- Free energy --------------------------------
//-----------------------------------------------------------------------------
//=============================================================================

double experimental_calculate_Uef1(double ulr, double clr, double hlr, double Eef02, double Uef0, double beta){
    //return hlr * ulr * beta;
    return hlr*hlr/2*beta;
}









#endif
