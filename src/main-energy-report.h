// =============================================================================
// =============================================================================
// ============ Energy reports and the calculation of free energies ============
// =============================================================================
// =============================================================================
class IET_Report {
  public:
  // fundamental parameters
    double mass, mass_mol, density;
  public:
  // basic energies
    double lj, coulsr, coullr, entropy, N, N0, dN, Ng, dNg;
    double Uef0, Uef1; // electric field energy
  // total correlations
    double cuv, clr, chuv, chlr;
  // free energy
    double excess_chem[3];  // [0]: GF, [1]: RISM, [2]: hybrid
    double zeta[4];         // [0]: molecular, [1]: atomic, [2]: molecular with closure, [3]: atomic with closure
  public:
    void operator += (IET_Report & o){
      // basic energies
        lj += o.lj; coulsr += o.coulsr; coullr += o.coullr;
        entropy += o.entropy * o.mass / o.mass_mol;
        N += o.N * o.mass / o.mass_mol; N0 += o.N0 * o.mass / o.mass_mol; dN = N - N0;
        Ng += o.Ng * o.mass / o.mass_mol; dNg = Ng - N0;
        Uef0 += o.Uef0;
        Uef1 += o.Uef1;
      // total correlations
        cuv += o.cuv;
        clr += o.clr;
        chuv += o.chuv;
        chlr += o.chlr;
      // free energy
        for (int i=0; i<sizeof(excess_chem)/sizeof(excess_chem[0]); i++) excess_chem[i] += o.excess_chem[i];
        zeta[0] += o.zeta[0] * o.mass / o.mass_mol;
        zeta[1] += o.zeta[1];
        zeta[2] += o.zeta[2] * o.mass / o.mass_mol;
        zeta[3] += o.zeta[3];
    }
    void operator *= (double scaling){
        mass *= scaling; mass_mol *= scaling; density *= scaling;
        lj *= scaling; coulsr *= scaling; coullr *= scaling; entropy *= scaling;
        N *= scaling; N0 *= scaling; Ng *= scaling; dNg *= scaling;
        cuv *= scaling; clr *= scaling; chuv *= scaling; chlr *= scaling;
        Uef0 *= scaling; Uef1 *= scaling;
        for (int i=0; i<sizeof(excess_chem)/sizeof(excess_chem[0]); i++) excess_chem[i] *= scaling;
        for (int i=0; i<sizeof(zeta)/sizeof(zeta[0]); i++) zeta[i] *= scaling;
    }
};

double calculate_zeta_by_chuv(int closure, double factor, double uuv, double huv, double cuv, double hlr){
    double chuv = huv - cuv; double t_over_ch = 0; double s, t;
    switch (closure) {
        case CLOSURE_HNC            : t_over_ch = 1; break;
        case CLOSURE_PLHNC          : t_over_ch = 1; break;
        case CLOSURE_MHNC           : t_over_ch = 1; break; // not defined
        case CLOSURE_MSA            : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_KGK            : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_PY             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : ln(1+chuv)/chuv; break;
        //case CLOSURE_D2             : t_over_ch = 1 - chuv/2; break;
        case CLOSURE_D2             : t_over_ch = 1; break;
        case CLOSURE_HNCB           : t_over_ch = 1; break;  // undefined
        //case CLOSURE_KH             : t_over_ch = -uuv+chuv<=0? 1 : chuv==0? 1 : (ln(1-uuv+chuv) + uuv)/chuv; break;
        case CLOSURE_KH             : t_over_ch = huv<=0? 1 : chuv==0? 1 : (ln(1+huv) + uuv)/chuv; break; //t_over_ch = 1; break;
        case CLOSURE_PSE2           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2)/chuv; break;
        case CLOSURE_PSE3           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6)/chuv; break;
        case CLOSURE_PSE4           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24)/chuv; break;
        case CLOSURE_PSE5           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120)/chuv; break;
        case CLOSURE_PSE6           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720)/chuv; break;
        case CLOSURE_PSE7           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040)/chuv; break;
        case CLOSURE_PSE8           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320)/chuv; break;
        case CLOSURE_PSE9           : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/362880)/chuv; break;
        case CLOSURE_PSE10          : t_over_ch = 1; break; // t_over_ch = chuv<=0? 1 : ln(1+chuv+chuv*chuv/2+chuv*chuv*chuv/6+chuv*chuv*chuv*chuv/24+chuv*chuv*chuv*chuv*chuv/120+chuv*chuv*chuv*chuv*chuv*chuv/720+chuv*chuv*chuv*chuv*chuv*chuv*chuv/5040+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/40320+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/362880+chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv*chuv/3628800)/chuv; break;
        case CLOSURE_MS             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : sqrt(1+2*chuv)/chuv; break;
        case CLOSURE_MSHNC          : t_over_ch = chuv<=0? 1 : sqrt(1+2*chuv)/chuv; break;
        case CLOSURE_BPGGHNC        : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : pow(1+factor*chuv, 1.0/factor)/chuv; break;
        case CLOSURE_VM             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : 1 - chuv/2/(1+factor*chuv); break;
        case CLOSURE_MP             : t_over_ch = 1; break; // t_over_ch = chuv==0? 1 : ((1+factor)*exp(chuv/(1+factor)) - factor)/chuv; break;
        default:
          #ifdef _EXPERIMENTAL_
            t_over_ch = experimental_calculate_t_over_ch_by_chuv(closure, factor, uuv, huv, cuv, hlr);
          #else
            t_over_ch = 1;
          #endif
            break;
    }
    return t_over_ch * cuv;
}
double calculate_Chandler_Guv(int closure, double closure_factor, double ccutoff, double uuv, double huv, double hlr, double cuv, double clr, double dd){
    double excessive_GF = (- (0.5*huv + 1) * cuv) * dd;
    double t = huv-cuv;
    if (closure<0) return excessive_GF;
    switch (closure) {
        case CLOSURE_HNC            : return excessive_GF + 0.5*huv*huv;
        case CLOSURE_PLHNC          : return excessive_GF + (huv<exp(ccutoff)? 0.5*huv*huv : 0);
        case CLOSURE_MSA            : return excessive_GF;
        case CLOSURE_KGK            : return excessive_GF;
        case CLOSURE_PY             : return excessive_GF + 0;
        case CLOSURE_D2             : return excessive_GF + 0.5*huv*huv - 0.5*t*t - huv*t*t/3;
        case CLOSURE_HNCB           : return excessive_GF + 0.5*huv*huv; // not finished
        case CLOSURE_KH             : return excessive_GF + (huv<0? 0.5*huv*huv : 0);
        case CLOSURE_PSE2           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 3 ) / 6;
        case CLOSURE_PSE3           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 4 ) / 24;
        case CLOSURE_PSE4           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 5 ) / 120;
        case CLOSURE_PSE5           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 6 ) / 720;
        case CLOSURE_PSE6           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 7 ) / 5040;
        case CLOSURE_PSE7           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 8 ) / 40320;
        case CLOSURE_PSE8           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 9 ) / 362880;
        case CLOSURE_PSE9           : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 10) / 3628800;
        case CLOSURE_PSE10          : return excessive_GF + 0.5*huv*huv - ((-uuv+t)>0?1:0) * pow((-uuv+t), 11) / 39916800;
        case CLOSURE_MHNC           :
        case CLOSURE_MS             :
        case CLOSURE_MSHNC          :
        case CLOSURE_BPGGHNC        :
        case CLOSURE_VM             :
        case CLOSURE_MP             : return excessive_GF;
    }
    return excessive_GF;
}
