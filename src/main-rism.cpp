namespace RISMHI3D_RISMNS {
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>  General Closure Treatment  >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    void perform_closure(IET_Param * sys, size_t N3, size_t i3begin, size_t i3end, __REAL__ **** auuv, __REAL__ **** aulr, __REAL__ **** acuv, __REAL__ **** ahuv, __REAL__ **** aextra, __REAL__ **** ares, __REAL__ **** ahlr, __REAL__ **** aclr, __REAL__ **** add){ // extra: extra content for e.g. HNCB
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * uuv = &auuv[iv][0][0][0];
            __REAL__ * ulr = &aulr[iv][0][0][0];
            __REAL__ * cuv = &acuv[iv][0][0][0];
            __REAL__ * huv = &ahuv[iv][0][0][0];
            __REAL__ * res = &ares[iv][0][0][0];
            __REAL__ * hlr = ahlr? &ahlr[iv][0][0][0] : nullptr;
            __REAL__ * clr = aclr? &aclr[iv][0][0][0] : nullptr;
            __REAL__ * extra = aextra? &aextra[iv][0][0][0] : nullptr;
            __REAL__ * dd = add? &add[sys->nvm][0][0][0] : nullptr;
            double factor = sys->closure_factors[iv];
            double cutoff = sys->ccutoff;
            double t = 0;
            int closure = sys->closures[iv];
            //if (i3begin<0) i3begin=0; if (i3end>N3) i3end = N3;
            if (i3end>N3) i3end = N3;

            if (closure==CLOSURE_NONE){
                for (size_t i3=i3begin; i3<i3end; i3++) res[i3] = 0;
            } else if (closure==CLOSURE_HNC){   // HNC, can be scaled with -cf
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + factor*(huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_HARDSPHERE){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    if (uuv[i3]>cutoff){
                        res[i3] = -1 - huv[i3]; huv[i3] = -1;
                    } else res[i3] = 0;
                }
            } else if (closure==CLOSURE_MSA){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    if (uuv[i3]>cutoff){
                        //t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                        t = -1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        res[i3] = - uuv[i3] - cuv[i3];
                    }
                }
            } else if (closure==CLOSURE_KGK){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = t>-1? t : -1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    double t = huv[i3] - cuv[i3];
                    double E_mhnc = - t*t/(2 + 2*0.8*t);
                    t = exp(-uuv[i3] + (huv[i3] - cuv[i3]) + E_mhnc) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PY){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3])*(1 + huv[i3] - cuv[i3]) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PLHNC){
                double expC = exp(cutoff);
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = (t>cutoff? t+expC-cutoff : exp(t)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_KH){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = (t>0? t+1 : exp(t)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_D2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + t - t*t/2) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_HNCB){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    double hh = extra? extra[i3] : 0;
                    t = exp(-uuv[i3] + t - hh*hh/2) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_HARDSPHERE){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = uuv[i3]>cutoff? uuv[i3] : 0;
                    t = exp(-t + (huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }

            } else if (closure==CLOSURE_PSE2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE3){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE4){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE5){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE6){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE7){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE8){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE9){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320 + t*t*t*t*t*t*t*t*t/362880;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE10){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320 + t*t*t*t*t*t*t*t*t/362880 + t*t*t*t*t*t*t*t*t*t/3628800;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_BPGGHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    if (t<=0){
                        t = exp(-uuv[i3] + t) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        t = exp(-uuv[i3] + pow(1+factor*t, 1.0/factor)-1) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    }
                }
            } else if (closure==CLOSURE_VM){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + t - t*t/2/(1+factor*t)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MS){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + sqrt(1+2*t)-1) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MSHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    if (t<=0){
                        t = exp(-uuv[i3] + t) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        t = exp(-uuv[i3] + sqrt(1+2*t)-1) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    }
                }
            } else if (closure==CLOSURE_MP){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3])*((1+factor)*exp(t/(1+factor)) - factor) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER1){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER3){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else { // default: HNC (no scaling) or experimental
              #ifdef _EXPERIMENTAL_
                if (!experimental_closure(closure, uuv, huv, hlr, cuv, clr, dd, res, i3begin, i3end, factor, cutoff)){
              #endif
                  for (size_t i3=i3begin; i3<i3end; i3++){
                      t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                      res[i3] = t - huv[i3]; huv[i3] = t;
                  }
              #ifdef _EXPERIMENTAL_
                }
              #endif
            }

        // closure enhancement
            if (sys->closure_enhance_level==2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    //double g = 1 + huv[i3]; res[i3] /= sqrt(1 + g*g);
                    //res[i3] /= sqrt(1 + huv[i3]*huv[i3]);
                    res[i3] /= (1 + huv[i3]*huv[i3]);
                }
            } else if (sys->closure_enhance_level==1){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    res[i3] /= sqrt(1 + huv[i3]*huv[i3]);
                }
            } else if (sys->closure_enhance_level>MACHINE_REASONABLE_ERROR){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    res[i3] /= pow(1 + huv[i3]*huv[i3], sys->closure_enhance_level/2);
                }
            }
        }
    }

    void perform_closure(IET_Param * sys, IET_arrays * arr, int id){
        size_t N3 = arr->nx*arr->ny*arr->nz;
      #ifdef _LOCALPARALLEL_
        size_t idjamin = N3 / sys->nt * id; size_t idjamax = N3 / sys->nt * (id+1); if (id+1 >= sys->nt) idjamax = N3;
      #else
        size_t idjamin = 0; size_t idjamax = N3;
      #endif
        perform_closure(sys, N3, idjamin, idjamax, arr->uuv, arr->ulr, arr->cuv, arr->huv, arr->rismhi_cache[2], arr->res, arr->hlr, arr->clr, arr->dd);
    }
    void perform_closure(IET_Param * sys, IET_arrays * arr){
      #ifdef _LOCALPARALLEL_
        for (int i=1; i<sys->nt; i++) sys->mp_tasks[i] = MPTASK_RISM_CLOSURE;
        perform_closure(sys, arr, 0);
        wait_subroutines(sys);
      #else
        size_t N3 = arr->nx*arr->ny*arr->nz;
        perform_closure(sys, N3, 0, N3, arr->uuv, arr->ulr, arr->cuv, arr->huv, arr->rismhi_cache[2], arr->res, arr->hlr, arr->clr, arr->dd);
      #endif
    }


    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                            >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>     Original RISM Main     >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                            >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void perform_rism_equation(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ **** dd, int nx, int ny, int nz, int nv, Vector box, __REAL__ *** wvv, double dk_wvv, int n_wvv, __REAL__ *** nhkvv, double dk_nhkvv, int n_nhkvv, double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb, double * nhkvv_scaling){
            size_t N3 = nx * ny * nz; size_t N4 = N3 * nv;
          // 1. C(r) * w(r) -> huv
            perform_3rx1k_convolution(&arr->fftw_mp, cuv, nx, ny, nz, box, nv, nv, wvv, dk_wvv, sys->xvv_k_shift, n_wvv, huv, fftin, fftout, planf, planb, true);
          // 2. C(r)n(r) -> res, C(r)n(r) * (w(r) + nhkvv(r)) -> huv
            if (sys->use_homogeneous_rism){
                perform_3rx1k_convolution(&arr->fftw_mp, cuv, nx, ny, nz, box, nv, nv, nhkvv, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, huv, fftin, fftout, planf, planb, false);
            } else {
                for (int iv=0; iv<nv; iv++){
                    int ivm = sys->av[iv].iaa; //fprintf(stderr, "%d -(mol)-> %d : %s\n", iv, ivm, sys->av[iv].name);
                    double nhkvv_scaling_factor = nhkvv_scaling?nhkvv_scaling[iv] : 1;
                    if (dd){
                        for (size_t i3=0; i3<N3; i3++) res[iv][0][0][i3] = cuv[iv][0][0][i3] * dd[ivm][0][0][i3] * nhkvv_scaling_factor;
                    } else {
                        for (size_t i3=0; i3<N3; i3++) res[iv][0][0][i3] = cuv[iv][0][0][i3] * sys->nbulk[ivm] * nhkvv_scaling_factor;
                    }
                }
                perform_3rx1k_convolution(&arr->fftw_mp, res, nx, ny, nz, box, nv, nv, nhkvv, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, huv, fftin, fftout, planf, planb, false);
            }
            lap_timer_fftw();
    }
    void perform_rism_equation(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ *** wvv=nullptr, __REAL__ *** nhkvv=nullptr, double * nhkvv_scaling=nullptr){
        if (!wvv) wvv = arr->convolution_wvv; if (!nhkvv) nhkvv = arr->convolution_nhkvv;
        perform_rism_equation(sys, arr, cuv, huv, res, arr->dd, arr->nx, arr->ny, arr->nz, arr->nv, arr->box, wvv, arr->dk_wvv, arr->n_wvv, nhkvv, arr->dk_nhkvv, arr->n_nhkvv, arr->fftin, arr->fftout, arr->planf, arr->planb, nhkvv_scaling);
    }
    void perform_rism_equation_without_hi(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ *** wvv=nullptr, __REAL__ *** nhkvv=nullptr, double * nhkvv_scaling=nullptr){
        if (!wvv) wvv = arr->convolution_wvv; if (!nhkvv) nhkvv = arr->convolution_nhkvv;
        perform_rism_equation(sys, arr, cuv, huv, res, nullptr, arr->nx, arr->ny, arr->nz, arr->nv, arr->box, wvv, arr->dk_wvv, arr->n_wvv, nhkvv, arr->dk_nhkvv, arr->n_nhkvv, arr->fftin, arr->fftout, arr->planf, arr->planb, nhkvv_scaling);
    }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Closure preparation       >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void prepare_closure(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; int N4 = N3 * sys->nv;
        int n_prepare = 0; int prepare_cmds[MAX_SOL]; memset(prepare_cmds, 0, sizeof(prepare_cmds));

        for (int i=0; i<sys->nv; i++){
            if (sys->closures[i] == CLOSURE_HNCB) prepare_cmds[n_prepare++] = CLOSURE_HNCB;
            //else if (sys->closures[i] == CLOSURE_RCCA2) prepare_cmds[n_prepare++] = CLOSURE_RCCA2;
        }
        if (n_prepare>0){
            __REAL__ **** extra = arr->rismhi_cache[2];
            __REAL__ **** temp = arr->rismhi_cache[3];
            clear_tensor4d(extra, N4);
            for (int ip=0; ip<n_prepare; ip++){
                if (prepare_cmds[ip] == CLOSURE_HNCB){
                    perform_rism_equation(sys, arr, arr->huv, temp, arr->res, arr->convolution_wvv, arr->convolution_nhkvv, sys->renorm_lj_scaling);
                /*} else if (prepare_cmds[ip] == CLOSURE_RCCA2){
                    double beta = sys->default_temperature / sys->temperature;
                    for (int iv=0; iv<sys->nv; iv++){
                        double charge = sys->av[iv].charge;
                        double dipole = sys->dipole_mv[sys->av[iv].iaa];
                        for (size_t i3=0; i3<N3; i3++){
                            temp[iv][0][0][i3] = beta * dipole * Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]).mod();
                        }
                    }*/
                }
                for (int iv=0; iv<sys->nv; iv++) if (sys->closures[iv]==prepare_cmds[ip]){
                    cp_tensor3d(temp[iv], extra[iv], N3);
                }
            }
        }
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Unrenormalized OZ Main    >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void perform_3drism_loop_ssoz(IET_Param * sys, IET_arrays * arr){
      // 1. RISM equation
        perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res, arr->convolution_wvv, arr->convolution_nhkvv, sys->renorm_lj_scaling);
      // 2. closure: ΔH -> res
        prepare_closure(sys, arr);
        perform_closure(sys, arr);
      // x. end
        lap_timer_rism();

    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Renormalized RISM Main    >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void perform_3drism_loop_rrism(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
      // 1. hsr -> arr->huv
        perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res, arr->convolution_wvv, arr->convolution_nhkvv, sys->renorm_lj_scaling);
      // 2. h = hsr + hlr
        for (size_t i4=0; i4<N4; i4++){
            arr->huv[0][0][0][i4] += arr->hlr[0][0][0][i4];
        }
      // 3. sr closure
        prepare_closure(sys, arr);
        perform_closure(sys, arr);
        lap_timer_rism();
    }

    void perform_3drism_loop_vrism(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
      // 1. hsr -> arr->huv
        perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res, arr->convolution_wvv, arr->convolution_nhkvv, sys->renorm_lj_scaling);
      // 2. h = hsr + hlr
        for (size_t i4=0; i4<N4; i4++){
            arr->huv[0][0][0][i4] += arr->hlr[0][0][0][i4];
        }
      // 3. sr closure
        prepare_closure(sys, arr);
        perform_closure(sys, arr);
        lap_timer_rism();
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>  General Renormalization Scheme  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void prepare_3d_rism(IET_Param * sys, IET_arrays * arr){
      // long range total correlation: arr->hlr = (β U_{CoulLR+SR} q / εr) ρv * χvv
        size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
        double beta = sys->default_temperature / sys->temperature;
      // 1. arr->clr = clr = - ulr
        for (size_t i4=0; i4<N4; i4++) arr->clr[0][0][0][i4] = - arr->ulr[0][0][0][i4];
      // 2. arr->hlr = clr ρv * χvv
        if (sys->hlr_no_hi){
            perform_rism_equation_without_hi(sys, arr, arr->clr, arr->hlr, arr->res, arr->wvv_hlr, arr->nhkvv_hlr, sys->renorm_coulomb_scaling);
        } else {
            perform_rism_equation(sys, arr, arr->clr, arr->hlr, arr->res, arr->wvv_hlr, arr->nhkvv_hlr, sys->renorm_coulomb_scaling);
        }
    }
    void prepare_3d_vrism(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
        for (int iv=0; iv<sys->nv; iv++) if (sys->av[iv].reverse_rism_factor!=1){
            double reverse_rism_factor = sys->av[iv].reverse_rism_factor;
            for (size_t i3=0; i3<N3; i3++) arr->ulr[iv][0][0][i3] = reverse_rism_factor * arr->ulr[iv][0][0][i3];
        }
        prepare_3d_rism(sys, arr);
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>      Main RISM Control      >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bool perform_3drism(FILE * flog, IET_Param * sys, IET_arrays * arr, bool guess_init_n=true){
        bool success = false; double err = 1;
        double beta = sys->default_temperature / sys->temperature;
        size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = sys->nv * N3;
        bool is_out_tty = isatty(fileno(flog));
        bool multiple_closure = false; for (int iv=1; iv<sys->nv; iv++) if (sys->closures[iv] != sys->closures[0]) multiple_closure = true;

      // step 1: other preparations
        if (sys->ietal==IETAL_RRISM){
            prepare_3d_rism(sys, arr);
        } else if (sys->ietal==IETAL_VRISM){
            prepare_3d_vrism(sys, arr);
        }
      // step 2: RISM loop
        //size_t flog_pointer = 0;
        for (int istep=0; istep<sys->stepmax_rism; istep++){
            while (sys->suspend_calculation){ usleep(100); sys->is_suspend_calculation = true; continue; }
            sys->is_suspend_calculation = false;
            //if (flog && flog!=stderr && flog!=stdout) flog_pointer = ftell(flog);
          // 1. Perform RISM calculation
            if (sys->ietal==IETAL_SSOZ){            // SSOZ
                perform_3drism_loop_ssoz(sys, arr);
            } else if (sys->ietal==IETAL_VRISM){    // VRISM
                perform_3drism_loop_vrism(sys, arr);
            } else if (sys->ietal==IETAL_RRISM){    // RRISM
                perform_3drism_loop_rrism(sys, arr);
            } else {
                return false;
            }
            lap_timer_rism();
          // 2. DIIS step in
            if (sys->ndiis_rism<=1){
                for (size_t i4=0; i4<N4; i4++){
                    double res = arr->res[0][0][0][i4];
                    err += res * res;
                    arr->cuv[0][0][0][i4] += res * sys->delrism;
                }
                err = sqrt(fabs(err)/N4);
            } else {
                arr->diis_current = &arr->diis_rism;
                err = arr->diis_rism.advance(sys, &arr->cuv[0][0][0][0], &arr->res[0][0][0][0], sys->delrism, true);
            }
          // 3. report and loop control
            bool stop_loop = false;
                if (err <= sys->errtolrism){ success = true; stop_loop = true; }
                if (err < 0 || err > 1e7){ success = false; stop_loop = true; }
                if (istep+1 >= sys->stepmax_rism) stop_loop = true;
            //if (flog && (sys->detail_level>=1 || istep+1>=sys->stepmax_rism || stop_loop)){
            if (flog && (sys->detail_level>=1 || stop_loop)){
                fprintf(flog, "  %s-%s%s step %d,", IETAL_name[sys->ietal], CLOSURE_name[sys->closures[0]], multiple_closure?"...":"", istep+1);
                fprintf(flog, " stdev:"); fprintf(flog, (err<1e-3?" %.1e":" %.4f"), err); //fprintf(flog, (" %.2g"), err);
                if (sys->ndiis_rism>1) fprintf(flog, " (DIISx%d)", arr->diis_rism.ndiis);
                if (sys->debug_level>=3||sys->detail_level>=3){ char time_buffer[64]; fprintf(flog, " (%s)", get_current_time_text(time_buffer)); }
                fprintf(flog, "  %s", stop_loop? "\n" : (istep<sys->stepmax_rism&&is_out_tty&&sys->detail_level==1)?"  \r":"  \n");
                fflush(flog);
            }
            lap_timer_diis();
            fflush(flog);
            //if (!stop_loop && sys->detail_level==1 && flog && flog!=stderr && flog!=stdout) fseek(flog, flog_pointer, SEEK_SET);
            if (stop_loop) break;
        }

        return success;
    }

    bool perform_inverse_3drism(FILE * flog, IET_Param * sys, IET_arrays * arr){
        bool success = false; double err = 1;
        size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = sys->nv * N3;
        bool is_out_tty = isatty(fileno(flog));
        __REAL__ **** huv_save = arr->rismhi_cache[4];
        __REAL__ **** cuv_save = arr->rismhi_cache[5];

        /*

        __REAL__ **** cache1 = arr->rismhi_cache[2];

        for (int iv=0; iv<sys->nv; iv++){
            int ivm = sys->av[iv].iaa;
            //if (arr->dd){
            //    for (size_t i3=0; i3<N3; i3++) cache1[iv][0][0][i3] = (1+huv_save[iv][0][0][i3]) * arr->dd[ivm][0][0][i3] - 1;
            //} else {
            //    for (size_t i3=0; i3<N3; i3++) cache1[iv][0][0][i3] = (1+huv_save[iv][0][0][i3]) * sys->nbulk[ivm] - 1;
            //}
            for (size_t i3=0; i3<N3; i3++) cache1[iv][0][0][i3] = huv_save[iv][0][0][i3];
        }

        perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, sys->nv, sys->nv, arr->xvv_inverse, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, arr->cuv, arr->fftin, arr->fftout, arr->planf, arr->planb, true);

        perform_3rx1k_convolution(&arr->fftw_mp, arr->cuv, arr->nx, arr->ny, arr->nz, arr->box, sys->nv, sys->nv, arr->wvv, arr->dk_wvv, sys->xvv_k_shift, arr->n_wvv, arr->huv, arr->fftin, arr->fftout, arr->planf, arr->planb, true);

        perform_3rx1k_convolution(&arr->fftw_mp, arr->cuv, arr->nx, arr->ny, arr->nz, arr->box, sys->nv, sys->nv, arr->nhkvv, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, arr->huv, arr->fftin, arr->fftout, arr->planf, arr->planb, false);

        for (size_t i4=0; i4<N4; i4++){
            if (1+huv_save[0][0][0][i4]<MACHINE_REASONABLE_ERROR) arr->huv[0][0][0][i4] = -1;
        }

        return success;

        //*/

        //cp_tensor4d<__REAL__>(arr->huv, huv_save, N4); //cp_tensor4d<__REAL__>(arr->huv, arr->cuv, N4);
        //clear_tensor4d(arr->hlr, N4); clear_tensor4d(arr->clr, N4);
        for (size_t i4=0; i4<N4; i4++){
            arr->huv[0][0][0][i4] = arr->cuv[0][0][0][i4] = huv_save[0][0][0][i4];
        }

        for (int istep=0; istep<sys->stepmax_rism; istep++){
            while (sys->suspend_calculation){ usleep(100); sys->is_suspend_calculation = true; continue; }
            sys->is_suspend_calculation = false;
          // SCF iteration
            perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res); // no scaling scheme at here
            for (size_t i4=0; i4<N4; i4++){
                if (1+huv_save[0][0][0][i4]<MACHINE_REASONABLE_ERROR) arr->huv[0][0][0][i4] = -1;
                arr->res[0][0][0][i4] = huv_save[0][0][0][i4] - arr->huv[0][0][0][i4];

            }
          // step in
            if (sys->ndiis_rism<=1){
                for (size_t i4=0; i4<N4; i4++){
                    double res = arr->res[0][0][0][i4];
                    err += res * res;
                    arr->cuv[0][0][0][i4] += res * sys->delrism;
                }
                err = sqrt(fabs(err)/N4);
            } else {
                arr->diis_current = &arr->diis_rism;
                err = arr->diis_rism.advance(sys, &arr->cuv[0][0][0][0], &arr->res[0][0][0][0], sys->delrism, true);
            }
          // report and loop control
            bool stop_loop = false;
                if (err <= sys->errtolrism){ success = true; stop_loop = true; }
                if (err < 0 || err > 1e7){ success = false; stop_loop = true; }
                if (istep+1 >= sys->stepmax_rism) stop_loop = true;
            //if (flog && (sys->detail_level>=1 || istep+1>=sys->stepmax_rism || stop_loop)){
            if (flog && (sys->detail_level>=1 || stop_loop)){
                fprintf(flog, "  %s step %d,", IETAL_name[sys->ietal], istep+1);
                fprintf(flog, " stdev:"); fprintf(flog, (err<1e-3?" %.1e":" %.4f"), err);
                if (sys->ndiis_rism>1) fprintf(flog, " (%dx%d DIIS)", arr->diis_rism.ndiis, arr->diis_rism.ndiis);
                if (sys->debug_level>=3||sys->detail_level>=3){ char time_buffer[64]; fprintf(flog, " (%s)", get_current_time_text(time_buffer)); }
                fprintf(flog, "  %s", stop_loop? "\n" : (istep<sys->stepmax_rism&&is_out_tty&&sys->detail_level==1)?"  \r":"  \n");
                fflush(flog);
            }

        }

        return success;
    }

}

namespace RISMHI3D_RISMNS {
    bool perform_3d_iet(FILE * flog, IET_Param * sys, IET_arrays * arr, bool guess_init_n=true){
        if (sys->ietal==IETAL_SSOZ){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_3drism(algorithm=SSOZ)\n");
            return perform_3drism(flog, sys, arr, guess_init_n);
        } else if (sys->ietal==IETAL_RRISM || sys->ietal==IETAL_VRISM){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_3drism(algorithm=RISM)\n");
            return perform_3drism(flog, sys, arr, guess_init_n);
        } else if (sys->ietal==IETAL_IRISM){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_inverse_3drism()\n");
            return perform_inverse_3drism(flog, sys, arr);
        } else {
            return false;
        }
    }


}
