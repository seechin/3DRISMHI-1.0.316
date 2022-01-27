//=============================================================================
//-----------------------------------------------------------------------------
//----------------------------- general functions -----------------------------
//-----------------------------------------------------------------------------
//=============================================================================
#ifdef _EXPERIMENTAL_
    double gas_density_1bar(double temperature){  // return: density of gas in nm^-3
        return 0.02462732363*298/temperature;  // 298K: 40.8946 mol/m^3, 1 mol/L = 0.602214086 nm^-3
    }
    double gas_dd_1bar(double temperature, double bulk_density){ return gas_density_1bar(temperature) / bulk_density; }
    double gas_dd_1bar(IET_Param * sys){ return gas_density_1bar(sys->temperature) / sys->density_hi; }
#endif

//=============================================================================
//-----------------------------------------------------------------------------
//------------------------ main-analysis-parameters.cpp -----------------------
//------------------------              and             -----------------------
//------------------------    main-preprocessing.cpp    -----------------------
//-----------------------------------------------------------------------------
//=============================================================================
const char * szHelpExperimental = "\
  Experimental features: (all optional, cautious!)\n\
    -zeta-forbid-missing    [zeta] section must define zeta of every pair\n\
      -zeta-allow-missing   [zeta] section allows missing terms\n\
    -lse_lambda, -lsl       gas/vapor density, default: -lsl 1\n\
    -theta-expand           FFT expanding of theta, default: 0\n\
    -ld-expand              radius factor for local density kernel, default: 1\n\
    -do hi,factor/ih=...    ih/is/oh/os, scale HI inside/outside hard/soft core\n\
    -do post-hi:keep/discard/complement keep/discard/fill HI density\n\
";
#ifdef _EXPERIMENTAL_
    int analysis_experinemtal_parameter_line(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line, const char * script_path){
        int ret = 0; int i = *argi; bool analysis_script = !script_name? false : (!script_name[0]? false : true);
        if (!script_name || !script_name[0]) script_name = (char*)"args\0\0";
        StringNS::string key = argv[i];

      // RISM related

      // HI related
        if (key=="-lse" || key=="--lse" || key=="lse"){
            if (i+1<argc && argv[i+1][0]!='-'){
                i++; StringNS::string key = argv[i];
                if (key=="default"){
                    sys->ex.lse_equation = EOS_LIQUID_LAMBDA;
                } else if (key=="liquid"){
                    sys->ex.lse_equation = EOS_LIQUID;
                } else if (key=="lambda" || key=="liquid-lambda" || key=="liquid_lambda" || key=="liquid-extend-to-gas" || key=="liquid_extend_to_gas"){
                    sys->ex.lse_equation = EOS_LIQUID_LAMBDA;
                } else if (key=="gas-liquid" || key=="gas_liquid" || key=="gas-or-liquid" || key=="gas_or_liquid"){
                    sys->ex.lse_equation = EOS_GAS_OR_LIQUID;
                } else if (key=="liquid-gas" || key=="liquid_gas" || key=="liquid-or-gas" || key=="liquid_or_gas"){
                    sys->ex.lse_equation = EOS_GAS_OR_LIQUID;
                } else {
                    //fprintf(sys->log(), "%s%s : %s[%d][%ld] : error : undefined LSE identifier \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+argv[i]-line+first_char_offset, argv[i], sys->is_log_tty?color_string_end:"");
                    fprintf(sys->log(), "%s%s : %s[%d] : unknown LSE identifier \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(script_name), script_line, argv[i], sys->is_log_tty?color_string_end:"");
                    ret = -1;
                }
            }
        } else if (key=="-lsl" || key=="--lsl" || key=="lsl" || key=="-lse_lambda" || key=="--lse_lambda" || key=="lse_lambda" || key=="-lse-lambda" || key=="--lse-lambda" || key=="lse-lambda" || key=="-lse_lambda" || key=="--lse_lambda" || key=="lse_lambda" || key=="-lse-lambda" || key=="--lse-lambda" || key=="lse-lambda"){
            if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ex.lse_lambda = atof(argv[++i]);
            if (sys->ex.lse_lambda < MACHINE_REASONABLE_ERROR) sys->ex.lse_lambda = MACHINE_REASONABLE_ERROR;
        /*} else if (key=="-lsls" || key=="--lsls" || key=="lsls" || key=="-lslls" || key=="--lslls" || key=="lslls" || key=="-lse-ln-lambda-surplus" || key=="--lse-ln-lambda-surplus" || key=="lse-ln-lambda-surplus" || key=="-lse_ln_lambda_surplus" || key=="--lse_ln_lambda_surplus" || key=="lse_ln_lambda_surplus" || key=="-lse-ln-lambda-surplus" || key=="--lse-ln-lambda-surplus" || key=="lse-ln-lambda-surplus" || key=="-lse_ln_lambda_surplus" || key=="--lse_ln_lambda_surplus" || key=="lse_ln_lambda_surplus"){
            if ((i+1<argc && StringNS::is_string_number(argv[i+1]))){
                i++; double v = atof(argv[i]); for (int j=0; j<MAX_SOL; j++) sys->ex.lse_lnlambda_surplus[j] = v; sys->ex.n_lse_lnlambda_surplus = MAX_SOL;
                if ((i+1<argc && StringNS::is_string_number(argv[i+1]))){
                    sys->ex.n_lse_lnlambda_surplus = 1;
                    while (i+1<argc && StringNS::is_string_number(argv[i+1])){
                        i++; sys->ex.lse_lnlambda_surplus[sys->ex.n_lse_lnlambda_surplus++] = atof(argv[i]); if (sys->ex.n_lse_lnlambda_surplus >= MAX_SOL) break;
                    }
                }
                //printf("lsls(%d):\n", sys->ex.n_lse_lnlambda_surplus); for (int i=0; i<sys->ex.n_lse_lnlambda_surplus; i++) printf("  %d : %g\n", i, sys->ex.lse_lnlambda_surplus[i]);
            }*/
        } else if (key=="-theta-expand" || key=="--theta-expand" || key=="theta-expand" || key=="-theta_expand" || key=="--theta_expand" || key=="theta_expand" || key=="-theta-expanding" || key=="--theta-expanding" || key=="theta-expanding" || key=="-theta_expanding" || key=="--theta_expanding" || key=="theta_expanding" || key=="-t-e" || key=="--t-e" || key=="t-e" || key=="-t_e" || key=="--t_e" || key=="t_e"){
            if ((i+1<argc && StringNS::is_string_number(argv[i+1]))){ i++;
                sys->ex.theta_expand = atof(argv[i]);
            }
      // no match
        } else {
            ret = 1;
        }
        *argi = i;
        return ret;
    }
    bool set_sys_options_experimental(IET_Param * sys, StringNS::string * sl, int i, bool flag){
        bool ret = false;

        return ret;
    }
#endif
#ifdef _EXPERIMENTAL_
    bool calculate_lse_ab_eperimental_automatically(IET_Param * sys, IET_arrays * arr){
      // calculate lse_b
        if (sys->calc_ab_automatically && arr->n_zeta>0){
            double lse_abn = 0; double average_density_hi = 0;
            for (int i=0; i<sys->nvm; i++){
                lse_abn += sys->llambda[i] * sys->nbulk[i];
                for (int j=0; j<sys->nvm; j++) lse_abn += - sys->nbulk[i] * sys->nbulk[j] * arr->zeta[i][j][0] * sys->density_hi * sys->density_hi;
                average_density_hi += sys->density_hi / sys->nvm;
            }
            if (sys->calc_ab_automatically==1){
                if (sys->ex.lse_equation==EOS_LIQUID_LAMBDA){
                    sys->lse_a = (lse_abn/average_density_hi + ln(sys->ex.lse_lambda))/sys->lse_b;
                } else {
                    sys->lse_a = (lse_abn/average_density_hi)/sys->lse_b;
                }
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : %s : lse_a autogened, A=%g, b=%g, lambda=%g\n", software_name, EOS_equation_names[sys->ex.lse_equation], sys->lse_a, sys->lse_b, sys->ex.lse_lambda);
            } else if (sys->calc_ab_automatically==2){
                if (sys->ex.lse_equation==EOS_LIQUID_LAMBDA){
                    sys->lse_b = (lse_abn/average_density_hi + ln(sys->ex.lse_lambda))/sys->lse_a;
                } else {
                    sys->lse_b = (lse_abn/average_density_hi)/sys->lse_a;
                }
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : %s : lse_b autogened, a=%g, B=%g, lambda=%g\n", software_name, EOS_equation_names[sys->ex.lse_equation], sys->lse_a, sys->lse_b, sys->ex.lse_lambda);
            }
            //printf("LSEANB: %12f %12f %12f\n", lse_abn, lse_abn/sys->lse_a/average_density_hi, average_density_hi);
        }
        return true;
    }
    void perform_experimental_cmd_hi_solver(IET_Param * sys, IET_arrays * arr, double * params, int nw){
        if (nw<2) return ;
        sys->lse_a = params[0]; sys->lse_b = params[1];
        if (nw>=3) sys->ex.lse_lambda = params[2];
        arr->solver_hi.set_param(sys->lse_a, sys->lse_b, &sys->ex, 0, 2, true);
    }
    void experimental_read_IETCMD_LSE(IET_Param * sys, IET_command * cmd, int ic){
        if (cmd[ic].command==IETCMD_LSE && cmd[ic].step>=1){
            double lse_a = sys->lse_a; bool lse_a_set = false;
            double lse_b = sys->lse_b; bool lse_b_set = false;
            double lse_l = sys->ex.lse_lambda; bool lse_l_set = false;
            for (int i=0; i<MAX_CMD_PARAMS; i++) if (cmd[ic].command_params_int[i]!=0){
                if (cmd[ic].command_params_int[i]==1){
                    lse_a = cmd[ic].command_params_double[i]; lse_a_set = true;
                } else if (cmd[ic].command_params_int[i]==2){
                    lse_b = cmd[ic].command_params_double[i]; lse_b_set = true;
                } else if (cmd[ic].command_params_int[i]==12){
                    lse_l = cmd[ic].command_params_double[i]; lse_l_set = true;
                }
            }
            if (lse_a_set) sys->lse_a = lse_a;
            if (lse_b_set) sys->lse_b = lse_b;
            if (lse_a_set&&!lse_b_set){
                sys->calc_ab_automatically = 2;
            } else if (!lse_a_set&&lse_b_set){
                sys->calc_ab_automatically = 1;
            }
            if (lse_l_set) sys->ex.lse_lambda = lse_l;
            //fprintf(sys->log(), "%s : solver_hi.set_param(A=%g, B=%g, lambda=%g)\n", software_name, sys->lse_a, sys->lse_b, sys->ex.lse_lambda);
        }
    }
#endif

//=============================================================================
//-----------------------------------------------------------------------------
//--------------------- main-rism.cpp : closure treatment ---------------------
//-----------------------------------------------------------------------------
//=============================================================================
#ifdef _EXPERIMENTAL_
    double inline rism_cch_from_h(double h, double nbulk_rism);
    bool experimental_closure(int closure, __REAL__ * uuv, __REAL__ * ulr, __REAL__ * huv, __REAL__ * hlr, __REAL__ * cuv, __REAL__ * clr, __REAL__ * dd, __REAL__ * res, size_t i3begin, size_t i3end, double factor, double cutoff, double nbulk_rism){
        double t;
        double h_nbulk_rism = nbulk_rism - 1;
        if (closure==CLOSURE_LR){
            for (size_t i3=i3begin; i3<i3end; i3++){
                t = exp(-uuv[i3] + 0.5*(huv[i3] - cuv[i3])) - 1;
                res[i3] = t - huv[i3]; huv[i3] = t;
            }
        } else if (closure==CLOSURE_PLLR){
            double expC = exp(cutoff);
            for (size_t i3=i3begin; i3<i3end; i3++){
                double hhere = rism_cch_from_h(huv[i3], nbulk_rism);
                t = -uuv[i3] + 0.5*(hhere - cuv[i3]);
                t = (t>cutoff? t+expC-cutoff : exp(t)) * nbulk_rism - 1;
                res[i3] = t - huv[i3]; huv[i3] = t;
            }
        } else if (closure==CLOSURE_D2HNC){
            for (size_t i3=i3begin; i3<i3end; i3++){
                double hhere = rism_cch_from_h(huv[i3], nbulk_rism);
                t = hhere - cuv[i3];
                if (t<=0){
                    t = exp(-uuv[i3] + t)*nbulk_rism - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                } else {
                    //t = exp(-uuv[i3] + (sqrt(1+4*t)-1)/2) - 1;
                    t = exp(-uuv[i3] + t - t*t/2)*nbulk_rism - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            }
        } else if (closure==CLOSURE_PYHNC){
            for (size_t i3=i3begin; i3<i3end; i3++){
                double hhere = rism_cch_from_h(huv[i3], nbulk_rism);
                t = hhere - cuv[i3];
                if (t<=0){
                    t = exp(-uuv[i3] + t)*nbulk_rism - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                } else {
                    t = exp(-uuv[i3])*(1 + t)*nbulk_rism - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            }
        } else if (closure==CLOSURE_YBGOZ){
            for (size_t i3=i3begin; i3<i3end; i3++){
                double hhere = rism_cch_from_h(huv[i3], nbulk_rism);
                double hh = (hhere<-1?-1:hhere);
                //t = exp(-uuv[i3] + factor*(hhere - cuv[i3])) * (1 - hh/2/(1+hh/2)) * nbulk_rism - 1;
                t = exp(-uuv[i3] + factor*(hhere - cuv[i3])) * (1 - hh/2/(1+hh/2)) * nbulk_rism - 1;
                res[i3] = t - huv[i3]; huv[i3] = t;
            }
        } else return false;
        return true;
    }
#endif

//=============================================================================
//-----------------------------------------------------------------------------
//--------------- main-build-ff.cpp : experimental uuv methods ----------------
//-----------------------------------------------------------------------------
//=============================================================================
#ifdef _EXPERIMENTAL_
size_t append_save_data_immediately(IET_Param * sys, IET_arrays * arr, FILE * flog, __REAL__ **** data4, int nx, int ny, int nz, int nv, FILE ** pfout, const char * _filename, const char * title=nullptr, const char * text=nullptr, double time_stamp=0, int * filter_array=nullptr, int filter_size=0){
    if (!sys || !arr || !flog || !data4 || !_filename) return 0;
    char filename[MAX_PATH]; memset(filename, 0, sizeof(filename)); strncpy(filename, _filename, sizeof(filename));
    return append_save_data(pfout, filename, flog, title, text, nx, ny, nz, nv, &data4[0][0][0][0], time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size);
}
#endif

#ifdef _EXPERIMENTAL_

    /*double ld_kernel_function(double k, int si, int sj, double kernel_params[2*MAX_SOL]){
        if (si!=sj) return 0;
        double ksigma = k*kernel_params[si];
        return k==0? 1 : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k)*kernel_params[MAX_SOL+si];
    }*/

    void perform_theta_expansion_guv_hshi(IET_Param * sys, IET_arrays * arr, int * hi_param_indicator, double * hi_param, int n_hi_param){
        double theta_expand = sys->ex.theta_expand;
        double kernel_expand = sys->ld_kernel_expand;
        for (int i=0; i<n_hi_param; i++){
            //if (hi_param_indicator[i]) printf("[%d] : %x : %g\n", i, hi_param_indicator[i], hi_param[i]);
            if (hi_param_indicator[i]==0xB34349) theta_expand = hi_param[i];  // "theta expansion"
            else if (hi_param_indicator[i]==0xB59E5B) kernel_expand = hi_param[i];  // "kernel expansion"
        }

        size_t N3 = arr->nx*arr->ny*arr->nz; double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double grid = 2*1.732050808 / Vector(arr->nx, arr->ny, arr->nz).mod();

        /*if (sys->hial==HIAL_H2HI){
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++){
                //if (arr->dd[ivm][0][0][i3]/sys->nbulk[ivm] < 0.5) arr->theta[ivm][0][0][i3] = 0; else arr->theta[ivm][0][0][i3] = 1;
                arr->theta[ivm][0][0][i3] = arr->dd[ivm][0][0][i3]/sys->nbulk[ivm];
                if (arr->theta[ivm][0][0][i3]>1) arr->theta[ivm][0][0][i3] = 1; else if (arr->theta[ivm][0][0][i3]<0) arr->theta[ivm][0][0][i3] = 0;
            }
        }*/

        if (theta_expand>0){
            if (sys->debug_level>=2){
                if (kernel_expand!=sys->ld_kernel_expand){
                    fprintf(sys->log(), "DEBUG:: post_prepare_theta(expansion=%g,kernel_expansion=%g)\n", theta_expand, kernel_expand);
                } else {
                    fprintf(sys->log(), "DEBUG:: post_prepare_theta(expansion=%g)\n", theta_expand);
                }
            }

            if (kernel_expand!=sys->ld_kernel_expand) build_ld_kernel(sys, arr->ld_kernel, sys->nvm, arr->n_gvv, theta_expand);

            __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++) cache2[ivm][0][0][i3] = arr->theta[ivm][0][0][i3] * (1+grid) - grid;
            perform_3rx1k_convolution(&arr->fftw_mp, cache2, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true); //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density_2"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache1[0][0][0][0], 1, sys); if (file_out) fclose(file_out);

          // volume expanding
            double theta_expanding_cutoff = theta_expand>1? 0 : (1-theta_expand)*(1-theta_expand)*(1+theta_expand/2);
            //if (theta_expanding_cutoff<MACHINE_REASONABLE_ERROR) theta_expanding_cutoff = MACHINE_REASONABLE_ERROR;
            for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++) if (arr->theta[iv][0][0][i3]<1){
                arr->theta[iv][0][0][i3] = cache1[iv][0][0][i3]>=theta_expanding_cutoff? 1 : 0;
            }

            if (kernel_expand!=sys->ld_kernel_expand) build_ld_kernel(sys, arr->ld_kernel, sys->nvm, arr->n_gvv, sys->ld_kernel_expand);
        }
    }

    void perform_theta_expansion_guv_eehi(IET_Param * sys, IET_arrays * arr, int * hi_param_indicator, double * hi_param, int n_hi_param){
        double theta_expand = sys->ex.theta_expand;
        double kernel_expand = sys->ld_kernel_expand;
        double dielect_hi = sys->dielect_hi;
        for (int i=0; i<n_hi_param; i++){
            if (hi_param_indicator[i]==0xB34349) theta_expand = hi_param[i];  // "theta expansion"
            else if (hi_param_indicator[i]==0xB59E5B) kernel_expand = hi_param[i];  // "kernel expansion"
            else if (hi_param_indicator[i]==0x724E6D) dielect_hi = hi_param[i];  // "dielect hi"
        }

        size_t N3 = arr->nx*arr->ny*arr->nz; double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double grid = 2*1.732050808 / Vector(arr->nx, arr->ny, arr->nz).mod();

        /*if (sys->hial == HIAL_E2HI){
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++){
                //if (arr->dd[ivm][0][0][i3]/sys->nbulk[ivm] < 0.5) arr->theta[ivm][0][0][i3] = 0; else arr->theta[ivm][0][0][i3] = 1;
                arr->theta[ivm][0][0][i3] = arr->dd[ivm][0][0][i3]/sys->nbulk[ivm];
                if (arr->theta[ivm][0][0][i3]>1) arr->theta[ivm][0][0][i3] = 1; else if (arr->theta[ivm][0][0][i3]<0) arr->theta[ivm][0][0][i3] = 0;
            }
        }*/

        if (arr->zeta && theta_expand>0){
            if (sys->debug_level>=2){
                if (kernel_expand!=sys->ld_kernel_expand){
                    fprintf(sys->log(), "DEBUG:: post_prepare_theta(expansion=%g,kernel_expansion=%g,dielect_hi=%g)\n", theta_expand, kernel_expand, dielect_hi);
                } else {
                    fprintf(sys->log(), "DEBUG:: post_prepare_theta(expansion=%g,dielect_hi=%g)\n", theta_expand, dielect_hi);
                }
            }

            if (kernel_expand!=sys->ld_kernel_expand) build_ld_kernel(sys, arr->ld_kernel, sys->nvm, arr->n_gvv, theta_expand);

            double beta = sys->default_temperature / sys->temperature;
            __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];

          // 1. compute the distribution of p•E/<ζ>
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++) if (arr->theta[ivm][0][0][i3]>=1-MACHINE_REASONABLE_ERROR){
                //double pE = beta * fabs(sys->dipole_mv[ivm] * Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]).mod()) / dielect_hi;
                double pE = sys->dipole_mv[ivm] * Vector(arr->Ecoul0[0][0][0][i3], arr->Ecoul0[1][0][0][i3], arr->Ecoul0[2][0][0][i3]).mod() / dielect_hi;
                cache2[ivm][0][0][i3] = fabs(pE / arr->zeta[ivm][ivm][0]);
            }
          // 2. expand by convolution
            perform_3rx1k_convolution(&arr->fftw_mp, cache2, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
            //save_data_immediately(sys, arr, &cache3[0][0][0][0],  arr->nx, arr->ny, arr->nz, arr->nvm, "local_density_2", "ld");

          // 3. volume expanding
            double theta_expanding_cutoff_global = theta_expand>1? 0 : (1-theta_expand)*(1-theta_expand)*(1+theta_expand/2);
            //if (theta_expanding_cutoff<MACHINE_REASONABLE_ERROR) theta_expanding_cutoff = MACHINE_REASONABLE_ERROR;
            for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++) if (arr->theta[iv][0][0][i3]<1){
                arr->theta[iv][0][0][i3] = cache1[iv][0][0][i3]>=theta_expanding_cutoff_global? 1 : 0;
            }

            if (kernel_expand!=sys->ld_kernel_expand) build_ld_kernel(sys, arr->ld_kernel, sys->nvm, arr->n_gvv, sys->ld_kernel_expand);
        }
    }

    void post_prepare_theta(IET_Param * sys, IET_arrays * arr, int * hi_param_indicator, double * hi_param, int n_hi_param){
        if (sys->hial==HIAL_HSHI){
            perform_theta_expansion_guv_hshi(sys, arr, hi_param_indicator, hi_param, n_hi_param);
        /*} else if (sys->hial==HIAL_H2HI){
                perform_theta_expansion_guv_hshi(sys, arr, hi_param_indicator, hi_param, n_hi_param);
        } else if (sys->hial==HIAL_EEHI || sys->hial==HIAL_E2HI){
            if (arr->Ecoul0){
                perform_theta_expansion_guv_eehi(sys, arr, hi_param_indicator, hi_param, n_hi_param);
            } else {
                fprintf(sys->log(), "%s : perform_theta_expansion_guv_hshi() instead of perform_theta_expansion_guv_eehi() as EF field is not generated\n", software_name);
                perform_theta_expansion_guv_hshi(sys, arr, hi_param_indicator, hi_param, n_hi_param);
            }
        */
        }

        /*if (sys->guvmal == GUVMAL_THETA && sys->ex.theta_expand>0 && sys->hial == HIAL_HSHI){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: post_prepare_theta(theta,expansion=%g)\n", sys->ex.theta_expand);
            size_t N3 = arr->nx*arr->ny*arr->nz; double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double grid = 2*1.732050808 / Vector(arr->nx, arr->ny, arr->nz).mod();
            __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];
            for (int ivm=0; ivm<sys->nvm; ivm++) for (size_t i3=0; i3<N3; i3++) cache2[ivm][0][0][i3] = arr->guvm[ivm][0][0][i3] * (1+grid) - grid;
            perform_3rx1k_convolution(&arr->fftw_mp, cache2, arr->nx, arr->ny, arr->nz, arr->box, sys->nvm, sys->nvm, arr->ld_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache1, arr->fftin, arr->fftout, arr->planf, arr->planb, true); //FILE * file_out = nullptr; char fn[MAX_PATH]; strcpy(fn, "local_density_2"); append_save_data(&file_out, fn, sys->log(), "ld", "", arr->nx, arr->ny, arr->nz, arr->nvm, &cache1[0][0][0][0], 1, sys); if (file_out) fclose(file_out);


          // volume expanding
            double theta_expanding_cutoff = sys->ex.theta_expand>1? 0 : (1-sys->ex.theta_expand)*(1-sys->ex.theta_expand)*(1+sys->ex.theta_expand/2);
            if (theta_expanding_cutoff<MACHINE_REASONABLE_ERROR) theta_expanding_cutoff = MACHINE_REASONABLE_ERROR;
            for (int iv=0; iv<sys->nvm; iv++) for (size_t i3=0; i3<N3; i3++) if (arr->guvm[iv][0][0][i3]<1){
                arr->guvm[iv][0][0][i3] = cache1[iv][0][0][i3]>=theta_expanding_cutoff? 1 : 0;
            }

        }*/
    }
#endif

#ifdef _EXPERIMENTAL_
    double experimental_calculate_solvation_free_energy(IET_Param * sys, IET_arrays * arr){
        double beta = sys->default_temperature / sys->temperature;
        double degree_of_freedoms = 3;
        return -(arr->report_total->zeta[0] + (arr->report_total->N-arr->report_total->N0)*degree_of_freedoms/2) / beta;
    }

  // new command specifier can be defined within ID of [50,90)
    #define IETCMD_POST_HI              50
    #define IETCMD_MERGE_THETA_DD       52
  // new parameter specifier can be defined within ID of [500,900)
    #define IETCMD_POST_HI_v_None       501
    #define IETCMD_POST_HI_v_Clear      502
    #define IETCMD_POST_HI_v_Complement 503

    bool experimental_analysis_command(StringNS::string key, int i_param_list, IET_command &cmd, int &cmd_type, FILE * flog, bool is_log_tty, const char * script_name, int script_line){
        bool success = true;
        if (key=="post-hi" || key=="post_hi" || key=="hi-post" || key=="hi_post"){
            cmd_type = IETCMD_POST_HI; cmd.command = IETCMD_POST_HI;
        } else if (key=="mtd" || key=="m-t-d" || key=="m_t_d" || key=="merge-theta-dd" || key=="merge_theta_dd"){
            cmd_type = IETCMD_MERGE_THETA_DD; cmd.command = IETCMD_MERGE_THETA_DD; cmd.step = 0;
        } else success = false;
        return success;
    }

    bool experimental_analysis_command_params(StringNS::string * sl, int nw, int &i, int &i_param_list, IET_command &cmd, int &cmd_type, FILE * flog, bool is_log_tty, const char * script_name, int script_line){
        bool success = true;
        if (cmd.command==IETCMD_POST_HI){
            int post_hi_method = IETCMD_POST_HI_v_None;
            if (sl[i]=="none" || sl[i]=="keep" || sl[i]=="hi"){
                post_hi_method = IETCMD_POST_HI_v_None;
            } else if (sl[i]=="remove" || sl[i]=="discard" || sl[i]=="constant" || sl[i]=="const"){
                post_hi_method = IETCMD_POST_HI_v_Clear;
            } else if (sl[i]=="complement" || sl[i]=="complementary"){
                post_hi_method = IETCMD_POST_HI_v_Complement;
            }
            if (i_param_list+1 < MAX_CMD_PARAMS){
                cmd.command_params_int[i_param_list++] = post_hi_method;
            } else {
                fprintf(flog, "%s%s : %s[%d] = post-hi : error : too many identifiers%s\n", is_log_tty?color_string_of_error:"", software_name, script_name, script_line, is_log_tty?color_string_end:""); success = false;
            }
            cmd.step = i_param_list;
        } else if (cmd.command==IETCMD_MERGE_THETA_DD){
        } else success = false;
        return success;
    }


    bool experimental_process_command(IET_Param * sys, IET_arrays * arr, IET_command * cmd, int ic){
        bool success = true; size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
        if (arr->dd && cmd[ic].command==IETCMD_POST_HI){
            if (sys->debug_level>=2){ fprintf(sys->log(), "DEBUG:: cmd[%d] = post_hi[", ic+1); for (int i=0; i<cmd[ic].step && i<MAX_CMD_PARAMS; i++) fprintf(sys->log(), i==0?"%s":",%s", cmd[ic].command_params_int[i]==IETCMD_POST_HI_v_None?"keep":cmd[ic].command_params_int[i]==IETCMD_POST_HI_v_Clear?"clear":cmd[ic].command_params_int[i]==IETCMD_POST_HI_v_Complement?"complement":"unknown"); fprintf(sys->log(), "]\n"); }
            for (int iv=0; iv<sys->nvm && iv<cmd[ic].step; iv++){
                if (cmd[ic].command_params_int[iv]==IETCMD_POST_HI_v_Clear){
                    size_t N3 = arr->nx * arr->ny * arr->nz;
                    //fprintf(sys->log(), "clear hi %d\n", ic+1);
                    for (size_t i3=0; i3<N3; i3++) arr->dd[iv][0][0][i3] = sys->nbulk[iv];
                }
            }
            for (int iv=0; iv<sys->nvm && iv<cmd[ic].step; iv++){
                if (cmd[ic].command_params_int[iv]==IETCMD_POST_HI_v_Complement){
                    size_t N3 = arr->nx * arr->ny * arr->nz;
                    //fprintf(sys->log(), "clear hi %d\n", ic+1);
                    for (size_t i3=0; i3<N3; i3++){
                        double dd_already = 0; double dd_weight = 0;
                        for (int jv=0; jv<sys->nvm; jv++){
                            if (cmd[ic].command_params_int[jv]==IETCMD_POST_HI_v_Complement){
                                dd_weight += sys->nbulk[jv];
                            } else {
                                dd_already += arr->dd[jv][0][0][i3];
                            }
                        }
                        arr->dd[iv][0][0][i3] = (1 - (dd_already>1?1:dd_already)) * sys->nbulk[iv] / dd_weight;
                    }
                }
            }
        } else if (cmd[ic].command==IETCMD_MERGE_THETA_DD){
            if (arr->dd){
                size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4m = N3 * sys->nvm;
                for (size_t i4=0; i4<N4m; i4++){
                    arr->dd[0][0][0][i4] *= arr->theta[0][0][0][i4];
                }
            }
        } else success = false;
        return success;
    }

#endif
