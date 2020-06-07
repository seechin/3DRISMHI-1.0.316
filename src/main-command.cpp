
int process_command_sequence(int time_of_run, IET_Param * sys, IET_arrays * arr, RDF_data * rdf, RDF_data * rdfs, FILE ** pfout, int nframe, double time_stamp, int ic_begin=0){
  // run_stage: 1: initiation; -1: finalization; 0: normal
    if (ic_begin<0 || ic_begin>=sys->ncmd) return ic_begin;
    int ret = 0; bool success = true;
    bool show_run_info = sys->detail_level>=1;
    FILE * flog = sys->log();
    bool istty = sys->is_log_tty;

    IET_Param sys_bak; memcpy(&sys_bak, sys, sizeof(IET_Param));

  // temporary variables
    int i_cmd_save_filter = -1;

  // run commands
    IET_command * cmd = sys->cmd;
    for (int ic=ic_begin; ic<sys->ncmd; ic++) if (cmd[ic].command){
        ret = ic;
        if (!success) break;
      // time of running filter
        if (time_of_run>0){ // at the beginning
            if (cmd[ic].time_to_run <= 0) continue;
        } else if (time_of_run<0){ // at the end
            if (cmd[ic].time_to_run >= 0) continue;
        } else {
            if (cmd[ic].time_to_run != 0) continue;
        }

        if (sys->debug_level>=2) print_command(sys->log(), ic+1, &cmd[ic], "DEBUG:: ", "");

      // handling
//printf("nframe=%d\n", nframe); printf("sys->cmd[%d]=%d: %d, ", ic, sys->cmd[ic].command, sys->cmd[ic].step); for (int i=0; i<MAX_CMD_PARAMS; i++) if (sys->cmd[ic].command_params_int[i] || sys->cmd[ic].command_params_double[i]) printf(" (%d,%g)", sys->cmd[ic].command_params_int[i], sys->cmd[ic].command_params_double[i]); printf("\n");

        if (cmd[ic].command == IETCMD_NOP){
        } else if (cmd[ic].command == IETCMD_END){
            ret = -1; break;
        } else if (cmd[ic].command == IETCMD_DONE){
            break;
      // settings commands
        } else if (cmd[ic].command==IETCMD_CLEAR  ){
            arr->reset_for_calculation(false, true, true, false);
            arr->is_energy_calculated = false;
            arr->is_rdf_calculated = false;
            reset_rdf(sys, rdf, sys->out_rdf_bins);
        } else if (cmd[ic].command==IETCMD_RESET  ){
            arr->reset_for_calculation(true, true, true, true);
            arr->is_energy_calculated = false;
            arr->is_rdf_calculated = false;
            reset_rdf(sys, rdf, sys->out_rdf_bins);
        } else if (cmd[ic].command==IETCMD_HI_SOLVER){
            if (cmd[ic].step>=2){
                sys->lse_a = cmd[ic].command_params_double[0];
                sys->lse_b = cmd[ic].command_params_double[1];
              #ifdef _EXPERIMENTAL_
                perform_experimental_cmd_hi_solver(sys, arr, cmd[ic].command_params_double, cmd[ic].step);
              #else
                arr->solver_hi.set_param(sys->lse_a, sys->lse_b, 0, 2, true);
                fprintf(sys->log(), "solver_hi.set_param(A=%g, B=%g)\n", sys->lse_a, sys->lse_b);
              #endif
            } else {
              fprintf(sys->log(), "%s : LSE: A=%g B=%g\n", software_name, sys->lse_a, sys->lse_b);
            }
        } else if (cmd[ic].command==IETCMD_LSE){
            if (cmd[ic].step>=1){
              #ifdef _EXPERIMENTAL_
                experimental_read_IETCMD_LSE(sys, cmd, ic);
              #else
                double lse_a = sys->lse_a; bool lse_a_set = false;
                double lse_b = sys->lse_b; bool lse_b_set = false;
                for (int i=0; i<MAX_CMD_PARAMS; i++) if (cmd[ic].command_params_int[i]!=0){
                    if (cmd[ic].command_params_int[i]==1){
                        lse_a = cmd[ic].command_params_double[i]; lse_a_set = true;
                    } else if (cmd[ic].command_params_int[i]==2){
                        lse_b = cmd[ic].command_params_double[i]; lse_b_set = true;
                    }
                }
                if (lse_a_set) sys->lse_a = lse_a;
                if (lse_b_set) sys->lse_b = lse_b;
                if (lse_a_set&&!lse_b_set){
                    sys->calc_ab_automatically = 2;
                } else if (!lse_a_set&&lse_b_set){
                    sys->calc_ab_automatically = 1;
                }
              #endif
                calculate_lse_ab_automatically(sys, arr);
                #ifdef _EXPERIMENTAL_
                    arr->solver_hi.set_param(sys->lse_a, sys->lse_b, &sys->ex, 0, 2, true);
                #else
                    arr->solver_hi.set_param(sys->lse_a, sys->lse_b, 0, 2, true);
                #endif
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: solver_hi.set_param(A=%g, B=%g)\n", sys->lse_a, sys->lse_b);
            } else {
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: LSE: A=%g B=%g\n", sys->lse_a, sys->lse_b);
            }

        /*
          #ifdef _EXPERIMENTAL_
            experimental_read_IETCMD_LSE(sys, cmd, ic);
          #else
            double lse_a = sys->lse_a; bool lse_a_set = false;
            double lse_b = sys->lse_b; bool lse_b_set = false;
            for (int i=0; i<MAX_CMD_PARAMS; i++) if (cmd[ic].command_params_int[i]!=0){
                if (cmd[ic].command_params_int[i]==1){
                    lse_a = cmd[ic].command_params_double[i]; lse_a_set = true;
                } else if (cmd[ic].command_params_int[i]==2){
                    lse_b = cmd[ic].command_params_double[i]; lse_b_set = true;
                }
            }
            if (lse_a_set&&!lse_b_set){
                sys->calc_ab_automatically = 2;
            } else if (!lse_a_set&&lse_b_set){
                sys->calc_ab_automatically = 1;
            }
          #endif
            calculate_lse_ab_automatically(sys, arr);
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: LSE: A=%g B=%g\n", sys->lse_a, sys->lse_b);
            #ifdef _EXPERIMENTAL_
                arr->solver_hi.set_param(sys->lse_a, sys->lse_b, &sys->ex, 0, 2, true);
            #else
                arr->solver_hi.set_param(sys->lse_a, sys->lse_b, 0, 2, true);
            #endif
            */
        } else if (cmd[ic].command==IETCMD_SET       ){
            for (int i=0; i<cmd[ic].step&&i<MAX_CMD_PARAMS; i++){
                if (cmd[ic].command_params_int[i] == IETCMD_v_rism_dielect){
                    apply_dielect_from_dipole(sys);
                } else if (cmd[ic].command_params_int[i] == IETCMD_v_uuv){
                    sys->scale_lj = cmd[ic].command_params_double[i];
                    sys->scale_coul = cmd[ic].command_params_double[i];
                } else if (cmd[ic].command_params_int[i] == IETCMD_v_ulj){
                    sys->scale_lj = cmd[ic].command_params_double[i];
                } else if (cmd[ic].command_params_int[i] == IETCMD_v_ucoul){
                    sys->scale_coul = cmd[ic].command_params_double[i];
                }
            }
      // IO commands
        } else if (cmd[ic].command==IETCMD_LOAD){
            bool force_field_loaded[3] = { false, false, false };
            if (!szfn_in[0]){
                if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : input file (-i) not specified and nothing is loaded.\n", software_name, ic+1);
            } else for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
                char filename[MAX_PATH+32];
                size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3*arr->nv; size_t N4m = arr->nvm;
                bool load_ret = true;
                if (cmd[ic].command_params_int[is]==IETCMD_v_ulj){
                    snprintf(filename, sizeof(filename), "%s.lj", szfn_in);
                    load_ret = load_array_from_file(&arr->ulj[0][0][0][0], sys->nv, arr->nz, arr->ny, arr->nx, filename, "lj", flog);
                    if (load_ret) force_field_loaded[0] = true;
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ucoul){
                    snprintf(filename, sizeof(filename), "%s.coulomb", szfn_in);
                    load_ret = load_array_from_file(&arr->ucoulsr[0][0][0], 1, arr->nz, arr->ny, arr->nx, filename, "coulomb", flog);
                    if (load_ret) clear_tensor3d(arr->ucoullr, arr->nz, arr->ny, arr->nx);
                    if (load_ret) force_field_loaded[1] = true;
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ucoul2 && arr->Ecoul0){
                    snprintf(filename, sizeof(filename), "%s.ef", szfn_in);
                    load_ret = load_array_from_file(&arr->Ecoul0[0][0][0][0], 3, arr->nz, arr->ny, arr->nx, filename, "electrostatic field", flog);
                    if (load_ret) force_field_loaded[2] = true;
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_cuv){
                    snprintf(filename, sizeof(filename), "%s.cuv", szfn_in);
                    load_ret = load_array_from_file(&arr->cuv[0][0][0][0], sys->nv, arr->nz, arr->ny, arr->nx, filename, "cuv", flog);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_csr){
                    snprintf(filename, sizeof(filename), "%s.cuv", szfn_in);
                    load_ret = load_array_from_file(&arr->cuv[0][0][0][0], sys->nv, arr->nz, arr->ny, arr->nx, filename, "csr", flog);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_huv){
                    snprintf(filename, sizeof(filename), "%s.huv", szfn_in);
                    load_ret = load_array_from_file(&arr->huv[0][0][0][0], sys->nv, arr->nz, arr->ny, arr->nx, filename, "huv", flog);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_hlr){
                    snprintf(filename, sizeof(filename), "%s.hlr", szfn_in);
                    load_ret = load_array_from_file(&arr->hlr[0][0][0][0], sys->nv, arr->nz, arr->ny, arr->nx, filename, "huv", flog);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_dd){
                    snprintf(filename, sizeof(filename), "%s.dd", szfn_in);
                    if (!arr->dd) arr->dd = arr->__dd;
                    load_ret = load_array_from_file(&arr->dd[0][0][0][0], sys->nmv, arr->nz, arr->ny, arr->nx, filename, "dd", flog);
                }
            }
            if (force_field_loaded[0] && force_field_loaded[1] && force_field_loaded[2]){
                arr->uuv_is_ready = true; arr->frame_stamp = nframe;
            }
        } else if (cmd[ic].command==IETCMD_SAVE_FILTER){
            i_cmd_save_filter = ic;
            if (sys->debug_level>=2){
                fprintf(sys->log(), "DEBUG:: cmd[%d] = %s:", ic+1, (cmd[ic].step==1||sys->nv==1)?"pick-site":"pick-sites");
                if (cmd[ic].step<=0){
                    fprintf(sys->log(), "all\n");
                } else {
                    for (int i=0; i<cmd[ic].step&&i<MAX_CMD_PARAMS; i++){
                        if (cmd[ic].command_params_int[i]<0) fprintf(sys->log(), "all"); else fprintf(sys->log(), "%d", cmd[ic].command_params_int[i]);
                        if (i+1<cmd[ic].step&&i+1<MAX_CMD_PARAMS) fprintf(sys->log(), ",");
                    }
                    fprintf(sys->log(), "\n");
                }
            }
        //} else if (cmd[ic].command==IETCMD_SAVE_EXIST && sys->b_output_filename_autogenerated){
            // prevent saving when -o not explicitly specified
        } else if (cmd[ic].command==IETCMD_SAVE || cmd[ic].command==IETCMD_SAVE_EXIST){
            char save_iterms_buffer[128]; memset(save_iterms_buffer, 0, sizeof(save_iterms_buffer)); long int total_size = 0; long int original_size = 0;
            for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
                const char * save_item = ""; switch (cmd[ic].command_params_int[is]) {
                    case IETCMD_v_cmd: save_item = "cmd";   break;
                    case IETCMD_v_ulj: save_item = "lj";   break;
                    case IETCMD_v_ucoul: save_item = "coul";  break;
                    case IETCMD_v_ucoul2: save_item = "ef";   break;
                    case IETCMD_v_rmin: save_item = "rmin";   break;
                    case IETCMD_v_uuv: save_item = "uuv";     break;
                    case IETCMD_v_ulr: save_item = "ulr";     break;
                    case IETCMD_v_cuv: save_item = "cuv";     break;
                    case IETCMD_v_clr: save_item = "clr";     break;
                    case IETCMD_v_huv: save_item = "huv";     break;
                    case IETCMD_v_hlr: save_item = "hlr";     break;
                    case IETCMD_v_dd:  save_item = "dd";      break;
                    case IETCMD_v_ddp: save_item = "ddp";     break;
                    case IETCMD_v_theta: save_item = "theta"; break;
                    case -IETCMD_v_dd: save_item = "nphi";    break;
                    case IETCMD_v_guv: save_item = "guv";     break;
                    case IETCMD_v_ld : save_item = "ld";      break;
                }
                if (save_item[0]) snprintf(&save_iterms_buffer[strlen(save_iterms_buffer)], sizeof(save_iterms_buffer)-strlen(save_iterms_buffer), "%s%s", (!save_iterms_buffer[0])?"":",", save_item);
            };
          lap_timer_others();
            if (!szfn_out[0]){
                if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : output file (-o) not specified and nothing is saved.\n", software_name, ic+1);
            } else for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
                char filename[MAX_PATH+32];
                char save_info_text[40960]; print_save_extra_info(sys, save_info_text, sizeof(save_info_text));
                size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3*arr->nv; size_t N4m = arr->nvm;
                size_t write_size = 0; size_t write_original_size = 0; int * filter_array = i_cmd_save_filter>0? cmd[i_cmd_save_filter].command_params_int : nullptr; int filter_size = i_cmd_save_filter>0? cmd[i_cmd_save_filter].step : 0;
                if (cmd[ic].command_params_int[is]==IETCMD_v_cmd){
                    size_t len = strnlen(save_info_text, sizeof(save_info_text));
                    snprintf(&save_info_text[len], sizeof(save_info_text)-len, "\n%s@%s:%s$", username, hostname, szfn_path);
                    for (int i=0; i<sys->argc; i++){
                        size_t len = strnlen(save_info_text, sizeof(save_info_text));
                        snprintf(&save_info_text[len], sizeof(save_info_text)-len, " %s", sys->argv[i]);
                    }
                    //printf("\33[31m[%s]\n\33[0m", save_info_text);
                    arr->res[0][0][0][0] = 0;
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "cmd", save_info_text, 1, 1, 1, 1, arr->res, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ulj){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "lj", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->ulj, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "lj", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->ulj[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //snprintf(filename, sizeof(filename), "%s.lj", szfn_out);
                    //compress_data(filename, "LJ", flog, arr->nv, arr->nx, arr->ny, arr->nz, &arr->ulj[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ucoul){
                    for (size_t i3=0; i3<N3; i3++) arr->res[0][0][0][i3] = arr->ucoulsr[0][0][i3] + arr->ucoullr[0][0][i3];
                    total_size += append_save_data(pfout, szfn_out, flog, "coul", save_info_text, arr->nx, arr->ny, arr->nz, 1, &arr->res[0][0][0][0], time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size);
                    original_size += 1*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //snprintf(filename, sizeof(filename), "%s.coulomb", szfn_out);
                    //compress_data(filename, "Coulomb potential", flog, 1, arr->nx, arr->ny, arr->nz, &arr->res[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ucoul2 && arr->Ecoul0){
                    total_size += append_save_data(pfout, szfn_out, flog, "ef", save_info_text, arr->nx, arr->ny, arr->nz, 3, &arr->Ecoul0[0][0][0][0], time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size);
                    original_size += 3*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //snprintf(filename, sizeof(filename), "%s.ef", szfn_out);
                    //compress_data(filename, "Electrostatic field", flog, 3, arr->nx, arr->ny, arr->nz, &arr->Ecoul0[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_uuv){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "uuv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->uuv, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "uuv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->uuv[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ulr){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "ulr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->ulr, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "ulr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->ulr[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_cuv){
                    for (size_t i4=0; i4<N4; i4++) arr->res[0][0][0][i4] = arr->cuv[0][0][0][i4] + arr->clr[0][0][0][i4];
                    total_size += select_append_save_data(filter_array, filter_size, arr->rismhi_cache[2], &write_size, &write_original_size, pfout, szfn_out, flog, "cuv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->res, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_csr){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "csr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->cuv, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "cuv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->cuv[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //snprintf(filename, sizeof(filename), "%s.cuv", szfn_out);
                    //compress_data(filename, "cuv", flog, arr->nv, arr->nx, arr->ny, arr->nz, &arr->cuv[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_clr){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "clr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->clr, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "clr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->clr[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_huv){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "huv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->huv, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "huv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->huv[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //snprintf(filename, sizeof(filename), "%s.huv", szfn_out);
                    //compress_data(filename, "huv", flog, arr->nv, arr->nx, arr->ny, arr->nz, &arr->huv[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_hlr){
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "hlr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->hlr, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "hlr", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->hlr[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_dd){
                    if (arr->dd){
                        total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "dd", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, arr->dd, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                        //total_size += append_save_data(pfout, szfn_out, flog, "dd", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, &arr->dd[0][0][0][0], time_stamp, sys);
                        //original_size += arr->nvm*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                        //snprintf(filename, sizeof(filename), "%s.dd", szfn_out);
                        //compress_data(filename, "dd", flog, arr->nvm, arr->nx, arr->ny, arr->nz, &arr->dd[0][0][0][0], sys->err_output_uv);
                    } else if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : HI not performed and dd will not be saved.\n", software_name, ic+1);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ddp){
                    if (arr->dd){
                        total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "ddp", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, arr->ddpot_hi, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                        //total_size += append_save_data(pfout, szfn_out, flog, "dd", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, &arr->dd[0][0][0][0], time_stamp, sys);
                        //original_size += arr->nvm*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                        //snprintf(filename, sizeof(filename), "%s.dd", szfn_out);
                        //compress_data(filename, "dd", flog, arr->nvm, arr->nx, arr->ny, arr->nz, &arr->dd[0][0][0][0], sys->err_output_uv);
                    } else if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : HI not performed and ddp will not be saved.\n", software_name, ic+1);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_theta){
                    if (arr->dd){
                        total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "ddp", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, arr->theta, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    } else if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : HI not performed and theta will not be saved.\n", software_name, ic+1);
                } else if (cmd[ic].command_params_int[is]==-IETCMD_v_dd){
                    if (arr->nphi){
                        total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "nphi", save_info_text, arr->nx, arr->ny, arr->nz, arr->nvm, arr->nphi, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    } else if (show_run_info) fprintf(flog, "%s : cmd[%d] : waning : HI not performed and nphi will not be saved.\n", software_name, ic+1);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_guv){
                    //snprintf(filename, sizeof(filename), "%s.guv", szfn_out);
                    double beta = sys->default_temperature / sys->temperature;
                    for (int iv=0; iv<sys->nv; iv++){
                        int ivm = sys->av[iv].iaa; //fprintf(sys->log(), "%d -(mol)-> %d : %s\n", iv, ivm, sys->av[iv].name);
                        for (size_t i3=0; i3<N3; i3++){
                            arr->res[iv][0][0][i3] = (1+arr->huv[iv][0][0][i3]) * (arr->dd? arr->dd[ivm][0][0][i3] : sys->nbulk[ivm]);
                            if (beta*arr->ulj[iv][0][0][i3]<sys->ucutoff_hs && arr->res[iv][0][0][i3]<MACHINE_REASONABLE_ERROR) arr->res[iv][0][0][i3] = MACHINE_REASONABLE_ERROR; // non hardsphere correction
                        }
                    }
                    total_size += select_append_save_data(filter_array, filter_size, arr->rismhi_cache[2], &write_size, &write_original_size, pfout, szfn_out, flog, "guv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->res, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "guv", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, &arr->res[0][0][0][0], time_stamp, sys);
                    //original_size += arr->nv*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                    //compress_data(filename, "density profile", flog, arr->nv, arr->nx, arr->ny, arr->nz, &arr->res[0][0][0][0], sys->err_output_uv);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_rmin && arr->r2uvmin){
                    double rc = sys->rvdw; if (rc>sys->rcoul) rc = sys->rcoul; double rc2 = rc*rc;
                    for (size_t i3=0; i3<N3; i3++){
                        if (arr->r2uvmin[0][0][i3]<rc2 && arr->r2uvmin[0][0][i3]>=0){
                            arr->res[0][0][0][i3] = sqrt(arr->r2uvmin[0][0][i3]);
                        } else {
                            arr->res[0][0][0][i3] = -1;
                        }
                    }
                    total_size += select_append_save_data(filter_array, filter_size, arr->rismhi_cache[2], &write_size, &write_original_size, pfout, szfn_out, flog, "rmin", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->res, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                    //total_size += append_save_data(pfout, szfn_out, flog, "rmin", save_info_text, arr->nx, arr->ny, arr->nz, 1, &arr->res[0][0][0][0], time_stamp, sys);
                    //original_size += 1*N3*sizeof(__REAL__) + sizeof(CompressPageHeader);
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){ char filename_without_ext[MAX_PATH]; file_without_extension(szfn_rdf, filename_without_ext);
                    double rcutoff = sys->rvdw>sys->rcoul?sys->rvdw:sys->rcoul;
                    snprintf(filename, sizeof(filename), "%s.rdf", filename_without_ext);
                    bool allow_to_save = true;
                    if (sys->output_override==0 && access(filename, F_OK) != -1){ // file already exist and not appendable
                        if (!get_nonoverwrite_filename(filename, flog)) allow_to_save = false;
                    } else if (sys->output_override==-1 && access(filename, F_OK) != -1){
                        fprintf(flog, "%s : overide output file %s\n", software_name, filename);
                    }
                    if (allow_to_save){
                        if (!arr->is_rdf_calculated){ calculate_rdf(sys, arr, rdf, rcutoff, sys->out_rdf_bins); arr->is_rdf_calculated = true; }
                        if (cmd[ic].time_to_run>=0){
                            print_rdf(filename, sys, arr, rdf, rcutoff, sys->out_rdf_bins);
                        } else {
                            print_rdf(filename, sys, arr, rdfs, rcutoff, sys->out_rdf_bins);
                        }
                    }
                } else if (cmd[ic].command_params_int[is]==IETCMD_v_ld){
                    arr->calculate_local_density_from_guv(sys, arr->res, arr->rismhi_cache[2]);
                    total_size += select_append_save_data(filter_array, filter_size, arr->res, &write_size, &write_original_size, pfout, szfn_out, flog, "rmin", save_info_text, arr->nx, arr->ny, arr->nz, arr->nv, arr->res, time_stamp, sys, arr->compress_buffer, arr->compress_buffer_size); original_size += write_original_size;
                }
            }
          lap_timer_io();
            //if ((sys->debug_level>0 || sys->detail_level>0) && total_size>0){
            //if (total_size>0){
            if (save_iterms_buffer[0]){
                char save_items_mem_buffer[2][64];
                fprintf(flog, "  save %s to %s%s%s, totally %s", save_iterms_buffer, sys->is_log_tty?prompt_path_prefix:"\"", szfn_out, sys->is_log_tty?prompt_path_suffix:"\"", print_memory_value(save_items_mem_buffer[0], sizeof(save_items_mem_buffer[0]), total_size));
                if (original_size>total_size) fprintf(flog, " compressed from %s", print_memory_value(save_items_mem_buffer[1], sizeof(save_items_mem_buffer[1]), original_size));
                fprintf(flog, "\n");
            }
      // display results
        } else if (cmd[ic].command==IETCMD_DISPLAY || cmd[ic].command==IETCMD_REPORT){
          // calculate energy
            //if (!arr->is_energy_calculated){ recalculate_energy(sys, arr); arr->is_energy_calculated = true; }
            recalculate_energy(sys, arr); arr->is_energy_calculated = true;
          // display energy
            if (cmd[ic].command==IETCMD_DISPLAY){
                display_IET_report(sys, arr, sys->cmd, ic);
            } else {
                print_IET_report(sys, arr, sys->cmd, ic);
            }
          // display other terms, e.g. rdf
            double rcutoff = sys->rvdw>sys->rcoul?sys->rvdw:sys->rcoul;
            bool need_print_rdf = false;
            for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
                if (cmd[ic].command_params_int[is]==IETCMD_v_rdf) need_print_rdf = true;
            }
            if (need_print_rdf){
                if (!arr->is_rdf_calculated){ calculate_rdf(sys, arr, rdf, rcutoff, sys->out_rdf_bins); arr->is_rdf_calculated = true; }
                if (cmd[ic].time_to_run>=0){
                    print_rdf(flog, sys, arr, rdf, rcutoff, sys->out_rdf_bins);
                } else {
                    print_rdf(flog, sys, arr, rdfs, rcutoff, sys->out_rdf_bins);
                }
            }
      // scaling values for arrays
        } else if (cmd[ic].command==IETCMD_SCALE     ){
            for (int i=0; i<MAX_CMD_PARAMS; i++){
                double scale = cmd[ic].command_params_double[i];
                size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3*arr->nv; size_t N4m = arr->nvm;
                if (cmd[ic].command_params_int[i]==IETCMD_v_uuv){
                    for (size_t i4=0; i4<N4; i4++){ arr->ulj[0][0][0][i4] *= scale; }
                    for (size_t i3=0; i3<N3; i3++){ arr->ucoulsr[0][0][i3] *= scale; arr->ucoullr[0][0][i3] *= scale; }
                } else if (cmd[ic].command_params_int[i]==IETCMD_v_ulj){
                    for (size_t i4=0; i4<N4; i4++){ arr->ulj[0][0][0][i4] *= scale; }
                } else if (cmd[ic].command_params_int[i]==IETCMD_v_ucoul){
                    for (size_t i3=0; i3<N3; i3++){ arr->ucoulsr[0][0][i3] *= scale; arr->ucoullr[0][0][i3] *= scale; }
                }
            }
      // setting values for arrays
        } else if (cmd[ic].command==IETCMD_CLOSURE   ){
            for (int i=0; i<sys->nav; i++){
                int iaa = sys->av[i].iaa;
                if (iaa>=0 && iaa<cmd[ic].step && iaa<MAX_CMD_PARAMS) sys->closures[i] = cmd[ic].command_params_int[iaa];
                else if (cmd[ic].step<=1) sys->closures[i] = cmd[ic].command_params_int[0];
            }
        } else if (cmd[ic].command==IETCMD_CLOSURE_A ){
            for (int i=0; i<sys->nav; i++){
                if (i>=0 && i<cmd[ic].step && i<MAX_CMD_PARAMS) sys->closures[i] = cmd[ic].command_params_int[i];
                else if (cmd[ic].step<=1) sys->closures[i] = cmd[ic].command_params_int[0];
            }
        } else if (cmd[ic].command==IETCMD_TEST || cmd[ic].command==IETCMD_TEST_SAVE){
            perform_sub_test(sys, arr, &cmd[ic]);
        } else if (cmd[ic].command==IETCMD_CF        ){
            for (int i=0; i<sys->nav; i++){
                int iaa = sys->av[i].iaa;
                if (iaa>=0 && iaa<cmd[ic].step && iaa<MAX_CMD_PARAMS) sys->closure_factors[i] = cmd[ic].command_params_double[iaa];
                else if (cmd[ic].step<=1) sys->closure_factors[i] = cmd[ic].command_params_double[0];
            }
        } else if (cmd[ic].command==IETCMD_CF_A      ){
            for (int i=0; i<sys->nav; i++){
                if (i>=0 && i<cmd[ic].step && i<MAX_CMD_PARAMS) sys->closure_factors[i] = cmd[ic].command_params_double[i];
                else if (cmd[ic].step<=1) sys->closure_factors[i] = cmd[ic].command_params_double[0];
            }
        } else if (cmd[ic].command==IETCMD_DIELECT   ){
            for (int i=0; i<sys->nav; i++){
                int iaa = sys->av[i].iaa;
                if (iaa>=0 && iaa<cmd[ic].step && iaa<MAX_CMD_PARAMS) sys->dielect[i] = cmd[ic].command_params_double[iaa];
            }
            for (int i=0; i<cmd[ic].step && i<MAX_CMD_PARAMS && i<MAX_SOL; i++) sys->dielect_mol[i] = cmd[ic].command_params_double[i];
            double dielecti = 0; double density_dielect = 0; for (int ivm=0; ivm<sys->nvm; ivm++){ dielecti += sys->density_mv[ivm] / sys->dielect_mol[ivm]; density_dielect += sys->density_mv[ivm]; } sys->mean_dielect = density_dielect / dielecti;
            //if (sys->esal==CoulAL_Coulomb) sys->mean_dielect = 1;
        } else if (cmd[ic].command==IETCMD_DENSITY   ){
            for (int i=0; i<sys->nav; i++){
                int iaa = sys->av[i].iaa;
                if (iaa>=0 && iaa<cmd[ic].step && iaa<MAX_CMD_PARAMS) sys->density_av[i] = cmd[ic].command_params_double[iaa];
            }
            for (int i=0; i<cmd[ic].step && i<MAX_CMD_PARAMS && i<MAX_SOL; i++) sys->density_mv[i] = cmd[ic].command_params_double[i];
      // rebuild force field
        } else if (cmd[ic].command==-IETCMD_BUILD_FF){
            arr->frame_stamp = nframe;
            arr->uuv_is_ready = true;
            build_uuv_base_on_force_field(sys, arr);
        } else if (cmd[ic].command==IETCMD_BUILD_FF){
            arr->uuv_is_ready = false;
            build_force_field_auto(sys, arr, flog, nframe);
            build_uuv_base_on_force_field(sys, arr);
        } else if (cmd[ic].command==IETCMD_BUILD_UUV){
            int esal_save = sys->esal;
            double mean_dielect_save = sys->mean_dielect;
            double yukawa_dielect_save = sys->dielect_yukawa;
            double yukawa_rc_save = sys->rc_yukawafft;
            sys->esal = cmd[ic].command_params_int[0];
            sys->mean_dielect = sys->dielect_yukawa = cmd[ic].command_params_double[0];
            sys->rc_yukawafft = cmd[ic].command_params_double[1]<=0? 1/MACHINE_REASONABLE_ERROR : cmd[ic].command_params_double[1];
            build_uuv_base_on_force_field(sys, arr);
            sys->esal = esal_save;
            sys->mean_dielect = mean_dielect_save;
            sys->dielect_yukawa = yukawa_dielect_save;
            sys->rc_yukawafft = yukawa_rc_save;
        } else if (cmd[ic].command==IETCMD_HOLD){
            if (cmd[ic].step <=0){
                sys->_n_hold_list = 0;
            } else {
                sys->_n_hold_list = cmd[ic].step;
                for (int i=0; i<sys->nv && i<cmd[ic].step; i++) sys->_hold_list[i] = cmd[ic].command_params_int[i];
            }
      // running HI
        } else if (cmd[ic].command>=2000 && cmd[ic].command<2500){  // HI
            int hial = cmd[ic].command - 2000;

            build_force_field_auto(sys, arr, flog, nframe);
            build_uuv_base_on_force_field(sys, arr);

            sys->hial = hial;
            sys->stepmax_hi = cmd[ic].step;
            if (!arr->dd) arr->dd = arr->__dd;

            if (!arr->zeta && hial>=HIAL_HI){
                fprintf(sys->log(), "%s%s : error : cannot perform %s: -zeta or [zeta] not defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, HIAL_name[hial], sys->is_log_tty?color_string_end:"");
                success = false;
            } else {
                if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_hi()\n");
                bool hi_converged = perform_hi(sys, arr, cmd[ic].command_params_int, cmd[ic].command_params_double, MAX_CMD_PARAMS);
                //if (sys->detail_level>=1&&!hi_converged) fprintf(flog, "  %s not converged\n", HIAL_name[sys->hial]);
                if (arr->dd&&(sys->debug_level>=3||sys->debug_show_crc)) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(dd)         = %08X\n", check_real_crc(&arr->dd[0][0][0][0], arr->nx*arr->ny*arr->nz));
            }

            /*if (!arr->zeta && hial>=HIAL_HI){
                fprintf(sys->log(), "%s%s : error : cannot perform %s: -zeta or [zeta] not defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, HIAL_name[hial], sys->is_log_tty?color_string_end:"");
                success = false;
            } else {
                build_force_field_auto(sys, arr, flog, nframe);
                build_uuv_base_on_force_field(sys, arr);

                sys->hial = hial;
                sys->stepmax_hi = cmd[ic].step;
                if (!arr->dd) arr->dd = arr->__dd;
                if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_hi()\n");
                bool hi_converged = perform_hi(sys, arr, cmd[ic].command_params_int, cmd[ic].command_params_double, MAX_CMD_PARAMS);
                //if (sys->detail_level>=1&&!hi_converged) fprintf(flog, "  %s not converged\n", HIAL_name[sys->hial]);
                if (arr->dd&&(sys->debug_level>=3||sys->debug_show_crc)) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(dd)         = %08X\n", check_real_crc(&arr->dd[0][0][0][0], arr->nx*arr->ny*arr->nz));
            }*/
      // running SSOZ/RISM
        } else if (cmd[ic].command==1500){ //RISM_NONE
            build_force_field_auto(sys, arr, flog, nframe);
            build_uuv_base_on_force_field(sys, arr);
            size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3 * sys->nv;
            for (size_t i4=0; i4<N4; i4++){
                arr->huv[0][0][0][i4] = arr->uuv[0][0][0][i4]>sys->ucutoff_hs? -1 : 0;
                arr->cuv[0][0][0][i4] = arr->uuv[0][0][0][i4]>sys->ucutoff_hs? -sys->ucutoff_hs : -arr->uuv[0][0][0][i4];
            }
        } else if (cmd[ic].command>1500 && cmd[ic].command<2000){  // RISM
            int ietal = cmd[ic].command - 1500; // SSOZ, RISM, etc

            //clear_tensor4d(arr->huv, N4); clear_tensor4d(arr->cuv, N4);
                // this is crucially important for some closures.
                // As long as RISM is SCF method, a smart guess of initial field does little help

            if (!arr->wvv || !arr->nhkvv){
                if (sys->gvv_specification==0) fprintf(sys->log(), "%s%s : error : cannot perform %s: -gvv not specified%s\n", sys->is_log_tty?color_string_of_error:"", software_name, IETAL_name[ietal], sys->is_log_tty?color_string_end:""); else fprintf(sys->log(), "%s%s : error : cannot perform %s: cannot use %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, IETAL_name[ietal], sys->gvv_specification>0?"-gvv":"-hvv", sys->is_log_tty?color_string_end:"");
                success = false;
            } else {
                sys->ietal = ietal;
                sys->stepmax_rism = cmd[ic].step;
                if (cmd[ic].command_params_int[3]>0) for (int i=0; i<MAX_SOL; i++) sys->closures[i] = cmd[ic].command_params_int[3];
                if (cmd[ic].command_params_int[4]>0) sys->esal = cmd[ic].command_params_int[4]-1;
                if (cmd[ic].command_params_int[5]>0) apply_dielect_from_dipole(sys);

                build_force_field_auto(sys, arr, flog, nframe);
                build_uuv_base_on_force_field(sys, arr);

                //if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_3d_iet()\n");
                arr->diis_rism.ndiis = 0;
                bool rism_converged = RISMHI3D_RISMNS::perform_3d_iet(flog, sys, arr, true);
                //if (sys->detail_level>=1&&!rism_converged) fprintf(flog, "  %s not converged\n", IETAL_name[sys->ietal]);
                if (sys->debug_level>=3||sys->debug_show_crc) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(huv)        = %08X\n", check_real_crc(&arr->huv[0][0][0][0], sys->nv*arr->nx*arr->ny*arr->nz));
            }

      // block commands
        } else if (cmd[ic].command==IETCMD_TI){
            if (cmd[ic].step<2){
                fprintf(flog, " %s: warning : will not perform TI because too few steps (%d)\n", software_name, cmd[ic].step);
            } else {
                size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3*arr->nv; size_t N4m = arr->nvm;
                IET_Report TI_energy; memset(&TI_energy, 0, sizeof(IET_Report));
              // 1. generate forcefield and save to safe buffer
                build_force_field_auto(sys, arr, flog, nframe);
              // 2-6. TI windows
                for (int iti=0; iti<cmd[ic].step; iti++){ double scale = (iti+1.0) / cmd[ic].step;
                    arr->is_energy_calculated = false;
                    if (sys->detail_level>=1){
                        if (flog) fprintf(flog, " TI step %d of %d, scaling=%g:\n", iti+1, cmd[ic].step, scale);
                    } else {
                        if (flog && istty && iti==0){ char buffer[53]; buffer[0] = '['; buffer[51] = ']'; buffer[52] = 0; for (int i=1; i<=50; i++) buffer[i] = ' '; fprintf(flog, " TI progress %s%s", buffer, iti+1>=cmd[ic].step?"\n":"\r");
                        }
                    }
                  // 2. scale potential
                    if (scale!=1){
                        for (size_t i4=0; i4<N4; i4++){ arr->ulj[0][0][0][i4] *= scale; }
                        for (size_t i3=0; i3<N3; i3++){
                            arr->ucoulsr[0][0][i3] *= scale; arr->ucoullr[0][0][i3] *= scale;
                            if (arr->Ecoul0){ arr->Ecoul0[0][0][0][i3] *= scale; arr->Ecoul0[1][0][0][i3] *= scale; arr->Ecoul0[2][0][0][i3] *= scale; }
                        }
                    }
                  // 3. process remaining commands until "done"
                    arr->reset_for_calculation(false, true, true, true);
                    int ic_next = process_command_sequence(time_of_run, sys, arr, rdf,rdfs, pfout, nframe, time_stamp, ic+1);
                    if (ic_next<0){ ret = ic_next; break; }
                    if (sys->detail_level<1 && flog && istty){
                        char buffer[53]; buffer[0] = '['; buffer[51] = ']'; buffer[52] = 0;
                        for (int i=1; i<=50; i++) buffer[i] = i<=(iti+1)/(double)cmd[ic].step*50? '=' : ' ';
                        fprintf(flog, " TI progress %s%s", buffer, iti+1>=cmd[ic].step?"\n":"\r");
                    }
                  // 4. calculate energy
                    //arr->solvation_energy(sys, sys->ccutoff); double total_mass = arr->total_energy.mass;
                    if (iti+1==cmd[ic].step) print_default_IET_report(sys, arr); //arr->display_solvation_energy_full(sys, flog, nullptr);
                    double total_mass = arr->report_total->mass;
                  // 5. calculate energy and add to TI statistics, display at the last frame
                    (*arr->report_total) *= 1.0/cmd[ic].step; TI_energy += *arr->report_total; TI_energy.mass = total_mass;
                    if (iti+1==cmd[ic].step){
                        memcpy(&arr->report_total, &TI_energy, sizeof(TI_energy));
                        //arr->display_solvation_energy_full(sys, flog, nullptr, "", "totalTI", false, false, true);
                        print_default_IET_report(sys, arr, "totalTI");
                        ic = ic_next;
                    }
                  // 6. restore potential
                    if (scale!=1){
                        for (size_t i4=0; i4<N4; i4++){ arr->ulj[0][0][0][i4] /= scale; }
                        for (size_t i3=0; i3<N3; i3++){
                            arr->ucoulsr[0][0][i3] /= scale; arr->ucoullr[0][0][i3] /= scale;
                            if (arr->Ecoul0){ arr->Ecoul0[0][0][0][i3] /= scale; arr->Ecoul0[1][0][0][i3] /= scale; arr->Ecoul0[2][0][0][i3] /= scale; }
                        }
                    }
                }
            }
      #ifdef _EXPERIMENTAL_
        } else if (experimental_process_command(sys, arr, cmd, ic)){
      #endif
        }
      // no more command and will end here
    }

    memcpy(sys, &sys_bak, sizeof(IET_Param));
    return ret;
}
