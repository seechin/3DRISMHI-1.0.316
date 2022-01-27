
const char * get_command_name(int cmd){
    switch (cmd) {
        case IETCMD_NOP         : return "nop";
        case IETCMD_END         : return "end";
        case IETCMD_DONE        : return "done";
        case IETCMD_CLEAR       : return "clear";
        case IETCMD_RESET       : return "reset";
        case IETCMD_SET         : return "set";
        case IETCMD_SCALE       : return "scale";
        case IETCMD_LOAD        : return "load";
        case IETCMD_SAVE        : return "save";
        case IETCMD_SAVE_EXIST  : return "savee";
        case IETCMD_SAVE_FILTER : return "save-filter";
        case IETCMD_DISPLAY     : return "display";
        case IETCMD_REPORT      : return "report";
        case IETCMD_HI_SOLVER   : return "hi-solver (internal command)";
        case IETCMD_LSE         : return "lse";
        case IETCMD_CLOSURE     : return "closure";
        case IETCMD_CLOSURE_A   : return "closure-a";
        case IETCMD_CF          : return "cf";
        case IETCMD_CF_A        : return "cf-a";
        case IETCMD_DIELECT     : return "dielect";
        case IETCMD_DENSITY     : return "density";
        case IETCMD_BUILD_FF    : return "build-ff";
        case -IETCMD_BUILD_FF   : return "skip-ff";
        case IETCMD_BUILD_UUV   : return "build-uuv";
        case IETCMD_TI          : return "TI";
        case IETCMD_HOLD        : return "hold";
        case IETCMD_RDF_CONTENT : return "rdf-content";
        case IETCMD_TEMPERATURE : return "temperature";
        case IETCMD_TEST        : return "test";
        case IETCMD_TEST_SAVE   : return "test-and-save";
        default:
            if (cmd>=2000 && cmd<2000+sizeof(HIAL_name)/sizeof(HIAL_name[0])){
                return HIAL_name[cmd-2000];
            } else if (cmd>=1500 && cmd<1500+sizeof(IETAL_name)/sizeof(IETAL_name[0])){
                return IETAL_name[cmd-1500];
            } else return "(unknown)";
    }
}

const char * get_variable_name(int var){
    switch (var) {
        case IETCMD_v_temperature : return "temperature";
        case IETCMD_v_Coulomb     : return "Coulomb";
        case IETCMD_v_dielect_y   : return "dielect-y";
        case IETCMD_v_rbohr       : return "rbohr";
        case IETCMD_v_cmd         : return "cmd";
        case IETCMD_v_mass        : return "mass";
        case IETCMD_v_uuv         : return "uuv";
        case IETCMD_v_ulr         : return "ulr";
        case IETCMD_v_ulj         : return "lj";
        case IETCMD_v_ucoul       : return "coul";
        case IETCMD_v_ucoul2      : return "ef";
        case IETCMD_v_ucoulsr     : return "coulsr";
        case IETCMD_v_ucoullr     : return "coullr";
        case IETCMD_v_dd          : return "dd";
        case IETCMD_v_huv         : return "huv";
        case IETCMD_v_hlr         : return "hlr";
        case IETCMD_v_guv         : return "guv";
        case IETCMD_v_cuv         : return "cuv";
        case IETCMD_v_clr         : return "clr";
        case IETCMD_v_rmin        : return "rmin";
        case IETCMD_v_rdf         : return "rdf";
        case IETCMD_v_Euv         : return "Euv";
        case IETCMD_v_Ef          : return "EF";
        case IETCMD_v_DeltaN      : return "dN";
        case IETCMD_v_DeltaN0     : return "dN0";
        case IETCMD_v_TS          : return "TS";
        case IETCMD_v_rism_dielect: return "rism-dielect";
        case IETCMD_v_HFE         : return "HFE";
        case IETCMD_v_ddp         : return "ddp";
        case IETCMD_v_ld          : return "ld";
        case IETCMD_v_Yukawa      : return "Yukawa";
        case IETCMD_v_LocalCoulomb: return "local-coul";
        case IETCMD_v_PMV:          return "Volume";
        case IETCMD_v_Ef1:          return "Hef1";
        case IETCMD_v_csr:          return "csr";
        case -IETCMD_v_cuv:         return "chuv";
        case -IETCMD_v_clr:         return "chlr";
        case IETCMD_v_excess_GF:    return "exGF";
        case IETCMD_v_excess_RISM:  return "excess";
        case IETCMD_v_excess_hyb:   return "exHyb";
        case IETCMD_v_zeta_hnc:     return "zetaHNC";
        case IETCMD_v_zeta_closure: return "zetaClos";
        default: return "(unknown)";
    }
}
void print_command_vname_and_vvalues(FILE * out, IET_command * cmd){ for (int i=0; cmd->command_params_int[i]&&i<MAX_CMD_PARAMS; i++) fprintf(out, i==0?"%s=%g":",%s=%g", get_variable_name(cmd->command_params_int[i]), cmd->command_params_double[i]); }
void print_command_vnames(FILE * out, IET_command * cmd){ for (int i=0; cmd->command_params_int[i]&&i<MAX_CMD_PARAMS; i++) fprintf(out, i==0?"%s":",%s", get_variable_name(cmd->command_params_int[i])); }
void print_command_vvalues(FILE * out, IET_command * cmd){ for (int i=0; cmd->command_params_int[i]&&i<MAX_CMD_PARAMS; i++) fprintf(out, i==0?"%g":",%g", cmd->command_params_double[i]); }
void print_command_ivalues(FILE * out, IET_command * cmd){ for (int i=0; cmd->command_params_int[i]&&i<MAX_CMD_PARAMS; i++) fprintf(out, i==0?"%d":",%d", cmd->command_params_int[i]); }
const char * print_command(FILE * out, int command_id, IET_command * cmd, const char * prefix, const char * surfix){
    const char * command_name = get_command_name(cmd->command);
    fprintf(out, "%scmd[%d] = %s", prefix, command_id, command_name);
    switch (cmd->command) {
        case IETCMD_SET: case IETCMD_SCALE:
            fprintf(out, "("); print_command_vname_and_vvalues(out, cmd); fprintf(out, ")"); break;
        case IETCMD_LOAD: case IETCMD_SAVE: case IETCMD_SAVE_EXIST: case IETCMD_SAVE_FILTER: case IETCMD_DISPLAY: case IETCMD_REPORT:
            fprintf(out, "("); print_command_vnames(out, cmd); fprintf(out, ")"); break;
        case IETCMD_CLOSURE: case IETCMD_CLOSURE_A:
            fprintf(out, "("); for (int i=0; i<cmd->step&&i<MAX_CMD_PARAMS; i++) fprintf(out, i==0?"%s":",%s", CLOSURE_name[cmd->command_params_int[i]]); fprintf(out, ")");
            break;
        case IETCMD_CF: case IETCMD_CF_A: case IETCMD_DIELECT: case IETCMD_DENSITY:
            fprintf(out, "("); print_command_vvalues(out, cmd); fprintf(out, ")"); break;
        case IETCMD_HOLD:
            fprintf(out, "("); print_command_ivalues(out, cmd); fprintf(out, ")"); break;
        case IETCMD_TI:
            fprintf(out, ",step=%d", cmd->step);
    }
    if (cmd->command>=1500 && cmd->command<1500+sizeof(IETAL_name)/sizeof(IETAL_name[0])){  // RISM
        fprintf(out, ",step=%d", cmd->step);
    } else if (cmd->command>=2000 && cmd->command<2000+sizeof(HIAL_name)/sizeof(HIAL_name[0])){ // HI
        fprintf(out, ",step=%d", cmd->step);
    }
    fprintf(out, "%s\n", surfix);
    return command_name;
}


int maximum_default_processors(){
  #ifdef _LOCALPARALLEL_
    return sysconf(_SC_NPROCESSORS_ONLN);
  #else
    return 1;
  #endif
}
const char * get_short_fn(const char * filename){
    size_t i = strlen(filename)-1;
    while (i>0 && filename[i] && filename[i]!='/') i--;
    if (filename[i]=='/') i++;
    return &filename[i>0? i : 0];
}

void list_sys_files(IET_Param * sys, FILE * o, char * lineprefix=(char*)""){
    //fprintf(flog, "-b -e -s -f -log\n");
    bool istty = isatty(fileno(o));
    if (istty) fprintf(o, "%s", prompt_comment_prefix);
    fprintf(o, "%s%s %s\n", lineprefix, software_name, software_version);
    fprintf(o, "%spwd          %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_path), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-param       %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(info_file_name), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-solute      %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_solute), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    fprintf(o, "%s-traj        %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", get_second_fn(szfn_xtc), istty?prompt_path_suffix:"", istty?"\33[37m":"");
    if (!istty && !(!szfn_log[0] || StringNS::string(szfn_log)=="con"||StringNS::string(szfn_log)=="screen"||StringNS::string(szfn_log)=="stdout"||StringNS::string(szfn_log)=="stderr")){
        fprintf(o, "%s-log         %s\n", lineprefix, szfn_log[0]?get_second_fn(szfn_log):"con");
    } else {
        fprintf(o, "%s-log         %s%s%s%s\n", lineprefix, istty?prompt_path_prefix:"", szfn_log[0]?get_second_fn(szfn_log):"con", istty?prompt_path_suffix:"", istty?"\33[37m":"");
    }
    fprintf(o, "%s-be          %.0f,%.0f\n", lineprefix, sys->time_begin, sys->time_end);
    fprintf(o, "%s-dt          %g\n", lineprefix, sys->time_step);
    fprintf(o, "%s-nice        %d\n", lineprefix, sys->nice_level);
    if (istty) fprintf(o, "%s", prompt_comment_suffix);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void main_generate_default_gvv_map(IET_Param * sys, IET_arrays * arr, FILE * flog){
    int icol = 0; int iargv = 0; char gvv_map_buffer[128]; memset(gvv_map_buffer, 0, sizeof(gvv_map_buffer));
    char * pgvv_map[3] = { &gvv_map_buffer[0], &gvv_map_buffer[30], &gvv_map_buffer[60] };

    for (int iv=0; iv<sys->nv; iv++){
        for (int jv=iv; jv<sys->nv; jv++){
            snprintf(pgvv_map[0], 28, "%d", iv+1);
            snprintf(pgvv_map[1], 28, "%d", jv+1);
            snprintf(pgvv_map[2], 28, "%d", icol+1);

            //fprintf(flog, "    %4s %4s %4s\n", pgvv_map[0], pgvv_map[1], pgvv_map[2]);
            analysis_gvv_map(sys, pgvv_map, &iargv, 3, (char*)"inner", icol+1);
            icol ++;
        }
    }
    if (sys->gvv_map){
        if (sys->mode_test) fprintf(flog, "%s : [gvv_map] set to default%s\n", software_name, sys->mode_test?":":".");
        if (sys->mode_test){ fprintf(flog, "    %4s", ""); for (int iv=0; iv<sys->nv; iv++) fprintf(flog, " %4s", sys->av[iv].name); fprintf(flog, "\n"); }
        for (int iv=0; iv<sys->nv; iv++){
            if (sys->mode_test) fprintf(flog, "    %4s", sys->av[iv].name);
            for (int jv=0; jv<sys->nv; jv++){
                int icol = -1;
                for (int i=0; i<sys->n_gvv_map; i++) if ((sys->gvv_map[i].grpi==iv+1 && sys->gvv_map[i].grpj==jv+1)||(sys->gvv_map[i].grpi==jv+1 && sys->gvv_map[i].grpj==iv+1)) icol = sys->gvv_map[i].col;
                if (sys->mode_test) fprintf(flog, " %4d", icol);
            }
            if (sys->mode_test) fprintf(flog, "\n");
        }
    }
}

void build_force_field_auto(IET_Param * sys, IET_arrays * arr, FILE * flog, int nframe){
    if (arr->frame_stamp != nframe) arr->uuv_is_ready = false;
    if (!arr->uuv_is_ready){
        arr->reset_for_calculation(true, false, false, false);
     // short range, as well as short range part of coulomb_p2
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_sr()\n");
        build_force_field_sr(sys, arr);
      // long range
        if (sys->perform_pme){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_lr()\n");
            build_force_field_lr(sys, arr, true);
        } else clear_tensor3d(arr->ucoullr, arr->nx * arr->ny * arr->nz);
      // external field
        if (sys->do_external_electrostatic_field){
            int nx = sys->nr[0]; int ny = sys->nr[1]; int nz = sys->nr[2];
            for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) for (int iz=0; iz<nz; iz++){
                double delta_x = arr->box.x/nx*(ix-nx/2);
                double delta_y = arr->box.y/ny*(iy-ny/2);
                double delta_z = arr->box.z/nz*(iz-nz/2);
                arr->ucoullr[iz][iy][ix] += delta_x*sys->external_electrostatic_field.x + delta_y*sys->external_electrostatic_field.y + delta_z*sys->external_electrostatic_field.z;
            }
        }

      // done
        lap_timer_uuv();
        arr->frame_stamp = nframe;
        arr->uuv_is_ready = true;
    }
}
void build_uuv_base_on_force_field(IET_Param * sys, IET_arrays * arr, int es_algorithm=-1){
    if (es_algorithm<0) es_algorithm = sys->esal;
    if (es_algorithm==CoulAL_YukawaFFT){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_YukawaFFT(kappa=%g)\n", 1/sys->rc_yukawafft);
        build_force_field_uuv_YukawaFFT(sys, arr);
    } else if (es_algorithm==CoulAL_Dielect){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur(dielect=%g)\n", sys->mean_dielect);
        build_force_field_uuv_ur(sys, arr, sys->mean_dielect);
    } else {
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_uuv_ur()\n");
        build_force_field_uuv_ur(sys, arr);
    }

    //if (sys->ietal==IETAL_SSOZ){
    if (!sys->rism_coulomb_renormalization){
        if (sys->debug_level>=2){
            fprintf(sys->log(), "DEBUG:: electric field renormalization disabled for RISM\n");
            fprintf(sys->log(), "DEBUG:: build_force_field_uuv_post_merge()\n");
        }
        build_force_field_uuv_post_merge(sys, arr);
    } else {
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: electric field renormalization enabled for RISM\n");
    }


    if (sys->debug_level>=3||sys->debug_show_crc) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(uuv)        = %08X\n", check_real_crc(&arr->uuv[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));
    if (sys->debug_level>=3||sys->debug_show_crc) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(ulr)        = %08X\n", check_real_crc(&arr->ulr[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));
    if (sys->debug_level>=3||sys->debug_show_crc) fprintf(sys->log(), "DEBUG:: (debug only) check_real_crc(hlr)        = %08X\n", check_real_crc(&arr->hlr[0][0][0][0], arr->nv*arr->nx*arr->ny*arr->nz));

    //if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: build_force_field_coulp2(E)\n");
    //build_force_field_coulp2(sys, arr, arr->ulpbe, nullptr, arr->Ecoul[0], arr->Ecoul[1], arr->Ecoul[2]);


}

void reset_rdf(IET_Param * sys, RDF_data * rdf, int n_rdf_bins){
    for (int i1=0; i1<n_rdf_bins; i1++) for (int ig=0; ig<sys->n_rdf_grps; ig++){
        rdf[ig].n[i1] = 0;
        rdf[ig].g[i1] = 0;
    }
}

void calculate_rdf(IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    if (!rdf || sys->n_rdf_grps<=0) return;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3;
    if (sys->rdf_content==IETCMD_v_dd){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * dd1 = arr->dd? &arr->dd[sys->av[iv].iaa][0][0][0] : nullptr;
            double nbulk = sys->nbulk[sys->av[iv].iaa];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (dd1? dd1[i3] : nbulk);
        }
    } else if (sys->rdf_content==-IETCMD_v_dd){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * dd1 = &arr->nphi[sys->av[iv].iaa][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (dd1? dd1[i3] : 0);
        }
    } else if (sys->rdf_content==IETCMD_v_ucoul){
        for (int iv=0; iv<sys->nv; iv++){
            double charge = sys->av[iv].charge;
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = &arr->ucoulsr[0][0][0];
            __REAL__ * src2 = &arr->ucoullr[0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (src1[i3] + src2[i3]) * charge * sys->scale_coul;
        }
    } else if (sys->rdf_content==IETCMD_v_Ef){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = &arr->Ecoul0[0][0][0][0];
            __REAL__ * src2 = &arr->Ecoul0[1][0][0][0];
            __REAL__ * src3 = &arr->Ecoul0[2][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (0.5/(4*PI*COULCOOEF)) * Vector(src1[i3], src2[i3], src3[i3]).pow2() * sys->scale_coul;
        }
    } else if (sys->rdf_content==IETCMD_v_cuv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * cuv1 = &arr->cuv[iv][0][0][0];
            __REAL__ * ulr1 = &arr->ulr[iv][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = cuv1[i3] - ulr1[i3];
        }
    } else if (sys->rdf_content==-IETCMD_v_cuv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * ulr1 = &arr->ulr[iv][0][0][0];
            __REAL__ * hlr1 = &arr->hlr[iv][0][0][0];
            __REAL__ * huv1 = &arr->huv[iv][0][0][0];
            __REAL__ * cuv1 = &arr->cuv[iv][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = huv1[i3] - cuv1[i3] + hlr1[i3] + ulr1[i3];
        }
    } else if (sys->rdf_content==IETCMD_v_ulj){
        for (size_t i4=0; i4<N4; i4++) arr->res[0][0][0][i4] = arr->ulj[0][0][0][i4]>sys->ccutoff? 0 : arr->ulj[0][0][0][i4];
    } else if (sys->rdf_content==IETCMD_v_guv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * huv1 = &arr->huv[iv][0][0][0];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = huv1[i3] + 1;
        }
    } else if (sys->rdf_content==IETCMD_v_uuv || sys->rdf_content==IETCMD_v_huv || sys->rdf_content==-IETCMD_v_huv){
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * src1 = nullptr;
            double scaling = 1;
            if (sys->rdf_content==IETCMD_v_ulj){        src1 = &arr->ulj[iv][0][0][0];    scaling = sys->scale_lj;
            } else if (sys->rdf_content==IETCMD_v_uuv){ src1 = &arr->uuv[iv][0][0][0];    scaling = sys->scale_coul;
            } else if (sys->rdf_content==IETCMD_v_huv){ src1 = &arr->huv[iv][0][0][0];
            } else if (sys->rdf_content==-IETCMD_v_huv){ src1 = &arr->hlr[iv][0][0][0];
            } else {                                    src1 = &arr->huv[iv][0][0][0];
            }
            for (size_t i3=0; i3<N3; i3++) res1[i3] = src1[i3] * scaling;
        }
    } else { // guv
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * res1 = &arr->res[iv][0][0][0];
            __REAL__ * huv1 = &arr->huv[iv][0][0][0];
            __REAL__ * dd1 = arr->dd? &arr->dd[sys->av[iv].iaa][0][0][0] : nullptr;
            double nbulk = sys->nbulk[sys->av[iv].iaa];
            double nbulk_rism = sys->nbulk_rism[sys->av[iv].iaa];
            for (size_t i3=0; i3<N3; i3++) res1[i3] = (1+huv1[i3]) * (dd1? dd1[i3] : nbulk) / nbulk_rism;
        }
    }

    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        memset(rdf[ig].g, 0, sizeof(double)*N1); memset(rdf[ig].n, 0, sizeof(double)*N1);
        int is = rdf[ig].is; int iv = rdf[ig].iv;
          if (iv<0 || iv>=sys->nv) continue;
          if (is>=sys->nas) continue;
        if (is<0 && arr->r2uvmin){
            for (int iz=0; iz<sys->nr[2]; iz++) for (int iy=0; iy<sys->nr[1]; iy++) for (int ix=0; ix<sys->nr[0]; ix++){
                double r = sqrt(fabs(arr->r2uvmin[iz][iy][ix]));
                int ithis = (int)(floor(r/dr)); if (ithis>=0 && ithis<N1){
                    rdf[ig].g[ithis] += arr->res[iv][iz][iy][ix];
                    rdf[ig].n[ithis] ++;
                }
            }
        } else {
            Vector center = sys->traj.atom[is].r; //printf("calculate rdf for pair: %d - %d, ref_center: %g %g %g\n", rdf[ig].is, rdf[ig].iv, center.x, center.y, center.z);
            for (int iz=0; iz<sys->nr[2]; iz++) for (int iy=0; iy<sys->nr[1]; iy++) for (int ix=0; ix<sys->nr[0]; ix++){
                //double r = (center - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz)).mod();

                //double cx = center.x; if (sys->pbc_x) cx = ff_pbc_d(cx, arr->box.x);
                //double cy = center.y; if (sys->pbc_y) cy = ff_pbc_d(cy, arr->box.y);
                //double cz = center.z; if (sys->pbc_z) cz = ff_pbc_d(cz, arr->box.z);
                //double r = (Vector(cx, cy, cz) - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz)).mod();

                Vector dv = Vector(center - Vector(ix*arr->dx, iy*arr->dy, iz*arr->dz));
                if (sys->pbc_x) dv.x = ff_pbc_d(dv.x+arr->box.x/2, arr->box.x)-arr->box.x/2;
                if (sys->pbc_y) dv.y = ff_pbc_d(dv.y+arr->box.y/2, arr->box.y)-arr->box.y/2;
                if (sys->pbc_z) dv.z = ff_pbc_d(dv.z+arr->box.z/2, arr->box.z)-arr->box.z/2;
                double r = dv.mod();

                int ithis = (int)(floor(r/dr)); if (ithis>=0 && ithis<N1){
                    rdf[ig].g[ithis] += arr->res[iv][iz][iy][ix];
                    rdf[ig].n[ithis] ++;
                }
            }
        }
    }
}
void add_rdf_to_sum(IET_Param * sys, IET_arrays * arr, RDF_data * rdf, RDF_data * rdfs, int n_rdf_bins){
    int N1 = n_rdf_bins;
    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        for (int i=0; i<N1; i++){
            rdfs[ig].g[i] += rdf[ig].g[i];
            rdfs[ig].n[i] += rdf[ig].n[i];
        }
    }
}

bool print_rdf(FILE * file, IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    if (!rdf || sys->n_rdf_grps<=0) return false;
    if (!file) return false;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    unsigned long mask64 = generate_double_sig_dig_mask(sys->output_significant_digits);
      if (sys->output_significant_digits==-4) mask64 = generate_double_sig_dig_mask(7);
      if (sys->output_significant_digits==-8) mask64 = generate_double_sig_dig_mask(15);

    int grp_min = sys->rdf_grps[0].grp; int grp_max = sys->rdf_grps[0].grp;
    for (int i=1; i<sys->n_rdf_grps; i++){ if (sys->rdf_grps[i].grp<grp_min) grp_min = sys->rdf_grps[i].grp; if (sys->rdf_grps[i].grp>grp_max) grp_max = sys->rdf_grps[i].grp; }

  // 1. print titles
    fprintf(file, sys->output_title_format, "r");
    for (int igrp=grp_min; igrp<=grp_max; igrp++){
        char buffer[MAX_NAME];
        for (int ig=0; ig<sys->n_rdf_grps; ig++) if (sys->rdf_grps[ig].grp==igrp){
            if (rdf[ig].is<0){
                snprintf(buffer, sizeof(buffer), "all-%d/%s.%s", rdf[ig].iv+1, sys->av[rdf[ig].iv].mole, sys->av[rdf[ig].iv].name);
            } else {
                snprintf(buffer, sizeof(buffer), "%d/%s.%s-%d/%s.%s", rdf[ig].is+1, sys->as[rdf[ig].is].mole, sys->as[rdf[ig].is].name, rdf[ig].iv+1, sys->av[rdf[ig].iv].mole, sys->av[rdf[ig].iv].name);
            }
            fprintf(file, sys->output_title_format, buffer);
            break;
        }
    }; fprintf(file, "\n");
  // 2. print the RDF data
    for (int i1=0; i1<N1; i1++){
        fprintf(file, sys->output_realnumber_format, (i1+0.5)*dr);
        for (int igrp=grp_min; igrp<=grp_max; igrp++){
            double value = 0; double n_value = 0; int encounter_times = 0;
            for (int ig=0; ig<sys->n_rdf_grps; ig++) if (sys->rdf_grps[ig].grp==igrp){
                value += rdf[ig].g[i1];
                n_value += rdf[ig].n[i1];
                encounter_times ++;

                //value += rdf[ig].n[i1]!=0? rdf[ig].g[i1]/rdf[ig].n[i1] : rdf[ig].g[i1];
                //n_value ++;
            }
            if (encounter_times>0) fprintf(file, sys->output_realnumber_format, n_value!=0? value/n_value : value);
        }
        fprintf(file, "\n");
    }

  /*
    //fprintf(file, "%11s", "r");
    fprintf(file, sys->output_title_format, "r");
    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        char buffer[MAX_NAME];
        if (rdf[ig].is<0){
            snprintf(buffer, sizeof(buffer), "all-%d%s", rdf[ig].iv+1, sys->av[rdf[ig].iv].name);
        } else {
            snprintf(buffer, sizeof(buffer), "%d%s-%d%s", rdf[ig].is+1, sys->as[rdf[ig].is].name, rdf[ig].iv+1, sys->av[rdf[ig].iv].nele);
        }
        //fprintf(file, " %11s", buffer);
        fprintf(file, sys->output_title_format, buffer);
    }; fprintf(file, "\n");
    for (int i1=0; i1<N1; i1++){
        //fprintf(file, "%11.4f", (i1+0.5)*dr);
        fprintf(file, sys->output_realnumber_format, (i1+0.5)*dr);
        for (int ig=0; ig<sys->n_rdf_grps; ig++){
            double value = rdf[ig].n[i1]!=0? rdf[ig].g[i1]/rdf[ig].n[i1] : rdf[ig].g[i1];
            //bit_and(&value, &value, &mask64, 8);
            //fprintf(file, " %11.4f", value);
            fprintf(file, sys->output_realnumber_format, value);
        }
        fprintf(file, "\n");
    }
  */

    return true;
}

bool print_rdf(const char * filename, IET_Param * sys, IET_arrays * arr, RDF_data * rdf, double r_max_rdf, int n_rdf_bins){
    bool fp_by_fopen = false;
    if (!rdf || sys->n_rdf_grps<=0) return false;
    if (!filename) return false;
    double dr = r_max_rdf / n_rdf_bins; int N1 = n_rdf_bins;
    FILE * file = stdout; StringNS::string fn = filename;
    if (fn=="stdout" || fn=="screen" || fn=="con"){ file = stdout;
    } else if (fn=="stderr"){   file = stderr;
    } else {    file = fopen(filename, "w"); fp_by_fopen = true;
    }
    if (!file){ fprintf(sys->log(), "%s%s : error : cannot write RDF to %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, filename, sys->is_log_tty?color_string_end:""); return false; }
    bool ret = print_rdf(file, sys, arr, rdf, r_max_rdf, n_rdf_bins);
    if (fp_by_fopen){
        size_t rdf_file_size = ftell(file);
        fclose(file);

        char save_items_mem_buffer[64];
        fprintf(sys->log(), "  save %s to %s%s%s, totally %s\n", "rdf", sys->is_log_tty?prompt_path_prefix:"\"", filename, sys->is_log_tty?prompt_path_suffix:"\"", print_memory_value(save_items_mem_buffer, sizeof(save_items_mem_buffer), rdf_file_size));
    }
    return ret;
}


const char * get_report_title(int print_item){
    switch (print_item) {
        case IETCMD_v_mass:         return "mass";
        case IETCMD_v_ulj :         return "LJSR";
        case IETCMD_v_ucoul:        return "Coul";
        case IETCMD_v_TS:           return "-TS";
        case IETCMD_v_DeltaN:       return "DN";
        case IETCMD_v_DeltaN0:      return "DN_vac";
        case IETCMD_v_PMV:          return "Volume";
        case IETCMD_v_Ef:           return "Hef0";
        case IETCMD_v_Ef1:          return "Hef1";
        case IETCMD_v_cuv:          return "cuv";
        case IETCMD_v_csr:          return "csr";
        case IETCMD_v_clr:          return "clr";
        case -IETCMD_v_cuv:         return "chuv";
        case -IETCMD_v_clr:         return "chlr";
        case IETCMD_v_excess_GF:    return "exGF";
        case IETCMD_v_excess_RISM:  return "excess";
        case IETCMD_v_excess_hyb:   return "exHyb";
        case IETCMD_v_zeta_hnc:     return "zetaHNC";
        case IETCMD_v_zeta_closure: return "zetaClos";
    }
    return "";
}
double get_report_value(IET_Report * report, int print_item){
    switch (print_item) {
        case IETCMD_v_mass:         return report->mass;
        case IETCMD_v_ulj :         return report->lj;
        case IETCMD_v_ucoul:        return report->coulsr + report->coullr;
        case IETCMD_v_TS:           return report->entropy;
        case IETCMD_v_DeltaN:       return report->dN;
        case IETCMD_v_DeltaN0:      return report->dNg;
        case IETCMD_v_PMV:          return report->dN;
        case IETCMD_v_Ef:           return report->Uef0;
        case IETCMD_v_Ef1:          return report->Uef1;
        case IETCMD_v_cuv:          return report->cuv + report->clr;
        case IETCMD_v_csr:          return report->cuv;
        case IETCMD_v_clr:          return report->clr;
        case -IETCMD_v_cuv:         return report->chuv + report->chlr;
        case -IETCMD_v_csr:         return report->chuv;
        case -IETCMD_v_clr:         return report->chlr;
        case IETCMD_v_excess_GF:    return report->excess_chem[0];
        case IETCMD_v_excess_RISM:  return report->excess_chem[1];
        case IETCMD_v_excess_hyb:   return report->excess_chem[2];
        case IETCMD_v_zeta_hnc:     return report->zeta[0];
        case IETCMD_v_zeta_closure: return report->zeta[2];
    }
    return 0;
}

void display_IET_report(IET_Param * sys, IET_arrays * arr, IET_command * cmd, int ic, const char * total_title="total"){
    FILE * flog = sys->log();
    bool ex_display_rdf = false;

    double average_mass = 0; double total_volume = 0;
    for (int iv=0; iv<sys->nv; iv++){
        average_mass += sys->av[iv].mass * sys->av[iv].multi * sys->nbulk[sys->av[iv].iaa];
        total_volume -= arr->report_sites[iv].dN / sys->bulk_density_mv[sys->av[iv].iaa] * arr->report_sites[iv].mass / arr->report_sites[iv].mass_mol;
    }
  // print title
    int n_print_count = 0; for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++) if (cmd[ic].command_params_int[is]==IETCMD_v_rdf); else n_print_count ++;
    if (sys->detail_level>=1 && n_print_count>0){
        fprintf(flog, "%9s", "#item");
        fprintf(flog, sys->output_title_format, total_title);
        for (int iv=0; iv<sys->nv; iv++) fprintf(flog, sys->output_title_format, sys->av[iv].name);
        fprintf(flog, "\n");
    }
    for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
      // print item name
        if (get_report_title(cmd[ic].command_params_int[is])[0]) fprintf(flog, "%9s", get_report_title(cmd[ic].command_params_int[is]));
      // first display total
        if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){
            ex_display_rdf = true;
            continue;
        } else if (cmd[ic].command_params_int[is]==IETCMD_v_mass){
            fprintf(flog, sys->output_realnumber_format, average_mass);
        } else if (cmd[ic].command_params_int[is]==IETCMD_v_PMV){
            fprintf(flog, sys->output_realnumber_format, total_volume);
        } else {
            double value = get_report_value(arr->report_total, cmd[ic].command_params_int[is]);
            fprintf(flog, sys->output_realnumber_format, value);
        }
      // then display all sites
        for (int iv=0; iv<sys->nv; iv++){
            if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){
                ex_display_rdf = true;
            } else if (cmd[ic].command_params_int[is]==IETCMD_v_mass){
                char mass_string[32];
                snprintf(mass_string, sizeof(mass_string), "%.0fx%d", sys->av[iv].mass, sys->av[iv].multi);
                fprintf(flog, sys->output_realnumber_format, mass_string);
                average_mass += sys->av[iv].mass * sys->av[iv].multi * sys->nbulk[sys->av[iv].iaa];
            } else if (cmd[ic].command_params_int[is]==IETCMD_v_PMV){
                double this_volume = - arr->report_sites[iv].dN / sys->bulk_density_mv[sys->av[iv].iaa];
                fprintf(flog, sys->output_realnumber_format, this_volume);
                total_volume += this_volume * arr->report_sites[iv].mass / arr->report_sites[iv].mass_mol;
            } else {
                double value = get_report_value(&arr->report_sites[iv], cmd[ic].command_params_int[is]);
                fprintf(flog, sys->output_realnumber_format, value);
            }
        }
      // finally
        fprintf(flog, "\n");
    }
}

void print_IET_report(IET_Param * sys, IET_arrays * arr, IET_command * cmd, int ic, const char * total_title="total"){
    FILE * flog = sys->log();
    bool ex_display_rdf = false;

    double average_mass = 0; double total_volume = 0;
    for (int iv=0; iv<sys->nv; iv++){
        average_mass += sys->av[iv].mass * sys->av[iv].multi * sys->nbulk[sys->av[iv].iaa];
        total_volume -= arr->report_sites[iv].dN / sys->bulk_density_mv[sys->av[iv].iaa] * arr->report_sites[iv].mass / arr->report_sites[iv].mass_mol;
    }

    int n_print_count = 0; for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++) if (cmd[ic].command_params_int[is]==IETCMD_v_rdf); else n_print_count ++;

    if (n_print_count>0){
        fprintf(flog, "%7s", "Atom");
        for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
            if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){
                ex_display_rdf = true;
            } else if (cmd[ic].command_params_int[is]==IETCMD_v_mass){
                fprintf(flog, " %7s", "mass");
            } else {
                fprintf(flog, sys->output_title_format, get_report_title(cmd[ic].command_params_int[is]));
            }
        }
        fprintf(flog, "\n");
    }
  // print each sites
    for (int iv=0; iv<sys->nv; iv++){
        if (n_print_count>0) fprintf(flog, "%7s", sys->av[iv].name);
        for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
            if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){
                ex_display_rdf = true;
            } else if (cmd[ic].command_params_int[is]==IETCMD_v_mass){
                char mass_string[32];
                snprintf(mass_string, sizeof(mass_string), "%.0fx%d", sys->av[iv].mass, sys->av[iv].multi);
                fprintf(flog, " %7s", mass_string);
            } else if (cmd[ic].command_params_int[is]==IETCMD_v_PMV){
                double this_volume = - arr->report_sites[iv].dN / sys->bulk_density_mv[sys->av[iv].iaa];
                fprintf(flog, sys->output_realnumber_format, this_volume);
            } else {
                double value = get_report_value(&arr->report_sites[iv], cmd[ic].command_params_int[is]);
                fprintf(flog, sys->output_realnumber_format, value);
            }
        }
        if (n_print_count>0) fprintf(flog, "\n");
    }
  // print total
    if (n_print_count>0) fprintf(flog, "%7s", total_title);
    for (int is=0; is<cmd[ic].step && is<MAX_CMD_PARAMS; is++){
        if (cmd[ic].command_params_int[is]==IETCMD_v_rdf){
            ex_display_rdf = true;
        } else if (cmd[ic].command_params_int[is]==IETCMD_v_mass){
            fprintf(flog, " %7g", average_mass);
        } else if (cmd[ic].command_params_int[is]==IETCMD_v_PMV){
            fprintf(flog, sys->output_realnumber_format, total_volume);
        } else {
            double value = get_report_value(arr->report_total, cmd[ic].command_params_int[is]);
            fprintf(flog, sys->output_realnumber_format, value);
        }
    }
    if (n_print_count>0) fprintf(flog, "\n");
}

void print_default_IET_report(IET_Param * sys, IET_arrays * arr, const char * total_title="total"){
    FILE * flog = sys->log();
    IET_command cmd;
    cmd.command = IETCMD_REPORT; cmd.step = 0;
    cmd.command_params_int[cmd.step++] = IETCMD_v_mass;
    cmd.command_params_int[cmd.step++] = IETCMD_v_DeltaN;
    cmd.command_params_int[cmd.step++] = IETCMD_v_DeltaN0;
    cmd.command_params_int[cmd.step++] = IETCMD_v_TS;
    cmd.command_params_int[cmd.step++] = IETCMD_v_ulj;
    cmd.command_params_int[cmd.step++] = IETCMD_v_ucoul;
    cmd.command_params_int[cmd.step++] = IETCMD_v_PMV;
    print_IET_report(sys, arr, &cmd, 0, total_title);
}

void apply_dielect_from_dipole(IET_Param * sys){
    sys->dielect_from_dipole = true;
    double beta = sys->default_temperature / sys->temperature;
    for (int ivm=0; ivm<sys->nvm; ivm++){
        sys->dielect_mol[ivm] = 1 + 4*PI/3*COULCOOEF * beta * sys->density_mv[ivm] * sys->dipole_mv[ivm]*sys->dipole_mv[ivm];
    }
    for (int i=0; i<sys->nav; i++){ int iaa = sys->av[i].iaa; if (iaa>=0 && iaa<MAX_CMD_PARAMS) sys->dielect[i] = sys->dielect_mol[iaa]; }
    double dielecti = 0; double density_dielect = 0; for (int ivm=0; ivm<sys->nvm; ivm++){ dielecti += sys->density_mv[ivm] / sys->dielect_mol[ivm]; density_dielect += sys->density_mv[ivm]; } sys->mean_dielect = density_dielect / dielecti;
    //if (sys->esal==CoulAL_Coulomb) sys->mean_dielect = 1;
}

void generate_default_output_filename(IET_Param * sys, char * filename, int size_filename, const char * traj_filename, int size_traj_filename){
    int ixtc = 0; int jxtc = 0;
    int len_traj_filename = strnlen(traj_filename, size_traj_filename);
    for (ixtc=len_traj_filename; ixtc>=0; ixtc--) if (traj_filename[ixtc]=='/') break;
    while (ixtc<size_traj_filename && traj_filename[ixtc]=='/') ixtc++;
    for (jxtc=len_traj_filename; jxtc>=0; jxtc--) if (traj_filename[jxtc]=='.') break;
    if (ixtc<0||ixtc>=len_traj_filename) ixtc = 0; if (jxtc<=ixtc) jxtc = len_traj_filename;
    char trajname[MAX_PATH]; memset(trajname, 0, sizeof(trajname)); memcpy(trajname, &traj_filename[ixtc], jxtc-ixtc);

    if (!filename[0]){
        time_t rawtime; time(&rawtime); struct tm * tmtime = localtime(&rawtime);
        char time_buffer[40]; snprintf(time_buffer, sizeof(time_buffer), "%02d%02d%02d_%02d%02d", tmtime->tm_year-100, tmtime->tm_mon+1, tmtime->tm_mday, tmtime->tm_hour, tmtime->tm_min);
        int nmsolute = 1; for (int i=1; i<sys->nas; i++) if (StringNS::string(sys->as[i].mole) != sys->as[i-1].mole) nmsolute ++;
        int nmsolvent = 1; for (int i=1; i<sys->nav; i++) if (StringNS::string(sys->av[i].mole) != sys->av[i-1].mole) nmsolvent ++;
        snprintf(filename, size_filename, "%s.%s%s%s.%s.ts4s", trajname, sys->av[0].mole, nmsolvent>1?"_":"", nmsolvent>2?"etc":nmsolvent==2?sys->av[sys->nav-1].mole:"", time_buffer);

        sys->b_output_filename_autogenerated = true;
    }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void main_print_header(IET_Param * sys, IET_arrays * arr, FILE * flog, int argc, char * argv[]){
    bool islogtty = isatty(fileno(flog));
        fprintf(flog, "============================================================================ \n");
        fprintf(flog, "%s %s\n", software_name, software_version);
        char time_buffer[40];
        fprintf(flog, "%s begins at %s\n%s@%s: %s%s%s\n", software_name, get_current_time_text(time_buffer), username, hostname, islogtty?prompt_path_prefix:"", szfn_path, islogtty?prompt_path_suffix:"");
        for (int i=0; i<argc; i++) fprintf(flog, "%s ", argv[i]); fprintf(flog, "\n");
        fprintf(flog, "---------------------------------------------------------------------------- \n");
}
void main_print_tailer(IET_Param * sys, IET_arrays * arr, bool success=true){
    FILE * flog = sys->log();
    char time_buffer[2][40];
    if (success){
      #ifdef _LOCALPARALLEL_
        if (sys->nt>1){
            fprintf(flog, "%s ends at %s (%d %s, %s s)\n", software_name, get_current_time_text(time_buffer[0]), sys->nt, sys->mp_by_fork?"forks":"threads", display_time(__total_timer, time_buffer[1]));
        } else {
            fprintf(flog, "%s ends at %s (%s s)\n", software_name, get_current_time_text(time_buffer[0]), display_time(__total_timer, time_buffer[1]));
        }
      #else
        fprintf(flog, "%s ends at %s (%s s)\n", software_name, get_current_time_text(time_buffer[0]), display_time(__total_timer, time_buffer[1]));
      #endif
    }
    fprintf(flog, "============================================================================ \n");
}

bool main_prepair_rdf_grps_to_pairs(IET_Param * sys, IET_arrays * arr, RDF_data ** prdf, int * pn_rdf_pairs){
    bool success = true;

    RDF_data * rdf = nullptr;
    int n_rdf_pairs = 0;

    int npr = map_rdf_grps_to_pairs(sys, nullptr, sys->detail_level>=1); //printf("pair number %d, bins %d\n", npr, sys->out_rdf_bins);
    if (npr>0){
        rdf = (RDF_data*) memalloc( (sizeof(RDF_data)+2*(sizeof(double)*(1+sys->out_rdf_bins))) * npr ); n_rdf_pairs = npr;
        if (rdf) for (int i=0; i<npr; i++){
            rdf[i].is = rdf[i].iv = 0;
            rdf[i].g = (double*) (((char *)rdf) + sizeof(RDF_data)*npr + 2*i*(sizeof(double)*(1+sys->out_rdf_bins)));
            rdf[i].n = (double*) (((char *)rdf) + sizeof(RDF_data)*npr + (2*i+1)*(sizeof(double)*(1+sys->out_rdf_bins)));
            for (int j=0; j<sys->out_rdf_bins; j++) rdf[i].g[j] = rdf[i].n[j] = 0;
        }
        int npr_mapped = map_rdf_grps_to_pairs(sys, rdf, false); //for (int i=0; i<n_rdf_pairs; i++) printf("rdf_pair[%d] = %d(%s.%s)-%d(%s.%s)\n", i+1, rdf[i].is+1, sys->as[rdf[i].is].mole,sys->as[rdf[i].is].name, rdf[i].iv+1, sys->av[rdf[i].iv].mole, sys->av[rdf[i].iv].name);
        if (npr_mapped<0) success = false;
    } else success = false;
    if (!success){
        fprintf(sys->log(), "%s : please use the following indicators to specify RDF groups: %s\n", software_name, (sys->debug_level>=3)?"":"(-debug 3 for full lists)");
        int max_display_items = 10;
        if (sys->debug_level>=3) max_display_items = sys->nas;
        fprintf(sys->log(), "  solutes:"); for (int i=0; i<sys->nas && i<max_display_items; i++) fprintf(sys->log(), " %d/%s:%s", i+1, sys->as[i].mole, sys->as[i].name); fprintf(sys->log(), sys->nas>max_display_items?" ...\n":"\n");
        if (sys->debug_level>=3) max_display_items = sys->nv;
        fprintf(sys->log(), "  solvents:"); for (int i=0; i<sys->nv && i<max_display_items; i++) fprintf(sys->log(), " %d/%s:%s", i+1, sys->av[i].mole, sys->av[i].name); fprintf(sys->log(), sys->nv>max_display_items?" ...\n":"\n");
    }

    *prdf = rdf;
    *pn_rdf_pairs = n_rdf_pairs;

    return success;
}
void main_print_solute_list(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    int max_display_solute_num = 10;
    //if (sys->detail_level>=2 && sys->debug_level>=2) max_display_solute_num = 99999999;
    fprintf(flog, "%ssolutes: %d atoms in %s: ", prefix, sys->nas, get_short_fn(szfn_xtc));
      int iaa_count = 0;
      for (int i=0; i<sys->nas; i++) if (i==0||StringNS::string(sys->as[i].mole)!=sys->as[i-1].mole){
          iaa_count++;
          if (iaa_count<max_display_solute_num) fprintf(flog, i==0?"%s":", %s", sys->as[i].mole);
          else if (iaa_count==max_display_solute_num) fprintf(flog, " ...");
      };
    if (iaa_count>=max_display_solute_num)fprintf(flog, " (%d residues)", iaa_count); fprintf(flog, "\n");
}
void main_print_solvent_list(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    fprintf(flog, "%ssolvents: %d %s in: ", prefix, sys->nav, sys->nav>1?"sites":"site");
    for (int i=0; i<sys->nav; i++) if (i==0||StringNS::string(sys->av[i].mole)!=sys->av[i-1].mole){
        fprintf(flog, i==0?"%s":", %s", sys->av[i].mole);
    }
    fprintf(flog, "\n");
    fprintf(flog, "%ssolvent densities: ", prefix);
    for (int i=0; i<sys->nav; i++) if (i==0||StringNS::string(sys->av[i].mole)!=sys->av[i-1].mole){
        int iaa = sys->av[i].iaa;
        fprintf(flog, iaa==0?"%g/%g" : ", %g/%g", sys->density_mv[iaa], sys->bulk_density_mv[iaa]);
    }
    fprintf(flog, "\n");
}
void main_print_forcefield_info(IET_Param * sys, IET_arrays * arr, FILE * flog, const char * prefix){
    const char * ffstring = "(undefined)";
    if (!sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 502.97)<=MACHINE_REASONABLE_ERROR) ffstring = "amber";
    if (!sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 120.27)<=MACHINE_REASONABLE_ERROR) ffstring = "gaff";
    if (sys->mixingrule_sigma_geometric && fabs(sys->default_temperature- 120.27)<=MACHINE_REASONABLE_ERROR) ffstring = "oplsaa";
    fprintf(flog, "%sforcefield %s, rvdw %g, rcoul %g, temperature %g\n", prefix, ffstring, sys->rvdw, sys->rcoul, sys->temperature);
}

void print_save_extra_info(IET_Param * sys, char * buffer, size_t size){
    char time_buffer[40];
    //const char * short_info = get_short_fn(get_second_fn(info_file_name));
    //const char * short_solu = get_short_fn(get_second_fn(szfn_solute));
    //const char * short_traj = get_short_fn(get_second_fn(szfn_xtc));
    //snprintf(buffer, size, "%s %s %s -p %s -s %s -f %s -pwd %s", software_name, software_version, get_current_time_text(time_buffer), short_info, short_solu, short_traj, szfn_path);
    snprintf(buffer, size, "%s %s %s", software_name, software_version, get_current_time_text(time_buffer));
}

size_t select_append_save_data(int * filter, int nfilter, __REAL__ **** temp, size_t * total_size, size_t * original_size, FILE ** pfile, char filename[MAX_PATH], FILE * flog, const char * title, const char * text, int nx, int ny, int nz, int nv, __REAL__ **** data, double time_stamp, IET_Param * sys, __REAL__ * _compressBuf, size_t _allocate_memory_size){

    bool save_all = false; bool save_none = true; int N3 = nx * ny * nz;
    if (!filter || nfilter<=0){
        save_all = true; save_none = false;
    } else{
        for (int i=0; i<nfilter; i++){
            if (filter[i]<0){ save_all = true; save_none = false; }
            else if (filter[i]>0) save_none = false;
        }
    }
    if (save_none){
//fprintf(stderr, "\33[31msave: none\n\33[0m");
        *total_size = *original_size = 0;
    } else if (save_all){
//fprintf(stderr, "\33[31msave: all\n\33[0m");
        *total_size += append_save_data(pfile, filename, flog, title, text, nx, ny, nz, nv, &data[0][0][0][0], time_stamp, sys, _compressBuf, _allocate_memory_size);
        *original_size = nx*ny*nz*nv*sizeof(__REAL__) + sizeof(CompressPageHeader);
    } else {
        int n_save_sites = 0;
        for (int i=0; i<nfilter; i++) if (filter[i]>0 && filter[i]<=nv && n_save_sites<nv){
//fprintf(stderr, "\33[31msave: item[%d]=%d (of %d)\n\33[0m", i, filter[i]-1, nfilter);
            cp_tensor3d(data[filter[i]-1], temp[n_save_sites++], N3);

            //int N3 = nx*ny*nz;
            //__REAL__ * src = &data[n_save_sites][0][0][0]; __REAL__ * dst = &temp[filter[i]-1][0][0][0];
            //for (int j=0; j<N3; j++) dst[j] = src[j];
            //n_save_sites ++;
        }
        *total_size = append_save_data(pfile, filename, flog, title, text, nx, ny, nz, n_save_sites, &temp[0][0][0][0], time_stamp, sys, _compressBuf, _allocate_memory_size);
        *original_size = nx*ny*nz*n_save_sites*sizeof(__REAL__) + sizeof(CompressPageHeader);
    }
    return *total_size;
}

void perform_sub_test(IET_Param * sys, IET_arrays * arr, IET_command * cmd){
    if (!cmd) return;
    PDBAtom * ia = sys->traj.atom; SoluteAtomSite * as = sys->as;
    size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = N3 * sys->nv;
    __REAL__ **** cache1 = arr->rismhi_cache[2]; __REAL__ **** cache2 = arr->rismhi_cache[3];

    for (int iparam=0; iparam<cmd->step; iparam++){
        if (cmd->command_params_int[iparam]==IETCMD_v_ucoul || cmd->command_params_int[iparam]==IETCMD_v_hlr || cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_sub_test(%s%s)\n", cmd->command_params_int[iparam]==IETCMD_v_ucoul?"coulomb":cmd->command_params_int[iparam]==IETCMD_v_hlr?"hlr":cmd->command_params_int[iparam]==IETCMD_v_Yukawa?"Yukawa":"unknown", cmd->command==IETCMD_TEST_SAVE?", save":"");

            if (cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
                build_force_field_uuv_YukawaFFT(sys, arr);
            } else {
                build_uuv_base_on_force_field(sys, arr);
            }
            RISMHI3D_RISMNS::prepare_3d_rism(sys, arr);
            /*
            for (size_t i3=0; i3<N3; i3++) cache1[0][0][0][i3] = 1;
            perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, 1, 1, arr->local_coulomb_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
            double convolution_mean = 0; for (size_t i3=0; i3<N3; i3++) convolution_mean += cache2[0][0][0][i3]; convolution_mean /= N3;
            if (sys->debug_level>=2) fprintf(sys->log(), "  1*Heaviside(%g-r) = %g\n", sys->rlocal_coul?sys->rcoul:sys->rlocal_coul, convolution_mean);

            for (size_t i3=0; i3<N3; i3++) cache1[0][0][0][i3] = arr->ucoulsr[0][0][i3] + arr->ucoullr[0][0][i3];
            perform_3rx1k_convolution(&arr->fftw_mp, cache1, arr->nx, arr->ny, arr->nz, arr->box, 1, 1, arr->local_coulomb_kernel, arr->dk_nhkvv, sys->xvv_k_shift, arr->n_nhkvv, cache2, arr->fftin, arr->fftout, arr->planf, arr->planb, true);
            for (size_t i3=0; i3<N3; i3++) cache2[0][0][0][i3] = cache1[0][0][0][i3] - cache2[0][0][0][i3]/convolution_mean;
            */

            if (cmd->command==IETCMD_TEST){
                for (int is=0; is<sys->nas; is++){
                    int ix = (sys->pbc_x? ff_pbc_d(ia[is].r.x, arr->box.x) : ia[is].r.x) / arr->box.x * arr->nx;
                    int iy = (sys->pbc_y? ff_pbc_d(ia[is].r.y, arr->box.y) : ia[is].r.y) / arr->box.y * arr->ny;
                    int iz = (sys->pbc_z? ff_pbc_d(ia[is].r.z, arr->box.z) : ia[is].r.z) / arr->box.z * arr->nz;
                    if (ix>=0&&ix<arr->nx && iy>=0&&iy<arr->ny && iz>=0&&iz<arr->nz){
                        if (cmd->command_params_int[iparam]==IETCMD_v_ucoul || cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
                            fprintf(sys->log(), "#%6d %6s %6s q=%8.3g (%4d %4d %4d) coul=%12.2f\n", is+1, as[is].mole, as[is].name, as[is].charge, ix, iy, iz, arr->ucoulsr[iz][iy][ix]+arr->ucoullr[iz][iy][ix]);
                        } else {
                            fprintf(sys->log(), "#%6d %6s %6s q=%8.3g (%4d %4d %4d) coul=%12.2f hlr=%12.2f\n", is+1, as[is].mole, as[is].name, as[is].charge, ix, iy, iz, arr->ucoulsr[iz][iy][ix]+arr->ucoullr[iz][iy][ix], arr->hlr[0][iz][iy][ix]);
                        }

                    }
                }
            } else if (cmd->command==IETCMD_TEST_SAVE){
                if (cmd->command_params_int[iparam]==IETCMD_v_ucoul){
                  // save sys->(ucoulsr + ucoullr)
fprintf(sys->log(), "\33[31msave Coulomb\n\33[0m");
                } else if (cmd->command_params_int[iparam]==IETCMD_v_Yukawa){
fprintf(sys->log(), "\33[31msave Yukawa\n\33[0m");
                  // save sys->ulr
                } else if (cmd->command_params_int[iparam]==IETCMD_v_hlr){
fprintf(sys->log(), "\33[31msave hlr\n\33[0m");
                  // save sys->hlr
                }
            }
        } else if (cmd->command_params_int[iparam]==IETCMD_v_LocalCoulomb){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_sub_test(local-Coulomb)\n");

        } else {
            fprintf(sys->log(), " %s%s : error : perform_sub_test(%d): unrecognized command%s\n", sys->is_log_tty?color_string_of_error:"", software_name, cmd->command_params_int[iparam], sys->is_log_tty?color_string_end:"");
        }
    }

    return ;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void generate_batch_command_1(char * buffer, int len_buffer, IET_Param * sys, IET_arrays * arr, int argc, char * argv[]){

}

void search_solutes_and_run(IET_Param * sys, IET_arrays * arr, int argc, char * argv[]){
  // szfn_solute and szfn_xtc: both exist and confirmed to be directories
  // scan the solute folder
    struct dirent ** file_list_entry;
    int n_files_in_solute = scandir(szfn_solute, &file_list_entry, NULL, alphasort);

    if (n_files_in_solute<=0){
        if (file_list_entry) free(file_list_entry); return ;
    }

  // allocate memory for command line
    int command_line_buffer_length = MAX_PATH*2;
    for (int i=0; i<argc; i++) command_line_buffer_length += strlen(argv[i]) + 1;
    sys->len_batch_cmd = command_line_buffer_length;
    sys->len_batch_out = 4*1024*1024; // 4 MB screen buffer at the most

  #ifdef _LOCALPARALLEL_
    for (int i=0; i<sys->ntb; i++){
        sys->batch_cmd[i] = (char *)memalloc(sys->len_batch_cmd);
        sys->batch_out[i] = (char *)memalloc(sys->len_batch_out);
    }
    char * command_line_buffer = sys->batch_cmd[0];
  #else
    char * command_line_buffer = (char *)memalloc(command_line_buffer_length);
  #endif
    //printf("command_line_buffer_length = %d\n", command_line_buffer_length);

  // create threads or forks
    int sys_nt = sys->nt; sys->nt = sys->ntb;
  #ifdef _LOCALPARALLEL_
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: create_subroutines(%s=%d)\n", sys->mp_by_fork?"np":"nt", sys->nt);
    create_subroutines(sys, arr);
    if (sys->debug_level>=0){
        if (sys->nt<=1){
            fprintf(sys->log(), "%s : master process %d (nice %d)\n", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]));
        } else if (sys->mp_by_fork){
            fprintf(sys->log(), "%s : master process %d (nice %d), %d process%s:", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]), sys->nt, sys->nt>1?"es":"");
            for (int i=0; i<sys->nt; i++) fprintf(sys->log(), " %d", __forkpid[i]);
            fprintf(sys->log(), "\n");
        } else {
            fprintf(sys->log(), "%s : master process %d (nice %d), %d thread%s\n", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]), sys->nt, sys->nt>1?"s":"");
        }
    }
  #else
    fprintf(sys->log(), "%s : master process %d (nice %d)\n", software_name, sys->pid, getpriority(PRIO_PROCESS, sys->pid));
  #endif
    sys->nt = sys_nt;

  // run all tasks

    bool mp_memory_cleared[MAX_THREADS]; for (int i=0; i<MAX_THREADS; i++) mp_memory_cleared[i] = true;

    for (int i=0; i<n_files_in_solute; i++){
        StringNS::string fn = file_list_entry[i]->d_name;
        StringNS::string ext = file_extension(fn);

        if (ext=="solute" || ext=="prmtop"){
            StringNS::string fn_no_ext = file_without_extension(fn);
                char filename_no_ext[MAX_PATH]; memset(filename_no_ext, 0, sizeof(filename_no_ext)); memcpy(filename_no_ext, fn_no_ext.text, fn_no_ext.length);

            char solute_name[MAX_PATH]; memset(solute_name, 0, sizeof(solute_name));
            bool append_solute_path = !szfn_solute[0]? false : StringNS::string(szfn_solute)=="."? false : true;
            snprintf(solute_name, sizeof(solute_name), "%s%s%s.%s", append_solute_path?get_second_fn(szfn_solute):"", append_solute_path?"/":"", filename_no_ext, ext.text);

            char xtc_name[MAX_PATH]; memset(xtc_name, 0, sizeof(xtc_name));
            bool append_xtc_path = !szfn_xtc[0]? false : StringNS::string(szfn_xtc)=="."? false : true;

            if (!xtc_name[0]) snprintf(xtc_name, sizeof(xtc_name), "%s%s%s.pdb", append_xtc_path?get_second_fn(szfn_xtc):"", append_xtc_path?"/":"", filename_no_ext);
            if (access(xtc_name, F_OK) != -1){
            } else xtc_name[0] = 0;

            if (!xtc_name[0]) snprintf(xtc_name, sizeof(xtc_name), "%s%s%s.gro", append_xtc_path?get_second_fn(szfn_xtc):"", append_xtc_path?"/":"", filename_no_ext);
            if (access(xtc_name, F_OK) != -1){
            } else xtc_name[0] = 0;

            if (!xtc_name[0]) snprintf(xtc_name, sizeof(xtc_name), "%s%s%s.xtc", append_xtc_path?get_second_fn(szfn_xtc):"", append_xtc_path?"/":"", filename_no_ext);
            if (access(xtc_name, F_OK) != -1){
            } else xtc_name[0] = 0;

            if (xtc_name[0] && !is_dir(solute_name) && !is_dir(xtc_name)){
                memset(command_line_buffer, 0, command_line_buffer_length);
                // strcpy(command_line_buffer, "echo ");
                for (int j=0; j<argc; j++){
                    int len_buffer = strnlen(command_line_buffer, command_line_buffer_length);
                    char * current_buffer = &command_line_buffer[len_buffer];

                    StringNS::string key = argv[j];
                    if (key == "-s" || key == "--s" || key == "-solute" || key == "--solute" || key == "solute"){
                        strcat(current_buffer, "-s ");
                        strcat(current_buffer, solute_name);
                        strcat(current_buffer, " ");
                        j++;
                    } else if (key == "-f" || key == "--f" || key == "traj" || key == "-traj" || key == "conf" || key == "-conf" || key == "conformation" || key == "-conformation" || key == "conformations" || key == "-conformations"){
                        strcat(current_buffer, "-f ");
                        strcat(current_buffer, xtc_name);
                        strcat(current_buffer, " ");
                        j++;
                    } else if (key == "-log" || key == "--log" || key == "log"){
                        j++;
                    } else {
                        strcat(current_buffer, key.text);
                        strcat(current_buffer, " ");
                    }
                }
                // printf("command: %s\n", command_line_buffer);

              // run the new command
              // followings are substitution of: system(command_line_buffer);
                #ifdef _LOCALPARALLEL_
                    if (sys->ntb<=1){
                        FILE * fd = popen(command_line_buffer, "r"); if (fd){
                            while (fgets(sys->batch_out[0], sys->len_batch_out, fd)){
                                fprintf(sys->log(), "%s", sys->batch_out[0]);
                                fflush(sys->log());
                            }
                            pclose(fd);
                        }
                    } else {
                      // assign jobs to thread: thread 1 ~ max
                        bool job_assigned = false;
                        if (i+1<n_files_in_solute){     // the last job will be assigned to the main thread
                            for (int j=1; j<sys->ntb && !job_assigned; j++){
                                if (__mp_tasks[j] == MPTASK_NONE && mp_memory_cleared[j]){
                                    strncpy(sys->batch_cmd[j], command_line_buffer, sys->len_batch_cmd);
                                    mp_memory_cleared[j] = false;
                                    __mp_tasks[j] = MPTASK_RUN;
                                    job_assigned = true;
                                }
                            }
                        }
                        if (job_assigned) continue;
                      // assign job to main thread: thread 0
                        if (!job_assigned){ int id = 0;
                            FILE * fd = popen(sys->batch_cmd[id], "r"); if (fd){
                                memset(sys->batch_out[id], 0, sys->len_batch_out);
                                // sprintf(sys->batch_out[id], "id %d :", id);
                                while (true){
                                    size_t len_now = strnlen(sys->batch_out[id], sys->len_batch_out);
                                    if (!fgets(&sys->batch_out[id][len_now], sys->len_batch_out-len_now-1, fd)) break;
                                }
                                //printf("%s", sys->batch_out[id]);
                                pclose(fd);
                            }
                        }
                      // wait for everyone
                        while (true){
                            int nbusy = 0; for (int j=1; j<sys->ntb; j++) if (__mp_tasks[j] != MPTASK_NONE) nbusy++;
                            if (nbusy>0){
                                usleep(100);
                            } else {
                                for (int i=1; i<sys->ntb && !mp_memory_cleared[i]; i++){
                                    fprintf(sys->log(), "%s", sys->batch_out[i]);
                                    mp_memory_cleared[i] = true;
                                }
                                fprintf(sys->log(), "%s", sys->batch_out[0]);
                                fflush(sys->log());
                                break;
                            }
                        }
                    }
                #else
                    FILE * fd = popen(command_line_buffer, "r"); if (fd){
                        while (fgets(sys->batch_out[0], sys->len_batch_out, fd)){
                            fprintf(sys->log(), "%s", sys->batch_out[0]);
                        }
                        pclose(fd);
                    }
                #endif

            }
        }
    }

  // end all subroutines
      sys_nt = sys->nt; sys->nt = sys->ntb;
    #ifdef _LOCALPARALLEL_
      if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: send stop signal to subroutines\n");
      for (int i=1; i<MAX_THREADS; i++) __mp_tasks[i] = MPTASK_TERMINATE;
      //wait_subroutines_end(sys);
      wait_subroutines(sys);
      if (sys->debug_level>=2){ char buffer[1024];
          fprintf(sys->log(), "DEBUG:: %s done in %s sec\n", sys->mp_by_fork?"process":"thread[0]", display_time(__total_timer, buffer));
      }
    #endif
      sys->nt = sys_nt;

  // release memory
    if (file_list_entry) free(file_list_entry);

  // clean up everything and exit
    FILE * flog =sys->log();
    sys_nt = sys->nt; sys->nt = sys->ntb;
    display_memory_cost(flog, _memory_blk_total, _memory_total);
    main_print_tailer(sys, arr, true);
    sys->nt = sys_nt;
  #ifdef _FFTWMPPARALLEL_
    fftw_cleanup_threads();
  #endif
    mem_dispose_all(); lap_timer_alloc_memory();
    if (flog && flog!=stdout && flog!=stderr){ FILE * flog_close = flog; flog = nullptr; fclose(flog_close); }
    if (file_out && file_out!=stdout && file_out!=stderr) fclose(file_out);

    exit(0);
}
