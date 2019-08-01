//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool analysis_topology_repeat_line(char * filename, int nline, int nas, StringNS::string * sl, int nw, int * prepeat_lines, int * prepeat_times, FILE * flog){
    if ((nw>1&&sl[0]=="#repeat") || (nw>2&&sl[0]=="#"&&sl[1]=="repeat")){
        int ibegin = 1; if (sl[0]=="#") ibegin = 2; bool analysis_success = false;
        if (ibegin+3<nw && ((sl[ibegin+1]=="lines"||sl[ibegin+1]=="atoms"||sl[ibegin+1]=="line"||sl[ibegin+1]=="atom") && (sl[ibegin+3]=="times"))){
            int repeat_lines = translate_string_to_number(sl[ibegin+0], -1);
            int repeat_times = translate_string_to_number(sl[ibegin+2], -1);
            if (repeat_lines<=0 || repeat_times<=0){ analysis_success = false;
            } else if (repeat_lines>nas){
                fprintf(flog, "%s : %s[%d] : error: not enough atoms to repeat\n", software_name, get_second_fn(filename), nline);
                return false;
            } else {
               *prepeat_lines = repeat_lines; *prepeat_times = repeat_times;
               return true;
            }
        }
        if (analysis_success){
        } else {
            fprintf(flog, "%s : %s[%d] : syntex error: \"", software_name, get_second_fn(filename), nline);
            for (int i=0; i<nw; i++) fprintf(flog, i==0?"%s":" %s", sl[i].text);
            fprintf(flog, "\"\n");
            return false;
        }
        return true;
    } else return false;
}
void read_solute_ff_expand_memory(IET_Param * sys){
    if (sys->nas >= sys->nasmax){
        sys->nasmax += INITIAL_SOLUTE_ATOMS; SoluteAtomSite* asold = sys->as;
        sys->as = (SoluteAtomSite*) malloc(sizeof(AtomSite) * sys->nasmax);
        memcpy(sys->as, asold, sizeof(AtomSite)*sys->nas);
        free(asold);
    }
}
int read_solute_ff(IET_Param * sys, char * filename, bool allow = 0){
    bool allow_compile = allow;
    if (!filename || !filename[0]){ fprintf(sys->log(), "%s%s : error : no -solute (-s) file%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); return -1; }
    FILE * file = fopen(filename, "r");
    if (!file){ fprintf(sys->log(), "%s%s : error : cannot open top: %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, filename, sys->is_log_tty?color_string_end:""); return -1; }
    StringNS::string sl[MAX_WORD]; char input[4096];

    if (!sys->as){ sys->nasmax = INITIAL_SOLUTE_ATOMS; sys->as = (SoluteAtomSite*) malloc(sizeof(AtomSite) * sys->nasmax); }

    fseek(file, 0, SEEK_SET); int iline = 0; while (fgets(input, sizeof(input), file)){ iline ++;
        int nw = StringNS::analysis_line(input, sl, MAX_WORD, true); if (nw<1) continue;
      // preprocessing
        if (sl[0].text[0] == '#'){
            if ((nw>1&&sl[0]=="#repeat") || (nw>2&&sl[0]=="#"&&sl[1]=="repeat")){
                int repeat_lines = 0; int repeat_times = 0;
                if (!analysis_topology_repeat_line(filename, iline, sys->nas, sl, nw, &repeat_lines, &repeat_times, sys->log())){
                    sys->nas = 0; break;
                } else if (repeat_lines<=0 || repeat_times<=0) {
                    sys->nas = 0; break;
                } else {
                    int repeat_begin = sys->nas - repeat_lines; int repeat_end = sys->nas;
                    for (int it=0; it<repeat_times; it++){
                        for (int ir=repeat_begin; ir<repeat_end; ir++){
                            //sys->as[sys->nas++].init(sys->nas, 0, sys->as[ir].mole, 0, sys->as[ir].name, sys->as[ir].mass, sys->as[ir].charge, sys->as[ir].sigma, sys->as[ir].epsilon);
                            //sys->as[sys->nas++].init(sys->nas+1, sys->as[ir].name, sys->as[ir].iaa, sys->as[ir].mole, sys->as[ir].mass, sys->as[ir].charge, sys->as[ir].sigma, sys->as[ir].epsilon);
                            memcpy(&sys->as[sys->nas], &sys->as[ir], sizeof(sys->as[sys->nas]));
                            sys->as[sys->nas].index = sys->as[ir].index + repeat_lines*(it+1);
                            sys->as[sys->nas].iaa = sys->as[ir].iaa + (sys->as[repeat_end-1].iaa-sys->as[repeat_begin].iaa+1)*(it+1);
                            sys->nas++;
                            read_solute_ff_expand_memory(sys);
                        }
                    }
                    if (sys->debug_level>=1) fprintf(sys->log(), "%s : %s[%d] : repeat %d line%s %d times\n", software_name, get_second_fn(filename), iline, repeat_lines, repeat_lines>1?"s":"", repeat_times);
                        //fprintf(sys->log(), "REPEAT %d lines %d times, range: %d %d\n", repeat_lines, repeat_times, repeat_begin, repeat_end);
                }
            }
            continue;
      // comments
        } else if (sl[0].text[0] == ';' || (sl[0].text[0] == '/' && sl[0].text[1] == '/')){
            continue;
      // section selector
        } else if ((sl[0].text[0] == '[' && sl[1] == "solute") || sl[0] == "[solute]"){
            allow_compile = true;
        } else if ((sl[0].text[0] == '[' && sl[1] == "atom") || sl[0] == "[atom]"){
            allow_compile = true;
        } else if (sl[0].text[0] == '['){
            allow_compile = false;
      // atom/bond content reading
        } else if (allow_compile){
            int nwa = 0; for (nwa=1; nwa<nw; nwa++) if (sl[nwa].text[0]=='#' || sl[nwa].text[0]==';' || (sl[nwa].text[0]=='/'&&sl[nwa].text[1]=='/')) break;
          // atom with index and iaa: 8 or more words
            if (nwa>=8){ // line: index atom iaa mole mass charge sigma epsilon [bond]
                SoluteAtomSite * sas = &sys->as[sys->nas];
                sas->init(atoi(sl[0].text), sl[1].text, atoi(sl[2].text), sl[3].text, atof(sl[4].text), atof(sl[5].text), atof(sl[6].text), atof(sl[7].text));
                for (int ipm=8; ipm<nwa; ipm++){
                    if (sl[ipm].length>5&&StringNS::string(sl[ipm].text, 5)=="bond:"){
                        StringNS::string bline = StringNS::string(&sl[ipm].text[5], sl[ipm].length-5);
                        StringNS::string bl[MAX_BONDS_PER_ATOM]; int nbl = StringNS::analysis_csv_line(bline, bl, MAX_BONDS_PER_ATOM, true);
                        for (sas->nbond=0; sas->nbond<nbl; sas->nbond++) sas->ibond[sas->nbond] = atoi(bl[sas->nbond].text);
                    }
                }
                sys->nas++;
                //sys->as[sys->nas++].init(atoi(sl[1].text), 0, sl[2].text, atoi(sl[3].text), sl[0].text, atof(sl[4].text), atof(sl[5].text), atof(sl[6].text), atof(sl[7].text));
                read_solute_ff_expand_memory(sys);
          // basic atom line: only 6 words (fixed)
            } else if (nwa==6){ // line: atom mole mass charge sigma epsilon
                sys->as[sys->nas++].init(sys->nas+1, sl[0].text, 1, sl[1].text, atof(sl[2].text), atof(sl[3].text), atof(sl[4].text), atof(sl[5].text));
                //sys->as[sys->nas++].init(sys->nas, 0, sl[1].text, 0, sl[0].text, atof(sl[2].text), atof(sl[3].text), atof(sl[4].text), atof(sl[5].text));
                read_solute_ff_expand_memory(sys);
            }
        }
    }

//for (int i=0; i<sys->nas; i++) printf("%5d : %s(%d).%s: %12f %12f %12f\n", i, sys->as[i].mole, sys->as[i].iaa, sys->as[i].name, sys->as[i].charge, sys->as[i].sigma, sys->as[i].epsilon);
    if (sys->nas<=0) fprintf(sys->log(), "%s : no atom is loaded from: %s\n", software_name, filename);
    fclose(file); return sys->nas;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bool read_solvent_gvv(IET_Param * sys, IET_arrays * arr, char * filename, int nv, int mv){
    bool success = true; FILE * flog = sys->log();
    if (!filename || !filename[0]) return false;
    FILE * file = fopen(filename, "r"); if (!file){ fprintf(flog, "%s : cannot open %s: %s\n", software_name, sys->gvv_specification>0?"gvv":"hvv", get_second_fn(filename)); return false; }
    int max_gvv_map_col = 0; for (int i=0; i<sys->n_gvv_map; i++) if (sys->gvv_map[i].col > max_gvv_map_col) max_gvv_map_col = sys->gvv_map[i].col;

//printf("gvv_specification= %d\n", sys->gvv_specification);
  // read through to calculate line number
    char input[4096]; StringNS::String sl[MAX_SOL*2];
    int nline = 0; int nlinex = 0; int iline = 0;
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
//printf("gvv: %d: %d %d %d %s\n", nline, nw, sys->n_gvv_map, max_gvv_map_col, filename);
        if (nw < max_gvv_map_col) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        nline ++;
    };
    if (nline<=0){ fprintf(flog, "%s : %s file %s does not contain any data\n", software_name, sys->gvv_specification>0?"gvv":"hvv", get_second_fn(filename)); return false; }
    nlinex = nline * (sys->xvv_extend>1? sys->xvv_extend : 1);
    if (arr->gvv_data) free(arr->gvv_data); arr->gvv_data = init_tensor3d<__REAL__>(1, sys->n_gvv_map+1, nlinex, 0);

  // read all gvv data
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
        if (nw < max_gvv_map_col) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        for (int i=0; i<sys->n_gvv_map; i++) if (sys->gvv_map[i].col-1 >= 0 && sys->gvv_map[i].col-1 < nw){
            arr->gvv_data[0][i][iline] = atof(sl[sys->gvv_map[i].col-1].text) + (sys->gvv_specification<0? 1 : 0);
        }
        iline ++;
    };
    for (int i=0; i<sys->n_gvv_map; i++) for (int j=nline; j<nlinex; j++) arr->gvv_data[0][i][j] = 1;
    arr->n_gvv = nlinex;
//for (int i=0; i<sys->n_gvv_map; i++) printf("map: gvv_data[%d] = data[%d]\n", i, sys->gvv_map[i].col-1); printf("gvv data:"); for (int i=0; i<iline; i++){ for(int j=0; j<sys->n_gvv_map+1; j++) printf(j<sys->n_gvv_map?"%12f ":"%12f ", arr->gvv_data[0][j][i]); printf("\n"); }

    fclose(file);
  // assign data to gvv
    arr->gvv = init_tensor3d_pointer(nv, mv, nlinex);
    for (int i=0; i<nv; i++) for (int j=0; j<mv; j++) arr->gvv[i][j] = &arr->gvv_data[0][sys->n_gvv_map][0];
    for (int i=0; i<sys->n_gvv_map; i++){
        if (sys->gvv_map[i].grpi-1<nv && sys->gvv_map[i].grpj-1<mv) arr->gvv[sys->gvv_map[i].grpi-1][sys->gvv_map[i].grpj-1] = &arr->gvv_data[0][i][0];
        if (sys->gvv_map[i].grpj-1<nv && sys->gvv_map[i].grpi-1<mv) arr->gvv[sys->gvv_map[i].grpj-1][sys->gvv_map[i].grpi-1] = &arr->gvv_data[0][i][0];
        //printf("%d assign to %d,%d\n", i, sys->gvv_map[i].grpi, sys->gvv_map[i].grpj);
    }
    //for (int i=0; i<nv; i++) for (int j=0; j<mv; j++) if (arr->gvv[i][j] == &arr->gvv_data[0][sys->n_gvv_map][0]) printf("not assigned: %d,%d\n", sys->gvv_map[i].grpi, sys->gvv_map[i].grpj);
    arr->dk_gvv = arr->drrism>0&&arr->n_gvv>0? PI / (1+arr->n_gvv) / arr->drrism : 1;
    /*for (int i=0; i<nv; i++) for (int j=i; j<mv; j++) if (arr->gvv[i][j] == &arr->gvv_data[0][sys->n_gvv_map][0]){
        fprintf(sys->log(), "%s : error : gvv_map[%s.%s][%s.%s] not defined\n", software_name, sys->av[i].mole, sys->av[i].name, sys->av[j].mole, sys->av[j].name);
        success = false;
    }*/
//for (int k=0;k<iline; k++){ for (int i=0; i<nv; i++) for(int j=0; j<mv; j++) printf("%9.6f ", arr->gvv[i][j][k]); printf("\n"); }
    return success;
}
bool generate_solvent_xvv(IET_Param * sys, IET_arrays * arr, char * filename, int nv, int mv, int nline, double xvv_extend, bool b_g_wvv=true, bool b_g_nhkvv=true){
  // prepare
    FILE * flog = sys->log();
    bool istty = sys->is_log_tty;
    double * fftin  = (double*) malloc(sizeof(double)*nline);
    double * fftout = (double*) malloc(sizeof(double)*nline);
    fftw_plan plan = fftw_plan_r2r_1d(nline, fftin, fftout, (fftw_r2r_kind)FFTW_RODFT00, FFTW_ESTIMATE);
    double dr = sys->drrism;
    double dk = PI / (1+nline) / dr;   //printf("drrism: %g\n", sys->drrism);
    double factor1 = 2 * PI * dr;

    char report_buffer[256]; memset(report_buffer, 0, sizeof(report_buffer));
  // calculate wvv according to bond list
    if (b_g_wvv){
        int xvv_length = nline;
        arr->n_wvv = xvv_length; arr->dk_wvv = dk;
        __REAL__ *** xvv = init_tensor3d<__REAL__>(nv, nv, xvv_length+2);
        clear_tensor3d(xvv, nv*nv*(xvv_length+2));
        arr->wvv  = xvv; arr->wvv_hlr = arr->wvv;
        for (int iv=0; iv<nv && iv<mv; iv++) for (int i=0; i<xvv_length+2; i++) xvv[iv][iv][i] = 1;
        for (int ib=0; ib<sys->n_bond_list; ib++){
            int addri = -1; int addrj = -1;
            for (int i=0; i<nv; i++) if (sys->av[i].grp == sys->bond_list[ib].grpi){ addri = i; break; };
            for (int i=0; i<nv; i++) if (sys->av[i].grp == sys->bond_list[ib].grpj){ addrj = i; break; };
            double bond = sys->bond_list[ib].bond;      //printf("bond: %d %d %g\n", addri, addrj, bond);
            double bond_sigma = sys->bond_list[ib].bond_stdev;
            if (addri>=0 && addrj>=0 && addri<nv && addrj<nv){
                for (int i=0; i<xvv_length+2; i++){
                    double k = i<xvv_length? (i+1) * dk : 0; double ksigma = k * bond_sigma;
                    xvv[addri][addrj][i] += ((k*bond == 0 ? 1 : (sin(k*bond) / (k*bond))) / sys->av[addri].multi) * exp(-ksigma*ksigma/2);
                    xvv[addrj][addri][i] += ((k*bond == 0 ? 1 : (sin(k*bond) / (k*bond))) / sys->av[addrj].multi) * exp(-ksigma*ksigma/2);
                }
            }
        };
        strcat(report_buffer, " wvv");
        if (sys->n_bond_list<=0 && sys->nv>1) fprintf(flog, "%s : warning : [bond] section undefined.\n", software_name);
        //if (sys->detail_level>=1) fprintf(flog, "%s : wvv generated from %s%s%s\n", software_name, istty?prompt_path_prefix:"", get_second_fn(filename), istty?"\33[0m":"");
        //for (int k=0;k<xvv_length; k++){ for (int i=0; i<nv; i++) for(int j=0; j<mv; j++) printf("%9.6f ", xvv[i][j][k]); printf("\n"); }
    } else arr->n_wvv = 0;

  // calculate nhkvv according to gvv
    if (b_g_nhkvv){
        int xvv_length = nline;
        arr->n_nhkvv = xvv_length; arr->dk_nhkvv = dk;
        arr->nhkvv = init_tensor3d<__REAL__>(nv, nv, xvv_length+2); arr->nhkvv_hlr = arr->nhkvv;
        clear_tensor3d(arr->nhkvv, nv*nv*(xvv_length+2));
        for (int iv=0; iv<nv; iv++) for (int jv=0; jv<mv; jv++){
            if (arr->gvv[iv][jv] == &arr->gvv_data[0][sys->n_gvv_map][0]){
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : gvv_map[%s.%s][%s.%s] set to zero\n", software_name, sys->av[iv].mole, sys->av[iv].name, sys->av[jv].mole, sys->av[jv].name);
            } else {
                double density = sys->bulk_density_mv[sys->av[iv].iaa] * sys->av[jv].multi;
                //printf("density %d %d : %12f\n", iv, jv, sys->density_av[iv]);
                arr->nhkvv[iv][jv][xvv_length+1] = 0;
                for (int i=0; i<nline; i++){
                    fftin[i] = (arr->gvv[iv][jv][i] - 1) * (i+1)*dr;
                    arr->nhkvv[iv][jv][xvv_length+1] += 2 * fftin[i] * (i+1)*dr * factor1 * density;
                        // additional factor 2 : refer to "FFTW 3.3.8: What FFTW Really Computes / 4.8.4 1d Real-odd DFTs (DSTs)"
                }
                arr->nhkvv[iv][jv][xvv_length] = arr->nhkvv[iv][jv][xvv_length+1];
                fftw_execute(plan);
                for (int i=0; i<xvv_length; i++){
                    arr->nhkvv[iv][jv][i] += fftout[i] * factor1 / ((i+1)*dk) * density;
                }
            }
        }
        strcat(report_buffer, " nhkvv");
        //if (sys->detail_level>=1) fprintf(flog, "%s : nhkvv generated from %s%s%s\n", software_name, istty?prompt_path_prefix:"", get_second_fn(filename), istty?"\33[0m":"");
        //for (int k=0;k<nline; k++){ printf("%3d ", k+1); for (int i=0; i<nv; i++) for(int j=0; j<mv; j++) printf("%9.6f ", arr->nhkvv[i][j][k]); printf("\n"); }
    } else arr->n_nhkvv = 0;

    //printf("wvv[0]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->wvv[i][j][arr->n_wvv+1]); printf("\n"); } printf("wvv[1]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->wvv[i][j][0]); printf("\n"); } printf("wvv[2]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->wvv[i][j][1]); printf("\n"); } printf("nhk[0]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->nhkvv[i][j][arr->n_nhkvv+1]); printf("\n"); } printf("nhk[1]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->nhkvv[i][j][0]); printf("\n"); } printf("nhk[2]:\n"); for (int i=0; i<nv; i++){ for (int j=0; j<nv; j++) printf(" %9.3g", arr->nhkvv[i][j][1]); printf("\n"); }

  // build Yukawa FFT kernel
    if (arr->n_nhkvv>0){
        int xvv_length = arr->n_nhkvv;

        arr->yukawa_kernel = init_tensor3d<__REAL__>(1, 1, xvv_length);
        double alpha = 1/sys->rc_yukawafft; double alpha2 = alpha*alpha;
        for (int i=0; i<xvv_length; i++){
            double k = (i)*dk; double k2 = k*k;
            arr->yukawa_kernel[0][0][i] = k2 / (k2 + alpha2);

            //double k = (i+1)*dk; double k2 = k*k;
            //arr->yukawa_kernel[0][0][i] = k2 / (k2 + alpha2);
        }

        arr->local_coulomb_kernel = init_tensor3d<__REAL__>(1, 1, xvv_length);
        double rlocal = sys->rlocal_coul<=0? sys->rcoul : sys->rlocal_coul;
        double local_coulomb_kernel_V = (4*PI/3*rlocal*rlocal*rlocal);
        for (int i=0; i<xvv_length; i++){
            double k = (i+0)*dk;
            double ksigma = k*rlocal;
            arr->local_coulomb_kernel[0][0][i] = k==0? 1 : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k) / local_coulomb_kernel_V;
        }

        arr->ld_kernel = init_tensor3d<__REAL__>(sys->nvm, sys->nvm, xvv_length);
        clear_tensor3d(arr->ld_kernel, sys->nvm*sys->nvm*xvv_length);
        for (int ivm=0; ivm<sys->nvm; ivm++){
            //double density = sys->bulk_density_mv[sys->av[ivm].iaa];
            double density = sys->bulk_density_mv[ivm];
            double rmol = pow(fabs(1.0/sys->bulk_density_mv[ivm]/4.0*3.0/PI), 1.0/3);

            /*
            for (int i=0; i<xvv_length; i++) fftin[i] = (i*dr<=rmol? 1 : 0) * (i+1)*dr;
            fftw_execute(plan);
            for (int i=0; i<xvv_length-1; i++) arr->ld_kernel[ivm][ivm][i+1] = fftout[i] * factor1 / ((i+1)*dk) * density;
            arr->ld_kernel[ivm][ivm][0] = 1;
            //*/

            ///*
            for (int i=0; i<xvv_length; i++){
                //double k = (i+1)*dk; if (k<MACHINE_REASONABLE_ERROR) k = MACHINE_REASONABLE_ERROR;
                double k = i*dk;
                double ksigma = k*rmol;
                arr->ld_kernel[ivm][ivm][i] = k==0? 1 : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k)*density;
            }
            //*/
        }
    }

  // finalization
    fftw_destroy_plan(plan); free(fftin); free(fftout);
    if (strlen(report_buffer)>0 && sys->detail_level>=2) fprintf(flog, "%s :%s generate from %s%s%s\n", software_name, report_buffer, istty?prompt_path_prefix:"", get_second_fn(filename), istty?prompt_path_suffix:"");

    return true;
}
void backup_solvent_xvv_with_shift_only(IET_Param * sys, IET_arrays * arr){
    arr->wvv_save = init_tensor3d<__REAL__>(sys->nv, sys->nv, arr->n_wvv+2);
    for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++){
        arr->wvv_save[iv][jv][0] = arr->wvv[iv][jv][arr->n_wvv+1];
        for (int k=0; k<arr->n_wvv; k++) arr->wvv_save[iv][jv][k+1] = arr->wvv[iv][jv][k];
    }
    arr->nhkvv_save = init_tensor3d<__REAL__>(sys->nv, sys->nv, arr->n_nhkvv+2);
    for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++){
        arr->nhkvv_save[iv][jv][0] = arr->nhkvv[iv][jv][arr->n_nhkvv+1];
        for (int k=0; k<arr->n_nhkvv; k++) arr->nhkvv_save[iv][jv][k+1] = arr->nhkvv[iv][jv][k];
    }

    /*MatrixNS::Matrix A, Ai, T; A.init(sys->nv); Ai.init(sys->nv); T.init(sys->nv);
    arr->xvv_inverse = init_tensor3d<__REAL__>(sys->nv, sys->nv, arr->n_nhkvv);
    for (int k=0; k<arr->n_nhkvv; k++){
        for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++) *A.e(iv, jv) = arr->wvv_save[iv][jv][k] + arr->nhkvv_save[iv][jv][k];
        A.inverse(Ai, T);
        for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++) arr->xvv_inverse[iv][jv][k] = *Ai.e(iv, jv);
    }
    A.dispose(); Ai.dispose(); T.dispose();

    FILE * flog = sys->log();
    int n_nhkvv = arr->n_nhkvv; int n_wvv = arr->n_wvv;
    for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++){
        if (sys->debug_level>=1){
            fprintf(flog, "debug:: (xvv_inverse)[%d][%d] = (%.3f)", iv+1, jv+1, arr->xvv_inverse[iv][jv][0]);
            if (n_wvv>1 && n_nhkvv>1) fprintf(flog, ", (%.3f)", arr->xvv_inverse[iv][jv][1]);
            if (n_wvv>2 && n_nhkvv>2) fprintf(flog, ", (%.3f)", arr->xvv_inverse[iv][jv][2]);
            if (n_wvv>3 && n_nhkvv>3) fprintf(flog, ", (%.3f)", arr->xvv_inverse[iv][jv][3]);
            if (n_wvv>4 && n_nhkvv>4) fprintf(flog, ", (%.3f)", arr->xvv_inverse[iv][jv][4]);
            fprintf(flog, (n_wvv>4&&n_nhkvv>4)? " ... [%s.%s][%s.%s]\n":" [%s.%s][%s.%s]\n", sys->av[iv].mole, sys->av[iv].name, sys->av[jv].mole, sys->av[jv].name);
        }
    }*/
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void debug_display_solvent_xvv(IET_Param * sys, IET_arrays * arr, int nv){
    double dr = sys->drrism; FILE * flog = sys->log();
  // display at deep debug mode
        fprintf(sys->log(), "%s : DEBUG : wvv =\n", software_name);
        for (int i=0; i<arr->n_wvv; i++){
            fprintf(sys->log(), "  %9.4g", arr->dk_wvv*(i)/PI);
            for (int iv=0; iv<nv; iv++) for (int jv=0; jv<nv; jv++) fprintf(sys->log(), " %9.4f", arr->wvv[iv][jv][i]);
            fprintf(sys->log(), "\n");
        }
  // display at deep debug mode
        fprintf(sys->log(), "%s : DEBUG : nhkvv =\n", software_name);
        for (int i=0; i<arr->n_nhkvv; i++){
            fprintf(sys->log(), "  %9.4g", arr->dk_nhkvv*(i)/PI);
            for (int iv=0; iv<nv; iv++) for (int jv=0; jv<nv; jv++) fprintf(sys->log(), " %9.4f", arr->nhkvv[iv][jv][i]);
            fprintf(sys->log(), "\n");
        }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int __read_solvent_wvv(IET_Param * sys, char * filename, IET_arrays * arr, int nv){
    if (!filename || !filename[0]) return -1;
    FILE * file = fopen(filename, "r"); if (!file) return -2;

    char input[4096]; StringNS::String sl[MAX_SOL*2];
    int nline = 0; int iline = 0;
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
//printf("wvv: %d %d\n", nw, nv);
        if (nw<nv*nv) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        nline ++;
    };
    arr->wvv = init_tensor3d<__REAL__>(nv, nv, nline);
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
        if (nw<nv*nv) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        if (sys->transpose_vv){
            for (int i=0; i<nw && i<nv*nv && iline<nline; i++) arr->wvv[i%nv][i/nv][iline] = atof(sl[i].text);
        } else {
            for (int i=0; i<nw && i<nv*nv && iline<nline; i++) arr->wvv[i/nv][i%nv][iline] = atof(sl[i].text);
        }
        iline ++;
    }; arr->n_wvv = nline;

    fclose(file); return iline;
}
int __read_solvent_nhkvv(IET_Param * sys, char * filename, IET_arrays * arr, int nv, int nv2){
  // nv: sys->nv; nv2: already existing lattice
    if (!filename || !filename[0]) return -1;
    FILE * file = fopen(filename, "r"); if (!file) return -2;

    char input[4096]; StringNS::String sl[MAX_SOL*2];
    int nline = 0; int iline = 0;
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
        if (nw<nv*nv2) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        nline ++;
    };
    arr->nhkvv = init_tensor3d<__REAL__>(nv2, nv, nline);
    fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
        int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
        if (nw<nv*nv2) continue;
        if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
        if (sys->transpose_vv){
            for (int i=0; i<nw && i<nv*nv2 && iline<nline; i++) arr->nhkvv[i%nv][i/nv][iline] = atof(sl[i].text);
        } else {
            for (int i=0; i<nw && i<nv*nv2 && iline<nline; i++) arr->nhkvv[i/nv][i%nv][iline] = atof(sl[i].text);
        }
        iline ++;
    }; arr->n_nhkvv = nline;

    fclose(file); return iline;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool generate_solvent_zeta(int & n_zeta, double & dk_zeta, IET_Param * sys, IET_arrays * arr, int nvm, double beta = 1){
    bool success = true;

    n_zeta = arr->n_nhkvv; dk_zeta = arr->dk_nhkvv;
    double density_hi = sys->density_hi; double * bulk_density_mv = sys->bulk_density_mv;

    arr->zeta = init_tensor3d<__REAL__>(nvm, nvm, n_zeta);
    for (int ivm=0; ivm<nvm; ivm++) for (int jvm=0; jvm<nvm; jvm++){
        int ilist = -1; for (int iilist=0; iilist<sys->n_zeta_list; iilist++) if (sys->zeta_list[iilist].iaai==ivm && sys->zeta_list[iilist].iaaj==jvm) ilist = iilist;
        if (ilist>=0){
            double rc_zeta0 = sys->zeta_list[ilist].rc_zeta0;
            double local_kernel_V = (4*PI/3*rc_zeta0*rc_zeta0*rc_zeta0);
            double zeta0 = sys->zeta_list[ilist].zeta0 * beta;
            for (int ik=0; ik<n_zeta; ik++){
                double k = ik * dk_zeta;
                double ksigma = k*rc_zeta0;
                arr->zeta[ivm][jvm][ik] = (k==0? local_kernel_V : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k)) * zeta0;
            }

            if (sys->debug_level>=1){
                double kT_kJ_mol = sys->temperature / 120.27;
                double kT_kcal_mol = sys->temperature / 120.27 / 4.184;
                fprintf(sys->log(), "debug:: rho[%d]x<zeta[%d][%d]> = %.4g x %.4g = %.4g kJ/mol = %.4g kcal/mol\n", jvm, ivm, jvm, sys->bulk_density_mv[jvm], arr->zeta[ivm][jvm][0], arr->zeta[ivm][jvm][0]*sys->bulk_density_mv[jvm]*kT_kJ_mol, arr->zeta[ivm][jvm][0]*sys->bulk_density_mv[jvm]*kT_kcal_mol);
            }

        } else {
            if (sys->zeta_list_allow_missing) for (int ik=0; ik<n_zeta; ik++) arr->zeta[ivm][jvm][ik] = 0;
            int indexi = -1; int indexj = -1; for (int i=0; i<sys->nav && (indexi<0||indexj<0); i++){ if (sys->av[i].iaa==ivm) indexi = i; if (sys->av[i].iaa==jvm) indexj = i; }
            if (sys->zeta_list_allow_missing){
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : zeta[%s][%s] set to zero\n", software_name, indexi>=0?sys->av[indexi].mole:"mol_unknown", indexj>=0?sys->av[indexj].mole:"mol_unknown");
            } else {
                fprintf(sys->log(), "%s%s : error : zeta[%s][%s] undefined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, indexi>=0?sys->av[indexi].mole:"mol_unknown", indexj>=0?sys->av[indexj].mole:"mol_unknown", sys->is_log_tty?color_string_end:"");
                success = false;
            }
        }
    }

  // renorm zeta with density
    for (int ivm=0; ivm<nvm; ivm++) for (int jvm=0; jvm<nvm; jvm++){
        double old_mean_zetak = arr->zeta[ivm][jvm][0];
        //double factor = bulk_density_mv[ivm] * bulk_density_mv[jvm] / density_hi / density_hi;
        double factor = bulk_density_mv[jvm] / density_hi;
        for (int i=0; i<n_zeta; i++) arr->zeta[ivm][jvm][i] *= factor;
    }

    return success;
}
bool generate_solvent_zeta(IET_Param * sys, IET_arrays * arr){
    return generate_solvent_zeta(arr->n_zeta, arr->dk_zeta, sys, arr, sys->nvm, sys->default_temperature / sys->temperature);
}
bool read_solvent_zeta_fft(__REAL__ *** zeta, int & n_zeta, double & dk_zeta, IET_Param * sys, IET_arrays * arr, char * filename, int nvm, int mvm, int nline, double beta = 1){
    //for (int il=0; il<10; il++) for (int iv=0; iv<nvm; iv++){ for (int jv=0; jv<mvm; jv++){ fprintf(sys->log(), " %11g", zeta[iv][jv][il]); }; fprintf(sys->log(), "\n"); };
  // prepare
    //double beta = sys->default_temperature / sys->temperature;
    FILE * flog = sys->log();
    bool istty = sys->is_log_tty;
    double * fftin  = (double*) malloc(sizeof(double)*nline);
    double * fftout = (double*) malloc(sizeof(double)*nline);
    fftw_plan plan = fftw_plan_r2r_1d(nline, fftin, fftout, (fftw_r2r_kind)FFTW_RODFT00, FFTW_ESTIMATE);
    double dr = sys->drhi;
    double dk = PI / (1+nline) / dr;   //printf("drrism: %g\n", sys->drrism);
    double factor1 = 2 * PI * dr;

  // calculate zkvv according to zeta
    //arr->n_zeta = nline; arr->dk_zeta = dk;
    n_zeta = nline; dk_zeta = dk;
    for (int iv=0; iv<nvm; iv++) for (int jv=0; jv<mvm; jv++){
        for (int i=0; i<nline; i++) fftin[i] = zeta[iv][jv][i] * (i+1)*dr;
        fftw_execute(plan);
        for (int i=0; i<nline; i++) zeta[iv][jv][i] = fftout[i] * factor1 / ((i+1)*dk) * beta;
    }

  // finalization
    fftw_destroy_plan(plan); free(fftin); free(fftout);
    //if (sys->detail_level>=1) fprintf(flog, "%s : zkvv generated from %s%s%s\n", software_name, istty?prompt_path_prefix:"", get_second_fn(filename), istty?"\33[0m":"");
    //for (int il=0; il<10; il++) for (int iv=0; iv<nvm; iv++){ for (int jv=0; jv<mvm; jv++){ fprintf(sys->log(), " %11g", zeta[iv][jv][il]); }; fprintf(sys->log(), "\n"); }; return false;

    return true;
}
bool read_solvent_zeta(IET_Param * sys, IET_arrays * arr, char * fnzeta, int nvm){
    bool success_livm = false;
    arr->drhi = sys->drhi;
  // HI: zeta
    if (fnzeta && fnzeta[0]){
        char * filename = fnzeta; success_livm = true;
        double * bulk_density_mv = sys->bulk_density_mv; double density_hi = sys->density_hi; SoluteAtomSite * as = sys->as;

        FILE * file = fopen(filename, "r"); if (!file){
            fprintf(sys->log(), "%s : cannot open -zeta %s\n", software_name, szfn_zeta); return false;
        }

        char input[4096]; StringNS::String sl[MAX_SOL*2];
        int nline = 0; int iline = 0;
        fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
            int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
            if (nw<nvm*nvm) continue;
            if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
            nline ++;
        };

        arr->zeta = init_tensor3d<__REAL__>(nvm, nvm, nline);

        fseek(file, 0, SEEK_SET); while (fgets(input, sizeof(input), file)){
            int nw = analysis_line_params(input, sl, sizeof(sl)/sizeof(sl[0]), true);
            if (nw<nvm*nvm) continue;
            if (sl[0].text[0]=='#' || sl[0].text[0]=='[' || sl[0].text[0]==';') continue;
            if (sys->transpose_vv){
                for (int i=0; i<nw && i<nvm*nvm && iline<nline; i++) arr->zeta[i%nvm][i/nvm][iline] = atof(sl[i].text);
            } else {
                for (int i=0; i<nw && i<nvm*nvm && iline<nline; i++) arr->zeta[i/nvm][i%nvm][iline] = atof(sl[i].text);
            }
            iline ++;
        }; arr->n_zeta = nline;
      // scale zeta
        if (sys->n_zeta_scaling_factor>0 && sys->n_zeta_scaling_factor>=sys->nvm){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: scale_zeta:");
            for (int ivm=0; ivm<sys->nvm; ivm++) for (int jvm=0; jvm<sys->nvm; jvm++) if (sys->zeta_scaling_factor[ivm]>0 && sys->zeta_scaling_factor[jvm]>0){
                double factor = sqrt(sys->zeta_scaling_factor[ivm] * sys->zeta_scaling_factor[jvm]);
                if (sys->debug_level>=2) fprintf(sys->log(), " %g->", arr->zeta[ivm][jvm][0]);
                for (int i=0; i<arr->n_zeta; i++) arr->zeta[ivm][jvm][i] *= factor;
                if (sys->debug_level>=2) fprintf(sys->log(), "%g", arr->zeta[ivm][jvm][0]);
            }
            if (sys->debug_level>=2) fprintf(sys->log(), "\n");
        }
        if (sys->n_zeta_scaling_factor>0 && sys->n_zeta_scaling_factor<sys->nvm){
            fprintf(sys->log(), "%s%s : error : too few (%d) scaling factors of zeta (%d required).%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->n_zeta_scaling_factor, sys->nvm, sys->is_log_tty?color_string_end:"");
            return false;
        }

      // do fft if not in k-space
        read_solvent_zeta_fft(arr->zeta, arr->n_zeta, arr->dk_zeta, sys, arr, szfn_zeta, sys->nvm, sys->nvm, arr->n_zeta, sys->default_temperature / sys->temperature);

      // show value of zeta integral
        if (sys->debug_level>=1 && bulk_density_mv && as){
            double kT_kJ_mol = sys->temperature / 120.27;
            double kT_kcal_mol = sys->temperature / 120.27 / 4.184;
            //for (int ivm=0; ivm<nvm; ivm++) for (int jvm=0; jvm<nvm; jvm++) fprintf(flog, "%s : debug : rho[%d]x<zeta[%d][%d]> = %.4g x %.4g = %.4g kJ/mol = %.4g kcal/mol\n", software_name, jvm, ivm, jvm, sys->density_mv[jvm], arr->zeta[ivm][jvm][0], arr->zeta[ivm][jvm][0]*sys->density_mv[jvm]*kT_kJ_mol, arr->zeta[ivm][jvm][0]*sys->density_mv[jvm]*kT_kcal_mol);
            for (int ivm=0; ivm<nvm; ivm++) for (int jvm=0; jvm<nvm; jvm++) fprintf(sys->log(), "debug:: rho[%d]x<zeta[%d][%d]> = %.4g x %.4g = %.4g kJ/mol = %.4g kcal/mol\n", jvm, ivm, jvm, sys->bulk_density_mv[jvm], arr->zeta[ivm][jvm][0], arr->zeta[ivm][jvm][0]*sys->bulk_density_mv[jvm]*kT_kJ_mol, arr->zeta[ivm][jvm][0]*sys->bulk_density_mv[jvm]*kT_kcal_mol);
        }

      // renorm zeta with density
        for (int ivm=0; ivm<nvm; ivm++) for (int jvm=0; jvm<nvm; jvm++){
            double old_mean_zetak = arr->zeta[ivm][jvm][0];
    //printf("solvent %d %d: density %12f %12f , density_HI %12f %12f, zeta[%s,%s][0] %12f -> ", ivm, jvm, bulk_density_mv[ivm], bulk_density_mv[jvm], density_hi, density_hi, as[ivm].mole, as[jvm].mole, arr->zeta[ivm][jvm][0]);
            double factor = bulk_density_mv[ivm] * bulk_density_mv[jvm] / density_hi / density_hi;
            for (int i=0; i<nline; i++) arr->zeta[ivm][jvm][i] *= factor;
            //if (sys->debug_level>=1 && bulk_density_mv && as) fprintf(flog, "%s : debug:: <zetak[%s][%s]> = %g <-- %g\n", software_name, as[ivm].mole, as[jvm].mole, arr->zeta[ivm][jvm][0], old_mean_zetak);
        }
//fprintf(stderr, "\33[31m arr->zeta[0][0][0] = %g\n\33[0m", arr->zeta[0][0][0]);

        fclose(file);

        if (iline<=0){
            fprintf(sys->log(), "%s : no valid zeta data in %s (need %d cols)\n", software_name, szfn_zeta, sys->nvm*sys->nvm); success_livm = false;
        }
    } else if (sys->hial>=HIAL_HI){
        fprintf(sys->log(), "%s%s : error : HI is on but zeta is not specified%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:"");
        success_livm = false;
    } else success_livm = true;
    // calculate dk
    arr->dk_zeta = arr->drhi>0&&arr->n_zeta>0? PI / (1+arr->n_zeta) / arr->drhi : 1;
    return success_livm;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// other advanced process for post handling of inputs
void debug_show_rism_xvv_matrix(IET_Param * sys, IET_arrays * arr, int nv, int nvm){
  // display xvv mean value in debug mode
    fprintf(sys->log(), "debug:: wvv[,][0] = \n");
    for (int iv=0; iv<nv; iv++){ fprintf(sys->log(), "  "); for (int jv=0; jv<nv; jv++) fprintf(sys->log(), " %10.4g", arr->wvv[iv][jv][0]); fprintf(sys->log(), "\n"); }
    fprintf(sys->log(), "debug:: nhkvv[,][0] = \n");
    for (int iv=0; iv<nv; iv++){ fprintf(sys->log(), "  "); for (int jv=0; jv<nv; jv++) fprintf(sys->log(), " %10.4g", arr->nhkvv[iv][jv][0]); fprintf(sys->log(), "\n"); }
    fprintf(sys->log(), "debug:: xvv[,][0] = \n");
    for (int iv=0; iv<nv; iv++){ fprintf(sys->log(), "  "); for (int jv=0; jv<nv; jv++) fprintf(sys->log(), " %10.4g", arr->wvv[iv][jv][0] + arr->nhkvv[iv][jv][0]); fprintf(sys->log(), "\n"); }
    if (arr->zeta){
        fprintf(sys->log(), "debug:: zetavv[,][0] = \n");
        for (int iv=0; iv<nvm; iv++){ fprintf(sys->log(), "  "); for (int jv=0; jv<nvm; jv++) fprintf(sys->log(), " %10.4g", arr->zeta[iv][jv][0]); fprintf(sys->log(), "\n"); }
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool calculate_lse_ab_automatically(IET_Param * sys, IET_arrays * arr, bool b_read_solvent_xvv, bool b_read_solvent_zeta){
  #ifdef _EXPERIMENTAL_
    return calculate_lse_ab_eperimental_automatically(sys, arr, b_read_solvent_xvv, b_read_solvent_zeta);
  #else
  // Post Handling
    if (b_read_solvent_zeta){
      // calculate lse_b
        if (sys->calc_ab_automatically && arr->n_zeta>0){
            double lse_abn = 0; double average_density_hi = 0;
            for (int i=0; i<sys->nvm; i++){
                lse_abn += sys->llambda[i] * sys->nbulk[i];
                for (int j=0; j<sys->nvm; j++) lse_abn += - sys->nbulk[i] * sys->nbulk[j] * arr->zeta[i][j][0] * sys->density_hi * sys->density_hi;
                average_density_hi += sys->density_hi / sys->nvm;
            }
            if (sys->calc_ab_automatically==1){
                sys->lse_a = (lse_abn/average_density_hi)/sys->lse_b;
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : lse_a autogened, A=%g, B=%g\n", software_name, sys->lse_a, sys->lse_b);
            } else if (sys->calc_ab_automatically==2){
                sys->lse_b = (lse_abn/average_density_hi)/sys->lse_a;
                if (sys->detail_level>=2) fprintf(sys->log(), "%s : lse_b autogened, A=%g, B=%g\n", software_name, sys->lse_a, sys->lse_b);
            }
            //printf("LSEANB: %12f %12f %12f\n", lse_abn, lse_abn/sys->lse_a/average_density_hi, average_density_hi);
        }
    }
    return true;
//printf("wvv: %d %d nhkvv: %d %d\n", rdvv, arr->n_wvv, rdvv, arr->n_nhkvv); for (int i=0; i<10&&i<arr->n_wvv; i++){ printf(" wvv"); for (int j=0; j<sys->nv; j++) for (int k=0; k<sys->nv; k++) printf(" %7.2f", arr->wvv[j][k][i]); printf(" nhkvv"); for (int j=0; j<sys->nv; j++) for (int k=0; k<sys->nv; k++) printf(" %7.2f", arr->nhkvv[j][k][i]); printf("\n"); }
//printf("wvv: %d %d nhkvv: %d %d\n", rdvv, arr->n_wvv, rdvv, arr->n_nhkvv);
//for (int i=arr->n_wvv-10; i>=0&&i<arr->n_wvv; i++){ printf(" wvv"); for (int j=0; j<sys->nv; j++) for (int k=0; k<sys->nv; k++) printf(" %7.2f", arr->wvv[j][k][i]); printf(" nhkvv"); for (int j=0; j<sys->nv; j++) for (int k=0; k<sys->nv; k++) printf(" %7.2f", arr->nhkvv[j][k][i]); printf("\n"); }
  #endif
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int map_rdf_grps_to_pairs(IET_Param * sys, RDF_data * buffer=nullptr, bool echo_warning=false){
    int npr = 0; bool success = true;

    for (int ig=0; ig<sys->n_rdf_grps; ig++){
        //printf("  RDF GROUP[%d]: %d (%s.%s) - %d (%s.%s)\n", ig, sys->rdf_grps[ig].is, sys->rdf_grps[ig].ms, sys->rdf_grps[ig].as, sys->rdf_grps[ig].iv, sys->rdf_grps[ig].mv, sys->rdf_grps[ig].av);
      // prepare
        StringNS::string ms, as, mv, av;
          ms = sys->rdf_grps[ig].ms; if (ms=="*") ms = "\0";
          as = sys->rdf_grps[ig].as; if (as=="*") as = "\0";
          mv = sys->rdf_grps[ig].mv; if (mv=="*") mv = "\0";
          av = sys->rdf_grps[ig].av; if (av=="*") av = "\0";
      // mark all the solutes selected
        bool mark_all_solute_as_whole = false;
        if ((ms=="all" && as.length==0) || (as=="all" && ms.length==0) || (as=="all" && ms=="all")){
            //fprintf(sys->log(), "%s.%s - %s.%s\n", ms.text, as.text, mv.text, av.text);
            mark_all_solute_as_whole = true;
        } else for (int i=0; i<sys->nas; i++){
            //printf("now mark solute %d\n", i);
            sys->as[i].reserved = 0;
            if (sys->rdf_grps[ig].is==0 || sys->rdf_grps[ig].is==i+1){ // match index
                sys->as[i].reserved = 1;
                //printf("solute %d: mark with ID\n", i);
            } else if (sys->rdf_grps[ig].is<0){ // match name
                bool mole_name_matched = ms.length<1? true : ms==sys->as[i].mole? true : false;
                bool atom_name_matched = as.length<1? true : as==sys->as[i].name? true : false;
                if (mole_name_matched && atom_name_matched) sys->as[i].reserved = 1;
                //printf("solute %d: %s with name: %d, %d (%s %s)\n", i, sys->as[i].reserved==1?"mark":"not marked", mole_name_matched, atom_name_matched, sys->as[i].name, as.text);
            } //printf("solute %d: not marked\n", i);
        }
      // mark all the solvents selected
        int atom_selection[MAX_SOL];
        for (int i=0; i<sys->nv; i++){
            atom_selection[i] = 0;
            if (sys->rdf_grps[ig].iv==0 || sys->rdf_grps[ig].iv==i+1){ // match index
                atom_selection[i] = 1;
            } else if (sys->rdf_grps[ig].iv<0){ // match name
                bool mole_name_matched = mv.length<1? true : mv==sys->av[i].mole? true : false;
                bool atom_name_matched = av.length<1? true : av==sys->av[i].name? true : false;
                if (mole_name_matched && atom_name_matched) atom_selection[i] = 1;
            }
        }
      // gathering all the marked atoms and form pairs
        int n_pair_this_grp = 0;
        if (mark_all_solute_as_whole){
            for (int j=0; j<sys->nv; j++) if (atom_selection[j]){
                //printf("RDF pair: %d (%s.%s) - %d (%s.%s)\n", i+1, sys->as[i].mole, sys->as[i].name, j+1, sys->av[j].mole, sys->av[j].name);
                if (buffer){ buffer[npr].is = -1; buffer[npr].iv = j; }
                n_pair_this_grp ++; npr ++;
            }
        } else {
            for (int i=0; i<sys->nas; i++) for (int j=0; j<sys->nv; j++) if (sys->as[i].reserved && atom_selection[j]){
                //printf("RDF pair: %d (%s.%s) - %d (%s.%s)\n", i+1, sys->as[i].mole, sys->as[i].name, j+1, sys->av[j].mole, sys->av[j].name);
                if (buffer){ buffer[npr].is = i; buffer[npr].iv = j; }
                n_pair_this_grp ++; npr ++;
            }
        }
        //if (n_pair_this_grp<1 && echo_warning){
        if (n_pair_this_grp<1){
            char tmp[2][96]; tmp[0][0] = tmp[1][0] = 0;
            if (sys->rdf_grps[ig].is==0) strncpy(tmp[0], "*", 96); else if (sys->rdf_grps[ig].is>0) snprintf(tmp[0], 96, "%d", sys->rdf_grps[ig].is); else snprintf(tmp[0], 96, "%s%s%s", ms.text, ms.length<1?"":".", as.text);
            if (sys->rdf_grps[ig].iv==0) strncpy(tmp[1], "*", 96); else if (sys->rdf_grps[ig].iv>0) snprintf(tmp[1], 96, "%d", sys->rdf_grps[ig].iv); else snprintf(tmp[1], 96, "%s%s%s", mv.text, mv.length<1?"":".", av.text);
            fprintf(sys->log(), "%s%s : error : no atom pair matches rdf_group[%d]=%s-%s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, ig+1, tmp[0], tmp[1], sys->is_log_tty?color_string_end:"");
            success = false;
        }

      // clear everything
        for (int i=0; i<sys->nas; i++) sys->as[i].reserved = 0;
    }
    return success? npr : -1;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef _EXPERIMENTAL_
    void shift_xvv_by_1(IET_Param * sys, __REAL__ *** wvv, int n_wvv, __REAL__ *** nhkvv, int n_nhkvv){
        for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++){
            for (int i=n_wvv; i>0; i--) wvv[iv][jv][i] = wvv[iv][jv][i-1];
            wvv[iv][jv][0] = wvv[iv][jv][n_wvv+1];
            for (int i=n_nhkvv; i>0; i--) nhkvv[iv][jv][i] = nhkvv[iv][jv][i-1];
            nhkvv[iv][jv][0] = nhkvv[iv][jv][n_nhkvv+1];
        }
    }
    void insert_xvv_at_0(IET_Param * sys, __REAL__ *** wvv, __REAL__ *** nhkvv){
        for (int iv=0; iv<sys->nv; iv++) for (int jv=0; jv<sys->nv; jv++){
            double nhkvv_inter_coupling = - 1 + sys->density_mv[sys->av[jv].iaa] * (1/sys->bulk_density_mv[sys->av[jv].iaa] - 1/sys->bulk_density_mv[sys->av[iv].iaa]);
            double enhance_factor_wvv = wvv[iv][jv][0]==0? 1 : sys->av[jv].multi / wvv[iv][jv][0];
            wvv[iv][jv][0] *= enhance_factor_wvv;
            nhkvv[iv][jv][0] = -wvv[iv][jv][0];
        }
    }
#endif

void perform_xvv_enhancement(IET_Param * sys, IET_arrays * arr){
  #ifdef _EXPERIMENTAL_
    perform_xvv_eenhancement(sys, arr);
  #else
    shift_xvv_by_1(sys, arr->wvv, arr->n_wvv, arr->nhkvv, arr->n_nhkvv);
    insert_xvv_at_0(sys, arr->wvv, arr->nhkvv);
  #endif
}
