const char * szHelpInteractive = "\
  Interactive mode command:\n\
    help, ha/-ha/-h/-help           show help\n\
    end, complete, finalize         end the calculation\n\
    q, exit, quit, ^D               directly exit program\n\
    r, run, continue, resume, g, go continue the calculation\n\
    ouuv, auuv, ohuv, o, ohi, ocuv  save files (see main help)\n\
    truv                            set the compression precision\n\
    delvv, delrism, delhi           change step factors\n\
    errtol, errtolrism, errtolhi    change error tolerance\n\
    stepmaxrism, stepmaxhi          change RISM/HI max steps\n\
    cceil, ccutoff                  change the HNC cutoff\n\
    closure-enhance                 change the closure enhancement factor\n\
    lap, mem, reports               reports of the time and memory cost till now\n\
    Guv, energy                     report the solvation energy till now\n\
    [other commands]                execute OS command via system(...)\n\
";
void subroutine_interactive(void * _id){
    IET_Param * sys = global_sys; IET_arrays * arr = global_arr;
    FILE * fin = stdin;
    FILE * fout = stderr;
    bool run_first_time = true;
    char *userName = getenv("LOGNAME");
    int iszfn_path = 0; for (iszfn_path=strlen(szfn_path)-1; iszfn_path>0; iszfn_path--) if (szfn_path[iszfn_path]=='/') break; if (szfn_path[iszfn_path]=='/') iszfn_path ++;
    char input[4096];
    while (true){
        bool interupt_main_loop = false;
        if (!sys->suspend_calculation){
            if (!fgets(input, sizeof(input), fin)) break;
            sys->suspend_calculation = true;
            interupt_main_loop = true;
        }
        while (!sys->is_suspend_calculation) usleep(100);
        lap_timer_none();

        if (run_first_time){ run_first_time = false; fprintf(fout, "%s", szHelpInteractive); }
        if (interupt_main_loop){
            recalculate_energy(sys, arr); //arr->display_solvation_energy_full(sys, fout, nullptr);
            print_default_IET_report(sys, arr);
        }

        time_t tt = time(nullptr); struct tm * tm_struct = localtime(&tt);
        fprintf(fout, "%s@%s:%s%s%s %02d:%02d:%02d> ", userName, software_name, sys->is_log_tty?prompt_path_prefix:"", &szfn_path[iszfn_path], sys->is_log_tty?prompt_path_suffix:"", tm_struct->tm_hour, tm_struct->tm_min, tm_struct->tm_sec);
        if (!fgets(input, sizeof(input), fin)){ fprintf(fout, "\n"); break; }

        StringNS::string sl[20]; int nw = analysis_line(input, sl, 20, true);
        if (sl[0]=="help"){
            fprintf(fout, "%s", szHelpInteractive);
        } else if (sl[0]=="-h" || sl[0]=="--h" || sl[0]=="-help" || sl[0]=="--help" || sl[0]=="-ha" || sl[0]=="--ha" || sl[0]=="ha"){
            char * help_search_str = nullptr;
            if (nw<=1 || sl[1]=="all"){
                fprintf(fout, "%s %s\n", software_name, software_version);
                fprintf(fout, "%s", szHelpMore);
                printf("%s%s%s%s%s", szHelp1, szHelpXTC, szHelpMP, szHelp2, szHelpLibZ);
                fprintf(fout, szHelpSecRISM3D_RISM, check_error_tol(0), check_error_tol(0));
                fprintf(fout, szHelpSecRISM3D_HI, check_error_tol(0), check_error_tol(0));
                fprintf(fout, "%s", szHelpSecATOM);
                fprintf(fout, "%s", szHelpSecRISM3DS);
                fprintf(fout, "%s%s", szHelpAdvanced, szHelpAdvancedOthers);
                fprintf(fout, "%s", szHelpNote);
            } else if (sl[1]=="sects"||sl[1]=="sections"){
                fprintf(fout, szHelpSecRISM3D_RISM, check_error_tol(0), check_error_tol(0));
                fprintf(fout, szHelpSecRISM3D_HI, check_error_tol(0), check_error_tol(0));
                fprintf(fout, "%s", szHelpSecRISM3DS);
            } else if (sl[1]=="sect"||sl[1]=="section"){
                if (nw>2){
                    if (sl[2]=="RISMHI" || sl[2]=="RISM") fprintf(fout, "%s", szHelpSecRISM3D_RISM);
                    if (sl[2]=="RISMHI" || sl[2]=="HI") fprintf(fout, "%s", szHelpSecRISM3D_HI);
                    if (sl[2]=="ATOM") fprintf(fout, "%s", szHelpSecATOM);
                }
            } else for (int iw=1; iw<nw; iw++){
                help_search_str = sl[iw].text;
                const char * helps[] = { szHelp1, szHelpXTC, szHelpMP, szHelp2, szHelpLibZ, szHelpMore, szHelpAdvanced, szHelpAdvancedOthers, szHelpSecRISM3D_RISM, szHelpSecRISM3D_HI, szHelpSecRISM3DS, szHelpSecATOM, szHelpNote };
                bool found = false; StringNS::string hkey = help_search_str;
                for (int ih=0; ih<sizeof(helps)/sizeof(helps[0]); ih++){
                    int im = 0; StringNS::string templ = helps[ih];
                    while (im<templ.length-hkey.length){
                        if (StringNS::string(&templ.text[im], hkey.length) == hkey){
                            found = true;
                            int begin = im; int end = im;
                            while (begin>0 && helps[ih][begin-1]!='\n') begin --;
                            while (end+1<templ.length && helps[ih][end+1]!='\r' && helps[ih][end+1]!='\n') end ++;
                            if (isatty(fileno(stdout))){
                                int out1=begin; int out2=begin;
                                while (out2<end-hkey.length){
                                    if (StringNS::string(&templ.text[out2], hkey.length) == hkey){
                                        fwrite(&helps[ih][out1], 1, out2-out1, stdout);
                                        printf("%s", prompt_highlight_prefix);
                                        fwrite(&helps[ih][out2], 1, hkey.length, stdout);
                                        printf("%s", prompt_highlight_suffix);
                                        out2 += hkey.length; out1 = out2;
                                    } else out2 ++;
                                }
                                if (out1<end) fwrite(&helps[ih][out1], 1, end-out1+1, stdout);
                                printf("\n");
                            } else {
                                fwrite(&helps[ih][begin], 1, end-begin+1, stdout);
                                printf("\n");
                            }
                            im = end+1;
                        }
                        im ++;
                    }
                }
                if (!found) printf("%s : no help entry for %s\n", software_name, help_search_str);
            }
        } else if (sl[0]=="end" || sl[0]=="complete" || sl[0]=="finalize" || sl[0]=="final" || sl[0]=="fine" || sl[0]=="fin"){
            sys->stepmax_rism = sys->stepmax_hi = 0; sys->suspend_calculation = false;
        } else if (sl[0]=="q" || sl[0]=="exit" || sl[0]=="quit"){
            break;
        } else if (sl[0] == "g" || sl[0] == "go" || sl[0] == "r" || sl[0]=="continue" || sl[0] == "resume" || sl[0]=="run"){
            sys->suspend_calculation = false;
        } else if (sl[0]=="delvv"){
            if (nw>1) sys->delrism = sys->delhi = atof(sl[1].text);
        } else if (sl[0]=="delrism"){
            if (nw>1) sys->delrism = atof(sl[1].text);
        } else if (sl[0]=="delhi"){
            if (nw>1) sys->delhi = atof(sl[1].text);

        } else if (sl[0]=="errtol"){
            if (nw>1) sys->errtolrism = sys->errtolhi = atof(sl[1].text);
        } else if (sl[0]=="errtolrism" || sl[0]=="errtol_rism"){
            if (nw>1) sys->errtolrism = atof(sl[1].text);
        } else if (sl[0]=="errtolhi" || sl[0]=="errtol_hi"){
            if (nw>1) sys->errtolhi = atof(sl[1].text);

        } else if (sl[0]=="stepmax" || sl[0]=="maxstep"){
            if (nw>1) sys->stepmax_rism = sys->stepmax_hi = atoi(sl[1].text);
        } else if (sl[0]=="stepmaxrism" || sl[0]=="stepmax_rism" || sl[0]=="maxsteprism" || sl[0]=="maxstep_rism"){
            if (nw>1) sys->stepmax_rism = atoi(sl[1].text);
        } else if (sl[0]=="stepmaxhi" || sl[0]=="stepmax_hi" || sl[0]=="maxstephi" || sl[0]=="maxstep_hi"){
            if (nw>1) sys->stepmax_hi = atoi(sl[1].text);
        } else if (sl[0]=="energy" || (sl[0]=="Guv" && sl[0].text[0]=='G')){
            recalculate_energy(sys, arr); //arr->display_solvation_energy_full(sys, fout, nullptr); arr->display_solvation_correlations(sys, fout, nullptr);
            print_default_IET_report(sys, arr);
        } else if (sl[0]=="cceil" || sl[0]=="gceil"){
            if (nw>1){
                double exp_cutoff = atof(sl[1].text);
                if (exp_cutoff>0) sys->ccutoff = log(exp_cutoff);
                else { fprintf(stderr, "%s : interactive : invalid HNC ceil %s\n", software_name, sl[1].text); }
            }
        } else if (sl[0]=="ccutoff" || sl[0]=="gcutoff"){
            if (nw>1){
                sys->ccutoff = log(atof(sl[1].text));
            }
        } else if (sl[0]=="closure-enhance" || sl[0]=="closure_enhance"){
            if (nw>1){
                sys->closure_enhance_level = atof(sl[1].text);
            }
        } else if (sl[0]=="lap-report" || sl[0]=="lap" || sl[0]=="reports"){
            lap_display_timers(fout);
        } else if (sl[0]=="memory-report" || sl[0]=="mem-report" || sl[0]=="memo-report" || sl[0]=="memory" || sl[0]=="memo" || sl[0]=="mem" || sl[0]=="reports"){
            fprintf(fout, "memory total assigned : %d blocks, ", _memory_blk_total);
              if (_memory_total>1024*1024*1024) fprintf(fout, "%.1f GB\n", _memory_total/(1024.0*1024*1024));
              else if (_memory_total>1024*1024) fprintf(fout, "%.1f MB\n", _memory_total/(1024.0*1024));
              else if (_memory_total>1024) fprintf(fout, "%.1f KB\n", _memory_total/(1024.0));
        } else {
            for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = ' ';
            unsigned char sys_ret = system(sl[0].text);
        }
        lap_timer_none();
        //if (nw>0) fprintf(fout, "%s@%s $ ", software_name, szfn_path);
    }
    #ifdef _LOCALPARALLEL_
        for (int i=0; i<sys->nt; i++) __mp_tasks[i] = MPTASK_TERMINATE;
    #endif
    exit(0);
}




void allocate_interactive(){
    pthread_t t;
    pthread_create(&t, nullptr, (void *(*)(void *))&subroutine_interactive, nullptr);
}
