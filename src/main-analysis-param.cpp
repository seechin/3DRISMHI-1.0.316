inline double atomised_molecule_uv_potential(IET_Param * sys, IET_arrays * arr, int iv, int ix, int iy, int iz);
bool analysis_command(IET_Param * sys, char * line, const char * line_orig, const char * script_name, int script_line, int first_char_offset);
int analysis_parameter_line(IET_Param * sys, const char * input, char * script_name, int script_line, const char * script_path=nullptr);
int analysis_parameter_line(IET_Param * sys, const char * in, char * script_name, int script_line, const char * script_path);
int analysis_params_file(IET_Param * sys, char * filename, int recursive_layer=0);
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int __iaanow = 0;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
const char * szHelp1 = "\
  Basic options:\n\
    -p, -solvent            iet parameter file, can be screen/con\n\
    $IETLIB (env string)    the folder containing iet parameter files\n\
    -s, -solute             solute file (solute file, prmtop file, or a folder)\n\
                            batch mode enabled when both -s and -f are folders\n\
";
#ifdef _GROMACS_
  const char * szHelpXTC = "\
    -f/traj/conf[ormation]  trajectory file: pdb, gro, xtc, or a folder\n\
     -b, -e, -be            xtc time range (ps), default all frames: -be 0,0\n\
     -dt                    xtc time step (ps), default all frames: -dt 0\n\
";
#else
  const char * szHelpXTC = "\
    -f, traj                trajectory file: pdb, gro, or a folder\n\
";
#endif
#if defined(_LOCALPARALLEL_) && defined(_LOCALPARALLEL_PTHREAD_)
 #if defined(_FFTWMPPARALLEL_)
 const char * szHelpMP = "\
   -nt/thread, -np/fork    thread or process number for FF/IET calculations\n\
   -j[obs], -ntb, -npb     thread or process number for batch mode\n\
   -ntf, -nt-fft           thread number for multithread fftw\n\
";
 #else
  const char * szHelpMP = "\
    -nt/thread, -np/fork   thread or process number for FF/IET calculations\n\
   -j[obs], -ntb, -npb     thread or process number for batch mode\n\
";
 #endif
#elif defined(_LOCALPARALLEL_)
const char * szHelpMP = "\
    -np/fork                process number for FF/IET calculations\n\
    -npb                    process number for batch mode\n\
";
#else
  const char * szHelpMP = "";
#endif
const char * szHelp2 = "\
    -nice[-level] 0         nice level of this process\n\
    -nr 100x100x100         grid number(s). Default: 100, or 100 100 100 \n\
    -rc 1                   interaction cutoff for LJ and Coulomb\n\
    -pwd, -cwd              directory for output, default: current folder\n\
    -log                    log file, could be stdout/con/screen/stderr\n\
    -v[bose] [0/1/2]        verbose level: 0=-silent, 1=-brief, 2=-detailed/-v\n\
    -debug[-level] 0/1/2/3  debug level: 0=-nodebug, 1=-debug, 2=-debugging\n\
      -[no-]debug-crc       calculate & display CRC32 of important data blocks\n\
      -[no-]debug-xvv       display some elements of wvv & nhkvv\n\
";
#ifdef _LIBZ_
  const char * szHelpLibZ = "\
    -o[ut][0/1/2]           output file name. Compress: 0:none, 1:fast, 2:best\n\
     -ov[erwrite][0/1/2]    output overwrites existing (see -o for compression)\n\
    -a[ppend][0/1/2]        output appends existing (see -o for compression)\n\
";
#else
  const char * szHelpLibZ = "\
    -o[ut]                  output file name (in TS4S)\n\
     -ov[eride]             output overwrites existing\n\
    -a[ppend]               output appends existing\n\
";
#endif
const char * szHelp3 = "\
    -run                    run script, can specify anything except -p or -h\n\
    -cmd, -do               add a command to the queue\n\
    -do-rism-kh             perform 3DRISM-KH\n\
    -do-rismhi-d2/kh        perform 3DRISM-HI-D2MSA/KH\n\
    -[default-]box          default size of box, applied if missing in trajs\n\
    -set-box                set box size to all frames\n\
";

const char * szHelpCommands = "\
  Commands:\n\
    clear/reset             clear memory for IET calculations\n\
    end/stop                end command queue for current frame\n\
    done/break              end command queue for current block\n\
    build-force-field       force to (re-)build ff\n\
    rism/hi[:step=100]      run rism/hi with given steps\n\
      rism: rism/ssoz, hi: nohi/hi/hshi\n\
    closure[-a/m]:...       20 (or less) closures, or one closure for all\n\
      closures: HNC/HNCB, PY/PH, MSA/KGK, KH/PSE2/.../PSE10\n\
      other closures: D2/DH/MS, MHNC/BPGGHNC/VM/MP\n\
    cf[-a/m]=1,1,...        20 (or less) closure factors\n\
      closure/cf surfix: -a atomic, -m molecular, or specify one for all\n\
      closure = closure-m\n\
    ti[_begin],step=10      do TI for following commands till \"done/break\"\n\
    hold:none,1,2,3,...     the sites that being held in rism/ssoz\n\
    density,dielect         20 (or less) items, each for a molecule\n\
    temperaute              set the temperature of the whole system\n\
    scale,param=val,...     scale values of: ff/uuv, lj, coul\n\
    display/calc[ulate]     calculate or display quantities\n\
      display/calc: energy, Euv/Ef/LJ/Coul, Cuv, dN, TS, SFE/Guv, rdf, all\n\
    report                  generate summary of Euv/Ef, Cuv or all\n\
    save, save-at-the-end, save-at-each-frame\n\
      save: ff, lj, coul; iet, cuv, huv, dd; all; guv, rdf\n\
      pick/save-sites: the sites picked for output, default: all\n\
    rdf-content             the same as -rdf-content but in -cmd/-do\n\
    Command running time specifier:\n\
      @b[egin], @e[nd]      run only at beginning/end (before/after all frames)\n\
      @number (e.g. @1, @5) insert command to specified location, begin with 1\n\
";
const char * szHelpAdvanced = "\
  Shortcuts:\n\
    -Lorentz-Berthelot      = -arith-sigma\n\
    -Yukawa/YukawaFFT       = -Coulomb YukawaFFT\n\
  Advanced settings: (all optional)\n\
    [#]include filename    (in -p file) include another file\n\
    -%14.7g                 output data format\n\
    -temperature, -T 298    temperature in Kelvin\n\
    -ff opls[aa]/amber/gaff\n\
      -amberff/-ffamber     = -Tdef 502.97 -arith-sigma\n\
      -oplsff/-ffopls[aa]   = -Tdef 120.27 -geo-sigma\n\
      -gaff/-ffgaff         = -Tdef 120.27 -arith-sigma\n\
    -coul[omb], -es         Electricstatic field algorithm, can be:\n\
                            Coulomb, dielect/dm, Yukawa/YukawaFFT\n\
    -dielect 1              dielect const of each solvent molecule\n\
    -rdf-grps               pairs of RDFs: solute_atom-solvent_atom[=grp],...\n\
    -rdf-content rdf        (-rdf-for) rdf/h/dd/c/ch, lj/coul/ef/uuv/ff\n\
    -rdf-bins 50            number of rdf bins\n\
    -arith[metic]-sigma     sigma combining rule: arithmetic rule\n\
    -geo[metric]-sigma      sigma combining rule: geometric rule\n\
    -list, -ls              list settings and parameters (no calculation)\n\
    -test                   test settings but do not actually run calculations\n\
    -rb 0.052911            (-Bohr-radius) minimal hardsphere radius\n\
    -rq 0                   (-rcharge) partial charge bead size (0: sigma/2)\n\
    -significant-digits     (-sd), digits for IETS output, can be float/double\n\
    -dynamic-delvv 1        (-closure-enhance) dynamic factor for convergence\n\
    -xvv-extend 0           number of times to extend wvv&nhkvv\n\
    -xvv-scheme so          scheme for xvv, can be so=amber/no/sc/sx/s1\n\
    -xvv-scale              factors to scale hvv, order: columns in [gvv_map]\n\
    -allow/-forbid          turn on/off options:\n\
        pme                 allow/forbid to perform PME, default on\n\
        [Coulomb-]renorm[alization] allow/forbid Coulomb renormalization in RISM\n\
        missing-xvv         allow/forbid to skip undefeined xvv, default on\n\
        missing-zeta        allow/forbid to skip undefeined zeta, default on\n\
        exceed-ram          allow/forbid to exceed physical memory, default off\n\
        alloc-rmin/Ef       * allow/forbid to generate rmin/Ef\n\
    -page[-]size 4096       page size in compressed TS4S file\n\
    -zeta-line/term/item    define a zeta term, see [zeta] for details\n\
    -dielect-hi             dielect const for Coulomb based HI theories\n\
    -ccutoff 5              cutoff for hybrid closures\n\
    -cceil 148              (-gceil) = exp(ccutoff), cutoff for hybrid closures\n\
    -dielect-y[ukawa] 1     dielectric const for Yukawa\n\
    -gcutoff-liquid         Liquid density threshold for reported correlations\n\
    -gcutoff-ef             Liquid density threshold for electric field energy\n\
    -Tdef 120.27            (-default_temperature) T of forcefield energy unit\n\
    -external-electrofield  external electric field in eV/nm, default: 0 0 0\n\
";

#ifdef _INTERACTIVE_
  const char * szHelpInteractive = "\
    -[no-]interactive       (-intervene) allow press enter to intervene\n\
";
#else
  const char * szHelpInteractive = "";
#endif

const char * szHelpInternal = "\
  Internal features: (all optional, cautious!)\n\
    -dielect-from-dipole    automatically calculate dielects from dipoles\n\
      -dielect-use-original disable any recalculation of dielect constants\n\
    -dipole-from-dielect    automatically calculate dipoles from dielects\n\
      -dipole-use-original  disable any recalculation of dipoles\n\
    -[in]homo-iet           homogeneous rism, default: off\n\
    -pme-gamma              * long range PME decreasor, default value: 2/rcoul\n\
    -gvv-align[ment]-left   the input gvv is aligned left (default)\n\
    -gvv-align[ment]-right  the input gvv is aligned right\n\
";

const char * szHelpSecRISM3D_RISM = "\
  [iet] or [solvent] (also [iet-hi]) defines RISM parameters in -p/-i:\n\
    -gvv, -hvv              solvent-solvent RDF (in r space) file and dr\n\
    -ndiisrism              (-ndiis) DIIS steps for RISM, default: 0 or 5\n\
    -errtolrism             (-errtol) RISM error tolerance, default: 1e-12\n\
    -density                densities for each solvent molecule\n\
    -bulk[-]density         bulk densities for each solvent molecule\n\
    -rlj, -rcoul            LJ/Coulomb cutoff, -rc for both, default: 1\n\
";
//-dielect-rism-on/off    automatically set (or not) dielect for RISM.\n
const char * szHelpSecRISM3D_HI = "\
  [hi] (also [iet-hi]) defines HI parameters in -p/-i:\n\
    -ndiishi                (-ndiis) DIIS steps for HI, default: 0 or 5\n\
    -errtolhi               (-errtol) HI error tolerance, default: 1e-12\n\
    -dipole                 dipoles for all solvent molecules\n\
    -zvv, -zeta             zeta file in r-space and dr, nvm^2 cols and -ff unit\n\
    -scale.zeta             scaling factors for each molecule\n\
    -ccutoff-[hs]hi         potential cutoff for theta function in HSHI\n\
    -lsa, -lsb              default: -lsa 0.3 -lsb auto\n\
                            Auto value of B = (lnλ - <ab|zeta|ab>/<c|c>)/Aρ\n\
    -theta                  theta cutoff of potential (in kT), default: -theta 5\n\
    -phi-cutoff             lowerbound of phi in HI equation, default: -20\n\
    -lambda                 a list of λ. llambda for lnλ, default: all 1 (lnλ=0)\n\
";
const char * szHelpSecRISM3DS = "\
  [solute] defines solute in -s:\n\
    atom_name mole_name mass charge sigma epsilon\n\
    repeat above atoms: \"#repeat 3/three atoms 106/one-hundred-and-six times\"\n\
";
const char * szHelpSecATOM = "\
  [atom] defines solvent molecules in -p/-i:\n\
    7 cols: atom_name mole_name index group_index charge sigma epsilon\n\
  [bond] or [pair] defines pairs for solvent atoms in -p/-i:\n\
    3 cols: atom_name_or_index_1 atom_name_or_index_2 bond_length\n\
    4 cols: atom_name_or_index_1 atom_name_or_index_2 bond_length bond_stdev\n\
  [gvv_map] defines the mapping of the solvent-solvent RDF (for gvv)\n\
    3 cols: atom_name_or_index_1 atom_name_or_index_2 gvv_col_num\n\
  [zeta] defines the information of zeta functions:\n\
    4 cols: molecule_1 molecule_2 zeta_cutoff_r zeta_at_origin [zeta_unit]\n\
";
const char * szHelpNote = "Note: \n\
  1. Default units: [L]=nm, [E]=kJ/mol, [q]=e\n\
  2. All the paremeters could be defined in -p or command line arguments.\n\
  3. Command line arguments will overwrite settings in -p.\n\
  4. Files included in -p will be searched from the path of -p file.\n\
  5. -s and -f are searched from the current working directory.\n\
  6. -o and -log are searched from the current or -pwd directory.\n\
  6. Use '--' to recognize filenames with leading '-'\n\
";
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int closure_from_string(StringNS::string closure_str, STKeywordTableUnit * alias, int count){
    for (int i=0; i<count; i++) if (alias[i].name){
        if (closure_str == alias[i].name) return alias[i].id;
    }
    return -1;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int analysis_solvent_atom_additional_word(IET_Param * sys, AtomSite * atom, StringNS::string sline, const char * script_name, int script_line){
    if (sline.text[0]=='#' || (sline.text[0]=='/'&&sline.text[1]=='/')) return -1; // stop further analysis
    if (sline=="done" || sline=="break") return -1;

    if (sline.sub(0, 7)=="charge:" || sline.sub(0, 7)=="charge="){
        atom->charge = atof(&sline.text[7]);
    } else if (sline.sub(0, 2)=="q:" || sline.sub(0, 2)=="q="){
        atom->charge = atof(&sline.text[2]);
    } else if (sline.sub(0, 4)=="esp:" || sline.sub(0, 4)=="esp="){
        atom->charge_esp = atof(&sline.text[4]);
    } else if (sline.sub(0, 11)=="charge-esp:" || sline.sub(0, 11)=="charge-esp=" || sline.sub(0, 11)=="charge_esp:" || sline.sub(0, 11)=="charge_esp="){
        atom->charge_esp = atof(&sline.text[11]);
    } else if (sline.sub(0, 6)=="q-esp:" || sline.sub(0, 6)=="q-esp=" || sline.sub(0, 6)=="q_esp:" || sline.sub(0, 6)=="q_esp="){
        atom->charge_esp = atof(&sline.text[6]);

    } else if (sline.sub(0, 3)=="dc:" || sline.sub(0, 3)=="dc=" || sline.sub(0, 18)=="dipole-correction:" || sline.sub(0, 18)=="dipole-correction=" || sline.sub(0, 18)=="dipole_correction:" || sline.sub(0, 18)=="dipole_correction="){
        StringNS::string this_sline = sline.sub(18);
            if (sline.sub(0, 3)=="dc:" || sline.sub(0, 3)=="dc=") this_sline = sline.sub(3);
        StringNS::string sl[4]; int nw = analysis_csv_line(sline.sub(3), sl, 4);
        if (nw>=2){
            atom->ren_charge = atof(sl[0].text);
            atom->ren_bond = atof(sl[1].text);
            atom->ren_dielect = 1;
            if (nw>=3) atom->ren_dielect = atof(sl[2].text);
            //if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: %s[%d] : solvent %s.%s defined with dipole correction: q=%g, bond=%g, dielect=%g\n", get_second_fn(script_name), script_line, atom->mole, atom->name, atom->ren_charge, atom->ren_bond, atom->ren_dielect);
        } else {
            fprintf(sys->log(), "%s%s : %s[%d] : error : not enough arguments for dc: %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(script_name), script_line, sline.text, sys->is_log_tty?color_string_end:"");
            return 0;
        }

//printf("atom %s:%s has dc %g %g %g\n", atom->mole, atom->name, atom->ren_charge, atom->ren_bond, atom->ren_dielect);

    } else {
        fprintf(sys->log(), "%s%s : %s[%d] : error : unrecognizable string %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(script_name), script_line, sline.text, sys->is_log_tty?color_string_end:"");
        return 0;
    }

    return 1;
}
bool analysis_solvent_atom_7(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line, int nargs){
    bool ret = true;
    //if (argc<7){ fprintf(sys->log(), "%s : %s[%d] : incomplete atom line\n", software_name, script_name, script_line); return false; }
    if (sys->nav >= MAX_SOL){ fprintf(sys->log(), "%s : %s[%d] : too many solvents\n", software_name, get_second_fn(script_name), script_line); return false; }
    int ifound = -1; int grp = atoi(argv[3]); for (int i=0; i<sys->nav && ifound<0; i++) if (sys->av[i].grp == grp) ifound = i;
    int id = atoi(argv[2]);
    ElementNS::Element * element = ElementNS::get_element_parameter(ElementNS::get_atom_element(argv[0]));
    if (!element) fprintf(sys->log(), "%s%s : %s[%d] : warning : no element found for atom %s%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, get_second_fn(script_name), script_line, argv[0], sys->is_log_tty?color_string_end:"");
    double mass = element? element->m : 0;
    if (ifound<0 || ifound>=sys->nav){
      // preparing
        bool is_this_atom_key = false;
        if (sys->nav<=0){
            is_this_atom_key = true;
        } else if (StringNS::string(argv[1]) != sys->av[sys->nav-1].mole){
            __iaanow ++; is_this_atom_key = true;
        } else is_this_atom_key = false;
      // fill atom infromation
        sys->av[sys->nav].init(id, grp, argv[1], __iaanow, argv[0], mass, atof(argv[4]), atof(argv[5]), atof(argv[6]));
      // fill additional information of atom
        for (int iarg=7; iarg<nargs; iarg++){
            int ret_additional = analysis_solvent_atom_additional_word(sys, &sys->av[sys->nav], argv[iarg], script_name, script_line);
            if (ret_additional==0) ret = false; else if (ret_additional<0) break;
        }
        /*if (nargs>7){
            sys->av[sys->nav].charge = atof(argv[7]);
            sys->av[sys->nav].dipole = atof(argv[8]);
        }*/
        sys->av[sys->nav].is_key = is_this_atom_key;
        sys->nav ++;
    } else {
        sys->av[ifound].multi ++;
    }
  // register to the atom list
    if (sys->n_atom_list+1 >= sys->nmax_atom_list || !sys->atom_list){
        sys->nmax_atom_list += 100;
        SYSITEM_AtomNameList * al = (SYSITEM_AtomNameList*) malloc(sizeof(SYSITEM_AtomNameList) * sys->nmax_atom_list);
        if (sys->atom_list){
            memcpy(al, sys->atom_list, sizeof(SYSITEM_AtomNameList)*sys->nmax_atom_list); free(sys->atom_list);
        }
        sys->atom_list = al;
        if (!sys->atom_list){ fprintf(sys->log(), "%s : malloc failure\n", software_name); exit(-1); }
    }
    if (sys->atom_list){
        memset(&sys->atom_list[sys->n_atom_list], 0, sizeof(SYSITEM_AtomNameList));
        StringNS::string st_argv0 = argv[0]; StringNS::string st_argv1 = argv[1];
        memcpy(sys->atom_list[sys->n_atom_list].name, argv[0], st_argv0.length>31?31:st_argv0.length);
        memcpy(sys->atom_list[sys->n_atom_list].mole, argv[1], st_argv1.length>31?31:st_argv1.length);
        sys->atom_list[sys->n_atom_list].index = id;
        sys->atom_list[sys->n_atom_list].grp = grp;
        sys->atom_list[sys->n_atom_list].iaa = __iaanow;
        sys->n_atom_list ++;
    }
  // end
    return ret;
}
bool analysis_solvent_atom(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line){
    int nargs = 0;
    for (; nargs<argc; nargs++) if (argv[nargs][0]=='#' || argv[nargs][0]==';') break;
    if (nargs<7){
        fprintf(sys->log(), "%s : %s[%d] : incomplete atom line\n", software_name, get_second_fn(script_name), script_line); return false;
    } else {
        return analysis_solvent_atom_7(sys, argv, argi, nargs, script_name, script_line, nargs);
    }
}
void analysis_format_string(IET_Param * sys, const char * fmt){
    int len = strlen(fmt);
    int i = 0;
    bool left_aligned = false; if (i<len && fmt[i]=='-'){ i++; left_aligned = true; }
    int fmt_len = 0; {
        int begin = i; while (i<len&&fmt[i]>='0'&&fmt[i]<='9') i++; int end = i;
        if (end>begin) fmt_len = atoi(&fmt[begin]);
    }
    int fmt_float = 0;
    if (i<len && fmt[i]=='.'){ i++;
        int begin = i; while (i<len&&fmt[i]>='0'&&fmt[i]<='9') i++; int end = i;
        if (end>begin) fmt_float = atoi(&fmt[begin]);
    }
    char fmt_indicator[2]; strncpy(fmt_indicator, "f", sizeof(fmt_indicator));
    if (i<len && (fmt[i]=='e'||fmt[i]=='E'||fmt[i]=='f'||fmt[i]=='F'||fmt[i]=='g'||fmt[i]=='G')){ fmt_indicator[0] = fmt[i]; i++; }

    if (fmt_len<7) fmt_len = 7;

    if (fmt_len<=0 && fmt_float<=0){
        snprintf(sys->output_realnumber_format, sizeof(sys->output_realnumber_format), " %%%s%s", left_aligned?"-":"", fmt_indicator);
    } else if (fmt_float<=0){
        snprintf(sys->output_realnumber_format, sizeof(sys->output_realnumber_format), " %%%s%d%s", left_aligned?"-":"", fmt_len, fmt_indicator);
    } else if (fmt_len<=0){
        snprintf(sys->output_realnumber_format, sizeof(sys->output_realnumber_format), " %%%s.%d%s", left_aligned?"-":"", fmt_float, fmt_indicator);
    } else {
        snprintf(sys->output_realnumber_format, sizeof(sys->output_realnumber_format), " %%%s%d.%d%s", left_aligned?"-":"", fmt_len, fmt_float, fmt_indicator);
    }
    if (fmt_len>0){
        sys->output_realnumber_length = fmt_len;
        snprintf(sys->output_title_format, sizeof(sys->output_realnumber_format), " %%%s%ds", left_aligned?"-":"", fmt_len);
    } else {
        snprintf(sys->output_title_format, sizeof(sys->output_realnumber_format), " %%%s%ds", left_aligned?"-":"", 7);
    }

    //printf(" align=%s length=%d float=%d indicator=[%s] fmt=[%s][%s]\n", left_aligned?"left":"right", fmt_len, fmt_float, fmt_indicator, sys->output_realnumber_format, sys->output_title_format);

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool analysis_solvent_bond(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line, const char * session_name){
    int nargs = 0;
    for (; nargs<argc; nargs++) if (argv[nargs][0]=='#' || argv[nargs][0]==';') break;
    if (nargs<3){ fprintf(sys->log(), "%s : %s[%d] : incomplete %s line\n", software_name, get_second_fn(script_name), script_line, session_name); return false; }
    if (nargs>=3){
        int atomi = -1; int atomj = -1; StringNS::string st_argv0 = argv[0]; StringNS::string st_argv1 = argv[1];
        if (StringNS::is_string_number(st_argv0)) atomi = atoi(argv[0]); else if (sys->atom_list) atomi = search_atom_list_grp(sys->atom_list, 0, sys->n_atom_list, argv[0], -1);
        if (StringNS::is_string_number(st_argv1)) atomj = atoi(argv[1]); else if (sys->atom_list) atomj = search_atom_list_grp(sys->atom_list, 0, sys->n_atom_list, argv[1], -1);
        double bond = atof(argv[2]);
        double bond_stdev = nargs>=4? atof(argv[3]) : 0;
        double bond_weight = nargs>=5? atof(argv[4]) : 1;
        if (bond<=0){
            fprintf(sys->log(), "%s%s : warning : ignoring bond[%s-%s]=%g%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, argv[0], argv[1], bond, sys->is_log_tty?color_string_end:"");
        } else {
            if (atomi>=0 && atomj>=0){
                if (sys->n_bond_list+1 <= sys->nmax_bond_list || ! sys->bond_list){
                    sys->nmax_bond_list += 100;
                    SYSITEM_BondList * bl = (SYSITEM_BondList*) malloc(sizeof(SYSITEM_BondList) * sys->nmax_bond_list);
                    if (sys->bond_list){
                        memcpy(bl, sys->bond_list, sizeof(SYSITEM_BondList)*sys->n_bond_list);
                        free(sys->bond_list);
                    }
                    sys->bond_list = bl;
                }
                if (!sys->bond_list){ fprintf(sys->log(), "%s : malloc failure\n", software_name); return false; }
                if (atomi>atomj){ int t = atomi; atomi = atomj; atomj = t; }
                sys->bond_list[sys->n_bond_list].grpi = atomi;
                sys->bond_list[sys->n_bond_list].grpj = atomj;
                sys->bond_list[sys->n_bond_list].bond = bond;
                sys->bond_list[sys->n_bond_list].bond_stdev = bond_stdev;
                sys->bond_list[sys->n_bond_list].weight = (bond_weight>1||bond_weight<=0)? 1 : bond_weight;
                sys->n_bond_list ++;
            } else {
                fprintf(sys->log(), "%s%s : error : undefined", sys->is_log_tty?color_string_of_error:"", software_name);
                if (atomi<0) fprintf(sys->log(), " %s %s", StringNS::is_string_number(st_argv0)?"atom group":"atom", argv[0]);
                if (atomi<0 && atomj<0) fprintf(sys->log(), " and");
                if (atomj<0) fprintf(sys->log(), " %s %s", StringNS::is_string_number(st_argv1)?"atom group":"atom", argv[1]);
                fprintf(sys->log(), "%s\n", sys->is_log_tty?color_string_end:"");
                return false;
            }
        }
        return true;
    } else return false;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool analysis_gvv_map(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line){
    int nargs = 0;
    for (; nargs<argc; nargs++) if (argv[nargs][0]=='#' || argv[nargs][0]==';') break;
    if (nargs<3){ fprintf(sys->log(), "%s : %s[%d] : incomplete gvv-map line\n", software_name, get_second_fn(script_name), script_line); return false; }
    if (nargs>=3){
        int atomi = -1; int atomj = -1; StringNS::string st_argv0 = argv[0]; StringNS::string st_argv1 = argv[1];
        if (StringNS::is_string_number(st_argv0)) atomi = atoi(argv[0]); else if (sys->atom_list) atomi = search_atom_list_grp(sys->atom_list, 0, sys->n_atom_list, argv[0], -1);
        if (StringNS::is_string_number(st_argv1)) atomj = atoi(argv[1]); else if (sys->atom_list) atomj = search_atom_list_grp(sys->atom_list, 0, sys->n_atom_list, argv[1], -1);
        int icol = atoi(argv[2]);
        if (icol<0){
            fprintf(sys->log(), "%s%s : warning : ignoring gvv-map[%s-%s]=%d%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, argv[0], argv[1], icol, sys->is_log_tty?color_string_end:"");
        } else {
            if (atomi>=0 && atomj>=0){
                if (sys->n_gvv_map+1 <= sys->nmax_gvv_map || ! sys->gvv_map){
                    sys->nmax_gvv_map += MAX_SOL;
                    SYSITEM_PairMapping * bl = (SYSITEM_PairMapping*) malloc(sizeof(SYSITEM_PairMapping) * sys->nmax_gvv_map);
                    if (sys->gvv_map){
                        memcpy(bl, sys->gvv_map, sizeof(SYSITEM_PairMapping)*sys->n_gvv_map);
                        free(sys->gvv_map);
                    }
                    sys->gvv_map = bl;
                }
                if (!sys->gvv_map){ fprintf(sys->log(), "%s : malloc failure\n", software_name); return false; }
                if (atomi>atomj){ int t = atomi; atomi = atomj; atomj = t; }
                sys->gvv_map[sys->n_gvv_map].grpi = atomi;
                sys->gvv_map[sys->n_gvv_map].grpj = atomj;
                sys->gvv_map[sys->n_gvv_map].col  = icol;
                sys->n_gvv_map ++;
            } else {
                fprintf(sys->log(), "%s%s : error : undefined", sys->is_log_tty?color_string_of_error:"", software_name);
                if (atomi<0) fprintf(sys->log(), " %s %s", StringNS::is_string_number(st_argv0)?"atom group":"atom", argv[0]);
                if (atomi<0 && atomj<0) fprintf(sys->log(), " and");
                if (atomj<0) fprintf(sys->log(), " %s %s", StringNS::is_string_number(st_argv1)?"atom group":"atom", argv[1]);
                fprintf(sys->log(), "%s\n", sys->is_log_tty?color_string_end:"");
                return false;
            }
        }
        return true;
    } else return false;
}

bool analysis_zeta_line(IET_Param * sys, char * argv[], int * argi, int argc, const char * script_name, int script_line){
    int nargs = 0;
    for (; nargs<argc; nargs++) if (argv[nargs][0]=='#' || argv[nargs][0]==';') break;
    if (nargs<=0) return true;
    if (nargs<4){ fprintf(sys->log(), "%s%s : %s[%d] : incomplete zeta line%s\n", sys->is_log_tty?color_string_of_synerr:"", software_name, script_name?get_second_fn(script_name):"args", script_line, sys->is_log_tty?color_string_end:""); return false; }

    int iaai = -1; int iaaj = -1; StringNS::string st_argv0 = argv[0]; StringNS::string st_argv1 = argv[1];
    if (StringNS::is_string_number(st_argv0)) iaai = atoi(argv[0]); else if (sys->atom_list) iaai = search_mole_list(sys->atom_list, 0, sys->n_atom_list, argv[0], -1);
    if (StringNS::is_string_number(st_argv1)) iaaj = atoi(argv[1]); else if (sys->atom_list) iaaj = search_mole_list(sys->atom_list, 0, sys->n_atom_list, argv[1], -1);
    double rc_zeta0 = atof(argv[2]);
    double deltaG = atof(argv[3]); double deltaG_defining_temperature = sys->default_temperature;
    if (argc>4 && !StringNS::is_string_number(argv[4])){ StringNS::string energy_unit_here = argv[4];
        if (energy_unit_here=="kJ/mol" || energy_unit_here=="(kJ/mol)" || energy_unit_here=="kJ_mol" || energy_unit_here=="(kJ_mol)"){
            deltaG_defining_temperature = 120.27;
        } else if (energy_unit_here=="kcal/mol" || energy_unit_here=="(kcal/mol)" || energy_unit_here=="kcal_mol" || energy_unit_here=="(kcal_mol)"){
            deltaG_defining_temperature = 502.97;
        } else if (energy_unit_here=="KT" || energy_unit_here=="(KT)"){
            deltaG_defining_temperature = sys->temperature;
        } else {
            fprintf(sys->log(), "%s%s : error : unrecognizable energy unit %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, argv[4], sys->is_log_tty?color_string_end:"");
            return false;
        }
    }

    int i_match = -1;
    if (sys->n_zeta_list>0 && sys->n_zeta_list){
        for (int i=0; i<sys->nmax_zeta_list; i++) if (sys->zeta_list[i].iaai==iaai && sys->zeta_list[i].iaaj==iaaj) i_match = i;
    }

    if (i_match>=0 && i_match<sys->n_zeta_list){
        sys->zeta_list[i_match].iaai = iaai;
        sys->zeta_list[i_match].iaaj = iaaj;
        sys->zeta_list[i_match].deltaG_in_Kelvin = deltaG * deltaG_defining_temperature;
        sys->zeta_list[i_match].rc_zeta0 = rc_zeta0;
        //printf("zeta[%d] = %d %d %g %g (%g)\n", i_match+1, iaai, iaaj, rc_zeta0, deltaG, deltaG_defining_temperature);
    } else {
        if (iaai>=0 && iaaj>=0 && iaai<sys->nmv && iaaj<sys->nmv){
            if (sys->n_zeta_list+1 <= sys->nmax_zeta_list || ! sys->zeta_list){
                sys->nmax_zeta_list += MAX_SOL;
                SYSITEM_ZetaList * zl = (SYSITEM_ZetaList*) malloc(sizeof(SYSITEM_ZetaList) * sys->nmax_zeta_list);
                if (sys->zeta_list){
                    memcpy(zl, sys->zeta_list, sizeof(SYSITEM_ZetaList)*sys->n_zeta_list);
                    free(sys->zeta_list);
                }
                sys->zeta_list = zl;
            }
            if (!sys->zeta_list){ fprintf(sys->log(), "%s : malloc failure\n", software_name); return false; }
            sys->zeta_list[sys->n_zeta_list].iaai = iaai;
            sys->zeta_list[sys->n_zeta_list].iaaj = iaaj;
            sys->zeta_list[sys->n_zeta_list].deltaG_in_Kelvin = deltaG * deltaG_defining_temperature;
            sys->zeta_list[sys->n_zeta_list].rc_zeta0 = rc_zeta0;
            sys->n_zeta_list ++;
        } else {
            fprintf(sys->log(), "%s%s : error : undefined molecule %s and %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, argv[0], argv[1], sys->is_log_tty?color_string_end:"");
            return false;
        }
        //printf("zeta[%d] = %d %d %g %g (%g)\n", sys->n_zeta_list, iaai, iaaj, rc_zeta0, deltaG, deltaG_defining_temperature);
    }
    return true;
}

bool set_sys_options(IET_Param * sys, const char * sline, bool flag, const char * script_name, int script_line){
    bool ret = true; StringNS::string sl[MAX_WORD]; int nw = analysis_csv_line(sline, sl, MAX_WORD);
    for (int i=0; i<nw; i++){
        if (sl[i]=="pme"){
            sys->perform_pme = flag;                //printf("PME: %s\n", flag?"true":"false");
        } else if (sl[i]=="renorm" || sl[i]=="renormalization" || sl[i]=="Coulomb-renorm" || sl[i]=="Coulomb_renorm" || sl[i]=="Coulomb-renormalization" || sl[i]=="Coulomb_renormalization"){
            sys->rism_coulomb_renormalization = flag;
        } else if (sl[i]=="missing-xvv" || sl[i]=="missing_xvv" || sl[i]=="xvv-missing" || sl[i]=="xvv_missing" || sl[i]=="skip-missing-xvv" || sl[i]=="skip_missing_xvv" || sl[i]=="skip-xvv-missing" || sl[i]=="skip_xvv_missing"){
            sys->check_zero_xvv = flag;             //printf("check_zero_xvv: %s\n", flag?"true":"false");
        } else if (sl[i]=="exceed-memory-capacity"||sl[i]=="exceed_memory_capacity"||sl[i]=="exceed-ram"||sl[i]=="exceed_ram"){
            sys->ignore_memory_capacity = flag;
            _ignore_memory_capacity = flag;         //printf("exceed_memory_capacity: %s\n", flag?"true":"false");
        } else if (sl[i]=="missing-zeta" || sl[i]=="missing_zeta" || sl[i]=="zeta-missing" || sl[i]=="zeta_missing" || sl[i]=="skip-missing-zeta" || sl[i]=="skip_missing_zeta" || sl[i]=="skip-zeta-missing" || sl[i]=="skip_zeta_missing"){
            sys->zeta_list_allow_missing = flag;    //printf("zeta_allow_missing: %s\n", flag?"true":"false");
        } else if (sl[i]=="rmin" || sl[i]=="alloc-rmin" || sl[i]=="alloc_rmin"){
            sys->b_allow_r2uvmin = flag;
        } else if (sl[i]=="Ef" || sl[i]=="alloc-Ef" || sl[i]=="alloc_Ef"){
            sys->b_allow_Ecoul0 = flag;
      #ifdef _EXPERIMENTAL_
        } else if (!set_sys_options_experimental(sys, sl, i, flag)){
      #else
        } else {
      #endif
            char buffer[1024]; strncpy(buffer, sl[i].text, sizeof(buffer));
            fprintf(sys->log(), "%s : %s[%d] : unknown flag \"%s\"\n", software_name, get_second_fn(script_name), script_line, buffer); ret = false;
        }
    }
    return ret;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void pre_analysis_params(int argc, char * argv[], int * debug_level, int * detail_level, char * log_filename, bool * help_out){
    *help_out = false; if (argc<=1) *help_out = true;
    for (int i=0; i<argc; i++){
        if (argv[i][0]=='-'){ StringNS::string key = argv[i];
            if (key == "-h" || key == "--h" || key == "-help" || key == "--help" || key == "-ha" || key == "--ha" || key=="-version" || key=="--version"){
                *help_out = true;
            } else if (key=="-debug-level" || key=="debug-level"  || key=="--debug-level" || key=="-debug_level" || key=="debug_level"  || key=="--debug_level"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; *debug_level = atoi(argv[i]); }
            } else if (key=="-debugging" || key=="debugging"  || key=="--debugging" || key=="-debuging" || key=="debuging"  || key=="--debuging"){
                *debug_level = 2;
            } else if (key=="-debug" || key=="debug"  || key=="--debug"){
                *debug_level = 1; if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; *debug_level = atoi(argv[i]); }
            } else if (key=="-v" || key=="--v" || key=="-verbose" || key=="verbose" || key=="--verbose"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; *detail_level = atoi(argv[i]); }
                else *detail_level = 2;
            } else if (key=="-detail-level" || key=="detail-level" || key=="--detail-level" || key=="-detail_level" || key=="detail_level" || key=="--detail_level"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; *detail_level = atoi(argv[i]); }
            } else if (key=="-detail" || key=="detail" || key=="--detail"){
                *detail_level = 2; if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; *detail_level = atoi(argv[i]); }
            } else if (key=="-brief" || key=="brief"  || key=="--brief"){
                *detail_level = 1;
            } else if (key=="-silent" || key=="silent" || key=="--silent"){
                *detail_level = 0;
            } else if (key=="-detailed" || key=="detailed" || key=="--detailed" || key=="-details" || key=="details" || key=="--details"){
                *detail_level = 2;
            } else if (key == "-log" || key == "log" || key == "--log"){
                if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(log_filename, argv[i]); }
            } else if (key == "-pwd" || key == "pwd" || key == "--pwd" || key == "-cwd" || key == "cwd" || key == "--cwd"){
                if (i+1<argc){ i++;
                    if (!cp_file_with_path(szfn_path, argv[i], nullptr)) szfn_path[0] = 0;   //printf("szfn_path = [%s]\n", szfn_path);
                }
            }
        }
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int analysis_parameter_line_pre(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line){
    int ret = 0; int i = *argi; bool analysis_script = !script_name? false : (!script_name[0]? false : true);
    if (!script_name || !script_name[0]) script_name = (char*)"args\0\0";
    StringNS::string key = argv[i];
    if (!analysis_script && (key=="-version" || key=="--version")){
        ret = 2;
    } if (!analysis_script && (key == "-h" || key == "--h" || key == "-help" || key == "--help" || key == "-ha" || key == "--ha" || key == "-ha" || key == "---help")){
        ret = 2 + 4;    // version string + brief help
        //if (key != "-h" && key != "--h") ret |= 8;
        if (key=="-ha"||key=="--ha"||key=="-help"||key=="--help") ret |= 8 + 16; // detailed help
        if (i+1<argc && argv[i+1][0] != '-'){
            i++; help_search_str = argv[i]; ret = 2 + 512; // version string + search for help strings
        }
    } else if (!analysis_script && (key == "-p" || key == "--p" || key == "-param" || key == "param" || key == "--param" || key=="-solvent" || key=="--solvent" || key=="solvent")){
        if (i+1<argc && argv[i+1][0]!='#' && argv[i+1][0]!=';' && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){
            if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++;
            i++; bool istty = sys->is_log_tty;
            bool cp_success = cp_file_with_path(info_file_name, argv[i], nullptr);
            if (!cp_success && sys->library_path){
                char sz_fn_tmp[MAX_PATH]; snprintf(sz_fn_tmp, sizeof(sz_fn_tmp), "%s%s%s", sys->library_path, sys->library_path[strlen(sys->library_path)-1]=='/'?"":"/", argv[i]);
                cp_success = cp_file_with_path(info_file_name, sz_fn_tmp, nullptr);
            }
            if (!cp_success){
                fprintf(sys->log(), "%s : %s[%d] : cannot open %s %s%s%s\n", software_name, get_second_fn(script_name), script_line, key.text, istty?prompt_path_prefix:"", argv[i], istty?prompt_path_suffix:""); ret = 1;
            } else if (!is_a_text_file(info_file_name)){
                fprintf(sys->log(), "%s : %s[%d] : not a text file: %s %s%s%s\n", software_name, get_second_fn(script_name), script_line, key.text, istty?prompt_path_prefix:"", argv[i], istty?prompt_path_suffix:""); ret = 1;
            }
            //if (StringNS::string(argv[i])=="help") ret = 8 + 16; else info_file_name = argv[i];
        }
    } else i++;
    *argi = i - 1;
    return ret;
}
int analysis_parameter_line(IET_Param * sys, char * argv[], int * argi, int argc, char * script_name, int script_line, const char * script_path=nullptr){
    int ret = 0; int i = *argi; bool analysis_script = (!script_name || !script_name[0])? false : (!script_name[0]? false : true);
    if (!analysis_script) script_name = (char*)"args\0\0";
    StringNS::string key = argv[i];
    bool istty = sys->is_log_tty;
    if (key == "--"){
   // ---------------------------------------------------------------------------------
   // The basic parameters for software -----------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if ((key == "-h" || key == "--h" || key == "-help" || key == "--help" || key == "help" || key == "-ha" || key == "--ha" || key == "ha")){
        if (analysis_script){
            fprintf(sys->log(), "%s : %s[%d] : \"%s\" ignored. Only valid as command arguments.\n", software_name, script_name, script_line, argv[i]);
        } else {
            ret = 2; if (key != "-h" && key != "--h") ret |= 4 + 8 + 16 + 32 + 64; if (key=="-help"||key=="--help"||key=="-ha"||key=="--ha") ret |= 1024+2048;
            if (i+1<argc && argv[i+1][0] != '-'){ i++; StringNS::string key2 = argv[i]; ret = 2;
                if (key2 == "all") ret = 2 + 4 + 8 + 16 + 32 + 64;
                else if (key2 == "none") ret = 2;
                else if (key2 == "advanced") ret = 2 + 4;
                else if (key2 == "iet" || key2 == "rism") ret = 8 + 1024;
                else if (key2 == "hi") ret = 8 + 2048;
                else if (key2 == "iethi" || key2 == "rismhi") ret = 8 + 1024 + 2048;
                else if (key2 == "atom") ret = 16;
                //else { fprintf(sys->log(), "%s : warning : no entry of %s\n%s\n", software_name, argv[i], szHelpMore); ret = 1; }
                else { help_search_str = argv[i]; ret = 512; }
            }
        }
    } else if (!analysis_script && key.text[0] != '-'){
        fprintf(sys->log(), "%s%s : %s[%d] : unknown string \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, script_name, script_line, argv[i], istty?color_string_end:""); ret = 1;
  #ifdef _LOCALPARALLEL_
   #ifdef _LOCALPARALLEL_PTHREAD_
    } else if (key == "-nt" || key == "--nt" || key == "nt" || key == "-thread" || key == "--thread" || key == "thread"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->nt = atoi(argv[++i]); if (sys->nt<1) sys->nt = 1; else if (sys->nt>MAX_THREADS) sys->nt = MAX_THREADS;
        #ifdef _LOCALPARALLEL_PTHREAD_
            sys->mp_by_fork = false;
        #else
            sys->mp_by_fork = true;
        #endif
    } else if (key=="-j" || key=="--j" || key=="-jobs" || key=="--jobs" ||key == "-ntb" || key == "--ntb" || key == "ntb"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ntb = atoi(argv[++i]); if (sys->ntb<1) sys->ntb = 1; else if (sys->ntb>MAX_THREADS) sys->ntb = MAX_THREADS;
        #ifdef _LOCALPARALLEL_PTHREAD_
            sys->mp_by_fork = false;
        #else
            sys->mp_by_fork = true;
        #endif
   #endif
   #ifdef _FFTWMPPARALLEL_
    } else if (key == "-ntf" || key == "--ntf" || key == "ntf" || key == "-nt-fft" || key == "--nt-fft" || key == "nt-fft" || key == "-nt_fft" || key == "--nt_fft" || key == "nt_fft"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ntf = atoi(argv[++i]); if (sys->ntf<1) sys->ntf = 1; else if (sys->ntf>MAX_THREADS) sys->ntf = MAX_THREADS;
   #endif
    } else if (key == "-np" || key == "--np" || key == "np" || key == "-fork" || key == "--fork" || key == "fork"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->nt = atoi(argv[++i]); if (sys->nt<1) sys->nt = 1; else if (sys->nt>MAX_THREADS) sys->nt = MAX_THREADS;
        sys->mp_by_fork = true;
    } else if (key == "-npb" || key == "--npb" || key == "npb"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ntb = atoi(argv[++i]); if (sys->ntb<1) sys->ntb = 1; else if (sys->ntb>MAX_THREADS) sys->ntb = MAX_THREADS;
        sys->mp_by_fork = true;
  #endif
    } else if (key == "-nice" || key == "--nice" || key == "nice" || key == "-nice-level" || key == "--nice-level" || key == "nice-level" || key == "-nice_level" || key == "--nice_level" || key == "nice_level"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->nice_level = atoi(argv[++i]);
    } else if (key == "-p" || key == "--p" || key == "-param" || key == "param" || key == "--param" || key == "-solvent" || key == "solvent" || key == "--solvent"){
        if (analysis_script){
            fprintf(sys->log(), "%s : %s[%d] : cannot specify \"%s\". Only valid as command arguments.\n", software_name, script_name, script_line, argv[i]); ret = 1;
        } else {
            if (i+1<argc && argv[i+1][0]!='#' && argv[i+1][0]!=';' && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){
                if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++;
                i++;
                bool cp_success = cp_file_with_path(info_file_name, argv[i], nullptr);
                if (!cp_success && sys->library_path){
                    char sz_fn_tmp[MAX_PATH]; snprintf(sz_fn_tmp, sizeof(sz_fn_tmp), "%s%s%s", sys->library_path, sys->library_path[strlen(sys->library_path)-1]=='/'?"":"/", argv[i]);
                    cp_success = cp_file_with_path(info_file_name, sz_fn_tmp, nullptr);
                }
                if (!cp_success){
                    fprintf(sys->log(), "%s : %s[%d] : cannot open %s %s%s%s\n", software_name, get_second_fn(script_name), script_line, key.text, istty?prompt_path_prefix:"", argv[i], istty?prompt_path_suffix:""); ret = 1;
                } else if (!is_a_text_file(info_file_name)){
                    fprintf(sys->log(), "%s : %s[%d] : not a text file: %s %s%s%s\n", software_name, get_second_fn(script_name), script_line, key.text, istty?prompt_path_prefix:"", argv[i], istty?prompt_path_suffix:""); ret = 1;
                }
            }
        }
    } else if (StringNS::string(key.text, 2)=="-%" || StringNS::string(key.text, 2)=="--%"){   // format string
        analysis_format_string(sys, StringNS::string(key.text, 2)=="--%"? &key.text[3] : &key.text[2]);
   // ---------------------------------------------------------------------------------
   // The input/output control parameters ---------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key == "-s" || key == "--s" || key == "-solute" || key == "--solute" || key == "solute"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; if (!cp_file_with_path(szfn_solute, argv[i], script_path)){ fprintf(sys->log(), "%s : %s[%d] : cannot open -s %s%s%s\n", software_name, get_second_fn(script_name), script_line, istty?prompt_path_prefix:"", argv[i], istty?prompt_path_suffix:""); ret = 1; }; }
    } else if (key == "-b" || key == "--b"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->time_begin = atof(argv[++i]);
    } else if (key == "-e" || key == "--e"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->time_end = atof(argv[++i]);
    } else if (key == "-dt" || key == "--dt"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->time_step = atof(argv[++i]);
    } else if (key == "-be" || key == "--be"){
        if (i+1<argc && argv[i+1][0] != '-'){ i++;
            StringNS::string slp[2]; char sep[4] = { '~', ',', 0, 0 };
            int nw = StringNS::analysis_general_line(sep, argv[i], slp, 2);
            if (nw>=2){
                char slpt[2][64]; memset(slpt, 0, sizeof(slpt));
                for (int i=0; i<2; i++) memcpy(slpt[i], slp[i].text, slp[i].length>sizeof(slpt[i])-1?sizeof(slpt[i])-1:slp[i].length);
                sys->time_begin = atof(slpt[0]); sys->time_end = atof(slpt[1]);
            } else {
                fprintf(sys->log(), "%s%s : %s[%d] : incorrect -be: %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
            }
        }
    } else if (key == "-f" || key == "--f" || key == "traj" || key == "-traj" || key == "conf" || key == "-conf" || key == "conformation" || key == "-conformation" || key == "conformations" || key == "-conformations"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_xtc, argv[i]); }
    } else if (key == "-pwd" || key == "pwd"){
        if (i+1<argc){ i++;
            if (!szfn_path[0]){
                fprintf(sys->log(), "%s%s : %s[%d] : cannot find directory: -pwd %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
            }
        }
        /*if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){
            if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++;
            i++; int chdir_ret = chdir(argv[i]);
            char * p_szfn_path = getcwd(szfn_path, sizeof(szfn_path));
        }*/
    } else if (key == "-log" || key == "--log" || key == "log"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_log, argv[i]); }
    } else if (key == "-i" || key == "--i" || key == "-in" || key == "--in" || key == "in" ){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_in, argv[i]); }
    } else if (key == "-o" || key == "--o" || key == "-out" || key == "--out" || key == "out"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 0; sys->output_compress_level = 1;
    } else if (key == "-o0" || key == "--o0" || key == "-out0" || key == "--out0" || key == "out0"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 0; sys->output_compress_level = 0;
    } else if (key == "-o1" || key == "--o1" || key == "-out1" || key == "--out1" || key == "out1"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 0; sys->output_compress_level = 1;
    } else if (key == "-o2" || key == "--o2" || key == "-out2" || key == "--out2" || key == "out2"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 0; sys->output_compress_level = 2;
    } else if (key == "-a" || key == "--a" || key == "-append" || key == "--append" || key == "append" ){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 1; sys->output_compress_level = 1;
    } else if (key == "-a0" || key == "--a0" || key == "-append0" || key == "--append0" || key == "append0" ){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 1; sys->output_compress_level = 0;
    } else if (key == "-a1" || key == "--a1" || key == "-append1" || key == "--append1" || key == "append1" ){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 1; sys->output_compress_level = 1;
    } else if (key == "-a2" || key == "--a2" || key == "-append2" || key == "--append2" || key == "append2" ){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = 1; sys->output_compress_level = 2;
    } else if (key == "-ov" || key == "--ov" || key == "-overwrite" || key == "--overwrite" || key == "overwrite"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = -1; sys->output_compress_level = 1;
    } else if (key == "-ov0" || key == "--ov0" || key == "-overwrite0" || key == "--overwrite0" || key == "overwrite0"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = -1; sys->output_compress_level = 0;
    } else if (key == "-ov1" || key == "--ov1" || key == "-overwrite1" || key == "--overwrite1" || key == "overwrite1"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = -1; sys->output_compress_level = 1;
    } else if (key == "-ov2" || key == "--ov2" || key == "-overwrite2" || key == "--overwrite2" || key == "overwrite2"){
        if (i+1<argc && (argv[i+1][0]!='-'||(i+2<argc&&StringNS::string(argv[i+1])=="--"))){ if (i+2<argc&&StringNS::string(argv[i+1])=="--") i++; i++; strcpy(szfn_out, argv[i]); };
        sys->output_override = -1; sys->output_compress_level = 2;

    } else if (key == "-rdf-bins" || key == "--rdf-bins" || key == "rdf-bins" || key == "-rdf_bins" || key == "--rdf_bins" || key == "rdf_bins"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->out_rdf_bins = atoi(argv[i]); }
    } else if (key == "-rdf-grps" || key == "--rdf-grps" || key == "rdf-grps" || key == "-rdf_grps" || key == "--rdf_grps" || key == "rdf_grps" || key == "-rdf-grp" || key == "--rdf-grp" || key == "rdf-grp" || key == "-rdf_grp" || key == "--rdf_grp" || key == "rdf_grp"){
        if (i+1<argc && argv[i+1][0]!='-'){ i++;
            StringNS::string sl[MAX_WORD]; StringNS::string sls[4]; StringNS::string sld[4]; int grp = -1;
            int nsl = analysis_csv_line(argv[i], sl, MAX_WORD, false, true);
            char sep_sv[4] = { '-', '-', '=', '=' }; char sep_ma[4] = { '.', ':', 0, 0 };
            for (int ist=0; ist<nsl; ist++){
                int id[2] = { 0, 0 }; StringNS::string mol[2] = { "", "" }; StringNS::string atom[2] = { "", "" };
                sls[0] = sls[1] = "\0";
                int nsls = analysis_general_line(sep_sv, sl[ist], sls, 4, false, true);
                for (int isls=0; isls<2; isls++){
                    if (isls>=nsls || sls[isls]=="*" || sls[isls]=="*.*" || sls[isls]=="*:*" || sls[isls]==".*" || sls[isls]=="*." || sls[isls]=="." || sls[isls]=="all"){
                        id[isls] = -1; // not specified, matches anything
                        sys->b_allow_r2uvmin = true;
                    } else if (StringNS::is_string_number(sls[isls])){
                        id[isls] = sls[isls].ToInt();
                    } else {
                        id[isls] = -1; // negative index meaning skip it and match names
                        int nsld = analysis_general_line(sep_ma, sls[isls], sld, 4, false, true);
                        if (nsld<1){
                            mol[isls] = atom[isls] = "*";
                        } else if (nsld>1){
                            mol[isls] = sld[0]; atom[isls] = sld[1];
                        } else {
                            mol[isls] = "*"; atom[isls] = sld[0];
                        }
                    }
                }
                if (nsls>2 && StringNS::is_string_number(sls[2])) grp = atoi(sls[2].text); else grp = -1;
                if (sys->n_rdf_grps < MAX_RDF_GRPS){
                    sys->rdf_grps[sys->n_rdf_grps].init(id[0], id[1], mol[0], atom[0], mol[1], atom[1], grp);
                    sys->n_rdf_grps ++;
                }
            }
        }
        int max_grp = sys->rdf_grps[0].grp; for (int i=1; i<sys->n_rdf_grps; i++) if (sys->rdf_grps[i].grp>max_grp) max_grp = sys->rdf_grps[i].grp;
        for (int i=0; i<=sys->n_rdf_grps; i++) if (sys->rdf_grps[i].grp<0) sys->rdf_grps[i].grp = ++max_grp;
//for (int i=0; i<sys->n_rdf_grps; i++) printf("%d %d.%s.%s - %d.%s.%s = %d\n", i+1, sys->rdf_grps[i].is, sys->rdf_grps[i].ms, sys->rdf_grps[i].as, sys->rdf_grps[i].iv, sys->rdf_grps[i].mv, sys->rdf_grps[i].av, sys->rdf_grps[i].grp);
    } else if (key == "-gcutoff-liquid" || key == "gcutoff-liquid" || key == "--gcutoff-liquid" || key == "-gcutoff_liquid" || key == "gcutoff_liquid" || key == "--gcutoff_liquid" || key == "-gcutoff" || key == "gcutoff" || key == "--gcutoff"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->gcutoff_liquid_occupation = atof(argv[++i]);
        }
    } else if (key == "-gcutoff-ef" || key == "gcutoff-ef" || key == "--gcutoff-ef" || key == "-gcutoff_ef" || key == "gcutoff_ef" || key == "--gcutoff_ef"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->gcutoff_ef_occupation = atof(argv[++i]);
        }
    } else if (key == "-rdf-content" || key == "rdf-content" || key == "--rdf-content" || key == "-rdf_content" || key == "rdf_content" || key == "--rdf_content" || key == "-rdf-for" || key == "rdf-for" || key == "--rdf-for" || key == "-rdf_for" || key == "rdf_for" || key == "--rdf_for"){
        if (i+1<argc && argv[i+1][0]!='-'){ i++; StringNS::string rdf_content = argv[i];
            if (rdf_content == "rdf" || rdf_content == "density"){
                sys->rdf_content = IETCMD_v_rdf;
            } else if (rdf_content == "g" || rdf_content == "guv"){
                sys->rdf_content = IETCMD_v_guv;
            } else if (rdf_content == "h" || rdf_content == "huv"){
                sys->rdf_content = IETCMD_v_huv;
            } else if (rdf_content == "hlr"){
                sys->rdf_content = -IETCMD_v_huv;
            } else if (rdf_content == "dd"){
                sys->rdf_content = IETCMD_v_dd;
            } else if (rdf_content == "nphi"){
                sys->rdf_content = -IETCMD_v_dd;
            } else if (rdf_content == "c" || rdf_content == "cuv"){
                sys->rdf_content = IETCMD_v_cuv;
            } else if (rdf_content == "ch" || rdf_content == "chuv"){
                sys->rdf_content = -IETCMD_v_cuv;
            } else if (rdf_content == "lj"){
                sys->rdf_content = IETCMD_v_ulj;
            } else if (rdf_content == "coul" || rdf_content == "coulomb"){
                sys->rdf_content = IETCMD_v_ucoul;
            } else if (rdf_content == "ff" || rdf_content == "uuv"){
                sys->rdf_content = IETCMD_v_uuv;
            } else if (rdf_content == "ef"){
                sys->rdf_content = IETCMD_v_Ef;
                sys->b_allow_Ecoul0 = true;
            } else {
                fprintf(sys->log(), "%s%s : %s[%d] : unrecognizable rdf content %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
            }
        }
    } else if (key == "-sd" || key == "sd" || key == "--sd" || key == "-significant-digits" || key == "significant-digits" || key == "--significant-digits" || key == "-significant_digits" || key == "significant_digits" || key == "--significant_digits"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->output_significant_digits = atoi(argv[++i]);
        } else if (i+1<argc && StringNS::string(argv[i+1])=="float"){
            i++; sys->output_significant_digits = 7;
        } else if (i+1<argc && StringNS::string(argv[i+1])=="double"){
            i++; sys->output_significant_digits = 18;
        } else if (i+1<argc) { i++;
            fprintf(sys->log(), "%s%s : %s[%d] : unrecognizable option for output precision: %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
        }
    } else if (key == "-pagesize" || key == "pagesize" || key == "--pagesize" || key == "-page-size" || key == "page-size" || key == "--page-size" || key == "-page_size" || key == "page_size" || key == "--page_size"){
        size_t compress_page_size = sys->compress_page_size;
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) compress_page_size = atol(argv[++i]);
        if (compress_page_size<1024){
            fprintf(sys->log(), "%s%s : warning : size of compress page too small (%lu), use 1KB instead%s\n", istty?color_string_of_synwarn:"", software_name, compress_page_size, istty?color_string_end:"");
            compress_page_size = 1024;
        } else if (compress_page_size>1024*1024*500){
            fprintf(sys->log(), "%s%s : warning : size of compress page too large (%lu), use 500MB instead%s\n", istty?color_string_of_synwarn:"", software_name, compress_page_size, istty?color_string_end:"");
            compress_page_size = 1024*1024*500;
        }
        sys->compress_page_size = compress_page_size;
    } else if (key == "-gvv" || key == "--gvv" || key == "gvv" || key == "-hvv" || key == "--hvv" || key == "hvv"){
        if (i+1<argc){
            if (key == "-gvv" || key == "--gvv" || key == "gvv") sys->gvv_specification = 1; else sys->gvv_specification = -1;
            sys->drrism = 0; n_szfn_gvvs = 0;
            while (i+1<argc){
                if (StringNS::is_string_number(argv[i+1])){
                    i++; double new_drrism = atof(argv[i]);
                    if (sys->drrism<=0) sys->drrism = new_drrism;
                    if (fabs(sys->drrism - new_drrism)>MACHINE_REASONABLE_ERROR){
                        fprintf(sys->log(), "%s%s : %s[%d] : dr of gvv redefined : %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:"");
                        ret = 1;
                    }
                } else if (StringNS::string(argv[i+1]) == "--" || argv[i+1][0] != '-'){
                    i++; if (StringNS::string(argv[i]) == "--") i++;
                    if (i<argc){
                        if (n_szfn_gvvs+1>=MAX_GVV_FILES){
                        } else {
                            if (!cp_file_with_path(&szfn_gvvs[n_szfn_gvvs*MAX_PATH], argv[i], script_path)){
                                fprintf(sys->log(), "%s%s : %s[%d] : file not found : %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:"");
                                ret = 1;
                            }
                            n_szfn_gvvs ++;
                        }

                    }
                } else break;
            }
            if (n_szfn_gvvs<1) sys->gvv_specification = 0;
        }
    } else if (key=="-gvv-align-left" || key=="--gvv-align-left" || key=="gvv-align-left" || key=="-gvv_align_left" || key=="--gvv_align_left" || key=="gvv_align_left"){
        sys->gvv_is_left_aligned = true;
    } else if (key=="-gvv-alignment-left" || key=="--gvv-alignment-left" || key=="gvv-alignment-left" || key=="-gvv_alignment_left" || key=="--gvv_alignment_left" || key=="gvv_alignment_left"){
        sys->gvv_is_left_aligned = true;
    } else if (key=="-gvv-align-right" || key=="--gvv-align-right" || key=="gvv-align-right" || key=="-gvv_align_right" || key=="--gvv_align_right" || key=="gvv_align_right"){
        sys->gvv_is_left_aligned = false;
    } else if (key=="-gvv-alignment-right" || key=="--gvv-alignment-right" || key=="gvv-alignment-right" || key=="-gvv_alignment_right" || key=="--gvv_alignment_right" || key=="gvv_alignment_right"){
        sys->gvv_is_left_aligned = false;
    } else if (key == "-zvv" || key == "--zvv" || key == "zvv" || key == "-zeta" || key == "--zeta" || key == "zeta"){
        const char * error_message = nullptr;
        int nw=0; StringNS::string sl[2]; if (i+2<argc && argv[i+1][0]!='-' && argv[i+2][0]!='-'){ sl[0] = argv[i+1]; sl[1] = argv[i+2]; nw = 2; i+=2; } else if (i+1<argc && argv[i+1][0]!='-'){ sl[0] = argv[i+1]; nw = 1; i++; }
        if (nw==2){
            if (StringNS::is_string_number(sl[0])){
                sys->drhi = atof(sl[0].text);
                if (!cp_file_with_path(szfn_zeta, sl[1].text, script_path)) error_message = "file not exist";
            } else if (StringNS::is_string_number(sl[1])){
                sys->drhi = atof(sl[1].text);
                if (!cp_file_with_path(szfn_zeta, sl[0].text, script_path)) error_message = "file not exist";
            } else {
                error_message = "please specify only one file";
            }
        } else if (nw==1){
            if (!cp_file_with_path(szfn_zeta, sl[0].text, script_path)) error_message = "file not exist";
        }
        if (error_message){ fprintf(sys->log(), "%s%s : %s[%d] : %s : %s%s %s %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, error_message, key.text[0]=='-'?"":"-", key.text, nw>=1?sl[0].text:"", nw>=2?sl[1].text:"", istty?color_string_end:""); ret = 1; }
    } else if (key == "-zeta-line" || key == "--zeta-line" || key == "-zeta_line" || key == "-zeta_line" || key == "-zeta-item" || key == "--zeta-item" || key == "-zeta_item" || key == "-zeta_item" || key == "-zeta-term" || key == "--zeta-term" || key == "-zeta_term" || key == "-zeta_term"){
        if (i+1<argc && argv[i+1][0]!='-'){
            i++; int ibegin = i;
            while (i+1<argc && (StringNS::is_string_number(argv[i+1]) || argv[i+1][0]!='-')) i++;
            int iend = i+1<argc?i+1:argc; int idx = 0;
            if (!analysis_zeta_line(sys, &argv[ibegin], &idx, iend-ibegin, script_path, ibegin)) ret = 1;
        }
   // ---------------------------------------------------------------------------------
   // The basic parameters for algorithm ----------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key == "-dr" || key == "--dr" || key == "dr"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            if (sys->drrism>0 && sys->drrism!=atof(argv[i+1])){ fprintf(sys->log(), "%s%s : %s[%d] : drrism changed to %s from %g%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i+1], sys->drrism, istty?color_string_end:""); ret = 1; }
            if (sys->drhi>0 && sys->drhi!=atof(argv[i+1])){ fprintf(sys->log(), "%s%s : %s[%d] : drhi changed to %s from %g%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i+1], sys->drhi, istty?color_string_end:""); ret = 1; }
            sys->drrism = sys->drhi = atof(argv[++i]);
        }
    } else if (key == "-drxvv" || key == "--drxvv" || key == "drxvv" || key == "-drrism" || key == "--drrism" || key == "drrism"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            if (sys->drrism>0 && sys->drrism!=atof(argv[i+1])){ fprintf(sys->log(), "%s%s : %s[%d] : %sdrrism%s changed to %s from %g%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, istty?prompt_highlight_prefix:"", istty?prompt_highlight_suffix:"", argv[i+1], sys->drrism, istty?color_string_end:""); ret = 1; }
            sys->drrism = atof(argv[++i]);
        }
    } else if (key == "-drzeta" || key == "--drzeta" || key == "drzeta" || key == "-drhi" || key == "--drhi" || key == "drhi"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            if (sys->drhi>0 && sys->drhi!=atof(argv[i+1])){ fprintf(sys->log(), "%s%s : %s[%d] : %sdrhi%s changed to %s from %g%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, istty?prompt_highlight_prefix:"", istty?prompt_highlight_suffix:"", argv[i+1], sys->drhi, istty?color_string_end:""); ret = 1; }
            sys->drhi = atof(argv[++i]);
        }
    } else if (key == "-nr" || key == "nr" || key == "--nr"){
        if (i+1<argc && argv[i+1][0]!='-'){
            bool analysis_compact = false;
            for (int j=0; argv[i+1][j] && !analysis_compact; j++) if (argv[i+1][j]==',' || argv[i+1][j]=='x' || argv[i+1][j]=='X') analysis_compact = true;
            if (analysis_compact){ i++;
                StringNS::string slp[3]; char sep[4] = { ',', 'x', 'X', 0 };
                int nwp = StringNS::analysis_general_line(sep, argv[i], slp, 3);
                if (nwp<3){
                    fprintf(sys->log(), "%s%s : %s[%d] : incorrect -nr: %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
                } else {
                    char slpt[3][64]; memset(slpt, 0, sizeof(slpt));
                    for (int i_=0; i_<3; i_++) memcpy(slpt[i_], slp[i_].text, slp[i_].length>sizeof(slpt[i_])-1?sizeof(slpt[i_])-1:slp[i_].length);
                    for (int j=0; j<3; j++) sys->nr[j] = atoi(slpt[j]);
                }
            } else {
                if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->nr[0] = sys->nr[1] = sys->nr[2] = atoi(argv[++i]);
                if (i+2<argc && StringNS::is_string_number(argv[i+1]) && StringNS::is_string_number(argv[i+2])){ sys->nr[1] = atoi(argv[++i]); sys->nr[2] = atoi(argv[++i]); }
            }
        }
    } else if (key == "-box" || key == "box" || key == "--box" || key == "-default-box" || key == "default-box" || key == "--default-box" || key == "-default_box" || key == "default_box" || key == "--default_box" || key == "-set-box" || key == "set-box" || key == "--set-box" || key == "-set_box" || key == "set_box" || key == "--set_box"){
        if (i+1<argc && argv[i+1][0]!='-'){
            if (key == "-set-box" || key == "set-box" || key == "--set-box" || key == "-set_box" || key == "set_box" || key == "--set_box"){
                sys->force_box = true;
            } else {
                sys->force_box = false;
            }
            bool analysis_compact = false;
            for (int j=0; argv[i+1][j] && !analysis_compact; j++) if (argv[i+1][j]==',' || argv[i+1][j]=='x' || argv[i+1][j]=='X') analysis_compact = true;
            if (analysis_compact){ i++;
                StringNS::string slp[3]; char sep[4] = { ',', 'x', 'X', 0 };
                int nwp = StringNS::analysis_general_line(sep, argv[i], slp, 3);
                if (nwp<3){
                    fprintf(sys->log(), "%s%s : %s[%d] : incorrect -nr: %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
                } else {
                    char slpt[3][64]; memset(slpt, 0, sizeof(slpt));
                    for (int i_=0; i_<3; i_++) memcpy(slpt[i_], slp[i_].text, slp[i_].length>sizeof(slpt[i_])-1?sizeof(slpt[i_])-1:slp[i_].length);
                    sys->box = Vector(atof(slpt[0]), atof(slpt[1]), atof(slpt[2]));
                }
            } else {
                if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->box.x = sys->box.y = sys->box.z = atof(argv[++i]);
                if (i+2<argc && StringNS::is_string_number(argv[i+1]) && StringNS::is_string_number(argv[i+2])){ sys->box.y = atof(argv[++i]); sys->box.z = atof(argv[++i]); }
            }
        }
    } else if (key == "-dielect-from-dipole" || key == "--dielect-from-dipole" || key == "dielect-from-dipole" || key == "-dielect_from_dipole" || key == "--dielect_from_dipole" || key == "dielect_from_dipole"        || key == "-dielect-rism-on" || key == "--dielect-rism-on" || key == "dielect-rism-on"){
        sys->dielect_from_dipole = true;
    } else if (key == "-dielect-use-original" || key == "--dielect-use-original" || key == "dielect-use-original" || key == "-dielect_use_original" || key == "--dielect_use_original" || key == "dielect_use_original"        || key == "-dielect-rism-off" || key == "--dielect-rism-off" || key == "dielect-rism-off"){
        sys->dielect_from_dipole = false;
    } else if (key == "-dipole-from-dielect" || key == "--dipole-from-dielect" || key == "dipole-from-dielect" || key == "-dipole_from_dielect" || key == "--dipole_from_dielect" || key == "dipole_from_dielect"){
        sys->dipole_from_dielect = true;
    } else if (key == "-dipole-use-original" || key == "--dipole-use-original" || key == "dipole-use-original" || key == "-dipole_use_original" || key == "--dipole_use_original" || key == "dipole_use_original"){
        sys->dipole_from_dielect = false;
    } else if (key == "-errtol" || key == "--errtol" || key == "errtol"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->errtolhi = sys->errtolrism = check_error_tol(atof(argv[++i]));
    } else if (key == "-errtolrism" || key == "errtolrism"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->errtolrism = check_error_tol(atof(argv[++i]));
    } else if (key == "-errtolhi" || key == "errtolhi"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->errtolhi = check_error_tol(atof(argv[++i]));
    } else if (key == "-del" || key == "del" || key == "-delvv" || key == "delvv"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->delhi = sys->delrism = atof(argv[++i]);
    } else if (key == "-delrism" || key == "delrism"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->delrism = atof(argv[++i]);
    } else if (key == "-delhi" || key == "delhi"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->delhi = atof(argv[++i]);
    } else if (key == "-ndiis" || key == "ndiis"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ndiis_rism = sys->ndiis_hi = atoi(argv[++i]);
    } else if (key == "-ndiisrism" || key == "ndiisrism"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ndiis_rism = atoi(argv[++i]);
    } else if (key == "-ndiishi" || key == "ndiishi"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ndiis_hi = atoi(argv[++i]);
   // ---------------------------------------------------------------------------------
   // The force field parameters ------------------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key == "-coulomb" || key == "--coulomb" || key == "coulomb" || key == "-coul" || key == "--coul" || key == "coul" || key == "-es" || key == "--es" || key == "es"){
        if (i+1<argc && argv[i+1][0]!='-'){ i++; StringNS::string key2 = argv[i];
            int esal_here = -1;
            for (int j=0; esal_here<0 && j<sizeof(CoulAL_names)/sizeof(CoulAL_names[0]); j++) if (key2==CoulAL_names[j]) esal_here = j;
            for (int j=0; esal_here<0 && j<sizeof(CoulAL_alias)/sizeof(CoulAL_alias[0]); j++) if (key2==CoulAL_alias[j].name) esal_here = CoulAL_alias[j].id;
            if (esal_here>=0){
                sys->esal = esal_here;
                if (sys->esal==CoulAL_Coulomb){
                    sys->esal = CoulAL_Coulomb;
                } else if (sys->esal==CoulAL_YukawaFFT){
                    if (i+1<argc && StringNS::is_string_number(argv[i+1])){
                        i++; sys->rc_yukawafft = atof(argv[i]);
                    }
                }
            } else {
                fprintf(sys->log(), "%s%s : %s[%d] : unknown Coulomb algorithm \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:"");
                ret = 1;
            }
        } else {
            if (key == "-coulomb" || key == "--coulomb" || key == "coulomb" || key == "-coul" || key == "--coul" || key == "coul"){
                sys->esal = CoulAL_Coulomb;
            }
        }
    } else if (key == "-YukawaFFT" || key == "--YukawaFFT" || key == "YukawaFFT" || key == "-YukawaFT" || key == "--YukawaFT" || key == "YukawaFT" || key == "-YukawaDFT" || key == "--YukawaDFT" || key == "YukawaDFT" || key == "-Yukawa" || key == "--Yukawa" || key == "Yukawa"){
        char command_buffer[128];
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; snprintf(command_buffer, sizeof(command_buffer), "-Coulomb Yukawa %g", atof(argv[i]));
        } else {
            snprintf(command_buffer, sizeof(command_buffer), "-Coulomb Yukawa");
        }
        ret = analysis_parameter_line(sys, command_buffer, script_name, script_line, script_path);
    //} else if (key == "-Yukawa" || key == "--Yukawa" || key == "Yukawa"){
    //    ret = analysis_parameter_line(sys, "-Coulomb Yukawa", script_name, script_line, script_path);
    //} else if (key == "-Debye-Huckel" || key == "--Debye-Huckel" || key == "Debye-Huckel" || key == "-Debye_Huckel" || key == "--Debye_Huckel" || key == "Debye_Huckel"){
    //    ret = analysis_parameter_line(sys, "-Coulomb Debye-Huckel", script_name, script_line, script_path);
    } else if (key == "-dielect-mixing" || key == "--dielect-mixing" || key == "dielect-mixing" || key == "-dielect_mixing" || key == "--dielect_mixing" || key == "dielect_mixing"){
        ret = analysis_parameter_line(sys, "-Coulomb dielect-mixing", script_name, script_line, script_path);
    } else if (key == "-dielectric-mixing" || key == "--dielectric-mixing" || key == "dielectric-mixing" || key == "-dielectric_mixing" || key == "--dielectric_mixing" || key == "dielectric_mixing"){
        ret = analysis_parameter_line(sys, "-Coulomb dielect-mixing", script_name, script_line, script_path);
    } else if (key == "-dielect-hi" || key == "--dielect-hi" || key == "dielect-hi" || key == "-dielect_hi" || key == "--dielect_hi" || key == "dielect_hi" || key == "-hi-dielect" || key == "--hi-dielect" || key == "hi-dielect" || key == "-hi_dielect" || key == "--hi_dielect" || key == "hi_dielect"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->dielect_hi = atof(argv[++i]);
    } else if (key == "-dielect" || key == "--dielect" || key == "dielect"){
        if ((i+1<argc && StringNS::is_string_number(argv[i+1]))){
            i++; double v = atof(argv[i]); for (int j=0; j<MAX_SOL; j++) sys->dielect_mol[j] = v; sys->n_dielect_mol = MAX_SOL;
            if ((i+1<argc && StringNS::is_string_number(argv[i+1]))){
                sys->n_dielect_mol = 1;
                while (i+1<argc && StringNS::is_string_number(argv[i+1])){
                    i++; sys->dielect_mol[sys->n_dielect_mol++] = atof(argv[i]); if (sys->n_dielect_mol >= MAX_SOL) break;
                }
            }
        }
    } else if (key == "-dielect-yukawa" || key == "--dielect-yukawa" || key == "dielect-yukawa" || key == "-dielect-y" || key == "--dielect-y" || key == "dielect-y" || key == "-dielect_yukawa" || key == "--dielect_yukawa" || key == "dielect_yukawa" || key == "-dielect_y" || key == "--dielect_y" || key == "dielect_y"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->dielect_yukawa = atof(argv[++i]);
            if (sys->dielect_yukawa<0) fprintf(sys->log(), "%s%s : %s[%d] : warning: invalid dielect_Yukawa %s%s\n", istty?color_string_of_synwarn:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:"");
        }
    } else if (key=="-ld-expand" || key=="--ld-expand" || key=="ld-expand" || key=="-ld_expand" || key=="--ld_expand"  || key=="ld_expand" || key=="-ld-kernel-expand" || key=="--ld-kernel-expand" || key=="ld-kernel-expand" || key=="-ld_kernel_expand" || key=="--ld_kernel_expand"  || key=="ld_kernel_expand" || key=="-ld-kernel-expanding" || key=="--ld-kernel-expanding" || key=="ld-kernel-expanding" || key=="-ld_kernel_expanding" || key=="--ld_kernel_expanding"  || key=="ld_kernel_expanding"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->ld_kernel_expand = atof(argv[++i]); if (sys->ld_kernel_expand<MACHINE_REASONABLE_ERROR) sys->ld_kernel_expand = MACHINE_REASONABLE_ERROR;
        }
    } else if (key == "-T" || key == "--T" || key == "temperature" || key == "-temperature" || key == "--temperature"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->temperature = atof(argv[++i]);
    } else if (key == "-Tdef" || key == "--Tdef" || key == "default-temperature" || key == "-default-temperature" || key == "--default-temperature" || key == "default_temperature" || key == "-default_temperature" || key == "--default_temperature"){ sys->default_temperature = atof(argv[++i]);

    } else if (key=="-external-elecfield" || key=="--external-elecfield" || key=="external-elecfield" || key=="-external_elecfield" || key=="--external_elecfield" || key=="external_elecfield"){
        if (i+3<argc && StringNS::is_string_number(argv[i+1]) && StringNS::is_string_number(argv[i+2]) && StringNS::is_string_number(argv[i+3])){
            sys->do_external_electrostatic_field = true;
            sys->external_electrostatic_field.x = atof(argv[i+1]);
            sys->external_electrostatic_field.y = atof(argv[i+2]);
            sys->external_electrostatic_field.z = atof(argv[i+3]);
            i += 3;
        }
    } else if (key=="-gcutoff" || key=="--gcutoff" || key=="gcutoff" || key=="-ccutoff" || key=="--ccutoff" || key=="ccutoff"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->ccutoff = atof(argv[++i]);
    } else if (key=="-gceil" || key=="--gceil" || key=="gceil" || key=="-cceil" || key=="--cceil" || key=="cceil" || key=="-closure-ceil" || key=="--closure-ceil" || key=="closure-ceil" || key=="-closure_ceil" || key=="--closure_ceil" || key=="closure_ceil"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            double exp_cutoff = atof(argv[++i]);
            if (exp_cutoff>0) sys->ccutoff = log(exp_cutoff);
            else { fprintf(sys->log(), "%s%s : %s[%d] : invalid HNC ceil %s%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1; }
        }
    } else if (key=="-dynamic-delvv" || key=="--dynamic-delvv" || key=="dynamic-delvv" || key=="-dynamic_delvv" || key=="--dynamic_delvv" || key=="dynamic_delvv" || key=="-enhance-closure" || key=="--enhance-closure" || key=="enhance-closure" || key=="-enhance_closure" || key=="--enhance_closure" || key=="enhance_closure" || key=="-closure-enhance" || key=="--closure-enhance" || key=="closure-enhance" || key=="-closure_enhance" || key=="--closure_enhance" || key=="closure_enhance" || key=="-closure-enhancement" || key=="--closure-enhancement" || key=="closure-enhancement" || key=="-closure_enhancement" || key=="--closure_enhancement" || key=="closure_enhancement"){
        sys->closure_enhance_level = 1;
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->closure_enhance_level = atof(argv[i]);
        }
    } else if (key=="-no-dynamic-delvv" || key=="--no-dynamic-delvv" || key=="no-dynamic-delvv" || key=="-no_dynamic_delvv" || key=="--no_dynamic_delvv" || key=="no_dynamic_delvv" || key=="-no-enhance-closure" || key=="--no-enhance-closure" || key=="no-enhance-closure" || key=="-no_enhance_closure" || key=="--no_enhance_closure" || key=="no_enhance_closure" || key=="-no-closure-enhance" || key=="--no-closure-enhance" || key=="no-closure-enhance" || key=="-no_closure_enhance" || key=="--no_closure_enhance" || key=="no_closure_enhance" || key=="-no-closure-enhancement" || key=="--no-closure-enhancement" || key=="no-closure-enhancement" || key=="-no_closure_enhancement" || key=="--no_closure_enhancement" || key=="no_closure_enhancement"){
        sys->closure_enhance_level = 0;
    } else if (key=="-hlr-no-hi" || key=="--hlr-no-hi" || key=="hlr-no-hi" || key=="-hlr_no_hi" || key=="--hlr_no_hi" || key=="-hlr_no_hi"){
        sys->hlr_no_hi = true;
    } else if (key=="-hlr-with-hi" || key=="--hlr-with-hi" || key=="hlr-with-hi" || key=="-hlr_with_hi" || key=="--hlr_with_hi" || key=="-hlr_with_hi"){
        sys->hlr_no_hi = false;
    } else if (key=="-homo-iet" || key=="--homo-iet" || key=="homo-iet" || key=="-homo_iet" || key=="--homo_iet" || key=="homo_iet"){
        sys->use_homogeneous_rism = true;
    } else if (key=="-inhomo-iet" || key=="--inhomo-iet" || key=="inhomo-iet" || key=="-inhomo_iet" || key=="--inhomo_iet" || key=="inhomo_iet"){
        sys->use_homogeneous_rism = false;
    } else if (key=="-rb" || key=="rb" || key=="--rb" || key=="-rbohr" || key=="rbohr" || key=="--rbohr" || key=="-Bohr-radius" || key=="--Bohr-radius" || key=="Bohr-radius" || key=="-Bohr_radius" || key=="--Bohr_radius" || key=="Bohr_radius" ){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->rbohr = atof(argv[++i]);
        }
    } else if (key=="-rq" || key=="rq" || key=="--rq" || key=="-rcharge" || key=="rcharge" || key=="--rcharge"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            sys->rbcharge = atof(argv[++i]);
        }
    } else if (key=="-xvv-extend" || key=="--xvv-extend" || key=="xvv-extend" || key=="-xvv_extend" || key=="--xvv_extend" || key=="xvv_extend"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->xvv_extend = atof(argv[i]); }
    } else if (key=="-xvv-scheme" || key=="--xvv-scheme" || key=="xvv-scheme"){
        if (i+1<argc && argv[i+1][0]!='-'){
            i++; sys->xvv_scheme = argv[i];
        }
    } else if (key=="-hvv-scale" || key=="--hvv-scale" || key=="hvv-scale" || key=="-xvv-scale" || key=="--xvv-scale" || key=="xvv-scale" || key=="-scale-hvv" || key=="--scale-hvv" || key=="scale-hvv" || key=="-scale-xvv" || key=="--scale-xvv" || key=="scale-xvv"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; int i_begin = i; int i_end = i+1; while (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; i_end = i+1; }
            sys->n_xvv_scale = i_end-i_begin; if (sys->xvv_scale) free(sys->xvv_scale); sys->xvv_scale = (double*) malloc(sizeof(double)*(sys->n_xvv_scale+1)); if (!sys->xvv_scale){ fprintf(sys->log(), "%s%s : malloc failure%s\n", istty?color_string_of_synerr:"", software_name, istty?color_string_end:""); exit(-1); }
            for (int j=0; j<sys->n_xvv_scale; j++) sys->xvv_scale[j] = atof(argv[i_begin+j]);
        }
  #ifdef _INTERACTIVE_
    } else if (key=="-interactive" || key=="--interactive" || key=="interactive" || key=="-intervene" || key=="--intervene" || key=="intervene"){
        sys->allow_interactive = true;
    } else if (key=="-no-interactive" || key=="--no-interactive" || key=="no-interactive" || key=="-no_interactive" || key=="--no_interactive" || key=="no_interactive" || key=="-no-intervene" || key=="--no-intervene" || key=="no-intervene" || key=="-no_intervene" || key=="--no_intervene" || key=="no_intervene"){
        sys->allow_interactive = false;
  #endif
    } else if (key=="-pbc" || key=="--pbc" || key=="pbc"){
        if (i+1<argc && argv[i+1][0]!='-'){ i++; StringNS::string pbc_str = argv[i];
            if (pbc_str == "no" || pbc_str == "none"){
                sys->pbc_x = sys->pbc_y = sys->pbc_z = false;
            } else if (pbc_str == "all"){
                sys->pbc_x = sys->pbc_y = sys->pbc_z = true;
            } else {
                sys->pbc_x = sys->pbc_y = sys->pbc_z = false;
                for (int i=0; i<pbc_str.length; i++){
                  if (pbc_str.text[i]=='x' || pbc_str.text[i]=='X') sys->pbc_x = true;
                  else if (pbc_str.text[i]=='y' || pbc_str.text[i]=='Y') sys->pbc_y = true;
                  else if (pbc_str.text[i]=='z' || pbc_str.text[i]=='Z') sys->pbc_z = true;
                  else { fprintf(sys->log(), "%s%s : %s[%d] : unknown type of pbc \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, pbc_str.text, istty?color_string_end:""); ret = 1; }
                }
            }
        }
    } else if (key=="-nopbc" || key=="--nopbc" || key=="nopbc"){      sys->pbc_x = sys->pbc_y = sys->pbc_z = false;
    } else if (key=="-density" || key=="--density" || key=="density"){ sys->nmv = 0;
        while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->density_mv[sys->nmv++] = atof(argv[i]); if (sys->nmv >= MAX_SOL) break;
        }
    } else if (key=="-bulkdensity" || key=="--bulkdensity" || key=="bulkdensity" || key=="-bulk-density" || key=="--bulk-density" || key=="bulk-density" || key=="-bulk_density" || key=="--bulk_density" || key=="bulk_density"){
        sys->nmvb = 0;
        while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->bulk_density_mv[sys->nmvb++] = atof(argv[i]); if (sys->nmvb >= MAX_SOL) break;
        }
    } else if (key == "-rlj" || key == "--rlj" || key == "rlj" || key == "-rvdw" || key == "--rvdw" || key == "rvdw"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->rvdw = atof(argv[++i]);
    } else if (key == "-rcoul" || key == "--rcoul" || key == "rcoul"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->rcoul = atof(argv[++i]);
    } else if (key == "-rlocal-coul" || key == "--rlocal-coul" || key == "rlocal-coul" || key == "-rlocal_coul" || key == "--rlocal_coul" || key == "rlocal_coul"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->rlocal_coul = atof(argv[++i]);
    } else if (key == "-rc" || key == "--rc" || key == "rc"){ if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->rvdw = sys->rcoul = atof(argv[++i]);
    } else if (key == "-pme-gamma" || key == "--pme-gamma" || key == "pme-gamma" || key == "-pme_gamma" || key == "--pme_gamma" || key == "pme_gamma"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])) sys->gamma_erf = atof(argv[++i]); sys->gamma_auto_generate = false;
    } else if (key == "-arithmetic-sigma" || key == "arithmetic-sigma" || key == "--arithmetic-sigma" || key == "-arith-sigma" || key == "arith-sigma" || key == "--arith-sigma" || key == "-arithmetic_sigma" || key == "arithmetic_sigma" || key == "--arithmetic_sigma" || key == "-arith_sigma" || key == "arith_sigma" || key == "--arith_sigma"){
        sys->mixingrule_sigma_geometric = false; sys->mixingrule_sigma_geometric_specified = true;
    } else if (key == "-geometric-sigma" || key == "geometric-sigma" || key == "--geometric-sigma" || key == "-geo-sigma" || key == "geo-sigma" || key == "--geo-sigma" || key == "-geometric_sigma" || key == "geometric_sigma" || key == "--geometric_sigma" || key == "-geo_sigma" || key == "geo_sigma" || key == "--geo_sigma"){
        sys->mixingrule_sigma_geometric = true; sys->mixingrule_sigma_geometric_specified = true;
    } else if (key == "-allow"||key == "--allow"||key == "allow" || key == "-enable"||key == "--enable"||key == "enable"){
        if (i+1<argc && argv[i+1][0]!='-') if (!set_sys_options(sys, argv[++i], true, script_name, script_line)) ret = 1;
    } else if (key == "-forbid"||key == "--forbid"||key == "forbid" || key == "-disable"||key == "--disable"||key == "disable"){
        if (i+1<argc && argv[i+1][0]!='-') if (!set_sys_options(sys, argv[++i], false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,7) == "-allow-" || key.sub(0,7) == "-allow_"){
        if (!set_sys_options(sys, key.sub(7).text, true, script_name, script_line)) ret = 1;
    } else if (key.sub(0,8) == "--allow-" || key.sub(0,8) == "--allow_"){
        if (!set_sys_options(sys, key.sub(8).text, true, script_name, script_line)) ret = 1;
    } else if (key.sub(0,8) == "-enable-" || key.sub(0,8) == "-enable_"){
        if (!set_sys_options(sys, key.sub(8).text, true, script_name, script_line)) ret = 1;
    } else if (key.sub(0,9) == "--enable-" || key.sub(0,9) == "--enable_"){
        if (!set_sys_options(sys, key.sub(9).text, true, script_name, script_line)) ret = 1;
    } else if (key.sub(0,8) == "-forbid-" || key.sub(0,8) == "-forbid_"){
        if (!set_sys_options(sys, key.sub(8).text, false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,9) == "--forbid-" || key.sub(0,9) == "--forbid_"){
        if (!set_sys_options(sys, key.sub(9).text, false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,9) == "-disable-" || key.sub(0,9) == "-disable_"){
        if (!set_sys_options(sys, key.sub(9).text, false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,10) == "--disable-" || key.sub(0,10) == "--disable_"){
        if (!set_sys_options(sys, key.sub(10).text, false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,6) == "-skip-" || key.sub(0,6) == "-skip_"){
        if (!set_sys_options(sys, key.sub(6).text, false, script_name, script_line)) ret = 1;
    } else if (key.sub(0,7) == "--skip-" || key.sub(0,7) == "--skip_"){
        if (!set_sys_options(sys, key.sub(0,7).text, false, script_name, script_line)) ret = 1;
   // ---------------------------------------------------------------------------------
   // The iteration control parameters ------------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key == "-stepmax" || key == "--stepmax" || key == "stepmax" || key == "-maxstep" || key == "--maxstep" || key == "maxstep"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->stepmax_rism = sys->stepmax_hi = atoi(argv[i]);
        }
    } else if (key == "-step" || key == "--step" || key == "step" || key == "-steps" || key == "--steps" || key == "steps"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->stepmax_hi = sys->stepmax_rism = atoi(argv[i]);
            if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->stepmax_rism = atoi(argv[i]); }
        }
    } else if (key=="-dof" || key=="--dof" || key=="dof" || key=="-degree-of-freedom" || key=="-degree-of-freedom" || key=="degree-of-freedom" || key=="-degree_of_freedom" || key=="--degree_of_freedom" || key=="degree_of_freedom"){
        sys->n_degree_of_freedom = 0; while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; if (sys->n_degree_of_freedom<MAX_SOL) sys->degree_of_freedom[sys->n_degree_of_freedom++] = atof(argv[i]);
        }
    } else if (key=="-scale.zeta" || key=="--scale.zeta" || key=="scale.zeta"){
        sys->n_zeta_scaling_factor = 0; while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; if (sys->n_zeta_scaling_factor<MAX_SOL) sys->zeta_scaling_factor[sys->n_zeta_scaling_factor++] = atof(argv[i]);
        }
    } else if (key=="-ucutoff-hi" || key=="--ucutoff-hi" || key=="ucutoff-hi" || key=="-ucutoff_hi" || key=="--ucutoff_hi" || key=="ucutoff_hi" || key=="-ucutoff-hshi" || key=="--ucutoff-hshi" || key=="ucutoff-hshi" || key=="-ucutoff_hshi" || key=="--ucutoff_hshi" || key=="ucutoff_hshi"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->ucutoff_hs = atof(argv[i]);
        }
    } else if (key=="-ccutoff-hi" || key=="--ccutoff-hi" || key=="ccutoff-hi" || key=="-ccutoff_hi" || key=="--ccutoff_hi" || key=="ccutoff_hi" || key=="-ccutoff-hshi" || key=="--ccutoff-hshi" || key=="ccutoff-hshi" || key=="-ccutoff_hshi" || key=="--ccutoff_hshi" || key=="ccutoff_hshi" || key=="-ccutoff-theta" || key=="--ccutoff-theta" || key=="ccutoff-theta" || key=="-ccutoff_theta" || key=="--ccutoff_theta" || key=="ccutoff_theta"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->ccutoff_hi = atof(argv[i]);
        }
    } else if (key=="-dipole" || key=="--dipole" || key=="dipole"){ sys->ndipole = 0;
        while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->dipole_mv[sys->ndipole++] = atof(argv[i]); if (sys->ndipole >= MAX_SOL) break;
        }
    } else if (key == "-do" || key == "--do" || key == "do" || key == "-cmd" || key == "--cmd" || key == "cmd" || key == "-command" || key == "--command" || key == "command"){
        const char * line_start = argv[i];
        while (i+1<argc && argv[i+1][0]!='-'){
            i++;
            char argv_copy[4096]; memset(argv_copy, 0, sizeof(argv_copy)); strncpy(argv_copy, argv[i], sizeof(argv_copy)-1);
            bool analysis_success = analysis_command(sys, argv_copy, line_start, script_name, script_line, &argv[i][0]-line_start);
            if (!analysis_success) ret = 1;
        }
    } else if (key == "-run" || key == "--run" || key == "run"){
        if (i+1<argc){ i++; char filename[MAX_PATH]; memset(filename, 0, sizeof(filename)); strncpy(filename, argv[i], sizeof(filename)-2);
            ret = analysis_params_file(sys, filename);
        }
   // ---------------------------------------------------------------------------------
   // Advanced settings ---------------------------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key=="-v" || key=="--v" || key=="-verbose" || key=="verbose" || key=="--verbose"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->detail_level = atoi(argv[i]); }
        else sys->detail_level = 2;
    } else if (key=="-detail-level" || key=="detail-level" || key=="--detail-level" || key=="-detail_level" || key=="detail_level" || key=="--detail_level"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->detail_level = atoi(argv[i]); }
    } else if (key=="-detail" || key=="detail" || key=="--detail"){
        sys->detail_level = 2; if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->detail_level = atoi(argv[i]); }
    } else if (key=="-brief" || key=="brief"  || key=="--brief"){
        sys->detail_level = 1;
    } else if (key=="-silent" || key=="silent" || key=="--silent"){
        sys->detail_level = 0;
    } else if (key=="-detailed" || key=="detailed" || key=="--detailed" || key=="-details" || key=="details" || key=="--details"){
        sys->detail_level = 2;
    } else if (key=="-debug-level" || key=="debug-level"  || key=="--debug-level" || key=="-debug_level" || key=="debug_level"  || key=="--debug_level"){
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->debug_level = atoi(argv[i]); }
    } else if (key=="-debugging" || key=="debugging"  || key=="--debugging" || key=="-debuging" || key=="debuging"  || key=="--debuging"){
        sys->debug_level = 2;
    } else if (key=="-debug" || key=="debug"  || key=="--debug"){
        sys->debug_level = 1; if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; sys->debug_level = atoi(argv[i]); }
    } else if (key=="-debug-crc" || key=="debug-crc" || key=="--debug-crc" || key=="-debug_crc" || key=="debug_crc" || key=="--debug_crc" || key=="-check-crc" || key=="check-crc" || key=="--check-crc" || key=="-check_crc" || key=="check_crc" || key=="--check_crc"){
        sys->debug_show_crc = true;
    } else if (key=="-no-debug-crc" || key=="no-debug-crc" || key=="--no-debug-crc" || key=="-no_debug_crc" || key=="no_debug_crc" || key=="--no_debug_crc"){
        sys->debug_show_crc = false;
    } else if (key=="-debug-xvv" || key=="debug-xvv" || key=="--debug-xvv" || key=="-debug_xvv" || key=="debug_xvv" || key=="--debug_xvv" || key=="-check-xvv" || key=="check-xvv" || key=="--check-xvv" || key=="-check_xvv" || key=="check_xvv" || key=="--check_xvv"){
        sys->debug_xvv = true;
        if (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->debug_xvv_terms = atoi(argv[i]);
        }
    } else if (key=="-no-debug-xvv" || key=="no-debug-xvv" || key=="--no-debug-xvv" || key=="-no_debug_xvv" || key=="no_debug_xvv" || key=="--no_debug_xvv"){
        sys->debug_xvv = false;
    } else if (key=="-test" || key=="test"  || key=="--test"){              sys->mode_test = true;
    } else if (key=="-nodebug" || key=="nodebug"  || key=="--nodebug"){     sys->debug_level = 0;
    } else if (key == "-list" || key == "--list" || key == "list" || key == "-listonly" || key == "--listonly" || key == "listonly" || key == "-list-only" || key == "--list-only" || key == "list-only" || key == "-list_only" || key == "--list_only" || key == "list_only"){
        sys->listonly = true; sys->listall = true;
    } else if (key == "-ls" || key == "--ls" || key == "ls"){
        sys->listonly = true; sys->listall = false;
    } else if (key=="-lambda" || key=="--lambda" || key=="lambda" || key=="-llambda" || key=="--llambda" || key=="llambda"){
        bool is_lambda = (key=="-lambda" || key=="--lambda" || key=="lambda"); double logEE = log(EE);
        while (i+1<argc && StringNS::is_string_number(argv[i+1])){
            i++; sys->llambda[sys->nllambda] = atof(argv[i]);
            if (is_lambda) sys->llambda[sys->nllambda] = log(fabs(sys->llambda[sys->nllambda]) + 1e-120)/logEE;
            sys->nllambda ++; if (sys->nllambda >= MAX_SOL) break;
        }
    } else if (key=="-lsa" || key=="--lsa" || key=="lsa" || key=="-lse_a" || key=="--lse_a" || key=="lse_a" || key=="-lse-a" || key=="--lse-a" || key=="lse-a" || key=="-les_a" || key=="--les_a" || key=="les_a" || key=="-les-a" || key=="--les-a" || key=="les-a"){
        if (i+1<argc && argv[i+1][0] != '-'){ i++;
            if (StringNS::is_string_number(argv[i])){
                sys->lse_a = atof(argv[i]); if (sys->calc_ab_automatically==1) sys->calc_ab_automatically = 0;
            } else if (StringNS::string(argv[i]) == "auto"){
                sys->lse_a = 0.3; sys->calc_ab_automatically = 1;
            }
        }
    } else if (key=="-lsb" || key=="-lsb" || key=="lsb" || key=="-lse_b" || key=="-lse_b" || key=="lse_b" || key=="-lse-b" || key=="--lse-b" || key=="lse-b" || key=="-les_b" || key=="-les_b" || key=="les_b" || key=="-les-b" || key=="--les-b" || key=="les-b"){
        if (i+1<argc && argv[i+1][0] != '-'){ i++;
            if (StringNS::is_string_number(argv[i])){
                sys->lse_b = atof(argv[i]); if (sys->calc_ab_automatically==2) sys->calc_ab_automatically = 0;
            } else if (StringNS::string(argv[i]) == "auto"){
                sys->lse_b = 50; sys->calc_ab_automatically = 2;
            }
        }
   // ---------------------------------------------------------------------------------
   // Shortcuts -----------------------------------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key == "-ff" || key == "--ff" || key == "ff"){
        if (i+1<argc && argv[i+1][0]!='-'){ i++; StringNS::string key2 = argv[i];
            if (key2 == "amber"){
                sys->forcefield_prefix = FFPREFIX_AMBER;
                ret = analysis_parameter_line(sys, "-Tdef 502.97", script_name, script_line, script_path);
                ret = analysis_parameter_line(sys, "-arith-sigma", script_name, script_line, script_path);
            } else if (key2 == "opls" || key2 == "oplsaa"){
                sys->forcefield_prefix = FFPREFIX_OPLS;
                ret = analysis_parameter_line(sys, "-Tdef 120.27", script_name, script_line, script_path);
                ret = analysis_parameter_line(sys, "-geo-sigma", script_name, script_line, script_path);
            } else if (key2 == "gaff"){
                sys->forcefield_prefix = FFPREFIX_GAFF;
                ret = analysis_parameter_line(sys, "-Tdef 120.27", script_name, script_line, script_path);
                ret = analysis_parameter_line(sys, "-arith-sigma", script_name, script_line, script_path);
            } else {
                fprintf(sys->log(), "%s%s : %s[%d] : unknown forcefield specification: \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:"");
                ret = 1;
            }
        }
    } else if (key == "-ffamber" || key == "--ffamber" || key == "ffamber" || key == "-amberff" || key == "--amberff" || key == "amberff"){
        ret = analysis_parameter_line(sys, "-Tdef 502.97", script_name, script_line, script_path);
        ret = analysis_parameter_line(sys, "-arith-sigma", script_name, script_line, script_path);
    } else if (key == "-ffopls" || key == "--ffopls" || key == "ffopls" || key == "-oplsff" || key == "--oplsff" || key == "oplsff" || key == "-ffoplsaa" || key == "--ffoplsaa" || key == "ffoplsaa" || key == "-oplsaaff" || key == "--oplsaaff" || key == "oplsaaff"){
        ret = analysis_parameter_line(sys, "-Tdef 120.27", script_name, script_line, script_path);
        ret = analysis_parameter_line(sys, "-geo-sigma", script_name, script_line, script_path);
    } else if (key == "-ffgmaff" || key == "--ffgmaff" || key == "ffgmaff" || key == "-gaff" || key == "--gaff" || key == "gaff"){
        ret = analysis_parameter_line(sys, "-Tdef 120.27", script_name, script_line, script_path);
        ret = analysis_parameter_line(sys, "-arith-sigma", script_name, script_line, script_path);
    } else if (key == "-Lorentz-Berthelot" || key == "--Lorentz-Berthelot" || key == "Lorentz-Berthelot" || key == "-Lorentz_Berthelot" || key == "--Lorentz_Berthelot" || key == "Lorentz_Berthelot"){
        ret = analysis_parameter_line(sys, "-arith-sigma", script_name, script_line, script_path);
    } else if (key == "-do-rism-kh" || key == "-do-kh-rism" || key == "-do-khrism"){
        ret = analysis_parameter_line(sys, "-cmd closure=kh rism report:energy,rdf savee:cmd,guv,cuv", script_name, script_line, script_path);
    } else if (key == "-do-rismhi-d2" || key == "-do-d2-rismhi" || key == "-do-d2rismhi"){
        ret = analysis_parameter_line(sys, "-cmd hshi closure=d2 rism report:energy,rdf savee:cmd,guv,cuv", script_name, script_line, script_path);
  // ---------------------------------------------------------------------------------
  // Shortcuts for old versions ------------------------------------------------------
  // ---------------------------------------------------------------------------------
    } else if (key=="-skip-zero-xvv" || key=="--skip-zero-xvv" || key=="skip-zero-xvv" || key=="-skip_zero_xvv" || key=="--skip_zero_xvv" || key=="skip_zero_xvv" || key=="-skip-missing-xvv" || key=="--skip-missing-xvv" || key=="skip-missing-xvv" || key=="-skip_missing_xvv" || key=="--skip_missing_xvv" || key=="skip_missing_xvv" || key=="-xvv-skip-missing" || key=="--xvv-skip-missing" || key=="xvv-skip-missing" || key=="--xvv_skip_missing" || key=="--xvv_skip_missing" || key=="xvv_skip_missing"){
        ret = analysis_parameter_line(sys, "-allow missing-xvv", script_name, script_line, script_path);
    } else if (key=="-zeta-forbid-missing" || key=="--zeta-forbid-missing" || key=="zeta-forbid-missing" || key=="-forbid-zeta-missing" || key=="--forbid-zeta-missing" || key=="forbid-zeta-missing" || key=="-zeta_forbid_missing" || key=="--zeta_forbid_missing" || key=="zeta_forbid_missing" || key=="-forbid_zeta_missing" || key=="--forbid_zeta_missing" || key=="forbid_zeta_missing"){
        ret = analysis_parameter_line(sys, "-forbid missing-zeta", script_name, script_line, script_path);
    } else if (key == "-ignore-ram" || key == "--ignore-ram" || key == "ignore-ram" || key == "-ignore_ram" || key == "--ignore_ram" || key == "ignore_ram" || key == "-ignore-memory-capacity" || key == "--ignore-memory-capacity" || key == "ignore-memory-capacity" || key == "-ignore_memory_capacity" || key == "--ignore_memory_capacity" || key == "ignore_memory_capacity"){
        ret = analysis_parameter_line(sys, "-allow exceed-memory-capacity", script_name, script_line, script_path);
    } else if (key == "-bound-to-ram" || key == "--bound-to-ram" || key == "bound-to-ram" || key == "-bound_to_ram" || key == "--bound_to_ram" || key == "bound_to_ram" || key == "-check-memory-capacity" || key == "--check-memory-capacity" || key == "check-memory-capacity" || key == "-check_memory_capacity" || key == "--check_memory_capacity" || key == "check_memory_capacity"){
        ret = analysis_parameter_line(sys, "-forbid exceed-memory-capacity", script_name, script_line, script_path);
   // ---------------------------------------------------------------------------------
   // Done for all known parameters ---------------------------------------------------
   // ---------------------------------------------------------------------------------
    } else if (key.text[0] == '-'){
        #ifdef _EXPERIMENTAL_
            if (analysis_experinemtal_parameter_line(sys, argv, &i, argc, script_name, script_line, script_path)>0){
                fprintf(sys->log(), "%s%s : %s[%d] : unknown option \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
            }
        #else
            fprintf(sys->log(), "%s%s : %s[%d] : unknown option \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
        #endif
    } else {
        fprintf(sys->log(), "%s%s : %s[%d] : unknown string \"%s\"%s\n", istty?color_string_of_synerr:"", software_name, get_second_fn(script_name), script_line, argv[i], istty?color_string_end:""); ret = 1;
    }
    *argi = i;
    return ret;
}
int analysis_parameter_line(IET_Param * sys, const char * in, char * script_name, int script_line, const char * script_path){
    char input[4096]; strcpy(input, in);
    StringNS::string sl[MAX_WORD];
    int nw = analysis_line_params(input, sl, MAX_WORD, true); if (nw<=0) return 0;
    char * sltxt[MAX_WORD];
    for (int i=0; i<nw; i++){ sl[i].text[sl[i].length] = 0; sltxt[i] = sl[i].text; }
    int idx = 0;
    return analysis_parameter_line(sys, sltxt, &idx, nw, script_name, script_line, script_path);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int analysis_params_file(IET_Param * sys, char * filename, int recursive_layer){
//    if (recursive_layer > MAX_INCLUDE_RECURSIVE) return -1;
    int error = 0; bool success = true; bool istty = sys->is_log_tty;
    if (filename && filename[0]){
        int idx = 0; int nline = 0; int cmd_compile = 1;
        StringNS::string sl[MAX_WORD]; char * sltxt[MAX_WORD]; char input[4096]; memset(input, 0, sizeof(input));
        FILE * finfo = (StringNS::string(filename)=="screen" || StringNS::string(filename)=="con" || StringNS::string(filename)=="stdin")? stdin : fopen(filename, "r");
        if (finfo){
            if (sys->debug_level>=1 && sys->log()) fprintf(sys->log(), "debug:: analysis_params_file(%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", filename, sys->is_log_tty?prompt_path_suffix:"\"");
            char info_file_path[MAX_PATH]; info_file_path[0] = 0; if (finfo!=stdin) get_file_path(filename, info_file_path); bool in_block_comment = false;
            while (fgets(input, sizeof(input)-1, finfo)){ nline++;

                StringNS::string string_line = input;
                if (in_block_comment){
                    string_line.length = 0;
                    for (int i=0; i<sizeof(input)-1&&input[i]; i++){
                        if (input[i]=='*'&&input[i+1]=='/'){
                            in_block_comment = false; string_line = &input[i+2]; break;
                        }
                    }
                } else {
                    for (int i=0; i<sizeof(input)-1&&input[i]; i++){
                        if (input[i]=='/'&&input[i+1]=='*'){
                            in_block_comment = true; string_line.length = i; break;
                        }
                    }
                }
                int nw = analysis_line_params(string_line, sl, MAX_WORD, true); if (nw<=0) continue ; //int nw = analysis_line_params(input, sl, MAX_WORD, true); if (nw<=0) continue ;
                for (int i=0; i<nw; i++){ sl[i].text[sl[i].length] = 0; sltxt[i] = sl[i].text; } idx = 0;
                if (sl[0].text[0] == ';'){ continue;
                } else if ((sl[0]=="#echo"||sl[0]=="@echo") || (nw>1 && (sl[0]=="#"||sl[0]=="@") && sl[1]=="echo")){
                    if (istty) fprintf(sys->log(), "%s", color_string_of_echo);
                    for (int i=((sl[0]=="#"||sl[0]=="@")?2:1); i<nw; i++) fprintf(sys->log(), "%s%s", sl[i].text, i==nw-1?"":" ");
                    if (istty) fprintf(stdout, "%s\n", color_string_end); else fprintf(sys->log(), "\n");
                } else if ((sl[0]=="#warning"||sl[0]=="@warning") || (nw>1 && (sl[0]=="#"||sl[0]=="@") && sl[1]=="warning")){
                    fprintf(sys->log(), "%s : %s[%d]: WARNING: %s", software_name, get_second_fn(filename), nline, istty?color_string_of_warning:"");
                    for (int i=((sl[0]=="#"||sl[0]=="@")?2:1); i<nw; i++) fprintf(sys->log(), "%s%s", sl[i].text, i==nw-1?"":" ");
                    fprintf(sys->log(), "%s\n", istty?color_string_end:"");
                } else if ((sl[0]=="#error"||sl[0]=="@error") || (nw>1 && (sl[0]=="#"||sl[0]=="@") && sl[1]=="error")){
                    fprintf(sys->log(), "%s : %s[%d]: ERROR: %s", software_name, get_second_fn(filename), nline, istty?color_string_of_error:"");
                    for (int i=((sl[0]=="#"||sl[0]=="@")?2:1); i<nw; i++) fprintf(sys->log(), "%s%s", sl[i].text, i==nw-1?"":" ");
                    fprintf(sys->log(), "%s\n", istty?color_string_end:"");
                    error = -1;
                } else if ((nw>1 && (sl[0]=="include" || sl[0]=="#include" || sl[0]=="@include")) || (nw>2 && (sl[0]=="#"||sl[0]=="@") && sl[1]=="include")){
                    char filenameinclude[MAX_PATH]; filenameinclude [0] = 0;
                    if (sl[0]=="#" || sl[0]=="@") strcpy(filenameinclude, sl[2].text); else strcpy(filenameinclude, sl[1].text);
                    if (recursive_layer >= MAX_INCLUDE_RECURSIVE){
                        fprintf(sys->log(), "%s : %s[%d]: error: too many includes\n", software_name, get_second_fn(filename), nline);
                        error = -1;
                    } else {
                        FILE * ftmp = fopen(filenameinclude, "r");
                        if (ftmp){ fclose(ftmp);
                            error = analysis_params_file(sys, filenameinclude, recursive_layer+1);
                        } else {
                            fprintf(sys->log(), "%s%s : %s[%d]: error: cannot open include file %s%s%s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(filename), nline, istty?prompt_path_prefix:"", filenameinclude, istty?prompt_path_suffix:"", sys->is_log_tty?color_string_end:"");
                            error = -1;
                        }
                    }
                } else if (sl[0].text[0] == '#' || sl[0].text[0] == ';'){
                    continue;
                } else if (sl[0].text[0] == '[' || sl[0].text[sl[0].length-1] == ':'){
                    if (nw>1 && sl[0] == "["){
                        if (sl[1] == "iet"){ cmd_compile = 1;
                        } else if (sl[1] == "solvent"){ cmd_compile = 1;
                        } else if (sl[1] == "rism"){ cmd_compile = 1;
                        } else if (sl[1] == "rism3d"){ cmd_compile = 1;
                        } else if (sl[1] == "iethi"){ cmd_compile = 1;
                        } else if (sl[1] == "iet-hi"){ cmd_compile = 1;
                        } else if (sl[1] == "iet_hi"){ cmd_compile = 1;
                        } else if (sl[1] == "rismhi"){ cmd_compile = 1;
                        } else if (sl[1] == "rismhi3d"){ cmd_compile = 1;
                        } else if (sl[1] == "hi"){ cmd_compile = 1;
                        } else if (sl[1] == "atom"){ cmd_compile = 2;
                        } else if (sl[1] == "bond"){ cmd_compile = 3;
                        } else if (sl[1] == "pair"){ cmd_compile = 4;
                        } else if (sl[1] == "gvv-map"){ cmd_compile = 5;
                        } else if (sl[1] == "gvv_map"){ cmd_compile = 5;
                        } else {
                            cmd_compile = 0; fprintf(sys->log(), "%s : %s[%d]: warning: skipping section [%s]\n", software_name, get_second_fn(filename), nline, sl[1].text);
                        }
                    } else if (sl[0]=="[solvent]" || sl[0]=="solvent:"){ cmd_compile = 1;
                    } else if (sl[0]=="[iet]"     || sl[0]=="iet:"    ){ cmd_compile = 1;
                    } else if (sl[0]=="[hi]"      || sl[0]=="hi:"     ){ cmd_compile = 1;
                    } else if (sl[0]=="[iethi]"   || sl[0]=="iethi:"  ){ cmd_compile = 1;
                    } else if (sl[0]=="[rism]"    || sl[0]=="rism:"   ){ cmd_compile = 1;
                    } else if (sl[0]=="[rismhi]"  || sl[0]=="rismhi:" ){ cmd_compile = 1;
                    } else if (sl[0]=="[atom]"    || sl[0]=="atom:"   ){ cmd_compile = 2;
                    } else if (sl[0]=="[bond]"    || sl[0]=="bond:"   ){ cmd_compile = 3;
                    } else if (sl[0]=="[pair]"    || sl[0]=="pair:"   ){ cmd_compile = 4;
                    } else if (sl[0]=="[gvv-map]" || sl[0]=="gvv-map:"){ cmd_compile = 5;
                    } else if (sl[0]=="[gvv_map]" || sl[0]=="gvv_map:"){ cmd_compile = 5;
                    } else if (sl[0]=="[zeta]"    || sl[0]=="zeta:"   ){ cmd_compile = 6;
                    } else {
                        cmd_compile = 0; fprintf(sys->log(), "%s : %s[%d]: warning: skipping section %s\n", software_name, get_second_fn(filename), nline, sl[0].text);
                    }
                } else if (cmd_compile == 1) {
                    int _error = analysis_parameter_line(sys, sltxt, &idx, nw, filename, nline, info_file_path[0]?info_file_path:nullptr);
                    if (_error){ success = false; error |= _error; }
                } else if (cmd_compile == 2) {  // [atom]
                    success = analysis_solvent_atom(sys, sltxt, &idx, nw, filename, nline);
                    if (!success) error = -1;
                } else if (cmd_compile == 3 || cmd_compile == 4) {  // [bond]
                    success = analysis_solvent_bond(sys, sltxt, &idx, nw, filename, nline, cmd_compile==3?"bond":"pair");
                    if (!success) error = -1;
                } else if (cmd_compile == 5){   // [gvv-map]
                    success = analysis_gvv_map(sys, sltxt, &idx, nw, filename, nline);
                    if (!success) error = -1;
                } else if (cmd_compile == 6){   // [zeta]
                    success = analysis_zeta_line(sys, sltxt, &idx, nw, filename, nline);
                    if (!success) error = -1;
                }
            }
            fclose(finfo);
        } else {
            fprintf(sys->log(), "%s%s : error : info file \"%s\" not found%s\n", sys->is_log_tty?color_string_of_error:"", software_name, filename, sys->is_log_tty?color_string_end:""); success = false; error = 1;
        }
    }
    return error;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int analysis_params(IET_Param * sys, int argc, char * argv[]){
    bool success = true; int error = 0;
    strcpy(szfn_log, "stdout");
  // analysis command line params: pre compiling only handles with -help and -p/-s
    __iaanow = 0;
    if (argc<2){ printf("%s %s\n", software_name, software_version); return 1; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '#' || argv[i][0] == ';'){
        } else {
            int _error = analysis_parameter_line_pre(sys, argv, &i, argc, nullptr, i); if (_error){ success = false; error |= _error; }
        }
        if (error >= 2 && error < 20480){
            if (error==2){
                printf("%s\n", software_version);
            } else {
                printf("%s %s %s\n", software_name, software_version, copyright_string);
            }
            if (error & 4){ // basic help information
                printf("%s%s%s%s%s%s", szHelp1, szHelpXTC, szHelpMP, szHelp2, szHelpLibZ, szHelp3);
            }
            if (error & 8){ // more help information
                printf("%s", szHelpCommands);
                printf("%s%s%s%s", szHelpSecRISM3D_RISM, szHelpSecRISM3D_HI, szHelpSecATOM, szHelpSecRISM3DS);
            }
            if (error & 16){
                printf("%s%s%s", szHelpAdvanced, szHelpInteractive, szHelpInternal);
              #ifdef _EXPERIMENTAL_
                printf("%s", szHelpExperimental);
              #endif
                printf("%s", szHelpNote);
            }
            if ((error & 512) && help_search_str){
                const char * helps[] = { szHelp1, szHelpXTC, szHelpMP, szHelp2, szHelpLibZ, szHelp3, szHelpCommands, szHelpAdvanced, szHelpInteractive, szHelpInternal,
                  #ifdef _EXPERIMENTAL_
                    szHelpExperimental,
                  #endif
                    szHelpSecRISM3D_RISM, szHelpSecRISM3D_HI, szHelpSecRISM3DS, szHelpSecATOM, szHelpNote };
                int line_num = 1; bool is_stdout_tty = isatty(fileno(stdout));
                bool found = false; StringNS::string hkey = help_search_str;
                for (int ih=0; ih<sizeof(helps)/sizeof(helps[0]); ih++){
                    int im = 0; int last_im = 0; StringNS::string templ = helps[ih];
                    while (im<templ.length-hkey.length){
                        if (StringNS::string(&templ.text[im], hkey.length) == hkey){
                            found = true;
                            int begin = im; int end = im;
                            while (begin>0 && helps[ih][begin-1]!='\n') begin --;
                            while (end+1<templ.length && helps[ih][end+1]!='\r' && helps[ih][end+1]!='\n') end ++;
                            //while (end<templ.length){ end ++; if (helps[ih][end]=='\r' || helps[ih][end]=='\n'){ end--; break; } }
                            printf("%s%3d%s ", is_stdout_tty?prompt_comment_prefix:"", line_num, is_stdout_tty?prompt_comment_suffix:"");
                            if (is_stdout_tty){
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

                        if (im>=templ.length-hkey.length) im = templ.length;
                        for (int i=last_im; i<im; i++) if (helps[ih][i]=='\n') line_num ++;
                        last_im = im;
                    }
                }
                if (!found) printf("%s : no help entry for %s\n", software_name, help_search_str);
            }
            return 1;
        }
    }
    if (!success) return error;
  // analysis script params
    __iaanow = 0;
    error = analysis_params_file(sys, info_file_name);
    if (error) return error;
  // analysis command line again to override the parameters
    __iaanow = 0;
    if (argc<2){ printf("%s %s\n", software_name, software_version); return 1; }
    for (int i=1; i<argc; i++){
        int _error = analysis_parameter_line(sys, argv, &i, argc, nullptr, i); if (_error){
            success = false; error |= _error;
        }
    }
    return error;
}
int analysis_params_post(IET_Param * sys){
    int error = 0;
  // post handling
  // above all: check whether solvent file is defined
    if (!info_file_name[0]){
        fprintf(sys->log(), "%s%s : error : solvent file (-p or -solvent) not defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:"");
        error = 1;
    }
  // basic: parameter auto filling based on settings
    sys->nv = sys->nav;
    if (sys->gamma_auto_generate) sys->gamma_erf = 2 / (sys->rcoul<=0? 1e-12 : sys->rcoul);
    if (sys->ndiis_hi>MAX_DIIS) sys->ndiis_hi = MAX_DIIS;
    if (sys->ndiis_rism>MAX_DIIS) sys->ndiis_rism = MAX_DIIS;
  // output filename extension
    if (szfn_out[0]){
        StringNS::string ext = file_extension(szfn_out);
        if (ext!="ts4s" && ext!="ts32" && ext!="ts64") strncat(szfn_out, ".ts4s", MAX_PATH-strnlen(szfn_out, MAX_PATH));
    }
  // forcefield prefix generating (forcefield prefix only for output)
    if (fabs(sys->default_temperature-502.97)<=1 && !sys->mixingrule_sigma_geometric){
        sys->forcefield_prefix = FFPREFIX_AMBER;
    } else if (fabs(sys->default_temperature-120.27)<=1 && sys->mixingrule_sigma_geometric){
        sys->forcefield_prefix = FFPREFIX_GAFF;
    } else if (fabs(sys->default_temperature-120.27)<=1 && !sys->mixingrule_sigma_geometric){
        sys->forcefield_prefix = FFPREFIX_OPLS;
    } else sys->forcefield_prefix = 0;
  // molecular atom mapping
    sys->nvm = 0; for (int i=0; i<sys->nv; i++) if (sys->nvm<=sys->av[i].iaa) sys->nvm = sys->av[i].iaa+1;
  // RISM dielect rebuilding if necessary
    if (!error && (sys->dielect_from_dipole || sys->dipole_from_dielect)){
        if (sys->dielect_from_dipole && sys->dipole_from_dielect){
            fprintf(sys->log(), "%s%s : error : cannot do both dielect-from-dipole and dipole-from-dielect%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); error = 1;
        } else if (sys->dielect_from_dipole){
            double beta = sys->default_temperature / sys->temperature;
            for (int ivm=0; ivm<sys->nvm; ivm++){
                sys->dielect_mol[ivm] = 1 + 4*PI/3*COULCOOEF * beta * sys->density_mv[ivm] * sys->dipole_mv[ivm]*sys->dipole_mv[ivm];
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: dielect[%d] set to %g\n", ivm+1, sys->dielect_mol[ivm]);
            }
        } else if (sys->dipole_from_dielect){
            double beta = sys->default_temperature / sys->temperature;
            for (int ivm=0; ivm<sys->nvm; ivm++){
                sys->dipole_mv[ivm] = sqrt(fabs((sys->dielect_mol[ivm]-1) / (4*PI/3*COULCOOEF * beta * sys->density_mv[ivm])));
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: dipole[%d] set to %g\n", ivm+1, sys->dipole_mv[ivm]);
            }
        }
    }
  // add density/dielectric to atoms
    for (int i=0; i<sys->nv; i++) sys->density_av[i] = sys->density_mv[sys->av[i].iaa];
    for (int i=0; i<sys->nv; i++) sys->dielect[i] = sys->dielect_mol[sys->av[i].iaa];
//printf("number of mols: %d\n", sys->nvm);
    int nvmmap = 0;
    for (int i=0; i<sys->nvm; i++){
        sys->vmmapi[i][0] = sys->vmmapi[i][1] = nvmmap;
        for (int j=0; j<sys->nv; j++) if (sys->av[j].iaa == i){ sys->vmmap[nvmmap] = j; sys->vmmapi[i][1] = ++nvmmap; }
    }
//for (int i=0; i<sys->nvm; i++){ printf("mol %d : %d atoms : density %12f\n", i, sys->vmmapi[i][1] - sys->vmmapi[i][0], sys->density_mv[i]); for (int j=sys->vmmapi[i][0]; j<sys->vmmapi[i][1]; j++) printf("  %d:%s:%s", j,sys->av[sys->vmmap[j]].mole,sys->av[sys->vmmap[j]].name); printf("\n"); }
    if (!error){
        if (sys->stepmax_rism<0) sys->stepmax_rism = 0; else if (sys->stepmax_rism>1000000000) sys->stepmax_rism = 1000000000;
    }
  // check sigma
    if (!error){
        for (int iv=0; iv<sys->nv; iv++){
            if (fabs(sys->av[iv].sigma)<=MACHINE_REASONABLE_ERROR || fabs(sys->av[iv].epsilon)<=MACHINE_REASONABLE_ERROR){
                for (int ib=0; ib<sys->n_bond_list; ib++) if (sys->bond_list[ib].grpi==iv || sys->bond_list[ib].grpj==iv){
                }
            }
        }
    }
  // other check
    if (!error) if (!sys->mixingrule_sigma_geometric_specified){
        fprintf(sys->log(), "%s%s : error : -ff : force field type not defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); error = 1;
    }
    //if (!error) if (sys->gvv_specification==0){
    //    fprintf(sys->log(), "%s : error : solvent correlation (-gvv) not specified.\n", software_name); error = 1;
    //}
    if (!error) if (sys->nav<1){
        fprintf(sys->log(), "%s%s : error : no solvent defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); error = 1;
    }
    if (!error) if (sys->nmv < sys->nvm){
        if (sys->nmv<1){
            fprintf(sys->log(), "%s%s : error : -density must be defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); error = 1;
        } else {
            fprintf(sys->log(), "%s%s : error : -density : %d item%s defined, but %d required%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nmv, sys->nmv>1?"s":"", sys->nvm, sys->is_log_tty?color_string_end:""); error = 1;
        }
    }
    if (!error) if (sys->nmvb < sys->nmv){
        if (sys->nmv<=1){
            sys->nmvb = sys->nmv; sys->bulk_density_mv[0] = sys->density_mv[0];
        } else if (sys->nmvb<1){
            fprintf(sys->log(), "%s%s : error : -bulk_density must be defined for mixture%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); error = 1;
        } else {
            fprintf(sys->log(), "%s%s : error : -bulk_density : %d item%s defined, but %d required%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nmvb, sys->nmvb>1?"s":"", sys->nmv, sys->is_log_tty?color_string_end:""); error = 1;
        }
    }
    if (!error) if (sys->n_dielect_mol < sys->nvm){
        fprintf(sys->log(), "%s%s : error : -dielectric : %d item%s defined, but %d required%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->n_dielect_mol, sys->n_dielect_mol>1?"s":"", sys->nvm, sys->is_log_tty?color_string_end:""); error = 1;
    }
    if (!error){
        if (sys->nmvam <= 0){
            /*for (int iv=0; iv<sys->nvm; iv++){
                int na_in_m = 0; const char * mole = ""; for (int i=0; i<sys->nv; i++) if (sys->av[i].iaa==iv && sys->av[i].mass>=4){ na_in_m += sys->av[i].multi; mole = sys->av[i].mole;
                }
                sys->density_hi[iv] = na_in_m * sys->density_mv[iv];
                if (mole[0]) fprintf(sys->log(), "%s : density_hi[%s] = %g <-- %d x %g\n", software_name, mole, sys->density_hi[iv], na_in_m, sys->density_mv[iv]); else fprintf(sys->log(), "%s : density_hi[%d] = %g <-- %d x %g\n", software_name, iv+1, sys->density_hi[iv], na_in_m, sys->density_mv[iv]);
            }*/
            sys->nmvam = sys->nvm;
        }
        if (sys->nmvam == 1){
            sys->nmvam = sys->nvm;
        } else if (sys->nmvam < sys->nvm){
            fprintf(sys->log(), "%s%s : error : too few atomic density (-density-hi) : (%d vs %d)%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nmvam, sys->nvm, sys->is_log_tty?color_string_end:""); error = 1;
        }
    }
    if (!error){
        double total_nbulk = 0;
        for (int i=0; i<sys->nvm; i++){
            sys->nbulk[i] = sys->density_mv[i] / sys->bulk_density_mv[i];
            sys->nbulk_rism[i] = sys->enable_nbulk_rism? sys->nbulk[i] : 1;
            total_nbulk += sys->nbulk[i];
        }
        for (int i=0; i<sys->nvm; i++) sys->nbulk[i] /= total_nbulk;
        total_nbulk = 1;

        if (sys->debug_level>=1){
            fprintf(sys->log(), "%s : nbulk autogenerated :", software_name); for (int i=0; i<sys->nvm; i++) fprintf(sys->log(), " %g%%", sys->nbulk[i]*100); fprintf(sys->log(), "\n");
        }
    }
    if (!error){
        if (sys->nllambda==0){
            for (int i=0; i<sys->nvm; i++) sys->llambda[i] = 0;
        } else if (sys->nllambda < sys->nvm){
            fprintf(sys->log(), "%s%s : error : too few llambda (%d vs %d)%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nllambda, sys->nvm, sys->is_log_tty?color_string_end:""); error = 1;
        }
    }
    if (!error){
      // serial mixing rule of dielectric constance
        double dielecti = 0; double density_dielect = 0;
        for (int ivm=0; ivm<sys->nvm; ivm++){
            dielecti += sys->density_mv[ivm] / sys->dielect_mol[ivm]; density_dielect += sys->density_mv[ivm];
        }
        sys->mean_dielect = density_dielect / dielecti;
        //if (sys->esal==CoulAL_Coulomb) sys->mean_dielect = 1;
    }
    if (!error){
      // calculate Debye wave number
        /*if (sys->n_debye_rate < sys->nvm){
            fprintf(sys->log(), "%s : error : %d Debye decay rates required, but only %d given\n", software_name, sys->n_debye_rate, sys->nvm); return false;
        }*/
        if (sys->n_debye_rate > sys->nvm) sys->n_debye_rate = sys->nvm;
        for (int i=0; i<MAX_DEBYE_TERMS; i++) sys->debye_kappa[i] = 0;
        for (int ivm=0; ivm<sys->n_debye_rate; ivm++)/* if (sys->debye_rate[ivm] <= 0)*/{
            const char * mole = ""; for (int i=0; i<sys->nv; i++) if (sys->av[i].iaa==ivm) mole = sys->av[i].mole;
            double beta = sys->default_temperature / sys->temperature;
            double qq = 0; double qqq = 0; double qqqq = 0; double q5 = 0; double q6 = 0;
            for (int ia=0; ia<sys->nv; ia++) if (sys->av[ia].iaa == ivm){
                qq   += sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].multi;
                qqq  += sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].multi;
                qqqq += sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].multi;
                q5 += sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].multi;
                q6 += sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].charge * sys->av[ia].multi;
            }
            double kappa2 = 4*PI* COULCOOEF * beta * sys->density_mv[ivm] / sys->mean_dielect * qq;
            sys->debye_rate[ivm] = sqrt(fabs(kappa2));
            sys->debye_kappa[1] += kappa2;
            sys->debye_kappa[2] += 4*PI* COULCOOEF * beta * beta * sys->density_mv[ivm] / sys->mean_dielect * qqq / 2;
            sys->debye_kappa[3] += 4*PI* COULCOOEF * beta*beta*beta * sys->density_mv[ivm] / sys->mean_dielect * qqqq / 6;
            sys->debye_kappa[4] += 4*PI* COULCOOEF * beta*beta*beta*beta * sys->density_mv[ivm] / sys->mean_dielect * q5 / 24;
            sys->debye_kappa[5] += 4*PI* COULCOOEF * beta*beta*beta*beta*beta * sys->density_mv[ivm] / sys->mean_dielect * q6 / 120;


            char mol_number[64]; snprintf(mol_number, sizeof(mol_number), "#%d", ivm+1);
            if (!mole[0]) mole = mol_number;

            if (sys->esal==CoulAL_YukawaFFT && sys->debug_level>=1) fprintf(sys->log(), "debug:: Debye_k[%s] = %.4g <-- dielectric %.4g, q^2 %g\n", mole, sys->debye_rate[ivm], sys->dielect_mol[ivm], qq);


            //if (mole[0]) fprintf(sys->log(), "%s : debye_rate[%s] = %g <-- dielectric %g\n", software_name, mole, sys->debye_rate[ivm], sys->dielect_mol[ivm]); else fprintf(sys->log(), "%s : debye_rate[%d] = %g <-- dielectric %g\n", software_name, ivm+1, sys->debye_rate[ivm], sys->dielect_mol[ivm]);
//fprintf(sys->log(), "debye_rate[%d] = %f\n", ivm, sys->debye_rate[ivm]);
        }
        double debye_kappa_1 = sys->debye_kappa[1]; for (int i=0; i<MAX_DEBYE_TERMS; i++) sys->debye_kappa[i] = debye_kappa_1<MACHINE_REASONABLE_ERROR? 0 : sys->debye_kappa[i]/debye_kappa_1;
        //if (sys->debug_level>=1/* && sys->esal==CoulAL_Yukawa*/) fprintf(sys->log(), "%s : system debye expansion: %g %g %g %g %g\n", software_name, sys->debye_kappa[1], sys->debye_kappa[2], sys->debye_kappa[3], sys->debye_kappa[4], sys->debye_kappa[5]);

        if (sys->esal==CoulAL_YukawaFFT){
        //if (sys->esal==CoulAL_Yukawa || sys->esal==CoulAL_DH){
            double yukawa_alpha = 0; for (int i=0; i<sys->nvm; i++) yukawa_alpha += sys->debye_rate[i] * sys->debye_rate[i] * sys->nbulk[i]; yukawa_alpha = sqrt(fabs(yukawa_alpha));
            sys->yukawa_alpha = yukawa_alpha;
            //if (yukawa_alpha * sys->rcoul < 2){
            //    fprintf(sys->log(), "%s : warning : -rcoul too short, use %g (suggested by Debye-rate %g) or use LPBI instead\n", software_name, 10/yukawa_alpha, yukawa_alpha);
                //error ++;
            //} else if (yukawa_alpha * sys->rcoul < 5){
            //    fprintf(sys->log(), "%s : warning : -rcoul too short, suggest %g (by Debye-rate %g)\n", software_name, 10/yukawa_alpha, yukawa_alpha);
            //}
        }
      // check Debye wave length of YukawaFFT
        if (sys->rc_yukawafft<0) sys->rc_yukawafft = 1 / sys->yukawa_alpha;
    }
    if (sys->cmd_flag_rdf_ever_display && sys->n_rdf_grps<=0){
        fprintf(sys->log(), "%s%s : warning : rdf will not show as -rdf-grps is missing%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->is_log_tty?color_string_end:"");
    }
//printf("density: "); for (int i=0; i<sys->nmv; i++) printf("%11f ", sys->density_mv[i]); printf("\n");
//printf("density_hi: "); for (int i=0; i<sys->nmv; i++) printf("%11f ", sys->density_hi[i]); printf("\n");
//printf("%d atoms: corre dr: %12f, grids: %dx%dx%d\n", sys->nav, sys->dr, sys->nr[0], sys->nr[1], sys->nr[2]); for (int i=0; i<sys->nv; i++) printf("    %9s(%d).%-6s : %12f %12f %12f %12f\n", sys->av[i].mole, sys->av[i].iaa, sys->av[i].name, sys->bulk_density_av[i], sys->av[i].charge, sys->av[i].sigma, sys->av[i].epsilon);
    return error;
}
int analysis_params_post_commands(IET_Param * sys){
    int error = 0; if (!sys->cmd || sys->ncmd<=0) return 0;

  // check command queue
    for (int ic=0; ic<sys->ncmd; ic++){
        /*if (sys->cmd[ic].command>1500&&sys->cmd[ic].command<2000){
            sys->cmd_flag_rism_ever_performed = true;
        } else if (sys->cmd[ic].command==IETCMD_DISPLAY || sys->ncmd==IETCMD_REPORT){
            for (int i=0; i<sys->cmd[ic].step; i++){
                if (sys->cmd[ic].command_params_int[i] == IETCMD_v_rdf){
                    sys->cmd_flag_rdf_ever_display = true;
                } else {
                    sys->cmd_flag_energy_ever_display = true;
                }
            }
        } else if (sys->cmd[ic].command==IETCMD_TI){
            sys->cmd_flag_energy_ever_display = true;
        } else */
        if (sys->cmd[ic].command==IETCMD_CLOSURE || sys->ncmd==IETCMD_CLOSURE_A){
            for (int i=0; i<sys->cmd[ic].step; i++){
                if (sys->cmd[ic].command_params_int[i] == CLOSURE_RBC_HNC){
                    sys->cmd_rbc_ljr_allowed = true;
                } else if (sys->cmd[ic].command_params_int[i] == CLOSURE_RBC_KH){
                    sys->cmd_rbc_ljr_allowed = true;
                }
            }
        }
    }

  // check global settings
    for (int i=0; i<sys->nv; i++){
        if (sys->closures[i]==CLOSURE_RBC_HNC){
            sys->cmd_rbc_ljr_allowed = true;
        } else if (sys->closures[i]==CLOSURE_RBC_KH){
            sys->cmd_rbc_ljr_allowed = true;
        }
    }

    return error;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



int check_hi_diis_default_steps(IET_Param * sys){ if (sys->ndiis_hi<=0) sys->ndiis_hi = 5; return sys->ndiis_hi; }
int check_rism_diis_default_steps(IET_Param * sys){ if (sys->ndiis_rism<=0) sys->ndiis_rism = 5; return sys->ndiis_rism; }

bool analysis_command(IET_Param * sys, char * line, const char * line_orig, const char * script_name, int script_line, int first_char_offset){
    bool success = true; int cmd_type = -1;
    IET_command cmd; memset(&cmd, 0, sizeof(IET_command)); int cmd_location = sys->ncmd + 1;
    cmd_type = 0; cmd.command = 0;
    int last_cmd_type = cmd_type; int last_param_num = cmd.step;

    StringNS::string sl[200];
    int nw = analysis_word_line(line, sl, 200);

    int i_param_list = 0;

  // analysis commands
    if (nw>0){
        if (sl[0]=="nop" || sl[0]=="noop"){
            cmd_type = 0; cmd.command = IETCMD_NOP; cmd.step = 1;
        } else if (sl[0]=="norism"){
            cmd_type = 0; cmd.command = IETCMD_NOP; cmd.step = 1;
      // cmd: I/O operations
        } else if (sl[0]=="clear"){ cmd_type = 0; cmd.command = IETCMD_CLEAR; cmd.step = 1;
        } else if (sl[0]=="reset"){ cmd_type = 0; cmd.command = IETCMD_RESET; cmd.step = 1;
        } else if (sl[0]=="end" || sl[0]=="stop" || sl[0]=="exit"){
            cmd_type = 0; cmd.command = IETCMD_END; cmd.step = 1;
        } else if (sl[0]=="done" || sl[0]=="break" || sl[0]=="return"){
            cmd_type = 0; cmd.command = IETCMD_DONE; cmd.step = 1;
        } else if (sl[0]=="set"){       cmd_type = 1; cmd.command = IETCMD_SET;
        } else if (sl[0]=="let"){       cmd_type = 1; cmd.command = IETCMD_SET;
        } else if (sl[0]=="scale"){     cmd_type = 2; cmd.command = IETCMD_SCALE;
        } else if (sl[0]=="load"){
            cmd_type = 5; cmd.command = IETCMD_LOAD; cmd.time_to_run = 0;
            //fprintf(sys->log(), "%s : %s[%d][%ld] : warning : \"%s\" redirected to \"load-at-the-beginning\"\n", software_name, script_name, script_line, 1+sl[0].text-line+first_char_offset, sl[0].text);
        } else if (sl[0]=="save-filter" || sl[0]=="save_filter" || sl[0]=="save-site" || sl[0]=="save-sites" || sl[0]=="save_site" || sl[0]=="save_sites" || sl[0]=="pick" || sl[0]=="pick-site" || sl[0]=="pick-sites" || sl[0]=="pick_site" || sl[0]=="pick_sites"){
            cmd_type = 61; cmd.command = IETCMD_SAVE_FILTER; cmd.time_to_run = 0;
        } else if (sl[0]=="save"){
            cmd_type = 6; cmd.command = IETCMD_SAVE; cmd.time_to_run = 0;
            //fprintf(sys->log(), "%s : %s[%d][%ld] : warning : \"%s\" redirected to \"save-at-the-end\"\n", software_name, script_name, script_line, 1+sl[0].text-line+first_char_offset, sl[0].text);
        } else if (sl[0]=="savee" || sl[0]=="save-exist"){
            cmd_type = 6; cmd.command = IETCMD_SAVE_EXIST; cmd.time_to_run = 0;
        } else if (sl[0]=="calculate"){ cmd_type = 7; cmd.command = IETCMD_DISPLAY;
        } else if (sl[0]=="calc"){      cmd_type = 7; cmd.command = IETCMD_DISPLAY;
        } else if (sl[0]=="display"){   cmd_type = 7; cmd.command = IETCMD_DISPLAY;
        } else if (sl[0]=="print"){     cmd_type = 7; cmd.command = IETCMD_DISPLAY;
        } else if (sl[0]=="report"){    cmd_type = 7; cmd.command = IETCMD_REPORT;
        } else if (sl[0]=="hi-solver"){ cmd_type = 0; cmd.command = IETCMD_HI_SOLVER;
        } else if (sl[0]=="lse"){       cmd_type = 0; cmd.command = IETCMD_LSE;
      // cmd: kernel calculation: HI, IET
        } else if (sl[0]=="build-force-field" || sl[0]=="build_force_field" || sl[0]=="build-ff" || sl[0]=="build_ff" || sl[0]=="rebuild-force-field" || sl[0]=="rebuild_force_field" || sl[0]=="rebuild-ff" || sl[0]=="rebuild_ff"){
            cmd_type = 13; cmd.command = IETCMD_BUILD_FF; cmd.step = 1;
        } else if (sl[0]=="no-build-force-field" || sl[0]=="no_build_force_field" || sl[0]=="no-build-ff" || sl[0]=="no_build_ff" || sl[0]=="skip-force-field" || sl[0]=="skip_force_field" || sl[0]=="skip-ff" || sl[0]=="skip_ff"){
            cmd_type = 13; cmd.command = -IETCMD_BUILD_FF; cmd.step = 1;
        } else if (sl[0]=="build-uuv" || sl[0]=="build_uuv" || sl[0]=="rebuild-uuv" || sl[0]=="rebuild_uuv"){
            cmd_type = 13; cmd.command = IETCMD_BUILD_UUV; cmd.step = 1;
            cmd.command_params_int[0] = CoulAL_Coulomb;
            cmd.command_params_double[0] = 1;
            cmd.command_params_double[1] = 0;
        } else if (sl[0]=="hi"){        cmd_type = 10; cmd.command = 2000+HIAL_HSHI;      cmd.step = sys->stepmax_hi;
            check_hi_diis_default_steps(sys);
            //fprintf(sys->log(), "%s : %s[%d][%ld] : warning : \"%s\" redirected to \"hshi\"\n", software_name, script_name, script_line, 1+sl[0].text-line+first_char_offset, sl[0].text);
        } else if (sl[0]=="nohi"){      cmd_type = 10; cmd.command = 2000+HIAL_NOHI;      cmd.step = sys->stepmax_hi;
        } else if (sl[0]=="hshi"){      cmd_type = 10; cmd.command = 2000+HIAL_HSHI;      cmd.step = sys->stepmax_hi;
            check_hi_diis_default_steps(sys);
        //} else if (sl[0]=="h2hi"){      cmd_type = 10; cmd.command = 2000+HIAL_H2HI;      cmd.step = sys->stepmax_hi;
        //    check_hi_diis_default_steps(sys);
        //} else if (sl[0]=="eehi"){      cmd_type = 10; cmd.command = 2000+HIAL_EEHI;      cmd.step = sys->stepmax_hi;
        //    check_hi_diis_default_steps(sys);
        //    sys->b_allow_Ecoul0 = true;
        //} else if (sl[0]=="e2hi"){      cmd_type = 10; cmd.command = 2000+HIAL_E2HI;      cmd.step = sys->stepmax_hi;
        //    check_hi_diis_default_steps(sys);
        //    sys->b_allow_Ecoul0 = true;
        //} else if (sl[0]=="dphi"){      cmd_type = 10; cmd.command = 2000+HIAL_DPHI;      cmd.step = sys->stepmax_hi;
        //    check_hi_diis_default_steps(sys);
        //    sys->b_allow_Ecoul0 = true;
        //} else if (sl[0]=="dnhi"){      cmd_type = 10; cmd.command = 2000+HIAL_DNHI;      cmd.step = sys->stepmax_hi;
        //    check_hi_diis_default_steps(sys);
        //    sys->b_allow_Ecoul0 = true;
        } else if (sl[0]=="init-rism" || sl[0]=="norism"){
            cmd_type = 12; cmd.command = 1500+IETAL_NONE;
        } else if (sl[0]=="ssoz"){      cmd_type = 12; cmd.command = 1500+IETAL_SSOZ;   cmd.step = sys->stepmax_rism;
            sys->rism_coulomb_renormalization = false;
            check_rism_diis_default_steps(sys);
            sys->cmd_flag_rism_ever_performed = true;
        } else if (sl[0]=="rism"){      cmd_type = 12; cmd.command = 1500+IETAL_RRISM;  cmd.step = sys->stepmax_rism;
            sys->rism_coulomb_renormalization = true;
            check_rism_diis_default_steps(sys);
            sys->cmd_flag_rism_ever_performed = true;
        } else if (sl[0]=="ti" || sl[0]=="tibegin" || sl[0]=="ti_begin" || sl[0]=="ti-begin"){
            cmd_type = 14; cmd.command = IETCMD_TI; cmd.step = 10;
            sys->b_allow_Ecoul0 = true;
        } else if (sl[0]=="hold"){
            cmd_type = IETCMD_HOLD; cmd.command = IETCMD_HOLD; cmd.step = 0;
      // cmd: set values of arrays
        } else if (sl[0]=="closure"){   cmd_type = 21; cmd.command = IETCMD_CLOSURE;
        } else if (sl[0]=="closure-a" || sl[0]=="closure_a"){ cmd_type = 21; cmd.command = IETCMD_CLOSURE_A;
        } else if (sl[0]=="closure-m" || sl[0]=="closure_m"){ cmd_type = 21; cmd.command = IETCMD_CLOSURE;
        } else if (sl[0]=="cf"){        cmd_type = 22; cmd.command = IETCMD_CF;
        } else if (sl[0]=="closure-factor" || sl[0]=="closure_factor"){ cmd_type = 22; cmd.command = IETCMD_CF;
        } else if (sl[0]=="closure-factor-a" || sl[0]=="closure_factor_a"){ cmd_type = 22; cmd.command = IETCMD_CF_A;
        } else if (sl[0]=="cf-a" || sl[0]=="cf_a"){ cmd_type = 22; cmd.command = IETCMD_CF_A;
        } else if (sl[0]=="closure-factor-m" || sl[0]=="closure_factor_m"){ cmd_type = 22; cmd.command = IETCMD_CF;
        } else if (sl[0]=="cf-m" || sl[0]=="cf_m"){ cmd_type = 22; cmd.command = IETCMD_CF;
        } else if (sl[0]=="density"){   cmd_type = 23; cmd.command = IETCMD_DENSITY;
        } else if (sl[0]=="dielect"){   cmd_type = 24; cmd.command = IETCMD_DIELECT;
        } else if (sl[0]=="test"){ cmd_type = IETCMD_TEST; cmd.command = IETCMD_TEST;
        } else if (sl[0]=="test-and-save"){ cmd_type = IETCMD_TEST_SAVE; cmd.command = IETCMD_TEST_SAVE;
        } else if (sl[0]=="rdf-content"){ cmd_type = IETCMD_RDF_CONTENT; cmd.command = IETCMD_RDF_CONTENT;
        } else if (sl[0]=="temperature"){ cmd_type = IETCMD_TEMPERATURE; cmd.command = IETCMD_TEMPERATURE;
      // cmd: experimental
      #ifdef _EXPERIMENTAL_
        } else if (experimental_analysis_command(sl[0], i_param_list, cmd, cmd_type, sys->log(), sys->is_log_tty, script_name, script_line)){
      #endif
      // cmd: done
        } else {
            char buffer[512]; memset(buffer, 0, sizeof(buffer)); memcpy(buffer, sl[0].text, sl[0].length>sizeof(buffer)-1?sizeof(buffer)-1:sl[0].length);
            fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined command \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[0].text-line+first_char_offset, buffer, sys->is_log_tty?color_string_end:"");
            success = false;
        }
    }


  // analysis parameters of commands
    if (success && nw>1){
        //printf(" (%d, type %d) parameter list:", cmd.command, cmd_type); for (int i=1; i<nw; i++){ char buffer[512]; memset(buffer, 0, sizeof(buffer)); memcpy(buffer, sl[i].text, sl[i].length>sizeof(buffer)-1?sizeof(buffer)-1:sl[i].length); printf(" [%s]", buffer); } printf ("\n");
        for (int i=1; i<nw; i++) if (sl[i].length>0){
            if (sl[i].text[-1]=='@'){
                StringNS::string at_location; at_location.text = &sl[i].text[0]; at_location.length = sl[i].length - 0;
                if (at_location.length>0){
                    if (StringNS::is_string_number(at_location)){
                        cmd_location = atoi(&sl[i].text[1]);
                    } else if (at_location=="b" || at_location=="begin" || at_location=="beginning"){
                        cmd.time_to_run = 1;
                    } else if (at_location=="e" || at_location=="end" || at_location=="ending" || at_location=="final"){
                        cmd.time_to_run = -1;
                    } else {
                        char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                        fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : unrecognizable location: \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                        success = false; sl[i].text[sl[i].length] = cc;
                    }
                }
            } else if (cmd_type==10){ // HI and CI
                if ((sl[i]=="step" || sl[i]=="steps")){
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        cmd.step = atoi(sl[i+1].text);
                        i++;
                    }
                } else if (sl[i]=="factor"){ // format: hi,factor=1,2,3...
                    while (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        double value = atof(sl[i+1].text); i++;
                        int i_available=-1; for (int ic=0; ic<MAX_CMD_PARAMS; ic++) if (cmd.command_params_int[ic]==0){ i_available = ic; break; }
                        if (i_available>=0 && i_available<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_available] = 1;
                            cmd.command_params_double[i_available] = value;
                        }
                    }
                } else if (sl[i]=="te" || sl[i]=="theta-expand" || sl[i]=="theta_expand" || sl[i]=="theta-expansion" || sl[i]=="theta_expansion"){
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        double value = atof(sl[i+1].text); i++;
                        int i_available=-1; for (int ic=0; ic<MAX_CMD_PARAMS; ic++) if (cmd.command_params_int[ic]==0){ i_available = ic; break; }
                        if (i_available>=0 && i_available<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_available] = 0xB34349; // last six digitals of md5sum("theta expansion")
                            cmd.command_params_double[i_available] = value;
                        }
                    }
                } else if (sl[i]=="ke" || sl[i]=="kernel-expand" || sl[i]=="kernel_expand" || sl[i]=="kernel-expansion" || sl[i]=="kernel_expansion" || sl[i]=="ld" || sl[i]=="ld-expand" || sl[i]=="ld_expand" || sl[i]=="ld-expansion" || sl[i]=="ld_expansion"){
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        double value = atof(sl[i+1].text); i++;
                        int i_available=-1; for (int ic=0; ic<MAX_CMD_PARAMS; ic++) if (cmd.command_params_int[ic]==0){ i_available = ic; break; }
                        if (i_available>=0 && i_available<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_available] = 0xB59E5B; // last six digitals of md5sum("kernel expansion")
                            cmd.command_params_double[i_available] = value;
                        }
                    }
                } else if (sl[i]=="dielect" || sl[i]=="dielect-hi" || sl[i]=="dielect_hi"){
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        double value = atof(sl[i+1].text); i++;
                        int i_available=-1; for (int ic=0; ic<MAX_CMD_PARAMS; ic++) if (cmd.command_params_int[ic]==0){ i_available = ic; break; }
                        if (i_available>=0 && i_available<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_available] = 0x724E6D; // last six digitals of md5sum("dielect hi")
                            cmd.command_params_double[i_available] = value;
                        }
                    }
                } else if (StringNS::is_string_number(sl[i])){ // format: hi,1,2,3...
                    double value =  atof(sl[i].text);
                    int i_available=-1; for (int ic=0; ic<MAX_CMD_PARAMS; ic++) if (cmd.command_params_int[ic]==0){ i_available = ic; break; }
                    if (i_available>=0 && i_available<MAX_CMD_PARAMS){
                        cmd.command_params_int[i_available] = 1;
                        cmd.command_params_double[i_available] = value;
                    }
                /*} else if (sl[i]=="factor"||sl[i]=="ih" || sl[i]=="is" || sl[i]=="oh" || sl[i]=="os"){
                    memset(cmd.command_params_int, 0, sizeof(cmd.command_params_int));
                    int factor_identifier = 1;
                        if (sl[i]=="oh") factor_identifier = -1;
                        else if (sl[i]=="is") factor_identifier = 2;
                        else if (sl[i]=="os") factor_identifier = -2;
                    int icpi = 0; while (icpi<MAX_CMD_PARAMS && i+1<nw && StringNS::is_string_number(sl[i+1])){
                        cmd.command_params_int[icpi] = factor_identifier;
                        cmd.command_params_double[icpi] = atof(sl[i+1].text);
                        i++; icpi ++;
                    }*/
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined variable \"%s\" for %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, HIAL_name[cmd.command-2000], sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
            } else if (cmd_type==12){ // RISM
                if (sl[i]=="step" || sl[i]=="steps"){
                    if (i+1<nw&&StringNS::is_string_number(sl[i+1])){ if (StringNS::is_string_number(sl[i+1])) cmd.step = atoi(sl[i+1].text); i++; }
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined variable \"%s\" for RISM%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
            } else if (cmd_type==14){ // TI
                if (sl[i]=="step" || sl[i]=="steps"){
                    if (i+1<nw&&StringNS::is_string_number(sl[i+1])){ if (StringNS::is_string_number(sl[i+1])) cmd.step = atoi(sl[i+1].text); i++; }
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined variable \"%s\" for TI%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                // parameter line in general format
            } else if (cmd.command==IETCMD_SET){
                if (sl[i]=="rism-dielect"){     if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_rism_dielect;
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined variable \"%s\" to set%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_LOAD){
                if (sl[i]=="lj"){               if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                } else if (sl[i]=="coul"){      if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                    if (i_param_list+2<MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul2;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save liquid structures : error : too many things to load%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="cuv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                } else if (sl[i]=="huv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_huv;
                } else if (sl[i]=="hlr"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_hlr;
                } else if (sl[i]=="dd"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_dd;
                } else if (sl[i]=="nphi"){      if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = -IETCMD_v_dd;
                } else if (sl[i]=="ff"){
                    if (i_param_list+3<MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul2;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save liquid structures : error : too many things to load%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="iet"){
                    if (i_param_list+3<MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_huv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_dd;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save liquid structures : error : too many things to load%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : \"%s\" undefined in loading%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_SAVE_FILTER){    // 0: unintialized; <0: everything
                if (sl[i]=="all"){
                    i_param_list = 0; cmd.command_params_int[i_param_list++] = -1;
                } else if (StringNS::is_string_number(sl[i])){
                    if (i_param_list+1 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = atoi(sl[i].text);
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save-filter : error : too many savings%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined saving term \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_SAVE || cmd.command==IETCMD_SAVE_EXIST){
                if (sl[i]=="nothing"){
                } else if (sl[i]=="all"||sl[i]=="everything"||sl[i]=="anything"||sl[i]=="checkpoint"){
                    if (i_param_list+13 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cmd;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul2;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_rmin; sys->b_allow_r2uvmin = true;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_uuv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulr;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_clr;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_huv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_hlr;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_dd;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_guv;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save all : error : too many things to save%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="potential"||sl[i]=="potentials"||sl[i]=="force-field"||sl[i]=="force_field"||sl[i]=="ff"){
                    if (i_param_list+3 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul2;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save potentials : error : too many things to save%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="iet"||sl[i]=="liquid"||sl[i]=="liquid-structure"||sl[i]=="liquid_structure"||sl[i]=="liquid-struct"||sl[i]=="liquid_struct"){
                    if (i_param_list+4 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_huv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_dd;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_guv;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = save liquid structures : error : too many things to save%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="cmd"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cmd;
                } else if (sl[i]=="command"){   if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cmd;
                } else if (sl[i]=="lj"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                } else if (sl[i]=="coul" || sl[i]=="coulomb"){
                    if (i_param_list < MAX_CMD_PARAMS)  cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                } else if (sl[i]=="electrostatic-field" || sl[i]=="electrostatic_field" || sl[i]=="electric-field" || sl[i]=="electric_field" || sl[i]=="ef"){
                    if (i_param_list < MAX_CMD_PARAMS)  cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul2;
                } else if (sl[i]=="cuv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                } else if (sl[i]=="huv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_huv;
                } else if (sl[i]=="hlr"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_hlr;
                } else if (sl[i]=="dd"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_dd;
                } else if (sl[i]=="ddp"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ddp;
                } else if (sl[i]=="theta"){     if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_theta;
                } else if (sl[i]=="nphi"){      if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = -IETCMD_v_dd;
                } else if (sl[i]=="density"){   if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_guv;
                } else if (sl[i]=="guv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_guv;
                } else if (sl[i]=="rdf"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_rdf;
                    sys->cmd_flag_rdf_ever_display = true;
                } else if (sl[i]=="rmin" || sl[i]=="min-dist" || sl[i]=="min_dist" || sl[i]=="mindist"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_rmin;
                    sys->b_allow_r2uvmin = true;
                } else if (sl[i]=="ld" || sl[i]=="local-density" || sl[i]=="local_density"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ld;
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : \"%s\" undefined in saving%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_SCALE){
                int this_command_params_int = -1;
                char sli_text[20]; strncpy(sli_text, sl[i].text, sizeof(sli_text)); if (sl[i].length<sizeof(sli_text)) sli_text[sl[i].length] = 0;
                if (sl[i]=="lj") this_command_params_int = IETCMD_v_ulj;
                else if (sl[i]=="coul") this_command_params_int = IETCMD_v_ucoul;
                else if (sl[i]=="ff") this_command_params_int = IETCMD_v_uuv;
                else if (sl[i]=="uuv") this_command_params_int = IETCMD_v_uuv;

                if (this_command_params_int>0){
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){
                        if (i_param_list<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_param_list] = this_command_params_int;
                            cmd.command_params_double[i_param_list] = atof(sl[i+1].text);
                            i++;
                            i_param_list ++;
                        }
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d][%ld] : warning : scaling factor of %s not defined.%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sli_text, sys->is_log_tty?color_string_end:"");
                    }
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : \"%s\" undefined for scaling%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if ((cmd.command==IETCMD_REPORT || cmd.command==IETCMD_DISPLAY)){
                if (StringNS::is_string_number(sl[i])){
                    cmd.command_params_double[i_param_list>0?i_param_list-1:i_param_list] = atof(sl[i].text);
                } else if (i_param_list+1<MAX_CMD_PARAMS){
                    if (sl[i]=="all"){
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_mass;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Ef; sys->b_allow_Ecoul0 = true;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_PMV;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_excess_GF;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_excess_RISM;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="energy"){
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_mass;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Ef; sys->b_allow_Ecoul0 = true;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_PMV;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="basic"){
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_mass;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        if (i_param_list+1<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_PMV;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="mass"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_mass;
                    } else if (sl[i]=="lj"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="coul"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="entropy" || sl[i]=="TS"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="N" || sl[i]=="dN"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="N0" || sl[i]=="dN0" || sl[i]=="dN_vac"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="volume" || sl[i]=="V" || sl[i]=="pmv"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_PMV;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="ef0" || sl[i]=="ef" || sl[i]=="Uef0" || sl[i]=="Hef0"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_Ef;
                        sys->cmd_flag_energy_ever_display = true;
                        sys->b_allow_Ecoul0 = true;
                    } else if (sl[i]=="ef1"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_Ef1;
                        sys->cmd_flag_energy_ever_display = true;
                        sys->b_allow_Ecoul0 = true;
                    } else if (sl[i]=="cuv" || sl[i]=="c" || sl[i]=="direct-correlation" || sl[i]=="direct_correlation"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="csr"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_csr;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="clr"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_clr;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="chuv" || sl[i]=="tuv"){
                        cmd.command_params_int[i_param_list++] = -IETCMD_v_cuv;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="chsr"){
                        cmd.command_params_int[i_param_list++] = -IETCMD_v_csr;
                    } else if (sl[i]=="chlr"){
                        cmd.command_params_int[i_param_list++] = -IETCMD_v_clr;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="excessGF" || sl[i]=="exGF" || sl[i]=="excess-GF" || sl[i]=="excess_GF" || sl[i]=="excessChemGF" || sl[i]=="excess-chem-GF" || sl[i]=="excessive-chemical-potential-GF" || sl[i]=="excess_chem_GF" || sl[i]=="excessive_chemical_potential_GF"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_excess_GF;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="chandler" || sl[i]=="excess" || sl[i]=="excessChem" || sl[i]=="excess-chem" || sl[i]=="excessive-chemical-potential" || sl[i]=="excess_chem" || sl[i]=="excessive_chemical_potential"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_excess_RISM;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="excessHyb" || sl[i]=="exHyb" || sl[i]=="excess-Hyb" || sl[i]=="excess_Hyb" || sl[i]=="excessChemHyb" || sl[i]=="excess-chem-hyb" || sl[i]=="excessive-chemical-potential-hybid" || sl[i]=="excess_chem_hyb" || sl[i]=="excessive_chemical_potential_hybrid"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_excess_hyb;
                    } else if (sl[i]=="zeta" || sl[i]=="zeta-hnc" || sl[i]=="zeta_hnc" || sl[i]=="zetaHNC"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_zeta_hnc;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="zeta-closure" || sl[i]=="zeta_closure"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_zeta_closure;
                        sys->cmd_flag_energy_ever_display = true;
                    } else if (sl[i]=="rdf"){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_rdf;
                        sys->cmd_flag_rdf_ever_display = true;
                    } else {
                        char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                        fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined %s item \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, cmd.command==IETCMD_REPORT?"report":"display", sl[i].text, sys->is_log_tty?color_string_end:"");
                        success = false; sl[i].text[sl[i].length] = cc;
                    }
                }
                cmd.step = i_param_list;
            /*} else if (cmd.command==IETCMD_DISPLAY){
                if (StringNS::is_string_number(sl[i])){
                    cmd.command_params_double[i_param_list>0?i_param_list-1:i_param_list] = atof(sl[i].text);
                } else if (sl[i]=="all"||sl[i]=="everything"||sl[i]=="anything"){
                    if (i_param_list+10 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_Euv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_Ef;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_Chandler_G;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_HFE;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = display all : error : too many things to display%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="rdf"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_rdf;
                    sys->cmd_flag_rdf_ever_display = true;
                } else if (sl[i]=="lj"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ulj;
                } else if (sl[i]=="coul"){      if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                } else if (sl[i]=="Euv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Euv;
                } else if (sl[i]=="Ef"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Ef;
                } else if (sl[i]=="energy"){    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_EuvDetail;
                } else if (sl[i]=="Cuv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                } else if (sl[i]=="dN"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN;
                } else if (sl[i]=="dN0"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_DeltaN0;
                } else if (sl[i]=="TS"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                } else if (sl[i]=="entropy"){   if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_TS;
                } else if (sl[i]=="HFE"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_HFE;
                } else if (sl[i]=="SFE"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_HFE;
                } else if (sl[i]=="GGF"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Chandler_G;
                } else if (sl[i]=="Chandler"){  if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Chandler_G;
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : \"%s\" undefined in displaying%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_REPORT){
                if (sl[i]=="all"||sl[i]=="everything"||sl[i]=="anything"){
                    if (i_param_list+2 < MAX_CMD_PARAMS){
                        cmd.command_params_int[i_param_list++] = IETCMD_v_EuvDetail;
                        sys->cmd_flag_energy_ever_display = true;
                        cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                    } else {
                        fprintf(sys->log(), "%s%s : %s[%d] = display all : error : too many things to display%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, sys->is_log_tty?color_string_end:""); success = false;
                    }
                } else if (sl[i]=="Euv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Euv;
                    sys->cmd_flag_energy_ever_display = true;
                } else if (sl[i]=="Ef"){        if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Ef;
                } else if (sl[i]=="energy"){    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_EuvDetail;
                    sys->cmd_flag_energy_ever_display = true;
                } else if (sl[i]=="Cuv"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_cuv;
                } else if (sl[i]=="rdf"){       if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_rdf;
                    sys->cmd_flag_rdf_ever_display = true;
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : \"%s\" undefined in displaying%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;*/
            } else if (cmd.command==IETCMD_HI_SOLVER){
                if (StringNS::is_string_number(sl[i]) && i_param_list<MAX_CMD_PARAMS){
                    cmd.command_params_int[i_param_list] = 1;
                    cmd.command_params_double[i_param_list++] = atof(sl[i].text);
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_LSE){
                if (i_param_list+1 < MAX_CMD_PARAMS){
                    int i_val = 0; double v_val = 0;
                    if (sl[i].text[0]>='A' && sl[i].text[0]<='Z') i_val = sl[i].text[0] - 'A' + 1;
                    if (sl[i].text[0]>='a' && sl[i].text[0]<='z') i_val = sl[i].text[0] - 'a' + 1;
                    if (i+1<nw && StringNS::is_string_number(sl[i+1])){ i++; v_val = atof(sl[i].text); }
                    cmd.command_params_int[i_param_list] = i_val;
                    cmd.command_params_double[i_param_list] = v_val;
                    i_param_list ++;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_BUILD_UUV){
                if (StringNS::is_string_number(sl[i])){
                    if (i_param_list<MAX_CMD_PARAMS){
                        cmd.command_params_double[i_param_list++] = atof(sl[i].text);
                    } else {
                        fprintf(sys->log(), "%s : %s[%d] : too many parameters for Coulomb algorithm \"%s\"\n", software_name, get_second_fn(script_name), script_line, sl[i].text);
                    }
                } else {
                    int esal_here = -1;
                    for (int j=0; esal_here<0 && j<sizeof(CoulAL_names)/sizeof(CoulAL_names[0]); j++) if (sl[i]==CoulAL_names[j]) esal_here = j;
                    for (int j=0; esal_here<0 && j<sizeof(CoulAL_alias)/sizeof(CoulAL_alias[0]); j++) if (sl[i]==CoulAL_alias[j].name) esal_here = j;
                    if (esal_here>=0){
                        cmd.command_params_int[0] = esal_here;
                        cmd.command_params_double[0] = 1;   // mean dielect
                        cmd.command_params_double[1] = 0;   // Yukawa decay length
                    } else {
                        fprintf(sys->log(), "%s : %s[%d] : unknown Coulomb algorithm \"%s\"\n", software_name, get_second_fn(script_name), script_line, sl[i].text);
                        success = false;
                    }
                }
            } else if (cmd.command==IETCMD_TEST || cmd.command==IETCMD_TEST_SAVE){ // test
                if (sl[i]=="coul" || sl[i]=="coulomb"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_ucoul;
                } else if (sl[i]=="hlr"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_hlr;
                } else if (sl[i]=="yukawa"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_Yukawa;
                } else if (sl[i]=="local-coul" || sl[i]=="local_coul" || sl[i]=="local-coulomb" || sl[i]=="local_coulomb"){
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = IETCMD_v_LocalCoulomb;
                } else {
                    char cc = sl[i].text[sl[i].length]; sl[i].text[sl[i].length] = 0;
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined variable \"%s\" for test%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false; sl[i].text[sl[i].length] = cc;
                }
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_CLOSURE || cmd.command==IETCMD_CLOSURE_A){
                int this_closure = closure_from_string(sl[i], CLOSURE_alias, n_CLOSURE_alias);
                if (this_closure<0){
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : undefined closure \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false;
                }
                if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = this_closure;
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_DENSITY || cmd.command==IETCMD_DIELECT || cmd.command==IETCMD_CF || cmd.command==IETCMD_CF_A){
                if (StringNS::is_string_number(sl[i])) if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_double[i_param_list++] = atof(sl[i].text);
                cmd.step = i_param_list;
            } else if (cmd.command==IETCMD_HOLD){
                if (sl[i]=="none" || sl[i]=="clear"){
                    cmd.step = 0; i_param_list = 0; memset(cmd.command_params_int, 0, sizeof(cmd.command_params_int));
                } else if (StringNS::is_string_number(sl[i])){
                    int this_hold = atoi(sl[i].text);
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_int[i_param_list++] = this_hold;
                    cmd.step = i_param_list;
                } else {
                  // search for atoms
                    char buffer[128]; memset(buffer, 0, sizeof(buffer)); memcpy(buffer, sl[i].text, sl[i].length>sizeof(buffer)?sizeof(buffer):sl[i].length);
                    int imole = search_mole_list(sys->atom_list, 0, sys->n_atom_list, buffer, -1);
                    int hold_grp[MAX_SOL]; int n_hold_grp = 0;
                    for (int i=0; i<sys->n_atom_list; i++) if (sys->atom_list[i].iaa==imole) hold_grp[n_hold_grp++] = sys->atom_list[i].grp;
                    for (int i=0; i<n_hold_grp; i++){ bool duplicate = false;
                        for (int j=0; j<i; j++) if (hold_grp[i]==hold_grp[j]) duplicate = true;
                        if (!duplicate && i_param_list<MAX_CMD_PARAMS){
                            cmd.command_params_int[i_param_list++] = hold_grp[i];
                            cmd.step = i_param_list;
                        }
                    }
                    if (n_hold_grp<=0){
                        fprintf(sys->log(), "%s%s : %s[%d][%ld] : -cmd hold: no such molecule: %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                        success = false;
                    }
                }
            } else if (cmd.command==IETCMD_RDF_CONTENT){ cmd.step = 0;
                if (sl[i] == "rdf" || sl[i] == "density"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_rdf;
                } else if (sl[i] == "g" || sl[i] == "guv"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_guv;
                } else if (sl[i] == "h" || sl[i] == "huv"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_huv;
                } else if (sl[i] == "hlr"){
                    cmd.step = 1; cmd.command_params_int[0] = -IETCMD_v_huv;
                } else if (sl[i] == "dd"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_dd;
                } else if (sl[i] == "nphi"){
                    cmd.step = 1; cmd.command_params_int[0] = -IETCMD_v_dd;
                } else if (sl[i] == "c" || sl[i] == "cuv"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_cuv;
                } else if (sl[i] == "ch" || sl[i] == "chuv"){
                    cmd.step = 1; cmd.command_params_int[0] = -IETCMD_v_cuv;
                } else if (sl[i] == "lj"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_ulj;
                } else if (sl[i] == "coul" || sl[i] == "coulomb"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_ucoul;
                } else if (sl[i] == "ff" || sl[i]=="uuv"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_uuv;
                } else if (sl[i] == "ef"){
                    cmd.step = 1; cmd.command_params_int[0] = IETCMD_v_Ef;
                    sys->b_allow_Ecoul0 = true;
                } else {
                    fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : unrecognizable rdf content \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, sl[i].text, sys->is_log_tty?color_string_end:"");
                    success = false;
                }
            } else if (cmd.command==IETCMD_TEMPERATURE){
                if (StringNS::is_string_number(sl[i])){
                    double factor = atof(sl[i].text);
                    if (i_param_list<MAX_CMD_PARAMS) cmd.command_params_double[i_param_list++] = factor;
                    cmd.step = i_param_list;
                }
          // experimental
          #ifdef _EXPERIMENTAL_
            } else if (experimental_analysis_command_params(sl, nw, i, i_param_list, cmd, cmd_type, sys->log(), sys->is_log_tty, script_name, script_line)){
          #endif
          // done
            } else {
                char buffer[512]; memset(buffer, 0, sizeof(buffer)); memcpy(buffer, sl[i].text, sl[i].length>sizeof(buffer)-1?sizeof(buffer)-1:sl[i].length);
                fprintf(sys->log(), "%s%s : %s[%d][%ld] : syntex error : cannot recognize \"%s\"%s\n", sys->is_log_tty?color_string_of_error:"", software_name, script_name, script_line, 1+sl[i].text-line+first_char_offset, buffer, sys->is_log_tty?color_string_end:"");
                success = false;
                cmd_type = -1;
            }
        }
    }

    if (nw>0 && success) sys->insert_command(&cmd, cmd_location-1); //int ic = sys->ncmd-1; printf("sys->cmd=%d, location=%d\n", sys->cmd[ic].command, ic); for (int i=0; i<MAX_CMD_PARAMS; i++) if (sys->cmd[ic].command_params_int[i] || sys->cmd[ic].command_params_double[i]) printf("     param[%d] = %12d %15g\n", i, sys->cmd[ic].command_params_int[i], sys->cmd[ic].command_params_double[i]);

    return success;
}
