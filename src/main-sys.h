//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//----------------------   System Param: IET command list   -----------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

class IET_command {
  public:
    unsigned int command;
      // RISM: ietal+1000
      // HI: hial+1500; other: commands
    int step;
      // RISM: stepmax_rism
      // HI: stepmax_rism
    int time_to_run; // 0: normal, everystep; 1: at the beginning; -1: at the end
    int command_params_int[MAX_CMD_PARAMS];
      // set: parameter id
      // closure: at most 10 molecular closures
      // RISM: ???, ???, ???, force_closure, force_Coulomb, force_rism_dielect
      // load: object index to load
      // save: object index to save
      // display: object index to display
    double command_params_double[MAX_CMD_PARAMS];
      // set: parameter value
      // closure: at most 10 molecular closure parameters
      // save: precision
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------   System Param: The system wide parameters -------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class IET_Param {
  #define MAX_DEBYE_TERMS     6
  #define MP_TASK_PARAMS      5
 // Fixed parameters
  public:   // task params
    FILE * flog; bool is_log_tty; int argc; char ** argv;
    const char * library_path;
    int nt, ntb, ntf; bool mp_by_fork; int nice_level; pid_t pid; //int mp_tasks[MAX_THREADS];
        char * batch_cmd[MAX_THREADS]; int len_batch_cmd;
        char * batch_out[MAX_THREADS]; int len_batch_out;
    bool suspend_calculation; bool is_suspend_calculation; bool listonly, listall;
  public:   // force field global params
    int forcefield_prefix; bool rism_coulomb_renormalization;
    bool mixingrule_sigma_geometric; bool mixingrule_sigma_geometric_specified;
    double rvdw, rcoul, rlocal_coul; double gamma_erf; bool gamma_auto_generate; double ccutoff;
    double rbohr, rbcharge; double epsilon_bohr;
    double mean_dielect; double dielect[MAX_SOL]; double dielect_mol[MAX_SOL]; int n_dielect_mol; //bool lj_is_eletric;
      double rc_yukawafft; double dielect_yukawa; double yukawa_alpha; double dielect_hi;
      bool dielect_from_dipole, dipole_from_dielect;
      double ld_kernel_expand;
    int esal; bool perform_pme;
      double debye_rate[MAX_SOL]; int n_debye_rate; double debye_kappa[MAX_DEBYE_TERMS];
    double temperature, default_temperature;
    double scale_lj, scale_coul, scale_hs;
  public:   // conformation and trajectory params
    bool pbc_x, pbc_y, pbc_z;
    __REAL__ drrism, drhi; int nr[3]; int nv/*site number in total*/;
    int output_significant_digits;
    double time_begin, time_end, time_step; bool handling_xtc;
    bool mode_test;
    int detail_level, debug_level; bool ignore_memory_capacity;
        bool debug_show_crc; bool debug_xvv; int debug_xvv_terms;
   #ifdef _INTERACTIVE_
    bool allow_interactive;
   #endif
  public:   // solvents
    int gvv_specification; bool gvv_is_left_aligned;
    double xvv_k_shift; double xvv_extend; const char * xvv_scheme; double * xvv_scale; int n_xvv_scale;
    AtomSite av[MAX_SOL]; int nav; double density_av[MAX_SOL];
    int nvm; int vmmap[MAX_SOL]; int vmmapi[MAX_SOL][2];
      int nmv; double density_mv[MAX_SOL]; int nmvam; double density_hi;
      int nmvb; double bulk_density_mv[MAX_SOL];
    double dipole_mv[MAX_SOL]; int ndipole;
    double zeta_scaling_factor[MAX_SOL]; int n_zeta_scaling_factor;
    bool check_zero_xvv;
    double degree_of_freedom[MAX_SOL]; int n_degree_of_freedom;
  public:   // solutes
    SoluteAtomSite * as; int nas; int nasmax; // solutes
    PDBAtomSet traj; // trajectory of solutes
    Vector box; bool force_box;  // default box
  public:   // HI params
    double llambda[MAX_SOL]; int nllambda; // ln\lambda is the integral constants
    double lse_a, lse_b; int calc_ab_automatically, calc_nbulk_automatically;
    double nbulk[MAX_SOL]; int nnbulk; // atomised mole fraction
        double nbulk_rism[MAX_SOL]; bool enable_nbulk_rism;
    double ccutoff_hi;
    #ifdef _EXPERIMENTAL_
      IET_Param_Exp ex; // experimental features
    #endif
  public:   // SCF iteration tolerances
    __REAL__ errtolhi, errtolrism; __REAL__ delhi, delrism; int ndiis_rism, ndiis_hi;
  public:   // output
    double gcutoff_liquid_occupation; double gcutoff_ef_occupation;
    int out_rdf_bins; RDFGroup rdf_grps[MAX_RDF_GRPS]; int n_rdf_grps; int rdf_content;
  public:   // experimental features
    bool hlr_no_hi;
    bool use_homogeneous_rism;
    double closure_enhance_level; double closure_enhance_cutoff;
    Vector external_electrostatic_field; bool do_external_electrostatic_field;
  public:   // name list and mappings
    SYSITEM_AtomNameList * atom_list; int n_atom_list; int nmax_atom_list;
    SYSITEM_BondList * bond_list; int n_bond_list; int nmax_bond_list;
    SYSITEM_PairMapping * gvv_map; int n_gvv_map; int nmax_gvv_map;
    SYSITEM_ZetaList * zeta_list; int n_zeta_list; int nmax_zeta_list; bool zeta_list_allow_missing;
    double * sigma_list; int n_sigma_list; int nmax_sigma_list;
    double * epsilon_list; int n_epsilon_list; int nmax_epsilon_list;
 // ------------------------------------------
 // SCF Iteration parameters and command queue
  public:   // RISM params
    int stepmax_rism;
    unsigned int ietal;
    int closures[MAX_SOL]; __REAL__ closure_factors[MAX_SOL];
  public:   // HI params
    int stepmax_hi;
    unsigned int hial; // HI algorithm
    double ucutoff_hs;
  public:   // command quene and command options
    IET_command * cmd; int ncmd, ncmd_max;
    bool cmd_flag_energy_ever_display;
    bool cmd_flag_rism_ever_performed;
    bool cmd_flag_rdf_ever_display;
    bool cmd_rbc_ljr_allowed;
 // ------------------------------------------
 // Other parameters
  public:   // output and compression
    int output_override; // -1: force override; 0: no override; 1: append
    int output_compress_level; // 0: no compression; 1: fast; 2: best; -1: default
    int compress_page_size;
    char output_realnumber_format[20]; char output_title_format[20]; int output_realnumber_length;
  public:   // other flags
    int ff_advanced_flags;
    bool b_output_filename_autogenerated;
    bool b_allow_r2uvmin, b_allow_Ecoul0;
  public:   // buffers for calculation
    int _hold_list[MAX_SOL]; int _n_hold_list;  // hold sites for HI and RISM
 // ------------------------------------------
 // ------------------------------------------
 // ------------------------------------------
 // ------------------------------------------
  public:
    void init(int _argc, char ** _argv){
        argc = _argc; argv = _argv;
        flog = stderr; is_log_tty = true;
        library_path = nullptr;
        nt = 1; ntb = 1; ntf = 1; nice_level = 0; pid = getpid(); // for (int i=0; i<MAX_THREADS; i++) mp_tasks[i] = 0;
        #ifdef _LOCALPARALLEL_PTHREAD_
            mp_by_fork = false;
        #else
            mp_by_fork = true;
        #endif
        memset(batch_cmd, 0, sizeof(batch_cmd)); len_batch_cmd = 0;
        memset(batch_out, 0, sizeof(batch_out)); len_batch_out = 0;

        forcefield_prefix = 0;
        rism_coulomb_renormalization = true;
        mixingrule_sigma_geometric = false; mixingrule_sigma_geometric_specified = false;
        suspend_calculation = false; is_suspend_calculation = false;
        listonly = false; listall = false;

        for (int i=0; i<MAX_SOL; i++){ closures[i] = CLOSURE_HNC; closure_factors[i] = 1; }
        drrism = drhi = 0; nr[0] = nr[1] = nr[2] = 100; nv = nvm = 0;
        //output_significant_digits = 160;
        output_significant_digits = 7;
        rbohr = 0.052911; rbcharge = 0; epsilon_bohr = 1/MACHINE_REASONABLE_ERROR;
        rc_yukawafft = -1; dielect_yukawa = mean_dielect = dielect_hi = 1; yukawa_alpha = 0;
          for (int i=0; i<MAX_SOL; i++) dielect[i] = dielect_mol[i] = 1; n_dielect_mol = 0; //lj_is_eletric = false;
          dielect_from_dipole = dipole_from_dielect = false;
          ld_kernel_expand = 1;
        esal = CoulAL_Coulomb; perform_pme = true;
            for (int i=0; i<MAX_SOL; i++) debye_rate[i] = 0; n_debye_rate = MAX_SOL; for (int i=0; i<MAX_DEBYE_TERMS; i++) debye_kappa[i] = 0;
        gvv_specification = 0; gvv_is_left_aligned = true; xvv_k_shift = 0; xvv_extend = 0; xvv_scheme = "so"; xvv_scale = nullptr; n_xvv_scale = 0;
        temperature = 298; default_temperature = 120.27;
        traj.count = 0; traj.atom = nullptr; traj.box = Vector(0,0,0); box = Vector(0, 0, 0); force_box = false;
        as = nullptr; nas = nasmax = nav = nmv = nmvb = 0; memset(av, 0, sizeof(av));
        rvdw = rcoul = 1; rlocal_coul = -1; gamma_erf = 2/rcoul; gamma_auto_generate = true;
        ccutoff = 3;
        pbc_x = pbc_y = pbc_z = true;
        scale_lj = scale_coul = scale_hs = 1;
        time_begin = time_end = time_step = 0; handling_xtc = false;
        detail_level = 1; debug_level = 0; mode_test = false; ignore_memory_capacity = false; debug_show_crc = false; debug_xvv = false; debug_xvv_terms = 4;
       #ifdef _INTERACTIVE_
        allow_interactive = false;
       #endif

        stepmax_rism = 100;
        errtolrism = check_error_tol(0); delrism = 0.7; ndiis_rism = 0;
        ietal = IETAL_NONE;

        stepmax_hi = 1000;
        errtolhi = check_error_tol(0); delhi = 1; ndiis_hi = 0;
        hial = HIAL_NONE; ucutoff_hs = 5;
        nmvam = 0; density_hi = 33.4;
        ndipole = 0;
        for (int i=0; i<MAX_SOL; i++) zeta_scaling_factor[i] = 1; n_zeta_scaling_factor = 0;
        check_zero_xvv = true;
        for (int i=0; i<MAX_SOL; i++) degree_of_freedom[i] = 3; n_degree_of_freedom = 0;
        for (int i=0; i<MAX_SOL; i++) llambda[i] = 0; nllambda = 0; nbulk[0] = 1; nnbulk = 0;
            enable_nbulk_rism = false; for (int i=0; i<MAX_SOL; i++) nbulk_rism[0] = 1;
        ccutoff_hi = 3;
        lse_a = 0.3; lse_b = 54; calc_ab_automatically = 2; calc_nbulk_automatically = false;
        #ifdef _EXPERIMENTAL_
            ex.init();
        #endif

        gcutoff_liquid_occupation = MACHINE_REASONABLE_ERROR;
        gcutoff_ef_occupation = MACHINE_REASONABLE_ERROR;
        out_rdf_bins = 50; n_rdf_grps = 0; memset(rdf_grps, 0, sizeof(rdf_grps)); rdf_content = 0;

        hlr_no_hi = false;
        use_homogeneous_rism = false;
        closure_enhance_level = 1; closure_enhance_cutoff = -1;
        external_electrostatic_field = Vector(0, 0, 0); do_external_electrostatic_field = false;

        atom_list = nullptr; n_atom_list = nmax_atom_list = 0;
        bond_list = nullptr; n_bond_list = nmax_bond_list = 0;
        gvv_map = nullptr; n_gvv_map = nmax_gvv_map = 0;
        zeta_list = nullptr; n_zeta_list = nmax_zeta_list = 0; zeta_list_allow_missing = true;
        sigma_list = nullptr; n_sigma_list = 0; nmax_sigma_list = 0;
        epsilon_list = nullptr; n_epsilon_list = 0; nmax_epsilon_list = 0;

       /*#ifdef _LOCALPARALLEL_
        advanced_paralling_ff_batch = true; advanced_paralling_fft_assign_node0 = true; advanced_mp_seq = false;
       #endif*/

        cmd = nullptr; ncmd = 0; ncmd_max = 0;
        cmd_flag_energy_ever_display = false;
        cmd_flag_rism_ever_performed  = false;
        cmd_flag_rdf_ever_display = false;
        cmd_rbc_ljr_allowed = false;

        output_override = 1; output_compress_level = 0; compress_page_size = 4096;
        strncpy(output_realnumber_format, " %11g", sizeof(output_realnumber_format));
        strncpy(output_title_format, " %11s", sizeof(output_realnumber_format));
        output_realnumber_length = 14;

        ff_advanced_flags = 0;
        b_output_filename_autogenerated = false;
        b_allow_r2uvmin = false;
        b_allow_Ecoul0 = false;

        _n_hold_list = 0;
    }
    FILE * log(){ return flog? flog : stderr; }
    bool prepare_traj_input(char * filename, FILE * flog){ bool success = true;
        char * ext = (char*)""; for (int i=strlen(filename)-1; i>0; i--) if (filename[i]=='.'){ ext = &filename[i]; break; }
        if (!ext[0]) { fprintf(flog, "%s : error : -f file has no extension\n", software_name); success = false; }
        else if (StringNS::string(ext) == ".pdb") success = prepare_input(filename, nullptr, nullptr);
        else if (StringNS::string(ext) == ".gro") success = prepare_input(nullptr, filename, nullptr);
        else if (StringNS::string(ext) == ".xtc") success = prepare_input(nullptr, nullptr, filename);
        else { fprintf(flog, "%s : error : unknown -f type %s\n", software_name, ext); success = false; }
        if (StringNS::string(ext) == ".xtc") handling_xtc = true; else handling_xtc = false;
        return success;
    }
    void insert_command(IET_command * cmdi, int index=-1){
        int ncmd_old = ncmd; int ncmd_max_old = ncmd_max; if (index<0) index = ncmd;
        if (ncmd<=index) ncmd = index+1;
        if (ncmd>=ncmd_max){
            ncmd_max = ncmd + 50;
            IET_command * croot = (IET_command *)malloc(sizeof(IET_command)*ncmd_max);
            if (croot){
                memset(croot, 0, sizeof(IET_command)*ncmd_max);
                if (cmd){
                    memcpy(croot, cmd, sizeof(IET_command)*ncmd_old);
                    free(cmd);
                }
                cmd = croot;
            } else ncmd_max = ncmd_max_old;
        }
        if (ncmd<ncmd_max){
            memcpy(&cmd[index], cmdi, sizeof(IET_command));
        } else { fprintf(stderr, "%s : malloc failure\n", software_name); }
    }
    void dump_text(FILE * o, const char * lineprefix=(char*)"", const char * headerprefix=nullptr){
        bool istty = isatty(fileno(o)); if (!headerprefix) headerprefix = lineprefix;
        fprintf(o, "%s%s[solvent]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
        fprintf(o, "%s# Global Settings%s\n", istty?"\33[37m":"", istty?"\33[0m":"");
          if (mixingrule_sigma_geometric_specified){
              fprintf(o, "%s%s\n", lineprefix, mixingrule_sigma_geometric?"geometric-sigma":"arithmetic-sigma");
          }
          fprintf(o, "%snr                  %dx%dx%d\n", lineprefix, nr[0], nr[1], nr[2]);
          fprintf(o, "%spbc                 ", lineprefix); if (pbc_x) fprintf(o, "x"); if (pbc_y) fprintf(o, "y"); if (pbc_z) fprintf(o, "z"); if (!pbc_x && !pbc_y && !pbc_z) fprintf(o, "none"); fprintf(o, "\n");
        fprintf(o, "%s# Force Field Setups%s\n", istty?"\33[37m":"", istty?"\33[0m":"");
          fprintf(o, "%scoulomb             %s\n", lineprefix, CoulAL_names[esal]);
          fprintf(o, "%srvdw                %g\n", lineprefix, rvdw);
          fprintf(o, "%srcoul               %g\n", lineprefix, rcoul);
          if (perform_pme) fprintf(o, "%s%sdo-PME%s\n", lineprefix, istty?"\33[37m":"", istty?"\33[0m":""); else fprintf(o, "%sno-PME\n", lineprefix);
          if (!gamma_auto_generate) fprintf(o, "%spme-gamma           %g\n", lineprefix, gamma_erf); else fprintf(o, "%s%s#gamma              %g%s\n", lineprefix, istty?"\33[37m":"", gamma_erf, istty?"\33[0m":"");
          fprintf(o, "%sccutoff             %g\n", lineprefix, ccutoff);
          //fprintf(o, "%sdielect             %.4e\n", lineprefix, dielect);
          //fprintf(o, "%s%-8s           ", lineprefix, lj_is_eletric?"dielect":"dielectx"); for (int i=0; i<nmv; i++) fprintf(o, " %f", dielect_mol[i]); fprintf(o, " \n");
          fprintf(o, "%sdielect            ", lineprefix); for (int i=0; i<nmv; i++) fprintf(o, " %g", dielect_mol[i]); fprintf(o, " \n");
          //fprintf(o, "%sdebye-decay-rate   ", lineprefix); for (int i=0; i<nmv; i++) fprintf(o, " %g", debye_rate[i]); fprintf(o, " \n");
          fprintf(o, "%stemperature         %g\n", lineprefix, temperature);
          fprintf(o, "%sdefault-temperature %g\n", lineprefix, default_temperature);
        fprintf(o, "%s# Advanced settings%s\n", istty?"\33[37m":"", istty?"\33[0m":"");
          fprintf(o, "%sdetail-level        %d\n", lineprefix, detail_level);
          fprintf(o, "%sdebug-level         %d\n", lineprefix, debug_level);
          #ifdef _INTERACTIVE_
            fprintf(o, "%s%s%s%s\n", lineprefix, istty?"\33[37m":"", allow_interactive?"interactive":"#interactive", istty?"\33[0m":"");
          #endif

        if (output_significant_digits==7){
            fprintf(o, "%ssignificant-digits  float\n", lineprefix);
        } else if (output_significant_digits==18){
            fprintf(o, "%ssignificant-digits  double\n", lineprefix);
        } else {
            fprintf(o, "%ssignificant-digits  %d\n", lineprefix, output_significant_digits);
        }

        fprintf(o, "%s%s[hi]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
          fprintf(o, "%serrtolhi            %g\n", lineprefix, errtolhi);
          fprintf(o, "%sdelhi               %g\n", lineprefix, delhi);
          fprintf(o, "%sndiishi             %d\n", lineprefix, ndiis_hi);
          fprintf(o, "%slse_a               %g\n", lineprefix, lse_a);
          fprintf(o, "%slse_b               %g\n", lineprefix, lse_b);
          if (nnbulk>0){
              fprintf (o, "%s%-6s             ", lineprefix, (calc_nbulk_automatically&&istty)?"\33[37m#nbulk":"nbulk");
              for (int i=0; i<nnbulk; i++) fprintf(o, " %g", nbulk[i]);
              fprintf(o, "%s\n", (calc_nbulk_automatically&&istty)?"\33[0m":"");
          }
          if (ndipole>0){ fprintf(o, "%sdipole             ", lineprefix); for (int i=0; i<ndipole; i++) fprintf(o, " %g", dipole_mv[i]); fprintf(o, "\n"); } else fprintf(o, "%s%s# no dipoles\n%s", lineprefix, istty?"\33[37m":"", istty?"\33[0m":"");
          fprintf(o, "%sdrhi                %g\n", lineprefix, drhi);
          if (szfn_zeta[0]) fprintf(o, "%szeta                %s%s%s\n", lineprefix, istty?prompt_path_prefix:"", (szfn_zeta), istty?prompt_path_suffix:"");
          if (n_zeta_scaling_factor>0){
              fprintf(o, "%sscale.zeta         ", lineprefix); for (int i=0; i<n_zeta_scaling_factor; i++) fprintf(o, " %g", zeta_scaling_factor[i]); fprintf(o, "\n");
            }
        fprintf(o, "%s%s[rism]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
          fprintf(o, "%serrtolrism          %g\n", lineprefix, errtolrism);
          fprintf(o, "%sdelrism             %g\n", lineprefix, delrism);
          fprintf(o, "%sndiisrism           %d\n", lineprefix, ndiis_rism);
          //fprintf(o, "%sclosure            ", lineprefix); for (int i=0; i<nv; i++) fprintf(o, " %s", CLOSURE_name[closures[i]]); fprintf(o, "\n");
          //fprintf(o, "%sclosure-factor     ", lineprefix); for (int i=0; i<nv; i++) fprintf(o, " %g", closure_factors[i]); fprintf(o, "\n");
          fprintf(o, "%sclosure-enhance     %g\n", lineprefix, closure_enhance_level);
          fprintf(o, "%srbohr               %f\n", lineprefix, rbohr);
          if (xvv_k_shift!=0) fprintf(o, "%sxvv-k-shift         %g\n", lineprefix, xvv_k_shift);
          if (xvv_extend!=0) fprintf(o, "%sxvv-extend          %g\n", lineprefix, xvv_extend);
          fprintf(o, "%sdensity            ", lineprefix); for (int i=0; i<nmv; i++) fprintf(o, " %g", density_mv[i]); fprintf(o, "\n");
        fprintf(o, "%sbulk-density       ", lineprefix); for (int i=0; i<nmv; i++) fprintf(o, " %g", bulk_density_mv[i]); fprintf(o, "\n");
          fprintf(o, "%sdrrism              %g\n", lineprefix, drrism);
          if (gvv_specification==0){
          } else if (gvv_specification>0 && n_szfn_gvvs>0){
              fprintf(o, "%sgvv                 %g", lineprefix, drrism);
              for (int j=0; j<n_szfn_gvvs; j++) fprintf(o, " %s%s%s", istty?prompt_path_prefix:"", &szfn_gvvs[j*MAX_PATH], istty?prompt_path_suffix:"");
              fprintf(o, "\n");
          } else if (gvv_specification<0){
              fprintf(o, "%sgvv                 %g", lineprefix, drrism);
              for (int j=0; j<n_szfn_gvvs; j++) fprintf(o, " %s%s%s", istty?prompt_path_prefix:"", &szfn_gvvs[j*MAX_PATH], istty?prompt_path_suffix:"");
              fprintf(o, "\n");
          }
          //fprintf(o, "%s  steps     max RISM: %d, max DIIS: %d\n", lineprefix, stepmax_rism, ndiis_rism);
      if (listall){
        if (nav>0){
            fprintf(o, "%s%s[atom]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
            int id_atom = 1;
            if (true){
                fprintf(o, "%s%s# %6s %6s %4s %4s %9s %9s %9s%s\n", lineprefix, istty?"\33[37m":"", "atom", "mole", "id", "grp", "charge", "sigma", "epsiln", istty?"\33[0m":"");
                for (int i=0; i<nav; i++){
                  for (int im=0; im<av[i].multi; im++){
                    char * atom_name = av[i].name; char * mole_name = av[i].mole;
                    if (atom_list) for (int i=0; i<n_atom_list; i++) if (atom_list[i].index == id_atom){ atom_name = atom_list[i].name; mole_name = atom_list[i].mole; }
                    fprintf(o, "%s  %s%6s %6s %4d %4d %9.4g %9.4g %9.4g%s", lineprefix, (istty&&av[i].is_key)?"\33[35m":"", atom_name, mole_name, /*av[i].id+im*/id_atom++, av[i].grp, av[i].charge, av[i].sigma, av[i].epsilon, (istty&&av[i].is_key)?"\33[0m":"");
                    fprintf(o, "\n");
                  }
                }
            }
            //for (int i=0; i<nav; i++) fprintf(o, "%s%6s %6s %9.4f %9.4f %9.4f %9.4f\n", lineprefix, av[i].mole, av[i].name, av[i].charge, av[i].sigma, av[i].epsilon, density_av[i]);
        }
        if (n_bond_list>0 && bond_list){
            fprintf(o, "%s%s[bond]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
            for (int ib=0; ib<n_bond_list; ib++){
                char * namei = (char*)""; char * namej = (char*)"";
                char * molei = (char*)""; char * molej = (char*)"";
                if (atom_list){
                    for (int i=0; i<n_atom_list; i++) if (atom_list[i].grp == bond_list[ib].grpi){ namei = atom_list[i].name; molei = atom_list[i].mole; break; }
                    for (int i=0; i<n_atom_list; i++) if (atom_list[i].grp == bond_list[ib].grpj){ namej = atom_list[i].name; molej = atom_list[i].mole; break; }
                }
                char buffer[2][64];
                if (namei[0]) snprintf(buffer[0], sizeof(buffer[0]), "%s.%s", molei[0]?molei:"*", namei); else snprintf(buffer[0], sizeof(buffer[0]), "%d", bond_list[ib].grpi);
                if (namej[0]) snprintf(buffer[1], sizeof(buffer[1]), "%s.%s", molej[0]?molej:"*", namej); else snprintf(buffer[1], sizeof(buffer[1]), "%d", bond_list[ib].grpj);
                fprintf(o, "%s %8s %8s %8g\n", lineprefix, buffer[0], buffer[1], bond_list[ib].bond);
            }
        }
        if (n_gvv_map>0 && gvv_map){
            fprintf(o, "%s%s[gvv_map]%s\n", istty?"\33[31m":"", headerprefix, istty?"\33[0m":"");
            for (int ib=0; ib<n_gvv_map; ib++){
                char * namei = (char*)""; char * namej = (char*)"";
                char * molei = (char*)""; char * molej = (char*)"";
                if (atom_list){
                    for (int i=0; i<n_atom_list; i++) if (atom_list[i].grp == gvv_map[ib].grpi){ namei = atom_list[i].name; molei = atom_list[i].mole; break; }
                    for (int i=0; i<n_atom_list; i++) if (atom_list[i].grp == gvv_map[ib].grpj){ namej = atom_list[i].name; molej = atom_list[i].mole; break; }
                }
                char buffer[2][64];
                if (namei[0]) snprintf(buffer[0], sizeof(buffer[0]), "%s.%s", molei[0]?molei:"*", namei); else snprintf(buffer[0], sizeof(buffer[0]), "%d", gvv_map[ib].grpi);
                if (namej[0]) snprintf(buffer[1], sizeof(buffer[1]), "%s.%s", molej[0]?molej:"*", namej); else snprintf(buffer[1], sizeof(buffer[1]), "%d", gvv_map[ib].grpj);
                fprintf(o, "%s %8s %8s %4d\n", lineprefix, buffer[0], buffer[1], gvv_map[ib].col);
            }
        }
      }
    }
};
