const char * software_name = "rismhi3d";
const char * software_version = "a.240.1465";
const char * copyright_string = "(c) Cao Siqin";

#define     __REAL__    double
#define     MACHINE_REASONABLE_ERROR    1e-12
#define _TTYPROMPTCOLOR_

#define MAX_SOL                     100     // Max atom site number
#define MAX_CMD_PARAMS              20      // Max parameter number for a command
#define MAX_THREADS                 100     // Max number of forks or threads
#define MAX_DIIS                    100     // Max DIIS steps
#define MAX_INCLUDE_RECURSIVE       20      // maximum include recursive levels

#include    "header.h"
#if defined(_GROMACS4_) || defined(_GROMACS5_) || defined(_GROMACS2016_) || defined(_GROMACS2018_)
  #define _GROMACS_
#endif
#include    <errno.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <stdint.h>
#include    <string.h>
#include    <math.h>
#include    <signal.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <unistd.h>
#include    <time.h>
#include    <libgen.h>
#include    <sys/time.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>
#if defined(_LOCALPARALLEL_PTHREAD_) || defined(_INTERACTIVE_)
  #include    <pthread.h>
#endif
#include    "fftw3.h"
#ifdef _LIBZ_
  #include  <zlib.h>
#endif
#ifndef nullptr
  #define nullptr NULL
#endif

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef _TTYPROMPTCOLOR_
    const char * prompt_path_prefix = "\33[4m";
    const char * prompt_path_suffix = "\33[0m";
    const char * prompt_comment_prefix = "\33[37m";
    const char * prompt_comment_suffix = "\33[0m";
    const char * prompt_highlight_prefix = "\33[1;31m";
    const char * prompt_highlight_suffix = "\33[0m";

    const char * color_string_of_echo    = "\33[7m"; //"\33[0;30;47m"
    const char * color_string_of_warning = "\33[0;30;103m"; // "\33[38;5;232;48;5;226m"; // "\33[0;44;38;5;228m";
    const char * color_string_of_error   = "\33[0;93;41m"; // "\33[0;41;38;5;228m"; //"\33[0;31;48;5;226m"
    const char * color_string_end        = "\33[0m";
#else
    const char * prompt_path_prefix = "";
    const char * prompt_path_suffix = "";
    const char * prompt_comment_prefix = "";
    const char * prompt_comment_suffix = "";
    const char * prompt_highlight_prefix = "";
    const char * prompt_highlight_suffix = "";

    const char * color_string_of_echo    = "";
    const char * color_string_of_warning = "";
    const char * color_string_of_error   = "";
    const char * color_string_end        = "";
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define PI  3.1415926535897932384626433832795
#define EE  2.7182818284590452353602874713527
#define COULCOOEF 138.9354846
#define MAX_PATH                1024
#define MAX_WORD                200
#define MAXTEXTFILEMB           16
#define MAX_RDF_GRPS            1000        // Max number of groups to output RDF
#define INITIAL_SOLUTE_ATOMS    500
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef _LOCALPARALLEL_
#define MPTASK_NONE             0
#define MPTASK_TERMINATE        -1
#define MPTASK_FFTW             11
#define MPTASK_MERGE_FFT_DATA   12
#define MPTASK_FFSR             13
#define MPTASK_MERGE_FF_DATA    14
#define MPTASK_RISM_CLOSURE     21
#define MPTASK_HI_SOLVER        22
#define MPTASK_HI_POTENTIAL     23
#define MPTASK_DIIS             24
#define MPTASK_DIIS_WEIGHT      25
#define MPTASK_DIIS_STEPIN      26
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include    "crc32_zlib.h"
#include    "Element.h"
#include    "String2.cpp"
#include    "Vector.cpp"
#include    "PDBAtom.cpp"
#include    "read_frame_abr.cpp"
#include    "common.cpp"
#include    "main-common.cpp"
#include    "main-matrix.cpp"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define FFPREFIX_NONE           0
#define FFPREFIX_AMBER          1
#define FFPREFIX_GAFF           2
#define FFPREFIX_OPLS           3
const char * FFPREFIX_names[] = { "amber", "gaff", "opls" };
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
struct STKeywordTableUnit { int id; const char * name; };
#define CoulAL_NONE             0
#define CoulAL_Coulomb          0
#define CoulAL_Dielect          1
#define CoulAL_YukawaFFT        2
#ifdef _EXPERIMENTAL_
  #define CoulAL_Dipole           6
  #define CoulAL_Local            7
  #define CoulAL_Local_Dielect    8
#endif
const char * CoulAL_names[] = { "Coulomb", "dielect", "YukawaFFT", "", "", "", "dipole", "local", "dlocal" };
STKeywordTableUnit CoulAL_alias [] = {
    { CoulAL_Coulomb            , "Coul" },
    { CoulAL_Dielect            , "dielectric" },
    { CoulAL_YukawaFFT          , "Yukawa" },
    { CoulAL_YukawaFFT          , "YukawaFT" },
    { CoulAL_YukawaFFT          , "YukawaDFT" },
  #ifdef _EXPERIMENTAL_
    { CoulAL_Dipole             , "dipole-ren" },
    { CoulAL_Dipole             , "dipole_ren" },
    { CoulAL_Dipole             , "dipole-renormalization" },
    { CoulAL_Dipole             , "dipole_renormalization" },
    { CoulAL_NONE               , "none" },
    { CoulAL_NONE               , "no" }
  #endif
};
#define CLOSURE_NONE            0
#define CLOSURE_HNC             1
#define CLOSURE_MSA             2
#define CLOSURE_KGK             3
#define CLOSURE_PY              4
#define CLOSURE_D2              5
#define CLOSURE_HNCB            6
#define CLOSURE_PLHNC           7
#define CLOSURE_HARDSPHERE      8
#define CLOSURE_KH              11
#define CLOSURE_PSE2            12
#define CLOSURE_PSE3            13
#define CLOSURE_PSE4            14
#define CLOSURE_PSE5            15
#define CLOSURE_PSE6            16
#define CLOSURE_PSE7            17
#define CLOSURE_PSE8            18
#define CLOSURE_PSE9            19
#define CLOSURE_PSE10           20
#define CLOSURE_MS              21
#define CLOSURE_BPGGHNC         22
#define CLOSURE_VM              23
#define CLOSURE_MP              24
#define CLOSURE_MHNC            25
const char * CLOSURE_name[100];
STKeywordTableUnit CLOSURE_alias[200] = {
  // key names here
    { CLOSURE_NONE              , "none" },
    { CLOSURE_HNC               , "HNC" },
    { CLOSURE_MSA               , "MSA" },
    { CLOSURE_KGK               , "KGK" },
    { CLOSURE_PY                , "PY" },
    { CLOSURE_D2                , "D2" },
    { CLOSURE_HNCB              , "HNCB" },
    { CLOSURE_PLHNC             , "PLHNC" },
    { CLOSURE_HARDSPHERE        , "HS" },
    { CLOSURE_KH                , "KH" },
    { CLOSURE_PSE2              , "PSE2" },
    { CLOSURE_PSE3              , "PSE3" },
    { CLOSURE_PSE4              , "PSE4" },
    { CLOSURE_PSE5              , "PSE5" },
    { CLOSURE_PSE6              , "PSE6" },
    { CLOSURE_PSE7              , "PSE7" },
    { CLOSURE_PSE8              , "PSE8" },
    { CLOSURE_PSE9              , "PSE9" },
    { CLOSURE_PSE10             , "PSE10" },
    { CLOSURE_MS                , "MS" },
    { CLOSURE_BPGGHNC           , "BPGGHNC" },
    { CLOSURE_VM                , "VM" },
    { CLOSURE_MP                , "Marucho-Pettitt" },
    { CLOSURE_MHNC              , "MHNC" },
  // alias in the following
    { CLOSURE_NONE              , "no" },
    { CLOSURE_HNC               , "Hyper-Netted-Chain" },
    { CLOSURE_MSA               , "Mean-Spherical-Approx" },
    { CLOSURE_MSA               , "Mean-Spherical-Approximation" },
    { CLOSURE_PY                , "Percus-Yevick" },
    { CLOSURE_D2                , "sD2" },
    { CLOSURE_HARDSPHERE        , "hardsphere" },
    { CLOSURE_KH                , "Kovalenko-Hirata" },
    { CLOSURE_KH                , "PLHNC" },
    { CLOSURE_KH                , "partial-linear-HNC" },
    { CLOSURE_KH                , "PSE1" },
    { CLOSURE_VM                , "Verlet" },
    { CLOSURE_VM                , "Verlet-Modified" },
    { CLOSURE_VM                , "Modified-Verlet" },
    { CLOSURE_MS                , "Martynov-Sarkisov" },
    { CLOSURE_BPGGHNC           , "BPGG" }
};
int n_CLOSURE_alias = sizeof(CLOSURE_alias)/sizeof(CLOSURE_alias[0]);
#define IETAL_NONE      0
#define IETAL_SSOZ      1
#define IETAL_RRISM     2
#define IETAL_VRISM     3
#define IETAL_IRISM     4
const char * IETAL_name[] = { "", "SSOZ", "RISM", "VRISM", "IRISM" };
#define HIAL_NONE       0
#define HIAL_CUTOFF     1
#define HIAL_ICUTOFF    2
#define HIAL_HI         4
#define HIAL_HSHI       4
#define HIAL_DPHI       5
#define HIAL_DNHI       6
#define HIAL_PLHI       7
#define HIAL_RES_PLHI   8
const char * HIAL_name[9] = { "", "CutoffHI", "iCutoffHI", "", "HSHI", "DPHI", "DNHI", "PLHI", "rPLHI" };
#define GUVMAL_THETA    1
#define GUVMAL_GUV      2
const char * GUVMAL_name[3] = { "", "Theta", "guv" };
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// for BG: integration pathway
#define PERFORM_3DBG_PATHWAY_LJ         1
#define PERFORM_3DBG_PATHWAY_LJ_2       2
#define PERFORM_3DBG_PATHWAY_Coul       3
#define PERFORM_3DBG_PATHWAY_Coul_2     4
#define PERFORM_3DBG_PATHWAY_LJ_Coul    5
#define PERFORM_3DBG_PATHWAY_HS         6
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define IETCMD_NOP          1
#define IETCMD_END          2
#define IETCMD_DONE         3
#define IETCMD_CLEAR        4
#define IETCMD_RESET        5
#define IETCMD_SET          10
#define IETCMD_SCALE        11
#define IETCMD_LOAD         12
#define IETCMD_SAVE         13
#define IETCMD_SAVE_FILTER  14
#define IETCMD_DISPLAY      15
#define IETCMD_REPORT       16
#define IETCMD_CLOSURE      21
#define IETCMD_CLOSURE_A    22
#define IETCMD_CF           23
#define IETCMD_CF_A         24
#define IETCMD_dielect      25
#define IETCMD_density      26
#define IETCMD_BUILD_FF     27
#define IETCMD_TI           28
#define IETCMD_TEST         98
#define IETCMD_TEST_SAVE    99
#define IETCMD_v_box            100
#define IETCMD_v_temperature    101
#define IETCMD_v_Coulomb        102
#define IETCMD_v_dion           103
#define IETCMD_v_dielect_y      104
#define IETCMD_v_rbohr          105
#define IETCMD_v_cmd            200
#define IETCMD_v_uuv            201
#define IETCMD_v_ulr            202
#define IETCMD_v_ulj            203
#define IETCMD_v_ucoul          204
#define IETCMD_v_ucoul2         205
#define IETCMD_v_ucoulsr        206
#define IETCMD_v_ucoullr        207
#define IETCMD_v_unh            208
#define IETCMD_v_unhlj          209
#define IETCMD_v_unhcoul        221
#define IETCMD_v_dd             222
#define IETCMD_v_huv            223
#define IETCMD_v_hlr            224
#define IETCMD_v_guv            225
#define IETCMD_v_cuv            226
#define IETCMD_v_clr            227
#define IETCMD_v_rmin           228
#define IETCMD_v_rdf            231
#define IETCMD_v_Euv            235
#define IETCMD_v_Ef             236
#define IETCMD_v_EuvDetail      237
#define IETCMD_v_DeltaN         238
#define IETCMD_v_DeltaN0        239
#define IETCMD_v_TS             240
#define IETCMD_v_rism_dielect   241
#define IETCMD_v_HFE            242
#define IETCMD_v_Chandler_G     243
#define IETCMD_v_EuvRep         244
#define	IETCMD_v_CuvRep         245
#define IETCMD_v_Mayer          246
#define IETCMD_v_ddp            247
// test commands
#define IETCMD_v_Yukawa         3001
#define IETCMD_v_LocalCoulomb   3002
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
char info_file_name[MAX_PATH];
const char * help_search_str = nullptr;
// path
char szfn_path[MAX_PATH];
char hostname[MAX_PATH], username[MAX_PATH];
// input file
char szfn_xtc[MAX_PATH];
char szfn_solute[MAX_PATH];
char szfn_gvv[MAX_PATH];
char szfn_zeta[MAX_PATH];
// log file
char szfn_log[MAX_PATH];
// TENSOR4D file
char szfn_in[MAX_PATH];
char szfn_out[MAX_PATH]; FILE * file_out = nullptr;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef _EXPERIMENTAL_
    #include "main-experimental.h"
#endif
#include    "main-sys.h"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef _LOCALPARALLEL_
  bool wait_subroutines(IET_Param * sys, int timeup_ms = -1);
#endif
#include    "main-diis.cpp"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include    "main-array.h"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include    "main-compress.cpp"
#ifdef _EXPERIMENTAL_
    #include "main-experimental.cpp"
#endif
#include    "main-analysis-param.cpp"
#include    "main-preprocessing.cpp"
#include    "main-build-ff.cpp"
#include    "main-hi.cpp"
#include    "main-rism.cpp"
#include    "main-mp.cpp"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include    "main-sub.cpp"
#include    "main-command.cpp"
IET_Param  * sys = nullptr;
IET_arrays * arr = nullptr;
RDF_data   * rdf = nullptr; int n_rdf_pairs = 0;    // RDF of current frame
RDF_data   * rdfs = nullptr; int n_rdfs_pairs = 0;  // RDF of all frames
#ifdef _INTERACTIVE_
  #include    "main-interactive.cpp"
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void init_software_constants(){
    for (int i=0; i<sizeof(CLOSURE_name)/sizeof(CLOSURE_name[0]); i++){ CLOSURE_name[i] = "";
        for (int j=0; j<sizeof(CLOSURE_alias)/sizeof(CLOSURE_alias[0]); j++) if (i==CLOSURE_alias[j].id){
            CLOSURE_name[i] = CLOSURE_alias[j].name; break;
        }
    }
    #ifdef _EXPERIMENTAL_
        init_experimental_contants();
    #endif
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------------------   main procedure   ------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
int main(int argc, char * argv[]){
    bool success = true; int error = 0; init_timer();
    bool help_out = false;
    init_software_constants();
  // determine the debug level, detail level, log filename and check -h
    int debug_level = -1; int detail_level = -1;
    pre_analysis_params(argc, argv, &debug_level, &detail_level, szfn_log, &help_out);
  // initialization of basic data structures
    memset(hostname, 0, sizeof(hostname)); gethostname(hostname, sizeof(hostname)-1);
    memset(username, 0, sizeof(username)); getlogin_r(username, sizeof(username)-1);
    sys = (IET_Param*) memalloc(sizeof(IET_Param)); if (sys) sys->init(argc, argv); else success = false;
      if (sys && debug_level>=0) sys->debug_level = debug_level;
      if (sys && detail_level>=0) sys->detail_level = detail_level;
      if (sys) sys->nt = maximum_default_processors();
    arr = (IET_arrays*) memalloc(sizeof(IET_arrays)); if (!arr) success = false;
    szfn_path[0] = 0; char * p_szfn_path = getcwd(szfn_path, sizeof(szfn_path));
    sys->library_path = getenv("IETLIB"); //printf("ietalIB=%s\n", sys->library_path);

  // prepare the log file
    FILE * flog = stdout; if (success){
        if (!szfn_log[0]){
        } else if (StringNS::string(szfn_log)=="screen" || StringNS::string(szfn_log)=="con" || StringNS::string(szfn_log)=="stdout"){ flog = stdout;
        } else if (StringNS::string(szfn_log)=="stderr"){ flog = stderr;
        } else { flog = fopen(szfn_log, "a"); if (!flog) flog = stderr;
        }
    }
    sys->flog = flog; sys->is_log_tty = isatty(fileno(flog));
    if (szfn_log[0] && (StringNS::string(szfn_log)=="screen" || StringNS::string(szfn_log)=="con")) sys->is_log_tty = false;
  // print the startup information
    if (success && flog && !sys->listonly && !help_out) main_print_header(sys, arr, flog, argc, argv);
  // initialize paralleling
  #ifdef _FFTWMPPARALLEL_
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: fftw_init_threads()\n");
    fftw_init_threads();
  #endif
  #ifdef _LOCALPARALLEL_
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: init_mp()\n");
    init_mp();
  #endif
  // analysis parameters
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: analysis_params()\n");
    error = analysis_params(sys, argc, argv); if (error) return error;
  // analysis parameters post
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: analysis_params_post()\n");
    error = analysis_params_post(sys); if (error){ if (flog) fclose(flog); return error; }
  // set priority
    if (sys->nice_level>0){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: set priority to %d\n", sys->nice_level);
        if (setpriority(PRIO_PROCESS, getpid(), sys->nice_level)!=0) fprintf(flog, "%s%s : warning : fail to change nice level%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->is_log_tty?color_string_end:"");
    }
  // read topology: all solute atom sites
    if (sys->debug_level>=1) fprintf(sys->log(), "debug:: read_solute_ff(%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", szfn_solute, sys->is_log_tty?prompt_path_suffix:"\"");
    if (read_solute_ff(sys, szfn_solute) <= 0){ success = false; return error; }
        //for (int i=0; i<sys->nas; i++){ printf("\33[37m Solute Atom %d%s:%d%s (%g, %g, %g, %g) has bond ", sys->as[i].iaa, sys->as[i].mole, sys->as[i].index, sys->as[i].name, sys->as[i].mass, sys->as[i].charge, sys->as[i].sigma, sys->as[i].epsilon); for (int j=0; j<sys->as[i].nbond; j++) printf("%d ", sys->as[i].ibond[j]); printf("\n\33[0m"); }
    if (success){
        sys->traj.count = sys->nas; sys->traj.box = Vector(0,0,0); sys->traj.atom = (PDBAtom*) memalloc(sizeof(PDBAtom) * (sys->nas+1));
    }
    PDBAtom * ia = sys->traj.atom;
    if (success){ // copy molecule and atom names to PDBAtom []
        for (int i=0; i<sys->traj.count; i++){ int len = 0;
            len = strlen(sys->as[i].mole); if (len>sizeof(ia[i].mole)-1) len = sizeof(ia[i].mole)-1; memcpy(ia[i].mole, sys->as[i].mole, len);
            len = strlen(sys->as[i].name); if (len>sizeof(ia[i].name)-1) len = sizeof(ia[i].name)-1; memcpy(ia[i].name, sys->as[i].name, len);
            ia[i].index = i+1;
        } //for (int i=0; i<sys->traj.count; i++){ char buf[120]; memset(buf, 0, sizeof(buf)); ia[i].export_pdb(buf); printf("%s", buf); }
    }
  lap_timer_analysis_param();
  // prepare RDF pairs
    if (sys->n_rdf_grps>0){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: rdf_grps_to_pairs()\n");
        main_prepair_rdf_grps_to_pairs(sys, arr, &rdf, &n_rdf_pairs);
        main_prepair_rdf_grps_to_pairs(sys, arr, &rdfs, &n_rdfs_pairs);
    }
  lap_timer_analysis_param();
  // checking and other preparation
    if (sys->nv>MAX_SOL){ fprintf(flog, "%s%s : error : too many solvent sites (%d > %d)%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nv, MAX_SOL, sys->is_log_tty?color_string_end:""); success = false; }
  // prepare memory and grids
    if (success){
        if (sys->nr[0]==0 && sys->nr[1]==0 && sys->nr[2]==0){
            fprintf(sys->log(), "%s%s : error : grid number (-nr) not defined%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->is_log_tty?color_string_end:""); success = false;
        } else if (sys->nr[0]<2 || sys->nr[1]<2 || sys->nr[2]<2){
            fprintf(sys->log(), "%s%s : error : grid number (-nr %dx%dx%d) too small%s\n", sys->is_log_tty?color_string_of_error:"", software_name, sys->nr[0], sys->nr[1], sys->nr[2], sys->is_log_tty?color_string_end:""); success = false;
        } else if (sys->nr[0]<10 || sys->nr[1]<10 || sys->nr[2]<10){
            fprintf(sys->log(), "%s%s : warning : is the grid number (-nr %dx%dx%d) too small?%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->nr[0], sys->nr[1], sys->nr[2], sys->is_log_tty?color_string_end:"");
        }
    }
    if (success){
        if (sys->debug_level>=1) fprintf(sys->log(), "debug:: arr->alloc(nv=%d, nvm=%d, nx=%d, ny=%d, nz=%d, n_solute=%d)\n", sys->nv, sys->nvm, sys->nr[0], sys->nr[1], sys->nr[2], sys->traj.count);
        arr->alloc(sys, sys->nv, sys->nvm, sys->nr[0], sys->nr[1], sys->nr[2], sys->traj.count);
        if (sys->debug_level>=2){ char buffer[64];
            fprintf(sys->log(), "DEBUG:: memory allocation: %d blocks, %s\n", _memory_blk_total, print_memory_value(buffer, sizeof(buffer), _memory_total));
        }
    }

  lap_timer_alloc_memory();
  // load, generate, and refine the solvent-solvent correlations: gvv, wvv, nhkvv, zetavv
    if (success && sys->gvv_specification!=0){
     // generate correlation map
        if (!sys->gvv_map && sys->nv>0){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: generate_default_gvv_map()\n");
            main_generate_default_gvv_map(sys, arr, flog);
        }
    }
    if (success){
      // load solvent self correlations, and generate K space correlations
        bool b_read_solvent_xvv = true;
        if (sys->gvv_specification!=0){
            if (sys->debug_level>=1) fprintf(sys->log(), "debug:: read_solvent_gvv(%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", szfn_gvv, sys->is_log_tty?prompt_path_suffix:"\"");
            b_read_solvent_xvv &= read_solvent_gvv(sys, arr, szfn_gvv, sys->nv, sys->nv);
            if (b_read_solvent_xvv && sys->debug_level>=2){
                if (sys->xvv_extend<=1) fprintf(sys->log(), "DEBUG:: generate_solvent_xvv(nv=%d, n_gvv=%d)\n", sys->nv, arr->n_gvv); else fprintf(sys->log(), "DEBUG:: generate_solvent_xvv(nv=%d, n_gvv=%gx%g)\n", sys->nv, arr->n_gvv/sys->xvv_extend, sys->xvv_extend);
            }
            if (b_read_solvent_xvv) b_read_solvent_xvv &= generate_solvent_xvv(sys, arr, szfn_gvv, sys->nv, sys->nv, arr->n_gvv, sys->xvv_extend);
        }
        if (sys->b_save_original_xvv && (arr->wvv && arr->nhkvv)){
            backup_solvent_xvv_with_shift_only(sys, arr);
        }

        bool b_read_solvent_zeta = true; bool b_analysis_post = true;
        if (sys->n_zeta_list>0 || szfn_zeta[0]){
            b_read_solvent_zeta = false; b_analysis_post = false;
            if (szfn_zeta[0]){
                if (sys->n_zeta_list>0) fprintf(sys->log(), "%s%s : warning : [zeta] section is ignored as -zeta is specified%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->is_log_tty?color_string_end:"");
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: read_solvent_zeta(%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", szfn_zeta, sys->is_log_tty?prompt_path_suffix:"\"");
                b_read_solvent_zeta = read_solvent_zeta(sys, arr, szfn_zeta, sys->nvm);
            } else {
                if (sys->debug_level>=1) fprintf(sys->log(), "debug:: generate_solvent_zeta(%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", szfn_zeta, sys->is_log_tty?prompt_path_suffix:"\"");
                b_read_solvent_zeta = generate_solvent_zeta(sys, arr);
            }
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: calculate_lse_ab_automatically()\n");
            b_analysis_post = calculate_lse_ab_automatically(sys, arr, b_read_solvent_xvv, b_read_solvent_zeta);
        }

        if (!b_read_solvent_xvv || !b_read_solvent_zeta || !b_analysis_post) success = false;
    }
    if (success && arr->nhkvv){
        perform_xvv_enhancement(sys, arr);
    }
    if (success && sys->debug_level>=3 && sys->listonly){
        debug_display_solvent_xvv(sys, arr, sys->nv);
    }
    if (success && sys->mode_test && arr->wvv && arr->nhkvv) debug_show_rism_xvv_matrix(sys, arr, sys->nv, sys->nvm);
  lap_timer_analysis_param();

  // create threads or forks
  #ifdef _LOCALPARALLEL_
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: create_subroutines(np=%d)\n", sys->nt);
    create_subroutines(sys, arr);
    if (sys->debug_level>=0){
        if (sys->nt<=1){
            fprintf(sys->log(), "%s : process %d (nice %d)\n", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]));
        } else if (sys->mp_by_fork){
            fprintf(sys->log(), "%s : process %d (nice %d), %d process%s:", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]), sys->nt, sys->nt>1?"es":"");
            for (int i=0; i<sys->nt; i++) fprintf(sys->log(), " %d", __forkpid[i]);
            fprintf(sys->log(), "\n");
        } else {
            fprintf(sys->log(), "%s : process %d (nice %d), %d thread%s\n", software_name, sys->pid, getpriority(PRIO_PROCESS, __forkpid[0]), sys->nt, sys->nt>1?"s":"");
        }
    }
  #else
    fprintf(sys->log(), "%s : process %d (nice %d)\n", software_name, sys->pid, getpriority(PRIO_PROCESS, sys->pid));
  #endif
  #ifdef _INTERACTIVE_
    if (success && sys->allow_interactive){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: allocate_interactive()\n");
        allocate_interactive();
    }
  #endif

  // open trajectory file
    if (success){
        if (szfn_xtc[0]){
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: prepare_traj_input(traj=%s%s%s)\n", sys->is_log_tty?prompt_path_prefix:"\"", szfn_xtc, sys->is_log_tty?prompt_path_suffix:"\"");
            success = sys->prepare_traj_input(szfn_xtc, flog);
        } else {
            fprintf(sys->log(), "%s%s : warning : trajectory (-f or -traj) not specified.%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->is_log_tty?color_string_end:"");
        }
    }
  lap_timer_io();

  // list only: display setups
    if (sys->listonly){
        bool islogtty = flog? isatty(fileno(flog)) : true;
        fprintf(flog, "%s# ==========================================================================%s\n", islogtty?prompt_comment_prefix:"", islogtty?prompt_comment_suffix:"");
        list_sys_files(sys, flog, (char*)"# ");
        fprintf(flog, "%s# ==========================================================================%s\n", islogtty?prompt_comment_prefix:"", islogtty?prompt_comment_suffix:"");
        sys->dump_text(flog, "  ", "");
        success = false;
    }

  // begin of run for no trajectory systems
    if (success && !szfn_xtc[0]){
        fprintf(flog, "%s%s : warning : trajectory (-f) not specified and nothing will be done.%s\n", sys->is_log_tty?color_string_of_warning:"", software_name, sys->is_log_tty?color_string_end:"");
        success = false;
    }

  // generate name of output file if not yet given
    if (!szfn_out[0]){
        generate_default_output_filename(sys, szfn_out, sizeof(szfn_out), szfn_xtc, sizeof(szfn_xtc));
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: generate_default_output_filename()=%s\n", szfn_out);
    }

    //if (sys->log()) fflush(sys->log());
  // read and do calculation of each frame
    if (success){
      // display trajectory, solute, solvent and force field information
        fprintf(flog, "%s trajectory: %s\n", software_name, szfn_xtc);
        main_print_solute_list(sys, arr, flog, "    ");
        main_print_solvent_list(sys, arr, flog, "    ");
        main_print_forcefield_info(sys, arr, flog, "    ");
    }

  // flush log before performing calculation
    if (sys->log()) fflush(sys->log());

    int frame = 0; int nframe = 0; double tframe_last = 0; TPAppendix tpa; memset(&tpa, 0, sizeof(tpa)); size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = sys->nv * N3;
    if (success) while (szfn_xtc[0]){
      // ====================================================================
      // ====================================================================
      // frame handling begins here =========================================
      // ====================================================================
      // ====================================================================
        arr->is_rdf_calculated = true; // block the calculation of RDF for skipped frames

      // read frame data
        if (sys->debug_level>=3) fprintf(sys->log(), "DEBUG:: read_frame(%d)\n", frame+1);
        int read_frame_ret = read_frame(&sys->traj, &frame, &tpa);
        if (read_frame_ret==0){ break;
        } else if (read_frame_ret==-1){
            fprintf(sys->log(), "%s%s : error : cannot read trajectory %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(szfn_xtc), sys->is_log_tty?color_string_end:""); break;
        } else if (read_frame_ret==-2){
            fprintf(sys->log(), "%s%s : error : incorrect number of atoms: %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(szfn_xtc), sys->is_log_tty?color_string_end:""); break;
        } else if (read_frame_ret<0){
            fprintf(sys->log(), "%s%s : error : fail for %s%s\n", sys->is_log_tty?color_string_of_error:"", software_name, get_second_fn(szfn_xtc), sys->is_log_tty?color_string_end:""); break;
        }
      // frame control: skip unwanted frames
        if (true){
            if (sys->time_begin>0 && tpa.time<sys->time_begin) continue;
            if (sys->time_end>0 && tpa.time>sys->time_end) break;
            if (sys->time_step>0 && tpa.time-tframe_last<sys->time_step) continue;
        };
        nframe ++;
        //if (sys->debug_level>=5){ for (int i=0; i<sys->nas; i++) fprintf(sys->log(), "DEBUG:::: %5s %5s : %12g %12g %12g\n", sys->as[i].mole, sys->as[i].name, sys->traj.atom[i].r.x, sys->traj.atom[i].r.y, sys->traj.atom[i].r.z);}
      // display frame information
        double drg[3]; drg[0] = sys->traj.box.x / sys->nr[0]; drg[1] = sys->traj.box.y / sys->nr[1]; drg[2] = sys->traj.box.z / sys->nr[2];
        if (sys->mode_test){
            if (sys->handling_xtc) fprintf(flog, "> testing: %.2f ps,", tpa.time); else fprintf(flog, "> testing: Frame %g,", tpa.time);
            fprintf(flog, " box=%gx%gx%g nm³, grid=%gx%gx%g nm³\n", sys->traj.box.x, sys->traj.box.y, sys->traj.box.z,
                    sys->traj.box.x/sys->nr[0], sys->traj.box.y/sys->nr[1], sys->traj.box.z/sys->nr[2]);
            continue;
        }
        if (flog){
            if (sys->handling_xtc) fprintf(flog, "> %g ps,", tpa.time); else fprintf(flog, "> Frame %g,", tpa.time);
            fprintf(flog, " box=%gx%gx%g nm³, pbc=%s%s%s, grid=%gx%gx%g nm³\n", sys->traj.box.x, sys->traj.box.y, sys->traj.box.z,
                (!sys->pbc_x&&!sys->pbc_y&&!sys->pbc_x)? "none" : sys->pbc_x?"x":"", sys->pbc_y?"y":"", sys->pbc_z?"z":"",
                sys->traj.box.x/sys->nr[0], sys->traj.box.y/sys->nr[1], sys->traj.box.z/sys->nr[2]);
        }
        tframe_last = tpa.time;
      lap_timer_io();

      // initialize buffers for calculation
        arr->reset_for_calculation();
      lap_timer_io();

      // set box and vecor scales, recalculating all dr and dk of 3D grids
        if (flog && !arr->set_scales(sys, sys->traj.box)){
            fprintf(flog, "%s%s : error : unable to init box, check wvv, nhkvv and frame box vector: %s", sys->is_log_tty?color_string_of_error:"", software_name, szfn_xtc);
            if (sys->handling_xtc) fprintf(flog, " %g ps%s\n", tpa.time, sys->is_log_tty?color_string_end:""); else fprintf(flog, " frame %d%s\n", nframe, sys->is_log_tty?color_string_end:"");
            continue;
        }

      // initialize flag variables for command processing
        arr->is_energy_calculated = false; arr->is_rdf_calculated = false;

      // run first frame commands for the first frame (commands with @begin)
        if (nframe<=1) if (process_command_sequence(1, sys, arr, rdf,rdfs, &file_out, nframe, tpa.time)<0) break;

      // run commands for each frame
        if (process_command_sequence(0, sys, arr, rdf,rdfs, &file_out, nframe, tpa.time)<0) break;

      // run implied commands
        if (sys->cmd_flag_rdf_ever_display){
            if (!arr->is_rdf_calculated){
                double rcutoff = sys->rvdw>sys->rcoul?sys->rvdw:sys->rcoul;
                if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: cmd[end-of-frame] = calculate_rdf(rc=%g, bins=%d)\n", rcutoff, sys->out_rdf_bins);
                calculate_rdf(sys, arr, rdf, rcutoff, sys->out_rdf_bins);
            }
            if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: cmd[end-of-frame] = add_rdf_to_sum(bins=%d)\n", sys->out_rdf_bins);
            add_rdf_to_sum(sys, arr, rdf, rdfs, sys->out_rdf_bins);
        }

      // ====================================================================
      // ====================================================================
      // frame handling ends here ===========================================
      // ====================================================================
      // ====================================================================
        if (sys->log()) fflush(sys->log());
    }

  // end of run: process the post commands (commands with @end)
    if (success){
        if (process_command_sequence(-1, sys, arr, rdf,rdfs, &file_out, nframe, tpa.time)<0) success = false;

        if (!sys->mode_test && !sys->cmd_flag_energy_ever_display && sys->cmd_flag_rism_ever_performed){
            if (!arr->is_energy_calculated){ recalculate_energy(sys, arr); arr->is_energy_calculated = true; }
            arr->display_solvation_energy_full(sys, flog, nullptr);
        }
    }

  // end all subroutines
  #ifdef _LOCALPARALLEL_
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: send stop signal to subroutines\n");
    for (int i=1; i<MAX_THREADS; i++) sys->mp_tasks[i] = MPTASK_TERMINATE;
    //wait_subroutines_end(sys);
    wait_subroutines(sys);
    if (sys->debug_level>=2){ char buffer[1024];
        fprintf(sys->log(), "DEBUG:: %s done in %s sec\n", sys->mp_by_fork?"process":"thread[0]", display_time(__total_timer, buffer));
    }
  #endif

  // display total CPU and memory usages
    lap_timer_others();
    if (flog){
        if (!sys->listonly && (sys->detail_level>=2 || sys->debug_level>=1)){
            //if (isatty(fileno(flog))) fprintf(flog, "\33[37m");
            lap_display_timers(flog);
            //display_memory_cost(flog, _memory_blk_total, _memory_total);
            //if (isatty(fileno(flog))) fprintf(flog, "\33[0m");
        }
        display_memory_cost(flog, _memory_blk_total, _memory_total);
    }
    //bool is_stderr_tty = isatty(fileno(stderr));
    if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: mem_dispose_all()\n");

  // display the end information
    if (success && flog) main_print_tailer(sys, arr);

  // dispose all
    mem_dispose_all(); lap_timer_alloc_memory();
    if (flog && flog!=stdout && flog!=stderr){ FILE * flog_close = flog; flog = nullptr; fclose(flog_close); }
    if (file_out && file_out!=stdout && file_out!=stderr) fclose(file_out);
    if (!success) return 1; return 0;
}
