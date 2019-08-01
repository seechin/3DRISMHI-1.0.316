const char * software_name = "gmxtop2solute";
const char * software_version = "0.239.1455";
const char * copyright_string = "(c) Cao Siqin";

#include    <errno.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <stdint.h>
#include    <string.h>
#include    <math.h>
#include    <signal.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <time.h>
#include    <sys/time.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>

#define     real    double

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define PI  3.1415926535897932384626433832795
#include    "header.h"
#include    "String2.cpp"
#include    "matrix.cpp"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int software_argc = 0;
char ** software_argv = nullptr;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*bool is_string_number(StringNS::string str){
    char * stop = nullptr;
    real v = strtod(str.text, &stop);
    if (stop && stop-str.text<str.length) return false; else return true;
}
int analysis_line(StringNS::string sline, StringNS::string * sl, int maxnw, bool separate_szstr = false){
    int nw = 0; int idx = 0;
    while (nw<maxnw){ sl[nw] = StringNS::seek_first_word(sline, &idx); if (sl[nw].length<=0) break; nw++; }
    if (separate_szstr) for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = 0;
    return nw;
}*/
int * init_int_vector(int n){ int * p = (int*) malloc(sizeof(int) * n); memset(p, 0, sizeof(int) * n); return p; }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// a vector
real * init_vector(int n){ real * p = (real*) malloc(sizeof(real) * n); memset(p, 0, sizeof(real) * n); return p; }
// a matrix
real ** init_matrix(int m, int n, int overflow_chars=0){
    int lenh = sizeof(real*) * m; int len = lenh + sizeof(real) * m * n + overflow_chars;
    char * p = (char*) malloc(len); memset(p, 0, len); real * d = (real*)(p + lenh);
    real ** a = (real**) p; for (int i=0; i<m; i++) a[i] = &d[i*n];
    return a;
}
void cp_matrix(real ** src, real ** dst, int m, int n){
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) dst[i][j] = src[i][j];
}
// a 3D tensor
real *** init_tensor3d(int nz, int ny, int nx, int overflow_chars=0){
    int lenhz = sizeof(real**) * nz; int lenhy = sizeof(real*) * nz*ny ;
    int len = lenhz + lenhy + sizeof(real) * nx * ny * nz + overflow_chars;
    char * p = (char*) malloc(len); memset(p, 0, len); real * d = (real*)(p + lenhz + lenhy);
    real *** a = (real***) p; real ** b = (real**) (p + lenhz);
    for (int i=0; i<nz; i++) a[i] = &b[i*ny];
    for (int i=0; i<nz; i++) for (int j=0; j<ny; j++) a[i][j] = &d[j*nx + i*nx*ny];
    return a;
}
void cp_tensor3d(real *** src, real *** dst, int nz, int ny, int nx){
    for (int i=0; i<nz; i++) for (int j=0; j<ny; j++) for (int k=0; k<nx; k++) dst[i][j][k] = src[i][j][k];
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define CLOSURE_D2      1
#define CLOSURE_HNC     2
#define CLOSURE_D2HNC   3
#define CLOSURE_KH      4
#define CLOSURE_D2KH    5
#define CLOSURE_PY      6
#define CLOSURE_PSE2    7
#define CLOSURE_P2E     8
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define MAX_SOL         100
#define MAX_PATH        512
#define MAX_PATHS       23
#define PRE_PATHS       3
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
const char * szHelp = "\
  Tool for lavender: translate gromacs topology to solute parameters\n\
  The input/output files:\n\
    -p, -top              topology file, TOP\n\
    -ffpath, -include     forcefield folder, multiple separated with \":\"\n\
    -o                    output file, default: screen\n\
    -debug                show debug info\n\
    -excl                 exclude group, default: SOL\n\
    -use-atom-name        (-an) use atom name, not atom type (default)\n\
    -use-atom-type        (-at) use atom type, not atom name\n\
    -solvent-format       (-for-gensolvent) output as solvent format\n\
    -bond, -no-bond[s]    show/hide bond information, default no bond\n\
    -no-index             (-ni) don't output atom and molecule index\n\
    -abbreviate           (-ab) allow #repeat commands, default off\n\
    -old-format           = -abbreviate -no-bond -no-index\n\
    -default              reset all options\n\
  The output format:\n\
    mole_name atom_name mass charge sigma epsilon\n\
";
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool debug = false; bool show_atom_spc_name = true; bool solvent_format = false;
    bool abbreviate_format = false; bool allow_bond = false; bool allow_index = true;
char * info_file_name = (char*)"";
char szfn_top[MAX_PATH];
char szfn_ffgmxt[PRE_PATHS][MAX_PATH];
char szfn_out[MAX_PATH];
char * szfn_ffpaths[MAX_PATHS]; int nszfn_ffpath = 0; char szfn_ffpath[MAX_PATH]; // search path for topology
char excl_grp[MAX_PATH];
int analysis_parameter_line(char * argv[], int * argi, int argc, char * script_name, int script_line){
    int ret = 0; int i = *argi; bool analysis_script = !script_name? false : (!script_name[0]? false : true);
    StringNS::string key = argv[i];
    if (!analysis_script && (key == "-h" || key == "-help" || key == "--h" || key == "--help")){ ret = 2;
    } else if (key == "-p" || key == "--p" || key == "-top" || key == "--top"){ if (i+1<argc){ i++; strcpy(szfn_top, argv[i]); }
    } else if (key == "-o" || key == "--o"){ if (i+1<argc){ i++; strcpy(szfn_out, argv[i]); }
    } else if (key == "-excl"){ if (i+1<argc){ i++; strcpy(excl_grp, argv[i]); }
    } else if (key == "-debug"){ debug = true;
    } else if (key=="-ffpath" || key=="--ffpath" || key=="-ffpath" || key=="-include" || key=="--include" || key=="include"){ if (i+1<argc){
        i++; strcpy(szfn_ffpath,argv[i]); for (int j=0; j<MAX_PATH && szfn_ffpath[j]; j++) if (szfn_ffpath[j]==':') szfn_ffpath[j] = 0;
      }
    } else if (key=="-an"||key=="--an"||key=="-use-atom-name"||key=="--use-atom-name"||key=="-use_atom_name"||key=="--use_atom_name"){
        show_atom_spc_name = true;
    } else if (key=="-at"||key=="--at"||key=="-use-atom-type"||key=="--use-atom-type"||key=="-use_atom_type"||key=="--use_atom_type"){
        show_atom_spc_name = false;
    } else if (key=="-solvent-format"||key=="--solvent-format"||key=="-solvent_format"||key=="--solvent_format"||key=="-for-gensolvent"||key=="--solvent-format"||key=="-for_gensolvent"||key=="--solvent_format"){
        solvent_format = true;
    } else if (key=="-ab"||key=="--ab"||key=="-abbreviate"||key=="--abbreviate"){
        abbreviate_format = true;
    } else if (key=="-nb"||key=="--nb"||key=="-no-bond"||key=="--no-bond"||key=="-no-bonds"||key=="--no-bonds"){
        allow_bond = false;
    } else if (key=="-bond"||key=="--bond"||key=="-bond"){
        allow_bond = true;
    } else if (key=="-ni"||key=="--ni"||key=="-no-index"||key=="--no-index"){
        allow_index = false;
    } else if (key=="-old"||key=="--old"||key=="-old-format"||key=="--old-format"){
        abbreviate_format = true; allow_bond = allow_index = false;
    } else if (key=="-default"||key=="--default"||key=="-default-format"||key=="--default-format"){
        show_atom_spc_name = true; solvent_format = false;
        abbreviate_format = false; allow_bond = false; allow_index = true;
    } else {
        strcpy(szfn_top, argv[i]);
    }
    *argi = i;
    return ret;
}
int analysis_params(int argc, char * argv[]){
    bool success = true; int error = 0;
    if (success){ szfn_ffgmxt[0][0] = szfn_ffgmxt[1][0] = szfn_ffgmxt[2][0] = 0;
        char * sz_env_ffgmxt = getenv ((char*)"GMXDATA");
        if (sz_env_ffgmxt){
            strncpy(szfn_ffgmxt[0], sz_env_ffgmxt, sizeof(szfn_ffgmxt[0]));
            strncpy(szfn_ffgmxt[1], sz_env_ffgmxt, sizeof(szfn_ffgmxt[1])-4); strncat(szfn_ffgmxt[1], "/top", MAX_PATH-1-strlen(szfn_ffgmxt[1]));
            strncpy(szfn_ffgmxt[2], sz_env_ffgmxt, sizeof(szfn_ffgmxt[2])-12); strncat(szfn_ffgmxt[2], "/gromacs/top", MAX_PATH-1-strlen(szfn_ffgmxt[2]));
        }
    }
  // analysis command line params
    if (argc<2){ printf("%s %s\n", software_name, software_version); return 0; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '#' || argv[i][0] == ';'){
        } else {
            int _error = analysis_parameter_line(argv, &i, argc, (char*)"", i);
            if (_error){ success = false; error |= _error; }
        }
    }
    if (error == 2){
        printf("%s %s %s\n", software_name, software_version, copyright_string);
        printf("%s", szHelp);
    }
    if (!success) return error;
  // prepare other params
    if (success){   // prepare the search path
        nszfn_ffpath = 1; for (int j=0; j<MAX_PATHS; j++) szfn_ffpaths[j] = nullptr; {
            for (int j=0; j<MAX_PATH && nszfn_ffpath<MAX_PATHS-PRE_PATHS;){
                if (szfn_ffpath[j]) szfn_ffpaths[nszfn_ffpath++] = &szfn_ffpath[j];
                j += strlen(&szfn_ffpath[j]) + 1; if (j>=MAX_PATH || szfn_ffpath[j]==0) break;
            }
        }; szfn_ffpaths[0] = (char*)"";
        if (szfn_ffgmxt[0][0]) szfn_ffpaths[nszfn_ffpath++] = szfn_ffgmxt[0];
        if (szfn_ffgmxt[1][0]) szfn_ffpaths[nszfn_ffpath++] = szfn_ffgmxt[1];
        if (szfn_ffgmxt[2][0]) szfn_ffpaths[nszfn_ffpath++] = szfn_ffgmxt[2];
      } // for (int j=0; j<nszfn_ffpath; j++) printf("include path: [%s]\n", szfn_ffpaths[j]);

    return 0;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class AtomType { public:
    char name[32], mole[32]; double mass; double charge; double sigma; double epsilon;
    void init(char * _name, char * _mole, double _mass, double _charge, double _sigma, double _epsilon){ int len = 0;
        memset(name, 0, sizeof(name)); len = (int)strlen(_name); if (len>31) len = 31; memcpy(name, _name, len);
        memset(mole, 0, sizeof(mole)); len = (int)strlen(_mole); if (len>31) len = 31; memcpy(mole, _mole, len);
        mass = _mass; charge = _charge; sigma = _sigma; epsilon = _epsilon;
        //printf("init %12s : %12f %12f %12f\n", name, charge, sigma, epsilon);
    }
};
AtomType * at = nullptr; int nat = 0; int natmax = 50000;
#define MAX_BONDS_PER_ATOM 12
class AtomTypeX : public AtomType { public:
    int index; int iaa; int nb; int ib[MAX_BONDS_PER_ATOM];
    AtomTypeX * next;
    void init(int _index, char * _name, char * _mole, double _mass, double _charge, double _sigma, double _epsilon){
        index = _index; nb = 0; iaa = 1;
        AtomType::init(_name, _mole, _mass, _charge, _sigma, _epsilon); next = nullptr;
    }
};
class MoleType { public:
    char name[32]; AtomTypeX * ar;
    void init(char * _name){
        memset(name, 0, sizeof(name)); int len = (int)strlen(_name); if (len>31) len = 31;
        memcpy(name, _name, len); ar = nullptr;
    }
};
MoleType * mt = nullptr; int nmt = 0; int nmtmax = 5000;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define MAX_STR 50
#define MAX_RECURS 50
FILE * fout = stdout;
int nrecursive = 0; StringNS::string sl[MAX_STR]; char input[4096];
bool analysis_top(char * filename, char * last_file_name, int last_line){
    bool success = true;
    bool ret = true; FILE * file = nullptr; int nline = 0; char fn[MAX_PATH];
    if (!file){
        strcpy(fn, last_file_name);
        for (int i=(int)strlen(fn)-1; i>=0; i--) if (i==0 || fn[i]=='/'){ fn[i+1] = 0; if (i==0) fn[i] = 0; break; }
        strncat(fn, filename, MAX_PATH-1-strlen(fn)); file = fopen(fn, "r");
    }
    if (!file) for (int i=1; i<nszfn_ffpath; i++){
        if (szfn_ffpaths[i][0]) { strcpy(fn, szfn_ffpaths[i]); strncat(fn, "/", MAX_PATH-1-strlen(fn)); } else fn[0] = 0;
        strncat(fn, filename, MAX_PATH-1-strlen(fn)); file = fopen(fn, "r"); if (file) break;
    }
    if (!file){ fprintf(stderr, "%s : %s[%d] : cannot open \"%s\"\n", software_name, last_file_name, last_line, filename); return false; }
    if (debug) fprintf(stderr, "%s : debug : handling %s\n", software_name, fn);

    int on_compile = 0; int imolnow = 0; int iindex = 0; int iimol = 0; int iaa_base = -1; int iaa_now = 0;
    while (success && fgets(input, sizeof(input), file)){ nline++;
        for (int i=0; i<sizeof(input) && input[i]; i++) if (input[i] == '\r' || input[i] =='\n') input[i] = 0;
        int nw = analysis_line(input, sl, MAX_STR, true); if (nw<=0) continue ;
        if (sl[0] == "#include" && nw>1){
            nrecursive ++; //printf("\033[31m%s include file: %s\033[0m\n", fn, sl[1].text);
            if (nrecursive >= MAX_RECURS){
                fprintf(stderr, "%s : %s[%d] : recursive overflow when opening \"%s\"\n", software_name, last_file_name, last_line, filename);
            } else {
                ret &= analysis_top(sl[1].text, fn, nline);
            }
            nrecursive --; //printf("\033[32mend recursive handling\033[0m\n");
        } else if (sl[0].text[0] == ';' || sl[0].text[0] == '#' || sl[0].text[0] == '*'){
        } else if (sl[0].text[0] == '['){
            if (sl[0] != "["){ sl[1] = &sl[0].text[1]; if (nw<2) nw = 2; }
            if (nw < 2) continue;
            if (sl[1] == "atomtypes"){ on_compile = 1;
            } else if (sl[1] == "moleculetype"){ on_compile = 2;
            } else if (sl[1] == "atoms"){ on_compile = 3;
            } else if (sl[1] == "molecules"){ on_compile = 4;
                if (solvent_format) fprintf(fout, "[atom]\n"); else fprintf(fout, "[solute]\n");
                fprintf(fout, "# generated by:"); for (int i=0; i<software_argc; i++) fprintf(fout, " %s", software_argv[i]); fprintf(fout, "\n");
            } else if (sl[1] == "system"){ on_compile = 5;
            } else if (sl[1] == "bonds"){ on_compile = 6;
            } else { on_compile = false;
            }
//if (on_compile) printf("section: %s[%d]: %s\n", fn, nline, sl[1].text);
        } else {
            if (on_compile==1){ // atomtypes
                for (int i=0; i<nw; i++) if (sl[i][0]=='#' || sl[i][0]==';') nw = i;
                if (nw>=7){
                    int icol_mass = 3; int icol_charge = 4; int icol_sigma = 6; int icol_epsilon = 7;
                    if (nw<8 || sl[4]=="A" || is_string_number(sl[5])){
                        icol_mass = 2; icol_charge = 3; icol_sigma = 5; icol_epsilon = 6;
                    }
                    int ist = nat;
                    for (int i=0; i<nat; i++) if (sl[0].Compare(at[i].name) == 0){ ist = i; break; }
                    if (debug && ist!=nat) fprintf(stderr, "%s : warning : overriding atomtype : %s\n", software_name, sl[0].text);
                    if (ist == nat) nat++; //at[nat++] = (AtomType*) malloc(sizeof(AtomType));
                    at[ist].init(sl[0].text, (char*)"*", atof(sl[icol_mass].text), atof(sl[icol_charge].text), atof(sl[icol_sigma].text), atof(sl[icol_epsilon].text));
                } else {
                    fprintf(stderr, "%s : %s[%d] : syntex error in atom type %s\n", software_name, fn, nline, sl[0].text); success = false;
                }
                if (nat >= natmax){ fprintf(stderr, "%s : error : too many atomtypes (max: %d)\n", software_name, natmax); success = false; }
            } else if (on_compile==2){ // moleculetype
                imolnow = -1;
                for (int i=0; i<nmt; i++) if (sl[0] == mt[i].name){ imolnow = i; break; }
                if (imolnow>=0) fprintf(stderr, "%s : warning : overriding moletype : %s\n", software_name, sl[0].text);
                if (imolnow<0){ imolnow = nmt; nmt++; }
                mt[imolnow].init(sl[0].text);
              if (nmt >= nmtmax){ fprintf(stderr, "%s : error : too many moletype (max: %d)\n", software_name, nmtmax); success = false; }
            } else if (on_compile==3){ // atoms
                if (imolnow<0 || imolnow>=nmtmax){ fprintf(stderr, "%s : warning : ignore atom %s without molecules\n", software_name, sl[0].text);
                } else {
                    int idef = -1; for (int i=0; i<nat; i++) if (sl[1].Equ(at[i].name)){ idef = i; break; }
                    if (idef<0){
                        fprintf(stderr, "%s : %s[%d] : error : atom %s not defined\n", software_name, fn, nline, sl[1].text); success = false;
                    } else {
                        //printf("Atom %s define: %s\n", sl[1].text, at[idef].name);
                        AtomTypeX * p = (AtomTypeX*) malloc(sizeof(AtomTypeX));
                        p -> init(atoi(sl[0].text), show_atom_spc_name?sl[4].text:at[idef].name, sl[3].text, at[idef].mass, at[idef].charge, at[idef].sigma, at[idef].epsilon);
                        p->iaa = atof(sl[2].text);
                        if (nw>=6) p->charge = atof(sl[6].text);
                        if (nw>=7) p->mass = atof(sl[7].text);
//printf("    %-12s %12f %12f %12f %12f\n", p->name, p->mass, p->charge, p->sigma, p->epsilon); printf("    %-12s %12f %12f %12f %12f\n", at[idef].name, at[idef].mass, at[idef].charge, at[idef].sigma, at[idef].epsilon);
                        if (!mt[imolnow].ar) mt[imolnow].ar = p; else {
                            AtomTypeX * q = mt[imolnow].ar; while (q->next) q = q->next; q->next = p;
                        }
                    }
                }
            } else if (on_compile==6 ){ // bond
                for (int i=0; i<nw; i++) if (sl[i][0]=='#' || sl[i][0]==';') nw = i;
                if (nw<2){
                    fprintf(stderr, "%s : warning : incomplete bond line \"%s\" ignored\n", software_name, sl[0].text);
                } else if (!(StringNS::is_string_number(sl[0]) && StringNS::is_string_number(sl[1]))){
                    fprintf(stderr, "%s : warning : incorrect bond line \"%s %s\" ignored\n", software_name, sl[0].text, sl[1].text);
                } else if (imolnow<0 || imolnow>=nmtmax){
                    fprintf(stderr, "%s : warning : ignore bond %s-%s without molecules\n", software_name, sl[0].text, sl[1].text);
                } else {
                    int bondi = atoi(sl[0].text); int bondj = atoi(sl[1].text);
                    AtomTypeX * ai = nullptr; AtomTypeX * aj = nullptr;
                    for (AtomTypeX * ax = mt[imolnow].ar; ax && (!ai || !aj); ax=ax->next){
                        if (ax->index==bondi) ai = ax; if (ax->index==bondj) aj = ax;
                    }
                    if (ai && aj){
                        if (ai->nb+1<MAX_BONDS_PER_ATOM) ai->ib[ai->nb++] = aj->index - ai->index;
                        if (aj->nb+1<MAX_BONDS_PER_ATOM) aj->ib[aj->nb++] = ai->index - aj->index;
                    }
                    //printf("define bond for %d : %d (%s:%s) - %d (%s:%s)\n", imolnow, bondi, ai?ai->mole:"nullptr", ai?ai->name:"nullptr", bondj, aj?aj->mole:"nullptr", aj?aj->name:"nullptr");
                }
            } else if (on_compile==4){ // molecules
                if (sl[0] == excl_grp){
                    fprintf(stderr, "# exclude: %s %d\n", sl[0].text, atoi(sl[1].text));
                    continue;
                }
                if (nw<2){
                    fprintf(stderr, "%s : %s[%d] : error : incomplete mole line\n", software_name, fn, nline); success = false; continue;
                }
                int imol = -1; for (int i=0; i<nmt; i++) if (sl[0] == mt[i].name){ imol = i; break; }
                if (imol<0){
                    fprintf(stderr, "%s : %s[%d] : error : molecule %s undefined\n", software_name, fn, nline, sl[0].text); success = false; continue;
                }
                int nm = atoi(sl[1].text);
                /*
                for (int im=0; im<nm; im++){
                    int ia = 0; iimol++;
                    for (AtomTypeX * p = mt[imol].ar; p; p=p->next){
                        //fprintf(fout, "%5d %6s %5d %6s %3d %12f %12f %12f %12f\n", ++iindex, sl[0].text, iimol, p->name, ++ia, p->mass, p->charge, p->sigma, p->epsilon);
                        ++iindex;
                        fprintf(fout, "%4s %4s %5.2f %12g %12g %12g\n", p->name, show_atom_spc_name?p->mole:sl[0].text, p->mass, p->charge, p->sigma, p->epsilon);
                    }
                }
                */
                fprintf(fout, "# molecule: %s\n", mt[imol].name);
                int n_atom_in_mol = 0;

                int count_display = nm; int count_abbreviate = 0;
                if (abbreviate_format){ count_display = 1; count_abbreviate = nm-1; }

                for (int ec=0; ec<count_display; ec++){
                    for (AtomTypeX * p = mt[imol].ar; p; p=p->next){
                        //fprintf(fout, "%5d %6s %5d %6s %3d %12f %12f %12f %12f\n", ++iindex, sl[0].text, iimol, p->name, ++ia, p->mass, p->charge, p->sigma, p->epsilon);
                        ++iindex; n_atom_in_mol ++;
                        if (solvent_format){
                            fprintf(fout, "%4s %4s %3d %3d %12g %12g %12g\n", p->name, show_atom_spc_name?p->mole:sl[0].text, iindex, iindex, p->charge, p->sigma, p->epsilon);
                        } else {
                            //fprintf(fout, "%4s %4s %5.2f %12g %12g %12g\n", p->name, show_atom_spc_name?p->mole:sl[0].text, p->mass, p->charge, p->sigma, p->epsilon);
                            if (allow_index){
                                if (iaa_base<0) iaa_base = 0; else if (p->iaa+iaa_base < iaa_now) iaa_base = iaa_now;
                                iaa_now = p->iaa + iaa_base;
                                fprintf(fout, "%6d %4s %6d %4s %5.2f %12g %12g %12g", iindex, p->name, iaa_now, show_atom_spc_name?p->mole:sl[0].text, p->mass, p->charge, p->sigma, p->epsilon);
                            } else {
                                fprintf(fout, "%4s %4s %5.2f %12g %12g %12g", p->name, show_atom_spc_name?p->mole:sl[0].text, p->mass, p->charge, p->sigma, p->epsilon);
                            }
                            if (allow_bond && p->nb>0){
                                fprintf(fout, "  bond:");
                                for (int i=0; i<p->nb; i++)fprintf(fout, "%d%s", p->ib[i], i+1==p->nb?"":",");
                            }
                            fprintf(fout, "\n");
                        }
                    }
                }
                if (count_abbreviate>1){ fprintf(fout, "#repeat %d atoms %d times\n", n_atom_in_mol, count_abbreviate); iindex += n_atom_in_mol*(count_abbreviate); }

                iimol+=nm;

            } else if (on_compile==5){ // system
            }
        }
    }

    fclose(file); return true;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char * argv[]){
    software_argc = argc; software_argv = argv;
    bool success = true; int error = 0; strcpy(excl_grp, ""); szfn_out[0] = 0;

    error = analysis_params(argc, argv); if (error) return error;
    if (error) success = false;

    at = (AtomType*) malloc(sizeof(AtomType) * natmax);
    mt = (MoleType*) malloc(sizeof(MoleType) * nmtmax);

    if (StringNS::string(szfn_out)=="con" || StringNS::string(szfn_out)=="stdout" || StringNS::string(szfn_out)=="screen"){
        fout = stdout;
    } else if (StringNS::string(szfn_out)=="stderr"){
        fout = stderr;
    } else if (szfn_out[0]) {
        fout = fopen(szfn_out, "w");
        if (!fout){ fprintf(stderr, "%s : error : cannot write to %s\n", software_name, szfn_out); success = false; }
    }


  // read topology
    if (success && szfn_top[0]) analysis_top(szfn_top, (char*)"", 1);
  // checking and other preparation

//for (int i=0; i<nat; i++) printf("AtomType %5d: %12s: %12f %12f %12f %12f\n", i, at[i].name, at[i].mass, at[i].charge, at[i].sigma, at[i].epsilon);
//for (int i=0; i<nmt; i++){ printf("MoleType %5d: %12s\n", i, mt[i].name); for (AtomTypeX * p = mt[i].ar; p; p=p->next) printf("    %-12s %12f %12f %12f %12f\n", p->name, p->mass, p->charge, p->sigma, p->epsilon); }

    if (!success) return 1;

    return 0;
}
