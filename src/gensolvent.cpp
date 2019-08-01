const char * software_name = "gensolvent";
const char * software_version = "0.231.1400";
const char * copyright_string = "(c) Cao Siqin";

#define     __REAL__    double
#define     MAX_SOL     100     // Max atom site number

#include    "header.h"
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
#ifndef nullptr
  #define nullptr NULL
#endif

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define PI  3.1415926535897932384626433832795
#define EE  2.7182818284590452353602874713527
#define COULCOOEF 138.9354846
#define MAX_PATH        1024
#define MAX_WORD        200
#define MAX_ITER_TIMES  20
#define MAX_ITEMS       50
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
const char * szHelp = "\
    -p                  forcefield file starts with [atom]\n\
                        can be screen/con\n\
    -s                  topology file, can be: pdb/gro\n\
    -nice               nice level of process, default 0\n\
    -f, -traj           trajectory file: pdb/gro/xtc\n\
    -o                  output file, default: screen/con\n\
      gvv (in -p file)  output solvent-RDF file, default: screen/con\n\
      zvv (in -p file)  output solvent-zeta file, default: screen/con\n\
    -gvv, -zeta         files to output gvv or zeta\n\
    -b, -e, -dt         time of begin, end and step of reading traj\n\
    -ff opls/amber/gaff force field type indicator\n\
    -arith/geo-sigma    LJ sigma combination rule: arithmetic or geometric\n\
    -Lorentz-Berthelot  = -arith-sigma\n\
    -rc                 (also -rvdw, -rcoul) cutoff of rdf calculation\n\
    -dr, -bins          dr (bins=rc/dr) of calculationg rdf and zeta\n\
    -force-override     override output files if exist, default not set\n\
    -detail             show more details of running\n\
    -rename-atom        rename atoms by element_name+index\n\
    -rename-mole        rename molecules by substitution: old:new\n\
    -rename[-all]       rename atoms and molecules, see -rename-mole\n\
";
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include    "Element.h"
#include    "String2.cpp"
#include    "Vector.cpp"
#include    "PDBAtom.cpp"
#include    "read_pdb_gro.cpp"
#include    "read_frame_abr.cpp"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class AtomSite { public:
    char name[32]; char mole[32]; char nmele[32];
    double mass; double charge;
    double sigma, sqrt_sigma; double epsilon, sqrt_epsilon; int iaa;
    int id, grp, multi; bool is_key;
    void init(int _id, int _grp, const char * _mole, int _iaa, const char * _name, const char * _nmele, double _mass, double _charge, double _sigma, double _epsilon){
        strncpy(name, _name, sizeof(name));
        strncpy(mole, _mole, sizeof(mole));
        strncpy(nmele, _nmele, sizeof(nmele));
        iaa = _iaa; mass = _mass; charge = _charge; sigma = _sigma; epsilon = _epsilon;
        id = _id; grp = _grp; multi = 1; is_key = false;
        sqrt_sigma = sqrt(fabs(sigma)); sqrt_epsilon = sqrt(fabs(epsilon));
        //printf("atom site %s.%s: %12f %12f %12f %12f\n", mole, name, charge, sigma, epsilon, sqrt_epsilon);
    }
};
class AtomIndex { public: int index, iaa, grp; bool is_center; };
class MoleculeIndex { public: char * mole; int iaa, ibegin, iend, icenter; };
class GroupInfo { public: int grp, ia, multi, nmole; };
class NameSubstitution { public: char name[2][64]; };
template <class DT> DT *** init_tensor3d(size_t nz, size_t ny, size_t nx, size_t overflow_chars=0){
    size_t lenhz = sizeof(DT**) * nz; size_t lenhy = sizeof(DT*) * nz*ny ;
    size_t len = lenhz + lenhy + sizeof(DT) * nx * ny * nz + overflow_chars;
    //printf("ALLOCATING: %12d + %12d + %d\n", lenhz+lenhy, sizeof(DT) * nx * ny * nz, overflow_chars);
    //printf("allocating %12f\n", sizeof(DT) * nx * ny * nz/(double)len); fflush(stdout);
    char * p = (char*) malloc(len); memset(p, 0, len); DT * d = (DT*)(p + lenhz + lenhy);
    DT *** a = (DT***) p; DT ** b = (DT**) (p + lenhz);
    for (size_t i=0; i<nz; i++) a[i] = &b[i*ny];
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) a[i][j] = &d[j*nx + i*nx*ny];
    return a;
}
template <class DT> void cp_tensor3d(DT *** src, DT *** dst, size_t nz, size_t ny, size_t nx){
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) for (size_t k=0; k<nx; k++) dst[i][j][k] = src[i][j][k];
}
template <class DT> void clear_tensor3d(DT *** dst, size_t nz, size_t ny, size_t nx){
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) for (size_t k=0; k<nx; k++) dst[i][j][k] = 0;
}
template <class DT> void cp_tensor3d(DT *** src, DT *** dst, size_t n3){
    DT * d = &dst[0][0][0]; DT * s = & src[0][0][0];
    for (size_t i3=0; i3<n3; i3++) d[i3] = s[i3];
}
template <class DT> void clear_tensor3d(DT *** dst, size_t n3){
    DT * d = &dst[0][0][0];
    for (size_t i3=0; i3<n3; i3++) d[i3] = 0;
}
double get_current_time_double(){
    struct timeval tv; gettimeofday(&tv, nullptr); return ((tv.tv_sec) + (tv.tv_usec) / 1000000.0);
}
int ff_pbc_i(int x, int box_x){
    if (x<0) x += floor(abs(x/box_x)+1)*box_x;
    return x%box_x;
}
double ff_pbc_d(double x, double box_x){
    if (x<0) x += floor(fabs(x/box_x)+1)*box_x;
    return fmod(x, box_x);
}
Vector pbc_v_sub(Vector r1, Vector r2, Vector box){
    double r1x = ff_pbc_d(r1.x, box.x); //if (r1x<-box.x/2) r1x += box.x; if (r1x>box.x/2) r1x -= box.x;
    double r1y = ff_pbc_d(r1.y, box.y); //if (r1y<-box.y/2) r1y += box.y; if (r1y>box.y/2) r1y -= box.y;
    double r1z = ff_pbc_d(r1.y, box.y); //if (r1z<-box.z/2) r1z += box.z; if (r1z>box.z/2) r1z -= box.z;

    double r2x = ff_pbc_d(r2.x, box.x); //if (r2x<-box.x/2) r2x += box.x; if (r2x>box.x/2) r2x -= box.x;
    double r2y = ff_pbc_d(r2.y, box.y); //if (r2y<-box.y/2) r2y += box.y; if (r2y>box.y/2) r2y -= box.y;
    double r2z = ff_pbc_d(r2.y, box.y); //if (r2z<-box.z/2) r2z += box.z; if (r2z>box.z/2) r2z -= box.z;

    double rx = r1x - r2x; double ry = r1y - r2y; double rz = r1z - r2z;

    //if (rx<-box.x/2) rx += box.x; if (rx>box.x/2) rx -= box.x;
    //if (ry<-box.y/2) ry += box.y; if (ry>box.y/2) ry -= box.y;
    //if (rz<-box.z/2) rz += box.z; if (rz>box.z/2) rz -= box.z;

    //return Vector(rx, ry, rz);
    return Vector(ff_pbc_d(rx, box.x), ff_pbc_d(ry, box.y), ff_pbc_d(rz, box.z));
}
Vector pbc_rm_v(Vector r, Vector box){
    return Vector(ff_pbc_d(r.x, box.x), ff_pbc_d(r.y, box.y), ff_pbc_d(r.z, box.z));
    //if (r.x<0) r.x += box.x; if (r.x>box.x) r.x -= box.x;
    //if (r.y<0) r.y += box.y; if (r.y>box.y) r.y -= box.y;
    //if (r.z<0) r.z += box.z; if (r.z>box.z) r.z -= box.z;
    //return r;
}
Vector pbc_v_sub_partial(Vector r1, Vector r2, Vector box){
    double rx = r1.x - r2.x; double ry = r1.y - r2.y; double rz = r1.z - r2.z;

    if (rx<-box.x/2) rx += box.x; if (rx>box.x/2) rx -= box.x;
    if (ry<-box.y/2) ry += box.y; if (ry>box.y/2) ry -= box.y;
    if (rz<-box.z/2) rz += box.z; if (rz>box.z/2) rz -= box.z;

    return Vector(rx, ry, rz);
    //return Vector(ff_pbc_d(rx, box.x), ff_pbc_d(ry, box.y), ff_pbc_d(rz, box.z));
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool show_warning = false; bool show_debug = false; bool force_override = false; bool force_rename_atom_names = false; int nice_level = 0;
bool arith_sigma = false; int ff_specifier = 0; // amber, opls, gaff
double time_begin = 0; double time_end = 0; double time_step = 0;
double rvdw = 1; double rcoul = 1; double rcutoff = 0;
double drrism = 0; double drhi = 0; int bins = 0; bool use_bins = false; int bins_rism = 0; int bins_hi = 0;
AtomSite av[MAX_SOL]; int av_i_in_mol[MAX_SOL]; int nav = 0; int nmv = 0; int ngv = 0;
GroupInfo gpi[MAX_SOL];
PDBAtomSet pdb_atom_set; AtomIndex * ai = nullptr; MoleculeIndex * mi = nullptr; int nm = 0;
double inverse_V = 0; int nmole_mv[MAX_SOL]; double dipole_mv[MAX_SOL]; double n_dipole_mv[MAX_SOL];
char szfp[MAX_PATH] = { 0, 0, 0, 0 };
char szfs[MAX_PATH] = { 0, 0, 0, 0 }; StringNS::string str_ext_fs = "";
char szff[MAX_PATH] = { 0, 0, 0, 0 }; StringNS::string str_ext_ff = "";
char szfo[MAX_PATH] = { 0, 0, 0, 0 };
char szfo_gvv[MAX_PATH] = { 0, 0, 0, 0 };
char szfo_zvv[MAX_PATH] = { 0, 0, 0, 0 };
char path_cwd[MAX_PATH]; char path_top[MAX_PATH];
NameSubstitution nss[MAX_ITEMS]; int n_nss = 0;
StringNS::string file_extension(StringNS::string fn){
    int ie = 0;
    for (ie=fn.length-1; ie>=0; ie--) if (fn.text[ie]=='.') break;
    if (ie==0||fn.text[ie]!='.') return StringNS::string((char*)"");
    while (ie<fn.length && fn.text[ie]=='.') ie++;
    return fn.Substring(ie, fn.length-ie);
}
void get_file_path(char * fn, char * pn){
    pn[0] = 0; int ifn = strlen(fn);
    for (; ifn>0 && fn[ifn]!='/'; ifn--);
    if (ifn>0 && fn[ifn]=='/'){
        char pt[MAX_PATH]; memset(pt, 0, MAX_PATH);
        for (int i=0; i<ifn; i++) pt[i] = fn[i]; if (!realpath(pt, pn)) return ;
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool analysis_args(int argc, char *argv[], int first_argv){
    bool exit_on_error = false;
    bool exit_on_help = false;

    if (argc > first_argv){
        for (int i = first_argv; i < argc; i ++){
            if (argv[i][0] == '-'){
                StringNS::string key = argv[i];
                if (key == "-?"){
                    printf("%s %s\n", software_name, software_version);
                    exit_on_help = true;
                } else if (key == "-h" || key == "-help" || key == "--help"){
                    printf("%s %s %s\n", software_name, software_version, copyright_string);
                    printf("%s", szHelp);
                    exit_on_help = true;
                } else if (key == "-p" || key == "-top" || key == "--top"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szfp, argv[i], sizeof(szfp));
                    }
                } else if (key == "-s"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szfs, argv[i], sizeof(szfs)); str_ext_fs = file_extension(argv[i]);
                        if (str_ext_fs!="pdb" && str_ext_fs!="gro"){
                            fprintf(stderr, "%s : argv[%d] : error : unknown type of -s %s. Support pdb/gro only\n", software_name, i, str_ext_fs.length==0?"(nullptr)":str_ext_fs.text);
                        }
                    }
                } else if (key == "-f" || key == "-traj" || key == "--traj"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szff, argv[i], sizeof(szff)); str_ext_ff = file_extension(argv[i]);
                        if (str_ext_ff!="pdb" && str_ext_ff!="gro" && str_ext_ff!="xtc"){
                            fprintf(stderr, "%s : argv[%d] : error : unknown type of -f %s. Support pdb/gro/xtc only.\n", software_name, i, str_ext_ff.length==0?"(nullptr)":str_ext_ff.text);
                        }
                    }
                } else if (key == "-o"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szfo, argv[i], sizeof(szff));
                    }
                } else if (key == "-gvv" || key == "--gvv"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szfo_gvv, argv[i], sizeof(szfo_gvv));
                    }
                } else if (key == "-zeta" || key == "--zeta" || key == "-zkvv" || key == "--zkvv"){
                    if (i+1<argc && argv[i+1][0]!='-'){
                        i++; strncpy(szfo_zvv, argv[i], sizeof(szfo_zvv));
                    }
                } else if (key == "-nice" || key == "--nice"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; nice_level = atoi(argv[i]); }
                } else if (key == "-b" || key == "-begin" || key == "--begin"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; time_begin = atof(argv[i]); }
                } else if (key == "-e" || key == "-end" || key == "--end"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; time_end = atof(argv[i]); }
                } else if (key == "-dt" || key == "--dt"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; time_step = atof(argv[i]); }
                } else if (key == "-dr" || key == "--dr"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; drrism = drhi = atof(argv[i]); use_bins = false; }
                } else if (key == "-bins" || key == "--bins"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; bins = bins_rism = bins_hi = atoi(argv[i]); use_bins = true; }
                } else if (key == "-rc" || key == "--rc" || key == "-rcutoff" || key == "--rcutoff"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; rvdw = rcoul = atof(argv[i]); }
                } else if (key == "-rvdw" || key == "--rvdw"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; rvdw = atof(argv[i]); }
                } else if (key == "-rcoul" || key == "--rcoul"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; rcoul = atof(argv[i]); }
                } else if (key == "-silent" || key == "--silent"){
                    show_warning = false;
                } else if (key == "-brief" || key == "--brief"){
                    show_warning = false;
                } else if (key == "-detail" || key == "--detail" || key == "-detailed" || key == "--detailed"){
                    show_warning = true;
                } else if (key == "-debug" || key == "--debug" || key == "-debuging" || key == "--debuging" || key == "-debugging" || key == "--debugging"){
                    show_warning = true; show_debug = true;
                } else if (key == "-arith-sigma" || key == "--arith-sigma" || key == "-arith_sigma" || key == "--arith_sigma"){
                    arith_sigma = true; ff_specifier = 3;
                } else if (key == "-Lorentz-Berthelot" || key == "--Lorentz-Berthelot" || key == "-Lorentz_Berthelot" || key == "--Lorentz_Berthelot"){
                    arith_sigma = true; ff_specifier = 3;
                } else if (key == "-geo-sigma" || key == "--geo-sigma" || key == "-geo_sigma" || key == "--geo_sigma"){
                    arith_sigma = false; ff_specifier = 2;
                } else if (key == "-ff" || key == "--ff" || key == "--forcefield" || key == "--forcefield"){
                    if (i+1<argc && argv[i+1][0]!='-'){ i++; StringNS::string key2 = StringNS::string(argv[i]);
                        if (key2=="opls" || key2=="oplsaa"){
                            arith_sigma = false; ff_specifier = 2;
                        } else if (key2=="gaff"){
                            arith_sigma = true; ff_specifier = 3;
                        } else if (key2=="amber"){
                            arith_sigma = true; ff_specifier = 1;
                        } else {
                            fprintf(stderr, "%s : argv[%d] : unkown forcefield type: %s\n", software_name, i, argv[i]);
                            exit_on_error = true;
                        }
                    }
                } else if (key == "-ffopls" || key == "--ffopls" || key == "-ffoplsaa" || key == "--ffoplsaa" || key == "-oplsff" || key == "--oplsff" || key == "-oplsaaff" || key == "--oplsaaff"){
                    arith_sigma = false; ff_specifier = 2;
                } else if (key == "-ffamber" || key == "--ffamber" || key == "-amberff" || key == "--amberff"){
                    arith_sigma = true; ff_specifier = 1;
                } else if (key == "-ffgaff" || key == "--ffgaff" || key == "-gaff" || key == "--gaff"){
                    arith_sigma = true; ff_specifier = 3;
                } else if (key == "-force-override" || key == "--force-override"){
                    force_override = true;
                } else if (key == "-rename-atom" || key == "--rename-atom" || key == "-rename_atom" || key == "--rename_atom"){
                    force_rename_atom_names = true;
                } else if (key == "-rename-mole" || key == "--rename-mole" || key == "-rename_mole" || key == "--rename_atom" || key == "-rename-all" || key == "--rename-all" || key == "-rename_all" || key == "--rename_all" || key == "-rename" || key == "--rename"){
                    n_nss = 0;
                    while (i+1<argc && argv[i+1][0]!='-'){
                        int number_of_seps = 0; for (int j=0; argv[i+1][j]; j++) if (argv[i+1][j]==':') number_of_seps ++;
                        if (number_of_seps!=1) break;
                        i++; if (n_nss<MAX_ITEMS){ memset(&nss[n_nss], 0, sizeof(nss[n_nss]));
                            int j = 0; for (; j<sizeof(nss[n_nss].name[0])-1 && argv[i][j] && argv[i][j]!=':'; j++) nss[n_nss].name[0][j] = argv[i][j];
                            j++; for (int k=0; k<sizeof(nss[n_nss].name[1])-1 && argv[i][j]; j++) nss[n_nss].name[1][k++] = argv[i][j];
                            n_nss ++;
                        }
                    }
                    if (!(key == "-rename-mole" || key == "--rename-mole" || key == "-rename_mole" || key == "--rename_atom")){
                        force_rename_atom_names = true;
                    }
                    //for (int i=0; i<n_nss; i++) printf("substitution: %s -> %s\n", nss[i].name[0], nss[i].name[1]);
                } else { fprintf(stderr, "%s : argv[%d] : error : unkown option %s. \n", software_name, i, argv[i]); exit_on_error = true; }
            } else { fprintf(stderr, "%s : argv[%d] : error : unkown string %s. \n", software_name, i, argv[i]); exit_on_error = true; }
        }
    } else { printf("%s %s\n", software_name, software_version); exit_on_help = true; }
    return !exit_on_help && !exit_on_error;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int __iaanow = 0; int __index = 0;
bool analysis_solvent_atom_7(char * argv[], int * argi, int argc, char * script_name, int script_line){
  //if (argc<7){ fprintf(stderr, "%s : %s[%d] : incomplete atom line\n", software_name, script_name, script_line); return false; }
    if (nav >= MAX_SOL){ fprintf(stderr, "%s : %s[%d] : too many solvents\n", software_name, (script_name), script_line); return false; }
    int id = atoi(argv[2]); int grp = atoi(argv[3]);
    ElementNS::Element * element = ElementNS::get_element_parameter(ElementNS::get_atom_element(argv[0]));
    if (!element) fprintf(stderr, "%s : %s[%d] : warning : no element found for atom %s\n", software_name, (script_name), script_line, argv[0]);
    double mass = element? element->m : 0;
    if (nav>0 && StringNS::string(argv[1]) != av[nav-1].mole) __iaanow ++;
    av[nav].init(++__index, grp, argv[1], __iaanow, argv[0], element?element->name:"?", mass, atof(argv[4]), atof(argv[5]), atof(argv[6]));
    nav ++;
  // end
    return true;
}
bool analysis_solvent_atom(char * argv[], int * argi, int argc, char * script_name, int script_line){
    int nargs = 0;
    for (; nargs<argc; nargs++) if (argv[nargs][0]=='#' || argv[nargs][0]==';') break;
    if (nargs<7){ fprintf(stderr, "%s : %s[%d] : incomplete atom line\n", software_name, (script_name), script_line); return false; }
    if (nargs>=7){
        return analysis_solvent_atom_7(argv, argi, nargs, script_name, script_line);
    } else return false;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int on_compile = 0;
bool read_solvent(char * filename, int iter_time){
    bool success = true; FILE * file = nullptr;
    if (!filename || !filename[0]){ fprintf(stderr, "%s : error : top file not specified (see -p)\n", software_name); return false; }
    if (StringNS::string(filename)=="screen"||StringNS::string(filename)=="con"||StringNS::string(filename)=="stdin") file = stdin; else file = fopen(filename, "r");
    if (!file){ fprintf(stderr, "%s : error : cannot open -top %s\n", software_name, filename&&filename[0]?filename:"(nullptr)"); return false; }

    if (iter_time==0){ get_file_path(filename, path_top); int __chdir = chdir(path_top); }

    StringNS::string sl[MAX_WORD]; int iline = 0; char input[4096];
    while (fgets(input, sizeof(input), file)){ iline ++;
        int nw = StringNS::analysis_line(input, sl, MAX_WORD, true); if (nw<1) continue;
        if (sl[0].text[0] == '#'){
            if (sl[0]=="#include" && nw>1){
                if (iter_time+1>MAX_ITER_TIMES){ fprintf(stderr, "%s : %s[%d] : too many includes\n", software_name, filename, iline); success = false; }
                else success &= read_solvent(sl[1].text, iter_time+1);
            } else if (sl[0]=="#" && sl[1]=="include" && nw>2){
                if (iter_time+1>MAX_ITER_TIMES){ fprintf(stderr, "%s : %s[%d] : too many includes\n", software_name, filename, iline); success = false; }
                else success &= read_solvent(sl[2].text, iter_time+1);
            } else continue;
        } else if (sl[0].text[0] == ';' || (sl[0].text[0] == '/'&&sl[0].text[1] == '/')){
            continue;
        } else if (sl[0].text[0]=='['){
            if (sl[0]=="[solvent]" || (nw>2&&sl[0]=="["&&sl[1]=="solvent"&&sl[2]=="]")) on_compile = 1;
            else if (sl[0]=="[rismhi]" || (nw>2&&sl[0]=="["&&sl[1]=="rismhi"&&sl[2]=="]")) on_compile = 1;
            else if (sl[0]=="[rism3d]" || (nw>2&&sl[0]=="["&&sl[1]=="rism3d"&&sl[2]=="]")) on_compile = 1;
            else if (sl[0]=="[rismhi3d]" || (nw>2&&sl[0]=="["&&sl[1]=="rismhi3d"&&sl[2]=="]")) on_compile = 1;
            else if (sl[0]=="[hi]" || (nw>2&&sl[0]=="["&&sl[1]=="hi"&&sl[2]=="]")) on_compile = 1;
            else if (sl[0]=="[atom]" || (nw>2&&sl[0]=="["&&sl[1]=="atom"&&sl[2]=="]")) on_compile = 2;
            else { if (show_warning) fprintf(stderr, "%s : %s[%d] : warning : section %s%s%s ignored\n", software_name, filename, iline, nw>2?"[":"", nw>2?sl[1].text:sl[0].text, nw>2?"]":""); on_compile = 0; }
        } else {
            if (on_compile==1){         // [solvent] section
              if (nw>1){
                if (sl[0]=="gvv"){
                    strncpy(szfo_gvv, sl[1].text, sizeof(szfo_gvv));
                } else if (sl[0]=="zvv" || sl[0]=="zeta"){
                    strncpy(szfo_zvv, sl[1].text, sizeof(szfo_zvv));
                } else if (sl[0]=="dr"){
                    drrism = drhi = atof(sl[1].text); use_bins = false;
                } else if (sl[0]=="drrism"){
                    drrism = atof(sl[1].text); use_bins = false;
                } else if (sl[0]=="drhi"){
                    drhi = atof(sl[1].text);
                } else if (sl[0]=="rcoul"){
                    rcoul = atof(sl[1].text);
                } else if (sl[0]=="rvdw"){
                    rvdw = atof(sl[1].text);
                } else if (sl[0]=="ffoplsaa" || sl[0]=="ffopls" || sl[0]=="oplsaaff" || sl[0]=="oplsff"){
                    arith_sigma = false; ff_specifier = 2;
                } else if (sl[0]=="ffamber" || sl[0]=="amberff"){
                    arith_sigma = true; ff_specifier = 1;
                } else if (sl[0]=="ffgaff" || sl[0]=="gaff"){
                    arith_sigma = true; ff_specifier = 3;
                } else if (sl[0]=="ff" || sl[0]=="forcefield"){
                    if (sl[1]=="opls" || sl[1]=="oplsaa"){
                        arith_sigma = false; ff_specifier = 2;
                    } else if (sl[1]=="gaff"){
                        arith_sigma = true; ff_specifier = 3;
                    } else if (sl[1]=="amber"){
                        arith_sigma = true; ff_specifier = 1;
                    } else {
                        fprintf(stderr, "%s : %s[%d] : unkown forcefield type: %s\n", software_name, filename, iline, sl[1].text);
                        success = false;
                    }
                }
              }
            } else if (on_compile==2){  // [atom] section
                char * args[MAX_WORD]; for (int i=0; i<nw; i++) args[i] = sl[i].text; int idx = 0;
                success = analysis_solvent_atom(args, &idx, nw, filename, iline);
            }
        }
    }

    if (iter_time==0){ int __chdir = chdir(path_cwd); }

    fclose(file);
    return success;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool establish_atom_indexing(AtomIndex * ai, AtomSite * as, int nas, int * as_i_in_mol, PDBAtomSet * pdb_atom_set){
    PDBAtom * a = pdb_atom_set->atom; int count = pdb_atom_set->count;
    bool success = true;

    int i_in_mol = 0; int iaa_last = -1;
    for (int i=0; i<count; i++){
        if (a[i].iaa!=iaa_last){ i_in_mol = 0; iaa_last = a[i].iaa; } else i_in_mol++;
        StringNS::string mole = a[i].mole; StringNS::string nmele = a[i].nmele;
        for (int j=0; j<nas; j++) if (as_i_in_mol[j]==i_in_mol) if (mole==as[j].mole && nmele==as[j].nmele){
            ai[i].index = as[j].id; ai[i].grp = as[j].grp; a[i].iaa = as[j].iaa;
            ai[i].is_center = false;
        }
    }

    return success;
}
double minimal_rcutoff(Vector box){
    double dc = box.x;
    if (dc>box.y) dc = box.y;
    if (dc<box.z) dc = box.z;
    return dc / 2;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char *argv[]){
  // Initialization
    bool success = true; bool allow_xvv_calculation = true;
    pdb_atom_set.count = 0; pdb_atom_set.atom = nullptr; pdb_atom_set.box = Vector(0, 0, 0);
    char * cwd = getcwd(path_cwd, sizeof(path_cwd)); strncpy(path_top, path_cwd, sizeof(path_top));

  // Andlysis command line and analysis top
    if (success) success = analysis_args(argc, &argv[0], 1);
    if (argc<=1) return 0;
    if (success) success = read_solvent(szfp, 0);
    if (success) success = analysis_args(argc, &argv[0], 1);
    if (szfo_gvv[0] || szfo_zvv[0]) allow_xvv_calculation = true; else allow_xvv_calculation = false;
    if (success){
        rcutoff = rvdw<rcoul? rcoul : rvdw;
        if (use_bins){
            drrism = rcutoff / bins;
            drhi = rcutoff / bins;
        } else {
            bins_rism = drrism<=0? 1 : rcutoff/drrism;
            bins_hi = drhi<=0? 1 : rcutoff/drhi;
            bins = bins_rism>bins_hi? bins_rism : bins_hi;
        }
        if (bins<=2 || rcutoff<0 || drrism<0 || drhi<0 || drrism>rcutoff*0.5 || drhi>rcutoff*0.5){
            if (allow_xvv_calculation){
                fprintf(stderr, "%s : error : incorrect dr(%g)/rc(%g) settings\n", software_name, drrism, rcutoff);
                success = false;
            }
        }
      // reassign group id
        int i_new_grp[MAX_SOL]; int now_i_new_grp = 1;
        for (int i=0; i<nav; i++){
            int i_same_grp_found = -1; for (int j=0; j<i&&i_same_grp_found<0; j++) if (av[i].grp==av[j].grp) i_same_grp_found = j;
            i_new_grp[i] = i_same_grp_found>=0? i_new_grp[i_same_grp_found] : now_i_new_grp++;
        }
        for (int i=0; i<nav; i++) av[i].grp = i_new_grp[i];
      // count the number of molecule types and groups
        for (int i=0; i<nav; i++) if (av[i].iaa+1>=nmv) nmv = av[i].iaa+1;
        for (int i=0; i<nav; i++) if (av[i].grp>=ngv) ngv = av[i].grp;
        //printf("total moles: %d, total grps: %d\n", nmv, ngv);
    }
    if (success){
        if (!szfo[0]) strcpy(szfo, "stdout");
        for (int i=0; i<nav; i++){
            int i_in_mol = 0; for (int j=0; j<i; j++) if (StringNS::string(av[i].mole)==av[j].mole) i_in_mol ++;
            av_i_in_mol[i] = i_in_mol;
        }
        //for (int i=0; i<nav; i++) printf("av_i_in_mol[%d, %s.%s, %d] = %d\n", av[i].id, av[i].mole, av[i].name, av[i].iaa, av_i_in_mol[i]);
        if (!ff_specifier){
            fprintf(stderr, "%s : error : forcefield type (-ff) not specified\n", software_name);
            success = false;
        }
    }
  // check name conflict
    if (success){
        for (int i=0; i<nav; i++) for (int j=0; j<i; j++){
            if (av[i].grp!=av[j].grp && StringNS::string(av[i].name)==av[j].name){
                fprintf(stderr, "%s : error : %s.%s(%d) of group %d and %s.%s(%d) of group %d should have different names\n", software_name, av[j].mole, av[j].name, av[j].id, av[j].grp, av[i].mole, av[i].name, av[i].id, av[i].grp);
                success = false;
            }
        }
    }
    //printf("paths:\n    %s\n    %s\n", path_cwd, path_top); printf("-o: %s %s %s\n", szfo, szfo_gvv, szfo_zvv);
    //for (int i=0; i<nav; i++) printf("atom[%d] of %d = %s.%s : %s, in grp %d\n", i, av[i].iaa, av[i].mole, av[i].name, av[i].nmele, av[i].grp);
  // set priority
    if (nice_level>0){ if (setpriority(PRIO_PROCESS, getpid(), nice_level)!=0) fprintf(stderr, "%s : warning : fail to change nice level\n", software_name); }
  // read -s file : establish conformation
    if (success){
        if (str_ext_fs=="pdb") success = read_system_file(&pdb_atom_set, szfs, nullptr);
        else if (str_ext_fs=="gro") success = read_system_file(&pdb_atom_set, nullptr, szfs);
        else success = false;
        if (success){ int iaa_last = -1;
          // reassign the iaa: starts from 1
            iaa_last = -1; int iaa_now = 0;
            for (int i=0; i<pdb_atom_set.count; i++){
                if (pdb_atom_set.atom[i].iaa!=iaa_last){
                    nm ++; iaa_last = pdb_atom_set.atom[i].iaa; pdb_atom_set.atom[i].iaa = ++iaa_now;
                } else pdb_atom_set.atom[i].iaa = iaa_now;
            }
            //printf("molecules: %d\n", nm);
          // generate the atom indexing
            iaa_last = -1;
            ai = (AtomIndex*) malloc(sizeof(AtomIndex) * pdb_atom_set.count); memset(ai, 0, sizeof(AtomIndex) * pdb_atom_set.count);
            for (int i=0; i<pdb_atom_set.count; i++){
                ElementNS::Element * element = ElementNS::get_element_parameter(ElementNS::get_atom_element(pdb_atom_set.atom[i].name));
                strncpy(pdb_atom_set.atom[i].nmele, element?element->name:"?", sizeof(pdb_atom_set.atom[i].nmele));
            }
            //for (int i=0; i<pdb_atom_set.count; i++) printf("%d.%d: %s.%s : %s\n", pdb_atom_set.atom[i].iaa, pdb_atom_set.atom[i].index, pdb_atom_set.atom[i].mole, pdb_atom_set.atom[i].name, pdb_atom_set.atom[i].nmele);
          // generate the molecule indexing
            mi = (MoleculeIndex*) malloc(sizeof(MoleculeIndex) * nm); memset(mi, 0, sizeof(MoleculeIndex) * nm);
            iaa_last = -1;
            for (int i=0; i<pdb_atom_set.count; i++){
                if (pdb_atom_set.atom[i].iaa!=iaa_last){ iaa_last = pdb_atom_set.atom[i].iaa; mi[iaa_last-1].mole = pdb_atom_set.atom[i].mole; mi[iaa_last-1].iaa = iaa_last; mi[iaa_last-1].ibegin = i; }
            }
            for (int i=0; i<nm; i++) mi[i].iend = i<nm-1? mi[i+1].ibegin : pdb_atom_set.count;
            //for (int i=0; i<nm; i++) printf("mole[%d]: %s, %d ~ %d\n", mi[i].iaa, mi[i].mole, mi[i].ibegin, mi[i].iend);
          // generate the group infomation
            for (int ig=0; ig<ngv; ig++){
                gpi[ig].grp = ig+1;
                gpi[ig].ia = -1; for (int i=0; i<nav && gpi[ig].ia<0; i++) if (av[i].grp == ig+1) gpi[ig].ia = i;
                gpi[ig].multi = 0; for (int i=0; i<nav; i++) if (av[i].grp == ig+1) gpi[ig].multi ++;
                gpi[ig].nmole = 0;
            }
            //for (int i=0; i<ngv; i++) printf("grp[%d] : %s.%s, multi=%d\n", gpi[i].grp, av[gpi[i].ia].mole, av[gpi[i].ia].name, gpi[i].multi);
        }
    }
  // specify the atom index
    if (success){
        success = establish_atom_indexing(ai, av, nav, av_i_in_mol, &pdb_atom_set);
        //for (int i=0; i<pdb_atom_set.count; i++) printf("%d.%d: %s.%s : %s : %d %d %d\n", pdb_atom_set.atom[i].iaa, pdb_atom_set.atom[i].index, pdb_atom_set.atom[i].mole, pdb_atom_set.atom[i].name, pdb_atom_set.atom[i].nmele, ai[i].index, ai[i].grp, ai[i].iaa);
        int n_non_matches = 0; for (int i=0; i<pdb_atom_set.count; i++){
            if (ai[i].index<=0 || ai[i].grp<=0 || ai[i].iaa<0){ n_non_matches ++; if (n_non_matches>100) break;
                success = false; fprintf(stderr, "%s : error : undefined atom (%d) %s.%s\n", software_name, pdb_atom_set.atom[i].index, pdb_atom_set.atom[i].mole, pdb_atom_set.atom[i].name);
            }
        }
        if (n_non_matches>100){ success = false; fprintf(stderr, "%s : error : too many undefined atoms\n", software_name); }
      // count atom numbers of each molecule
        for (int i=0; i<nm; i++) nmole_mv[av[ai[mi[i].ibegin].index-1].iaa] ++;
        //for (int i=0; i<nmv; i++) printf("molecule %d has %d molecules\n", i, nmole_mv[i]);
        for (int ia=0; ia<nav; ia++) gpi[av[ia].grp-1].nmole = nmole_mv[av[ia].iaa];
        //for (int i=0; i<ngv; i++) printf("grp[%d] : %s.%s, multi=%d, nmole=%d\n", gpi[i].grp, av[gpi[i].ia].mole, av[gpi[i].ia].name, gpi[i].multi, gpi[i].nmole);
      // search for center of each molecule, by the coordinate that given in -s file
        int center_polling[MAX_SOL]; memset(center_polling, 0, sizeof(center_polling));
//for (int i=0; i<nm; i++) printf("mole[%d]: %s, %d ~ %d\n", mi[i].iaa, mi[i].mole, mi[i].ibegin, mi[i].iend);
        for (int im=0; im<nm; im++){
            double mass = 0; Vector v = Vector(0, 0, 0);
            for (int i=mi[im].ibegin; i<mi[im].iend; i++){ int iv = ai[i].index-1;
                mass += av[iv].mass; v += pdb_atom_set.atom[i].r * av[iv].mass;
            }
            if (mass>0) v /= mass; double r_min = -1; int i_min = -1;
            for (int i=mi[im].ibegin; i<mi[im].iend; i++){
                double r = pbc_v_sub(pdb_atom_set.atom[i].r, v, pdb_atom_set.box).mod();
                if (i_min<0 || r<r_min){ r_min = r; i_min = i; }
            }
            if (i_min>=0 && i_min<pdb_atom_set.count){
                int iv = ai[i_min].index-1;
                center_polling[iv] ++;
            }
        }
        //for (int i=0; i<nav; i++) printf("center_polling[%d of %d] = %d\n", i, av[i].iaa, center_polling[i]);
        for (int im=0; im<nmv; im++){
            int i_polling = -1;
            for (int iv=0; iv<nav; iv++) if (av[iv].iaa == im){
                if (i_polling<0 || center_polling[iv]>center_polling[i_polling]) i_polling = iv;
            }
            mi[im].icenter = i_polling;
            //printf("mole %d (%s) select center atom: %d (%s.%s): mi[%d].icenter = %d\n", im, mi[im].mole, mi[im].icenter , av[mi[im].icenter ].mole, av[mi[im].icenter].name, im, mi[im].icenter);
        }
        for (int ia=0; ia<pdb_atom_set.count; ia++){
            ai[ia].is_center = ai[ia].index-1==mi[av[ai[ia].index-1].iaa].icenter? true : false;
            //printf("%satom %d.%d (%s.%s), iaa_v=%d, mi[iaa_v]=%d\n%s", ai[ia].is_center?"\33[31m":"", av[ai[ia].index-1].iaa, av[ai[ia].index-1].id, av[ai[ia].index-1].mole, av[ai[ia].index-1].name, av[ai[ia].index-1].iaa, mi[av[ai[ia].index-1].iaa].icenter, "\33[0m");

        }
    }
  // prepare -f file
    if (success){
        if (str_ext_ff=="pdb") success = prepare_input(szff, nullptr, nullptr);
        else if (str_ext_ff=="gro") success = prepare_input(nullptr, szff, nullptr);
        else if (str_ext_ff=="xtc") success = prepare_input(nullptr, nullptr, szff);
        else success = false;
    }
  // prepare other things
    double *** mol_pair_dist = nullptr; double *** mol_pair_dist_count = nullptr;
    double *** gvv = nullptr;
    double *** uvv = nullptr; double *** uvv_count = nullptr;
    double *** zeta = nullptr;
    if (success){
        memset(dipole_mv, 0, sizeof(dipole_mv)); memset(n_dipole_mv, 0, sizeof(n_dipole_mv));
        mol_pair_dist = init_tensor3d<double>(nav, nav, 2); mol_pair_dist_count = init_tensor3d<double>(nav, nav, 1);
        if (allow_xvv_calculation){
            gvv = init_tensor3d<double>(ngv, ngv, bins_rism);
            uvv = init_tensor3d<double>(nmv, nmv, bins_hi); uvv_count = init_tensor3d<double>(nmv, nmv, bins_hi);
            zeta = init_tensor3d<double>(nmv, nmv, bins_hi);
        }
        if (!mol_pair_dist || (allow_xvv_calculation && (!gvv || !zeta || !uvv || !uvv_count))){ fprintf(stderr, "%s : malloc failure\n", software_name); success = false; }
        else {
            clear_tensor3d<double>(mol_pair_dist, nav, nav, 2);
            clear_tensor3d<double>(mol_pair_dist_count, nav, nav, 2);
            if (gvv) clear_tensor3d<double>(gvv, ngv, ngv, bins_rism);
            if (uvv) clear_tensor3d<double>(uvv, nmv, nmv, bins_hi);
            if (uvv_count) clear_tensor3d<double>(uvv_count, nmv, nmv, bins_hi);
            if (zeta) clear_tensor3d<double>(zeta, nmv, nmv, bins_hi);
        }
    }

  // read the trajectory
    int iframe = 0; double nframe = 0; Vector minimal_box = Vector(-1, -1, -1);
    if (success){
        TPAppendix tpa; tpa.time = 0; double time_last = 0; double last_runout_time = get_current_time_double(); int report_count = 0;
        while (true){
            int read_frame_ret = read_frame(&pdb_atom_set, &iframe, &tpa);
            if (read_frame_ret==0){
                break;
            } else if (read_frame_ret==-1){
                fprintf(stderr, "%s : error : cannot read trajectory %s\n", software_name, szff);
                break;
            } else if (read_frame_ret==-2){
                fprintf(stderr, "%s : error : incorrect number of atoms: %s\n", software_name, szff);
                break;
            } else if (read_frame_ret<0){
                fprintf(stderr, "%s : error : fail for %s\n", software_name, szff);
                break;
            }
            bool skip_frame = false; if (time_begin>0 && tpa.time<time_begin) skip_frame = true; if (time_end>0 && tpa.time>time_end) break;
            if (tpa.time-time_last < time_step) continue;
            double runout_time = get_current_time_double(); if (runout_time-last_runout_time>0.1 || !skip_frame){ last_runout_time = runout_time; report_count++;
                if (skip_frame) fprintf(stderr, "%s : %12s : %11.1f : skipped     %s", software_name, szff, tpa.time, show_debug?"\n":"\r");
                else fprintf(stderr, "%s : %12s : %11.1f %8.3f,%8.3f,%8.3f     %s", software_name, szff, tpa.time, pdb_atom_set.box.x, pdb_atom_set.box.y, pdb_atom_set.box.z, show_debug?"\n":"\r");
            }
            if (skip_frame) continue; nframe ++; time_last = tpa.time;
          // ==========================================
            PDBAtom * a = pdb_atom_set.atom;
            if (minimal_box.x<0 || minimal_box.x<pdb_atom_set.box.x) minimal_box.x = pdb_atom_set.box.x;
            if (minimal_box.y<0 || minimal_box.y<pdb_atom_set.box.y) minimal_box.y = pdb_atom_set.box.y;
            if (minimal_box.z<0 || minimal_box.z<pdb_atom_set.box.z) minimal_box.z = pdb_atom_set.box.z;
           // PBC
            for (int i=0; i<pdb_atom_set.count; i++) pdb_atom_set.atom[i].r = pbc_rm_v(pdb_atom_set.atom[i].r, pdb_atom_set.box);
           // density
            inverse_V += 1 / (pdb_atom_set.box.x * pdb_atom_set.box.y * pdb_atom_set.box.z);
            for (int im=0; im<nm; im++){
                Vector dp = Vector(0,0,0); for (int i=mi[im].ibegin; i<mi[im].iend; i++){
                    dp = dp + pbc_v_sub_partial(a[i].r, a[mi[im].ibegin].r, pdb_atom_set.box) * av[ai[i].index-1].charge;
                }
                if (mi[im].iend > mi[im].ibegin){
                    dipole_mv[av[ai[mi[im].ibegin].index-1].iaa] += dp.mod(); n_dipole_mv[av[ai[mi[im].ibegin].index-1].iaa] += 1.0;
                }
            }
           // bond
            for (int im=0; im<nm; im++){
                for (int ia=mi[im].ibegin; ia<mi[im].iend; ia++) for (int ja=mi[im].ibegin; ja<mi[im].iend; ja++) if (ia!=ja){
                    double r = pbc_v_sub_partial(a[ia].r, a[ja].r, pdb_atom_set.box).mod();
                    int iv = ai[ia].index-1;
                    int jv = ai[ja].index-1;
                    mol_pair_dist[iv][jv][0] += r; mol_pair_dist[jv][iv][0] += r;
                    mol_pair_dist[iv][jv][1] += r*r; mol_pair_dist[jv][iv][1] += r*r;
                    mol_pair_dist_count[iv][jv][0] += 1; mol_pair_dist_count[jv][iv][0] += 1;
//if (show_debug) printf("  bond %2d %2d : %12g (%8.3f %8.3f %8.3f) - (%8.3f %8.3f %8.3f)\n", iv, jv, r, a[ia].r.x, a[ia].r.y, a[ia].r.z, a[ja].r.x, a[ja].r.y, a[ja].r.z);
                }
            }
           // gvv and zeta
            if (allow_xvv_calculation && (szfo_gvv[0] || szfo_zvv[0])) for (int im=0; im<nm; im++) for (int jm=im+1; jm<nm; jm++) if (im!=jm){
                Vector ri, rj; ri = rj = Vector(0, 0, 0); double Evdw = 0; double Ecoul = 0;
                for (int ia=mi[im].ibegin; ia<mi[im].iend; ia++) for (int ja=mi[jm].ibegin; ja<mi[jm].iend; ja++){
                    //double r = pbc_rm_v(a[ia].r - a[ja].r, pdb_atom_set.box).mod();
                    double r = pbc_v_sub_partial(a[ia].r, a[ja].r, pdb_atom_set.box).mod();
                  // gvv
                    if (szfo_gvv[0]){
                        int ir = (int) floor(r / drrism);
                        if (ir>=0 && ir<bins_rism){
                            double inc = 1.0 / (4*PI*(r*r*drrism));
                            //double ra = ir*drrism; double rb = (ir+1)*drrism;
                            //double inc = 1.0 / (4*PI*(rb*rb*rb - ra*ra*ra)/3);
                            gvv[ai[ia].grp-1][ai[ja].grp-1][ir] += inc;
                            gvv[ai[ja].grp-1][ai[ia].grp-1][ir] += inc;
                        }
                    }
                  // zeta
                    if (szfo_zvv[0]){
                        if (ai[ia].is_center && ai[ja].is_center){ ri = a[ia].r; rj = a[ja].r; }
                        int iv = ai[ia].index-1;
                        int jv = ai[ja].index-1;
                        double sigma = arith_sigma? (av[iv].sigma + av[jv].sigma)/2 : (av[iv].sqrt_sigma * av[jv].sqrt_sigma);
                        double epsilon = 4 * av[iv].sqrt_epsilon * av[jv].sqrt_epsilon;
                        double s = sigma / r; double s2 = s*s; double s6 = s2*s2*s2; double s12 = s6*s6;
                        Evdw += epsilon * (s12 - s6);
                        Ecoul += COULCOOEF * av[iv].charge * av[jv].charge / r;
//if (epsilon*(s12-s6)>10) printf("%s (%d) vs %s (%d) : sigma %f %f -> %f, epsilon %f %f -> %f, r=%f\n", av[iv].name, pdb_atom_set.atom[ia].index, av[jv].name, pdb_atom_set.atom[ja].index, av[iv].sigma, av[jv].sigma, sigma, av[iv].epsilon, av[jv].epsilon, epsilon, r);
                    }
                }
              // zeta
                if (szfo_zvv[0]){
                    double r = pbc_v_sub_partial(ri, rj, pdb_atom_set.box).mod();
                    int ir = (int) floor(r / drhi);
                    int imv = av[ai[mi[im].ibegin].index-1].iaa;
                    int jmv = av[ai[mi[jm].ibegin].index-1].iaa;
                    if (ir>=0 && ir<bins_hi){
                        uvv[imv][jmv][ir] += (Ecoul + Evdw);
                        uvv_count[imv][jmv][ir] += 1;
                        uvv_count[jmv][imv][ir] += 1;
                    }
                }
            }
        }
        if (report_count>0) fprintf(stderr, "\n");
    }

  // post handling
    if (success && nframe>0){
        inverse_V /= nframe;
        for (int i=0; i<nmv; i++) dipole_mv[i] /= n_dipole_mv[i];
        if (allow_xvv_calculation && szfo_gvv[0]){
            for (int ir=0; ir<bins_rism; ir++) for (int iv=0; iv<ngv; iv++) for (int jv=iv; jv<ngv; jv++) gvv[iv][jv][ir] /= nframe * gpi[iv].multi * gpi[jv].multi * gpi[iv].nmole * gpi[jv].nmole * inverse_V;
        }
        if (allow_xvv_calculation && szfo_zvv[0]){
            for (int ir=0; ir<bins_hi; ir++){
                for (int iv=0; iv<nmv; iv++) for (int jv=iv; jv<nmv; jv++){
                    if (uvv_count[iv][jv][ir]<=0) uvv[iv][jv][ir] = 0; else uvv[iv][jv][ir] /= uvv_count[iv][jv][ir];
                    double ra = ir*drhi; double rb = (ir+1)*drhi;
                    uvv_count[iv][jv][ir] /= nframe * nmole_mv[iv] * nmole_mv[jv] * inverse_V * (4*PI*(rb*rb*rb - ra*ra*ra)/3);
    //printf(" %11f %11f (%11f)", uvv[iv][jv][ir], uvv_count[iv][jv][ir], uvv[iv][jv][ir] * uvv_count[iv][jv][ir]);
                }
    //printf("\n");
            }
        }
        if (allow_xvv_calculation && szfo_zvv[0]){
            double rcutoff_max = minimal_rcutoff(minimal_box);
            int bins_hi_max = rcutoff_max / drhi; if (bins_hi_max > bins_hi) bins_hi_max = bins_hi;
            for (int iv=0; iv<nmv; iv++) for (int jv=iv; jv<nmv; jv++){
                double Integral = 0; for (int ir=0; ir<bins_hi_max; ir++){
                    Integral += (ir==0? (uvv[iv][jv][ir+1]-uvv[iv][jv][ir]) : ir==bins_hi_max-1? (uvv[iv][jv][ir] - uvv[iv][jv][ir-1]) : (uvv[iv][jv][ir+1]- uvv[iv][jv][ir-1])/2) * uvv_count[iv][jv][ir];
                    zeta[iv][jv][ir] = Integral;
                }
                Integral = zeta[iv][jv][bins_hi_max-1];
                for (int ir=0; ir<bins_hi_max; ir++) zeta[iv][jv][ir] -= Integral;
                for (int ir=bins_hi_max; ir<bins_hi; ir++) zeta[iv][jv][ir] = 0;
            }
        }
      // rename atoms and molecules
        if (force_rename_atom_names){
            for (int iv=0; iv<nav; iv++) snprintf(av[iv].name, sizeof(av[iv].name), "%s%d%s", av[iv].nmele, av[iv].id, av[iv].charge>0.95?"+":av[iv].charge<-0.95?"-":"");
        }
        if (n_nss>0){
            for (int iv=0; iv<nav; iv++){ int ifound = -1; StringNS::string key = av[iv].mole;
                for (int i=0; i<n_nss; i++) if (key == nss[i].name[0]) strncpy(av[iv].mole, nss[i].name[1], sizeof(av[iv].mole));
            }
        }
    }

  // report : -o
    if (success){
      // prepare file and path
        FILE * flog = nullptr;
        if (StringNS::string(szfo)=="con"||StringNS::string(szfo)=="screen"||StringNS::string(szfo)=="stdout"){ flog = stdout;
        } else if (StringNS::string(szfo)=="stderr"){ flog = stderr;
        } else { if (force_override) flog = fopen(szfo, "w"); else { flog = fopen(szfo, "r"); if (flog){ fclose(flog); flog = stdout; fprintf(stderr, "%s : warning : -o %s already exist, output to screen instead\n", software_name, szfo); } else flog = fopen(szfo, "w"); }
        }
        if (!flog){ flog = stdout; fprintf(stderr, "%s : warning : cannot write -o %s, output to screen instead\n", software_name, szfo); }
        bool istty = flog? isatty(fileno(flog)) : true;
      // output command line
        fprintf(flog, "# generated by %s %s\n# ", software_name, software_version); for (int i=0; i<argc; i++) fprintf(flog, " %s", argv[i]); fprintf(flog, "\n");
      // report [solvent]
        fprintf(flog, "[solvent]\n"); // fprintf(flog, "%s[solvent]\n%s", istty?"\33[31m":"", istty?"\33[30m":"");
        if (ff_specifier==1) fprintf(flog, " ff           amber\n");
        if (ff_specifier==2) fprintf(flog, " ff           oplsaa\n");
        if (ff_specifier==3) fprintf(flog, " ff           gaff\n");
        fprintf(flog, " rvdw         %g\n", rvdw);
        fprintf(flog, " rcoul        %g\n", rcoul);
        if (allow_xvv_calculation){
            fprintf(flog, " drrism       %g\n", drrism);
            if (StringNS::string(szfo_gvv)!="con"&&StringNS::string(szfo_gvv)!="screen"&&StringNS::string(szfo_gvv)!="stdout"&&StringNS::string(szfo_gvv)!="stderr") fprintf(flog, " gvv          %s\n", szfo_gvv); else fprintf(flog, " gvv          \n");
            fprintf(flog, " drhi         %g\n", drhi);
            if (StringNS::string(szfo_zvv)!="con"&&StringNS::string(szfo_zvv)!="screen"&&StringNS::string(szfo_zvv)!="stdout"&&StringNS::string(szfo_zvv)!="stderr") fprintf(flog, " zvv          %s\n", szfo_zvv); else fprintf(flog, " zvv          \n");
        }
        fprintf(flog, " density      "); for (int i=0; i<nmv; i++){ fprintf(flog, "%g ", nmole_mv[i] * inverse_V); }; fprintf(flog, "\n");
        fprintf(flog, " bulk-density "); for (int i=0; i<nmv; i++){ fprintf(flog, "%g ", nmole_mv[i] * inverse_V); }; fprintf(flog, "\n");
        fprintf(flog, " dipole       "); for (int i=0; i<nmv; i++) fprintf(flog, "%g ", dipole_mv[i]); fprintf(flog, "\n");
        double total_dp = 0; for (int i=0; i<nmv; i++) total_dp += nmole_mv[i] * dipole_mv[i];
        fprintf(flog, " dielect      "); for (int i=0; i<nmv; i++) fprintf(flog, "%g ", total_dp * inverse_V * 42.95537+1); fprintf(flog, "\n");
        //fprintf(flog, "\n");
      // report [atom]
        fprintf(flog, "[atom]\n"); // fprintf(flog, "%s[atom]\n%s", istty?"\33[31m":"", istty?"\33[30m":"");
        fprintf(flog, "%s#%6s %9s %3s %3s %11s %11s %11s\n%s", istty?"\33[37m":"", "name", "mole", "idx", "grp", "charge", "sigma", "epsilon", istty?"\33[0m":"");
        for (int i=0; i<nav; i++) fprintf(flog, " %6s %9s %3d %3d %11f %11f %11f\n", av[i].name, av[i].mole, av[i].id, av[i].grp, av[i].charge, av[i].sigma, av[i].epsilon);
        //fprintf(flog, "\n");
      // report [bond]
        fprintf(flog, "[bond]\n"); // fprintf(flog, "%s[bond]\n%s", istty?"\33[31m":"", istty?"\33[30m":"");
        fprintf(flog, "%s#%11s %11s %11s %11s\n%s", istty?"\33[37m":"", "atom1", "atom2", "bond(nm)", "stdev(nm)", istty?"\33[0m":"");
        for (int i=0; i<nav; i++) for (int j=i+1; j<nav; j++) if (av[i].iaa == av[j].iaa){
            mol_pair_dist[i][j][0] /= mol_pair_dist_count[i][j][0]; mol_pair_dist[i][j][1] /= mol_pair_dist_count[i][j][0];
            mol_pair_dist[i][j][1] = sqrt(fabs(mol_pair_dist[i][j][1] - mol_pair_dist[i][j][0]*mol_pair_dist[i][j][0]));
            char compound_names[2][72];
            snprintf(compound_names[0], sizeof(compound_names[0]), "%s.%s", av[i].mole, av[i].name);
            snprintf(compound_names[1], sizeof(compound_names[1]), "%s.%s", av[j].mole, av[j].name);
            fprintf(flog, " %11s %11s %11f %11f\n", compound_names[0], compound_names[1], mol_pair_dist[i][j][0], mol_pair_dist[i][j][1]);
        }
        //fprintf(flog, "\n");
      // finished
        if (flog && flog!=stdout && flog!=stderr) fclose(flog);
        fprintf(flog, "\n");
    }

  // report : gvv
    if (success && szfo_gvv[0]){
      // prepare file and path
        FILE * flog = nullptr;
        if (StringNS::string(szfo_gvv)=="con"||StringNS::string(szfo_gvv)=="screen"||StringNS::string(szfo_gvv)=="stdout"){ flog = stdout;
        } else if (StringNS::string(szfo_gvv)=="stderr"){ flog = stderr;
        } else { if (force_override) flog = fopen(szfo_gvv, "w"); else { flog = fopen(szfo_gvv, "r"); if (flog){ fclose(flog); flog = stdout; fprintf(stderr, "%s : warning : -o %s already exist, output to screen instead\n", software_name, szfo_gvv); } else flog = fopen(szfo_gvv, "w"); }
        }
        if (!flog){ flog = stdout; fprintf(stderr, "%s : warning : cannot write gvv file %s, output to screen instead\n", software_name, szfo_gvv); }
        bool istty = flog? isatty(fileno(flog)) : true;
      // output gvv
        if (flog==stdout || flog==stderr){
            if (StringNS::string(szfo_gvv)=="stdout" || StringNS::string(szfo_gvv)=="stderr" || StringNS::string(szfo_gvv)=="screen" || StringNS::string(szfo_gvv)=="con") fprintf(flog, "%sFile -gvv %s%s.%gnm.gvv :\n%s", istty?"\33[31m":"", av[0].mole, nmv>1?"-etc":"", drrism, istty?"\33[0m":"");
            else fprintf(flog, "%sFile -gvv %s :\n%s", istty?"\33[31m":"", szfo_gvv, istty?"\33[0m":"");
        }
      // output gvv
        double rcutoff_max = minimal_rcutoff(minimal_box);
        for (int ir=0; ir<bins_rism; ir++){
            if (ir*drrism >= rcutoff_max) break;
            for (int iv=0; iv<ngv; iv++) for (int jv=iv; jv<ngv; jv++){
                fprintf(flog, " %11f", gvv[iv][jv][ir]);
            }
            fprintf(flog, "\n");
        }
      // finished
        if (flog && flog!=stdout && flog!=stderr) fclose(flog);
    }

  // report : zvv
    if (success && szfo_zvv[0]){
      // prepare file and path
        FILE * flog = nullptr;
        if (StringNS::string(szfo_zvv)=="con"||StringNS::string(szfo_zvv)=="screen"||StringNS::string(szfo_zvv)=="stdout"){ flog = stdout;
        } else if (StringNS::string(szfo_zvv)=="stderr"){ flog = stderr;
        } else { if (force_override) flog = fopen(szfo_zvv, "w"); else { flog = fopen(szfo_zvv, "r"); if (flog){ fclose(flog); flog = stdout; fprintf(stderr, "%s : warning : -o %s already exist, output to screen instead\n", software_name, szfo_zvv); } else flog = fopen(szfo_zvv, "w"); }
        }
        if (!flog){ flog = stdout; fprintf(stderr, "%s : warning : cannot write gvv file %s, output to screen instead\n", software_name, szfo_zvv); }
        bool istty = flog? isatty(fileno(flog)) : true;
      // output gvv
        if (flog==stdout || flog==stderr){
            if (StringNS::string(szfo_zvv)=="stdout" || StringNS::string(szfo_zvv)=="stderr" || StringNS::string(szfo_zvv)=="screen" || StringNS::string(szfo_zvv)=="con") fprintf(flog, "%sFile -zvv %s%s.%gnm.zeta :\n%s", istty?"\33[31m":"", av[0].mole, nmv>1?"-etc":"", drrism, istty?"\33[0m":"");
            else fprintf(flog, "%sFile -zvv %s :\n%s", istty?"\33[31m":"", szfo_zvv, istty?"\33[0m":"");
        }
      // output zvv
        double rcutoff_max = minimal_rcutoff(minimal_box);
        for (int ir=0; ir<bins_rism; ir++){
            if (ir*drrism >= rcutoff_max) break;
            for (int iv=0; iv<nmv; iv++) for (int jv=iv; jv<nmv; jv++){
                fprintf(flog, " %11f", zeta[iv][jv][ir]);
            }
            fprintf(flog, "\n");
        }
      // finished
        if (flog && flog!=stdout && flog!=stderr) fclose(flog);
    }



    //if (ai) free(ai); if (mol_pair_dist) free(mol_pair_dist); if (mol_pair_dist_count) free(mol_pair_dist_count);
    //if (gvv) free(gvv); if (uvv) free (uvv); if (uvv_count) free (uvv_count); if (guvv) free(guvv); if (zeta) free(zeta);
    return 0;
}
