//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//---------------   System Param: Atom List (Solute and Solvent)   ----------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class AtomSite { public:
  // fundamental parameters
    char name[MAX_NAME]; char mole[MAX_NAME]; const char * nele;
    double mass; double charge; double charge_esp;  // charge for ES field, charge_esp for reporting ES energy
    double sigma, sqrt_sigma; double epsilon, sqrt_epsilon; int iaa;
    int id, grp, multi; bool is_key;
    unsigned int reserved;
  // advanced parameters
    double reverse_rism_factor;
    double ren_charge, ren_bond, ren_dielect;
  // initialization
    void init(int _id, int _grp, char * _mole, int _iaa, char * _name, double _mass, double _charge, double _sigma, double _epsilon){
        memset(name, 0, sizeof(name)); int len = strlen(_name); if (len>MAX_NAME-1) len = MAX_NAME-1; memcpy(name, _name, len);
        memset(mole, 0, sizeof(mole)); len = strlen(_mole); if (len>MAX_NAME-1) len = MAX_NAME-1; memcpy(mole, _mole, len);
        nele = ElementNS::get_atom_element(name);
        iaa = _iaa; mass = _mass; charge = _charge; charge_esp = _charge;
        sigma = _sigma; sqrt_sigma = sqrt(fabs(sigma));
        epsilon = _epsilon<0? -1e-15 : _epsilon; sqrt_epsilon = (_epsilon<0? -1 : 1) * sqrt(fabs(epsilon));
        id = _id; grp = _grp; multi = 1; is_key = false; reserved = 0;
        reverse_rism_factor = -1;
        ren_charge = 0; ren_bond = 1; ren_dielect = 1;
        //printf("atom site %s.%s: %12f %12f %12f %12f\n", mole, name, charge, sigma, epsilon, sqrt_epsilon);
    }
};
#define MAX_BONDS_PER_ATOM 6
class SoluteAtomSite { public:
    int index; char name[8]; int iaa; char mole[8];
    double mass, charge; double sigma, sqrt_sigma, epsilon, sqrt_epsilon;
    int nbond; int ibond[MAX_BONDS_PER_ATOM];
    int i_sigma_list; int n_sigma_list; int i_epsilon_list; int n_epsilon_list;
    unsigned int reserved;
    void init(int _index, const char * _name, int _iaa, const char * _mole, double _mass, double _charge, double _sigma, double _epsilon){
        index = _index;   memset(name, 0, sizeof(name)); strncpy(name, _name, sizeof(name)-1);
        iaa = _iaa; memset(mole, 0, sizeof(mole)); strncpy(mole, _mole, sizeof(mole)-1);
        mass = _mass; charge = _charge;
        sigma = _sigma; sqrt_sigma = sqrt(fabs(sigma));
        epsilon = _epsilon<0? -1e-15 : _epsilon; sqrt_epsilon = (_epsilon<0? -1 : 1) * sqrt(fabs(epsilon));
        i_sigma_list = -1; n_sigma_list = 0; i_epsilon_list = -1; n_epsilon_list = 0;
        reserved = 0;
    }
};
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class RDFGroup { public:
  // index: >0: index; ==0: mathces anything; <0: ignore.
  // mol and atom name: can be string; or leave blank ("") or use "*" to match anything
    int is, iv, grp; // index of solute and solvent. Begin with 1.
    char ms[MAX_NAME], mv[MAX_NAME], as[MAX_NAME], av[MAX_NAME]; // mol and atom name, ignored if
    void init(int _is, int _iv, StringNS::string _ms, StringNS::string _as, StringNS::string _mv, StringNS::string _av, int _grp=-1){
        memset(this, 0, sizeof(RDFGroup));
        is = _is; iv = _iv; grp = _grp;
        memcpy(ms, _ms.text, _ms.length>(MAX_NAME-1)?(MAX_NAME-1):_ms.length);
        memcpy(as, _as.text, _as.length>(MAX_NAME-1)?(MAX_NAME-1):_as.length);
        memcpy(mv, _mv.text, _mv.length>(MAX_NAME-1)?(MAX_NAME-1):_mv.length);
        memcpy(av, _av.text, _av.length>(MAX_NAME-1)?(MAX_NAME-1):_av.length);
    }
};
class RDF_data { public:
    int is, iv; // index of solute and solvent, must be positive
    double * g, * n;
};
class RDF_datas { public:
    int n_rdf_datas[2]; RDF_data * rdf_datas[2];
};
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-----------------------------   HI Equation Solver   ----------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class HIEquationSolver {
  private:
    bool trim;
    int count;
    double * xa, * ya;
  public:
    double lse_a, lse_b;
    #ifdef _EXPERIMENTAL_
        IET_Param_Exp ex;
    #endif
  public:
    void set_trim(){ trim = true; }
    void unset_trim(){ trim = false; }
    double f(double x){
      #ifdef _EXPERIMENTAL_
        return ex.HIEquationSolver_f_exp(x, lse_a, lse_b);
      #else
        double pre = lse_b * lse_a * exp((1-1/x)/lse_a); return ln(x+1e-15) + pre;
      #endif
    }
    void dispose(){
        if (count>0){ free(xa); free(ya); }
        count = 0; xa = ya = nullptr;
    }
    void prepare(double x_inf, double x_sup){
        for (int i=0; i<=count; i++){
            xa[i] = x_inf + (x_sup - x_inf) * i / (count);
            ya[i] = this->f(xa[i]);
//printf("%12g %12g\n", xa[i], ya[i]);
        }
    }
    void init(int N = 10000){
        lse_a = 0.3; lse_b = 51.8;
        count = N; trim = false;
        xa = (double*) memalloc(sizeof(double) * (count+1) * 2);
        ya = &xa[count+1];
        #ifdef _EXPERIMENTAL_
            ex.init();
        #endif
    }
    #ifdef _EXPERIMENTAL_
        void set_param(double _lse_a, double _lse_b, IET_Param_Exp * _ex, double _inf=0, double _sup=2, bool _set_trim=true){
            lse_a = _lse_a; lse_b = _lse_b; memcpy(&ex, _ex, sizeof(ex));
            trim = _set_trim; prepare(_inf, _sup);
        }
    #else
        void set_param(double _lse_a, double _lse_b, double _inf=0, double _sup=2, bool _set_trim=true){
            lse_a = _lse_a; lse_b = _lse_b;
            trim = _set_trim; prepare(_inf, _sup);
        }
    #endif
    double getx_ds(double y, double xinf, double xsup, double err){
        double inf = xinf; double sup = xsup;
        while (sup-inf>err){
            double xt = (inf + sup) / 2;
            double yt = f(xt);
            if (yt >= y) sup = xt;
            if (yt <= y) inf = xt;
        }
        return (inf + sup) / 2;
    }
    double getx_ds(double y, double err=1e-5){
        double inf = 0; double sup = 0;
        while (f(inf) > y) inf -= 1;
        while (f(sup) < y) sup += 1;
        return getx_ds(y, inf, sup, err);
    }
    double getx_bs(double y, double err=1e-5, bool further_ds=false){
        if (y>ya[count] || y<ya[0]){
            if (!trim) return getx_ds(y, err);
            else return y>ya[count]? xa[count] : xa[0];
        } else {
            int i = 0; int j = count;
            while (j-i>1){
                int k = (i+j)/2;
                if (ya[k]>=y) j = k;
                if (ya[k]<=y) i = k;
            }
            if (j==i) return xa[i];
            if (err > xa[j] - xa[i]){
                return xa[i];
            } else {
                if (further_ds){
                    return getx_ds(y, xa[i], xa[j], err);
                } else {
                    if (ya[i]==ya[j]) return (xa[i] + xa[j])/2;
                    return (xa[i]*(ya[j]-y) + xa[j]*(y-ya[i])) / (ya[j] - ya[i]); //
                }
            }
        }
    }
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-----------------------   System Param: Atom Name List   ------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class SYSITEM_AtomNameList {
  public:
    char name[MAX_NAME]; char mole[MAX_NAME]; int index; int grp; int iaa;
};
class SYSITEM_PairMapping { // mapping of gvv
  public:
    int grpi; int grpj; int col;
};
class SYSITEM_BondList {
  public:
    int grpi; int grpj; double bond; double bond_stdev;
    double weight; // weight must be (0,1], and all invalid numbers are recognized as 100%
};
class SYSITEM_ZetaList {
  public:
      //int iaai, iaaj; double zeta0, rc_zeta0;
      int iaai, iaaj; double deltaG_in_Kelvin, rc_zeta0;
};
//-------------------------------------
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, StringNS::string mole, StringNS::string atom, int default_ret){
    for (int i=begin; i<nal; i++){
        bool mole_match = mole=="*" || mole==al[i].mole;
        bool atom_match = atom=="*" || atom==al[i].name;
        if (mole_match && atom_match) return return_col==1? al[i].index : return_col==2? al[i].grp : return_col==3? al[i].iaa : i;
    }
    return default_ret;
}
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, const char * sz_mole, const char * sz_atom, int default_ret){
    StringNS::string mole = sz_mole? sz_mole : "*";
    StringNS::string atom = sz_atom? sz_atom : "*";
    return search_atom_list(return_col, al, begin, nal, mole, atom, default_ret);
}
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    char mole_name[MAX_NAME]; char atom_name[MAX_NAME]; memset(mole_name, 0, sizeof(mole_name)); memset(atom_name, 0, sizeof(atom_name));
    for (int i=0; compond[i]; i++){
        if (compond[i]=='.'){
            memcpy(mole_name, compond, i>(MAX_NAME-1)?(MAX_NAME-1):i); if (mole_name[0]==0) strcpy(mole_name, "*");
            const char * compond_next = &compond[i+1]; int compond_next_len = strlen(compond_next);
            memcpy(atom_name, compond_next, compond_next_len>(MAX_NAME-1)?(MAX_NAME-1):compond_next_len); if (atom_name[0]==0) strcpy(atom_name, "*");
            break;
        } else if (compond[i]==':' && compond[i+1]==':'){
            memcpy(mole_name, compond, i>(MAX_NAME-1)?(MAX_NAME-1):i); if (mole_name[0]==0) strcpy(mole_name, "*");
            const char * compond_next = &compond[i+2]; int compond_next_len = strlen(compond_next);
            memcpy(atom_name, compond_next, compond_next_len>(MAX_NAME-1)?(MAX_NAME-1):compond_next_len); if (atom_name[0]==0) strcpy(atom_name, "*");
            break;
        }
    }
    int compond_len = strlen(compond);
    if (!atom_name[0]){ strcpy(mole_name, "*"); memcpy(atom_name, compond, compond_len>(MAX_NAME-1)?(MAX_NAME-1):compond_len); }
    if (atom_name[0]=='*') return default_ret;
    return search_atom_list(return_col, al, begin, nal, mole_name, atom_name, default_ret);
}

/*int search_atom_list_2(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, StringNS::string compound, int default_ret){
    StringNS::string search_strings[2] = { "", "" };
    char sep_sv[4] = { '.', ':', '.', ':' };
    int nw = analysis_general_line(sep_sv, compound, search_strings, 2, false, true);
    if (nw==1){
        int imole = search_atom_list(return_col, al, begin, nal, search_strings[0], "*", -1);
        int iatom = search_atom_list(return_col, al, begin, nal, "*", search_strings[0], -1);
        if (iatom>=0) return iatom;
        else if (imole>=0) return imole;
        else return default_ret;
    } else if (nw==2){
        return search_atom_list(return_col, al, begin, nal, search_strings[0], search_strings[1], default_ret);
    }
    return default_ret;
}*/


int search_atom_list_index(SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    return search_atom_list(1, al, begin, nal, compond, default_ret);
}
int search_atom_list_grp(SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    return search_atom_list(2, al, begin, nal, compond, default_ret);
}
int search_mole_list(SYSITEM_AtomNameList * al, int begin, int nal, const char * sz_mole, int default_ret){
    if (StringNS::is_string_number(sz_mole)){
        int iaa = atoi(sz_mole);
        return iaa<=0 ? default_ret : iaa-1;
    } else {
        return search_atom_list(3, al, begin, nal, sz_mole, "*", default_ret);
    }
}
