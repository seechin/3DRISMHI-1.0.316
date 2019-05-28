
int analysis_line_words(StringNS::string line, StringNS::string * wl, int nwl_max){
    int nwl = 0; int idx = 0;
    while (idx < line.length && nwl < nwl_max){
        wl[nwl] = StringNS::seek_first_word(line, &idx);
        if (wl[nwl].length <= 0) break; else nwl ++;
    }
    return nwl;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Atom in pdb's format
class PDBAtom {
  public:
    int     index;
    char    name[6];
    char    mole[6];
    char    ichain[2];
    int     iaa;
    double  occupacy, Bfactor;
    char    nmele[4];
    Vector  r;
  public:
    PDBAtom(){ init(); }
    void init(){
        index = 0; iaa = 0; occupacy = Bfactor = 0;
        name[0] = mole[0] = ichain[0] = nmele[0] = 0;
        r = Vector(0,0,0);
    }
    void read_pdb_line(StringNS::string sline){
        index = string_to_int(seek_first_word(sline.Substring(6, 5)));
        StringNS::string t;
        memset(name, 0, sizeof(name));
          t = seek_first_word(sline.Substring(12, 4));
          if (t.length>0) memcpy(name, t.text, t.length);
        memset(mole, 0, sizeof(mole));
          t = seek_first_word(sline.Substring(17, 4));
          if (t.length>0) memcpy(mole, t.text, t.length);
        memset(ichain, 0, sizeof(ichain));
          t = seek_first_word(sline.Substring(21, 1));
          if (t.length>0) memcpy(ichain, t.text, t.length);
        iaa      = string_to_int(seek_first_word(sline.Substring(22, 4)));
        r.x      = string_to_double(seek_first_word(sline.Substring(30, 8))) / 10;
        r.y      = string_to_double(seek_first_word(sline.Substring(38, 8))) / 10;
        r.z      = string_to_double(seek_first_word(sline.Substring(46, 8))) / 10;
        occupacy = string_to_double(seek_first_word(sline.Substring(54, 6)));
        Bfactor  = string_to_double(seek_first_word(sline.Substring(60, 6)));
        memset(nmele, 0, sizeof(nmele));
          t = seek_first_word(sline.Substring(76, 2));
          if (t.length>0) memcpy(nmele, t.text, t.length);
    }
    void read_gro_line(StringNS::string sline){
        index = string_to_int(seek_first_word(sline.Substring(15, 5)));
        StringNS::string t;
        memset(name, 0, sizeof(name));
          t = seek_first_word(sline.Substring(10, 5));
          if (t.length>0) memcpy(name, t.text, t.length);
        memset(mole, 0, sizeof(mole));
          t = seek_first_word(sline.Substring(5, 5));
          if (t.length>0) memcpy(mole, t.text, t.length);
        ichain[0] = ichain[1] = 0;
        iaa = string_to_int(seek_first_word(sline.Substring(0, 5)));
        r.x      = string_to_double(seek_first_word(sline.Substring(20, 8)));
        r.y      = string_to_double(seek_first_word(sline.Substring(28, 8)));
        r.z      = string_to_double(seek_first_word(sline.Substring(36, 8)));
        occupacy = 0;
        Bfactor  = 0;
    }
    void export_pdb(char buffer[78]){
//ATOM      1  C1  FUL     1      10.270   8.200  31.670  1.00  0.00
        snprintf(buffer, 78, "ATOM  %5d %4s%4s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %s\n", \
          index%100000, name, mole, ichain, iaa, r.x*10, r.y*10, r.z*10, occupacy, Bfactor, nmele);
    }
    void export_gro(char buffer[45]){
        snprintf(buffer, 45, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", \
          iaa, mole, name, index, r.x, r.y, r.z);
    }
};
class PDBAtomSet {
  public:
    int      count;
    PDBAtom* atom;
    Vector   box;
};
