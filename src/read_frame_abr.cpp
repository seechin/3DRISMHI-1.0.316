// 2016.6.7: matrix boundary check
#ifndef TPAppendix
   class TPAppendix {
     public:
       double      time, prec;
       int         step;
   };
#endif
#ifndef MIN
  #define MIN(a,b) (a<b?a:b)
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// pdb reading
FILE* input_file_pdb=nullptr;
bool prepare_input_pdb(char* filename){
    if (!filename) return false;
    input_file_pdb = fopen(filename, "r"); if (!input_file_pdb) return false;
    fseek(input_file_pdb, 0, SEEK_SET); return true;
}
void end_input_pdb(){ if (input_file_pdb) fclose(input_file_pdb); }
int read_frame_pdb(PDBAtomSet* atoms, Vector * box, int * iframe){
    if (!input_file_pdb) return -1;
    int ia = 0; bool box_size_got = false; bool start_record = false;
    FILE* file = input_file_pdb; char input[4096];
    while (true){
        if (!fgets(input, sizeof(input), file)) return 0;
        int idx = 0; StringNS::string s = StringNS::seek_first_word(StringNS::string(input), &idx);
        if (StringNS::string(s.text, MIN(s.length, 4))=="ATOM" || StringNS::string(s.text, MIN(s.length, 6))=="HETATM"){
            PDBAtom atmp; atmp.read_pdb_line(StringNS::string(input));
            if (ia>=atoms->count) return -2;
            atoms->atom[ia++].r = atmp.r;
            //atoms->atom[ia++].read_pdb_line(StringNS::string(input));
        } else if (StringNS::string(s.text, MIN(s.length, 6))>="CRYST1"){
            box_size_got = true;
            StringNS::string swx = StringNS::seek_first_word(StringNS::string(input), &idx);
            StringNS::string swy = StringNS::seek_first_word(StringNS::string(input), &idx);
            StringNS::string swz = StringNS::seek_first_word(StringNS::string(input), &idx);
            box->x = string_to_double(swx) / 10.0;
            box->y = string_to_double(swy) / 10.0;
            box->z = string_to_double(swz) / 10.0;
        } else if (StringNS::string(s.text, MIN(s.length, 5))=="MODEL"){
            StringNS::string siframe = StringNS::seek_first_word(StringNS::string(input), &idx);
            *iframe = string_to_int(siframe);
            start_record = true;
        } else if (StringNS::string(s.text, MIN(s.length, 3))=="TER"){
            if (start_record){
                start_record = false; break;
            }
        } else if (StringNS::string(s.text, MIN(s.length, 6))=="ENDMDL"){
            if (start_record){
                start_record = false; break;
            }
        } else if (StringNS::string(s.text, MIN(s.length, 3))=="END"){
            if (start_record) break;
        }
    }
    if (ia<=0) return 0;
    if (ia!=atoms->count) return -2;
    return 1;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// gro reading
FILE* input_file_gro=nullptr;
bool prepare_input_gro(char* filename){
    if (!filename) return false;
    input_file_gro = fopen(filename, "r"); if (!input_file_gro) return false; else return true;
    fseek(input_file_gro, 0, SEEK_SET);
}
void end_input_gro(){ if (input_file_gro) fclose(input_file_gro); }
int read_frame_gro(PDBAtomSet* atoms, Vector * box, int * iframe){
    if (!input_file_gro) return -1;
    int ia = 0; bool box_size_got = false; bool start_record = false;
    FILE* file = input_file_gro; char input[4096];
  //2. read the leading 2 lines
    if (!fgets(input, sizeof(input), file)) return 0;
    if (!fgets(input, sizeof(input), file)) return 0;
    int count = string_to_int(StringNS::seek_first_word(StringNS::string(input)));
  //3. read the atoms
    for (int i=0; i<count; i++){
        if (!fgets(input, sizeof(input), file)) return 0;
        StringNS::string sline(input);
        if (sline.length<44) return -3;

        if (ia>=atoms->count) return -2;
        PDBAtom atmp; atmp.read_gro_line(sline);
        atoms->atom[ia++].r = atmp.r;
        //atoms->atom[i].read_gro_line(sline);
    }
  //4. read the box size
    if (!fgets(input, sizeof(input), file)) return 0;
    int idx = 0;
    StringNS::string swx = StringNS::seek_first_word(StringNS::string(input), &idx);
    StringNS::string swy = StringNS::seek_first_word(StringNS::string(input), &idx);
    StringNS::string swz = StringNS::seek_first_word(StringNS::string(input), &idx);
    box->x = string_to_double(swx);
    box->y = string_to_double(swy);
    box->z = string_to_double(swz);
  //x. end
    if (ia<=0) return 0;
    if (ia!=atoms->count) return -1;
    return 1;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// xtc reading
#ifdef _GROMACS_
  #ifdef _GROMACS4_
    #include    <gromacs/xtcio.h>
  #else
    #include    <gromacs/fileio/xtcio.h>
  #endif
    class XtcFileInfo {
      public:
        gmx_bool  valid;
        t_fileio* file;
        int     na, step;
        matrix  box;
        rvec    *x;
        real    time, prec;
        int     nframe;
    };
    XtcFileInfo input_file_xtc;
    bool prepare_input_xtc(char* filename){
        FILE * ft = fopen(filename, "r"); if (!ft) return false; fclose(ft);

        XtcFileInfo * x = &input_file_xtc;
        x->file = open_xtc(filename, "r");
        x->valid = true;
        x->nframe = 0;
        return x->valid;
    }
    void end_input_xtc(){
        if (input_file_xtc.valid){
            input_file_xtc.valid = false;
            close_xtc(input_file_xtc.file);
            input_file_xtc.file = 0;
            input_file_xtc.nframe = 0;
        }
    }
    int read_frame_xtc(PDBAtomSet* atoms, Vector * box, int * iframe, TPAppendix * tpa){
        XtcFileInfo * x = &input_file_xtc;

        if (!x->valid) return -1;

        bool ret = false; gmx_bool bok = false;
        if (x->nframe<=0){
            ret = read_first_xtc(x->file, &x->na, &x->step, &x->time, x->box, &x->x, &x->prec, &bok);
            x->nframe = 1;
        } else {
            for (int i=0; i<x->na && i<atoms->count; i++) x->x[i][0] = x->x[i][1] = x->x[i][2] = 0;
            ret = read_next_xtc (x->file, x->na, &x->step, &x->time, x->box, x->x, &x->prec, &bok);
            x->nframe ++;
        }
        *iframe = x->nframe;
        if (!ret || !bok) return 0;
        if (atoms->count != x->na) return -2;

        PDBAtom * a = atoms->atom;
        for (int i=0; i<x->na && i<atoms->count; i++){
            a[i].r.x = x->x[i][0];
            a[i].r.y = x->x[i][1];
            a[i].r.z = x->x[i][2];
        }
        box->x = x->box[0][0];
        box->y = x->box[1][1];
        box->z = x->box[2][2];
        tpa->step = x->step;
        tpa->time = x->time;
        tpa->prec = x->prec;
//printf("step %d, time %f, prec %f\n", x->step, x->time, x->prec);
        if (x->na!=atoms->count) return 0;
        return 1;
    }
#else
    bool prepare_input_xtc(char* filename){
        fprintf(stderr, "%s : error : xtc input not yet supported.\n", software_name);
        return false;
    }
    void end_input_xtc(){ }
    bool read_frame_xtc(PDBAtomSet* atoms, Vector * box, int * iframe, TPAppendix * tpa){ return false; }
#endif
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define FRAME_FMT_NUL   0
#define FRAME_FMT_PDB   1
#define FRAME_FMT_GRO   2
#define FRAME_FMT_XTC   3

int frame_format = FRAME_FMT_NUL;

//prepare for input
bool prepare_input(char* pdb, char* gro, char* xtc){
    if (pdb&&!pdb[0]) pdb = nullptr;
    if (gro&&!gro[0]) gro = nullptr;
    if (xtc&&!xtc[0]) xtc = nullptr;
    bool success = true;
    if (success){
        if (pdb){           success = prepare_input_pdb(pdb); frame_format = FRAME_FMT_PDB;
        } else if (gro){    success = prepare_input_gro(gro); frame_format = FRAME_FMT_GRO;
        } else if (xtc){    success = prepare_input_xtc(xtc); frame_format = FRAME_FMT_XTC;
        } else {
            return false;
        }
        if (!success){
            fprintf(stderr, "%s : \"-f\" : cannot open %s.\n", software_name, pdb?pdb:gro?gro:xtc);
            return false;
        }
    }
    return true;
}


int read_frame(PDBAtomSet * itp_atom_set, int * iframe, TPAppendix * tpa_=nullptr){
    TPAppendix tpa; int frame = 0;
    int ret = 0;
    if (frame_format==FRAME_FMT_PDB){           ret = read_frame_pdb(itp_atom_set, &itp_atom_set->box, &frame); (tpa_?tpa_:&tpa)->time = 1 + *iframe;
    } else if (frame_format==FRAME_FMT_GRO){    ret = read_frame_gro(itp_atom_set, &itp_atom_set->box, &frame); (tpa_?tpa_:&tpa)->time = 1 + *iframe;
    } else if (frame_format==FRAME_FMT_XTC){    ret = read_frame_xtc(itp_atom_set, &itp_atom_set->box, &frame, tpa_?tpa_:&tpa);
    } else ret = -1;
    (*iframe) ++;
    return ret;
}
