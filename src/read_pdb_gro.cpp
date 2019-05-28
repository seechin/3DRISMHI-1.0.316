bool load_pdb(PDBAtomSet* atom, char* pdb){
    if (!pdb || !atom) return false;
    bool start_record = true; bool box_size_got = false; char input[4096];
  //1. open file
    FILE * file = fopen(pdb, "r"); if (!file) return false;
    if (atom->atom) free(atom->atom); atom->count = 0;
  //2. read through to find the atom number, box size, etc
    fseek(file, 0, SEEK_SET);
    start_record = true;
    while (fgets(input, sizeof(input), file)){
        int idx = 0; StringNS::string s = StringNS::seek_first_word(StringNS::string(input), &idx);
        if (s=="ATOM"){ if (start_record) atom->count ++;
        } else if (s>="CRYST1"){
            box_size_got = true;
            StringNS::string swx = StringNS::seek_first_word(StringNS::string(input), &idx);
            StringNS::string swy = StringNS::seek_first_word(StringNS::string(input), &idx);
            StringNS::string swz = StringNS::seek_first_word(StringNS::string(input), &idx);
            atom->box.x = string_to_double(swx) / 10.0;
            atom->box.y = string_to_double(swy) / 10.0;
            atom->box.z = string_to_double(swz) / 10.0;
        } else if (s=="MODEL"){  start_record = true;
        } else if (s=="TER"){    start_record = false; break;
        } else if (s=="ENDMDL"){ start_record = false; break;
        } else if (s=="END"){    start_record = false; break;
        }
    }
    if (!box_size_got) return false;
  //3. allocate memory, and read through to read the atoms
    atom->atom = (PDBAtom*)malloc(sizeof(PDBAtom)*atom->count);
    if (!atom->atom){ atom->count = 0; fclose(file); return false; }
    fseek(file, 0, SEEK_SET);
    int ia = 0; start_record = false;
    while (fgets(input, sizeof(input), file)){
        int idx = 0; StringNS::string s = StringNS::seek_first_word(StringNS::string(input), &idx);
        if (s=="ATOM"){
            atom->atom[ia].read_pdb_line(StringNS::string(input));
            ia++;
        } else if (s=="MODEL"){  start_record = true;
        } else if (s=="TER"){    break;
        } else if (s=="ENDMDL"){ start_record = false; break;
        } else if (s=="END"){    break;
        }
    }
  //x. end
    if (file) fclose(file); //printf("%d atoms. box: %f %f %f\n", atom->count, atom->box.x, atom->box.y, atom->box.z);
    return true;
}

bool load_gro(PDBAtomSet* atom, char* gro){
    if (!gro || !atom) return false;
    char input[4096];
  //1. open file
    FILE * file = fopen(gro, "r"); if (!file) return false;
    if (atom->atom) free(atom->atom); atom->count = 0;
  //2. read the leading 2 lines
    fseek(file, 0, SEEK_SET);
    if (!fgets(input, sizeof(input), file)) return false;
    if (!fgets(input, sizeof(input), file)) return false;
    atom->count = string_to_int(StringNS::seek_first_word(StringNS::string(input)));
    atom->atom = (PDBAtom*)malloc(sizeof(PDBAtom)*atom->count);
    if (!atom->atom){ atom->count = 0; fclose(file); return false; }
  //3. read the atoms
    for (int i=0; i<atom->count; i++){
        if (!fgets(input, sizeof(input), file)) return false;
        StringNS::string sline(input);
        if (sline.length<44) return false;
        atom->atom[i].read_gro_line(sline);
    }
  //4. read the box size
    if (!fgets(input, sizeof(input), file)) return false;
    int idx = 0;
    StringNS::string swx = StringNS::seek_first_word(StringNS::string(input), &idx);
    StringNS::string swy = StringNS::seek_first_word(StringNS::string(input), &idx);
    StringNS::string swz = StringNS::seek_first_word(StringNS::string(input), &idx);
    atom->box.x = string_to_double(swx);
    atom->box.y = string_to_double(swy);
    atom->box.z = string_to_double(swz);

  //x. end
    if (file) fclose(file); //printf("%d atoms. box: %f %f %f\n", atom->count, atom->box.x, atom->box.y, atom->box.z);
    return true;
}

bool read_system_file(PDBAtomSet* atom, char* pdb, char* gro){
    if (pdb&&!pdb[0]) pdb = nullptr;
    if (gro&&!gro[0]) gro = nullptr;
    bool atom_handled = false;
    if (!atom_handled && (pdb || gro)){
        if (pdb){
            atom_handled = load_pdb(atom, pdb);
        } else if (gro){
            atom_handled = load_gro(atom, gro);
        }
    } else atom_handled = true;

    if (!atom_handled){
        fprintf(stderr, "%s : \"-s\" : no atom is loaded.\n", software_name);
        return false;
    } else {
        return true;
    }
}
