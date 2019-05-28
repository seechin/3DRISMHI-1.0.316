//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // cat alphabet
    //A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
    // cat alphabet_l
    //a   b   c   d   e   f   g   h   i   j   k   l   m   n   o   p   q   r   s   t   u   v   w   x   y   z
    // generate of single letter elements
    // for A in `cat alphabet`; do Z=`cat ~/workspace/all-elements.txt|awk -v e=$A '{if($1==e){print $2}}'`; printf "%s %4d\n" $A $Z; done | awk '{printf("%d,",$2)}'
#ifndef __ElementXNS__
#define __ElementXNS__
namespace ElementNS {
    //                           A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z
    const int atom_tab_1[] = {   0,   5,   6,   0,   0,   9,   0,   1,  53,   0,  19,   0,   0,   7,   8,  15,   0,   0,  16,   0,  92,  23,  74,   0,  39,   0 };
    const int atom_tab_a[] = {   0,   0,  89,   0,   0,   0,  47,   0,   0,   0,   0,  13,  95,   0,   0,   0,   0,  18,  33,  85,  79,   0,   0,   0,   0,   0 };
    const int atom_tab_b[] = {  56,   0,   0,   0,   4,   0,   0, 107,  83,   0,  97,   0,   0,   0,   0,   0,   0,  35,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_c[] = {  20,   0,   0,  48,  58,  98,   0,   0,   0,   0,   0,  17,  96, 112,  27,   0,   0,  24,  55,   0,  29,   0,   0,   0,   0,   0 };
    const int atom_tab_d[] = {   0, 105,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 110,   0,   0,   0,   0,   0,  66,   0 };
    const int atom_tab_e[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  68,  99,   0,  63,   0,   0,   0,   0,   0 };
    const int atom_tab_f[] = {   0,   0,   0,   0,  26,   0,   0,   0,   0,   0,   0, 114, 100,   0,   0,   0,   0,  87,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_g[] = {  31,   0,   0,  64,  32,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_h[] = {   0,   0,   0,   0,   2,  72,  80,   0,   0,   0,   0,   0,   0,   0,  67,   0,   0,   0, 108,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_i[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  49,   0,   0,   0,  77,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_j[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_k[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  36,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_l[] = {  57,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0, 103,   0,   0,  71, 116,   0,   0,   0,   0 };
    const int atom_tab_m[] = {   0,   0, 115, 101,   0,   0,  12,   0,   0,   0,   0,   0,   0,  25,  42,   0,   0,   0,   0, 109,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_n[] = {  11,  41,   0,  60,  10,   0,   0, 113,  28,   0,   0,   0,   0,   0, 102,  93,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_o[] = {   0,   0,   0,   0,   0,   0, 118,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  76,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_p[] = {  91,  82,   0,  46,   0,   0,   0,   0,   0,   0,   0,   0,  61,   0,  84,   0,   0,  59,   0,  78,  94,   0,   0,   0,   0,   0 };
    const int atom_tab_q[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_r[] = {  88,  37,   0,   0,  75, 104, 111,  45,   0,   0,   0,   0,   0,  86,   0,   0,   0,   0,   0,   0,  44,   0,   0,   0,   0,   0 };
    const int atom_tab_s[] = {   0,  51,  21,   0,  34,   0, 106,   0,  14,   0,   0,   0,  62,  50,   0,   0,   0,  38,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_t[] = {  73,  65,  43,   0,  52,   0,   0,  90,  22,   0,   0,  81,  69,   0,   0,   0,   0,   0, 117,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_u[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_v[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_w[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_x[] = {   0,   0,   0,   0,  54,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_y[] = {   0,  70,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int atom_tab_z[] = {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  30,   0,   0,   0,  40,   0,   0,   0,   0,   0,   0,   0,   0 };
    const int * atom_tab_2[] = { atom_tab_a, atom_tab_b, atom_tab_c, atom_tab_d, atom_tab_e, atom_tab_f, atom_tab_g, atom_tab_h, atom_tab_i, atom_tab_j, atom_tab_k, atom_tab_l, atom_tab_m, atom_tab_n, atom_tab_o, atom_tab_p, atom_tab_q, atom_tab_r, atom_tab_s, atom_tab_t, atom_tab_u, atom_tab_v, atom_tab_w, atom_tab_x, atom_tab_y, atom_tab_z };
    class Element {
      public:
        const char * name;
        int an/*atom number*/;
        double m/*weight*/;
    };
    Element elements[] = {
          {    "H" ,     1,          1.008 },
          {    "He",     2,      4.0026022 },
          {    "Li",     3,           6.94 },
          {    "Be",     4,     9.01218315 },
          {    "B" ,     5,          10.81 },
          {    "C" ,     6,         12.011 },
          {    "N" ,     7,         14.007 },
          {    "O" ,     8,         15.999 },
          {    "F" ,     9,    18.99840316 },
          {    "Ne",    10,       20.17976 },
          {    "Na",    11,    22.98976928 },
          {    "Mg",    12,         24.305 },
          {    "Al",    13,    26.98153857 },
          {    "Si",    14,         28.085 },
          {    "P" ,    15,      30.973762 },
          {    "S" ,    16,          32.06 },
          {    "Cl",    17,          35.45 },
          {    "Ar",    18,        39.9481 },
          {    "K" ,    19,       39.09831 },
          {    "Ca",    20,        40.0784 },
          {    "Sc",    21,     44.9559085 },
          {    "Ti",    22,        47.8671 },
          {    "V" ,    23,       50.94151 },
          {    "Cr",    24,       51.99616 },
          {    "Mn",    25,     54.9380443 },
          {    "Fe",    26,        55.8452 },
          {    "Co",    27,     58.9331944 },
          {    "Ni",    28,       58.69344 },
          {    "Cu",    29,        63.5463 },
          {    "Zn",    30,         65.382 },
          {    "Ga",    31,        69.7231 },
          {    "Ge",    32,        72.6308 },
          {    "As",    33,     74.9215956 },
          {    "Se",    34,        78.9718 },
          {    "Br",    35,         79.904 },
          {    "Kr",    36,        83.7982 },
          {    "Rb",    37,       85.46783 },
          {    "Sr",    38,         87.621 },
          {    "Y" ,    39,      88.905842 },
          {    "Zr",    40,        91.2242 },
          {    "Nb",    41,      92.906372 },
          {    "Mo",    42,         95.951 },
          {    "Tc",    43,             98 },
          {    "Ru",    44,        101.072 },
          {    "Rh",    45,     102.905502 },
          {    "Pd",    46,        106.421 },
          {    "Ag",    47,      107.86822 },
          {    "Cd",    48,       112.4144 },
          {    "In",    49,       114.8181 },
          {    "Sn",    50,       118.7107 },
          {    "Sb",    51,       121.7601 },
          {    "Te",    52,        127.603 },
          {    "I" ,    53,     126.904473 },
          {    "Xe",    54,       131.2936 },
          {    "Cs",    55,     132.905452 },
          {    "Ba",    56,       137.3277 },
          {    "La",    57,     138.905477 },
          {    "Ce",    58,       140.1161 },
          {    "Pr",    59,     140.907662 },
          {    "Nd",    60,       144.2423 },
          {    "Pm",    61,            145 },
          {    "Sm",    62,        150.362 },
          {    "Eu",    63,       151.9641 },
          {    "Gd",    64,        157.253 },
          {    "Tb",    65,     158.925352 },
          {    "Dy",    66,       162.5001 },
          {    "Ho",    67,     164.930332 },
          {    "Er",    68,       167.2593 },
          {    "Tm",    69,     168.934222 },
          {    "Yb",    70,       173.0451 },
          {    "Lu",    71,      174.96681 },
          {    "Hf",    72,        178.492 },
          {    "Ta",    73,     180.947882 },
          {    "W" ,    74,        183.841 },
          {    "Re",    75,       186.2071 },
          {    "Os",    76,        190.233 },
          {    "Ir",    77,       192.2173 },
          {    "Pt",    78,       195.0849 },
          {    "Au",    79,    196.9665695 },
          {    "Hg",    80,       200.5923 },
          {    "Tl",    81,         204.38 },
          {    "Pb",    82,         207.21 },
          {    "Bi",    83,     208.980401 },
          {    "Po",    84,            209 },
          {    "At",    85,            210 },
          {    "Rn",    86,            222 },
          {    "Fr",    87,            223 },
          {    "Ra",    88,            226 },
          {    "Ac",    89,            227 },
          {    "Th",    90,      232.03774 },
          {    "Pa",    91,     231.035882 },
          {    "U" ,    92,     238.028913 },
          {    "Np",    93,            237 },
          {    "Pu",    94,            244 },
          {    "Am",    95,            243 },
          {    "Cm",    96,            247 },
          {    "Bk",    97,            247 },
          {    "Cf",    98,            251 },
          {    "Es",    99,            252 },
          {    "Fm",   100,            257 },
          {    "Md",   101,            258 },
          {    "No",   102,            259 },
          {    "Lr",   103,            266 },
          {    "Rf",   104,            267 },
          {    "Db",   105,            268 },
          {    "Sg",   106,            269 },
          {    "Bh",   107,            270 },
          {    "Hs",   108,            277 },
          {    "Mt",   109,            278 },
          {    "Ds",   110,            281 },
          {    "Rg",   111,            282 },
          {    "Cn",   112,            285 },
          {    "Nh",   113,            286 },
          {    "Fl",   114,            289 },
          {    "Mc",   115,            290 },
          {    "Lv",   116,            293 },
          {    "Ts",   117,            294 },
          {    "Og",   118,            294 }
    };
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Element * get_element_parameter(const char * name){
        int Z = 0;
        if (name[0]>='A'&&name[0]<='Z'){
            if (name[1]>='a'&&name[1]<='z') Z = atom_tab_2[name[0]-'A'][name[1]-'a'];
            if (Z<=0) Z = atom_tab_1[name[0]-'A'];
        }
        if (Z-1>=0 && Z-1<sizeof(elements)/sizeof(elements[0])) return &elements[Z-1];
        return nullptr;
    }
    const char * get_atom_element(const char * name/*>=2bytes*/){
        Element * e = get_element_parameter(name);
        if (e) return e->name; else return "";


        int Z = 0;
        if (name[0]>='A'&&name[0]<='Z'){
            if (name[1]>='a'&&name[1]<='z') Z = atom_tab_2[name[0]-'A'][name[1]-'a'];
            if (Z<=0) Z = atom_tab_1[name[0]-'A'];
        }
        if (Z-1>=0 && Z-1<sizeof(elements)/sizeof(elements[0])) return elements[Z-1].name;
        return "";
    }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//    Element * get_element_parameter_old(const char * elem){
//        for (int i=0; i<sizeof(elements)/sizeof(elements[0]); i++){
//            if (elem[0]==elements[i].name[0] && elem[1]==elements[i].name[1]){
//                return &elements[i];
//            }
//        }
//        return nullptr;
//    }
//    const char * get_atom_element_old(const char * name/*>=2bytes*/){
//        switch (name[0]){
//          case 'A': if (name[1]=='l') return "Al"; else if (name[1]=='r') return "Ar"; else if (name[1]=='s') return "As";
//               else if (name[1]=='g') return "Ag"; else if (name[1]=='u') return "Au"; else if (name[1]=='t') return "Au";
//               else if (name[1]=='c') return "Ac"; else if (name[1]=='m') return "Am"; break;
//          case 'B': if (name[1]=='e') return "Be"; else if (name[1]=='r') return "Br"; else if (name[1]=='a') return "Ba";
//               else if (name[1]=='i') return "Bi"; else if (name[1]=='k') return "Bk"; else if (name[1]=='b') return "Bb";
//               else return "B"; break;
//          case 'C': if (name[1]=='l') return "Cl"; else if (name[1]=='a') return "Ca"; else if (name[1]=='r') return "Cr";
//               else if (name[1]=='o') return "Co"; else if (name[1]=='u') return "Cu"; else if (name[1]=='d') return "Cd";
//               else if (name[1]=='s') return "Cs"; else if (name[1]=='e') return "Ce"; else if (name[1]=='m') return "Cm";
//               else if (name[1]=='f') return "Cf"; else if (name[1]=='n') return "Cn"; else return "C"; break;
//          case 'D': if (name[1]=='y') return "Dy"; else if (name[1]=='b') return "Db"; else if (name[1]=='s') return "Ds";
//               break;
//          case 'E': if (name[1]=='u') return "Eu"; else if (name[1]=='r') return "Er"; else if (name[1]=='s') return "Es";
//               break;
//          case 'F': if (name[1]=='e') return "Fe"; else if (name[1]=='r') return "Fr"; else if (name[1]=='m') return "Fm";
//               else if (name[1]=='l') return "Fl"; else return "F"; break;
//          case 'G': if (name[1]=='a') return "Ga"; else if (name[1]=='e') return "Ge"; else if (name[1]=='d') return "Gd";
//               break;
//          case 'H': if (name[1]=='e') return "He"; else if (name[1]=='o') return "Ho"; else if (name[1]=='f') return "Hf";
//               else if (name[1]=='g') return "Hg"; else if (name[1]=='s') return "Hs"; else return "H"; break;
//          case 'I': if (name[1]=='n') return "In"; else if (name[1]=='r') return "Ir"; else return "I"; break;
//          case 'K': if (name[1]=='r') return "Kr"; else return "K"; break;
//          case 'L': if (name[1]=='i') return "Li"; else if (name[1]=='a') return "La"; else if (name[1]=='u') return "Lu";
//               else if (name[1]=='r') return "Lr"; else if (name[1]=='v') return "Lv"; break;
//          case 'M': if (name[1]=='g') return "Mg"; else if (name[1]=='n') return "Mn"; else if (name[1]=='o') return "Mo";
//               else if (name[1]=='d') return "Md"; else if (name[1]=='t') return "Mt"; else if (name[1]=='c') return "Mc";
//               break;
//          case 'N': if (name[1]=='e') return "Ne"; else if (name[1]=='a') return "Na"; else if (name[1]=='i') return "Ni";
//               else if (name[1]=='b') return "Nb"; else if (name[1]=='d') return "Nd"; else if (name[1]=='p') return "Np";
//               else if (name[1]=='o') return "No"; else if (name[1]=='h') return "Nh"; else return "N"; break;
//          case 'O': if (name[1]=='s') return "Os"; else if (name[1]=='g') return "Og"; else return "O"; break;
//          case 'P': if (name[1]=='d') return "Pd"; else if (name[1]=='r') return "Pr"; else if (name[1]=='m') return "Pm";
//               else if (name[1]=='t') return "Pt"; else if (name[1]=='b') return "Pb"; else if (name[1]=='o') return "Po";
//               else if (name[1]=='a') return "Pa"; else if (name[1]=='u') return "Pu"; else return "P"; break;
//          case 'R': if (name[1]=='b') return "Rb"; else if (name[1]=='u') return "Ru"; else if (name[1]=='h') return "Rh";
//               else if (name[1]=='e') return "Re"; else if (name[1]=='n') return "Rn"; else if (name[1]=='a') return "Ra";
//               else if (name[1]=='f') return "Rf"; else if (name[1]=='g') return "Rg"; break;
//          case 'S': if (name[1]=='i') return "Si"; else if (name[1]=='c') return "Sc"; else if (name[1]=='e') return "Se";
//               else if (name[1]=='r') return "Sr"; else if (name[1]=='n') return "Sn"; else if (name[1]=='b') return "Sb";
//               else if (name[1]=='m') return "Sm"; else if (name[1]=='g') return "Sg"; else return "S"; break;
//          case 'T': if (name[1]=='i') return "Ti"; else if (name[1]=='c') return "Tc"; else if (name[1]=='e') return "Te";
//               else if (name[1]=='b') return "Tb"; else if (name[1]=='m') return "Tm"; else if (name[1]=='a') return "Ta";
//               else if (name[1]=='l') return "Tl"; else if (name[1]=='h') return "Th"; else if (name[1]=='s') return "Ts";
//               break;
//          case 'U': return "U"; break;
//          case 'V': return "V"; break;
//          case 'W': return "W"; break;
//          case 'X': if (name[1]=='e') return "Xe"; break;
//          case 'Y': if (name[1]=='b') return "Yb"; return "Y"; break;
//          case 'Z': if (name[1]=='n') return "Zn"; else if (name[1]=='r') return "Zr"; break;
//        }
//        return "";
//    }
}
#endif
