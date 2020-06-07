//StringX ============================================
// fix: 2018.02.26 : is_string_number remove unused variable v
// version : 2017.01.21 : ToLong, ToInt
// version : 2015.05.05 : seek_first_general_word, analysis_general_line
// version : 2015.04.27 : seek_first_csv_word, analysis_csv_line
// version : 2018.02.26 : fix size_t conversion warning
// version : 2013.10.17 : is_string_number, seek_first_word, analysis_line
// version : 2011.01.12 : basic functions
//StringX ============================================
#ifndef nullptr
    #define nullptr 0
#endif

#ifndef __StringNS__
#define __StringNS__

//string class
namespace StringNS {
    class string {
      public:
        char * text;
        int    length;
      public:
        string(){ text = nullptr; length = 0; }
        string(const char * lpsz){ text = (char*)lpsz; length = (int)strlen(text); }
        string(char * lpsz){ text = lpsz; length = (int)strlen(text); }
        string(char * lpsz, int len){ text = lpsz; length = len; }
        void FromChars(char * lpsz){ text = lpsz; length = (int)strlen(text); }
        void FromChars(char * lpsz, int len){ text = lpsz; length = len; }
        inline char operator [] (int i){ return text[i]; }
        string Substring(int s, int l){ return string(&text[s], l); }
        string sub(int s, int l){ return string(s>=length?&text[length]:&text[s], l>length-s?length-s:l); }
        string sub(int s){ return string(s>=length?&text[length]:&text[s], s>length?0:length-s); }
        static string * Alloc(int len){
          #ifdef _WINDOWS_
            ::LPVOID lpvoid = ::GlobalAlloc(GPTR, len + sizeof(string) + 4);
          #else
            void * lpvoid = malloc(len + sizeof(string) + 4);
          #endif
            ((string *)lpvoid)->text   = (char *)lpvoid + sizeof(string);
            ((string *)lpvoid)->length = len;
            return (string *)lpvoid;
        }
        void Dispose(){
          #ifdef _WINDOWS_
            ::GlobalFree(this);
          #else
            free(this);
          #endif
        }
        unsigned long GetHashCode(int hl = 0){
            unsigned long u = 0; unsigned long crc = 16777213;
            for (int i = 0; i < length; i ++){
                unsigned char c = text[i]; if (c >= 'a' && c <= 'z') c += 'A' - 'a';
                u = (u * 256 + c) % crc;
            }
            if (hl) u %= hl; return (long)u;
        }
        long ToLong(){
            char long_buffer[48]; memset(long_buffer, 0, sizeof(long_buffer));
            memcpy(long_buffer, text, length>sizeof(long_buffer)-1? sizeof(long_buffer)-1 : length);
            return atof(long_buffer);
        }
        int ToInt(){
            char int_buffer[48]; memset(int_buffer, 0, sizeof(int_buffer));
            memcpy(int_buffer, text, length>sizeof(int_buffer)-1? sizeof(int_buffer)-1 : length);
            return atoi(int_buffer);
        }
        double ToDouble();
        bool Equ(string s){
                if (s.length != length) return false;
                for (int i = 0; i < length; i++) if (s.text[i] != text[i]) return false;
                return true;
        }
        int CompareIgnoreCase(string s){
                int i; for (i = 0; i < length || i < s.length; i ++){
                        if (i >= length) return -1; else if (i >= s.length) return 1;
                        unsigned char c1 = text[i]; unsigned char c2 = s.text[i];
                        if (c1 <= 'z' && c1 >= 'a') c1 += 'A' - 'a';
                        if (c2 <= 'z' && c2 >= 'a') c2 += 'A' - 'a';
                        if (c1 > c2) return 2; else if (c1 < c2) return -2;
                }
                return 0;
        }
        int Compare(string s){
                int i; for (i = 0; i < length || i < s.length; i ++){
                        if (i >= length) return -1; else if (i >= s.length) return 1;
                        if (text[i] > s.text[i]) return 1; else if (text[i] < s.text[i]) return -1;
                }
                return 0;
        }
        int operator - (string s){
            return CompareIgnoreCase(s);
        }
        bool operator == (string s){
            return CompareIgnoreCase(s) == 0;
        }
        bool operator == (char * str){
            return CompareIgnoreCase(string(str)) == 0;
        }
        bool operator == (const char * str){
            return CompareIgnoreCase(string(str)) == 0;
        }
        bool operator >= (string s){
            int a = CompareIgnoreCase(s); return a==0||a==1;
        }
        bool operator >= (char * str){
            int a = CompareIgnoreCase(string(str)); return a==0||a==1;
        }
        bool operator >= (const char * str){
            int a = CompareIgnoreCase(string(str)); return a==0||a==1;
        }
        bool operator > (string s){
            int a = CompareIgnoreCase(s); return a==1;
        }
        bool operator > (char * str){
            int a = CompareIgnoreCase(string(str)); return a==1;
        }
        bool operator > (const char * str){
            int a = CompareIgnoreCase(string(str)); return a==1;
        }
        bool operator <= (string s){
            int a = CompareIgnoreCase(s); return a==0||a==-1;
        }
        bool operator <= (char * str){
            int a = CompareIgnoreCase(string(str)); return a==0||a==-1;
        }
        bool operator <= (const char * str){
            int a = CompareIgnoreCase(string(str)); return a==0||a==-1;
        }
        bool operator < (string s){
            int a = CompareIgnoreCase(s); return a==-1;
        }
        bool operator < (char * str){
            int a = CompareIgnoreCase(string(str)); return a==-1;
        }
        bool operator < (const char * str){
            int a = CompareIgnoreCase(string(str)); return a==-1;
        }

        bool operator != (string s){
            return CompareIgnoreCase(s) != 0;
        }
        bool operator != (char * str){
            return CompareIgnoreCase(string(str)) != 0;
        }
        bool operator != (const char * str){
            return CompareIgnoreCase(string(str)) != 0;
        }
        void KMPnext(int * next){
            if (!next) return ;
            next[0] = -1; int i = 0; int j = -1;
            while (i < length){
                if (j == -1 || text[i] == text[j]){ ++i; ++j; next[i] = j; }
                else j = next[j];
            }
        }
        int  KMPindex(string t, int * next){
            if (!next) return -1;
            int i = 0; int j = 0;
            while (i < t.length && j < length){
                if (j == -1 || t.text[i] == text[j]){ ++ i; ++ j; }
                else j = next[j];
            }
            if (j >= length) return i - length;
            else return -1;
        }
        void KMPnext_ignore_case(int * next){
            if (!next) return ;
            next[0] = -1; int i = 0; int j = -1;
            while (i < length){
                char c1 = text[i]; if (c1 >= 'a' && c1 <= 'z') c1 += 'A' - 'a';
                char c2 = text[j]; if (c2 >= 'a' && c2 <= 'z') c2 += 'A' - 'a';
                if (j == -1 || c1 == c2){ ++i; ++j; next[i] = j; }
                else j = next[j];
            }
        }
        int  KMPindex_ignore_case(string t, int * next){
            if (!next) return -1;
            int i = 0; int j = 0;
            while (i < t.length && j < length){
                char c1 = t.text[i]; if (c1 >= 'a' && c1 <= 'z') c1 += 'A' - 'a';
                char c2 = text[j]; if (c2 >= 'a' && c2 <= 'z') c2 += 'A' - 'a';
                if (j == -1 || c1 == c2){ ++ i; ++ j; }
                else j = next[j];
            }
            if (j >= length) return i - length;
            else return -1;
        }
    };
    class String : public string {
      public:
        String() : string() { }
        String(char * lpsz) : string(lpsz) { }
        String(char * lpsz, int len) : string(lpsz, len) { }
    };

    bool string_cmp_routine(char * s1, char * s2){ //return true;
        int i = 0;
        while (true){
            unsigned char c1 = s1[i]; unsigned char c2 = s2[i];
            if (c1 >= 'a' && c1 <= 'z') c1 += 'A' - 'a';
            if (c2 >= 'a' && c2 <= 'z') c2 += 'A' - 'a';
            if (c1 == c2){ if (c1 == 0) return true; else i ++; }
            else return false;
        }
    }
}

//string function
namespace StringNS {
    int string_resemble_length(char * s1, char * s2){ //return true;
        int i = 0;
        while (true){
            unsigned char c1 = s1[i]; unsigned char c2 = s2[i];
            if (c1 >= 'a' && c1 <= 'z') c1 += 'A' - 'a';
            if (c2 >= 'a' && c2 <= 'z') c2 += 'A' - 'a';
            if (c1 == c2){ if (c1 == 0) return i; else i ++; }
            else return i;
        }
    }
    string seek_first_word(string src, int * idx=nullptr){
        int i = idx==nullptr? 0 : *idx; int endl = i;
        while (i < src.length && src.text[i] && (src.text[i]==' '||src.text[i]=='\t'||src.text[i]=='\r'||src.text[i]=='\n')) i++;
        char c_skip[2]; c_skip[0] = ' '; c_skip[1] = '\t';
        if (src.text[i]=='\'' || src.text[i] == '"'){ c_skip[0] = c_skip[1] = src.text[i]; i ++; }
        int e = i;
        while (e < src.length && src.text[e] && !(src.text[e]==c_skip[0]||src.text[e]==c_skip[1]||src.text[e]=='\r'||src.text[e]=='\n')) e++;
        if (src.text[e] == '"' || src.text[e] == '\''){
            endl = e + 1;
        } else {
            endl = e;
        }
        if (idx) *idx = endl;
        return string(&src.text[i], e - i);
    }
    bool is_string_number(string str){
        char * stop = nullptr;
        strtod(str.text, &stop);
        if (stop && stop-str.text<str.length) return false; else return true;
    }
    int analysis_line(string sline, string * sl, int maxnw, bool separate_szstr = false){
        int nw = 0; int idx = 0;
        while (nw<maxnw){ sl[nw] = seek_first_word(sline, &idx); if (sl[nw].length<=0) break; nw++; }
        if (separate_szstr) for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = 0;
        return nw;
    }
}

namespace StringNS {
    StringNS::string seek_first_csv_word(StringNS::string src, int * idx=nullptr, bool trim_spaces=false){
        int i = idx==nullptr? 0 : *idx; int endl = i;
        while (i < src.length && src.text[i] && (src.text[i]=='\r'||src.text[i]=='\n')) i++;
        int e = i; int in_quote = 0;
        while (e < src.length && src.text[e]){
            if (src.text[e]==',' || src.text[e]=='\r'||src.text[e]=='\n'){
                if (!in_quote) break;
            } else if (src.text[e]=='\''){  in_quote ^= 1;
            } else if (src.text[e]=='"'){   in_quote ^= 2;
            } else if (src.text[e]=='\\'){
                if (e+1<src.length) e++;
            }
            e++;
        }
        endl = e;
        if (endl < src.length && src.text[endl]==',') endl++;
        if (idx) *idx = endl;
        if (trim_spaces){
            while (i<e && (src.text[i]==' '||src.text[i]=='\t')) i++;
            while (e>i && (src.text[e-1]==' '||src.text[e-1]=='\t')) e--;
        }
        return StringNS::string(&src.text[i], e - i);
    }
    int analysis_csv_line(StringNS::string sline, StringNS::string * sl, int maxnw, bool separate_szstr = false, bool trim_spaces=false){
        while (sline.length>0 && (sline.text[sline.length-1]=='\n' || sline.text[sline.length-1]=='\r')) sline.length --;
        int nw = 0; int idx = 0;
        while (nw<maxnw){
            if (idx>=sline.length) break;
            sl[nw] = seek_first_csv_word(sline, &idx, trim_spaces);
            nw++;
        }
        if (separate_szstr) for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = 0;
        return nw;
    }

}

namespace StringNS {
    StringNS::string seek_first_wordline_word(StringNS::string src, int * idx=nullptr){
        int i = idx==nullptr? 0 : *idx; int endl = i;
        while (i < src.length && src.text[i] && (src.text[i]=='\r'||src.text[i]=='\n' || src.text[i]==' '||src.text[i]=='\t')) i++;
        int e = i; int in_quote = 0; int break_flag = 0; bool word_encounter = false; bool space_encounter = false;
        while (e < src.length && src.text[e]){
            switch (src.text[e]) {
                case '\r': case '\n':
                    break_flag = 3; break;
                case '\t': case ' ':
                    space_encounter = true; break;
                case '!': case '@': case '#': case '$': case '%': case '^': case '&': case '*': case '(': case ')': case '=': case '+':
                case '[': case '{': case ']': case '}': case '|': case '<': case '>': case '?': case '/': case '~':
                case ',': case ':': case ';': // note: '-' and '_' are seen as part of text
                    if (!in_quote) break_flag = 1; break;
                case '\'': in_quote ^= 1; break;
                case '"': in_quote ^= 2; break;
                case '`': in_quote ^= 2; break;
                case '\\': if (e+1<src.length) e++; break;
                default: if (space_encounter && !in_quote) break_flag = 2; break;
            }
            if (break_flag) break; else e++;
        }
        endl = e;
        if (endl < src.length){
            if (break_flag==1) endl++;
        }
        if (idx) *idx = endl;
          while (i<e && (src.text[i]==' '||src.text[i]=='\t')) i++;
          while (e>i && (src.text[e-1]==' '||src.text[e-1]=='\t')) e--;
        return StringNS::string(&src.text[i], e - i);
    }
    int analysis_word_line(StringNS::string sline, StringNS::string * sl, int maxnw, bool separate_szstr = false){
        while (sline.length>0 && (sline.text[sline.length-1]=='\n' || sline.text[sline.length-1]=='\r')) sline.length --;
        int nw = 0; int idx = 0;
        while (nw<maxnw){
            if (idx>=sline.length) break;
            sl[nw] = seek_first_wordline_word(sline, &idx);
            nw++;
        }
        if (separate_szstr) for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = 0;
        return nw;
    }

}

namespace StringNS {
    StringNS::string seek_first_general_word(const char sep[4], StringNS::string src, int * idx=nullptr, bool trim_spaces=false){
        int i = idx==nullptr? 0 : *idx; int endl = i;
        while (i < src.length && src.text[i] && (src.text[i]=='\r'||src.text[i]=='\n')) i++;
        int e = i; int in_quote = 0;
        while (e < src.length && src.text[e]){
            if (src.text[e]==sep[0] || src.text[e]==sep[1] || src.text[e]==sep[2] || src.text[e]==sep[3] || src.text[e]=='\r'||src.text[e]=='\n'){
                if (!in_quote) break;
            } else if (src.text[e]=='\''){  in_quote ^= 1;
            } else if (src.text[e]=='"'){   in_quote ^= 2;
            } else if (src.text[e]=='\\'){
                if (e+1<src.length) e++;
            }
            e++;
        }
        endl = e;
        if (endl < src.length && (src.text[endl]==sep[0] || src.text[endl]==sep[1] || src.text[endl]==sep[2] || src.text[endl]==sep[3])) endl++;
        if (idx) *idx = endl;
        if (trim_spaces){
            while (i<e && (src.text[i]==' '||src.text[i]=='\t')) i++;
            while (e>i && (src.text[e-1]==' '||src.text[e-1]=='\t')) e--;
        }
        return StringNS::string(&src.text[i], e - i);
    }
    int analysis_general_line(char sep[4], StringNS::string sline, StringNS::string * sl, int maxnw, bool separate_szstr = false, bool trim_spaces=false){
        while (sline.length>0 && (sline.text[sline.length-1]=='\n' || sline.text[sline.length-1]=='\r')) sline.length --;
        int nw = 0; int idx = 0;
        while (nw<maxnw){
            if (idx>=sline.length) break;
            sl[nw] = seek_first_general_word(sep, sline, &idx, trim_spaces);
            nw++;
        }
        if (separate_szstr) for (int i=0; i<nw; i++) sl[i].text[sl[i].length] = 0;
        return nw;
    }

}

namespace __StringNS__ {
    int __TextToInt__(char * text,int width){
        int sign;
        int i,val;

        sign = 1; val = 0;
        for (i=0; i<width; i++){
            switch (text[i]){
            case '+':
                break;
            case '-':
                sign = -1;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                val = val * 10 + text[i] - '0';
                break;
            case '.':
                return val * sign;;
            default:
                return val * sign;;
            }
        }
        return val * sign;
    }
    int __TextToHex__(char * text,int length){
        int sign;
        int i,val;

        sign = 1; val = 0;
        for (i=0; i<length; i++){
            switch (text[i]){
            case '+':
                break;
            case '-':
                sign = -1;
                break;
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                val = val * 16 + text[i] - '0';
                break;
            case 'A':
            case 'B':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
                val = val * 16 + text[i] - 'A' + 10;
                break;
            case 'a':
            case 'b':
            case 'c':
            case 'd':
            case 'e':
            case 'f':
                val = val * 16 + text[i] - 'a' + 10;
                break;
            default:
                return val * sign;;
            }
        }
        return val * sign;
    }
    double  __ten_factor__[32] = { 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001,
            0.00000000001, 0.000000000001, 0.0000000000001, 0.00000000000001, 0.000000000000001, 0.0000000000000001,
            0.00000000000000001, 0.000000000000000001, 0.0000000000000000001, 0.00000000000000000001,
            0.000000000000000000001, 0.0000000000000000000001, 0.00000000000000000000001, 0.000000000000000000000001,
            0.0000000000000000000000001, 0.00000000000000000000000001, 0.000000000000000000000000001,
            0.0000000000000000000000000001, 0.00000000000000000000000000001, 0.000000000000000000000000000001,
            0.0000000000000000000000000000001, 0.00000000000000000000000000000001 };
    float   __float_ten__ = 10.0;
    bool __TextToReal__(char * text,double * f, int width){
        double s; int i,sign,p; char t; bool startp;
        s = 0; p = 0; sign = 1; startp = false;
        i = 0;
        if (text[0]=='+') i++;
        else if (text[0]=='-'){ sign = -1; i++; }
        for (; i<width && p<32; i++){
            switch (text[i]){
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                t = text[i] - '0';
                startp==false? s = s * 10.0 + t : s += __ten_factor__[p] * t;
                if (startp == true) p++;
                break;
            case '.':
                startp = true;
                break;
            default:
                *f = 0;
                return false;
                break;
            }
        }
        *f = s * sign;
        return true;
    }
}
namespace __StringNS__ {
    double string_to_double(StringNS::string s){
        double ret = 0; __TextToReal__(s.text, &ret, s.length); return ret;
    }
    int string_to_int(StringNS::string s){
        return __TextToInt__(s.text, s.length);
    }
    int string_to_hex(StringNS::string s){
        return __TextToHex__(s.text, s.length);
    }
}

#endif
