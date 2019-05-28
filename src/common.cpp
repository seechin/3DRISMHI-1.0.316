//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------------------   Time and clock   -------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

double get_current_time_double(){
    struct timeval tv; gettimeofday(&tv, nullptr);
    return ((tv.tv_sec) * 1000 + (tv.tv_usec) / 1000);
}
double get_current_time_us(){
    struct timeval tv; gettimeofday(&tv, nullptr);
    return ((tv.tv_sec) * 1000000.0 + (tv.tv_usec));
}
char * get_current_time_text(char time_buffer[20]){
    time_t rawtime; time(&rawtime); struct tm * tmtime = localtime(&rawtime);
    snprintf(time_buffer, 20, "%04d-%02d-%02d,%02d:%02d:%02d", tmtime->tm_year+1900, tmtime->tm_mon+1, tmtime->tm_mday, tmtime->tm_hour, tmtime->tm_min, tmtime->tm_sec);
    return time_buffer;
}
StringNS::string file_extension(StringNS::string fn){
    int ie = 0;
    for (ie=fn.length-1; ie>=0; ie--) if (fn.text[ie]=='.') break;
    if (ie==0||fn.text[ie]!='.') return StringNS::string((char*)"");
    while (ie<fn.length && fn.text[ie]=='.') ie++;
    return fn.Substring(ie, fn.length-ie);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------   Trim real numbers : float and double   --------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

#define FILTER_REAL128_HIGH 0b11111111111111110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#define FILTER_REAL80_HIGH  0b11111111111111111000000000000000000000000000000000000000000000000000000000000000
#define FILTER_REAL64_HIGH  0b1111111111110000000000000000000000000000000000000000000000000000
#define FILTER_REAL32_HIGH  0b11111111100000000000000000000000
#define BITS_REAL128_HIGH   16
#define BITS_REAL80_HIGH    17
#define BITS_REAL64_HIGH    12
#define BITS_REAL32_HIGH    9
#define SIGN_REAL128_HIGH   0b10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#define SIGN_REAL80_HIGH    0b10000000000000000000000000000000000000000000000000000000000000000000000000000000
#define SIGN_REAL64_HIGH    0b1000000000000000000000000000000000000000000000000000000000000000
#define SIGN_REAL32_HIGH    0b10000000000000000000000000000000

unsigned int generate_float_sig_dig_mask(int bits){
    unsigned int mask = 0;
    for (int i=0; i<bits; i++){ mask >>= 1; mask |= SIGN_REAL32_HIGH; }
    for (int i=0; i<BITS_REAL32_HIGH; i++){ mask >>= 1; mask |= SIGN_REAL32_HIGH; }
    return mask;
}

unsigned long generate_double_sig_dig_mask(int bits){
    unsigned long mask = 0;
    for (int i=0; i<bits; i++){ mask >>= 1; mask |= SIGN_REAL64_HIGH; }
    for (int i=0; i<BITS_REAL64_HIGH; i++){ mask >>= 1; mask |= SIGN_REAL64_HIGH; }
    return mask;
}

inline void bit_and(void * _out, void * _in1, void * _in2, int size){
    for (int i=0; i<size; i++) ((unsigned char *)_out)[i] = ((unsigned char *)_in1)[i] & ((unsigned char *)_in2)[i];
    //unsigned char * out = (unsigned char *)_out;
    //unsigned char * in1 = (unsigned char *)_in1;
    //unsigned char * in2 = (unsigned char *)_in2;
    //for (int i=0; i<size; i++) out[i] = in1[i] & in2[i];
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-----------------------   ln, interpolation, and CRC   --------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

static double ____loge = -1;
double ln(double x){ if (____loge<=0) ____loge = log(2.7182818284590452353602874713527); return log(x) / ____loge; }

double interpolation123(double * a, int n, double dr, double x, double def_ret = 0){
    int i1 = (int)floor(x/dr); if (i1>x) i1--; int i2 = i1 + 1;
    if (i1<0 || i2>=n) return def_ret;
    return a[i1] * (i2-x) + a[i2] * (x-i1);
}
inline double vec_pow2(double x, double y, double z){ return x*x + y*y + z*z; }
inline double sqrt_norm2(double x, double y, double z){ return sqrt(x*x + y*y + z*z); }
unsigned int check_crc32(unsigned char * s, int l){
    unsigned long crc = crc32_zlib(0L, nullptr, 0); crc = crc32_zlib(crc, s, l); return (unsigned int) crc;
}
unsigned int check_real_crc(__REAL__ * s, int n, double errtr=1e-12){
    return check_crc32((unsigned char*)s, n*sizeof(__REAL__));
}
//unsigned int check_real_crc_with_truncation(__REAL__ * s, int n, double errtr=1e-12){
//    unsigned int v=0; double log10 = log(10.0);
//    v = crc32_zlib(0L, nullptr, 0);
//    for (int i=0; i<n; i++){
//        if (fabs(s[i])<1e-12){
//        } else {
//            double sv = floor(s[i]/errtr) * errtr;
//            double ls = log(fabs(sv)) / log10;
//            int ord = (int)floor(ls);
//            double val = exp((s[i]>0?log10:-log10)*(ls - ord));
//            float valf = floor(val/1e-5) * 1e-5;
//            v = crc32_zlib(v, (unsigned char *)&valf, sizeof(valf));
//            v = crc32_zlib(v, (unsigned char *)&ord, sizeof(ord));
//            //fprintf(stderr, "%18.6e %12f %5d\n", s[i], val, ord);
//        }
//    }
//    return v;
//}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------------   Special number printing   ---------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
char * print_memory_value(char * buffer, int size, size_t value){
    if (value<1000) snprintf(buffer, size, "%lu B", value);
    else if (value<10000) snprintf(buffer, size, "%.3f KB", value/1024.0);
    else if (value<100000) snprintf(buffer, size, "%.2f KB", value/1024.0);
    else if (value<1000000) snprintf(buffer, size, "%.1f KB", value/1024.0);
    else if (value<10000000) snprintf(buffer, size, "%.3f MB", value/1048576.0);
    else if (value<100000000) snprintf(buffer, size, "%.2f MB", value/1048576.0);
    else if (value<1000000000) snprintf(buffer, size, "%.1f MB", value/1048576.0);
    else if (value<10000000000) snprintf(buffer, size, "%.3f GB", value/1048576.0/1024.0);
    else if (value<100000000000) snprintf(buffer, size, "%.2f GB", value/1048576.0/1024.0);
    else if (value<1000000000000) snprintf(buffer, size, "%.1f GB", value/1048576.0/1024.0);
    else if (value<10000000000000) snprintf(buffer, size, "%.3f TB", value/1048576.0/1048576.0);
    else if (value<100000000000000) snprintf(buffer, size, "%.2f TB", value/1048576.0/1048576.0);
    else if (value<1000000000000000) snprintf(buffer, size, "%.1f TB", value/1048576.0/1048576.0);
    else snprintf(buffer, size, "%.3f PB", value/1048576.0/1048576.0/1024);
    return buffer;
}
char * print_percentage_value(char * buffer, int size, double value){
    if (value==1) snprintf(buffer, size, "100%%");
    else if (value>1) snprintf(buffer, size, "%.0f%%", value*100.0);
    else if (value>0.1) snprintf(buffer, size, "%.1f%%", value*100.0);
    else if (value>0.01) snprintf(buffer, size, "%.2f%%", value*100.0);
    else if (value>0.001) snprintf(buffer, size, "%.3f%%", value*100.0);
    else if (value>0.0001) snprintf(buffer, size, "%.3f%%%%", value*10000.0);
    else snprintf(buffer, size, "%g", value);
    return buffer;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//--------------------   Period Boundary Condition Treatment  ---------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
int ff_pbc_i(int x, int box_x){
    if (x<0) x += floor(abs(x/box_x)+1)*box_x;
    return x%box_x;
}
double ff_pbc_d(double x, double box_x){
    if (x<0) x += floor(fabs(x/box_x)+1)*box_x;
    return fmod(x, box_x);
}
#define ff_bound(a,b,c) (a<b? b : a>c? c : a)
