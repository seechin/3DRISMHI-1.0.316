//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//--------------------   String Functions: Enhanced StringNS   --------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
int analysis_line_params(StringNS::string sline, StringNS::string * sl, int maxnw, bool separate_szstr = false){
    int nw = StringNS::analysis_line(sline, sl, maxnw, separate_szstr);
    for (int i=0; i<nw; i++){
        if (sl[i].length>0){
            if (sl[i].text[0]=='#'){
                if (sl[i]=="#include" || sl[i]=="#warning" || sl[i]=="#error" || sl[i]=="#echo"){
                } else if (sl[i]=="#" && i+1<nw && (sl[i+1]=="include"||sl[i+1]=="warning"||sl[i+1]=="error"||sl[i+1]=="echo")){
                } else nw = i;
            } else if (sl[i].text[0]==';' || sl[i].text[0]=='@' || sl[i].text[0]=='%'){
                nw = i;
            } else if (sl[i].text[0]=='/' && sl[i].length>1 && sl[i].text[1]=='/'){
                nw = i;
            }
        } else nw = i;
    }
    return nw;
}
int parse_sz_line(char * input, char * ss[], int max_ss){
    StringNS::string sl[MAX_WORD];
    int nw = analysis_line_params(input, sl, MAX_WORD, true); if (nw>max_ss) nw = max_ss;
    for (int i=0; i<nw; i++) ss[i] = sl[i].text;
    return nw;
}
int translate_string_to_number(StringNS::string st, int ret_default){
    int ret = ret_default;

    if (StringNS::is_string_number(st)){
        ret = atoi(st.text);
    } else {
        int nw=0; StringNS::string sl[MAX_WORD]; char input[4096]; strncpy(input, st.text, sizeof(input));
        sl[0].text = input; nw = 1;
        for (int i=1; i<sizeof(input) && input[i]; i++){
            if (input[i]=='-' || input[i]==' ' || input[i]=='\t'){
                while (i<sizeof(input) && input[i] && (input[i]=='-' || input[i]==' ' || input[i]=='\t')){ input[i] = 0; i++; }
                if (input[i] && nw<MAX_WORD) sl[nw++].text = &input[i];
            }
        }
        for (int i=0; i<nw; i++) sl[i] = StringNS::string(sl[i].text);
        ret = 0; int value = 0; int last_hold_value = 0;
        for (int i=0; i<nw; i++){
            if (sl[i]=="and"){
            } else if (sl[i]=="zero"){
            } else if (sl[i]=="one"){       last_hold_value += 1;
            } else if (sl[i]=="two"){       last_hold_value += 2;
            } else if (sl[i]=="three"){     last_hold_value += 3;
            } else if (sl[i]=="four"){      last_hold_value += 4;
            } else if (sl[i]=="five"){      last_hold_value += 5;
            } else if (sl[i]=="six"){       last_hold_value += 6;
            } else if (sl[i]=="seven"){     last_hold_value += 7;
            } else if (sl[i]=="eight"){     last_hold_value += 8;
            } else if (sl[i]=="nine"){      last_hold_value += 9;
            } else if (sl[i]=="ten"){       last_hold_value += 10;
            } else if (sl[i]=="eleven"){    last_hold_value += 11;
            } else if (sl[i]=="twelve"){    last_hold_value += 12;
            } else if (sl[i]=="thirteen"){  last_hold_value += 13;
            } else if (sl[i]=="fourteen"){  last_hold_value += 14;
            } else if (sl[i]=="sixteen"){   last_hold_value += 15;
            } else if (sl[i]=="seventeen"){ last_hold_value += 16;
            } else if (sl[i]=="seventeen"){ last_hold_value += 17;
            } else if (sl[i]=="eighteen"){  last_hold_value += 18;
            } else if (sl[i]=="nineteen"){  last_hold_value += 19;
            } else if (sl[i]=="twenty"){    last_hold_value += 20;
            } else if (sl[i]=="thirty"){    last_hold_value += 30;
            } else if (sl[i]=="fourty"){    last_hold_value += 40;
            } else if (sl[i]=="forty"){     last_hold_value += 40;
            } else if (sl[i]=="fifty"){     last_hold_value += 50;
            } else if (sl[i]=="sixty"){     last_hold_value += 60;
            } else if (sl[i]=="seventy"){   last_hold_value += 70;
            } else if (sl[i]=="eighty"){    last_hold_value += 80;
            } else if (sl[i]=="ninety"){    last_hold_value += 90;
            } else if (sl[i]=="hundred"){   value = last_hold_value*100; last_hold_value = 0;
            } else if (sl[i]=="thousand"){  value = last_hold_value*1000; last_hold_value = 0;
            } else if (sl[i]=="million"){   value = last_hold_value*1000000; last_hold_value = 0;
            } else if (sl[i]=="billion"){   value = last_hold_value*1000000000; last_hold_value = 0;
            } else { ret = -1; }
//printf("ANALYSISING WORD %s : %d, %d, %d\n", sl[i].text, ret, value, last_hold_value);
        }
        if (ret >= 0) ret = last_hold_value + value;
    }
    return ret;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//--------------------------   Text File preprocessing   --------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
bool get_file_path(char * fn, char * pn){
    pn[0] = 0; int ifn = strlen(fn);
    for (; ifn>0 && fn[ifn]!='/'; ifn--);
    if (ifn>0 && fn[ifn]=='/'){
        char pt[MAX_PATH]; memset(pt, 0, MAX_PATH);
        for (int i=0; i<ifn; i++) pt[i] = fn[i]; if (!realpath(pt, pn)) return false;
    }
    return true;
}
bool cp_file_with_path(char * dst, const char * src, const char * path, bool append_short_path=true){
    bool ret = true;
    StringNS::string ssrc = src; if (ssrc=="stdin"||ssrc=="stdout"||ssrc=="stderr"||ssrc=="con"||ssrc=="screen"){
        strcpy(dst, src); memcpy(&dst[strlen(dst)+1], "\0\0\0\0", 2); return ret;
    }
    if (!path || !path[0] || src[0]=='/'){
        if (!realpath(src, dst)) ret = false;
//printf("[%s] : FROM %s\n", dst, src);
    } else {
        char pt[2*MAX_PATH]; strcpy(pt, path); strncat(pt, "/", sizeof(pt)-1-strlen(pt)); strncat(pt, src, sizeof(pt)-1-strlen(pt));
        if (!realpath(pt, dst)) ret = false;
//printf("[%s] : from %s\n", dst, pt);
    }
    if (!ret) snprintf(dst, MAX_PATH, " %s", src);
    if (append_short_path) strcpy(&dst[strlen(dst)+1], src); else memcpy(&dst[strlen(dst)+1], "\0\0\0\0", 2);
    return ret;
}
const char * get_second_fn(const char * fn){
    if (!fn[strlen(fn)+1]) return fn; else return &fn[strlen(fn)+1];
}
/*bool file_too_big(char * filename, int size_ceiling_MB){
    FILE * fp = fopen(filename, "r"); if (!fp) return true;
    fseek(fp, 0, SEEK_END); size_t size = ftell(fp); fclose(fp);
    if (1024*1024*(long int)size_ceiling_MB < (long int)size) return true;
    return false;
}*/
bool is_a_text_file(char * filename){
    StringNS::string sfilename = filename;
    if (sfilename=="stdin"||sfilename=="stdout"||sfilename=="stderr"||sfilename=="con"||sfilename=="screen") return true;

    FILE * fp = fopen(filename, "rb"); if (!fp) return false;
    bool ret = true; bool stop = false;
    unsigned char input[4096+1];
    bool line_comment = false; bool block_comment = false;
    while (true){
        memset(input, 0, sizeof(input));
        int br = fread(input, 1, sizeof(input)-1, fp);
        if (br<=0) break;
        for (int i=0; i<br; i++){
            if (!input[i]) {
                stop = true;
            } else if (input[i]=='\n' || input[i]=='\r'){
                line_comment = false;
            } else if (input[i]=='#' || input[i]=='@'){
                line_comment = true;
            } else if (input[i]=='/'){
                if (i+1<br && input[i+1]=='/'){
                    i++; line_comment = true;
                } else if (i+1<br && input[i+1]=='*'){
                    block_comment = true;
                }
            } else if (input[i]=='*'){
                if (i+1<br && input[i+1]=='/'){
                    i++; block_comment = false;
                }
            } else if (input[i]&0x80 || (input[i]<0x20 && input[i]!='\t')){
                if (line_comment){
                } else {
                    ret = false; stop = true;

                }
            }
        }
        if (br<sizeof(input)) break;
        if (stop) break;
    }
    if (fp) fclose(fp);

    return ret;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------   Memory Management for Multi Processing   -------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
#define MAX_MEMORYS 20000
int _memory_blk_total = 0; size_t _memory_total = 0; bool _ignore_memory_capacity = false;
void * _memory_pointers[MAX_MEMORYS]; long int _memory_size[MAX_MEMORYS];
size_t get_total_physical_memory(){
    size_t pages = sysconf(_SC_PHYS_PAGES);
    size_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
void mem_dispose_all(){
  #ifdef _LOCALPARALLEL_
    for (int i=0; i<_memory_blk_total; i++) if (_memory_pointers[i]){
//fprintf(stderr, "\033[32m%s %p\n\033[0m", _memory_mmap[i]?"munmap":"free", _memory_pointers[i]);
        munmap(_memory_pointers[i], _memory_size[i]);
    }
  #else
    for (int i=0; i<_memory_blk_total; i++) free(_memory_pointers[i]);
  #endif
}
void * memalloc(size_t size){
    if (!_ignore_memory_capacity){
        size_t physical_memory = get_total_physical_memory();
        if (size + _memory_total > physical_memory){ char buffer[64];
            fprintf(stderr, "%s : exceeds physical memory capacity. Mapped memory: %d bolcks, %s\n", software_name, _memory_blk_total, print_memory_value(buffer, sizeof(buffer), _memory_total));
            mem_dispose_all(); exit(-1);
        }
    }

    //{ char buffer[64]; size_t physical_memory = get_total_physical_memory(); fprintf(stderr, "%s : %.4f%% allocated, allocating %.4f%%. Mapped memory: %d bolcks, %s\n", software_name, (double)_memory_total/(double)physical_memory*100, (double)size/(double)physical_memory*100, _memory_blk_total, print_memory_value(buffer, sizeof(buffer), _memory_total)); }

//fprintf(stderr, "\33[1;31m::MEMALLOC: %12lu (%d %d)\n\33[0m", (long int)size, sizeof(long int), sizeof(size_t)); fflush(stderr);
  #ifdef _LOCALPARALLEL_
    void * p = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
    if (p==MAP_FAILED){ char buffer[64];
        fprintf(stderr, "%s : mmap failure. Mapped memory: %d bolcks, %s\n", software_name, _memory_blk_total, print_memory_value(buffer, sizeof(buffer), _memory_total));
        mem_dispose_all(); exit(-1);
    }
//fprintf(stderr, "\033[34mmmap: %p\n\033[0m", p);
//printf("\33[31m    allocated pointer %d\33[0m\n", (char *)p);
    if (p && p!=MAP_FAILED){
  #else
    void * p = malloc(size);
    if (!p){ char buffer[64];
        fprintf(stderr, "%s : malloc failure. Totally %d bolcks, %s was allocated and mapped.\n", software_name, _memory_blk_total, print_memory_value(buffer, sizeof(buffer), _memory_total));
        mem_dispose_all(); exit(-1);
    }
    if (p){
  #endif
        if (_memory_blk_total<MAX_MEMORYS){
            _memory_pointers[_memory_blk_total] = p;
            _memory_size[_memory_blk_total] = size;
        }
        _memory_blk_total ++; _memory_total += size;
    }

    return p;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------   Time cost for each part of the software   ------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
double __last_clk;
double __timer[20]; double __total_timer = 0;
void init_timer(){ __last_clk = get_current_time_double(); for (int i=0; i<20; i++) __timer[i] = 0; __total_timer = 0; }
void set_timer(double * timer){
    double __current_clk = get_current_time_double();
    if (timer) *timer += (__current_clk - __last_clk) / 1000.0;
    if (timer) __total_timer += (__current_clk - __last_clk) / 1000.0;
    __last_clk = __current_clk;
}
void lap_timer_none(){ set_timer(nullptr); }
void lap_timer_others(){ set_timer(&__timer[0]); }
void lap_timer_analysis_param(){ set_timer(&__timer[1]); }
void lap_timer_alloc_memory(){ set_timer(&__timer[2]); }
void lap_timer_io(){ set_timer(&__timer[3]); }
void lap_timer_uuv(){ set_timer(&__timer[4]); }
void lap_timer_livm(){ set_timer(&__timer[5]); }
void lap_timer_rism(){ set_timer(&__timer[6]); }
void lap_timer_diis(){ set_timer(&__timer[7]); }
void lap_timer_fftw(){ set_timer(&__timer[8]); }
char * display_time(double t, char * buf){
    int day = (int)(t / 86400);
    int hour = (int)((t - day*86400) / 3600);
    int minute = (int)((t - day*86400 - hour*3600) / 60);
    int second = (int)(((int)floor(t)) % 60);
    int milisec = (int)((t - floor(t)) * 1000);
    if (day<=0){
      if (hour<=0){
        if (minute<=0){
            sprintf(buf, "%d.%02d", second, (int)floor(milisec/10));
        } else sprintf(buf, "%d:%d.%02d", minute, second, (int)floor(milisec/10));
      } else sprintf(buf, "%d:%d:%d.%1d", hour, minute, second, (int)floor(milisec/100));
    } else sprintf(buf, "%d, %d:%d:%d.%1d", day, hour, minute, second, (int)floor(milisec/100));
    return buf;
}
double lap_display_timers(FILE * fp, bool display=true){
    char buf[1024]; char vbuf[64];

    if (display){
        fprintf(fp, "Computational time:\n");
        #define DISPLAY_TIMECOST_CRI 0.01
        if (__timer[1]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  read param        %12s s, %s\n", display_time(__timer[1], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[1]/__total_timer));
        if (__timer[2]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  memory management %12s s, %s\n", display_time(__timer[2], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[2]/__total_timer));
        if (__timer[3]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  input/output      %12s s, %s\n", display_time(__timer[3], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[3]/__total_timer));
        if (__timer[4]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  potential         %12s s, %s\n", display_time(__timer[4], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[4]/__total_timer));
        if (__timer[6]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  %-10s        %12s s, %s\n", "IET", display_time(__timer[6], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[6]/__total_timer));
        if (__timer[5]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  %-10s        %12s s, %s\n", "HI", display_time(__timer[5], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[5]/__total_timer));
        if (__timer[7]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  DIIS              %12s s, %s\n", display_time(__timer[7], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[7]/__total_timer));
        if (__timer[8]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  Convolution       %12s s, %s\n", display_time(__timer[8], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[8]/__total_timer));
        if (__timer[0]>DISPLAY_TIMECOST_CRI) fprintf(fp, "  others            %12s s, %s\n", display_time(__timer[0], buf), print_percentage_value(vbuf, sizeof(vbuf), __timer[0]/__total_timer));

        double this_time = __timer[1] + __timer[2] + __timer[3] + __timer[4];
        if (this_time>DISPLAY_TIMECOST_CRI*__total_timer) fprintf(fp, "  preparation       %12s s\n", display_time(this_time, buf));
        this_time = __timer[5] + __timer[6] + __timer[7] + __timer[8];
        if (this_time>DISPLAY_TIMECOST_CRI*__total_timer) fprintf(fp, "  calculation       %12s s\n", display_time(this_time, buf));

        fprintf(fp, "  total             %12s s\n", display_time(__total_timer, buf));
    }
    return __total_timer;
}
void display_memory_cost(FILE * flog, int memory_blk, long int memory_total){ char buffer[64];
    fprintf(flog, "Memory usage: totally %d blocks assigned, %s\n", memory_blk, print_memory_value(buffer, sizeof(buffer), memory_total));
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
double check_error_tol(double errtol){
    return errtol<MACHINE_REASONABLE_ERROR? MACHINE_REASONABLE_ERROR : errtol>1? 1 : errtol;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------------   Multi-Process Paralleling   -------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
#ifdef _LOCALPARALLEL_
    int * __active_threads = nullptr; pid_t * __forkpid;
    pthread_t * __thread_list; void init_mp(){
        __active_threads = (int*) memalloc(sizeof(int)); * __active_threads = 1;

        __thread_list = (pthread_t*) memalloc(sizeof(pthread_t)*MAX_THREADS);
        memset(__thread_list, 0, sizeof(pthread_t)*MAX_THREADS);
        #ifdef _LOCALPARALLEL_PTHREAD_
            __thread_list[0] = pthread_self();
        #else
            __thread_list[0] = 0;
        #endif

        __forkpid = (pid_t*) memalloc(sizeof(pid_t)*MAX_THREADS);
        memset(__forkpid, 0, sizeof(pid_t)*MAX_THREADS);
        __forkpid[0] = getpid();
    }

    int fork_here(int nthreads){ // return : thread indicator, 0: main thread
        if (nthreads==1){ * __active_threads = 1; return 0; }
        for (int i=1; i<nthreads; i++){
            int pid = fork();
            if (pid == 0){ (* __active_threads) ++; return i;
            } else {
                __forkpid[i] = pid;
            }
        }
        //printf("//// ////////////////// active %d\n", * __active_threads);
        return 0;
    }
    void join_here(int id, bool mp_by_fork){
      if (id>0){
        (*__active_threads) --; if (mp_by_fork) _exit(0);
      }
    }
    void wait_here(int id, int timeout_ms=1000){
        if (id>0) return ;
        double time0 = get_current_time_double();
        while ((*__active_threads) > 1){
            //printf("//// active %d\n", * __active_threads);
            //int status; wait(&status);
            usleep(1000);
            if (timeout_ms>0 && get_current_time_double()-time0 > timeout_ms){
                for (int id=1; id<MAX_THREADS; id++) if (__forkpid[id]>0) kill(__forkpid[id], SIGKILL);
                return;
            }
        }
    }
#endif

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//---------------   Vector, Matrix, Tensor 3D and 4D preparation   ----------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
int * init_int_vector(size_t n){ int * p = (int*) memalloc(sizeof(int) * n); memset(p, 0, sizeof(int) * n); return p; }
// a vector
__REAL__ * init_vector(size_t n){
    __REAL__ * p = (__REAL__*) memalloc(sizeof(__REAL__) * n); memset(p, 0, sizeof(__REAL__) * n); return p;
}
// a matrix
__REAL__ ** init_matrix(size_t m, size_t n, int overflow_chars=0){
    size_t lenh = sizeof(__REAL__*) * m; size_t len = lenh + sizeof(__REAL__) * m * n + overflow_chars;
    char * p = (char*) memalloc(len); memset(p, 0, len); __REAL__ * d = (__REAL__*)(p + lenh);
    __REAL__ ** a = (__REAL__**) p; for (size_t i=0; i<m; i++) a[i] = &d[i*n];
    return a;
}
void cp_matrix(__REAL__ ** src, __REAL__ ** dst, size_t m, size_t n){
    for (size_t i=0; i<m; i++) for (size_t j=0; j<n; j++) dst[i][j] = src[i][j];
}
// a 3D tensor pointer
__REAL__ *** init_tensor3d_pointer(size_t nz, size_t ny, size_t nx, size_t overflow_chars=0){
    size_t lenhz = sizeof(__REAL__**) * nz; size_t lenhy = sizeof(__REAL__*) * nz*ny ;
    size_t len = lenhz + lenhy + overflow_chars;
    //printf("ALLOCATING: %12d + %12d + %d\n", lenhz+lenhy, sizeof(__REAL__) * nx * ny * nz, overflow_chars);
    //printf("allocating %12f\n", sizeof(__REAL__) * nx * ny * nz/(double)len); fflush(stdout);
    char * p = (char*) memalloc(len); memset(p, 0, len);
    __REAL__ *** a = (__REAL__***) p; __REAL__ ** b = (__REAL__**) (p + lenhz);
    for (size_t i=0; i<nz; i++) a[i] = &b[i*ny];
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) a[i][j] = nullptr;
    return a;
}

__REAL__ *** init_tensor3d_pointer(__REAL__ *** src, size_t nz, size_t ny, size_t overflow_chars=0){
    size_t lenhz = sizeof(__REAL__**) * nz; size_t lenhy = sizeof(__REAL__*) * nz*ny ;
    size_t len = lenhz + lenhy + overflow_chars;
    //printf("ALLOCATING: %12d + %12d + %d\n", lenhz+lenhy, sizeof(__REAL__) * nx * ny * nz, overflow_chars);
    //printf("allocating %12f\n", sizeof(__REAL__) * nx * ny * nz/(double)len); fflush(stdout);
    char * p = (char*) memalloc(len); memset(p, 0, len);
    __REAL__ *** a = (__REAL__***) p; __REAL__ ** b = (__REAL__**) (p + lenhz);
    for (size_t i=0; i<nz; i++) a[i] = &b[i*ny];
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) a[i][j] = &src[i][j][0];
    return a;
}

__REAL__ *** init_tensor3d_pointer(__REAL__ * data, size_t nz, size_t ny, size_t nx, size_t overflow_chars=0){
    size_t lenhz = sizeof(__REAL__**) * nz; size_t lenhy = sizeof(__REAL__*) * nz*ny ;
    size_t len = lenhz + lenhy + overflow_chars;
    //printf("ALLOCATING: %12d + %12d + %d\n", lenhz+lenhy, sizeof(__REAL__) * nx * ny * nz, overflow_chars);
    //printf("allocating %12f\n", sizeof(__REAL__) * nx * ny * nz/(double)len); fflush(stdout);
    char * p = (char*) memalloc(len); memset(p, 0, len);
    __REAL__ *** a = (__REAL__***) p; __REAL__ ** b = (__REAL__**) (p + lenhz);
    for (size_t i=0; i<nz; i++) a[i] = &b[i*ny];
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) a[i][j] = &data[j*nx + i*nx*ny];
    return a;
}

// a 3D tensor
template <class DT> DT *** init_tensor3d(size_t nz, size_t ny, size_t nx, size_t overflow_chars=0){
    size_t lenhz = sizeof(DT**) * nz; size_t lenhy = sizeof(DT*) * nz*ny ;
    size_t len = lenhz + lenhy + sizeof(DT) * nx * ny * nz + overflow_chars;
//printf("ALLOCATING: %12d + %12d + %d\n", lenhz+lenhy, sizeof(DT) * nx * ny * nz, overflow_chars);
//printf("allocating %12f\n", sizeof(DT) * nx * ny * nz/(double)len); fflush(stdout);
    char * p = (char*) memalloc(len); memset(p, 0, len); DT * d = (DT*)(p + lenhz + lenhy);
    DT *** a = (DT***) p; DT ** b = (DT**) (p + lenhz);
    for (size_t i=0; i<nz; i++) a[i] = &b[i*ny];
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) a[i][j] = &d[j*nx + i*nx*ny];
    return a;
}
void cp_tensor3d(__REAL__ *** src, __REAL__ *** dst, size_t nz, size_t ny, size_t nx){
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) for (size_t k=0; k<nx; k++) dst[i][j][k] = src[i][j][k];
}
__REAL__ *** duplate_tensor3d(__REAL__ *** src, size_t nz, size_t ny, size_t nx, size_t overflow_chars=0){
    __REAL__ *** ret = init_tensor3d<__REAL__>(nz, ny, nx, overflow_chars);
    if (ret) cp_tensor3d(src, ret, nz, ny, nx);
    return ret;
}
void clear_tensor3d(__REAL__ *** dst, size_t nz, size_t ny, size_t nx){
    for (size_t i=0; i<nz; i++) for (size_t j=0; j<ny; j++) for (size_t k=0; k<nx; k++) dst[i][j][k] = 0;
}
void cp_tensor3d(__REAL__ *** src, __REAL__ *** dst, size_t n3){
    __REAL__ * d = &dst[0][0][0]; __REAL__ * s = & src[0][0][0];
    for (size_t i3=0; i3<n3; i3++) d[i3] = s[i3];
}
void clear_tensor3d(__REAL__ *** dst, size_t n3){
    __REAL__ * d = &dst[0][0][0];
    for (size_t i3=0; i3<n3; i3++) d[i3] = 0;
}
void clear_tensor3d_double(double *** dst, size_t n3){
    double * d = &dst[0][0][0];
    for (size_t i3=0; i3<n3; i3++) d[i3] = 0;
}
__REAL__ **** init_tensor4d(size_t n4, size_t n3, size_t n2, size_t n1, size_t overflow_chars=0){
    size_t lenh4 = sizeof(__REAL__***) * n4; size_t lenh3 = sizeof(__REAL__**) * n4*n3; size_t lenh2 = sizeof(__REAL__*) * n4*n3*n2;
    size_t len1 = sizeof(__REAL__) * n4*n3*n2*n1;
    size_t len = lenh4 + lenh3 + lenh2 + len1 + overflow_chars;
    char * p = (char*) memalloc(len);
    __REAL__ **** a4 = (__REAL__****) p;
    __REAL__ *** a3 = (__REAL__***) (p + lenh4);
    __REAL__ ** a2 = (__REAL__**) (p + lenh4 + lenh3);
    __REAL__ * a1 = (__REAL__*) (p + lenh4 + lenh3 + lenh2); __REAL__ * data = a1;
//printf("ALLOCATING: %12d + %12d + %d\n", lenh4+lenh3+lenh2, len1, overflow_chars);
//printf("allocating %12f\n", len1/(double)len); fflush(stdout);
//printf("pointer layer 1, from %8ld to %8ld -> %8ld to %8ld\n", (char*)&a4[0]-p, (char*)&a4[n4]-p, (char*)&a3[0]-p, (char*)&a3[n4*n3]-p); fflush(stdout);
    for (size_t i4=0; i4<n4; i4++) a4[i4] = & a3[i4*n3];
//printf("pointer layer 2, from %8ld to %8ld -> %8ld to %8ld\n", (char*)&a4[0][0]-p, (char*)&a4[n4-1][n3]-p, (char*)&a2[0]-p, (char*)&a2[n4*n3*n2]-p); fflush(stdout);
    for (size_t i4=0; i4<n4; i4++) for (size_t i3=0; i3<n3; i3++) a4[i4][i3] = & a2[(i3+i4*n3)*n2];
//printf("pointer layer 3, from %8ld to %8ld -> %8ld to %8ld\n", (char*)&a4[0][0][0]-p, (char*)&a4[n4-1][n3-1][n2]-p, (char*)&a1[0]-p, (char*)&a1[n4*n3*n2*n1]-p); fflush(stdout);
    for (size_t i4=0; i4<n4; i4++) for (size_t i3=0; i3<n3; i3++) for (size_t i2=0; i2<n2; i2++) a4[i4][i3][i2] = & a1[(i2+(i3+i4*n3)*n2)*n1];
//printf("pointer layer 4, from %8ld to %8ld -> %8ld to %8ld\n", (char*)&a4[0][0][0][0]-p, (char*)&a4[n4-1][n3-1][n2-1][n1]-p, (char*)&a1[0]-p, len); fflush(stdout);
    for (size_t i4=0; i4<n4; i4++) for (size_t i3=0; i3<n3; i3++) for (size_t i2=0; i2<n2; i2++) for (size_t i1=0; i1<n1; i1++) a4[i4][i3][i2][i1] = 0;
//printf("init tensor4d done\n"); fflush(stdout);
    return a4;
}
void cp_tensor4d(__REAL__ **** src, __REAL__ **** dst, size_t n4){
    __REAL__ * d = &dst[0][0][0][0]; __REAL__ * s = & src[0][0][0][0];
    for (size_t i4=0; i4<n4; i4++) d[i4] = s[i4];
}
void clear_tensor4d(__REAL__ **** dst, size_t n4){
    __REAL__ * d = &dst[0][0][0][0];
    for (size_t i4=0; i4<n4; i4++) d[i4] = 0;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//---------------------------   Serial 3D Convolution   ---------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
inline double interpolate(double k, __REAL__ * fk, double xvv_k_shift, int nr, double dk){
    int i1 = (int)floor((k-xvv_k_shift)/dk); int i2 = i1 + 1;
    if (i2>=nr) return 0; if (i1<0) return fk[0];
    return (i2 - k/dk) * fk[i1] + (k/dk - i1) * fk[i2];
}

void perform_3rx1k_convolution(__REAL__ *** f3r[], int nx, int ny, int nz, Vector box, int ni, int nj, __REAL__ *** f1k, double dk, double xvv_k_shift, int nf1k, __REAL__ *** out[], double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb, bool clear_out=true){
    double dx = box.x/nx; double dy = box.y/ny; double dz = box.z/nz; // int flips[3] = { nx, ny, nz };
    double dkx = 2*PI/(nx * dx); double dky = 2*PI/(ny * dy); double dkz = 2*PI/(nz * dz);
    double convolution_factor = 1.0 / (nx * ny * nz);
    size_t N3 = nx * ny * nz;

    double * fftin1 = &fftin[0][0][0];
    double * fftout1 = &fftout[0][0][0];

    if (clear_out) clear_tensor4d(out, ni * N3);
    for (int si=0; si<ni; si++){
        __REAL__ * out1 = &out[si][0][0][0];
        for (int sj=0; sj<nj; sj++) if (f1k[si][sj]){
            __REAL__ * f3r1 = &f3r[sj][0][0][0];
          // f3k_j = FFT[f3r_j] -> fftout
            for (size_t i3=0; i3<N3; i3++) fftin1[i3] = f3r1[i3];
            fftw_execute(planf);
          // f1k_ij f3k_j -> fftout
            for (int iz=0; iz<=nz/2; iz++) for (int iy=0; iy<=ny/2; iy++) for (int ix=0; ix<=nx/2; ix++){
                double kx = dkx * (ix+0); double ky = dky * (iy+0); double kz = dkz * (iz+0);
                double k = sqrt(kx*kx + ky*ky + kz*kz);
                double intp = interpolate(k, f1k[si][sj], xvv_k_shift*dk, nf1k, dk);
                double intp_factor = intp;
                unsigned int mask = 0; if (ix==0||ix>=nx-ix) mask |= 1; if (iy==0||iy>=ny-iy) mask |= 2; if (iz==0||iz>=nz-iz) mask |= 4;
                fftout[iz][iy][ix] *= intp_factor;
                if (!(mask&1)) fftout[iz][iy][nx-ix] *= intp_factor;
                if (!(mask&2)) fftout[iz][ny-iy][ix] *= intp_factor;
                if (!(mask&3)) fftout[iz][ny-iy][nx-ix] *= intp_factor;
                if (!(mask&4)) fftout[nz-iz][iy][ix] *= intp_factor;
                if (!(mask&5)) fftout[nz-iz][iy][nx-ix] *= intp_factor;
                if (!(mask&6)) fftout[nz-iz][ny-iy][ix] *= intp_factor;
                if (!(mask&7)) fftout[nz-iz][ny-iy][nx-ix] *= intp_factor;
            }
          // FFTi[f1k_ij f3k_j] -> fftin
            fftw_execute(planb);
          // f1r_ij * f3k_j -(+)-> out_i
            for (size_t i3=0; i3<N3; i3++) out1[i3] += fftin1[i3] * convolution_factor;
        }
    }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------------   MP Parallel 3D Convolution   ------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class RISMHI3D_FFTW_MP_UNIT {
  public:
    __REAL__ *** fftin;
    __REAL__ *** fftout;
    fftw_plan planf, planb;
    int si, sj;
    bool task_done; __REAL__ * out1;
    void init(int nx, int ny, int nz){
        si = sj = 0;
        fftin = init_tensor3d<__REAL__>(nz, ny, nx, 0);
        fftout = init_tensor3d<__REAL__>(nz, ny, nx, 0);
        planf = fftw_plan_r2r_3d(nz, ny, nx, &fftin[0][0][0], &fftout[0][0][0], (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, (fftw_r2r_kind)FFTW_FORWARD, FFTW_ESTIMATE);
        planb = fftw_plan_r2r_3d(nz, ny, nx, &fftout[0][0][0], &fftin[0][0][0], (fftw_r2r_kind)FFTW_BACKWARD, (fftw_r2r_kind)FFTW_BACKWARD, (fftw_r2r_kind)FFTW_BACKWARD, FFTW_ESTIMATE);
    }
};
class RISMHI3D_FFTW_MP {
  public:
    RISMHI3D_FFTW_MP_UNIT fft[MAX_THREADS];
    int nx, ny, nz, np;
    double dx, dy, dz;
    double dkx, dky, dkz, dk, xvv_k_shift, convolution_factor;
    size_t N3; int nf1k;
    int * mp_tasks; int n_active_jobs;
  public:
    __REAL__ **** f3r; __REAL__ *** f1k; __REAL__ **** out;
  public:
    void set_scale(Vector box){
        dx = box.x/nx; dy = box.y/ny; dz = box.z/nz;
        dkx = 2*PI/(nx * dx); dky = 2*PI/(ny * dy); dkz = 2*PI/(nz * dz);
        convolution_factor = 1.0 / (nx * ny * nz);
    }
    void init(int _np, int * _mp_tasks, int _nx, int _ny, int _nz, double _xvv_k_shift){
        np = _np; mp_tasks = _mp_tasks; n_active_jobs = 0; xvv_k_shift = _xvv_k_shift;
        nx = _nx; ny = _ny; nz = _nz; N3 = nx * ny * nz;
        for (int i=0; i<np; i++) fft[i].init(nx, ny, nz);
    }
    void set_field(__REAL__ **** _f3r, __REAL__ *** _f1k, double _dk, int _nf1k, __REAL__ **** _out){
        f3r = _f3r; f1k = _f1k; out = _out; dk = _dk; nf1k = _nf1k;
    }
    void reset_for_calculation(){
        for (int i=0; i<np; i++){
            clear_tensor3d_double(fft[i].fftin, nx*ny*nz);
            clear_tensor3d_double(fft[i].fftout, nx*ny*nz);
        }
    }
};

void perform_3rx1k_convolution_r1(int id, RISMHI3D_FFTW_MP * param, int si, int sj){
    __REAL__ **** f3r = param->f3r;
    __REAL__ ***  f1k = param->f1k;
    //__REAL__ **** out = param->out;
    int nx = param->nx; int ny = param->ny; int nz = param->nz;
    double dkx = param->dkx; double dky = param->dky; double dkz = param->dkz; double dk = param->dk;
    size_t N3 = nx * ny * nz;
    double * fftin1 = &param->fft[id].fftin[0][0][0];
    double *** fftout = param->fft[id].fftout;

    //__REAL__ * out1 = &out[si][0][0][0];

    param->fft[id].out1 = &param->out[si][0][0][0];
    param->fft[id].task_done = false;

    if (!f1k[si][sj]){
        memset(fftin1, 0, N3*sizeof(__REAL__));
    } else {
        __REAL__ * f3r1 = &f3r[sj][0][0][0];
      // f3k_j = FFT[f3r_j] -> fftout
        for (size_t i3=0; i3<N3; i3++) fftin1[i3] = f3r1[i3];
        fftw_execute(param->fft[id].planf);
      // f1k_ij f3k_j -> fftout
        for (int iz=0; iz<=nz/2; iz++) for (int iy=0; iy<=ny/2; iy++) for (int ix=0; ix<=nx/2; ix++){
            double kx = dkx * (ix+0); double ky = dky * (iy+0); double kz = dkz * (iz+0);
            double k = sqrt(kx*kx + ky*ky + kz*kz);
            double intp = interpolate(k, f1k[si][sj], param->xvv_k_shift*dk, param->nf1k, dk);
            double intp_factor = intp * param->convolution_factor;
            unsigned int mask = 0; if (ix==0||ix>=nx-ix) mask |= 1; if (iy==0||iy>=ny-iy) mask |= 2; if (iz==0||iz>=nz-iz) mask |= 4;
            fftout[iz][iy][ix] *= intp_factor;
            if (!(mask&1)) fftout[iz][iy][nx-ix] *= intp_factor;
            if (!(mask&2)) fftout[iz][ny-iy][ix] *= intp_factor;
            if (!(mask&3)) fftout[iz][ny-iy][nx-ix] *= intp_factor;
            if (!(mask&4)) fftout[nz-iz][iy][ix] *= intp_factor;
            if (!(mask&5)) fftout[nz-iz][iy][nx-ix] *= intp_factor;
            if (!(mask&6)) fftout[nz-iz][ny-iy][ix] *= intp_factor;
            if (!(mask&7)) fftout[nz-iz][ny-iy][nx-ix] *= intp_factor;
        }
      // FFTi[f1k_ij f3k_j] -> fftin
        fftw_execute(param->fft[id].planb);
    }

  // f1r_ij * f3k_j -(+)-> out_i
    //for (size_t i3=0; i3<N3; i3++) out1[i3] += fftin1[i3];
    param->fft[id].task_done = true;
}

#ifdef _LOCALPARALLEL_
    /*
    void sumback_3rx1k_convolution_r1(RISMHI3D_FFTW_MP * fftw_mp, int id){
        int it = id;
        if (fftw_mp->fft[it].task_done){
            size_t N3 = fftw_mp->nx * fftw_mp->ny * fftw_mp->nz;
            __REAL__ * out1 = fftw_mp->fft[it].out1;
            double * in1 = &fftw_mp->fft[it].fftin[0][0][0];
            for (size_t i3=0; i3<N3; i3++) out1[i3] += in1[i3];
            fftw_mp->fft[it].task_done = false;
        }
    }
    void sumback_3rx1k_convolution_r1(RISMHI3D_FFTW_MP * fftw_mp){
        for (int it=0; it<fftw_mp->np; it++) sumback_3rx1k_convolution_r1(fftw_mp, it);
    }
     */

    void sumback_3rx1k_convolution_t1(RISMHI3D_FFTW_MP * fftw_mp, int id, int n_active_jobs){
        size_t N3 = fftw_mp->nx * fftw_mp->ny * fftw_mp->nz;
        size_t begin = N3 / fftw_mp->np * id; size_t end = N3 / fftw_mp->np * (id+1); if (id+1 >= fftw_mp->np) end = N3;
        for (int it=0; it<n_active_jobs; it++){
            __REAL__ * out1 = fftw_mp->fft[it].out1;
            double * in1 = &fftw_mp->fft[it].fftin[0][0][0];
            for (size_t i3=begin; i3<end; i3++) out1[i3] += in1[i3];
        }
    }

    void sumback_3rx1k_convolution(RISMHI3D_FFTW_MP * fftw_mp){
        if (fftw_mp->np>1){
            for (int it=1; it<fftw_mp->np; it++) fftw_mp->mp_tasks[it] = MPTASK_MERGE_FFT_DATA;
            sumback_3rx1k_convolution_t1(fftw_mp, 0, fftw_mp->n_active_jobs);
            while (true){ int n_busy = 0;
                for (int i=0; i<fftw_mp->np; i++) if (fftw_mp->mp_tasks[i]!=MPTASK_NONE) n_busy++;
                if (n_busy<=0) break; usleep(100);
            }
        } else {
            for (int it=0; it<fftw_mp->np; it++) sumback_3rx1k_convolution_t1(fftw_mp, it, fftw_mp->n_active_jobs);
        }
    }
#endif

void perform_3rx1k_convolution(RISMHI3D_FFTW_MP * fftw_mp, __REAL__ *** f3r[], int nx, int ny, int nz, Vector box, int ni, int nj, __REAL__ *** f1k, double dk, double xvv_k_shift, int nf1k, __REAL__ *** out[], double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb, bool clear_out){

  #if defined(_LOCALPARALLEL_) && defined(_LOCALPARALLEL_FFTW_)
    if (fftw_mp && fftw_mp->np>1){
        fftw_mp->set_scale(box);
        fftw_mp->set_field(f3r, f1k, dk, nf1k, out);
        if (clear_out) clear_tensor4d(out, ni * fftw_mp->nx*fftw_mp->ny*fftw_mp->nz);

        int ss = 0; while (ss < ni*nj){
          // 1. assign jobs
            for (int i=0; i<fftw_mp->np; i++) fftw_mp->mp_tasks[i] = MPTASK_NONE;
            int njobs = 0; for (int i=0; i<fftw_mp->np && ss<ni*nj; i++){
                int si = ss/nj; int sj = ss%nj;
                if (f1k[si][sj]){
                    fftw_mp->fft[i].si = si; fftw_mp->fft[i].sj = sj;
                    fftw_mp->mp_tasks[i] = MPTASK_FFTW;
                    ss ++; njobs = i+1;
                } else {
                    ss ++; i--;
                }
            }
            fftw_mp->n_active_jobs = njobs;
          // 2. work on task 0 and wait
            perform_3rx1k_convolution_r1(0, fftw_mp, fftw_mp->fft[0].si, fftw_mp->fft[0].sj);
            fftw_mp->mp_tasks[0] = MPTASK_NONE;
            while (true){ int n_busy = 0;
                for (int i=0; i<njobs; i++) if (fftw_mp->mp_tasks[i]!=MPTASK_NONE) n_busy++;
                if (n_busy<=0) break; usleep(100);
            }
          // 3. collect data
            //for (int i=0; i<njobs; i++) sumback_3rx1k_convolution_r1(fftw_mp, i);
            //for (int it=0; it<fftw_mp->np; it++) sumback_3rx1k_convolution_t1(fftw_mp, it, njobs);
            sumback_3rx1k_convolution(fftw_mp);
          // 4. others
            fftw_mp->n_active_jobs = 0;
        }

        /*
        for (int ss=0; ss<ni*nj; ss++){
            int si = ss/nj; int sj = ss%nj;
            if (fftw_mp->assign_task_to_node_0 && fftw_mp->np>1){
                int i_task_idle = 0;
                for (int i=1; i<fftw_mp->np; i++) if (fftw_mp->mp_tasks[i]==MPTASK_NONE){ i_task_idle = i; break; }
                if (i_task_idle>0){
                    sumback_3rx1k_convolution_r1(fftw_mp, i_task_idle);
                    fftw_mp->fft[i_task_idle].si = si; fftw_mp->fft[i_task_idle].sj = sj;
                    fftw_mp->mp_tasks[i_task_idle] = MPTASK_FFTW;
                } else {
                    perform_3rx1k_convolution_r1(0, fftw_mp, si, sj);
                    sumback_3rx1k_convolution_r1(fftw_mp);
                }
            } else {
                while (true){
                    int i_task_idle = 0;
                    for (int i=1; i<fftw_mp->np; i++) if (fftw_mp->mp_tasks[i]==MPTASK_NONE){ i_task_idle = i; break; }
                    if (i_task_idle>0){
                        sumback_3rx1k_convolution_r1(fftw_mp, i_task_idle);
                        fftw_mp->fft[i_task_idle].si = si; fftw_mp->fft[i_task_idle].sj = sj;
                        fftw_mp->mp_tasks[i_task_idle] = MPTASK_FFTW;
                        break;
                    } else {
                        usleep(1000); continue;
                    }
                }
            }
        }

        double time0 = get_current_time_double();
        double timeup_ms = -1;
        while (true){
            bool finished = true; for (int i=1; i<fftw_mp->np; i++) if (fftw_mp->mp_tasks[i]!=MPTASK_NONE){ finished = false; break; }
            if (finished) break;
            if (timeup_ms>0 && get_current_time_double()-time0 > timeup_ms) break;
            usleep(100);
        }
        sumback_3rx1k_convolution_r1(fftw_mp);
         */
    } else {
        perform_3rx1k_convolution(f3r, nx, ny, nz, box, ni, nj, f1k, dk, xvv_k_shift, nf1k, out, fftin, fftout, planf, planb, clear_out);
    }
  #else
    /*for (int si=0; si<ni; si++) for (int sj=0; sj<nj; sj++){
        perform_3rx1k_convolution_r1(0, fftw_mp, si, sj);
        sumback_3rx1k_convolution_r1(fftw_mp);
    }*/
    perform_3rx1k_convolution(f3r, nx, ny, nz, box, ni, nj, f1k, dk, xvv_k_shift, nf1k, out, fftin, fftout, planf, planb, clear_out);
  #endif

}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//--------------   MP Parallel of Short Range Potential Calculation   -------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class RISMHI3D_FFSR_MP {
  public:
    int irange[MAX_THREADS][3];
    __REAL__ **** lj[MAX_THREADS];
    __REAL__ **** others[MAX_THREADS];
    __REAL__ *** coulsr[MAX_THREADS];
    __REAL__ *** coulp2[MAX_THREADS][3];
    __REAL__ *** r2uvmin[MAX_THREADS];
    __REAL__ *** pseudoliquid_potential[MAX_THREADS];
    //RISMHI3D_FFTW_MP_UNIT fft[MAX_THREADS];
    int np, nx, ny, nz, nv, N3, N4;
    int * mp_tasks;
  public:
    bool param_b[10]; int param_i[10]; double param_d[10]; __REAL__ *** param_t3[10]; __REAL__ **** param_t4[10];
  public:
    void init(int _np, int * _mp_tasks, int _nx, int _ny, int _nz, int _nv){
        np = _np; mp_tasks = _mp_tasks;
        nx = _nx; ny = _ny; nz = _nz; nv = _nv;
        N3 = nx * ny * nz; N4 = nv * N3;
        for (int it=0; it<np; it++){
            lj[it] = init_tensor4d(nv, nz, ny, nx, 0);
            others[it] = init_tensor4d(7, nz, ny, nx, 0);
            coulsr[it] = others[it][0];
            r2uvmin[it] = others[it][1];
            coulp2[it][0] = others[it][3];
            coulp2[it][1] = others[it][4];
            coulp2[it][2] = others[it][5];
            pseudoliquid_potential[it] = others[it][6];
        }
    }
    void reset_for_calculation(bool _clj, bool _ccoul){
        for (int it=0; it<np; it++){
            irange[it][0] = irange[it][1] = irange[it][2] = 0;
            if (_clj) clear_tensor4d(lj[it], N4);
            if (_ccoul) clear_tensor4d(others[it], N3*7);
            if (_ccoul) for (size_t i3=0; i3<N3; i3++) r2uvmin[it][0][0][i3] = -1;
        }
        for (int i=0; i<10; i++){ param_b[i] = false; param_i[i] = 0; param_d[i] = 0; param_t3[i] = nullptr; param_t4[i] = nullptr; }
    }
};


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------------------   Other procedures   -----------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------












/*

typedef double convolution_kernel_ptr(double k, int si, int sj, void * params);

void perform_3rxf1k_convolution(__REAL__ *** f3r[], int nx, int ny, int nz, Vector box, int ni, int nj, convolution_kernel_ptr * kernel_function, void * kernel_param, __REAL__ *** out[], double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb, bool clear_out=true){
    double dx = box.x/nx; double dy = box.y/ny; double dz = box.z/nz; // int flips[3] = { nx, ny, nz };
    double dkx = 2*PI/(nx * dx); double dky = 2*PI/(ny * dy); double dkz = 2*PI/(nz * dz);
    double convolution_factor = 1.0 / (nx * ny * nz);
    size_t N3 = nx * ny * nz;

    double * fftin1 = &fftin[0][0][0];
    double * fftout1 = &fftout[0][0][0];

    if (clear_out) clear_tensor4d(out, ni * N3);
    for (int si=0; si<ni; si++){
        __REAL__ * out1 = &out[si][0][0][0];
        for (int sj=0; sj<nj; sj++){
            __REAL__ * f3r1 = &f3r[sj][0][0][0];
          // f3k_j = FFT[f3r_j] -> fftout
            for (size_t i3=0; i3<N3; i3++) fftin1[i3] = f3r1[i3];
            fftw_execute(planf);
          // f1k_ij f3k_j -> fftout
            for (int iz=0; iz<=nz/2; iz++) for (int iy=0; iy<=ny/2; iy++) for (int ix=0; ix<=nx/2; ix++){
                double kx = dkx * (ix+0); double ky = dky * (iy+0); double kz = dkz * (iz+0);
                double k = sqrt(kx*kx + ky*ky + kz*kz);
                double intp = kernel_function(k, si, sj, kernel_param);
                double intp_factor = intp;
                unsigned int mask = 0; if (ix==0||ix>=nx-ix) mask |= 1; if (iy==0||iy>=ny-iy) mask |= 2; if (iz==0||iz>=nz-iz) mask |= 4;
                fftout[iz][iy][ix] *= intp_factor;
                if (!(mask&1)) fftout[iz][iy][nx-ix] *= intp_factor;
                if (!(mask&2)) fftout[iz][ny-iy][ix] *= intp_factor;
                if (!(mask&3)) fftout[iz][ny-iy][nx-ix] *= intp_factor;
                if (!(mask&4)) fftout[nz-iz][iy][ix] *= intp_factor;
                if (!(mask&5)) fftout[nz-iz][iy][nx-ix] *= intp_factor;
                if (!(mask&6)) fftout[nz-iz][ny-iy][ix] *= intp_factor;
                if (!(mask&7)) fftout[nz-iz][ny-iy][nx-ix] *= intp_factor;
            }
          // FFTi[f1k_ij f3k_j] -> fftin
            fftw_execute(planb);
          // f1r_ij * f3k_j -(+)-> out_i
            for (size_t i3=0; i3<N3; i3++) out1[i3] += fftin1[i3] * convolution_factor;
        }
    }
}
*/

void perform_3rxf1k_convolution(__REAL__ *** f3r[], int nx, int ny, int nz, Vector box, int ni, int nj, __REAL__ *** f1k, double dk, double xvv_k_shift, int nf1k, __REAL__ *** out[], double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb, bool clear_out=true){
    double dx = box.x/nx; double dy = box.y/ny; double dz = box.z/nz; // int flips[3] = { nx, ny, nz };
    double dkx = 2*PI/(nx * dx); double dky = 2*PI/(ny * dy); double dkz = 2*PI/(nz * dz);
    double convolution_factor = 1.0 / (nx * ny * nz);
    size_t N3 = nx * ny * nz;

    double * fftin1 = &fftin[0][0][0];
    double * fftout1 = &fftout[0][0][0];

    if (clear_out) clear_tensor4d(out, ni * N3);
    for (int si=0; si<ni; si++){
        __REAL__ * out1 = &out[si][0][0][0];
        for (int sj=0; sj<nj; sj++){
            __REAL__ * f3r1 = &f3r[sj][0][0][0];
          // f3k_j = FFT[f3r_j] -> fftout
            for (size_t i3=0; i3<N3; i3++) fftin1[i3] = f3r1[i3];
            fftw_execute(planf);
          // f1k_ij f3k_j -> fftout
            for (int iz=0; iz<=nz/2; iz++) for (int iy=0; iy<=ny/2; iy++) for (int ix=0; ix<=nx/2; ix++){
                double kx = dkx * (ix+0); double ky = dky * (iy+0); double kz = dkz * (iz+0);
                double k = sqrt(kx*kx + ky*ky + kz*kz);
                //double intp = interpolate(k, f1k[si][sj], xvv_k_shift*dk, nf1k, dk);
                double density = 33.4; double ksigma = k * pow(fabs(1.0/density/4.0*3.0/PI), 1.0/3);
                double intp = k==0? 1 : (4*PI*(sin(ksigma) - ksigma*cos(ksigma))/k/k/k)*density;
                double intp_factor = intp;
                unsigned int mask = 0; if (ix==0||ix>=nx-ix) mask |= 1; if (iy==0||iy>=ny-iy) mask |= 2; if (iz==0||iz>=nz-iz) mask |= 4;
                fftout[iz][iy][ix] *= intp_factor;
                if (!(mask&1)) fftout[iz][iy][nx-ix] *= intp_factor;
                if (!(mask&2)) fftout[iz][ny-iy][ix] *= intp_factor;
                if (!(mask&3)) fftout[iz][ny-iy][nx-ix] *= intp_factor;
                if (!(mask&4)) fftout[nz-iz][iy][ix] *= intp_factor;
                if (!(mask&5)) fftout[nz-iz][iy][nx-ix] *= intp_factor;
                if (!(mask&6)) fftout[nz-iz][ny-iy][ix] *= intp_factor;
                if (!(mask&7)) fftout[nz-iz][ny-iy][nx-ix] *= intp_factor;
            }
          // FFTi[f1k_ij f3k_j] -> fftin
            fftw_execute(planb);
          // f1r_ij * f3k_j -(+)-> out_i
            for (size_t i3=0; i3<N3; i3++) out1[i3] += fftin1[i3] * convolution_factor;
        }
    }
}
