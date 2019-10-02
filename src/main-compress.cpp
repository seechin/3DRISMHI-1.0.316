#include "compress.cpp"

size_t append_save_data(FILE ** pfile, char filename[MAX_PATH], FILE * flog, const char * title, const char * text, int nx, int ny, int nz, int nv, __REAL__ * data, double time_stamp, IET_Param * sys, void * _compressBuf, size_t _allocate_memory_size){
  // Append data directly to file. If file not exist then this will create one.
  // *pfile == stdin, will skip everything. This is to stop creating output files.
    if (!flog) return 0; if (!pfile || *pfile==stdin) return 0;
    if (!*pfile && !filename){ fprintf(flog, "%s : error : output file name not specified.\n", software_name); return 0; }
    if (!*pfile){
        if (StringNS::string(filename)=="stdout" || StringNS::string(filename)=="con" || StringNS::string(filename)=="screen"){
            *pfile = stdout;
        } else if (StringNS::string(filename)=="stderr"){
            *pfile = stderr;
        } else {
            if (sys->output_override==0 && access(filename, F_OK) != -1){ // file already exist and not appendable
                char filename_new[MAX_PATH]; filename_new[0] = 0;
                bool file_name_available = false;
                int max_rename_count = 10000;
                for (int i=1; i<max_rename_count; i++){
                    const char * filename_dir  = dirname(filename);
                    const char * filename_base = basename(filename);
                    if (!filename_dir || filename_dir[0]==0 || StringNS::string(filename_dir)=="."){
                        snprintf(filename_new, sizeof(filename_new), "#%s.%d#", filename_base, i);
                    } else {
                        snprintf(filename_new, sizeof(filename_new), "%s/#%s.%d#", filename_dir, filename_base, i);
                    }

                    if (access(filename_new, F_OK) == -1){
                        strncpy(filename, filename_new, MAX_PATH); file_name_available = true;
                        /*
                        rename(filename, filename_new);
                        if (access(filename, F_OK) != -1){
                            fprintf(flog, "%s : error : %s already exist and cannot be renamed. No output hence.\n", software_name, filename);
                            *pfile = stdin;
                        } else {
                            fprintf(flog, "%s : file %s is renamed to %s\n", software_name, filename, filename_new);
                            file_name_available = true;
                        }
                        */
                        break;
                    }
                }
                if (!file_name_available){
                    fprintf(flog, "%s : error : too many backup files: #%s.1# ~ #%s.%d#. No output hence.\n", software_name, filename, filename, max_rename_count);
                    *pfile = stdin;
                }
            } else if (sys->output_override==-1 && access(filename, F_OK) != -1){
                fprintf(flog, "%s : overide output file %s\n", software_name, filename);
            }
            if (!*pfile){
                *pfile = fopen(filename, sys->output_override==1?"a+":"w");
                if (!*pfile){
                    fprintf(flog, "%s : error : cannot write output to %s. No output hence.\n", software_name, filename);
                    *pfile = stdin;
                }
            }
        }
    }
    if (!*pfile || *pfile==stdin) return 0;

    size_t N3 = nx*ny*nz; size_t N4 = nv * N3; size_t data_size = sizeof(__REAL__)*nx*ny*nz*nv;

    /*if (sys->output_significant_digits<0 && -sys->output_significant_digits<sizeof(__REAL__)){
        if (-sys->output_significant_digits==sizeof(float)){
            float * dst = (float*)&data[0];
            for (size_t i4=0; i4<N4; i4++) dst[i4] = data[i4];
            data_size = sizeof(float)*nx*ny*nz*nv;
        } else if (-sys->output_significant_digits==sizeof(double)){
            double * dst = (double*)&data[0];
            for (size_t i4=0; i4<N4; i4++) dst[i4] = data[i4];
            data_size = sizeof(double)*nx*ny*nz*nv;
        }
    } else if (sys->output_significant_digits>0) {
        int significant_digits_bit = sys->output_significant_digits*log(10.0)/log(2.0);
        if (sys->output_significant_digits<=0 && sizeof(__REAL__)>=sizeof(float)){
            float * dst = (float*)&data[0]; unsigned int mask32 = generate_float_sig_dig_mask(significant_digits_bit);
            //printf("significant digits: %d, mask: %u\n", sys->output_significant_digits, mask32);
            for (size_t i4=0; i4<N4; i4++){ float tmp = data[i4]; bit_and(&dst[i4], &tmp, &mask32, 4); }
            data_size = sizeof(float)*nx*ny*nz*nv;
        } else if (sys->output_significant_digits<=18 && sizeof(__REAL__)>=sizeof(double)){
            double * dst = (double*)&data[0]; unsigned long mask64 = generate_double_sig_dig_mask(significant_digits_bit);
            //printf("significant digits: %d, mask: %lu\n", sys->output_significant_digits, mask64);
            for (size_t i4=0; i4<N4; i4++){ double tmp = data[i4]; bit_and(&dst[i4], &tmp, &mask64, 8); }
            data_size = sizeof(double)*nx*ny*nz*nv;
        }
    }*/
    if (sys->output_significant_digits>0) {
        int significant_digits_bit = sys->output_significant_digits*log(10.0)/log(2.0);
        if (sys->output_significant_digits<=18 && sizeof(__REAL__)>=sizeof(double)){
            double * dst = (double*)&data[0]; unsigned long mask64 = generate_double_sig_dig_mask(significant_digits_bit);
            //printf("significant digits: %d, mask: %lu\n", sys->output_significant_digits, mask64);
            for (size_t i4=0; i4<N4; i4++){ double tmp = data[i4]; bit_and(&dst[i4], &tmp, &mask64, 8); }
            data_size = sizeof(double)*nx*ny*nz*nv;
        }
    }

    unsigned short dimension[4]; dimension[0] = (unsigned short)nx; dimension[1] = (unsigned short)ny; dimension[2] = (unsigned short)nz; dimension[3] = (unsigned short)nv;

    bool show_save_details = sys->debug_level>0 || sys->detail_level>1;
    return write_data_page(*pfile, filename, time_stamp, title, show_save_details?flog:nullptr, dimension, data, data_size, text?text:"", sys->output_compress_level, sys->compress_page_size, _compressBuf, _allocate_memory_size);
}











bool load_array_from_file(__REAL__ * buffer, int nv, int nz, int ny, int nx, char * filename, const char * title, FILE * flog){
    printf("%s : load function is not finished\n", software_name);
    return false;
}

bool uncompress_data_array(const char * filename, FILE * file, CompressPageHeader & header, unsigned char * data, int len_data, bool show_error=true){
    printf("%s : load function is not finished\n", software_name);
    return false;
}
