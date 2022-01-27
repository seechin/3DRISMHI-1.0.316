const char * software_name = "ts4sdump";
const char * software_version = "0.296.1749";
const char * copyright_string = "(c) Cao Siqin";

#define     __REAL__    double
#define     MACHINE_REASONABLE_ERROR    1e-12
#define _TTYPROMPTCOLOR_

#include    "header.h"
#include    "main-header.h"
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
#include    <dirent.h>
#include    <sys/time.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>
#ifdef _LIBZ_
    #include  <zlib.h>
#endif

#define PI  3.1415926535897932384626433832795
#define EE  2.7182818284590452353602874713527
#define COULCOOEF 138.9354846

#include    "String2.cpp"
#include    "crc32_zlib.h"
#include    "common.cpp"
#include    "compress.cpp"

const char * szHelp ="\
  options:\n\
    [-f]                    file to handle\n\
    -n                      list number of frames\n\
    -l, -list[-frames]      list all the frames\n\
    -dim[ension]            list dimenstions (x y z v) of frames\n\
    -check                  check data frame and validate CRC\n\
    -d[ata], -e[xtract]     extract data of specific frame\n\
    -ex, -ey, -ez           extract all cutview frames according to x/y/z\n\
    -pwd, -epwd             extracting folder (for -ex/ey/ez)\n\
    -%...                   data format\n\
    -no-check-when-extract  don't check CRC when extract data\n\
  frame identifier:\n\
      by frame number:  a given frame number\n\
      by frame content: title@time\n\
";

#define WORK_LIST_N_FRAMES      1
#define WORK_LIST_HEADERS_L     2
#define WORK_LIST_HEADERS       3
#define WORK_LIST_DIM           4
#define WORK_CHECK_DATA         5
#define WORK_EXTRACT_DATA       6
#define WORK_EXTRACT_EDATA      7
#define WORK_EXTRACT_EDATA_BY_X 701
#define WORK_EXTRACT_EDATA_BY_Y 702
#define WORK_EXTRACT_EDATA_BY_Z 703


void print_uncompress_page_error(FILE * fout, int uret, int nframe){
    if (uret==0){         fprintf(fout, "%s : %d at #%d : unknown error\n", software_name, uret, nframe);
    } else if (uret==-1){ fprintf(fout, "%s : %d at #%d : memory overflow\n", software_name, uret, nframe);
    } else if (uret==-2){ fprintf(fout, "%s : %d at #%d : unknown precision\n", software_name, uret, nframe);
    } else if (uret==-3){ fprintf(fout, "%s : %d at #%d : compressor unmatched\n", software_name, uret, nframe);
    } else if (uret==-4){ fprintf(fout, "%s : %d at #%d : uncompressor version too low\n", software_name, uret, nframe);
    } else if (uret==-5){ fprintf(fout, "%s : %d at #%d : unknown compressor\n", software_name, uret, nframe);
    }
}


int main(int argc, char * argv[]){
    bool success = true;
    int work = WORK_LIST_HEADERS; const char * identifier_string = nullptr;
    int identifier_frame = 0; double identifier_time = 0; StringNS::string identifier_title = "";
    char * filename = nullptr; char * output_format = nullptr;
    bool check_crc_on_extract = true; const char * extract_folder = nullptr;

    if (argc<2){ printf("%s %s\n", software_name, software_version); return 0; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '-'){
            StringNS::string key = argv[i];
            if (key=="-h" || key=="--h" || key=="-help" || key=="--help"){
                printf("%s %s %s\n%s", software_name, software_version, copyright_string, szHelp); success = false;
            } else if (key=="-version" || key=="--version"){
                printf("%s\n", software_version); success = false;
            } else if (key=="-no-check-when-extract" || key=="--no-check-when-extract" || key=="-no_check_when_extract" || key=="--no_check_when_extract" || key=="-no-check-on-extract" || key=="--no-check-on-extract" || key=="-no_check_on_extract" || key=="--no_check_on_extract"){
                check_crc_on_extract = false;
            } else if (key=="-do-check-when-extract" || key=="--do-check-when-extract" || key=="-do_check_when_extract" || key=="--do_check_when_extract" || key=="-do-check-on-extract" || key=="--do-check-on-extract" || key=="-do_check_on_extract" || key=="--do_check_on_extract"){
                check_crc_on_extract = true;
            } else if (key=="-pwd" || key=="--pwd" || key=="-epwd" || key=="--epwd"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; extract_folder = argv[i]; }
            } else if (key=="-l" || key=="--l"){
                work = WORK_LIST_HEADERS_L;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-list" || key=="--list" || key=="-list-frames" || key=="--list-frames" || key=="-list_frames" || key=="--list_frames"){
                work = WORK_LIST_HEADERS;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-n" || key=="--n" || key=="-number-of-frames" || key=="--number-of-frames" || key=="-number_of_frames" || key=="--number_of_frames"){
                work = WORK_LIST_N_FRAMES;
            } else if (key=="-dim" || key=="--dim" || key=="-dimension" || key=="--dimension"){
                work = WORK_LIST_DIM;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-check" || key=="--check"){
                work = WORK_CHECK_DATA;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-d" || key=="--d" || key=="-data" || key=="--data"){
                work = WORK_EXTRACT_DATA;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-e" || key=="--e" || key=="-extract" || key=="--extract"){
                work = WORK_EXTRACT_EDATA;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-ex" || key=="--ex" || key=="-extract-x" || key=="--extract-x" || key=="-extract_x" || key=="--extract_x"){
                work = WORK_EXTRACT_EDATA_BY_X;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-ey" || key=="--ey" || key=="-extract-y" || key=="--extract-y" || key=="-extract_y" || key=="--extract_y"){
                work = WORK_EXTRACT_EDATA_BY_Y;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-ez" || key=="--ez" || key=="-extract-z" || key=="--extract-z" || key=="-extract_z" || key=="--extract_z"){
                work = WORK_EXTRACT_EDATA_BY_Z;
                if (i+1<argc && argv[i+1][0]!='-'){ i++; identifier_string = argv[i]; }
            } else if (key=="-f" || key=="--f"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; filename = argv[i]; }
            } else if (key.text[0]=='-' && key.text[1]=='%'){
                output_format = &key.text[1];
            } else if (key.text[0]=='-' && key.text[1]=='-' && key.text[2]=='%'){
                output_format = &key.text[2];
            } else {
                fprintf(stderr, "%s : error : unrecognized parameter %s\n", software_name, argv[i]); success = false;
            }
        } else {
            filename = argv[i];
        }
    }
    if (!success) return 0;

    if (identifier_string){
        StringNS::string sline_identifier = identifier_string;
        if (StringNS::is_string_number(sline_identifier)){
            identifier_frame = atoi(sline_identifier.text);
        } else {
            identifier_title = identifier_string;
            for (int i=0; i<identifier_title.length; i++) if (identifier_title.text[i]=='@' || identifier_title.text[i]==',') identifier_title.length = i;
            if (identifier_title.text[identifier_title.length]=='@'||identifier_title.text[identifier_title.length]==','){
                if (StringNS::is_string_number(&identifier_title.text[identifier_title.length+1])){
                    identifier_time = atoi(&identifier_title.text[identifier_title.length+1]);
                }
            }
        }
    }
    //printf("identifier: %d, %s[%d], %g\n", identifier_frame, identifier_title.text, identifier_title.length, identifier_time);

    if (!filename){ printf("%s : please specify the file to handle\n", software_name); success = false; }

    FILE * file = nullptr;
    if (success){
        file = fopen(filename, "r"); if (!file){
            fprintf(stderr, "%s : error : cannot opent %s\n", software_name, filename); success = false;
        }
    }

    unsigned char * cdata = nullptr; size_t cdata_size = 0; size_t cdata_max_size = 0;  // compressed data
    unsigned char * data = nullptr; size_t data_size = 0; size_t data_max_size = 0; // uncompressed data
    CompressPageHeader header; memset(&header, 0, sizeof(header));

    unsigned short compressor_ver[5]; memset(&compressor_ver, 0, sizeof(compressor_ver)); compressor_ver[0] = 'u'+'n'*256;
  #ifdef _LIBZ_
    compressor_ver[0] = 'z' + 'l'*256;
    get_version_array_from_string(&compressor_ver[1], zlibVersion());
  #endif

    if (success){ int nframe = 0; bool data_alread_extracted = false;
        if (success){
            char file_type[4]; memset(file_type, 0, sizeof(file_type));
            fseek(file, 0, SEEK_SET); size_t br = fread(file_type, sizeof(file_type), 1, file); fseek(file, 0, SEEK_SET);
            if (memcmp(&file_type[0], szCompressPageHeader, 4)){
                fprintf(stderr, "%s : error : unknown file format \"%c%c%c%c\"\n", software_name, file_type[0], file_type[1], file_type[2], file_type[3]);
                success = false;
            }
        }
        if (work!=0){
            fseek(file, 0, SEEK_SET); while (success){ char text_buffer[4096]; memset(text_buffer, 0, sizeof(text_buffer));
                if (!fread(&header, sizeof(header), 1, file)) break;
                nframe ++;
                size_t text_size = header.data_offset - sizeof(header);
                if (text_size>0){
                    size_t read_size = text_size>sizeof(text_buffer)-1?sizeof(text_buffer)-1:text_size;
                    if (!fread(text_buffer, read_size, 1, file)) break;
                    fseek(file, -read_size, SEEK_CUR);
                }

                if (memcmp(&header.identifier[0], szCompressPageHeader, 4)){
                    fprintf(stderr, "%s : error : unknown format \"%c%c%c%c\" at frame %d\n", software_name, header.identifier[0], header.identifier[1], header.identifier[2], header.identifier[3], nframe);
                    success = false;
                }; if (!success) break;
                bool pass_filter = true;
                    if (identifier_frame>0 && nframe!=identifier_frame) pass_filter = false;
                    if (identifier_time>0 && fabs((float)(identifier_time-header.time_stamp)) > MACHINE_REASONABLE_ERROR) pass_filter = false;
                    if (identifier_title.text[0]){
                        //printf("identifier: %c%c%c%c vs %c%c%c%c\n", identifier_title.text[0], identifier_title.text[1], identifier_title.text[2], identifier_title.text[3], header.title[0], header.title[1], header.title[2], header.title[3]);
                        if (identifier_title!="any" && identifier_title!="all"){
                            for (int i=0; i<4 && header.title[i] && identifier_title.text[i]; i++) if (header.title[i] != identifier_title.text[i]) pass_filter = false;
                        }
                    }
              // do work
                bool data_read = false;
                if (pass_filter){
                    bool ts4sdump_can_decode = true; bool is_stdout_tty = isatty(fileno(stdout));
                    int bit_precision = header.data_original_length/(header.dimensions[0]*header.dimensions[1]*header.dimensions[2]*header.dimensions[3]);
                    if (work==WORK_LIST_HEADERS || work==WORK_LIST_HEADERS_L){
                      // display all headers
                        char print_buffer[4][128]; char title[5]; title[4] = 0; memcpy(title, header.title, 4);
                        if (header.dimensions[0]<=1 && header.dimensions[1]<=1 && header.dimensions[2]<=1 && header.dimensions[3]<=1){
                            printf("%d %4s@%g", nframe, title, header.time_stamp);
                        } else {
                            if (work==WORK_LIST_HEADERS_L){
                                    printf("%d %4s@%g, %s:%dx%dx%dx%d", nframe, title, header.time_stamp, bit_precision==sizeof(float)?"real4":bit_precision==sizeof(double)?"real8":"unknown precision", header.dimensions[0], header.dimensions[1], header.dimensions[2], header.dimensions[3]);
                            } else {
                                printf("%d %4s@%g, %s:%dx%dx%dx%d, CRC32: 0x%08X", nframe, title, header.time_stamp, bit_precision==sizeof(float)?"real4":bit_precision==sizeof(double)?"real8":"unknown precision", header.dimensions[0], header.dimensions[1], header.dimensions[2], header.dimensions[3], header.crc32);
                            }

                            if (work==WORK_LIST_HEADERS){
                                if ((header.compressor_version[0]&0xFF)=='u'){
                                    printf(", %s", print_memory_value(print_buffer[0], sizeof(print_buffer[0]), header.data_length));
                                } else {
                                    bool ts4sdump_can_decode = true;
                                    if (compressor_version_compare(header.compressor_version, compressor_ver)!=0) ts4sdump_can_decode = false;
                                    if (((header.compressor_version[0]&0xFF)!='u') && (compressor_ver[1]<header.compressor_version[1] || compressor_ver[2]<header.compressor_version[2] || compressor_ver[3]<header.compressor_version[3])) ts4sdump_can_decode = false;
                                    printf(", %s (%s<- %s in %g KB, %s%s)", print_memory_value(print_buffer[0], sizeof(print_buffer[0]), header.data_length), ts4sdump_can_decode?"":is_stdout_tty?"\33[31m":"", print_memory_value(print_buffer[1], sizeof(print_buffer[1]), header.data_original_length), header.compressor_page_size/1024.0, print_percentage_value(print_buffer[3], sizeof(print_buffer[3]), header.data_length/(double)header.data_original_length), ts4sdump_can_decode?"":is_stdout_tty?"\33[0m":"");
                                }
                            }
                        }
                        if (text_size>0){
                            //printf("%s  # %s%s\n", is_stdout_tty?"\33[37m":"", buffer, is_stdout_tty?"\33[0m":"");
                            #ifdef _TTYPROMPTCOLOR_
                                if (header.dimensions[0]<=1 && header.dimensions[1]<=1 && header.dimensions[2]<=1 && header.dimensions[3]<=1 && is_stdout_tty){
                                    printf(" # \33[1m%s\33[0m\n", text_buffer);
                                } else {
                                    printf(" # %s\n", text_buffer);
                                }
                            #else
                                printf(" # %s\n", text_buffer);
                            #endif
                        } else printf("\n");
                    } else if (work==WORK_LIST_DIM){
                        char title[5]; title[4] = 0; memcpy(title, header.title, 4);
                        printf("%d %4s@%g %4d %4d %4d %4d\n", nframe, title, header.time_stamp, header.dimensions[0], header.dimensions[1], header.dimensions[2], header.dimensions[3]);
                    } else if (work==WORK_CHECK_DATA||(!data_alread_extracted && (work==WORK_EXTRACT_DATA || work==WORK_EXTRACT_EDATA || work==WORK_EXTRACT_EDATA_BY_X || work==WORK_EXTRACT_EDATA_BY_Y || work==WORK_EXTRACT_EDATA_BY_Z))){
                      // relocate memory
                        if (data_max_size < header.data_original_length){
                            data_max_size = header.data_original_length;
                            if (data) free(data); data = (unsigned char*) malloc(data_max_size); if (!data){ fprintf(stderr, "%s : malloc failure\n", software_name); break; }
                        }
                        if (cdata_max_size < header.data_length){
                            cdata_max_size = header.data_length;
                            if (cdata) free(cdata); cdata = (unsigned char*) malloc(cdata_max_size); if (!cdata){ fprintf(stderr, "%s : malloc failure\n", software_name); break; }
                        }
                      // read file
                        if (header.data_offset-sizeof(header)>0) fseek(file, header.data_offset-sizeof(header), SEEK_CUR);
                        if (!fread(cdata, header.data_length, 1, file)){
                            fprintf(stderr, "%s : error reading #%d\n", software_name, nframe); data_read = false;
                            cdata_size = 0;
                        } else {
                            cdata_size = header.data_length;
                            data_read = true;
                        }
                      // do work
                        if (cdata_size>0){
                            if (work==WORK_CHECK_DATA){
                                char print_buffer[3][128]; char title[5]; title[4] = 0; memcpy(title, header.title, 4);
                                unsigned int crc32_this_frame = check_crc32(cdata, cdata_size);
                                printf("#%d %4s@%g, check CRC32: 0x%08X (%s0x%08X)\n", nframe, title, header.time_stamp, crc32_this_frame, crc32_this_frame==header.crc32?"==":"!=", crc32_this_frame);
                            } else if (work==WORK_EXTRACT_DATA || work==WORK_EXTRACT_EDATA){
                                unsigned int crc32_this_frame = header.crc32; if (check_crc_on_extract) crc32_this_frame = check_crc32(cdata, cdata_size);
                                if (crc32_this_frame != header.crc32){
                                    fprintf(stderr, "%s : CRC32 failure at #%d\n", software_name, nframe);
                                } else {
                                    int uret = uncompress_page(file, filename, nframe, &header, compressor_ver, cdata, cdata_size, data, data_max_size);
                                    if (uret>0){
                                        int nx = header.dimensions[0]; int ny = header.dimensions[1]; int nz = header.dimensions[2]; int nv = header.dimensions[3];
                                        float * data_pointer_float = (float *) data; double * data_pointer_double = (double *) data;
                                        for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++){
                                            if (work==WORK_EXTRACT_EDATA) printf("%5d %5d %5d ", ix+1, iy+1, iz+1);
                                            for (int iv=0; iv<nv; iv++){
                                                if (bit_precision==sizeof(float)){
                                                    printf(output_format?output_format:"%11g", data_pointer_float[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                } else if (bit_precision==sizeof(double)){
                                                    printf(output_format?output_format:"%11g", data_pointer_double[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                }
                                                printf(iv+1<nv?" ":"\n");
                                            }
                                        }
                                        //bit_precision
                                    } else {
                                        print_uncompress_page_error(stderr, uret, nframe);
                                    }; //printf("%d\n", uret);
                                    data_alread_extracted = true;
                                }
                            } else if (work==WORK_EXTRACT_EDATA_BY_X || work==WORK_EXTRACT_EDATA_BY_Y || work==WORK_EXTRACT_EDATA_BY_Z){
                                unsigned int crc32_this_frame = header.crc32; if (check_crc_on_extract) crc32_this_frame = check_crc32(cdata, cdata_size);
                                if (crc32_this_frame != header.crc32){
                                    fprintf(stderr, "%s : CRC32 failure at #%d\n", software_name, nframe);
                                } else {
                                    int uret = uncompress_page(file, filename, nframe, &header, compressor_ver, cdata, cdata_size, data, data_max_size);
                                    if (uret>0){
                                      // output file name
                                        const char * bn = basename(filename); StringNS::string bn_no_ext = file_without_extension(bn);
                                        char fn_prefix[MAX_PATH]; memset(fn_prefix, 0, sizeof(fn_prefix));
                                        bool extract_folder_end_with_slash = strlen(extract_folder)<1? false : extract_folder[strlen(extract_folder)-1]=='/'? true : false;
                                        if (extract_folder) snprintf(fn_prefix, sizeof(fn_prefix), extract_folder_end_with_slash?"%s":"%s/", extract_folder);
                                        memcpy(&fn_prefix[strnlen(fn_prefix, sizeof(fn_prefix))], bn_no_ext.text, bn_no_ext.length>MAX_PATH-1?MAX_PATH-1:bn_no_ext.length);
                                      // dimensions
                                        int nx = header.dimensions[0]; int ny = header.dimensions[1]; int nz = header.dimensions[2]; int nv = header.dimensions[3];
                                        float * data_pointer_float = (float *) data; double * data_pointer_double = (double *) data;
                                      // extract to each dimensions
                                        char fn_here[MAX_PATH];
                                        if (work==WORK_EXTRACT_EDATA_BY_X){
                                            for (int ix=0; ix<nx; ix++){
                                                snprintf(fn_here, sizeof(fn_here), nx<10?"%s_x%d.txt":nx<100?"%s_x%02d.txt":nx<1000?"%s_x%03d.txt":nx<10000?"%s_x%04d.txt":"%s_x%d.txt", fn_prefix, ix+1);
                                                FILE * file = fopen(fn_here, "w"); if (file){
                                                    for (int iz=0; iz<nz; iz++) for (int iy=0; iy<ny; iy++){
                                                        fprintf(file, "%d %d %d ", ix+1, iy+1, iz+1);
                                                        for (int iv=0; iv<nv; iv++){
                                                            if (bit_precision==sizeof(float)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_float[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            } else if (bit_precision==sizeof(double)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_double[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            }
                                                            fprintf(file, iv+1<nv?" ":"\n");
                                                        }
                                                    }
                                                    fclose(file);
                                                } else {
                                                    fprintf(stderr, "%s : error : cannot write to %s\n", software_name, fn_here);
                                                }
                                            }
                                        } else if (work==WORK_EXTRACT_EDATA_BY_Y){
                                            for (int iy=0; iy<ny; iy++){
                                                snprintf(fn_here, sizeof(fn_here), ny<10?"%s_y%d.txt":ny<100?"%s_y%02d.txt":ny<1000?"%s_y%03d.txt":ny<10000?"%s_y%04d.txt":"%s_y%d.txt", fn_prefix, iy+1);
                                                FILE * file = fopen(fn_here, "w"); if (file){
                                                    for (int iz=0; iz<nz; iz++) for (int ix=0; ix<nx; ix++){
                                                        fprintf(file, "%d %d %d ", ix+1, iy+1, iz+1);
                                                        for (int iv=0; iv<nv; iv++){
                                                            if (bit_precision==sizeof(float)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_float[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            } else if (bit_precision==sizeof(double)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_double[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            }
                                                            fprintf(file, iv+1<nv?" ":"\n");
                                                        }
                                                    }
                                                    fclose(file);
                                                } else {
                                                    fprintf(stderr, "%s : error : cannot write to %s\n", software_name, fn_here);
                                                }
                                            }
                                        } else if (work==WORK_EXTRACT_EDATA_BY_Z){
                                            for (int iz=0; iz<nz; iz++){
                                                snprintf(fn_here, sizeof(fn_here), nz<10?"%s_z%d.txt":nz<100?"%s_z%02d.txt":nz<1000?"%s_z%03d.txt":nz<10000?"%s_z%04d.txt":"%s_z%d.txt", fn_prefix, iz+1);
                                                FILE * file = fopen(fn_here, "w"); if (file){
                                                    for (int iy=0; iy<ny; iy++) for (int ix=0; ix<nx; ix++){
                                                        fprintf(file, "%d %d %d ", ix+1, iy+1, iz+1);
                                                        for (int iv=0; iv<nv; iv++){
                                                            if (bit_precision==sizeof(float)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_float[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            } else if (bit_precision==sizeof(double)){
                                                                fprintf(file, output_format?output_format:"%g", data_pointer_double[ix+(iy+(iz+iv*nz)*ny)*nx]);
                                                            }
                                                            fprintf(file, iv+1<nv?" ":"\n");
                                                        }
                                                    }
                                                    fclose(file);
                                                } else {
                                                    fprintf(stderr, "%s : error : cannot write to %s\n", software_name, fn_here);
                                                }
                                            }
                                        }
                                      //done
                                        printf("%s\n", fn_prefix);
                                    } else {
                                        print_uncompress_page_error(stderr, uret, nframe);
                                    };
                                    data_alread_extracted = true;
                                }
                            }
                        }
                    }
                }
              // next frame
                if (!data_read) fseek(file, header.data_offset + header.data_length - sizeof(header), SEEK_CUR);
            }
            if (work==WORK_LIST_N_FRAMES){
                printf("%d\n", nframe);
            }
        }
    }


    if (file) fclose(file); if (cdata) free(cdata); if (data) free(data);
    return success? 0 : -1;
}
