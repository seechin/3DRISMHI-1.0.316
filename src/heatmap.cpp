const char * software_name = "heatmap";
const char * software_version = "0.296.1746";
const char * copyright_string = "(c) Cao Siqin";

#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#include    <string.h>
#pragma pack (1)

#include    "String2.cpp"

double check_boundary(double value, double inf, double sup){
    if (value < inf) return inf;
    if (value > sup) return sup;
    return value;
}

//#include    "bmp.h"
#ifndef WORD
    #define WORD unsigned short
#endif
#ifndef DWORD
    #define DWORD unsigned int
#endif
#ifndef LONG
    #define LONG unsigned int
#endif

typedef struct tagBITMAPFILEHEADER {
  WORD    bfType;
  DWORD   bfSize;
  WORD    bfReserved1;
  WORD    bfReserved2;
  DWORD   bfOffBits;
} BITMAPFILEHEADER; //, FAR *LPBITMAPFILEHEADER, *PBITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  DWORD biSize;
  LONG  biWidth;
  LONG  biHeight;
  WORD  biPlanes;
  WORD  biBitCount;
  DWORD biCompression;
  DWORD biSizeImage;
  LONG  biXPelsPerMeter;
  LONG  biYPelsPerMeter;
  DWORD biClrUsed;
  DWORD biClrImportant;
} BITMAPINFOHEADER; //, FAR *LPBITMAPINFOHEADER, *PBITMAPINFOHEADER;

#ifndef __REAL__
    #define __REAL__ double
#endif
#define MAX_WORD 100
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
__REAL__ ** init_matrix(size_t m, size_t n, int overflow_chars=0){
    size_t lenh = sizeof(__REAL__*) * m; size_t len = lenh + sizeof(__REAL__) * m * n + overflow_chars;
    char * p = (char*) malloc(len); memset(p, 0, len); __REAL__ * d = (__REAL__*)(p + lenh);
    __REAL__ ** a = (__REAL__**) p; for (size_t i=0; i<m; i++) a[i] = &d[i*n];
    return a;
}
void cp_matrix(__REAL__ ** src, __REAL__ ** dst, size_t m, size_t n){
    for (size_t i=0; i<m; i++) for (size_t j=0; j<n; j++) dst[i][j] = src[i][j];
}

const char * szHelp = "\
  options:\n\
    -f              data file\n\
    -c[olor]        color map file\n\
    -n[r] 1000x1000 dimensions of input data\n\
    -[image-]size   image_width image_height\n\
    -col            column of data to plot\n\
    -o out.bmp      filename for output.\n\
";

typedef struct tagRGBColorLegend {
    double value; unsigned char r, g, b;
} RGBColorLegend;
#define MAX_COLOR_LEGEND 1000

int main (int argc, char * argv[]){
    bool success = true;
    int image_width = 0; int image_height = 0;
    int data_dim_x = 0; int data_dim_y = 0;
    int icol = 0;
    const char * filename = NULL; const char * filecolor = NULL;
    const char * ofilename = "out.bmp";
    RGBColorLegend legend[MAX_COLOR_LEGEND]; int nlegend = 0;
    int debug = 0;

    if (argc<2){ printf("%s %s\n", software_name, software_version); return 0; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '-'){
            StringNS::string key = argv[i];
            if (key=="-h" || key=="--h" || key=="-help" || key=="--help"){
                printf("%s %s %s\n%s", software_name, software_version, copyright_string, szHelp); success = false;
            } else if (key=="-version" || key=="--version"){
                printf("%s\n", software_version); success = false;
            } else if (key=="-debug" || key=="--debug"){
                debug = 1;
            } else if (key=="-f" || key=="--f"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; filename = argv[i]; }
            } else if (key=="-c" || key=="--c" || key=="-color" || key=="--color"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; filecolor = argv[i]; }
            } else if (key=="-size" || key=="--size" || key=="-image-size" || key=="--image-size" || key=="-image_size" || key=="--image_size"){
                if (i+1<argc && argv[i+1][0]!='-'){
                    bool analysis_compact = false;
                    for (int j=0; argv[i+1][j] && !analysis_compact; j++) if (argv[i+1][j]==',' || argv[i+1][j]=='x' || argv[i+1][j]=='X') analysis_compact = true;
                    if (analysis_compact){ i++;
                        StringNS::string slp[2]; char sep[4] = { ',', 'x', 'X', 0 };
                        int nwp = StringNS::analysis_general_line(sep, argv[i], slp, 3);
                        if (nwp<2){
                            fprintf(stderr, "%s : %s[%d] : incorrect -nr: %s\n", software_name, "argv", i, argv[i]); success = false;
                        } else {
                            char slpt[2][64]; memset(slpt, 0, sizeof(slpt));
                            for (int i_=0; i_<2; i_++) memcpy(slpt[i_], slp[i_].text, slp[i_].length>sizeof(slpt[i_])-1?sizeof(slpt[i_])-1:slp[i_].length);
                            image_width = atoi(slpt[0]);
                            image_height = atoi(slpt[1]);
                        }
                    } else {
                        if (i+1<argc && StringNS::is_string_number(argv[i+1])) image_width = image_height = atoi(argv[++i]);
                    }
                    //printf("size: %d %d\n", image_width, image_height);
                }
            } else if (key=="-n" || key=="--n" || key=="-nr" || key=="--nr"){
                if (i+1<argc && argv[i+1][0]!='-'){
                    bool analysis_compact = false;
                    for (int j=0; argv[i+1][j] && !analysis_compact; j++) if (argv[i+1][j]==',' || argv[i+1][j]=='x' || argv[i+1][j]=='X') analysis_compact = true;
                    if (analysis_compact){ i++;
                        StringNS::string slp[2]; char sep[4] = { ',', 'x', 'X', 0 };
                        int nwp = StringNS::analysis_general_line(sep, argv[i], slp, 3);
                        if (nwp<2){
                            fprintf(stderr, "%s : %s[%d] : incorrect -nr: %s\n", software_name, "argv", i, argv[i]); success = false;
                        } else {
                            char slpt[2][64]; memset(slpt, 0, sizeof(slpt));
                            for (int i_=0; i_<2; i_++) memcpy(slpt[i_], slp[i_].text, slp[i_].length>sizeof(slpt[i_])-1?sizeof(slpt[i_])-1:slp[i_].length);
                            data_dim_x = atoi(slpt[0]);
                            data_dim_y = atoi(slpt[1]);
                        }
                    } else {
                        if (i+1<argc && StringNS::is_string_number(argv[i+1])) data_dim_x = data_dim_y = atoi(argv[++i]);
                    }
                    //printf("dim: %d %d\n", data_dim_x, data_dim_y);
                }
            } else if (key=="-col" || key=="--col" || key=="-column" || key=="--column"){
                if (i+1<argc && StringNS::is_string_number(argv[i+1])){ i++; icol = atoi(argv[i])-1; }
            } else if (key=="-o" || key=="--o"){
                if (i+1<argc && argv[i+1][0]!='-'){ i++; ofilename = argv[i]; }
            } else {
                fprintf(stderr, "%s : error : unrecognized parameter %s\n", software_name, argv[i]); success = false;
            }
        } else {
            fprintf(stderr, "%s : error : unrecognized parameter %s\n", software_name, argv[i]); success = false;
        }
    }
    if (!success) return 0;

    if (!filecolor || !filecolor[0]){
        fprintf(stderr, "%s : error : color file (-c) not specified\n", software_name);
        success = false;
    }
    if (!filename || !filename[0]){
        fprintf(stderr, "%s : error : data file (-f) not specified\n", software_name);
        success = false;
    }
    if (data_dim_x<=0 || data_dim_y<=0){
        fprintf(stderr, "%s : error : invalid dimensions of data: -n %d %d\n", software_name, data_dim_x, data_dim_y);
        success = false;
    }
    if (icol<0){
        fprintf(stderr, "%s : error : invalid column: -col %d\n", software_name, icol+1);
        success = false;
    }
    if (image_width<=0 || image_height<=0){
        image_width = image_height = 1000;
        if (data_dim_x>0 && data_dim_y>0){
            if (image_width/(double)image_height < data_dim_x/(double)data_dim_y){
                image_width = image_height * data_dim_x / data_dim_y;
            } else {
                image_height = image_width * data_dim_y / data_dim_x;
            }
            if (debug) fprintf(stderr, "%s : image size %dx%d\n", software_name, image_width, image_height);
        }
    }

  // read data legend file

    if (debug) printf("%s : read color file\n", software_name);

    if (success){
        FILE * file = fopen(filecolor, "r");
        if (!file){
            fprintf(stderr, "%s : error : cannot open the color file %s\n", software_name, filecolor); success = false;
        } else {
            int nline = 0; char input[40960]; StringNS::string sl[MAX_WORD];
            while (fgets(input, sizeof(input), file)){ nline++;
                int nw = analysis_line_params(input, sl, MAX_WORD, true); if (nw<=0) continue ;
                if (sl[0].text[0]=='#' || sl[0].text[0]=='@' || (sl[0].text[0]=='/'&&sl[0].text[1]=='/')) continue;
                if (nw<2) continue;
                if (nlegend>=MAX_COLOR_LEGEND) break;
                if (sl[1].text[0]=='0'&&(sl[1].text[1]=='x'||sl[1].text[1]=='X')){
                    legend[nlegend].value = atof(sl[0].text);
                    unsigned int color_hex = __StringNS__::string_to_hex(StringNS::string(&sl[1].text[2], sl[1].length-2));
                    legend[nlegend].r = (color_hex / 65536) % 256;
                    legend[nlegend].g = (color_hex / 256) % 256;
                    legend[nlegend].b = color_hex % 256;
                } else if (nw>=4){
                    legend[nlegend].value = atof(sl[0].text);
                    legend[nlegend].r = atoi(sl[1].text);
                    legend[nlegend].g = atoi(sl[2].text);
                    legend[nlegend].b = atoi(sl[3].text);
                }
                nlegend ++;
            }
            fclose(file);
        }
        //for (int i=0; i<nlegend; i++) printf("color item [%d] : %12g (%d,%d,%d)\n", i, legend[i].value, legend[i].r, legend[i].g, legend[i].b);
    }


  // read data file

    if (debug) printf("%s : read data file\n", software_name);

    double ** input_data = NULL; double ** plot_data = NULL;
    if (success){
        input_data = init_matrix(data_dim_x, data_dim_y);
        plot_data = init_matrix(image_width, image_height);
        if (!input_data||!plot_data){ fprintf(stderr, "malloc failure\n"); success = false; }
    }
    if (success){
        FILE * file = fopen(filename, "r");
        if (!file){
            fprintf(stderr, "%s : error : cannot open the data file %s\n", software_name, filename); success = false;
        } else {
            int nline = 0; int idata = 0; char input[40960]; StringNS::string sl[MAX_WORD];
            while (fgets(input, sizeof(input), file)){ nline++;
                int nw = analysis_line_params(input, sl, MAX_WORD, true); if (nw<=0) continue ;
                if (sl[0].text[0]=='#' || sl[0].text[0]=='@' || (sl[0].text[0]=='/'&&sl[0].text[1]=='/')) continue;

                if (nw<=icol) continue;
                if (idata>=data_dim_x*data_dim_y) break;
                int ix = idata % data_dim_x; int iy = idata / data_dim_x;
                input_data[ix][iy] = atof(sl[icol].text);
                idata++;

            }
            fclose(file);
        }
    }

    //for (int y=0; y<data_dim_y/5; y++) for (int x=0; x<data_dim_x/5; x++) printf("read data [%d][%d] = %g\n", x, y, input_data[x][y]);

    if (debug) printf("%s : interporlate data\n", software_name);

  // interporlate the output grids of data
    if (success){
        for (int ix=0; ix<image_width; ix++) for (int iy=0; iy<image_height; iy++){
            double dx = ix / (double)image_width * data_dim_x;
            double dy = iy / (double)image_height * data_dim_y;
            int dx1 = floor(dx); int dx2 = dx>dx1? dx+1 : dx1;
              dx1 = check_boundary(dx1, 0, data_dim_x-1);
              dx2 = check_boundary(dx2, 0, data_dim_x-1);
            int dy1 = floor(dy); int dy2 = dy>dy1? dy+1 : dy1;
              dy1 = check_boundary(dy1, 0, data_dim_y-1);
              dy2 = check_boundary(dy2, 0, data_dim_y-1);

            double wx1 = dx2==dx1? 1 : ((double)dx2-dx)/((double)dx2-dx1);
            double wx2 = 1 - wx1;
            double wy1 = dy2==dy1? 1 : ((double)dy2-dy)/((double)dy2-dy1);
            double wy2 = 1 - wy1;

            double v11 = input_data[dx1][dy1];
            double v12 = input_data[dx1][dy2];
            double v21 = input_data[dx2][dy1];
            double v22 = input_data[dx2][dy2];

            plot_data[ix][iy] = (v11*wx1 + v21*wx2) * wy1 + (v12*wx1 + v22*wx2) * wy2;

            //printf("interporlate %d %d from %d %d %d %d: %g %g %g %g (%g %g %g %g) -> %g\n", ix, iy, dx1, dx2, dy1, dy2, v11, v12, v21, v22, wx1, wx2, wy1, wy2, plot_data[ix][iy]);
        }
    }

    //for (int y=0; y<image_width/5; y++) for (int x=0; x<image_height/5; x++) printf("plot data [%d][%d] = %g\n", x, y, plot_data[x][y]);

    if (debug) printf("%s : write BMP\n", software_name);

  // write RGB data to BMP file
    unsigned char * rgb_data = NULL;
    if (success){

        BITMAPFILEHEADER head;
        BITMAPINFOHEADER info;

        int width = image_width; int height = image_height; int rgb24_align_gap = 0;
        if ((width*3)%4 != 0) rgb24_align_gap = 4 - (width*3)%4;
        int stride24 = width*3 + rgb24_align_gap;
        int size24 = stride24 * height;
        rgb_data = (unsigned char *)malloc(sizeof(unsigned char)*size24);
        if (!rgb_data){ fprintf(stderr, "malloc failure\n"); success = false; }

        head.bfType = 'B'+'M'*256;
        head.bfSize = sizeof(head) + sizeof(info) + size24;
        head.bfReserved1 = head.bfReserved2 = 0;
        head.bfOffBits = sizeof(head) + sizeof(info);

        if (debug) printf("data size: %d, %dx%d\n", size24, width, height);

        //printf("size of WORD: %d\n", (int)(((char *)&head.bfSize)-((char*)&head.bfType)));
        //for (int i=0; i<sizeof(head); i++) printf(" %02X", ((unsigned char *)&head)[i]); printf("\n");

        info.biSize = sizeof(tagBITMAPINFOHEADER);
        info.biWidth = width;
        info.biHeight = height;
        info.biPlanes = 1;
        info.biBitCount = 24;
        info.biCompression = 0;
        info.biSizeImage = 0;
        info.biXPelsPerMeter = 5670; // 144 DPI
        info.biYPelsPerMeter = 5670; // 144 DPI
        info.biClrUsed = 0;
        info.biClrImportant = 0;

        if (success) for (int x=0; x<width; x++){ for (int y=0; y<height; y++){
            unsigned int base = y*stride24+x*3;


            if (plot_data[x][y]<=legend[0].value){
                rgb_data[base+2] = legend[0].r;
                rgb_data[base+1] = legend[0].g;
                rgb_data[base+0] = legend[0].b;
            } else if (plot_data[x][y]>=legend[nlegend-1].value){
                rgb_data[base+2] = legend[nlegend-1].r;
                rgb_data[base+1] = legend[nlegend-1].g;
                rgb_data[base+0] = legend[nlegend-1].b;
            } else {
                int icolor1 = 0; int icolor2 = 0; double w1 = 1; double w2 = 0;
                for (int i=0; i<nlegend-1; i++) if (plot_data[x][y]>=legend[i].value && plot_data[x][y]<=legend[i+1].value){
                    icolor1 = i;
                    icolor2 = plot_data[x][y]==legend[i].value? i : i+1;
                    if (legend[icolor1].value==legend[icolor2].value){
                        w1 = 1; w2 = 0;
                    } else {
                        w1 = (legend[icolor2].value - plot_data[x][y]) / (double)(legend[icolor2].value - legend[icolor1].value);
                        w2 = 1 - w1;
                    }
                    break;
                }
                rgb_data[base+2] = legend[icolor1].r * w1 + legend[icolor2].r * w2;
                rgb_data[base+1] = legend[icolor1].g * w1 + legend[icolor2].g * w2;
                rgb_data[base+0] = legend[icolor1].b * w1 + legend[icolor2].b * w2;
                //printf("%d %d : %12g -> %d,%d,%d: with %d,%d,%d (%g,w%g) and %d,%d,%d (%g,w%g)\n",x, y, plot_data[x][y], rgb_data[base+0],rgb_data[base+1],rgb_data[base+2], legend[icolor1].r,legend[icolor1].g,legend[icolor1].b,legend[icolor1].value,w1,  legend[icolor2].r,legend[icolor2].g,legend[icolor2].b,legend[icolor2].value,w2);
            }

            //rgb_data[base+0] = 0;
            //rgb_data[base+1] = 0;
            //rgb_data[base+2] = (plot_data[x][y]+1)*255;
            //printf("%d %d : %12g -> %02X %02X %02X\n", x, y, plot_data[x][y], rgb_data[base+0], rgb_data[base+1], rgb_data[base+2]);

        }}

        FILE * fo = fopen(ofilename, "w"); if (fo){
            fwrite(&head, sizeof(head), 1, fo);
            fwrite(&info, sizeof(info), 1, fo);
            fwrite(rgb_data, size24, 1, fo);
            fclose(fo);
        } else {
            fprintf(stderr, "%s : error : cannot write BITMAP to %s\n", software_name, ofilename); success = false;
        }
    }

    if (rgb_data) free(rgb_data); if (plot_data) free(plot_data);
    if (input_data) free(input_data);
    return 0;
}
