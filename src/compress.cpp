//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//---------------------------  class CompressPageHeader  --------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

// Page structure:
  // CompressPageHeader data_header
  // char[] extra_data_goes_here
  // unsighed char[] page_data_goes here

const char * szCompressPageHeader = "TS64";
class IETSPageHeader {
  public:
    char identifier[4]; // must be "TS64"
    char title[4]; // can be: g/guv, h/huv/hsr/hlr, c/cuv/csr/clr, lj/coul, uuv/ulr, ...
    unsigned int crc32; // CRC32 of stored data
    size_t data_offset; // begin of data from begin of this header
    size_t data_length; // length of data
    size_t data_original_length; // length of original data, data_length=data_original_length implies no compression
    unsigned short compressor_version[4]; // {"un",0,0,0} means uncompressed; {"zl",1,2,11} means zlib v 1.2.11
    unsigned int compressor_page_size;
};

class CompressPageHeader : public IETSPageHeader {
  public:
    unsigned int dimensions[4]; // nx, ny, nz, nv: from low to high dimension
    double time_stamp; // frame number of frame time
};

void get_version_array_from_string(unsigned short array[4], const char * string){
    char string_buffer[128]; memset(string_buffer, 0, sizeof(string_buffer));
    strncpy(string_buffer, string, sizeof(string_buffer));
    for (int i=0; i<sizeof(string_buffer)&&string_buffer[i]; i++) if (string_buffer[i]=='.') string_buffer[i] = ' ';
    StringNS::string sl[4]; int nw = StringNS::analysis_line(string_buffer, sl, 4, true);
    for (int i=0; i<4; i++) array[i] = 0;
    for (int i=0; i<4&&i<nw; i++) array[i] = atoi(sl[i].text);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-------------------------------  read uncompress  -------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
int compressor_version_compare(unsigned short hd[4]/*data*/, unsigned short hc[4]/*this code*/){
  // return: 0: can decode; 1: version too low; -1: compressor unmached
    //printf("%c%c %d.%d.%d vs %c%c %d.%d.%d\n", hc[0]%256, hc[0]/256, hc[1], hc[2], hc[3], hd[0]%256, hd[0]/256, hd[1], hd[2], hd[3]);
    if ((hd[0]&0xFF)=='u' || (hd[0]&0xFF)=='U') return 0; // no compression of data
    else if (hd[0]!=hc[0]) return -1; // compressor unsupport
    else if (hc[1]>=hd[1] && hc[2]>=hd[2] && hc[3]>=hd[3]) return 0; // can decode
    else return 1; // canont decode
}

int uncompress_page(FILE * file, const char * filename, int block_index, CompressPageHeader * header, unsigned short uncompressor_ver[4], unsigned char * compressed, size_t length_compressed, unsigned char * data, size_t length_max_data){
  // return:
    //  positive: length of uncompressed data; 0: not doing anything
    //  -1: memory overflow; -2: unknown precision; -3: compressor not match; -4: compressor version too low; -5: unknown compressor
    int bit_precision = header->data_original_length/(header->dimensions[0]*header->dimensions[1]*header->dimensions[2]*header->dimensions[3]); if (bit_precision!=sizeof(float) && bit_precision!=sizeof(double)) return -2;
    if ((header->compressor_version[0]&0xFF)!='u' && (header->compressor_version[0]&0xFF)!='U'){
        int compressor_check = compressor_version_compare(header->compressor_version, uncompressor_ver);
        if (compressor_check<0) return -3; else if (compressor_check>0) return -4;
    }
    if ((header->compressor_version[0]&0xFF)=='u' || (header->compressor_version[0]&0xFF)=='U'){
        if (length_max_data<length_compressed) return -1;
        memcpy(data, compressed, length_compressed); return length_compressed;
    } else if (header->compressor_version[0]=='z'+'l'*256 || header->compressor_version[0]=='Z'+'L'*256){
      #ifdef _LIBZ_
        unsigned char * pc = compressed; unsigned char * pd = data;
        while (pd<data+header->data_original_length && pc<compressed+header->data_length){
          // uncompress this page
            int size_c_page = *(int*)pc; pc += 4;
            if (pc+size_c_page > compressed+header->data_length) return -1;
            size_t src_length = size_c_page;
            size_t dst_length = pd+header->compressor_page_size>data+header->data_original_length? data+header->data_original_length-pd : header->compressor_page_size;
            //printf("uncompress2: %lu <- %d\n", dst_length, size_c_page);
            //uncompress2(pd, &dst_length, pc, &src_length);
            uncompress(pd, &dst_length, pc, src_length);
          // locate to next page
            pc += size_c_page; pd += header->compressor_page_size;
        }
        return header->data_original_length;
      #else
        return -5;
      #endif
    } else return -5;
    return 0;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//------------------------------  write append data  ------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

size_t write_data_page(FILE * fo, const char * filename, double time_stamp, const char * title, FILE * flog, unsigned short dimensions[4], void * data, size_t data_size, StringNS::string comment, int compress_level, size_t page_size, void * _compressBuf, size_t _allocate_memory_size){
    if (!title) title = "";
    if (!fo) return 0;
    //double current_time = get_current_time_double();

  // fille header common fields
    CompressPageHeader header; memset(&header, 0, sizeof(header));
    for (int i=0; i<4&&i<sizeof(szCompressPageHeader); i++) header.identifier[i] = szCompressPageHeader[i];
    for (int i=0; i<4&&title[i]; i++) header.title[i] = title[i];
    header.data_offset = sizeof(header) + comment.length;
    header.data_original_length = data_size;
    header.compressor_page_size = data_size;
    for (int i=0; i<4; i++) header.dimensions[i] = dimensions[i];
    header.time_stamp = time_stamp;

  // compress data and do CRC
    header.data_length = 0; unsigned char * compressBuf = nullptr;
  #ifdef _LIBZ_
    if (compress_level!=0){
        header.compressor_page_size = page_size;
        size_t allocate_memory_size = data_size + (data_size/page_size+1)*page_size*sizeof(unsigned int);

        unsigned short version_number[4]; get_version_array_from_string(version_number, zlibVersion());
        header.compressor_version[0] = 'z'+'l'*256; for (int i=1; i<4; i++) header.compressor_version[i] = version_number[i-1];
        if (_compressBuf && _allocate_memory_size>=allocate_memory_size){
            compressBuf = (unsigned char *)_compressBuf;
        } else {
            compressBuf = (unsigned char *)malloc(allocate_memory_size);
        }
          #ifdef MAX_MEMORYS
            if (!compressBuf){ fprintf(stderr, "%s : malloc failure, totally %d bolcks\n", software_name, _memory_blk_total); exit(-1); }
          #else
            if (!compressBuf){ fprintf(stderr, "%s : malloc failure.\n", software_name); exit(-1); }
          #endif
        memset(compressBuf, 0, data_size);

        int level = compress_level==-1? Z_DEFAULT_COMPRESSION : compress_level==2? Z_BEST_COMPRESSION : compress_level==1? Z_BEST_SPEED : Z_NO_COMPRESSION;

      // compress by blocks
        size_t compressRet = 0; size_t last_compressedlen = 0;
        for (size_t begin=0; begin<data_size; begin+=page_size){
            unsigned int * this_blk_size = (unsigned int *) &compressBuf[last_compressedlen]; last_compressedlen += sizeof(unsigned int); header.data_length += sizeof(unsigned int);
            unsigned long int originallen = begin+page_size<data_size? page_size : data_size-begin;
            if (originallen+last_compressedlen>allocate_memory_size){ header.data_length = 0; compressRet = -1; break; }
            unsigned long int compressedlen = originallen;
            //printf("compress %d, size %lu\n", begin/page_size, originallen);
            compressRet += compress2(&compressBuf[last_compressedlen], &compressedlen, &((unsigned char*)data)[begin], originallen, level);
            header.data_length += compressedlen; last_compressedlen += compressedlen;
            *this_blk_size = (unsigned int) compressedlen;
        }

      // compress all together
        //unsigned long int compressedlen = data_size;
        //int compressRet = compress2(compressBuf, &compressedlen, (unsigned char*) data, (size_t)data_size, level);
        //header.data_length = compressedlen;
      // compress done
        //printf("compress2() at level %d, return %d, compressedlen %d from %d\n", level, compressRet, header.data_length, header.data_original_length);
    }
  #endif
    if (header.data_length==0 || !compressBuf){
        header.compressor_version[0] = 'u'+'n'*256; for (int i=1; i<4; i++) header.compressor_version[i] = 0;
        //strncpy(header.compressor_identifier, "uncompressed", sizeof(header.compressor_identifier));
        compressBuf = (unsigned char *)data;
        header.data_length = header.data_original_length;
    }
    header.crc32 = (unsigned int)crc32_zlib(0L, nullptr, 0); header.crc32 = (unsigned int)crc32_zlib(header.crc32, compressBuf, header.data_length);

  // write data
    fwrite(&header, sizeof(header), 1, fo);
    if (comment.length>0) fwrite(comment.text, comment.length, 1, fo);
    fwrite(compressBuf, header.data_length, 1, fo);
    fflush(fo);

  // report
    //current_time = (get_current_time_double() - current_time)/1000;
    char membuff[64]; char compbuff[64]; char tbuff[64];
    if (compress_level!=0 && header.data_length!=header.data_original_length){
        if (flog) fprintf(flog, "  save: %4s -> \"%s\", %s (%s), CRC32: 0x%08X\n", title, filename, print_memory_value(membuff, sizeof(membuff), header.data_length), print_percentage_value(compbuff, sizeof(compbuff), header.data_length/(double)header.data_original_length), header.crc32);
    } else {
        if (flog) fprintf(flog, "  save: %4s -> \"%s\", %s, CRC32: 0x%08X\n", title, filename, print_memory_value(membuff, sizeof(membuff), header.data_length), header.crc32);
    }

  // done and dispose all
  #ifdef _LIBZ_
    if (compressBuf && compressBuf!=data && compressBuf!=_compressBuf) free(compressBuf);
  #endif

    return header.data_offset + header.data_length;;
}
