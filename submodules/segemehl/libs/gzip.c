
#include <stdlib.h>
#include <sys/types.h>
//#include <endian.h>
#include "portable_endian.h"
#include "gzip.h"
#include <assert.h>

#define CHECK_BIT(var,pos) !!((var) & (1<<(pos)))


// Header which all gzip files must have (no optional fields)
typedef struct  {
    uint8_t id1;
    uint8_t id2;
    uint8_t compressionMethod;
    uint8_t flags;
    uint32_t mtime;
    uint8_t extraFlags;
    uint8_t operatingSystem;
    uint16_t lenExtraBytes;
} CommonHeader;



typedef union  {
    CommonHeader gzipHeader;
    char buf[sizeof (CommonHeader)];
} u_CommonHeader;



/// Try to read sizeof(gzip_u_CommonHeader) from f (current offset!)
/// and fill param header appropriately.
/// file position of f may be changed.
/// @return 0 on success, 1 if not enough bytes read and 2 if not a gzip file.
static int readCommonHeader(FILE *f, CommonHeader *header){
    u_CommonHeader u_header;
    size_t bytesRead = fread(u_header.buf, 1, sizeof(u_CommonHeader), f );
    if(bytesRead != sizeof (u_CommonHeader)){
      fprintf(stderr, "common header has a length of %ld instead of %ld - not good.\n", bytesRead, sizeof(u_CommonHeader));
      if(ferror(f)) perror("file read error");
      return -1;
    }
    *header = u_header.gzipHeader;

    if(header->id1 != 31 || header->id2 != 139){
        return 1; // not a gzip file
    }
    // gzip stores int's as little endian -> transform to host's endianness

    //uint32_t save = header->mtime;
    uint32_t tmp = bl_readLittleEndian32(header->mtime);
    //assert(save == header->mtime);
    //header->mtime = le32toh(header->mtime);
    //assert(header->mtime == tmp);
    header->mtime = tmp;

    //save = header->lenExtraBytes;
    uint16_t tmp16 = bl_readLittleEndian16(header->lenExtraBytes);
    //header->lenExtraBytes = le16toh(header->lenExtraBytes);
    //assert(header->lenExtraBytes == tmp16);
    header->lenExtraBytes = tmp16;
    return 0;
}

/// Fills the extraByes in header.
/// @return 0 on success, -1 on error
static int readExtraBytes(FILE *f, const
                             uint16_t lenExtraBytes,
                             gzip_Header *header){
    if(lenExtraBytes <= MAX_XTRA_BYTES){
         if (fread(header->extraBytes, 1, lenExtraBytes, f ) != lenExtraBytes) {
             fprintf(stderr, "readExtraBytes: too few bytes read\n");
             return -1;
         }
    } else {
        // we are not interested in the bytes, if > MAX_XTRA_BYTES
        if(fseek(f, lenExtraBytes, SEEK_CUR) != 0){
            fprintf(stderr, "readExtraBytes: fseek failed\n");
            return -1;
        }
    }
    return 0;
}

/// @return countof read bytes or -1 on err
static int readUntilNullTerminator(FILE* f){
    int readBytes=0;
    int c;
    while ((c = fgetc(f)) != EOF) {
        readBytes++;
        if(c == '\0'){
            return readBytes;
        }
    }
    return -1;
}

gzip_Header gzip_Header_default()
{
    gzip_Header header;
    header.offsetInFile = -1;
    header.lenTotal = 0;
    header.compressionMethod = 255;
    header.lenExtraBytes=0;

    return header;
}

/// Read the whole header from f into param header.
/// Only a few of the read fields are actually relevant in this program.
/// @return 0 on success,
///         > 0, if not a properly formed gzip file
///         < 0, on error
int gzip_readHeader(FILE *f, gzip_Header *header)
{
    header->offsetInFile = ftell(f);
    if(header->offsetInFile == -1){
        fprintf(stderr, "gzip_readHeader: ftell failed\n");
    }
 
    //fseek(f, 0, SEEK_END);
    //off_t end=ftell(f);
    //fseek(f, header->offsetInFile, SEEK_SET);

   // fprintf(stderr, "reading header with offset %ld, end:%ld, eof:%d\n", header->offsetInFile, end, feof(f));
    CommonHeader commonHeader;
    int commonReturn = readCommonHeader(f, &commonHeader);
    if(commonReturn != 0){
        return commonReturn;
    }
    header->compressionMethod = commonHeader.compressionMethod;
    if(CHECK_BIT(commonHeader.flags, 2)){
        // extra bytes set
        header->lenExtraBytes = commonHeader.lenExtraBytes;
        int readExtraRet = readExtraBytes(f, commonHeader.lenExtraBytes, header);
        if(readExtraRet != 0){
            return readExtraRet;
        }
    } else {
        header->lenExtraBytes = 0;
    }

    header->lenTotal = sizeof (commonHeader) + header->lenExtraBytes;

    // consume fname, comment, crc16.

    if(CHECK_BIT(commonHeader.flags, 3)){
        int readBytes = readUntilNullTerminator(f);
        if(readBytes == -1 ){
            fprintf(stderr, "gzip_readHeader: reading filename failed\n");
            return -1;
        }
        header->lenTotal += readBytes;
    }
    if(CHECK_BIT(commonHeader.flags, 4)){
        int readBytes = readUntilNullTerminator(f);
        if(readBytes == -1 ){
            fprintf(stderr, "gzip_readHeader: reading comment failed\n");
            return -1;
        }
        header->lenTotal += readBytes;
    }

    if(CHECK_BIT(commonHeader.flags, 1)){
        const int countCrcBytes = 2;
        if (fread(header->extraBytes, 1, countCrcBytes, f ) != countCrcBytes) {
            fprintf(stderr, "gzip_readHeader: reading crcBytes failed\n");
            return -1;
        }
        header->lenTotal += countCrcBytes;
    }
    return 0;
}


