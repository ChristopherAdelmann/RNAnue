
#include <assert.h>

#include <string.h>
//#include <endian.h>
#include "portable_endian.h"
#include "bgzip.h"

typedef struct {
    // next three fields are fixed for a bgzip header extension
    uint8_t id1;
    uint8_t id2;
    uint16_t len;
    uint16_t bsizeMinusOne;

} BgzipRawHeader ;

typedef union  {
    BgzipRawHeader header;
    char buf[sizeof (BgzipRawHeader)];
} u_BgzipRawHeader;



/// @return 0, on success 1, if gzipHeader does not belong to a bgzip-file
int bgzip_extractBgzHeader(gzip_Header* gzipHeader, bgzip_Header *bgzipHeader)
{
   
    if(gzipHeader->lenExtraBytes != sizeof (BgzipRawHeader)){
        return 1;
    }
    u_BgzipRawHeader u_bgzheader;
    // might be optimizable...
    memcpy(u_bgzheader.buf, gzipHeader->extraBytes, sizeof (u_BgzipRawHeader));

    if(u_bgzheader.header.id1 != 66 || u_bgzheader.header.id2 != 67){
        return 1; // not a bgzip file
    }

    uint16_t tmp16 = bl_readLittleEndian16(u_bgzheader.header.len);
    //u_bgzheader.header.len = le16toh(u_bgzheader.header.len);
    //assert(u_bgzheader.header.len == tmp16);
    u_bgzheader.header.len = tmp16; 
    if(u_bgzheader.header.len != 2){
        return 1;
    }
  
    tmp16 = bl_readLittleEndian16(u_bgzheader.header.bsizeMinusOne);
    //u_bgzheader.header.bsizeMinusOne = le16toh(u_bgzheader.header.bsizeMinusOne);
    u_bgzheader.header.bsizeMinusOne = tmp16;
    
    bgzipHeader->lenCompressedData = u_bgzheader.header.bsizeMinusOne;

    bgzipHeader->lenCompressedData-= (gzipHeader->lenExtraBytes + 19) ;
    return 0;
}

bgzip_Header bgzip_Header_default()
{
    bgzip_Header header;
    header.lenCompressedData = 0;
    return header;
}

/// The uncrompressed data is at the end of a block, so after this function was
/// called successfully, f points to the next block.
/// @return length of the uncompressed data or -1 on error. A well-formed
/// bgzip-file should end with a special block, whose uncompressed
/// length is 0, so this value is returned in that case.
int64_t bgzip_findLenUncompressedData(FILE *f, gzip_Header *gzipHeader,
                                       bgzip_Header *bgzipHeader)
{
    assert(gzipHeader->offsetInFile != -1);

    long target = gzipHeader->offsetInFile
            + gzipHeader->lenTotal
            + bgzipHeader->lenCompressedData
            +  sizeof (uint32_t) // CRC32
            ;
    if(fseek(f, target, SEEK_SET) != 0){
        fprintf(stderr, "bgzip_findLenUncompressedData: fseek failed\n");
        return -1;
    }
    union u_unint32_t{
        uint32_t val;
        char buf[sizeof (uint32_t)];
    } isize;

    if (fread(isize.buf, 1, sizeof (uint32_t), f ) != sizeof (uint32_t)) {
        fprintf(stderr, "bgzip_findLenUncompressedData: too few bytes read\n");
        return -1;
    }
    uint32_t tmp = bl_readLittleEndian32(isize.val);
    //isize.val = le32toh(isize.val);
    isize.val = tmp;

    return isize.val;
}


