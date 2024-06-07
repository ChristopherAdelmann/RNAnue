#ifndef BGZIP_H
#define BGZIP_H

#include <inttypes.h>
#include "gzip.h"
typedef struct {
    uint32_t lenCompressedData;
} bgzip_Header ;


bgzip_Header bgzip_Header_default();

int bgzip_extractBgzHeader(gzip_Header* gzipHeader, bgzip_Header *bgzipHeader);

int64_t bgzip_findLenUncompressedData(FILE* f, gzip_Header* gzipHeader,
                                   bgzip_Header *bgzipHeader);


#endif // BGZIP_H
