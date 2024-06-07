#ifndef GZIP_H
#define GZIP_H

#include <stdint.h>
#include <sys/types.h>
#include <stdio.h>
#include "portable_endian.h"

// current we are only interested in xtra-fields
// of the bgzip file-format, which has 6 bytes in the extra field
#define MAX_XTRA_BYTES 6



// only contains header fields relevant for this program
typedef struct {
    off_t offsetInFile;
    size_t lenTotal;
    uint8_t compressionMethod;
    uint16_t lenExtraBytes;
    char extraBytes[MAX_XTRA_BYTES];

} gzip_Header;


/// Always use this function when generating a new gzip_Header.
gzip_Header gzip_Header_default();

int gzip_readHeader(FILE* f, gzip_Header* header);



#endif // GZIP_H
