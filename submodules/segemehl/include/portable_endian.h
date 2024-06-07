
/*
 *
 *	portable_endian.h
 *  
 * 
 *  @author , 
 *  @company  
 *  @date 09/03/2018 08:00:05 CEST  
 *
 */

#ifndef PORTABLE_ENDIAN_H_
#define PORTABLE_ENDIAN_H_
#include <inttypes.h>
#include <stdint.h>

uint16_t bl_endianReverse16(uint16_t value);
uint32_t bl_endianReverse32(uint32_t value);
uint16_t bl_letoh16(uint16_t value) ;
uint32_t bl_letoh32(uint32_t value) ;
uint16_t bl_readLittleEndian16(uint16_t value);
uint32_t bl_readLittleEndian32(uint32_t value);

#ifdef __APPLE__
#include <machine/endian.h>
#include <libkern/OSByteOrder.h>

#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)
#define _le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)
#define _le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)

#define __BIG_ENDIAN    BIG_ENDIAN
#define __LITTLE_ENDIAN LITTLE_ENDIAN
#define __BYTE_ORDER    BYTE_ORDER

#else
#include <endian.h>
#endif
#endif /* PORTABLE_ENDIAN_H_ */

