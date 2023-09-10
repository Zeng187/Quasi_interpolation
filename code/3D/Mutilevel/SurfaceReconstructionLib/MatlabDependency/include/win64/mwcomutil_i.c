

/* this ALWAYS GENERATED file contains the IIDs and CLSIDs */

/* link this file in with the server and any clients */


 /* File created by MIDL compiler version 8.00.0603 */
/* at Tue Feb 06 01:45:20 2018
 */
/* Compiler settings for win64\mwcomutil.idl:
    Oicf, W1, Zp8, env=Win64 (32b run), target_arch=IA64 8.00.0603 
    protocol : dce , ms_ext, c_ext, robust
    error checks: allocation ref bounds_check enum stub_data 
    VC __declspec() decoration level: 
         __declspec(uuid()), __declspec(selectany), __declspec(novtable)
         DECLSPEC_UUID(), MIDL_INTERFACE()
*/
/* @@MIDL_FILE_HEADING(  ) */

#pragma warning( disable: 4049 )  /* more than 64k source lines */


#ifdef __cplusplus
extern "C"{
#endif 


#include <rpc.h>
#include <rpcndr.h>

#ifdef _MIDL_USE_GUIDDEF_

#ifndef INITGUID
#define INITGUID
#include <guiddef.h>
#undef INITGUID
#else
#include <guiddef.h>
#endif

#define MIDL_DEFINE_GUID(type,name,l,w1,w2,b1,b2,b3,b4,b5,b6,b7,b8) \
        DEFINE_GUID(name,l,w1,w2,b1,b2,b3,b4,b5,b6,b7,b8)

#else // !_MIDL_USE_GUIDDEF_

#ifndef __IID_DEFINED__
#define __IID_DEFINED__

typedef struct _IID
{
    unsigned long x;
    unsigned short s1;
    unsigned short s2;
    unsigned char  c[8];
} IID;

#endif // __IID_DEFINED__

#ifndef CLSID_DEFINED
#define CLSID_DEFINED
typedef IID CLSID;
#endif // CLSID_DEFINED

#define MIDL_DEFINE_GUID(type,name,l,w1,w2,b1,b2,b3,b4,b5,b6,b7,b8) \
        const type name = {l,w1,w2,{b1,b2,b3,b4,b5,b6,b7,b8}}

#endif !_MIDL_USE_GUIDDEF_

MIDL_DEFINE_GUID(IID, IID_IMWUtil,0xC47EA90E,0x56D1,0x11d5,0xB1,0x59,0x00,0xD0,0xB7,0xBA,0x75,0x44);


MIDL_DEFINE_GUID(IID, LIBID_MWComUtil,0x8A4A8B94,0x0DE0,0x4A7B,0xA6,0xBC,0x65,0x48,0x5C,0xEE,0x44,0xD0);


MIDL_DEFINE_GUID(CLSID, CLSID_MWField,0x26C540D8,0x6613,0x4527,0x99,0xC0,0xB6,0xE2,0xCE,0xEE,0x70,0x06);


MIDL_DEFINE_GUID(CLSID, CLSID_MWStruct,0xA6177DCF,0x37CD,0x4806,0xB2,0xAD,0x9A,0x88,0xF3,0x37,0x3A,0xE0);


MIDL_DEFINE_GUID(CLSID, CLSID_MWComplex,0xD1E247D7,0xF263,0x4F05,0x87,0x2F,0xE6,0xD7,0x79,0x8F,0x39,0xCA);


MIDL_DEFINE_GUID(CLSID, CLSID_MWSparse,0x556E64DA,0x46C2,0x45DE,0x9F,0x50,0x9D,0x02,0xA3,0x03,0x2A,0xC6);


MIDL_DEFINE_GUID(CLSID, CLSID_MWArg,0xB9EF3AB4,0x9E97,0x495D,0xA5,0x18,0x17,0x2E,0x7B,0x88,0xA8,0x54);


MIDL_DEFINE_GUID(CLSID, CLSID_MWArrayFormatFlags,0xE3240CC4,0xFBB2,0x4F65,0x86,0x1F,0xAD,0xC6,0x78,0x63,0x98,0xA9);


MIDL_DEFINE_GUID(CLSID, CLSID_MWDataConversionFlags,0x9AD155E5,0xB420,0x443D,0x90,0xBB,0x99,0x6D,0xDA,0x6C,0xE8,0xC9);


MIDL_DEFINE_GUID(CLSID, CLSID_MWUtil,0x9320093B,0x4CB8,0x44D9,0xB7,0x6A,0x1A,0x17,0x0C,0x3B,0x3A,0xF4);


MIDL_DEFINE_GUID(CLSID, CLSID_MWFlags,0xFC094E60,0x5DBA,0x4CE5,0x81,0xE3,0xC4,0x4A,0x93,0x8B,0x57,0xC5);

#undef MIDL_DEFINE_GUID

#ifdef __cplusplus
}
#endif



