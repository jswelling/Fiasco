/************************************************************
 *                                                          *
 *  dicom_parser.h                                          *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2004 Department of Statistics,         *
 *                        Carnegie Mellon University        *
 *                                                          *
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *                                                          *
 *                                                          *
 ************************************************************/

#ifndef INCL_DICOM_PARSER_H
#define INCL_DICOM_PARSER_H

typedef enum type_enum {
  UL, UI, SH, US, AE, AT, LO, OB, CS, DS, OW, OX, DL, SQ, PN, ST, LT,
  TM, DA, SS, IS, SL, AS, FD, DT, FL, XS, UT, OF, UNKNOWN, SPECIAL
} DicomElementType;

struct DataElement_struct;
struct DicomParser_struct;

/* This macro defines the shape of the methods we'll apply to elements.
 * Note that the parameter names are included!.
 */
#define DCM_METHOD_PROTOTYPE( foo ) \
  void foo(struct DicomParser_struct* parser, KVHash* info, \
	   struct DataElement_struct* de, \
	   FILE* f, const char* key, const char* def, void* hook)
typedef DCM_METHOD_PROTOTYPE( (*DataElementMethod) );

/* This macro defines the shape of the methods we'll use to test for
 * breaking out of the parsing loop..
 * Note that the parameter names are included!.
 */
#define DCM_PARSER_BREAK_TEST_PROTOTYPE( foo ) \
  int foo( KVHash* info, const struct DataElement_struct* de, \
           long long offset, void* hook )
typedef DCM_PARSER_BREAK_TEST_PROTOTYPE( (*ParserBreakTest) );

typedef struct DicomElementMethodTableEntry_stuct {
  long group;
  long element;
  const char* key;
  const char* def;
  DataElementMethod method;
} DicomElementMethodTableEntry;

typedef struct DicomDictDataElement_struct {
  long group;
  long element;
  DicomElementType type;
  const char* name;
} DicomDictDataElement;

/* For use in recognizing DICOM files */
#define DCM_SOP_CLASS_GROUP (0x8)
#define DCM_SOP_CLASS_ELEMENT (0x16)
#define DCM_DICOM_UID_STRING "1.2.840.10008.5.1.4."

/* What transfer syntax should we use for reading the
 * meta information header?
 */
#define DCM_META_INFO_TS_NAME "ExplicitVRLittleEndian"

typedef struct DataElement_struct {
  long group;
  long element;
  long length;
  long long payloadOffset;
  DicomElementType type;
  DicomDictDataElement* dictEntry;
  DicomElementMethodTableEntry* methodEntry;
} DicomDataElement;

typedef struct DicomParser_struct {
  int debug;
  DicomElementMethodTableEntry* elementMethodTable;
  long elementMethodTableSize;
} DicomParser;

typedef struct uid_struct {
  const char* name;
  const char* uid;
} UID;

typedef struct transfer_syntax_struct {
  const char* name;
  const char* UID;
  int isLittleEndian;
  int isExplicitVR;
  int isEncapsulated;
  int isBZip2ed;
  int isDeflated;
} TransferSyntax;

/* These functions may raise EXCEPTION_DICOM */
DicomParser* dcm_createParser(DicomElementMethodTableEntry* table, 
			      long tableSize);
void dcm_destroyParser(DicomParser* parser);
void dcm_setDebug(DicomParser* parser, const int val);
int dcm_getDebug(DicomParser* parser);
void dcm_parseStream(DicomParser* parser, KVHash* info, FILE* f, 
		     long long firstOffset, 
		     const TransferSyntax* transferSyntax,
		     ParserBreakTest breakTest,
		     void* hook);
const char* dcm_getElementTypeName(DicomElementType type);
const UID* dcm_getUIDByUIDString(const char* string);
const UID* dcm_getUIDByName(const char* name);
const TransferSyntax* dcm_getTransferSyntaxByUIDString(DicomParser* parser,
						       const char* uid);
const TransferSyntax* dcm_getTransferSyntaxByName(DicomParser* parser,
						  const char* name);
int dcm_transferSyntaxIsSupported(DicomParser* parser,
				  const TransferSyntax* ts);
const TransferSyntax* dcm_guessTransferSyntax(DicomParser* parser,
					      FILE* fp, long long offset);


#endif /* INCL_DICOM_PARSER_H */
