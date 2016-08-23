/*
 * rdb.h - Byte offsets within RDB header of P-files
 *		Based upon the GE rdbm.h file when compiled
 *		on a Sun workstation.
 */

#define RDB_HDR_SIZE		39940	/* size of the header in bytes */

#define RDB_HDR_SCAN_DATE	16
#define RDB_HDR_SCAN_TIME	26

#define RDB_HDR_USER0		216
#define RDB_HDR_USER1		220
#define RDB_HDR_USER2		224
#define RDB_HDR_USER3		228
#define RDB_HDR_USER4		232
#define RDB_HDR_USER5		236
#define RDB_HDR_USER6		240
#define RDB_HDR_USER7		244
#define RDB_HDR_USER8		248
#define RDB_HDR_USER9		252
#define RDB_HDR_USER10		256
#define RDB_HDR_USER11		260
#define RDB_HDR_USER12		264
#define RDB_HDR_USER13		268
#define RDB_HDR_USER14		272
#define RDB_HDR_USER15		276
#define RDB_HDR_USER16		280
#define RDB_HDR_USER17		284
#define RDB_HDR_USER18		288
#define RDB_HDR_USER19		292

#define RDB_HDR_NSLICES		68
#define RDB_HDR_FRAME_SIZE	80
#define RDB_HDR_DAB_0_START_RCV	200
#define RDB_HDR_DAB_0_STOP_RCV	202
#define RDB_HDR_RAW_PASS_SIZE	116
#define RDB_HDR_GW_POINT1_0	10244		
#define RDB_HDR_GW_POINT2_0	10256
#define RDB_HDR_GW_POINT3_0	10268
#define RDB_HDR_GW_POINT_STEP	4	/* number of bytes for each entry
					   in the gw_point1, gw_point2, and
					   gw_point3 arrays */
#define RDB_HDR_A_STEP		40	/* number of bytes in each
					   RDB_SLICE_INFO_ENTRY */
#define RDB_HDR_ROTATION	176
#define RDB_HDR_TRANSPOSE	178

#define RDB_HDR_IM_SLTHICK 	38942	/* slice thickness (mm)*/
#define RDB_HDR_IM_SCANSPACING  39032 	/* spacing between scans (mm) */
#define RDB_HDR_IM_CTR_R	39046	/* center R coord of plane image */
#define RDB_HDR_IM_CTR_A	39050	/* center A coord of plane image */
#define RDB_HDR_IM_CTR_S	39054   /* center S coord of plane image */
#define RDB_HDR_IM_TLC_R	39070   /* R coord of top left corner */
#define RDB_HDR_IM_TLC_A	39074   /* A coord of top left corner */
#define RDB_HDR_IM_TLC_S	39078   /* S coord of top left corner */
#define RDB_HDR_IM_TRC_R	39082   /* R coord of top right corner */
#define RDB_HDR_IM_TRC_A	39086   /* A coord of top right corner */
#define RDB_HDR_IM_TRC_S	39090   /* S coord of top right corner */
#define RDB_HDR_IM_BRC_R	39094   /* R coord of bottom right corner */
#define RDB_HDR_IM_BRC_A	39098   /* A coord of bottom right corner */
#define RDB_HDR_IM_BRC_S	39102   /* S coord of bottom right corner */
#define RDB_HDR_IM_TR		39110	/* TR - pulse repetition time
					   (microseconds) (int32) */
#define RDB_HDR_IM_TE		39118	/* TE - pulse echo time
					   (microseconds) (int32) */
#define RDB_HDR_IM_NEX		39134	/* number of excitations */
#define RDB_HDR_IM_FLIP		39170	/* flip angle (degrees) (int16) */
#define RDB_HDR_IM_PSD		39224	/* pulse sequence number */
