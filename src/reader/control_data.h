#define IMG_HDR_VERSION 3
#define IMG_MAGIC 0x494d4746

typedef struct imghdr 
{
  int img_magic;
  int img_hdr_length;
  int img_width;
  int img_height;
  int img_depth;
  int img_compress; /* should be 1 */
  int img_dwindow;
  int img_dlevel;
  int img_bgshade;
  int img_ovrflow;
  int img_undflow;
  int img_top_offset;
  int img_bot_offset;
  short img_version; /* IMG_HDR_VERSION */
  unsigned short img_checksum;
  int img_p_id;
  int img_l_id;
  int img_p_unpack;
  int img_l_unpack;
  int img_p_compress;
  int img_l_compresss;
  int img_p_histo;
  int img_l_histo;
  int img_p_text;
  int img_l_text;
  int img_p_graphics;
  int img_l_graphics;
  int img_p_dbHdr;
  int img_l_dbHdr;
  int img_levelOffset;
  int img_p_user;
  int img_l_user;
  int img_p_suite;
  int img_l_suite;
  int img_p_exam;
  int img_l_exam;
  int img_p_series;
  int img_l_series;
  int img_p_image;
  int img_l_image;
}
ImgHdr;

  
  
