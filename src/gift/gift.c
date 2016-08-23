/*
 * gift.c - MRI image file format converter
 *
 * Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *
 * HISTORY
 *	1/96	- Written by Greg Hood
 */

/*
 *	To avoid the O(n^2) code complexity of converting from any of n formats
 *	to any other format, we use an intermediate format.  Thus, we convert
 *	the source into gift's internal format, then convert that into the
 *	destination format.  This reduces the code complexity to O(n).  In order
 *	to avoid an I/O penalty of actually writing the intermediate gift
 *	format out to disk, we keep it in memory.  Since this data may be large,
 *	we only keep the header and one image in memory at a time.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <dirent.h>
#include <unistd.h>
#include "bio.h"
#include "gift.h"
#include "fmri.h"
#include "fiasco.h"
#include "analyze.h"
#include "spiral_sl.h"
#include "pgh.h"
#include "stdcrg.h"

#define CMD_LENGTH	512
#define MAX_FORMATS	16

typedef struct Format {
  char name[32];
  Boolean (*check_format)();		/* returns TRUE if the files appear to be in the format */
  void (*start_reading)();		/* reads the header */
  void (*read_image)();			/* reads one image in */
  void (*end_reading)();		/* any cleanup from reading */
  void (*start_writing)();		/* usually writes the header */
  void (*write_image)();		/* writes one image out */
  void (*end_writing)();		/* any cleanup from writing */
} Format;

static char rcsid[] = "$Id: gift.c,v 1.12 2007/04/19 22:39:44 welling Exp $";

Format formats[MAX_FORMATS];		/* table of formats we can handle */
int n_formats;				/* number of formats */
FILE *input;				/* file for reading images */
FILE *output;				/* file for writing images */
int from_format;			/* format we are converting from (index into formats array) */
int to_format = -1;			/* format we are converting to   (index into formats array) */
Filename to_basename = "";		/* prefix for files we are writing */
int compress_limit = 10000000;		/* bytes to write out before we attempt to compress */
int uncompress_limit = 10000000;	/* bytes to try to uncompress at once */
int tmp_count = 1;			/* numbering of tmp files */
Boolean delete_input = FALSE;		/* if TRUE, delete the input files after conversion */
Boolean compress_output = FALSE;	/* if TRUE, compress the output after conversion */

Boolean big_endian_input;		/* if TRUE, assume input is big-endian */
Boolean big_endian_output;		/* if TRUE, write output in big-endian format */

int slice_start = -1;			/* which slices to convert */
int slice_end = -1;
int slice_stride = 1;

int time_start = -1;			/* which times to convert */
int time_end = -1;
int time_stride = 1;

/* the intermediate format */
GiftHeader fh;	/* header */
char *image;		/* one image */

extern Boolean analyze_write_split_times;
extern Boolean pgh_write_split_data;
extern int spiral_sl_read_coil;

/* FORWARD DECLARATIONS */
FileList ReadArguments (int argc, char **argv);
void ParseRange (int *start, int *end, int *stride, char *s);
void PrintUsage ();
void Convert (char *name);
FileList ExpandBasename (Filename basename);
int DetermineFormat (Filename basename, FileList files, int to_format);
void SortFiles (FileList *pHead, FileList *pTail);
int SortFunction (const void *p1, const void *p2);
void Cleanup ();
void Init ();
void AddFormat ();

int main (int argc, char **argv)
{
  FileList basenames, b;

  Init();
  basenames = ReadArguments(argc, argv);

  /* convert every basename given on the command line */
  for (b = basenames; b != NULL; b = b->next)
    Convert(b->name);

  return(0);
}

/* Init registers all the different format handlers */
void
Init ()
{
  AddFormat("FIASCO", FiascoCheckFormat,
	    FiascoStartReading, FiascoReadImage, FiascoEndReading,
	    FiascoStartWriting, FiascoWriteImage, FiascoEndWriting);
  AddFormat("PGH", PghCheckFormat,
	    PghStartReading, PghReadImage, PghEndReading,
	    PghStartWriting, PghWriteImage, PghEndWriting);
  AddFormat("ANALYZE", AnalyzeCheckFormat,
	    AnalyzeStartReading, AnalyzeReadImage, AnalyzeEndReading,
	    AnalyzeStartWriting, AnalyzeWriteImage, AnalyzeEndWriting);
  AddFormat("SPIRAL_SL", SpiralSLCheckFormat,
	    SpiralSLStartReading, SpiralSLReadImage, SpiralSLEndReading,
	    SpiralSLStartWriting, SpiralSLWriteImage, SpiralSLEndWriting);
}

void
AddFormat (char *name,
	   Boolean (*check_format_func)(),
	   void (*start_reading_func)(),
	   void (*read_image_func)(),
	   void (*end_reading_func)(),
	   void (*start_writing_func)(),
	   void (*write_image_func)(),
	   void (*end_writing_func)())
{
  strcpy(formats[n_formats].name, name);
  formats[n_formats].check_format = check_format_func;
  formats[n_formats].start_reading = start_reading_func;
  formats[n_formats].read_image = read_image_func;
  formats[n_formats].end_reading = end_reading_func;
  formats[n_formats].start_writing = start_writing_func;
  formats[n_formats].write_image = write_image_func;
  formats[n_formats].end_writing = end_writing_func;
  ++n_formats;
}

/* ReadArguments scans the command-line arguments, records what
   processing options the user has selected, and returns a list of
   the file basenames to be converted. */
FileList
ReadArguments (int argc,
	       char **argv)
{
  FileList fl_head, fl_tail;
  int a;
  int i;
  char *s;
  char slice_string[256]; /* no good way to handle length */
  char out_format_string[256]; /* no good way to handle length */
  char time_string[256]; /* still no good way to handle length */
  Filename in_fname; /* still no good way to handle length */
  
  big_endian_input = bio_big_endian_machine;
  big_endian_output = bio_big_endian_machine;
  fl_head = NULL;
  fl_tail = NULL;

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
        Help( "selecttopic" );
      else
        Help( argv[2] );
    }

 /*** Parse command line and get input params ***/
  cl_scan( argc, argv );

  big_endian_input= cl_present("bei");
  big_endian_output= cl_present("beo");
  if (cl_present("lei")) big_endian_input= FALSE;
  if (cl_present("leo")) big_endian_output= FALSE;
  compress_output= cl_present("compress");
  delete_input= cl_present("d|delete");
  cl_get( "coil", "%option %d", &spiral_sl_read_coil );
  cl_get( "o|outname", "%option %s", to_basename );
  if (cl_get( "slice", "%option %s", slice_string )) {
    ParseRange(&slice_start, &slice_end, &slice_stride, slice_string);
  }
  if (cl_present("single")) {
    analyze_write_split_times = FALSE;
    pgh_write_split_data = FALSE;
  }
  if (cl_present("split")) {
    analyze_write_split_times = TRUE;
    pgh_write_split_data = TRUE;
  }
  if (cl_get( "to", "%option %s", out_format_string )) {
    for (i = 0; i < n_formats; ++i)
      if (strcasecmp(out_format_string, formats[i].name) == 0)
	break;
    if (i >= n_formats)
      Abort("-to was followed by invalid format name\n");
    to_format = i;
  };
  if (cl_get( "time", "%option %s", time_string )) {
    ParseRange(&time_start, &time_end, &time_stride, time_string);
  }

  while (cl_get( "", "%s", &in_fname)) {
    AppendToFileList( &fl_head, &fl_tail, in_fname, 0 );
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  if (fl_head == NULL) {
    fprintf(stderr,"%s: error: required input file(s) not given\n",argv[0]);
    PrintUsage();
  }

  if (to_format < 0)
    Abort("Output format must be specified by \"-to <format>\"\n");

  return(fl_head);
}

/* ParseRange scans a start:end:stride argument on the command line */
void
ParseRange (int *start,
	    int *end,
	    int *stride,
	    char *string)
{
  char *s;
  int n;

  s = string;
  if (sscanf(s, "%d%n", start, &n) != 1)
    Abort("Range %s does not begin with a number\n", string);
  s = &s[n];
  if (*s == '\0')
    {
      *end = *start;
      *stride = 1;
      return;
    }
  if (*s != ':' && *s != ',' && *s != '-')
    Abort("Invalid range: %s\n", string);
  ++s;
  if (sscanf(s, "%d%n", end, &n) != 1)
    Abort("Range %s does not have a valid ending value.\n", string);
  s = &s[n];
  if (*s == '\0')
    {
      *stride = 1;
      return;
    }
  if (*s != ':' && *s != ',' && *s != '-')
    Abort("Invalid range: %s\n", string);
  ++s;
  if (sscanf(s, "%d", stride) != 1)
    Abort("Range %s does not have a valid stride value.\n", string);
}

/* PrintUsage prints a synopsis of the options available when invoking the program. */
void
PrintUsage ()
{
  Help("usage");

  exit(1);
}

/* Convert the set of files specified by basename into the desired format */
void
Convert (Filename basename)
{
  FileList files;
  int t;
  int s;

  /* determine what specific files match the basename */
  files = ExpandBasename(basename);
  if (files == NULL)
    Abort("No files match the given argument: %s\n", basename);

  /* figure out what input format we are dealing with */
  from_format = DetermineFormat(basename, files, to_format);
  if (from_format < 0)
    Abort("Input %s is not in a recognized format.\n",
	  basename);
  if (to_basename[0] == '\0')
    strcpy(to_basename, basename);

  /* start the conversion process */
  (*formats[from_format].start_reading)(basename, files);
  (*formats[to_format].start_writing)(to_basename);

  /* convert each image in the data set that we are interested in */
  for (t = fh.dim[4].min; t <= fh.dim[4].max; t += fh.dim[4].stride)
    for (s = fh.dim[3].min; s <= fh.dim[3].max; s += fh.dim[3].stride)
      {
	(*formats[from_format].read_image)(t, s);
	(*formats[to_format].write_image)(t, s);
      }

  /* finish up */
  (*formats[to_format].end_writing)();		/* we finish writing before reading because if the -delete
						   flag is given, we want to be sure valid output is present
						   before deleting the input */
  (*formats[from_format].end_reading)();
  Cleanup(files);
}

/* construct a sorted list of all filenames which begin with
   the specified basename */
FileList
ExpandBasename (Filename basename)
{
  char *p;
  DIR *d;
  struct dirent *e;
  Filename dir_name;
  Filename prefix;
  int prefix_len;
  FileList fl_head, fl_tail;
  Filename name;

  /* separate out the directory name from the file name prefix */
  strcpy(dir_name, basename);
  if ((p = strrchr(dir_name, '/')) != NULL)
    {
      *p = '\0';
      strcpy(prefix, p+1);
    }
  else
    {
      strcpy(dir_name, ".");
      strcpy(prefix, basename);
    }
  prefix_len = strlen(prefix);

  if ((d = opendir(dir_name)) == NULL)
    Abort("Cannot open directory %s\n", dir_name);

  fl_head = NULL;
  fl_tail = NULL;
  while ((e = readdir(d)) != NULL)
    if (strncmp(e->d_name, prefix, prefix_len) == 0)
      {
	sprintf(name, "%s/%s", dir_name, e->d_name);
	AppendToFileList(&fl_head, &fl_tail, name, 0);
      }

  closedir(d);
  SortFiles (&fl_head, &fl_tail);
  return(fl_head);
}

/* add a filename onto a list of files; size is optionally used to hold the size of the file in bytes */
void
AppendToFileList (FileList *pHead,
		  FileList *pTail,
		  char *name,
		  int size)
{
  FileList fl_new;

  fl_new = (FileList) malloc(sizeof(FileListElement));
  fl_new->next = NULL;
  StringCopy(fl_new->name, name, MAX_FILENAME_LENGTH);
  fl_new->size = size;
  if (*pTail != NULL)
    (*pTail)->next = fl_new;
  else
    (*pHead) = fl_new;
  (*pTail) = fl_new;
}

/* sort a file list into ASCII order */
void
SortFiles (FileList *pHead,
	   FileList *pTail)
{
  int i, count;
  FileListElement **sort_table, *fle;

  count = 0;
  for (fle = *pHead; fle != NULL; fle = fle->next)
    ++count;

  sort_table = (FileListElement **) malloc(count * sizeof(FileListElement *));
  i = 0;
  for (fle = *pHead; fle != NULL; fle = fle->next)
    sort_table[i++] = fle;
  qsort(sort_table, count, sizeof(FileListElement *), SortFunction);
  *pHead = NULL;
  *pTail = NULL;
  for (i = 0; i < count; ++i)
    {
      fle = sort_table[i];
      if (i == 0)
	*pHead = fle;
      if (i == (count-1))
	{
	  *pTail = fle;
	  fle->next = NULL;
	}
      else
	fle->next = sort_table[i+1];
    }
  free(sort_table);
}

int SortFunction (const void *p1,
		  const void *p2)
{
  FileListElement *fle1, *fle2;

  fle1 = *((FileListElement **) p1);
  fle2 = *((FileListElement **) p2);
  return(strcmp(fle1->name, fle2->name));
}

/* determine what format the files are already in */
int
DetermineFormat (Filename basename, FileList files, int to_format)
{
  int i;

  /* first check if the files are already in the desired output format */
  if (to_basename == '\0')
    {
      if (formats[to_format].check_format != NULL &&
	  (*formats[to_format].check_format)(basename, files))
	Abort("Input %s already appears to be in %s format!\n",
	      basename, formats[to_format].name);
      /* reverse the endian-ness and try again */
      big_endian_input = !big_endian_input;
      if (formats[to_format].check_format != NULL &&
	  (*formats[to_format].check_format)(basename, files))
	Abort("Input %s already appears to be in %s format!\n",
	      basename, formats[to_format].name);
      big_endian_input = !big_endian_input;
    }

  /* next try the other formats */
  for (i = 0; i < n_formats; ++i)
    if (formats[i].check_format != NULL &&
	(*formats[i].check_format)(basename, files))
      return(i);

  /* finally try the other formats assuming reversed-endian data */
  big_endian_input = !big_endian_input;
  for (i = 0; i < n_formats; ++i)
    if (formats[i].check_format != NULL &&
	(*formats[i].check_format)(basename, files))
      return(i);
  big_endian_input = !big_endian_input;

  return(-1);
}

/* HasExtension returns TRUE if string <name> ends in string <ext> */
Boolean
HasExtension (char *name,
	      char *ext)
{
  int name_len;
  int ext_len;

  name_len = strlen(name);
  ext_len = strlen(ext);
  if (name_len > ext_len && strcmp(&name[name_len - ext_len], ext) == 0)
    return(TRUE);
  return(FALSE);
}

/* ReadImage reads the data for one image into memory; it assumes
   the file pointer is positioned at the beginning of that image in the input file,
   and that the data is stored in xy order.  Any exceptions to this
   have to be handled outside of this routine. */
void
ReadImage (int time,
	   int slice)
{
  bio_error = 0;
  bio_big_endian_input = big_endian_input;
  switch (fh.data_type)
    {
    case GIFT_UINT8:
      FRdUInt8Array(input, (unsigned char *) image, fh.n_items_per_image);
      break;
    case GIFT_INT16:
      FRdInt16Array(input, (short *) image, fh.n_items_per_image);
      break;
    case GIFT_INT32:
      FRdInt32Array(input, (int *) image, fh.n_items_per_image);
      break;
    case GIFT_FLOAT32:
      FRdFloat32Array(input, (float *) image, fh.n_items_per_image);
      break;
    case GIFT_FLOAT64:
      FRdFloat64Array(input, (double *) image, fh.n_items_per_image);
      break;
    default:
      Abort("Invalid fh.data_type\n");
    }
  if (bio_error)
    Abort("Cannot read image (time %d, slice %d) in input file.\n",
	  time, slice);
}

/* WriteImage writes the data for one image out to disk; it assumes
   the file pointer is positioned correctly in the output file,
   and that the data should be stored in xy order.
   Any exceptions to this have to be handled outside of this routine. */
void
WriteImage (int time,
	    int slice)
{
  bio_error = 0;
  bio_big_endian_output = big_endian_output;
  switch (fh.data_type)
    {
    case GIFT_UINT8:
      FWrUInt8Array(output, (unsigned char *) image, fh.n_items_per_image);
      break;
    case GIFT_INT16:
      FWrInt16Array(output, (short *) image, fh.n_items_per_image);
      break;
    case GIFT_INT32:
      FWrInt32Array(output, (int *) image, fh.n_items_per_image);
      break;
    case GIFT_FLOAT32:
      FWrFloat32Array(output, (float *) image, fh.n_items_per_image);
      break;
    case GIFT_FLOAT64:
      FWrFloat64Array(output, (double *) image, fh.n_items_per_image);
      break;
    default:
      Abort("Invalid fh.data_type\n");
    }
  if (bio_error)
    Abort("Cannot write image (time %d, slice %d) to output file.\n",
	  time, slice);
}

/* Compress requests the compression of file <name> of size <size> bytes;
   this compression may not happen immediately, but may be buffered so
   that many files can be compressed at once.  Calling this routine
   with a name of "" forces all pending files to be compressed. */
void
Compress (Filename name,
	  int size)
{
  static int total_size = 0;
  static int compress_cmd_length = 11;
  static char compress_cmd[CMD_LENGTH+sizeof(Filename)+2] = "compress -f";

  if (name[0] != '\0')
    {
      compress_cmd[compress_cmd_length++] = ' ';
      strcpy(&compress_cmd[compress_cmd_length], name);
      compress_cmd_length += strlen(name);
      total_size += size;
    }
  if ((name[0] == '0' && compress_cmd_length > 11) || compress_cmd_length > CMD_LENGTH ||
      total_size > compress_limit)
    {
      system(compress_cmd);
      compress_cmd[11] = '\0';
      compress_cmd_length = 11;
      total_size = 0;
    }
}

/* Uncompress uncompresses file <name> immediately.  The name of the
   uncompressed file is returned in <unc_name>. */
void
Uncompress (Filename unc_name,
	    char *name)
{
  char com[MAX_FILENAME_LENGTH + 128];
  char* uncompress_cmd= NULL;
  int len;

  if (!(uncompress_cmd= getenv("F_UNCOMPRESS")))
    uncompress_cmd= "uncompress";

  sprintf(unc_name, "/tmp/gift.%d.%d", getpid(), tmp_count++);
  sprintf(com, "%s -c %s >%s", uncompress_cmd, name, unc_name);
  system(com);
}

/* UncompressBatch uncompresses at least one file and possibly many files.
   A list of files is given to it in *fl.  It uncompresses files from
   that list until the maximum byte count to uncompress is exceeded, or
   else it reaches the end of the list.  The uncompressed files are
   left concatenated together in a single file, whose name is returned
   in <unc_name>.  The *fl pointer is updated to point to the remaining
   files on the list that were not uncompressed.  The function returns
   the total size (in bytes) of the uncompressed files. */
int
UncompressBatch (Filename unc_name, FileList *fl)
{
  int total_size;
  int count;
  char uncompress_cmd[CMD_LENGTH+sizeof(Filename)+64];
  int uncompress_cmd_length;
  char* unc_cmd= NULL;

  if (!(unc_cmd=getenv("F_UNCOMPRESS")))
    unc_cmd= "uncompress";
  
  sprintf(unc_name, "/tmp/gift.%d.%d", getpid(), tmp_count++);
  sprintf(uncompress_cmd, "%s -f -c",unc_cmd);
  uncompress_cmd_length= strlen(uncompress_cmd);
  total_size = 0;
  count = 0;
  while (*fl != NULL && uncompress_cmd_length < CMD_LENGTH && total_size < uncompress_limit)
    {
      uncompress_cmd[uncompress_cmd_length++] = ' ';
      strcpy(&uncompress_cmd[uncompress_cmd_length], (*fl)->name);
      uncompress_cmd_length += strlen((*fl)->name);
      total_size += (*fl)->size;
      ++count;
      *fl = (*fl)->next;
    }
  if (count > 0)
    {
      uncompress_cmd[uncompress_cmd_length++] = ' ';
      uncompress_cmd[uncompress_cmd_length++] = '>';
      strcpy(&uncompress_cmd[uncompress_cmd_length], unc_name);
      if (system(uncompress_cmd) != 0)
	Abort("uncompress failed\n");
    }
  return(total_size);
}


/* Delete removes the named file from the filesystem. */
void
Delete (Filename name)
{
  unlink(name);
}

/* Cleanup cleans up data structures after a format conversion. */
void
Cleanup (FileList files)
{
  FileListElement *fle, *nfle;

  /* finish compressing output files */
  Compress("", 0);

  /* deallocate arrays */
  if (image != NULL)
    free(image);
  if (fh.corrupt != NULL)
    {
      free(fh.corrupt);
      fh.corrupt = NULL;
    }

  /* deallocate the list of filenames matching the basename */
  for (fle = files; fle != NULL; fle = nfle)
    {
      nfle = fle->next;
      free(fle);
    }
}
