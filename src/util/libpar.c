/*
 *	libpar.c	Parallel task manager
 *
 *	This module implements a single-master/multiple-slave style
 *	of parallel processing.  Both the master and worker programs
 *	(which typically are instances of the same code) should call
 *	ParallelProcess upon startup.  In this call they supply the
 *	routines to be used for task execution as well as communication.
 *
 *	The user-supplied (*master_task) procedure farms out individual
 *	tasks to workers by setting up a Task record and then calling
 *	par_delegate_task.  When a worker has completed a task, it sets up
 *	a Result record before returning from (*worker_task).
 *	(*master_result) is then automatically called on the master
 *	to process the values in the Result record.
 *
 *	
 *	Copyright (c) 1995,1996,1997,1999,2001,2004  Pittsburgh Supercomputing
 *                                                   Center
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
 *	11/95 - Written by Greg Hood
 *	6/96  - Added capability for workers to return results to master
 *              (ghood@psc.edu)
 *	12/96 - Added Wait call and made group names passed through environment
 *		variables rather than by arguments (ghood@psc.edu)
 *	11/99 - added par_broadcast_context() for sending large contexts
 *		(ghood@psc.edu)
 *	3/00  - fixed bug in HandleRestart (ghood@psc.edu)
 *      3/01  - Converted to use either PVM or MPI (ghood@psc.edu)
 *      2/04  - Fixed bug that caused ready/idle list corruption upon
 *              worker termination (SWang@psych.uic.edu, ghood@psc.edu,
 *              welling@stat.cmu.edu)
 *
 */

/*
 * POSSIBLE FUTURE OPTIMIZATIONS:
 *	The scheduling could be made more sophisticated so that tasks having
 *	    the same context would be preferentially assigned to the same
 *	    worker, in order to reduce context switching.
 *	The scheduling could be made predictive so that a task would not
 *	    be assigned to a historically slow worker near the end of the
 *	    run when a faster worker will likely soon be available.
 */

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <netdb.h>
#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "errors.h"
#include "errors.h"
#include "misc.h"
#include "acct.h"
#include "par.h"

#if defined(PVM)
#include "pvm3.h"
#elif defined(MPI)
#include "mpi.h"
#endif

#define WORKER_RECEIVE_TIMEOUT	300	/* seconds before a worker will timeout
					   waiting for a message from the
					   master */
#define RETRY_HOST_INTERVAL	300	/* seconds between retrying to start
					   workers on a host */
#define PENDING_WORKERS_TIMEOUT	60	/* seconds to wait for pending
					   workers to show up before exiting */
#define REPORT_INTERVAL		5	/* we will report every this
					   many seconds when waiting
					   for pending workers to show up */

#define REQUEST_MSG		1	/* worker -> master */
#define RESTART_MSG		2	/* worker -> master */
#define CONTEXT_MSG		3	/* master -> worker  */
#define TASK_MSG		4	/* master -> worker  */
#define TERMINATE_MSG		6	/* master -> worker  */
#define WORKER_EXIT_MSG		7	/* PVM daemon -> master */
#define HOST_ADDED_MSG		8	/* PVM daemon -> master,
					   master -> master */
#define HOST_DELETED_MSG	9	/* PVM daemon -> master */
#define BROADCAST_CONTEXT_MSG	10	/* master -> worker */
#define BROADCAST_ACK_MSG	11	/* worker -> master */

#undef MIN
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#undef MAX
#define MAX(x,y)	((x) > (y) ? (x) : (y))

static char rcsid[] = "$Id: libpar.c,v 1.22 2005/07/07 20:13:08 welling Exp $";

/* the following buffer type is only used if we aren't using PVM buffers */
typedef struct Buffer {
#if defined(PVM)
  int bufid;			/* the PVM buffer holding the context info */
#elif defined(MPI)
  unsigned char *buffer;	/* actual data in the buffer */
  int size;			/* size of the buffer in bytes */
  int position;			/* reading or writing position */
#else
  int dummy;			/* to keep some compilers happy */
#endif
} Buffer;

typedef struct Context {
  int use_count;		/* # of times this structure is referenced */
  Buffer buffer;		/* the context buffer */
} Context;

typedef struct Task {
  struct Task *next;		/* next on the task queue */
  struct Task *prev;		/* previous on the task queue */
  Context *context;	        /* context associated with the task */
  int number;			/* number of this task */
  Buffer buffer;		/* the task buffer */
} Task;

typedef struct WorkerState {
  int tid;			/* the TID of this worker; note that for
				   MPI-based compilations, this is equal
				   to the rank */
  Hostname host;		/* the name of the host this worker
				   is running on */
  Boolean ready;		/* TRUE if this worker is ready to accept
				   a task */
  Boolean idle;			/* TRUE if this worker is ready to accept
				   a task and has no tasks currently
				   assigned to it */
  int next_ready;		/* the next worker on the ready list
				   (-1 terminates list) */
  int next_idle;		/* the next worker on the idle list
				   (-1 terminates list) */
  Context* last_context;	/* the context that was last sent to the
				   worker */
  Task *first_task;		/* the first task (of at most 2) assigned to
				   this worker */
  Task *second_task;		/* the second task (of at most 2) assigned to
				   this worker */
} WorkerState;

static Context *current_context = NULL;

static Task *first_queued_task = NULL;
static Task *last_queued_task = NULL;

static WorkerState workers[PAR_MAX_WORKERS];

static Boolean par = FALSE;	/* TRUE if we are going to run in parallel */

static int prog_argc;		/* the # of arguments this program was
				   invoked with, less -PAR_ arguments */
static char **prog_argv;	/* the arguments this program was invoked
				   with, less -PAR_ arguments */
static int par_argc;            /* the # of -PAR_ arguments this program
				   was invoked with */
static char **par_argv;         /* the -PAR_ arguments this program was
				   invoked with */
static char prog_name[128];	/* our program name */
static char group_name[128];    /* the name identifying our pvm process
				   group */

static int rank = -1;		/* the rank (instance #) of this process within
				   the process group (0 indicates this is the
				   master) */

static int my_tid;		/* the TID of this process */
static int master_tid;		/* the TID of the master process */
static Hostname my_host_name;	/* the name of the host on which this process
				   is running */
static Boolean spawn;		/* TRUE if we must start up the worker
				   processes */
static Boolean closed_worker_set;
                                /* TRUE if no new workers may be added to
				   the set of active workers */

static Boolean context_changed = FALSE;	/* TRUE if the context has changed
					   since the last call to
					   par_delegate_task */

static Par_Task task_number = 0; /* the next task number to be assigned */
static int tasks_outstanding = 0;/* counts the # of tasks that have been
				    delegated but not yet finished by the
				    workers */

static int n_workers = 0;	/* actual number of workers that are running */
static int n_workers_pending = 0; /* number of workers that should be starting,
				     but we haven't heard from yet */
static int ready_workers = -1;	/* a list of ready workers (i.e., those workers
				   having 0 or 1 tasks currently assigned to
				   them); a -1 terminates this list */
static int idle_workers = -1;	/* a list of idle workers (i.e., those workers
				   having 0 tasks currently assigned to them);
				   a -1 terminates this list */

static int par_verbose = FALSE;	/* if TRUE, report on normal events such as
				   task delegation and receiving worker
				   results */
static int debug_workers = FALSE; /* if TRUE, invoke workers within a
				     debugger */

static unsigned char *task_completed = NULL;
                                /* a bit array that tells for each task
				   whether it has been completed (a 1
				   indicated it has) */
static int task_completed_size = 0;
                                /* the number of bytes in task_completed */
static int task_completed_first = 0;
                                /* the number of the first task represented
				   in the table; all tasks with ID's lower
				   than this have already completed */
static int task_completed_offset = 0;
                                /* the index into task_completed at which
				   the task_completed_first task can be
				   found */
static int broadcast_count = 0;	/* number of broadcast messsages that
				   have been sent */
static int broadcast_first_ack_count = -1;
                                /* count of the last acknowledge received from
				   the first process that this process relayed
				   the broadcast to */
static int broadcast_second_ack_count = -1;
                                /* count of the last acknowledge received from
				   the second process that this process relayed
				   the the broadcast to */
static int broadcast_first_forward_tid = -1; /* TID of the first process to
						forward a broadcast to */
static int broadcast_second_forward_tid = -1; /* TID of the second process to
						 forward a broadcast to */
static int broadcast_backward_tid = -1; /* TID of the process from which this
					   process receives broadcasts */
static int n_processes = 1;	/* number of entries in process_tids */
static int process_tids[PAR_MAX_WORKERS+1];
                                /* stores all TIDs for processes
				   involved in a broadcast */
static int my_process_index = 0; /* index into process_tids for my process */

static int n_hosts_pending = 0;	/* number of hosts in the host table which
				   for which max_workers > 0 but which we
				   have not yet heard from */

static struct timeval longest_task= {0,0}; 
                                /* run time of longest task so far */

#if defined(PVM)
typedef struct HostState {
  Hostname name;		/* official name of the host */
  Hostname pvm_name;		/* name to use for this host in PVM calls;
				   this may differ from name if the user
				   manually added it to the PVM configuration
				   under a synonym */
  int dtid;			/* TID of the PVM daemon on the host
				   (-1 indicates that the host is not in the
				   current PVM configuration) */
  int max_workers;		/* maximum number of workers we should
				   try to start up on the host; if this number
				   is negative, we do not automatically recruit
				   this host, but if the user manually adds it
				   to the PVM configuration, then we try to
				   start abs(max_workers) on it */
  int n_workers;		/* number of workers actually running on
				   the host */
} HostState;

static int n_hosts = 0;		/* number of hosts */
static HostState hosts[PAR_MAX_HOSTS];

static int add_host_time = 0;	/* time at which to try to add hosts */

#elif defined(MPI)
static unsigned char *out_buffer = NULL;  /* output buffer */
static int out_size = 0;	/* size of output buffer in bytes */
static int out_position = 0;	/* index into output buffer */

static unsigned char *in_buffer = NULL;
static int in_size = 0;		/* size of input buffer in bytes */
static int in_position = 0;	/* index into input buffer */

static int sizeof_byte;		/* size of packed length of individual types */
static int sizeof_short;
static int sizeof_int;
static int sizeof_long;
static int sizeof_float;
static int sizeof_double;

#endif

/* the functions which must be supplied by the user of the this module */
static void (*par_master_task)() = NULL;	/* required */
static void (*par_master_result)() = NULL;	/* optional */
static void (*par_worker_context)() = NULL;	/* optional */
static void (*par_worker_task)() = NULL;	/* required */
static void (*par_worker_finalize)() = NULL;	/* optional */
static void (*par_pack_context)() = NULL;	/* optional */
static void (*par_unpack_context)() = NULL;	/* optional */
static void (*par_pack_task)() = NULL;		/* optional */
static void (*par_unpack_task)() = NULL;	/* optional */
static void (*par_pack_result)() = NULL;	/* optional */
static void (*par_unpack_result)() = NULL;	/* optional */

/* the internal functions */
static void ScanEnvironment();
static void DispatchTasks();
static void DispatchTask();
static int  FindReadyWorker();
static void HandleMessage();
static void HandleRequest();
static void HandleBroadcastAck();
static void HandleWorkerExit();
static void QueueTask();
static void RequeueTask();
static Context *ReuseContext();
static void DisuseContext();
static void FreeBuffer ();
static void PutOnReadyList();
static void RemoveFromReadyList();
static void PutOnIdleList();
static void RemoveFromIdleList();
static void DeclareWorkerDead();
static void PerformWorkerTasks();
static void ComposeRequest();
static int  MasterReceiveMessage(float timeout);
static int  WorkerReceiveMessage();
static void PrepareToSend();
static void Send();
static void GetHostNameFromTID();

#if defined(PVM)
static void StartPVM();
static void HandleHostAdded();
static void HandleHostDeleted();
static void HandleRestart();
static void GetOfficialHostname();
static void AddHosts();
static void StartWorkers();

#elif defined(MPI)
static void StartMPI();
static void ExpandInBuffer(int size);
static void ExpandOutBuffer(int size);
#endif

void
par_process (int argc,
	     char **argv,
	     char **envp,
	     void (*master_task)(),
	     void (*master_result)(),
	     void (*worker_context)(),
	     void (*worker_task)(),
	     void (*worker_finalize)(),
	     void (*pack_context)(),
	     void (*unpack_context)(),
	     void (*pack_task)(),
	     void (*unpack_task)(),
	     void (*pack_result)(),
	     void (*unpack_result)())
{
  struct hostent *e;

  /* save the user-supplied functions */
  if (master_task == NULL)
    Abort("libpar: master_task cannot be NULL\n");
  par_master_task = master_task;
  par_master_result = master_result;
  par_worker_context = worker_context;
  if (worker_task == NULL)
    Abort("libpar: worker_task cannot be NULL\n");
  par_worker_task = worker_task;
  par_worker_finalize = worker_finalize;
  par_pack_context = pack_context;
  par_unpack_context = unpack_context;
  par_pack_task = pack_task;
  par_unpack_task = unpack_task;
  par_pack_result = pack_result;
  par_unpack_result = unpack_result;

  if (gethostname(my_host_name, MAX_HOSTNAME_LENGTH) < 0)
    Abort("gethostname failed\n");
  e = gethostbyname(my_host_name);
  if (e != NULL)
    StringCopy(my_host_name, e->h_name, MAX_HOSTNAME_LENGTH);
  else
    strcpy(my_host_name, "unknown_host");

#if defined(PVM)      /* this has been compiled for PVM */
  spawn = TRUE;
  closed_worker_set = FALSE;

  /* check the relevant environment variables */
  ScanEnvironment(argc, argv);

  if (par)
    /* and the parallelism flag is on, so start running in a PVM mode */
    StartPVM(argc, argv, envp);

  if (!par) /* may just have been turned off by StartPVM */
    /* and the parallelism flag is off, so just run the master task */
    (*master_task)(prog_argc, prog_argv, envp);

#elif defined(MPI)   /* this has been compiled for MPI */
  spawn = FALSE;
  closed_worker_set = FALSE;

  /* we have to enroll with MPI here in order for the command-line
     arguments to be valid */

  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    Abort("Could not initialize MPI\n");
  /* check the relevant environment variables */
  ScanEnvironment(argc, argv);

  if (spawn)
    {
      Error("Sorry, spawning is not allowed with MPI-based libpar\n");
      spawn = FALSE;
    }

  if (par) 
    {
      /* the parallelism flag is on, so start running in an MPI mode */
      StartMPI(argc, argv, envp);
    }

  if (!par) /* par may possibly have been turned off by StartMPI */
    {
      /* the parallelism flag is off, so just run the master task */
      (*master_task)(prog_argc, prog_argv, envp);
      if (worker_finalize != NULL)
	(*worker_finalize)();
    }

  if (par_verbose) Report("Tid %d at MPI_Finalize\n",my_tid);
  if (MPI_Finalize() != MPI_SUCCESS) {
    Abort("MPI_Finalize failed\n");
  }

#else		/* this hasn't been compiled for PVM or MPI */
  spawn = FALSE;
  closed_worker_set = TRUE;

  /* check the relevant environment variables */
  ScanEnvironment(argc, argv);

  if (par)
    {
      Error("Sorry, this has not been compiled to use PVM or MPI\n");
      /* turn parallelism off */
      par = FALSE;
    }

  /* just run the master task */
  (*master_task)(prog_argc, prog_argv, envp);
  if (worker_finalize != NULL)
    (*worker_finalize)();
#endif
}

void par_set_hosts (char *s)
{
#ifdef PVM
  char *t;
  char *p;
  Hostname hn;
  char *tmp;
  char *ptmp;

  /* make a local copy of s because strtok is destructive */
  tmp = (char *) malloc((strlen(s)+1)*sizeof(char));
  strcpy(tmp, s);
  ptmp = tmp;

  n_hosts = 0;
  while ((t = strtok(ptmp, ", \t\n")) != NULL)
    {
      StringCopy(hn, t, MAX_HOSTNAME_LENGTH);
      hosts[n_hosts].max_workers = 1;
      hosts[n_hosts].n_workers = 0;
      hosts[n_hosts].dtid = -1;
      if ((p = strchr(hn, ':')) != NULL)
	{
	  *p = '\0';
	  GetOfficialHostname(hosts[n_hosts].name, hn);
	  sscanf(p+1, "%d", &hosts[n_hosts].max_workers);
	  *p = ':';
	}
      else
	GetOfficialHostname(hosts[n_hosts].name, hn);
      if (hosts[n_hosts].max_workers > 0)
	++n_hosts_pending;
      ++n_hosts;
      ptmp = NULL;
    }
  free(tmp);

#else
  Warning(1,"par_set_hosts() has no effect on non-PVM-based libpar\n");
#endif
}

static void
ScanEnvironment (int argc, char **argv)
{
  char *p;
  int i;
  int v;
  char s[256];

  if ((p = strrchr(argv[0], '/')) == NULL)
    StringCopy(prog_name, argv[0], 128);
  else
    StringCopy(prog_name, p+1, 128);

  if ((p = getenv("PAR_ENABLE")) != NULL &&
      sscanf(p, "%d", &v) == 1)
    par = v != 0;
  if ((p = getenv("PAR_DEBUG")) != NULL &&
      sscanf(p, "%d", &v) == 1)
    debug_workers = v != 0;
  if ((p = getenv("PAR_GROUP")) != NULL)
    StringCopy(group_name, p, sizeof(group_name) - 1);
  else
    group_name[0] = '\0';
  if ((p = getenv("PAR_NOSPAWN")) != NULL &&
      sscanf(p, "%d", &v) == 1)
    spawn = (v == 0);
  if ((p = getenv("PAR_CWD")) != NULL)
    if (chdir(p) != 0)
      Error("Could not change to working directory %s\n", p);
  if ((p = getenv("PAR_HOSTS")) != NULL)
    par_set_hosts(p);
  if ((p = getenv("PAR_VERBOSE")) != NULL &&
      sscanf(p, "%d", &v) == 1)
    par_verbose = (v != 0);

  prog_argc = 0;
  prog_argv = (char **) malloc(argc * sizeof(char *));
  par_argc = 0;
  par_argv = (char **) malloc(argc * sizeof(char *));
  for (i = 0; i < argc; ++i)
    {
      par_argv[par_argc++] = argv[i];
      if (strncmp(argv[i], "-PAR_", 5) == 0) {
	if (strncmp(argv[i], "-PAR_ENABLE=", 12) == 0)
	  {
	    if (sscanf(&argv[i][12], "%d", &v) == 1)
	      par = v != 0;
	  }
	else if (strncmp(argv[i], "-PAR_DEBUG=", 11) == 0)
	  {
	    if (sscanf(&argv[i][11], "%d", &v) == 1)
	      debug_workers = v != 0;
	  }
	else if (strncmp(argv[i], "-PAR_GROUP=", 11) == 0)
	  {
	    if (sscanf(&argv[i][11], "%s", group_name) != 1)
	      group_name[0] = '\0';
	  }
	else if (strncmp(argv[i], "-PAR_NOSPAWN=", 13) == 0)
	  {
	    if (sscanf(&argv[i][13], "%d", &v) == 1)
	      spawn = (v == 0);
	  }
	else if (strncmp(argv[i], "-PAR_CWD=", 9) == 0)
	  {
	    if (sscanf(&argv[i][9], "%s", s) == 1)
	      {
		if (chdir(s) != 0)
		  Error("Could not change to working directory %s\n", s);
	      }
	  }
	else if (strncmp(argv[i], "-PAR_HOSTS=", 11) == 0)
	  {
	    if (sscanf(&argv[i][11], "%s", s) == 1)
	      par_set_hosts(s);
	  }
	else if (strncmp(argv[i], "-PAR_VERBOSE=", 13) == 0)
	  {
	    if (sscanf(&argv[i][13], "%d", &v) == 1)
	      par_verbose = (v != 0);
	  }
      }
      else {
	prog_argv[prog_argc++]= argv[i];
      }
    }
  if (par_verbose) verbose= 1; /* set global flag in libmisc */
}

Par_Task
par_delegate_task ()
{
  int i;
  int new_task_completed_size;
  int index, bit;
  Task* task;

  /* if we are not running in parallel, do the worker task ourself */
  if (!par)
    {
      if (par_verbose)
	Report("Performing task %d myself\n", task_number);
      (*par_worker_task)();
      if (par_master_result != NULL)
	(*par_master_result)(task_number);
      return(task_number++);
    }

  if (task_number - task_completed_first >= 8*task_completed_size)
    {
      /* double the size of the task_completed table */
      if (task_completed_size == 0)
	{
	  new_task_completed_size = 1024;
	  task_completed = (unsigned char *) malloc(new_task_completed_size *
						    sizeof(unsigned char));
        }
      else
	{
	  new_task_completed_size = 2*task_completed_size;
	  task_completed = realloc(task_completed,
				   new_task_completed_size *
				   sizeof(unsigned char));
	}
      for (i = task_completed_offset; i < task_completed_size; ++i)
	task_completed[task_completed_size + i] = task_completed[i];
      task_completed_offset += task_completed_size;
      task_completed_size = new_task_completed_size;
    }
  index = (((task_number - task_completed_first) >> 3) + task_completed_offset)
    % task_completed_size;
  bit = (task_number - task_completed_first) & 7;
  task_completed[index] &= ~(0xff << bit);

  ++tasks_outstanding;
  if (context_changed)
    {
      current_context = (Context *) malloc(sizeof(Context));
      fprintf(stderr, "Context allocated at %lx\n", (uintptr_t)current_context);
      current_context->use_count = 1;
#if defined(PVM)
      current_context->buffer.bufid = pvm_mkbuf(PvmDataDefault);
      if (pvm_setsbuf(current_context->buffer.bufid) < 0)
	Abort("Could not set active buffer.\n");
      if (par_pack_context != NULL)
	(*par_pack_context)();
      if (pvm_setsbuf(0) < 0)
	Abort("Could not reset active buffer.\n");
#elif defined(MPI)
      out_position = 0;
      if (par_pack_context != NULL)
	(*par_pack_context)();
      current_context->buffer.buffer = out_buffer;
      current_context->buffer.size = out_size;
      current_context->buffer.position = out_position;
      out_buffer = NULL;
      out_size = 0;
      out_position = 0;
#endif
      context_changed = FALSE;
    }
  task = (Task*) malloc(sizeof(Task));
  task->next = NULL;
  task->prev = NULL;
  task->context = ReuseContext(current_context);
  task->number = task_number;

#if defined(PVM)
  task->buffer.bufid = pvm_mkbuf(PvmDataDefault);
  if (pvm_setsbuf(task->buffer.bufid) < 0)
    Abort("Could not set active buffer.\n");
  if (pvm_pkint(&task_number, 1, 1) < 0)
    Abort("Could not pack task_number\n");
  if (par_pack_task != NULL)
    (*par_pack_task)();
  if (pvm_setsbuf(0) < 0)
    Abort("Could not reset active buffer.\n");

#elif defined(MPI)
  out_position = 0;
  par_pkint(task_number);
  if (par_pack_task != NULL)
    (*par_pack_task)();
  task->buffer.buffer = out_buffer;
  task->buffer.size = out_size;
  task->buffer.position = out_position;
  out_buffer = NULL;
  out_size = 0;
  out_position = 0;

#endif

  QueueTask(task);
  DispatchTasks();
  return(task_number++);
}

void
par_set_context ()
{
  /* if we are not running in parallel, call WorkerContext
     directly */
  if (!par)
    {
      if (par_worker_context != NULL)
	(*par_worker_context)();
      return;
    }

  if (current_context != NULL)
    DisuseContext(&current_context);
  context_changed = TRUE;
}

void
par_broadcast_context ()
{
  int total_time;
  int bufid;
  int i;
  int oldid;
  struct timeval current_time;

  /* if we are not running in parallel, call WorkerContext directly */
  if (!par)
    {
      if (par_worker_context != NULL)
	(*par_worker_context)();
      return;
    }

  if (par_verbose)
    Report("Master going to broadcast context %d\n", broadcast_count);

  /* wait for all tasks to complete */
  while (tasks_outstanding > 0)
    (void) MasterReceiveMessage(PAR_FOREVER);

  /* wait for any pending workers to start up */
#if defined(T3D) || defined(T3E)
  while (n_workers_pending > 0)
    MasterReceiveMessage(PAR_FOREVER);
#else
  total_time = 0;
  while ((n_hosts_pending > 0 || n_workers_pending > 0) &&
	 total_time < PENDING_WORKERS_TIMEOUT)
    {
      if (gettimeofday(&current_time, NULL) != 0)
	Abort("Master could not get time of day.\n");
      printf("pbc: %d %d %d %ld %ld\n",
	     n_hosts_pending, n_workers_pending,
	     total_time,
	     current_time.tv_sec, current_time.tv_usec);
      if (MasterReceiveMessage(REPORT_INTERVAL) == 0)
	total_time += REPORT_INTERVAL;
    }
  if (gettimeofday(&current_time, NULL) != 0)
    Abort("Master could not get time of day.\n");
  printf("pbc_cls: %d %d %d %ld %ld\n",
	 n_hosts_pending, n_workers_pending,
	 total_time,
	 current_time.tv_sec, current_time.tv_usec);
  fflush(stdout);
  /* give up on any that have not responded */
  n_hosts_pending = 0;
  n_workers_pending = 0;
#endif

  /* doing a broadcast effectively closes the set of workers since any
     new workers could not be initialized with the proper context */
  closed_worker_set = TRUE;

  /* if there are many broadcasts already in progress, we have to
     wait here for acknowledgements */
  while (broadcast_count > MIN(broadcast_first_ack_count,
			       broadcast_second_ack_count) +
	                   PAR_MAX_OVERLAPPED_BROADCASTS)
    MasterReceiveMessage(PAR_FOREVER);

  if (broadcast_count == 0)
    {
      if (n_workers <= 0)
	Abort("Broadcast operation not possible without any workers!\n");
      process_tids[0] = master_tid;
      for (i = 0; i < n_workers; ++i)
	process_tids[i+1] = workers[i].tid;
      n_processes = n_workers + 1;
      if (n_processes > 1)
	broadcast_first_forward_tid = process_tids[1];
      if (n_processes > 2)
	broadcast_second_forward_tid = process_tids[2];
    }
  if (par_verbose)
    Report("Master broadcasting context %d\n", broadcast_count);
  if (broadcast_first_forward_tid >= 0)
    {
      /* prepare message */
      PrepareToSend();
      par_pkint(broadcast_count);
      if (broadcast_count == 0)
	{
	  par_pkint(n_processes);
	  par_pkintarray(process_tids, n_processes);
	}
      if (par_pack_context != NULL)
	(*par_pack_context)();
      Send(broadcast_first_forward_tid, BROADCAST_CONTEXT_MSG);
    }
  else 
    broadcast_first_ack_count = broadcast_count;

  if (broadcast_second_forward_tid >= 0)
    {
      /* prepare message */
      PrepareToSend();
      par_pkint(broadcast_count);
      if (broadcast_count == 0)
	{
	  par_pkint(n_processes);
	  par_pkintarray(process_tids, n_processes);
	}
      if (par_pack_context != NULL)
	(*par_pack_context)();
      Send(broadcast_second_forward_tid, BROADCAST_CONTEXT_MSG);
    }
  else 
    broadcast_second_ack_count = broadcast_count;
  ++broadcast_count;
}

int
par_wait (float timeout)
{
#if defined(PVM) || defined(MPI)
  /* if we are not running in parallel, there is nothing to do */
  if (!par)
    return(0);
  return(MasterReceiveMessage(timeout));
#else
  return(0);
#endif
}

void
par_finish ()
{
  int i;
  int total_time;

  /* if we are not running in parallel, there is nothing to do */
  if (!par)
    return;

  while (tasks_outstanding > 0)
    (void) MasterReceiveMessage(PAR_FOREVER);

#if defined(T3D) || defined(T3E) || defined(MPI)
  /* on the T3D/E we wait for all workers to start up before shutting down */
  while (n_workers_pending > 0)
    MasterReceiveMessage(PAR_FOREVER);
#else
  total_time = 0;
  while (n_workers_pending > 0)
    if (MasterReceiveMessage(REPORT_INTERVAL) == 0)
      {
	total_time += REPORT_INTERVAL;
	if (total_time >= PENDING_WORKERS_TIMEOUT)
	  {
	    if (REPORT_INTERVAL < PENDING_WORKERS_TIMEOUT)
	      Error("\n");
	    break;
	  }
	if (total_time == REPORT_INTERVAL)
	  Error("Finished all tasks but still waiting for pending workers");
	Error(".");
	fflush(stdout);
      }
#endif

#if defined(PVM)
  if (pvm_initsend(PvmDataDefault) < 0)
    Abort("Master: pvm_initsend failed.\n");
  for (i = 0; i < n_workers; ++i)
    if (pvm_send(workers[i].tid, TERMINATE_MSG) < 0)
      DeclareWorkerDead(i--);  /* we need the -- here because
				 DeclareWorkerDead will remove
				 element i from the workers
				 array, and so we want to rerun
				 the loop with the same index i */
#elif defined(MPI)
  if (par_verbose)
    Report("Sending terminate messages to %d workers\n",n_workers);
  for (i = 0; i < n_workers; ++i)
    if (MPI_Send(out_buffer, 0, MPI_PACKED,
		   workers[i].tid, TERMINATE_MSG,
		   MPI_COMM_WORLD) != MPI_SUCCESS)
      Abort("Master cannot send TERMINATE message to worker.\n");
#endif

  n_workers = 0;   /* consider all workers dead */
}

int
par_finished (Par_Task id)
{
  int index, bit;

  if (!par)
    return(id < task_number);

  if (id < task_completed_first)
    return(TRUE);
  if (id >= task_number)
    return(FALSE);
  index = (((id - task_completed_first) >> 3) + task_completed_offset) %
    task_completed_size;
  bit = (id - task_completed_first) & 7;
  return((task_completed[index] >> bit) & 1);
}

int
par_finished_up_to (Par_Task id)
{
  int index, bit, test_bits;

  if (!par)
    return(id < task_number);

  if (id < task_completed_first)
    return(TRUE);
  if (id >= (task_completed_first+8) ||
      id >= task_number)
    return(FALSE);
  index = (((id - task_completed_first) >> 3) + task_completed_offset) %
    task_completed_size;
  bit = (id - task_completed_first) & 7;
  test_bits = 0xff >> (7-bit);
  return((task_completed[index] & test_bits) == test_bits);
}

int
par_tasks_outstanding ()
{
  if (!par)
    return(0);
  return(tasks_outstanding);
}

Par_Task par_next_task()
{
  return(task_number);
}

int
par_enabled ()
{
  return(par);
}

int
par_workers ()
{
  if (!par)
    return(1);
  return(n_workers);
}

int
par_instance ()
{
  if (!par)
    return(0);
  return(rank);
}

void
par_restart ()
{
  if (par_verbose)
    Report("Worker requesting restart\n");
#ifdef PVM
  PrepareToSend();
  Send(master_tid, RESTART_MSG);
#else
  Abort("par_restart() is only supported in PVM-based versions of libpar.\n");
#endif
}

const char* 
par_arch ()
{
  static char* arch= NULL;

#if defined(PVM)
  if (!arch)
    {
      char* s = getenv("PVM_ARCH");
      if (s)
	arch = strdup(s);
    }
#endif

  if (!arch)
    {
      FILE *f = popen("uname -s", "r");
      if (f != NULL)
	{
	  char s[1024];
	  if (fscanf(f, "%s", s) == 1)
	    arch = strdup(s);
	  else
	    arch = "DUMMY";
	  pclose(f);
	}
      else
	arch = "DUMMY";
    }

  return arch;
}

/* DispatchTasks does not return until all queued tasks have been assigned to
   workers for processing */
static void
DispatchTasks ()
{
  /* make sure all broadcasts have completed */
  while (broadcast_first_ack_count < (broadcast_count - 1) ||
	 broadcast_second_ack_count < (broadcast_count - 1))
    MasterReceiveMessage(PAR_FOREVER);

  /* dispatch all tasks */
  while (first_queued_task != NULL)
      DispatchTask();
}

static void
DispatchTask ()
{
  int n;
  int res;
  Task *task;

  /* take the first task off the queue */
  task = first_queued_task;
  first_queued_task = task->next;
  if (first_queued_task != NULL)
    first_queued_task->prev = NULL;
  else
    last_queued_task = NULL;
  task->next = NULL;

  /* find a worker to run it on */
  n = FindReadyWorker();

  /* check if the context needs to be sent */
  if (task->context != NULL &&
      task->context != workers[n].last_context)
    {
#if defined(PVM)
      if (pvm_setsbuf(task->context->buffer.bufid) < 0)
	Abort("Master cannot switch to context send buffer.\n");
      res = pvm_send(workers[n].tid, CONTEXT_MSG);
      if (res < 0)
	if (res == PvmBadParam || res == PvmSysErr)
	  {
	    DeclareWorkerDead(n);
	    return;
	  }
	else Abort("Error sending context to worker.\n");

#elif defined(MPI)
      if (MPI_Send(task->context->buffer.buffer, task->context->buffer.position,
		   MPI_PACKED, workers[n].tid, CONTEXT_MSG,
		   MPI_COMM_WORLD) != MPI_SUCCESS)
	Abort("Error sending context to worker.\n");
#endif

      DisuseContext(&workers[n].last_context);
      workers[n].last_context = ReuseContext(task->context);
    }

  /* send the task */
  if (par_verbose)
    Report("Master sending task %d to worker %d on %s\n",
	   task->number, n, workers[n].host);

#if defined(PVM)
  if (pvm_setsbuf(task->buffer.bufid) < 0)
    Abort("Master cannot switch to task send buffer.\n");
  res = pvm_send(workers[n].tid, TASK_MSG);
  if (res < 0)
    if (res == PvmBadParam || res == PvmSysErr)
      {
	DeclareWorkerDead(n);
	return;
      }
    else Abort("Error sending task to worker.\n");

#elif defined(MPI)
  if (MPI_Send(task->buffer.buffer, task->buffer.position, MPI_PACKED,
	       workers[n].tid, TASK_MSG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Error sending task to worker.\n");
#endif

  /* record that worker n is performing the task */
  if (workers[n].first_task == NULL)
    {
      workers[n].first_task = task;
      RemoveFromIdleList(n);
    }
  else if (workers[n].second_task == NULL)
    {
      workers[n].second_task = task;
      RemoveFromReadyList(n);
    }
  else Abort("Worker %d requested more than 2 tasks at a time.\n", n);
}

static int
FindReadyWorker ()
{
  for (;;)
    {
      /* if a worker is totally idle, choose it */
      if (idle_workers >= 0)
	return(idle_workers);

      /* if no more workers are expected soon, then
	 any ready worker will do */
      if (n_workers_pending <= 0 && ready_workers >= 0)
	return(ready_workers);

      /* wait for something to happen */
      (void) MasterReceiveMessage(PAR_FOREVER);
    }
}

static int
MasterReceiveMessage (float timeout)
{
#if defined(PVM)
  int bufid;
  int current_time;
  struct timeval tmout;
  int n_bytes;
  int msg_tag;
  int tid;

#ifdef ACCT
  Acct(MSGWAITING);
#endif

  if (!spawn || closed_worker_set)
    {
      if (timeout == 0.0)
	{
	  /* check if any message has arrived */
	  bufid = pvm_probe(-1, -1);
	  if (bufid < 0)
	    Abort("Master cannot receive messages.\n");
	  if (bufid == 0)
	    /* no message has arrived */
	    return(0);
	  /* a message is there; now receive it */
	  bufid = pvm_recv(-1, -1);
	}
      else if (timeout < 0.0)
	{
	  /* we can't spawn new workers, so we'll have
	     to be patient and wait for one to free up */
	  bufid = pvm_recv(-1, -1);
	}
      else
      {
	/* wait until a message arrives or the timeout has expired */
	tmout.tv_sec = floor(timeout);
	tmout.tv_usec = floor(1000000.0*(timeout - tmout.tv_sec));
	bufid = pvm_trecv(-1, -1, &tmout);
      }
    }
  else
    {
      /* check if it's time to try to add new hosts */
      current_time = time(0);
      if (current_time > add_host_time)
	{
	  /* give up on any workers that haven't started yet */
	  n_workers_pending = 0;
	  /* add the new hosts */
	  AddHosts();
	  /* reset timer */
	  add_host_time = current_time + RETRY_HOST_INTERVAL;
	}

      if (timeout == 0.0)
	{
	  /* check if any message has arrived */
	  bufid = pvm_probe(-1, -1);
	  if (bufid < 0)
	    Abort("Master cannot receive messages.\n");
	  if (bufid == 0)
	    /* no message has arrived */
	    return(0);
	  /* a message is there; now receive it */
	  bufid = pvm_recv(-1, -1);
	}
      else if (timeout < 0.0)
	{
	  /* wait until a message arrives or it is time to
	     try adding a new host */
	  tmout.tv_sec = add_host_time - current_time;
	  tmout.tv_usec = 0;
	  bufid = pvm_trecv(-1, -1, &tmout);
	}
      else
	{
	  /* wait until a message arrives or it is time to
	     try adding a new host or the timeout expires */
	  tmout.tv_sec = floor(timeout);
	  tmout.tv_usec = floor(1000000.0*(timeout - tmout.tv_sec));
	  if ((add_host_time - current_time) <= tmout.tv_sec)
	    {
	      tmout.tv_sec = add_host_time - current_time;
	      tmout.tv_usec = 0;
	    }
	  bufid = pvm_trecv(-1, -1, &tmout);
	}
    }
#ifdef ACCT
  Acct(PROCESSING);
#endif

  if (bufid < 0)
    Abort("Master cannot receive messages.\n");
  if (bufid > 0)
    {
      if (pvm_bufinfo(bufid, &n_bytes, &msg_tag, &tid) < 0)
	Abort("Master cannot get buffer information.\n");
      HandleMessage(msg_tag, tid);
      return(1);
    }
  return(0);

#elif defined(MPI)
  int flag;
  int len;
  MPI_Status status;
  struct timeval end_time, current_time;
  struct timespec duration;

#ifdef ACCT
  Acct(MSGWAITING);
#endif

  if (timeout == 0.0)
    {
      /* check if any message has arrived */
      if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		     &flag, &status) != MPI_SUCCESS)
	Abort("Could not probe for messages in MasterReceiveMessage()\n");

      if (!flag)
	/* no message has arrived */
	return(0);

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");
      if (len > in_size)
	ExpandInBuffer(len);

      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }
  else if (timeout < 0.0)
    {
      /* we can't spawn new workers, so we'll have
	 to be patient and wait for one to free up */
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		    &status) != MPI_SUCCESS)
	Abort("Could not probe for messages in MasterReceiveMessage()\n");

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");
      if (len > in_size)
	ExpandInBuffer(len);

      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }
  else
    {
      /* wait until a message arrives or the timeout has expired */
      /* since MPI does not support timeouts, and we don't want to use
	 multiple threads, we have to poll (yuck) */
      end_time.tv_sec = 0;
      if (gettimeofday(&end_time, NULL) != 0)
	Abort("Master could not get time of day.\n");
      end_time.tv_sec += floor(timeout);
      end_time.tv_usec += floor(1000000.0*(timeout - floor(timeout)));
      if (end_time.tv_usec >= 1000000)
	{
	  end_time.tv_sec++;
	  end_time.tv_usec -= 1000000;
	}
      duration.tv_sec = 0;
      duration.tv_nsec = 10000000;	/* 10 milliseconds */
      do {
	/* check if any message has arrived */
	if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		       &flag, &status) != MPI_SUCCESS)
	  Abort("Could not probe for messages in MasterReceiveMessage()\n");

	if (flag)
	  break;

	nanosleep(&duration, NULL);

	if (gettimeofday(&current_time, NULL) != 0)
	  Abort("Master could not get time of day.\n");
      } while (current_time.tv_sec < end_time.tv_sec ||
	       (current_time.tv_sec == end_time.tv_sec &&
		current_time.tv_usec < end_time.tv_usec));

      if (!flag)
	return(0);

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");
      if (len > in_size)
	ExpandInBuffer(len);
      
      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }

#ifdef ACCT
  Acct(PROCESSING);
#endif

  in_position = 0;
  HandleMessage(status.MPI_TAG, status.MPI_SOURCE);
  return(1);

#else
  return(0);

#endif
}

static void
HandleMessage (int msg_tag, int tid)
{
  switch (msg_tag)
    {
    case REQUEST_MSG:
      HandleRequest(tid);
      break;
	      
    case BROADCAST_ACK_MSG:
      HandleBroadcastAck(tid);
      break;

    case WORKER_EXIT_MSG:
      HandleWorkerExit();
      break;

#if defined(PVM)
    case RESTART_MSG:
      HandleRestart(tid);
      break;

    case HOST_ADDED_MSG:
      HandleHostAdded();
      break;

    case HOST_DELETED_MSG:
      HandleHostDeleted();
      break;
#endif

    default:
      Error("Master received unknown message: %d\n", msg_tag);
    }
}

static void
HandleRequest (int tid)
{
  int tc;
  int n;
  int result_appended;
  int index, bit;
  int current;
  Hostname worker_host_name;

  /* locate the worker in the worker table */
  tc = par_upkint();
  if (tc < 0)
    par_upkstr(worker_host_name);
  for (n = 0; n < n_workers; ++n)
    if (workers[n].tid == tid)
      break;	/* found it */
  if (n >= n_workers)
    {
      /* this is a new worker */
      if (closed_worker_set)
	return;	/* sorry -- too late */
      --n_workers_pending;
      /* check that integer value sent is -1 */
      if (tc != -1)
	Abort("Initial worker message was invalid.\n");
#if defined(PVM)
      if (spawn)
	if (pvm_notify(PvmTaskExit, WORKER_EXIT_MSG, 1, &tid) < 0)
	  return;	/* ignore this -- it may be a worker that
			   got killed almost immediately */
#endif
      /* add it to the table of workers */
      workers[n_workers].tid = tid;
      strcpy(workers[n_workers].host, worker_host_name);
      workers[n_workers].ready = FALSE;
      workers[n_workers].idle = FALSE;
      workers[n_workers].last_context = NULL;
      workers[n_workers].first_task = NULL;
      workers[n_workers].second_task = NULL;
      PutOnReadyList(n_workers);
      PutOnIdleList(n_workers);
      if (par_verbose)
	Report("Worker %d started on host %s\n", n_workers, workers[n_workers].host);
      ++n_workers;
      return;
    }

  /* this is an old worker */
  /* update the first_task and second_task fields */
  if (tc >= 0 && workers[n].first_task->number == tc)
    {
      /* we are done with the first task */
      if (par_verbose)
	Report("Master received result of task %d from worker %d on host %s\n",
	       tc, n, workers[n].host);
      result_appended = par_upkint();
      if (result_appended && par_master_result != NULL)
	{
	  if (par_unpack_result != NULL)
	    (*par_unpack_result)();
	  (*par_master_result)(tc);
	}
      --tasks_outstanding;
      DisuseContext(&workers[n].first_task->context);

      /* free up the task buffer */
      FreeBuffer(&workers[n].first_task->buffer);

      /* mark as finished */
      index = (((workers[n].first_task->number - task_completed_first) >> 3) +
	       task_completed_offset) % task_completed_size;
      bit = (workers[n].first_task->number - task_completed_first) & 7;
      task_completed[index] |= 1 << bit;
      if (index == task_completed_offset)
	{
	  current = (((task_number - task_completed_first) >> 3) +
		     task_completed_offset) % task_completed_size;
	  while (task_completed_offset != current)
	    {
	      if (task_completed[task_completed_offset] != 0xff)
		break;
	      task_completed_first += 8;
	      if (++task_completed_offset >= task_completed_size)
		task_completed_offset = 0;
	    }
	}

      free(workers[n].first_task);
      workers[n].first_task = NULL;

      /* the second task is now the first */
      if (workers[n].second_task != NULL)
	{
	  workers[n].first_task = workers[n].second_task;
	  workers[n].second_task = NULL;
	}
    }
  else
    if (tc != -1)
      Error("Worker %d on host %s completed task %d instead of %d\n",
	    n, workers[n].host, tc, workers[n].first_task->number);

  if (workers[n].first_task == NULL)
    PutOnIdleList(n);
  if (workers[n].second_task == NULL)
    PutOnReadyList(n);
}

static void
HandleBroadcastAck (int from_tid)
{
  int count;

  count = par_upkint();
  if (par_verbose)
    Report("Master received BROADCAST_ACK to %d from %d\n",
	   count, from_tid);
  if (from_tid == broadcast_first_forward_tid)
    broadcast_first_ack_count = count;
  else if (from_tid == broadcast_second_forward_tid)
    broadcast_second_ack_count = count;
  else
    Abort("Received BROADCAST_ACK_MESSAGE from unexpected source.\n");
}


static void
HandleWorkerExit ()
{
  int stid;
  int n;

  stid = par_upkint();
  for (n = 0; n < n_workers; ++n)
    if (workers[n].tid == stid)
      {
	Error("Worker %d on host %s exited.\n", n, workers[n].host);
	DeclareWorkerDead(n);
	return;
      }
  /* if we get here, we didn't have a record of the worker anyway,
     so ignore the message */
  if (par_verbose)
    Report ("Unknown Worker tid %d.\n", stid);
}

static void
QueueTask (Task* task)
{
  /* put the task at the end of the queued task list so that
     it will wait its turn before getting dispatched */
  task->next = NULL;
  task->prev = last_queued_task;
  if (last_queued_task != NULL)
    last_queued_task->next = task;
  else
    first_queued_task = task;
  last_queued_task = task;
}
     
static void
RequeueTask (Task *task)
{
  /* put the task at the front of the queued task list so that
     it will get dispatched as soon as possible */
  task->next = first_queued_task;
  task->prev = NULL;
  if (first_queued_task != NULL)
    first_queued_task->prev = task;
  else
    last_queued_task = task;
  first_queued_task = task;
}
     
static Context*
ReuseContext (Context *context)
{
  if (context != NULL)
    context->use_count++;
  return(context);
}

static void
DisuseContext (Context **pContext)
{
  Context *context;
  int i;
  Task* t;

  context = *pContext;
  *pContext = NULL;
  if (context == NULL || --context->use_count > 0)
    return;
  if (context->use_count < 0)
    fprintf(stderr, "CONTEXT_ERROR: use_count of context %lx is %d\n",
	    (uintptr_t)context, context->use_count);
  FreeBuffer(&context->buffer);
  free(context);
  fprintf(stderr, "Context freed at %lx\n", (uintptr_t)context);

  /* make sure that context is not used anywhere */
  if (current_context == context)
    fprintf(stderr,
	    "CONTEXT ERROR: current_context is a pointer to deleted context %lx\n",
	    (uintptr_t)context);
  for (i = 0; i < n_workers; ++i)
    {
      if (workers[i].last_context == context)
	fprintf(stderr,
		"CONTEXT ERROR: worker %d has a pointer to deleted context %lx\n",
		i, (uintptr_t)context);
      if (workers[i].first_task != NULL &&
	  workers[i].first_task->context == context)
	fprintf(stderr,
		"CONTEXT ERROR: the first_task of worker %d has a pointer to deleted context %lx\n",
		i, (uintptr_t)context);
      if (workers[i].second_task != NULL &&
	  workers[i].second_task->context == context)
	fprintf(stderr,
		"CONTEXT ERROR: the second_task of worker %d has a pointer to deleted context %lx\n",
		i, (uintptr_t)context);
    }
  for (t = first_queued_task; t != NULL; t = t->next)
    if (t->context == context)
      fprintf(stderr,
	      "CONTEXT ERROR: task %lx on task queue has a pointer to deleted context %lx\n",
	      (uintptr_t)t, (uintptr_t)context);
}

static void
FreeBuffer (Buffer *buffer)
{
#if defined(PVM)
  if (pvm_freebuf(buffer->bufid) < 0)
    Abort("Could not free task buffer.\n");
#elif defined(MPI)
  if (buffer->buffer != NULL)
    {
      free(buffer->buffer);
      buffer->buffer = NULL;
    }
#endif
}

static void
PutOnReadyList (int n)
{
  if (workers[n].ready)
    /* worker was already on ready list */
    return;

  /* put on the ready list */
  workers[n].ready = TRUE;
  workers[n].next_ready = ready_workers;
  ready_workers = n;
}

static void
RemoveFromReadyList (int n)
{
  int pr;
  int r;

  if (!workers[n].ready)
    /* worker wasn't on ready list */
    return;

  /* take this worker off the ready list */
  workers[n].ready = FALSE;
  pr = -1;
  r = ready_workers;
  while (r >= 0)
    {
      if (r == n)
	{
	  if (pr >= 0)
	    workers[pr].next_ready = workers[n].next_ready;
	  else
	    ready_workers = workers[n].next_ready;
	  return;
	}
      pr = r;
      r = workers[r].next_ready;
    }
  Abort("Corrupt ready list.\n");
}

static void
PutOnIdleList (int n)
{
  if (workers[n].idle)
    /* worker was already on idle list */
    return;

  /* put on the idle list */
  workers[n].idle = TRUE;
  workers[n].next_idle = idle_workers;
  idle_workers = n;
}

static void
RemoveFromIdleList (int n)
{
  int pi;
  int i;

  if (!workers[n].idle)
    /* worker wasn't on idle list */
    return;

  /* take this worker off the idle list */
  workers[n].idle = FALSE;
  pi = -1;
  i = idle_workers;
  while (i >= 0)
    {
      if (i == n)
	{
	  if (pi >= 0)
	    workers[pi].next_idle = workers[n].next_idle;
	  else
	    idle_workers = workers[n].next_idle;
	  return;
	}
      pi = i;
      i = workers[i].next_idle;
    }
  Abort("Corrupt idle list.\n");
}

static void
DeclareWorkerDead (int n)
{
  if (workers[n].second_task != NULL)
    RequeueTask(workers[n].second_task);
  if (workers[n].first_task != NULL)
    RequeueTask(workers[n].first_task);
  RemoveFromReadyList(n);
  RemoveFromIdleList(n);

  --n_workers;
  if (n != n_workers)
    {
      /* move the last worker into the slot formerly occupied by the
	 dead worker */
      workers[n] = workers[n_workers];
      if (workers[n_workers].idle)
	{
	  RemoveFromIdleList(n_workers);
	  workers[n].idle = FALSE;
	  PutOnIdleList(n);
	}
      if (workers[n_workers].ready)
	{
	  RemoveFromReadyList(n_workers);
	  workers[n].ready = FALSE;
	  PutOnReadyList(n);
	}
    }
}

static void
PerformWorkerTasks ()
{
  int bufid;
  int msg_type;
  int task_num;
  int from_tid;
  int i;
  int dest;
  int oldid;
  struct timeval task_start;
  struct timeval task_end;
  long task_sec;
  long task_usec;

#ifdef PTHREADS
  /* launch another thread that will actually execute the tasks */
#endif

  ComposeRequest(-1, FALSE);
  Send(master_tid, REQUEST_MSG);

  for (;;)
    {
      msg_type = WorkerReceiveMessage(&from_tid);

      switch (msg_type)
	{
	case CONTEXT_MSG:
	  if (par_unpack_context != NULL)
	    (*par_unpack_context)();
	  if (par_worker_context != NULL)
	    (*par_worker_context)();
	  break;
	case TASK_MSG:
	  task_num = par_upkint();
	  if (par_verbose)
	    Report("Worker received task %d\n", task_num);
	  if (par_unpack_task != NULL)
	    (*par_unpack_task)();
	  if (par_worker_task != NULL) {
	    if (gettimeofday(&task_start, NULL) != 0)
	      Abort("Worker could not get time of day.\n");
	    (*par_worker_task)();
	    if (gettimeofday(&task_end, NULL) != 0)
	      Abort("Worker could not get time of day.\n");
	    task_sec= task_end.tv_sec-task_start.tv_sec;
	    task_usec= task_end.tv_usec-task_start.tv_usec;
	    if (task_usec<0) 
	      {
		task_sec -= 1;
		task_usec += 1000000;
	      }
	    if (task_sec > longest_task.tv_sec 
		||(task_sec==longest_task.tv_sec 
		   && task_usec>longest_task.tv_usec))
	      {
		longest_task.tv_sec= task_sec;
		longest_task.tv_usec= task_usec;
	      }
	  }
	  if (par_master_result != NULL)
	    {
	      ComposeRequest(task_num, TRUE);
	      if (par_pack_result != NULL)
		(*par_pack_result)();
	      Send(master_tid, REQUEST_MSG);
	    }
	  else
	    {
	      ComposeRequest(task_num, FALSE);
	      Send(master_tid, REQUEST_MSG);
	    }
	  break;
	case BROADCAST_CONTEXT_MSG:
	  broadcast_count = par_upkint();
	  if (broadcast_count == 0)
	    {
	      n_processes = par_upkint();
	      par_upkintarray(process_tids, n_processes);
	      /* find my_process_index */
	      for (i = 0; i < n_processes; ++i)
		if (process_tids[i] == my_tid)
		  {
		    my_process_index = i;
		    break;
		  }
	      if (i >= n_processes)
		Abort("Could not find my process tid in process array.\n");
	      dest = (my_process_index << 1) + 1;
	      if (dest < n_processes)
		broadcast_first_forward_tid = process_tids[dest];
	      dest = (my_process_index << 1) + 2;
	      if (dest < n_processes)
		broadcast_second_forward_tid = process_tids[dest];
	      dest = (my_process_index - 1) >> 1;
	      broadcast_backward_tid = process_tids[dest];
	    }
	  if (par_verbose)
	    Report("Process %d received broadcast context %d\n",
		   my_process_index, broadcast_count);
	  /* handle the broadcast message locally */
	  if (par_unpack_context != NULL)
	    (*par_unpack_context)();
	  if (par_worker_context != NULL)
	    (*par_worker_context)();

	  if (broadcast_first_forward_tid >= 0)
	    {
	      if (par_verbose)
		Report("Process %d relaying context %d onto tid %d\n",
		       my_process_index, broadcast_count,
		       broadcast_first_forward_tid);

	      /* prepare message */
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      if (broadcast_count == 0)
		{
		  par_pkint(n_processes);
		  par_pkintarray(process_tids, n_processes);
		}
	      if (par_pack_context != NULL)
		(*par_pack_context)();
	      Send(broadcast_first_forward_tid, BROADCAST_CONTEXT_MSG);
	    }
	  else
	    broadcast_first_ack_count = broadcast_count;

	  if (broadcast_second_forward_tid >= 0)
	    {
	      if (par_verbose)
		Report("Process %d relaying context %d onto tid %d\n",
		       my_process_index, broadcast_count,
		       broadcast_second_forward_tid);

	      /* prepare message */
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      if (broadcast_count == 0)
		{
		  par_pkint(n_processes);
		  par_pkintarray(process_tids, n_processes);
		}
	      if (par_pack_context != NULL)
		(*par_pack_context)();
	      Send(broadcast_second_forward_tid, BROADCAST_CONTEXT_MSG);
	    }
	  else 
	    broadcast_second_ack_count = broadcast_count;

	  if (broadcast_first_ack_count == broadcast_count &&
	      broadcast_second_ack_count == broadcast_count)
	    {
	      /* we can acknowledge this immediately */
	      if (par_verbose)
		Report("Process %d sending immediate acknowledge of %d\n",
		       my_process_index, broadcast_count);
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
	    }
	  break;
	case BROADCAST_ACK_MSG:
	  broadcast_count = par_upkint();
	  if (par_verbose)
	    Report("Process %d received acknowledgement of %d from %d\n",
		   my_process_index, broadcast_count, from_tid);
	  if (from_tid == broadcast_first_forward_tid)
	    {
	      broadcast_first_ack_count = broadcast_count;
	      if (broadcast_first_ack_count <= broadcast_second_ack_count)
		{
		  /* relay acknowledgement */
		  if (par_verbose)
		    Report("Process %d relaying acknowledgement of %d\n",
			   my_process_index, broadcast_count);
		  PrepareToSend();
		  par_pkint(broadcast_count);
		  Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
		}
	    }
	  else if (from_tid == broadcast_second_forward_tid)
	    {
	      broadcast_second_ack_count = broadcast_count;
	      if (broadcast_second_ack_count <= broadcast_first_ack_count)
		{
		  /* relay acknowledgement */
		  if (par_verbose)
		    Report("Process %d relaying acknowledgement of %d\n",
			   my_process_index, broadcast_count);
		  PrepareToSend();
		  par_pkint(broadcast_count);
		  Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
		}
	    }
	  else
	    Abort("Received BROADCAST_ACK_MSG from unexpected source.\n");
	  break;
	case TERMINATE_MSG:
	  if (par_verbose)
	    Report("Worker received TERMINATE_MSG: tid=%d\n", my_tid);
	  return;
	default:
	  Abort("Unknown PVM message type received.\n");
	  break;
	}
    }
}

static void
ComposeRequest (int last_task_completed, int result_appended)
{
  if (par_verbose)
    Report("Worker requesting task\n");
  PrepareToSend();
  par_pkint(last_task_completed);
  if (last_task_completed < 0)
    par_pkstr(my_host_name);
  par_pkint(result_appended);
}

static int
WorkerReceiveMessage (int *pfrom_tid)
{
#if defined(PVM)
  int bufid;
  int n_bytes;
  int msg_tag;
  int tid;
  struct timeval tmout;
  int i;

#ifdef ACCT
  Acct(MSGWAITING);
#endif
  tmout.tv_sec = MAX(WORKER_RECEIVE_TIMEOUT,10*(longest_task.tv_sec+1));
  tmout.tv_usec = 0;
  bufid = pvm_trecv(-1, -1, &tmout);

#ifdef ACCT
  Acct(PROCESSING);
#endif
  if (bufid < 0)
    Abort("Cannot receive task from master.\n");
  if (bufid == 0)
    Abort("Worker timed out waiting for message.\n");
  if (pvm_bufinfo(bufid, &n_bytes, &msg_tag, &tid) < 0)
    Abort("Cannot get buffer information.\n");

  *pfrom_tid = tid;
  return(msg_tag);

#elif defined(MPI)
  int flag;
  int len;
  MPI_Status status;
  struct timeval end_time, current_time;
  struct timespec duration;

  /* wait until a message arrives or the timeout has expired */
  /* since MPI does not support timeouts, and we don't want to use
     multiple threads, we have to poll (yuck) */
  end_time.tv_sec = 0;
  if (gettimeofday(&end_time, NULL) != 0)
    Abort("Worker could not get time of day.\n");
  end_time.tv_sec += MAX(WORKER_RECEIVE_TIMEOUT,10*(longest_task.tv_sec+1));
  duration.tv_sec = 0;
  duration.tv_nsec = 10000000;	/* 10 milliseconds */
  do {
    /* check if any message has arrived */
    if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		   &flag, &status) != MPI_SUCCESS)
      Abort("Could not probe for messages in WorkerReceiveMessage()\n");

    if (flag)
      break;

    nanosleep(&duration, NULL);

    if (gettimeofday(&current_time, NULL) != 0)
      Abort("Worker could not get time of day.\n");
  } while (current_time.tv_sec < end_time.tv_sec ||
	   (current_time.tv_sec == end_time.tv_sec &&
	    current_time.tv_usec < end_time.tv_usec));

  if (!flag)
    {
      Report("Worker timed out waiting for message. Exiting...\n");
      PrepareToSend ();
      par_pkint (my_tid);
      Send (master_tid, WORKER_EXIT_MSG);
      /* Fake msg so finalize will be done before exit.*/
      return (TERMINATE_MSG); 
    }

  /* a message is there; now receive it */
  if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
    Abort("Could not obtain number of bytes in message\n");
  if (len > in_size)
    ExpandInBuffer(len);

  if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &status) != MPI_SUCCESS)
    Abort("Worker could not receive message.\n");

#ifdef ACCT
  Acct(PROCESSING);
#endif

  in_position = 0;

  *pfrom_tid = status.MPI_SOURCE;
  return(status.MPI_TAG);

#else
  return(-1);
#endif
}

static void
PrepareToSend ()
{
#if defined(PVM)
  if (pvm_initsend(PvmDataDefault) < 0)
    Abort("pvm_initsend failed\n");
#elif defined(MPI)
  out_position = 0;
#endif
}

static void
Send (int tid, int tag)
{
#if defined(PVM)
  if (pvm_send(tid, tag) < 0)
    Abort("Cannot send message (tid = %d tag = %d)\n", tid, tag);
#elif defined(MPI)
  if (MPI_Send(out_buffer, out_position, MPI_PACKED,
	       tid, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Cannot send message (tid = %d tag = %d)\n", tid, tag);
#endif
}

void
par_pkbyte(unsigned char v)
{
#if defined(PVM)
  if (pvm_pkbyte((char*)&v, 1, 1) < 0)
    Abort("Could not pack byte into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_byte > out_size)
    ExpandOutBuffer(out_position + sizeof_byte);
  if (MPI_Pack(&v, 1, MPI_BYTE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack byte into MPI buffer\n");
#endif
}

void
par_pkshort(short v)
{
#if defined(PVM)
  if (pvm_pkshort(&v, 1, 1) < 0)
    Abort("Could not pack short into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_short > out_size)
    ExpandOutBuffer(out_position + sizeof_short);
  if (MPI_Pack(&v, 1, MPI_SHORT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack short into MPI buffer\n");
#endif
}

void
par_pkint(int v)
{
#if defined(PVM)
  if (pvm_pkint(&v, 1, 1) < 0)
    Abort("Could not pack int into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_int > out_size)
    ExpandOutBuffer(out_position + sizeof_int);
  if (MPI_Pack(&v, 1, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack int into MPI buffer\n");
#endif
}

void
par_pklong (long v)
{
#if defined(PVM)
  if (pvm_pklong(&v, 1, 1) < 0)
    Abort("Could not pack long into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_long > out_size)
    ExpandOutBuffer(out_position + sizeof_long);
  if (MPI_Pack(&v, 1, MPI_LONG, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack long into MPI buffer\n");
#endif
}

void
par_pkfloat (float v)
{
#if defined(PVM)
  if (pvm_pkfloat(&v, 1, 1) < 0)
    Abort("Could not pack float into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_float > out_size)
    ExpandOutBuffer(out_position + sizeof_float);
  if (MPI_Pack(&v, 1, MPI_FLOAT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack float into MPI buffer\n");
#endif
}

void
par_pkdouble (double v)
{
#if defined(PVM)
  if (pvm_pkdouble(&v, 1, 1) < 0)
    Abort("Could not pack double into PVM buffer\n");
#elif defined(MPI)
  if (out_position + sizeof_double > out_size)
    ExpandOutBuffer(out_position + sizeof_double);
  if (MPI_Pack(&v, 1, MPI_DOUBLE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack double into MPI buffer\n");
#endif
}

void
par_pkstr (char *v)
{
#if defined(PVM)
  if (pvm_pkstr(v) < 0)
    Abort("Could not pack string into PVM buffer\n");
#elif defined(MPI)
  int size;
  int string_len = strlen(v);

  if (MPI_Pack_size(string_len, MPI_CHAR, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of string\n");
  if (out_position + sizeof_int + size > out_size)
    ExpandOutBuffer(out_position + sizeof_int + size);
  if (MPI_Pack(&string_len, 1, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Pack(v, string_len, MPI_CHAR, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack string into MPI buffer\n");
#endif
}

void
par_pkbytearray (unsigned char *p, int n)
{
#if defined(PVM)
  if (pvm_pkbyte((char*)p, n, 1) < 0)
    Abort("Could not pack byte array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_BYTE, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of byte array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_BYTE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack byte array into MPI buffer\n");
#endif
}

void
par_pkshortarray (short *p, int n)
{
#if defined(PVM)
  if (pvm_pkshort(p, n, 1) < 0)
    Abort("Could not pack short array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_SHORT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of short array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_SHORT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack short array into MPI buffer\n");
#endif
}

void
par_pkintarray (int *p, int n)
{
#if defined(PVM)
  if (pvm_pkint(p, n, 1) < 0)
    Abort("Could not pack int array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_INT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of int array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack int array into MPI buffer\n");
#endif
}

void
par_pklongarray (long *p, int n)
{
#if defined(PVM)
  if (pvm_pklong(p, n, 1) < 0)
    Abort("Could not pack long array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_LONG, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of long array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_LONG, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack long array into MPI buffer\n");
#endif
}

void
par_pkfloatarray (float *p, int n)
{
#if defined(PVM)
  if (pvm_pkfloat(p, n, 1) < 0)
    Abort("Could not pack float array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_FLOAT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of float array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_FLOAT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack float array into MPI buffer\n");
#endif
}

void
par_pkdoublearray (double *p, int n)
{
#if defined(PVM)
  if (pvm_pkdouble(p, n, 1) < 0)
    Abort("Could not pack double array into PVM buffer\n");
#elif defined(MPI)
  int size;
  if (MPI_Pack_size(n, MPI_DOUBLE, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of double array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_DOUBLE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack double array into MPI buffer\n");
#endif
}

unsigned char
par_upkbyte ()
{
  unsigned char v;
#if defined(PVM)
  if (pvm_upkbyte((char*)&v, 1, 1) < 0)
    Abort("Could not unpack byte from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_BYTE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack byte from MPI buffer\n");
#endif
  return(v);
}

short
par_upkshort ()
{
  short v;
#if defined(PVM)
  if (pvm_upkshort(&v, 1, 1) < 0)
    Abort("Could not unpack short from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_SHORT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack short from MPI buffer\n");
#endif
  return(v);
}

int
par_upkint ()
{
  int v;
#if defined(PVM)
  if (pvm_upkint(&v, 1, 1) < 0)
    Abort("Could not unpack int from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack int from MPI buffer\n");
#endif
  return(v);
}

long
par_upklong ()
{
  long v;
#if defined(PVM)
  if (pvm_upklong(&v, 1, 1) < 0)
    Abort("Could not unpack long from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_LONG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack long from MPI buffer\n");
#endif
  return(v);
}

float
par_upkfloat ()
{
  float v;
#if defined(PVM)
  if (pvm_upkfloat(&v, 1, 1) < 0)
    Abort("Could not unpack float from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack float from MPI buffer\n");
#endif
  return(v);
}

double
par_upkdouble ()
{
  double v;
#if defined(PVM)
  if (pvm_upkdouble(&v, 1, 1) < 0)
    Abort("Could not unpack double from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_DOUBLE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack double from MPI buffer\n");
#endif
  return(v);
}

void
par_upkstr (char *s)
{
#if defined(PVM)
  if (pvm_upkstr(s) < 0)
    Abort("Could not unpack string from PVM buffer\n");
#elif defined(MPI)
  int string_len;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &string_len, 1, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Unpack(in_buffer, in_size, &in_position,
		 s, string_len, MPI_CHAR, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack string from MPI buffer\n");
  s[string_len] = '\0';
#endif
}

void
par_upkbytearray (unsigned char *p, int n)
{
#if defined(PVM)
  if (pvm_upkbyte((char*)p, n, 1) < 0)
    Abort("Could not unpack byte array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_BYTE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack byte array from MPI buffer\n");
#endif
}

void
par_upkshortarray (short *p, int n)
{
#if defined(PVM)
  if (pvm_upkshort(p, n, 1) < 0)
    Abort("Could not unpack short array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_SHORT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack short array from MPI buffer\n");
#endif
}

void
par_upkintarray (int *p, int n)
{
#if defined(PVM)
  if (pvm_upkint(p, n, 1) < 0)
    Abort("Could not unpack int array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack int array from MPI buffer\n");
#endif
}

void
par_upklongarray (long *p, int n)
{
#if defined(PVM)
  if (pvm_upklong(p, n, 1) < 0)
    Abort("Could not unpack long array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_LONG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack long array from MPI buffer\n");
#endif
}

void
par_upkfloatarray (float *p, int n)
{
#if defined(PVM)
  if (pvm_upkfloat(p, n, 1) < 0)
    Abort("Could not unpack float array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack float array from MPI buffer\n");
#endif
}

void
par_upkdoublearray (double *p, int n)
{
#if defined(PVM)
  if (pvm_upkdouble(p, n, 1) < 0)
    Abort("Could not unpack double array from PVM buffer\n");
#elif defined(MPI)
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_DOUBLE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack double array from MPI buffer\n");
#endif
}


#if defined(PVM)
/*-------------the following functions are all PVM-specific--------------*/

static void
StartPVM (const int argc,
	  const char **argv,
	  const char **envp)
{	
  char *env;
  char * p;

  /* enroll in pvm */
  my_tid = pvm_mytid();
  
  /* Collect output to the master task */
  pvm_catchout(stdout);

#if defined(T3D) || defined(T3E)
  /* on the T3D/E we can get our rank without using a group */
  rank = pvm_get_PE(my_tid);
#else
  /* construct unique group name if none was specified in the environment */
  if (group_name[0] == '\0')
    sprintf(group_name, "%s_%d", prog_name, my_tid);

  /* join pvm group */
  rank = pvm_joingroup(group_name);
#endif

  if (rank == 0)
	{
	  /* I am the master */
	  master_tid = my_tid;

	  /* request notification as new hosts are added to PVM */
	  if (spawn && pvm_notify(PvmHostAdd, HOST_ADDED_MSG, -1, NULL) < 0)
	    Abort("Cannot request notification of added hosts.\n");

	  /* put the group name into the environment */
	  env = (char *) malloc(16+strlen(group_name));
	  sprintf(env, "PAR_GROUP=%s", group_name);
	  putenv(env);

	  /* put the current working directory into the environment */
	  env = (char *) malloc(256*sizeof(char));
	  strcpy(env, "PAR_CWD=");
	  if (getcwd(&env[8], 256-8) == NULL)
	    Abort("Cannot get current working directory!\n");
	  putenv(env);

	  /* pass the following environment variables onto the workers */
	  env = (char *) malloc(128*sizeof(char));
	  if ((p = getenv("PVM_EXPORT")) != NULL)
	    sprintf(env, "PVM_EXPORT=PAR_ENABLE:PAR_GROUP:PAR_NOSPAWN:PAR_DEBUG:PAR_CWD:%s", p);
	  else
	    sprintf(env, "PVM_EXPORT=PAR_ENABLE:PAR_GROUP:PAR_NOSPAWN:PAR_DEBUG:PAR_CWD");
	  putenv(env);

#if defined(T3D) || defined(T3E)
	  n_workers_pending = atoi(getenv("MPP_NPES")) - 1;
#endif
	  (*par_master_task)(prog_argc, prog_argv, envp);
	  /* make sure everything is finished in case par_master_task does not
	     call par_finish itself */
	  par_finish();
	  /*	  pvm_lvgroup(group_name); */
#ifndef ALPHAMP
	  pvm_exit();		/* alphamp prints errors */
#endif
	}
      else
	{
	  /* I am a worker */
	  sleep(1);

#ifdef T3D
	  /* the T3D allows PE numbers to be used in place of tid's */
	  master_tid = 0;
#else
	  master_tid = pvm_gettid(group_name, 0);
#endif
	  PerformWorkerTasks();

	  if (par_worker_finalize != NULL)
	    (*par_worker_finalize)();

	  /*	  pvm_lvgroup(group_name); THIS CALL CAUSES A HANG! */
#ifndef ALPHAMP
	  pvm_exit();		/* alphamp prints errors */
#endif
#ifdef ACCT
	  PrintAcct(argv[0], my_tid);
#endif
	  exit(0);
	}
}

static void
HandleHostAdded ()
{
  int i, j;
  int nhost;
  int narch;
  struct pvmhostinfo *hostp;
  Hostname hname;

  /* ignore this message if we cannot start workers */
  if (!spawn || closed_worker_set)
    return;

  /* forget trying to decode the HOST_ADDED_MSG;
     we will have to call pvm_config anyway to
     get the ASCII host name, so just use the
     new info from pvm_config to update the
     host table */

  (void) pvm_config(&nhost, &narch, &hostp);
  for (i = 0; i < nhost; ++i)
    {
      GetOfficialHostname(hname, hostp[i].hi_name);
      for (j = 0; j < n_hosts; ++j)
	if (strcmp(hosts[j].name, hname) == 0)
	  {
	    /* already in host table */
	    if (hosts[j].dtid < 0)
	      {
		/* host was just added to PVM configuration */
		if (par_verbose)
		  Report("Host %s added to PVM configuration.\n", hname);
		if (pvm_notify(PvmHostDelete, HOST_DELETED_MSG, 1, &hostp[i].hi_tid) < 0)
		  Error("Could not request HOST_DELETED_MSG for %s\n", hname);
		hosts[j].dtid = hostp[i].hi_tid;
		StringCopy(hosts[j].pvm_name, hostp[i].hi_name, MAX_HOSTNAME_LENGTH);
		if (hosts[j].max_workers != 0)
		  StartWorkers(j, abs(hosts[j].max_workers));
		if (hosts[j].max_workers > 0)
		  --n_hosts_pending;
	      }
	    break;
	  }
      if (j >= n_hosts)
	{
	  /* add new entry to host table */
	  if (n_hosts >= PAR_MAX_HOSTS)
	    Abort("Too many hosts added.\n");
	  if (par_verbose)
	    Report("Unlisted host %s added to PVM configuration.\n", hname);
	  if (pvm_notify(PvmHostDelete, HOST_DELETED_MSG, 1, &hostp[i].hi_tid) < 0)
	    Error("Could not request HOST_DELETED_MSG for %s\n", hname);
	  strcpy(hosts[n_hosts].name, hname);
	  hosts[n_hosts].max_workers = -1;
	  hosts[n_hosts].n_workers = 0;
	  hosts[n_hosts].dtid = hostp[i].hi_tid;
	  StringCopy(hosts[n_hosts].pvm_name, hostp[i].hi_name, MAX_HOSTNAME_LENGTH);
	  StartWorkers(n_hosts++, 1);
	}
    }
}

static void
HandleHostDeleted ()
{
  int dtid;
  int n;

  if (pvm_upkint(&dtid, 1, 1) < 0)
    Abort("Invalid HOST_DELETED_MSG\n");
  for (n = 0; n < n_hosts; ++n)
    if (hosts[n].dtid == dtid)
      {
	/* record that host is no longer
	   in the PVM configuration */
	hosts[n].dtid = -1;
	/* this is a recoverable error if there are other hosts left
	   to process the tasks */
	Error("Host %s deleted from PVM configuration.\n", hosts[n].name);
	return;
      }
  /* if we get here, we didn't have a record
     of the host anyway, so ignore the
     message */
}

static void
HandleRestart (int tid)
{
  int n;
  int i;

  for (n = 0; n < n_workers; ++n)
    if (workers[n].tid == tid)
      {
	if (par_verbose)
	  Report("Worker %d on host %s requested a restart.\n", n, workers[n].host);
	DeclareWorkerDead(n);

	if (!spawn)
	  {
	    Error("A worker called par_restart(), but PAR_NOSPAWN is true.\n");
	    return;
	  }
	if (closed_worker_set)
	  {
	    Error("A worker called par_restart(), but the set of workers has been closed\ndue to the use of par_broadcast_context()\n");
	    return;
	  }
	
	for (i = 0; i < n_hosts; ++i)
	  if (strcmp(hosts[i].name, workers[n].host) == 0)
	    break;
	if (i >= n_hosts)
	  {
	    Error("Restart could not find host %s in host table -- ignored\n",
		  workers[n].host);
	    return;
	  }
	StartWorkers(i, 1);
	return;
      }
  Error("Restart requested by unregistered worker -- ignored\n");
}

static void
GetOfficialHostname (char *official_name, char *name)
{
  struct hostent *e;
  char local_name[256];
  char *p_name;

  if (strcmp(name, "localhost") == 0)
    {
      gethostname(local_name, 255);
      p_name = local_name;
    }
  else
    p_name = name;
      
  if ((e = gethostbyname(p_name)) == NULL)
    {
      strcpy(official_name, "???");
      return;
    }
  StringCopy(official_name, e->h_name, MAX_HOSTNAME_LENGTH);
}

static void
AddHosts ()
{
  char *host_names[PAR_MAX_HOSTS];
  int host_index[PAR_MAX_HOSTS];
  int host_info[PAR_MAX_HOSTS];
  int i;
  int nh = 0;
  static Boolean first_time = TRUE;

  /* only try to add the hosts that are not in PVM and
     that we should automatically add to PVM */
  for (i = 0; i < n_hosts; ++i)
    if (hosts[i].dtid < 0 && hosts[i].max_workers > 0)
      {
	host_names[nh] = hosts[i].name;
	host_index[nh] = i;
	++nh;
      }
  if (nh > 0)
    (void) pvm_addhosts(host_names, nh, host_info);
  
  /* check for errors */
  for (i = 0; i < nh; ++i)
    if (host_info[i] < 0)
      switch (host_info[i])
	{
	case PvmDupHost:
	  /* no problem */
	  break;

	case PvmBadParam:
	  Error("Bad hostname syntax: %s\n", host_names[i]);
	  break;

	case PvmNoHost:
	  Error("No such host: %s\n", host_names[i]);
	  break;

	case PvmCantStart:
	  Error("Cannot start PVM on %s\n", host_names[i]);
	  break;

	case PvmBadVersion:
	  Error("Wrong PVM version running on %s\n",
		host_names[i]);
	  break;

	case PvmOutOfRes:
	  Error("PVM out of system resources on %s\n",
		host_names[i]);
	  break;

	default:
	  Error("Unknown addhosts result code (%d) for %s\n",
		host_info[i], host_names[i]);
	  break;
	}

  if (first_time)
    {
      /* the first time AddHosts is called, we fake
	 a HOST_ADDED_MSG from PVM to force an update
	 of the hosts table in this program */
      if (pvm_initsend(PvmDataDefault) < 0)
	Abort("Master: pvm_initsend failed.\n");
      if (pvm_send(my_tid, HOST_ADDED_MSG) < 0)
	Abort("Could not send initial HOST_ADDED_MSG\n");
      first_time = FALSE;
    }
}

static void
StartWorkers (const int host_index, const int count)
{
  char **worker_args;
  char arg1[128];
  int worker_tids[PAR_MAX_WORKERS_PER_HOST];
  int n_workers_added;
  int i;

  /* set up the arguments to use when spawning workers */
  worker_args = (char **) malloc((prog_argc+par_argc) * sizeof(char *));
  for (i = 1; i < prog_argc; ++i)
    worker_args[i-1] = prog_argv[i];
  for (i=0; i<par_argc; ++i)
    worker_args[prog_argc+i-1]= par_argv[i];
  worker_args[prog_argc+par_argc-1] = NULL;

  /* start up workers on that host */
  n_workers_added = pvm_spawn(prog_name, worker_args,
			     (debug_workers ? (PvmTaskHost | PvmTaskDebug) :
			      PvmTaskHost),
			     hosts[host_index].pvm_name,
			     count,
			     worker_tids);
  free(worker_args);
  if (n_workers_added <= 0)
    Error("Could not start workers on %s: reason %d\n", hosts[host_index].name,
	  worker_tids[0]);
  else
    {
      if (n_workers_added < count)
	Error("Could only start %d workers on %s\n",
	      n_workers_added, hosts[host_index].name);
      n_workers_pending += n_workers_added;
    }

  Report("Running %d worker%s on %s\n",
	 n_workers_added,
	 n_workers_added != 1 ? "s" : "",
	 hosts[host_index].name);
}


#elif defined(MPI)
/*-------------the following functions are all MPI-specific--------------*/

static void
StartMPI (int argc,
	  char **argv,
	  char **envp)
{	
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    Abort("Could not obtain rank\n");
  my_tid = rank;

  if (MPI_Pack_size(1, MPI_BYTE, MPI_COMM_WORLD,
		    &sizeof_byte) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_SHORT, MPI_COMM_WORLD,
		    &sizeof_short) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD,
		    &sizeof_int) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_LONG, MPI_COMM_WORLD,
		    &sizeof_long) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_FLOAT, MPI_COMM_WORLD,
		    &sizeof_float) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD,
		    &sizeof_double) != MPI_SUCCESS)
    Abort("Could not determine packing size of elementary types\n");

  if (rank == 0)
    {
      /* I am the master */
      master_tid = 0;
      
      if (MPI_Comm_size(MPI_COMM_WORLD, &n_workers_pending) != MPI_SUCCESS)
	Abort("Cannot obtain number of workers\n");
      
      --n_workers_pending;	/* because the master is not a worker */
      if (n_workers_pending == 0)
	{
	  par = FALSE;
	  Warning(1,
		  "Parallelism disabled since running on only 1 processor\n");
	  return;
	}

      (*par_master_task)(prog_argc, prog_argv, envp);

      /* make sure everything is finished in case par_master_task does not
	 call par_finish itself */
      par_finish();
    }
  else
    {
      /* I am a worker */
      master_tid = 0;

      PerformWorkerTasks();
      if (par_worker_finalize != NULL)
	(*par_worker_finalize)();

#ifdef ACCT
      PrintAcct(argv[0], my_tid);
#endif
    }
}

static void
ExpandInBuffer (int size)
{
  while (in_size < size)
    if (in_size == 0)
      in_size = 1024;
    else
      in_size *= 2;
  in_buffer = realloc(in_buffer, in_size);
  if (in_buffer == NULL)
    Abort("Could not expand in_buffer to %d bytes\n", in_size);
}

static void
ExpandOutBuffer (int size)
{
  while (out_size < size)
    if (out_size == 0)
      out_size = 1024;
    else
      out_size *= 2;
  out_buffer = realloc(out_buffer, out_size);
  if (out_buffer == NULL)
    Abort("Could not expand out_buffer to %d bytes\n", out_size);
}

#endif	
