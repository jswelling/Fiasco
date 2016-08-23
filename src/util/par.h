/*
 *	Parallel Task Manager Library
 *
 *	Copyright (c) 1995,1996,1997,2001 Pittsburgh Supercomputing Center
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
 *	History:
 *		3/96: Written by Greg Hood
 *		11/95 - Written by Greg Hood
 *		6/96  - Added capability for workers to return results to
 *			master (ghood@psc.edu)
 *		2/97  - cleaned up for release, added several new functions,
 *			and made group names passed through environment
 *			variables rather than by arguments (ghood@psc.edu)
 *		11/99 - added broadcast feature for sending large contexts
 *              3/01  - Converted to use either PVM or MPI (ghood@psc.edu)
 *
 *	This library implements a single-master/multiple-worker style
 *	of parallel processing.  Both the master and worker programs
 *	(which typically are instances of the same code) should call
 *	par_process upon startup.  In this call they supply the
 *	routines to be used for task execution as well as communication.
 *
 *	The user-supplied (*master_task) procedure farms out individual
 *	tasks to workers by setting up a Task record and then calling
 *	par_delegate_task.  When a worker has completed a task, it sets up
 *	a Result record before returning from (*worker_task).
 *	(*master_result) is then automatically called on the master
 *	to process the values in the Result record.
 *
 */

/* the following constants may be modified if necessary,
   but the defaults should work for most sites */
#define PAR_MAX_HOSTS		64	/* maximum # of distinct hosts that
					   may be used at once */
#define PAR_MAX_WORKERS_PER_HOST 256	/* maximum # of workers on a single
					   host */
#define PAR_MAX_WORKERS		256	/* maximum # of workers in the entire
					   collection of hosts */
#define PAR_MAX_OVERLAPPED_BROADCASTS	8   /* maximum # of broadcasts that may
					       be issued before waiting for
					       acknowledgements */


/*----------- nothing beyond this point----------------*/
/*----------- should have to be modified---------------*/
/*---------------for installation----------------------*/

#define PAR_FOREVER	(-1.0)	/* a negative waiting time indicates
				   that we should wait forever */

typedef int Par_Task;		/* identifies a specific task that was doled
				   out by the master */

/*-----------------------------------------------------------------------
 *	APPLICATION PROGRAM INTERFACE FUNCTIONS
 *-----------------------------------------------------------------------*/

/* par_process is the function that should be called in the user's main()
   program to register the handlers, and transfer control over to the
   libpar library */
extern void par_process (int argc, char **argv, char **envp,
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
			 void (*unpack_result)());

/* par_set_hosts tells the library what hosts to use and how many
   workers to start on each host; the format of s is
   "<hostname>:<worker count> <hostname>:<worker count>";
   if a ":<worker count>" portion is omitted, it is assumed to be 1;
   if the worker count is negative, this indicates that no workers should
   be automatically started on the host, but if the user manually adds it
   to the PVM configuration, then we should start abs(worker count) workers
   on it; this function is only usable on PVM-based compilations of libpar.c,
   in which case it may be called at any time, but it will have no effect
   on workers already started  */
extern void par_set_hosts (char *s);

/* par_set_context informs the library that subsequent tasks
   should be performed using the current context; the library
   immediately calls the user's (*pack_context)() to gather up
   the current context and then saves it in a PVM buffer */
extern void par_set_context ();

/* par_broadcast_context informs the library that it should broadcast
   the current context to the workers; this is useful when the context is
   so large that it must be sent in parts; the library will wait for any
   outstanding tasks to complete before issuing the broadcast */
extern void par_broadcast_context ();

/* par_delegate_task causes the current task to be packed up
   (via a call to the user's (*pack_task)()) and delegated to the
   first available worker */
extern Par_Task par_delegate_task ();

/* par_finish waits for all delegated tasks to finish */
extern void par_finish ();

/* par_wait blocks for at most timeout seconds; it returns when
   either some event has happened (such as a task completing) or
   the timeout has expired; the return value is 1 if an event happened
   or 0 if the function timed out */
extern int par_wait (float timeout);

/* par_restart can be called by a worker when it has painted
   itself into a corner and wants to gracefully bail out before
   exiting; the master will then automatically restart another
   worker on the same host to take it's place; this call only works
   on PVM-based compilations of libpar.c */
extern void par_restart ();

/* par_finished returns 1 if the specified task has been
   completed, otherwise 0 */
extern int par_finished (Par_Task id);

/* par_finished_up_to returns 1 if all tasks up to and including
   the specified task have been completed; otherwise it returns 0 */
extern int par_finished_up_to (Par_Task id);

/* par_tasks_outstanding returns the number of tasks that have
   been delegated but have not yet been completed */
extern int par_tasks_outstanding ();

/* par_next_task returns the id which will be assigned to the next
   delegated task */
extern Par_Task par_next_task();

/* par_enabled returns 1 if parallelism has been enabled for this run;
   otherwise 0 */
extern int par_enabled ();

/* par_workers returns the number of workers that are currently
   in the configuration; this call is only valid in the master process */
extern int par_workers ();

/* par_instance returns the instance number of the process and may be
   called by both master and worker processes; it will return 0 for the
   master, 1 through n_workers for the worker processes, and -1 if the
   application is not running in parallel and the current process is
   playing the role of both master and worker. */
extern int par_instance ();

/* par_arch returns an architecture-specific string; for PVM-based compilations
   this will be $PVM_ARCH, otherwise the value printed by "uname -s" */
extern const char* par_arch ();

/* the following routines pack individual values into a message; note that
   these parameters are passed by value, not by pointer as with the PVM
   pk* routines */
void par_pkbyte(unsigned char v);
void par_pkshort(short v);
void par_pkint(int v);
void par_pklong(long v);
void par_pkfloat(float v);
void par_pkdouble(double v);
void par_pkstr(char *v);

/* the following routines pack arrays into a message */
void par_pkbytearray(unsigned char *p, int n);
void par_pkshortarray(short *p, int n);
void par_pkintarray(int *p, int n);
void par_pklongarray(long *p, int n);
void par_pkfloatarray(float *p, int n);
void par_pkdoublearray(double *p, int n);

/* the following routines unpack individual values from a message */
unsigned char par_upkbyte();
short par_upkshort();
int par_upkint();
long par_upklong();
float par_upkfloat();
double par_upkdouble();
void par_upkstr(char *s);

/* the following routines unpack arrays from a message */
void par_upkbytearray(unsigned char *p, int n);
void par_upkshortarray(short *p, int n);
void par_upkintarray(int *p, int n);
void par_upklongarray(long *p, int n);
void par_upkfloatarray(float *p, int n);
void par_upkdoublearray(double *p, int n);
