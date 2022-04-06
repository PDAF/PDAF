/*BHEADER**********************************************************************

  Copyright (c) 1995-2009, Lawrence Livermore National Security,
  LLC. Produced at the Lawrence Livermore National Laboratory. Written
  by the Parflow Team (see the CONTRIBUTORS file)
  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.

  This file is part of Parflow. For details, see
  http://www.llnl.gov/casc/parflow

  Please read the COPYRIGHT file or Our Notice and the LICENSE file
  for the GNU Lesser General Public License.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License (as published
  by the Free Software Foundation) version 2.1 dated February 1999.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
  and conditions of the GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
**********************************************************************EHEADER*/

#include <stdio.h>
#include <sys/param.h>
#include <sys/param.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/times.h>
//>>TSMP-PDAF addition beginning (from parflow v3.4.0:amps/mpi1/amps_init.c:36)
#include  <inttypes.h>
//<<TSMP-PDAF addition end


#include "amps.h"

//>>TSMP-PDAF addition beginning
MPI_Comm dacomm; /* kuw */
//<<TSMP-PDAF addition end

int amps_mpi_initialized = FALSE;

#ifdef AMPS_MALLOC_DEBUG
char amps_malloclog[MAXPATHLEN];
#endif

#ifndef CRAY_TIME
long AMPS_CPU_TICKS_PER_SEC;
#endif

int amps_size;
int amps_rank;
//>>TSMP-PDAF addition beginning (from parflow v3.4.0:amps/mpi1/amps_init.c:48)
int amps_node_rank;
int amps_node_size;
int amps_write_rank;
int amps_write_size;
MPI_Comm nodeComm = MPI_COMM_NULL;
MPI_Comm writeComm = MPI_COMM_NULL;
//<<TSMP-PDAF addition end

#ifdef AMPS_F2CLIB_FIX
int MAIN__()
{
}
#endif

/*===========================================================================*/
/**

Every {\em AMPS} program must call this function to initialize the
message passing environment.  This must be done before any other
{\em AMPS} calls.  \Ref{amps_Init} does a synchronization on all the
nodes and the host (if it exists).  {\bf argc} and {\bf argv} should be
the {\bf argc} and {\bf argv} that were passed to {\bf main}.  On some
ports command line arguments can be used to control the underlying
message passing system.  For example, on {\em Chameleon} the
{\bf -trace} flag can be use to create a log of communication events.

{\large Example:}
\begin{verbatim}
int main( int argc, char *argv)
{
   amps_Init(argc, argv);
   
   amps_Printf("Hello World");

   amps_Finalize();
}
\end{verbatim}

{\large Notes:}

@memo Initialize AMPS
@param argc Command line argument count [IN/OUT]
@param argv Command line argument array [IN/OUT]
@return 
*/
//>>TSMP-PDAF addition beginning (from parflow v3.4.0:amps/mpi1/amps_init.c:92)
/*Adler32 function to calculate hash based on node name and length */
const int MOD_ADLER = 65521;

uint32_t Adler32(unsigned char *data, size_t len) /* where data is the location of the data in physical memory and 
len is the length of the data in bytes */
{
  uint32_t a = 1, b = 0;
  size_t index;

/* Process each byte of the data in order */
  for (index = 0; index < len; ++index)
    {
      a = (a + data[index]) % MOD_ADLER;
      b = (b + a) % MOD_ADLER;
    }

  return (b << 16) | a;
}
//<<TSMP-PDAF addition end

int amps_Init(int *argc, char **argv[])
{
#ifdef AMPS_MPI_SETHOME
   char *temp_path;
   int length;
#endif

   char processor_name[MPI_MAX_PROCESSOR_NAME];
//>>TSMP-PDAF addition beginning (from parflow v3.4.0:amps/mpi1/amps_init.c:118)
   unsigned char processor_Name[MPI_MAX_PROCESSOR_NAME];
//<<TSMP-PDAF addition end
   int namelen;

//>>TSMP-PDAF comment out beginning
  // MPI_Init(argc, argv);
//<<TSMP-PDAF addition end
   amps_mpi_initialized = TRUE;
//>>TSMP-PDAF addition beginning
   //dacomm = MPI_Comm_f2c(comm_model_pdaf);
//<<TSMP-PDAF addition end

//>>TSMP-PDAF change beginning: MPI_COMM_WORLD -> dacomm
   MPI_Comm_size(dacomm, &amps_size);
   MPI_Comm_rank(dacomm, &amps_rank);
//<<TSMP-PDAF change end

//>>TSMP-PDAF addition beginning (from parflow v3.4.0:amps/mpi1/amps_init.c:127)
  /*Split the node level communicator based on Adler32 hash keys*/
  MPI_Get_processor_name(processor_name,&namelen);
  uint32_t checkSum = Adler32(processor_name, namelen);
  MPI_Comm_split(dacomm, checkSum, amps_rank, &amps_CommNode);
  MPI_Comm_rank(amps_CommNode, &amps_node_rank);
  MPI_Comm_size(amps_CommNode, &amps_node_size);
  int color;
  if (amps_node_rank == 0)
  {
     color = 0;
  }
  else
  { 
    color = 1;
  }
  //MPI_Comm_split(MPI_COMM_WORLD, color, amps_rank, &amps_CommWrite);
  MPI_Comm_split(dacomm, color, amps_rank, &amps_CommWrite);
  if(amps_node_rank == 0)
  {
	MPI_Comm_size(amps_CommWrite, &amps_write_size);
  }
//<<TSMP-PDAF addition end


#ifdef AMPS_STDOUT_NOBUFF
   setbuf (stdout, NULL);
#endif

#ifdef CASC_HAVE_GETTIMEOFDAY
   amps_clock_init();
#endif

#ifdef AMPS_MPI_SETHOME
   if( !amps_rank )
   {
      if( (temp_path = getenv("AMPS_EXE_DIR")) == NULL)
      {
	 printf("AMPS Error: can't get AMPS_EXE_DIR envirnment variabl\n");
	 exit(1);
      }

      length = strlen(temp_path)+1;
   }

//>>TSMP-PDAF change beginning: MPI_COMM_WORLD -> dacomm
   MPI_Bcast(&length, 1, MPI_INT, 0, dacomm);
//<<TSMP-PDAF change end

   if(amps_rank)
   {
      temp_path = malloc(length);
   }

//>>TSMP-PDAF change beginning: MPI_COMM_WORLD -> dacomm
   MPI_Bcast(temp_path, length, MPI_CHAR, 0, dacomm);
//<<TSMP-PDAF change end

   if(chdir(temp_path))
      printf("AMPS Error: can't set working directory to %s", temp_path);

   if(amps_rank)
   {
      free(temp_path);
   }

#endif

#ifdef AMPS_MALLOC_DEBUG
    dmalloc_logpath = amps_malloclog;
    sprintf(dmalloc_logpath, "malloc.log.%04d", amps_Rank(amps_CommWorld));
#endif

#ifdef TIMING
#ifndef CRAY_TIME
   AMPS_CPU_TICKS_PER_SEC = sysconf(_SC_CLK_TCK);
#endif
#endif

#ifdef AMPS_PRINT_HOSTNAME
   MPI_Get_processor_name(processor_name,&namelen);

   printf("Process %d on %s\n", amps_rank, processor_name);
#endif

   return 0;
}  


/*===========================================================================*/
/**

Initialization when ParFlow is being invoked by another application.
This must be done before any other {\em AMPS} calls.

{\large Example:}
\begin{verbatim}
int main( int argc, char *argv)
{
   amps_EmbeddedInit();
   
   amps_Printf("Hello World");

   amps_Finalize();
}
\end{verbatim}

{\large Notes:}

@memo Initialize AMPS
@return 
*/
int amps_EmbeddedInit(void)
{
//>>TSMP-PDAF change beginning: MPI_COMM_WORLD -> dacomm
   MPI_Comm_size(dacomm, &amps_size);
   MPI_Comm_rank(dacomm, &amps_rank);
//<<TSMP-PDAF change end

#ifdef AMPS_STDOUT_NOBUFF
   setbuf (stdout, NULL);
#endif

#ifdef CASC_HAVE_GETTIMEOFDAY
   amps_clock_init();
#endif

#ifdef AMPS_MALLOC_DEBUG
    dmalloc_logpath = amps_malloclog;
    sprintf(dmalloc_logpath, "malloc.log.%04d", amps_Rank(amps_CommWorld));
#endif

#ifdef TIMING
#ifndef CRAY_TIME
   AMPS_CPU_TICKS_PER_SEC = sysconf(_SC_CLK_TCK);
#endif
#endif

   return 0;
}  

//>>TSMP-PDAF addition beginning
/* kuw */
int amps_EmbeddedInit_tsmp(MPI_Comm subcomm)
{
   MPI_Comm_size(subcomm, &amps_size);
   MPI_Comm_rank(subcomm, &amps_rank);

   dacomm = subcomm;
}
//<<TSMP-PDAF addition end
