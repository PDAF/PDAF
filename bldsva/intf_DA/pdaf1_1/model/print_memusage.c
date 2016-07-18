/*-----------------------------------------------------------------------------------------
Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)

This file is part of TerrSysMP-PDAF

TerrSysMP-PDAF is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TerrSysMP-PDAF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LesserGeneral Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------
print_memusage.c: Print memory usage on Juqueen
-------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "wrapper_tsmp.h"

#include <spi/include/kernel/memory.h>

void print_memusage(int ts){
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;
  int rank;
  char fn[100];
  FILE *fnout;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  strcpy(fn,"memcheck");
  sprintf(fn,"%s_%05d_%05d",fn,rank,ts);
  fnout = fopen(fn,"a");

    
  
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
/*
#if 0
  printf("Allocated heap: %.2f MB, avail. heap: %.2f MB\n", double(heap)/(1024*1024), double(heapavail)/(1024*1024));
  printf("Allocated stack: %.2f MB, avail. stack: %.2f MB\n", double(stack)/(1024*1024), double(stackavail)/(1024*1024));
  printf("Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f MB, mmap: %.2f MB\n", double(shared)/(1024*1024), double(persist)/(1024*1024), double(guard)/(1024*1024), double(mmap)/(1024*1024));
#else
  printf("MEMSIZE heap: %.2f/%.2f stack: %.2f/%.2f mmap: %.2f MB\n", (double)heap/(1024*1024), (double)heapavail/(1024*1024), (double)stack/(1024*1024), (double)stackavail/(1024*1024), (double)mmap/(1024*1024));
  printf("MEMSIZE shared: %.2f persist: %.2f guard: %.2f MB\n", (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024));
#endif
*/
  fprintf(fnout,"MEMSIZE heap: %.2f/%.2f stack: %.2f/%.2f mmap: %.2f MB\n", (double)heap/(1024*1024), (double)heapavail/(1024*1024), (double)stack/(1024*1024), (double)stackavail/(1024*1024), (double)mmap/(1024*1024));
  fprintf(fnout,"MEMSIZE shared: %.2f persist: %.2f guard: %.2f MB\n", (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024));
  fclose(fnout);
}
