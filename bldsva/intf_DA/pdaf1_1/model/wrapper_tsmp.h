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
wrapper_tsmp.h: Wrapper functions for TerrSysMP (header file)
-------------------------------------------------------------------------------------------*/

#ifndef _WRAPPER_TSMP_H_
#define _WRAPPER_TSMP_H_

/* integer variables */
extern int model;
extern int mype_model;

/* functions */
void initialize_tsmp();
void finalize_tsmp();
void integrate_tsmp();
void print_update_pfb();
void update_tsmp();
void print_memusage(int ts);

#endif
