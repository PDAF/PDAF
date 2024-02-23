/*-----------------------------------------------------------------------------------------
Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)

This file is part of TSMP-PDAF

TSMP-PDAF is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TSMP-PDAF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LesserGeneral Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------
problem_saturationtopressure.h: Functions for converting saturation to pressure in ParFlow
                                (header file)
-------------------------------------------------------------------------------------------*/

typedef void (*SaturationToPressureInvoke) (Vector *phase_saturation , Vector *phase_pressure , Vector *phase_density , double gravity , ProblemData *problem_data , int fcn , int type);
typedef PFModule *(*SaturationToPressureInitInstanceXtraInvoke) (Grid *grid , double *temp_data );

/* problem_saturation.c */
void SaturationToPressure (Vector *phase_saturation , Vector *phase_pressure , Vector *phase_density , double gravity , ProblemData *problem_data , int fcn ,int type);
PFModule *SaturationToPressureInitInstanceXtra (Grid *grid , double *temp_data );
void SaturationToPressureFreeInstanceXtra (void );
PFModule *SaturationToPressureNewPublicXtra (void );
void SaturationToPressureFreePublicXtra (void );
int SaturationToPressureSizeOfTempData (void );


