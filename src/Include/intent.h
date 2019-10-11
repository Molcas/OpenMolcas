/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/

/* This macro exists to properly document the intent of the programmer.
   Unfortunately intent(out) deallocates allocatables upon entering
   procedures.
   Use this intent value if you want to show that a procedure
   does not depend on existing values in an allocatable,
   but depends on the fact that it is allocated. */
#define _OUT_ inout
