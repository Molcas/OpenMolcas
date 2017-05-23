************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Per-Olof Widmark                                 *
************************************************************************
      Real*8 Function NucExp(Z,A)
************************************************************************
*                                                                      *
* Routine: NucExp.                                                     *
* Purpose: This routine computes a nuclear radius in the form of a     *
*          gaussian exponent. The exponent is a function of the        *
*          nuclear charge (Z) and the nuclear mass number (A).         *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden.                                    *
* Written: September 1999                                              *
* History: none.                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Z  - The nuclear charge for the nucleus.                             *
* A  - The nuclear mass number for the nucleus.                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Algorithm:                                                           *
*                                                                      *
************************************************************************
      Implicit None
c      Real*8 CONST_BOHR_RADIUS_IN_SI_
#include "constants.fh"
*----------------------------------------------------------------------*
* Parameter list.                                                      *
*----------------------------------------------------------------------*
      Integer Z,A
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Real*8 A3, R, Xi
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      A3    = 1.0d0*DBLE(A)**(1.0d0/3.0d0)
*
      R=0.836d0*A3+0.570d0             ! fm
      R = R * 1.0D-15                  ! m
      R = R / CONST_BOHR_RADIUS_IN_SI_ ! bohr
      Xi=1.5D0/R**2
      NucExp=Xi
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Z)
      End
