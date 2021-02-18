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
* Copyright (C) 2000, Per-Olof Widmark                                 *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Integer Function ixMostAbundantIsotope(Z,Rc,Opt)
#include "proputil.fh"
************************************************************************
*                                                                      *
* Routine: ixMostAbundantIsotope                                       *
* Purpose: The mass number for the most abundant isotope of the atom   *
*          with charge Z is returned. For radioactive atoms, the       *
*          most stable isotope is returned.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden.                                    *
* Written: March 2000                                                  *
* History: March 2017, use isotopes module, Ignacio Fdez. Galvan       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Algorithm: The nuclear mass number is tabulated for most charges.    *
*            For Z=0, A=1 is returned, corresponding to a neutron.     *
*            For untabulated Z, A=176+Z.                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Z   - The nuclear charge for which the nuclear mass number is        *
*       returned.                                                      *
* Rc  - Return code.                                                   *
* Opt - Options.                                                       *
*                                                                      *
************************************************************************
      Use isotopes, only: Initialize_Isotopes, ElementList, MaxAtomNum
      Implicit None
*----------------------------------------------------------------------*
* Parameters.                                                          *
*----------------------------------------------------------------------*
      Integer    StopOnError
      Parameter (StopOnError=_OPT_STOP_ON_ERROR_)
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Integer Z
      Integer Rc
      Integer Opt
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Integer A
*----------------------------------------------------------------------*
* Compute A.                                                           *
*----------------------------------------------------------------------*
      Call Initialize_Isotopes()
      If(Z.lt.0) Then
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** ixMostAbundantIsotope: error'
         Write(6,'(a)') '***    Charge less than zero!'
         Write(6,'(a)') '***'
         If(iAnd(Opt,StopOnError).ne.0) Call Quit_OnUserError
         A=1
      Else If(Z.eq.0) Then
         A=1
      Else If(Z.gt.MaxAtomNum) Then
         A=176+Z
      Else
         A=ElementList(Z)%Isotopes(1)%A
      End If
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      ixMostAbundantIsotope=A
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Rc)
      End
