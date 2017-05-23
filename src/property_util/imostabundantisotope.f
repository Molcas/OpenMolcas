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
************************************************************************
      Integer Function iMostAbundantIsotope(Z)
************************************************************************
*                                                                      *
* Routine: iMostAbundantIsotope                                        *
* Purpose: This is a wrapper for ixMostAbundantIsotope, to return the  *
*          mass number for the most abundant isotope for atom with     *
*          charge Z.                                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden.                                    *
* Written: March 2000                                                  *
* History: none.                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Z  - The nuclear charge for which the nuclear mass number is         *
*      returned.                                                       *
*                                                                      *
************************************************************************
      Implicit  None
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Integer   Z
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Integer   Rc
      Integer   Opt
      Integer   A
*----------------------------------------------------------------------*
* External references.                                                 *
*----------------------------------------------------------------------*
      Integer   ixMostAbundantIsotope
      external  ixMostAbundantIsotope
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Rc=0
      Opt=0
      A=ixMostAbundantIsotope(Z,Rc,Opt)
      If(Rc.ne.0) Then
      Call SysAbendMsg('imostabundantisotope','Fail to get mass',' ')
      End if
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      iMostAbundantIsotope=A
      Return
      End
