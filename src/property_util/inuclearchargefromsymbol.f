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
      Integer Function iNuclearChargeFromSymbol(Symbol)
************************************************************************
*                                                                      *
* Routine: iNuclearChargeFromSymbol                                    *
* Purpose: Wrapper for ixNuclearChargeFromSymbol, to give the nuclear  *
*          charge for atom with symbol 'Symbol'.                       *
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
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Character*(*) Symbol
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Integer  Z
      Integer  Opt
      Integer  Rc
*----------------------------------------------------------------------*
* External references                                                  *
*----------------------------------------------------------------------*
      Integer   ixNuclearChargeFromSymbol
      External  ixNuclearChargeFromSymbol
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Opt=0
      Rc=0
      Z=ixNuclearChargeFromSymbol(Symbol,Rc,Opt)
      If(Rc.ne.0) Then
      Call SysAbendMsg('inuclearchargefromsymbol',
     & 'Fail to get nuclear charge',' ')
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      iNuclearChargeFromSymbol=Z
      Return
      End
