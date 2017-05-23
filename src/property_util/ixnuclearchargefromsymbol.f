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
      Integer Function ixNuclearChargeFromSymbol(Symbol,Rc,Opt)
#include "proputil.fh"
************************************************************************
*                                                                      *
* This function returns the nuclear charge of an atom based on the     *
* chemical symbol.                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
* Written: Feb. 2000                                                   *
* Modified: March 2017, Ignacio Fdez. Galvan (use periodic_table.fh)   *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Parameters.                                                          *
*----------------------------------------------------------------------*
      Integer    StopOnError
      Parameter (StopOnError=_OPT_STOP_ON_ERROR_)
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Character*(*) Symbol
      Integer Rc
      Integer Opt
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
#include "periodic_table.fh"
      Character*2 Sym1,Sym2
      Integer Index
      Integer i
*----------------------------------------------------------------------*
* External references.                                                 *
*----------------------------------------------------------------------*
      External UpCase
*----------------------------------------------------------------------*
* Locate symbol in table.                                              *
*----------------------------------------------------------------------*
      Index=0
      Sym1=AdjustL(Symbol)
      Call UpCase(Sym1)
      Do i=1,Num_Elem
         Sym2=AdjustL(PTab(i))
         Call UpCase(Sym2)
         if(Sym1.eq.Sym2) Index=i
      End Do
*----------------------------------------------------------------------*
* Are we successful.                                                   *
*----------------------------------------------------------------------*
      If(Index.eq.0) Then
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** NuclearChargeBySymbol: error'
         Write(6,'(2a)') '***    unknown atom: ',Symbol
         Write(6,'(a)') '***'
         If(iAnd(Opt,StopOnError).ne.0) Call Quit_OnUserError
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      ixNuclearChargeFromSymbol=Index
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Rc)
      End
