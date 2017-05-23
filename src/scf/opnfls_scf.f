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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      Subroutine OpnFls_SCF
************************************************************************
*                                                                      *
*     purpose: Open files needed by SCF                                *
*                                                                      *
*     called from: SCF                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
#include "file.fh"

#include "mxdm.fh"
#include "infscf.fh"
*
*---- Define local variables
      Logical test
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
      Call qEnter('OpnFls')
#endif
*
*---  open two-electron integral file ---------------------------------*
      Call f_Inquire(FnOrd,test)
      Call DecideOnDirect(.True.,test,DSCF,DoCholesky)
      If (.Not.DSCF .And. .Not.DoCholesky) Then
*        InVec=0
         iRc=-1
         iOpt=0
         Call OpnOrd(iRC,iOpt,FnOrd,LuOrd)
         If (iRc.ne.0) Then
            Write (6,*) 'OpnFls: Error opening ORDINT'
            Call QTrace
            Call Abend()
         End If
      End If
*
*---  open DNSMAT, dVxcdR, TWOHAM and GRADIENT ------------------------*
      Call DAName(LuDSt,FnDSt)
      Call DAName(LuOSt,FnOSt)
      Call DAName(LuTSt,FnTSt)
      Call DAName(LuGrd,FnGrd)
*
*---  open 2nd order update files      --------------------------------*
      Call DAName(LuDGd,FnDGd)
      Call DAName(Lux,Fnx)
      Call DAName(LuDel,FnDel)
      Call DAName(Luy,Fny)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
      Call qExit('OpnFls')
#endif
      Return
      End
