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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine GetInf(Info,nInfo,DoRys,nDiff,icase)
************************************************************************
*                                                                      *
* Object: to read all input information on the file INFO.              *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OpnCom                                                  *
*              ClsCom                                                  *
*              RdCom                                                   *
*              SetUp_RW                                                *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
      use Real_Spherical
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "nq_info.fh"
#include "status.fh"
      Logical DoRys
      Integer iix(2)
      Real*8 rix(2)
      Logical Found
#include "SysDef.fh"
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
*     Prologue
*
*     Call qEnter('GetInf')
*
      ioffr=0
*
*     Load the static input area.
*
      Call Get_Info_Static(ioffr)
*
*     Load the dynamic input area.
*
      Call Get_Info_Dynamic(Info,nInfo,ioffr,icase)
*                                                                      *
************************************************************************
*                                                                      *
*     Reaction field parameters
*
      Len = ilLoc(lRFEnd)-ilLoc(lRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('RFlInfo',lRFStrt,Len)
*
      Len = idLoc(rRFEnd)-idLoc(rRFStrt)
      Len = (Len+nByte_r)/nByte_r
      Call Get_dArray('RFrInfo',rRFStrt,Len)
*
      Len = iiLoc(iRFEnd)-iiLoc(iRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('RFiInfo',iRFStrt,Len)
*
      Len = iiLoc(cRFEnd)-iiLoc(cRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('RFcInfo',cRFStrt,Len)
*                                                                      *
************************************************************************
*                                                                      *
*     Numerical integration information and parameters
*
      Len = idLoc(rQEnd)-idLoc(rQStrt)
      Len = (Len+nByte_r)/nByte_r
      Call Get_dArray('Quad_r',rQStrt,Len)
*
      Len = iiLoc(iQEnd)-iiLoc(iQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('Quad_i',iQStrt,Len)
*
      Len = iiLoc(cQEnd)-iiLoc(cQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('Quad_c',cQStrt,Len)
*                                                                      *
************************************************************************
*                                                                      *
*     Load the transformation matrices
*
      If (iAngMx-1.ge.lMax) Then
         Call qpg_dArray('SewTInfo',Found,Len2)
         If(.not.Found .or.  Len2.eq.0) Then
            Call SysAbendMsg('getinf','Did not find:','SewTInfo')
         End If
         Len=Len2
*        mod by M.Schuetz: LenSph used to broadcast transformation
*        matrices to servers (parallel distributed SCF)
*        LenSph is a member of IInfo common block
         LenSph=Len
         If (Allocated(RSph)) Then
             Call WarningMessage(2,'RSph already allocated!')
             Call Quit_OnUserError()
         End If
         Call mma_allocate(RSph,LenSph,label='RSph')
         Call mma_allocate( ipSph,[0,iAngMx],label='ipSph')
         ipSph(0)=1
         Do 2 iAng = 0, iAngMx-1
            ipSph(iAng+1)= ipSph(iAng) + (iAng*(iAng+1)/2 + iAng + 1)**2
 2       Continue
         Call Get_dArray('SewTInfo',RSph(ipSph(0)),Len)
      Else
         Call Sphere(lMax)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Set the highest number of differentiations which will be
*     applied to the basis functions. In this case 2 + 1 ( the
*     kinetic energy operator and a differentiaion with respect
*     to a nuclear coordinate.
*
      nPrp=Max(lMax,3)
*
*     Setup of tables for coefficients of the Rys roots and weights.
*
      If (iAngMx.eq.0) nDiff=2
      If (DKroll.and.nOrdEF.gt.0) nDiff=nDiff+nOrdEF
      If (.Not.Test) Call Setup_RW(DoRys,nDiff)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up for contracted calculation
*
      Call Flip_Flop(.False.)
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_EFP()
*                                                                      *
************************************************************************
*                                                                      *
*     Epilogue, end
*
*     Call qExit('GetInf')
      Return
      End
