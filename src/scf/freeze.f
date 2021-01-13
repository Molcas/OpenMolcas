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
      SubRoutine Freeze(TrMat,nTrMat,OneHam,mBT)
************************************************************************
*                                                                      *
*     purpose: Modify transformation matrix such atomic orbitals we    *
*              want to freeze are put in the first position            *
*                                                                      *
*     input:                                                           *
*       TrMat   : transformation matrix of length nTrMat               *
*                                                                      *
*     output:                                                          *
*       TrMat   : modified transformation matrix                       *
*                                                                      *
*     called from: TrGen                                               *
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
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 TrMat(nTrMat), OneHam(mBT)
      Real*8, Dimension(:), Allocatable:: Temp
*
*
*---- Define local variables
      Integer MapBas(MxBas,MxSym)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
*---- Allocate temporary spaces
      Call mma_allocate(Temp,nBT,Label='Temp')
*
      call dcopy_(nBT,OneHam,1,Temp,1)
*
*---- Form an array saying which atomic orbitals are frozen
      iStart = 0
      Do iSym = 1, nSym
         nOr = nOrb(iSym)
         Do iFro = 1, nFro(iSym)
            OEMin = 1.0d+06
            iStrt = iStart
            iSta=0
            Do iBas = 1, nOr
               iStrt = iStrt + iBas
               If (Temp(iStrt).lt.OEMin) Then
                  MapBas(iFro,iSym)=iBas
                  OEMin = Temp(iStrt)
                  iSta = iStrt
               End If
            End Do
            If (iSta.ne.0) Temp(iSta) = -Temp(iSta)
         End Do
         iStart = iStart + nOr*(nOr + 1)/2
         Do i = 1, nFro(iSym) - 1
            k = i
            Do j = i + 1, nFro(iSym)
               If (MapBas(j,iSym).lt.MapBas(k,iSym)) k = j
            End Do
            If (k.ne.i) Then
               iSwap          = MapBas(k,iSym)
               MapBas(k,iSym) = MapBas(i,iSym)
               MapBas(i,iSym) = iSwap
            End If
         End Do
      End Do
*
*---- Move frozen atomic orbitals to the first positions
      iCMO = 1
      Do iSym = 1, nSym
         Do iFro = 1, nFro(iSym)
            ind1 = iCMO + (iFro              - 1)*nBas(iSym)
            ind2 = iCMO + (MapBas(iFro,iSym) - 1)*nBas(iSym)
            Call DSwap_(nBas(iSym),TrMat(ind1),1,TrMat(ind2),1)
         End Do
         iCMO = iCMO + nBas(iSym)*nOrb(iSym)
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(Temp)
*
      End
