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
*               2011, Francesco Aquilante                              *
************************************************************************
      SubRoutine DOne_SCF(nSym,nBas,nOrb,nFro,CMO,nCMO,Occ,Dlt,
     &                    alpha_density)
************************************************************************
*                                                                      *
*     purpose: Compute density matrix in AO basis                      *
*                                                                      *
*     input:                                                           *
*       nSym    : number of symmetries                                 *
*       nBas(i) : number of basis functions (i = 1, nSym)              *
*       nOrb(i) : number of orbitals (i = 1, nSym)                     *
*       nFro(i) : number of frozen orbitals (i = 1, nSym)              *
*       CMO     : molecular orbitals                                   *
*       Occ     : occupation numbers                                   *
*                                                                      *
*       alpha_density: .true. iff alpha MOs were sent in               *
*                                                                      *
*     output:                                                          *
*       Dlt     : density matrix in triangular storrage                *
*                                                                      *
*     called from: DMat, final                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: modified by F. Aquilante, Apr 2011                      *
*                                                                      *
************************************************************************
*
#include "compiler_features.h"
#ifndef POINTER_REMAP
      Use, Intrinsic :: ISO_C_BINDING
#endif
      Implicit Real*8 (a-h,o-z)
*
      Real*8, Target:: CMO(nCMO), Occ(*), Dlt(*)
      Real*8, Pointer:: pCMO(:,:), pOcc(:), pDlt(:)
      Integer nBas(nSym),nOrb(nSym),nFro(nSym)
      Logical alpha_density
*
#include "real.fh"
#include "WrkSpc.fh"
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
      COMMON  / SPAVE_I  / ip_DSc

*
*---- Statement function for triangular storage
      Ind(i,j) = i*(i - 1)/2 + j
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call NrmClc(CMO,nCMO,'DOne_SCF','CMO')
      nOcc=0
      Do iSym = 1, nSym
         nOcc = nOcc + nOrb(iSym)
      End Do
      Call NrmClc(Occ,nOcc,'DOne_SCF','Occ')
#endif
*
      ipCMO  = 1
      ipDlt  = 1
      ipOcc  = 1
      Do iSym = 1, nSym
*
         nBs = nBas(iSym)
         nOr = nOrb(iSym)
         nFr = nFro(iSym)
*
         lth = nBs*(nBs + 1)/2
         Call FZero(Dlt(ipDlt),lth)
         If (nOr.eq.0) Go To 100
*
#ifdef POINTER_REMAP
         pCMO(1:nBs,1:nBs) => CMO(ipCMO:ipCMO+nBs**2-1)
#else
         Call C_F_POINTER(C_LOC(CMO(ipCMO)), pCMO, [nBs,nBs])
#endif
         pOcc => Occ(ipOcc:ipOcc+nOr-1)
         pDlt => Dlt(ipDlt:ipDlt+lth-1)
*
         Do iRow = 1, nBs
            Sum = Zero
            Do i = nFr + 1, nOr
               Sum = Sum + pOcc(i)*pCMO(iRow,i)**2
*              Sum = Sum + pOcc(i)*pCMO(iRow,i)*pCMO(iRow,i)
            End Do
            pDlt(Ind(iRow,iRow)) = Sum
*
            Do iCol = 1, iRow - 1
               Sum = Zero
               Do i = nFr + 1, nOr
                  Sum = Sum + pOcc(i)*pCMO(iRow,i)*pCMO(iCol,i)
               End Do
               pDlt(Ind(iRow,iCol)) = 2.0D0*Sum
            End Do
         End Do
#ifdef _DEBUGPRINT_
         Call NrmClc(pDlt,nBs*(nBs+1)/2,'DOne_SCF','Dlt')
*        Call RecPrt('CMO',' ',pCMO,nBs,nBs)
*        Call RecPrt('Occ',' ',pOcc,1,nOr)
*        Call TriPrt('Dlt',' ',pDlt,nBs)
#endif
*
 100     Continue
*
         Nullify(pCMO,pOcc,pDlt)
         ipCMO = ipCMO + nBs**2
         ipDlt = ipDlt + lth
         ipOcc = ipOcc + nOr
      End Do
*
      xsign=1.0d0
      If (Do_SpinAV) Then
*
         If (alpha_density) xsign=-1.0d0
         iOffD=1
         lOff=0
         Do iSym = 1, nSym
            lth=nBas(iSym)*(nBas(iSym)+1)/2
*
            pDlt => Dlt(iOffD:iOffD+lth-1)
*
            ipDScc=ip_DSc+lOff
            Do j=1,nBas(iSym)
               Do i=1,j-1
                  ji=j*(j-1)/2+i
                  iDSc=ipDScc-1+nBas(iSym)*(j-1)+i
                  pDlt(ji)=pDlt(ji)+xsign*2.0d0*Work(iDSc)
               End Do
               jj=j*(j+1)/2
               iDSc=ipDScc-1+nBas(iSym)*(j-1)+j
               pDlt(jj)=pDlt(jj)+xsign*Work(iDSc)
            End Do
            iOff=iOff+lth
            lOff=lOff+nBas(iSym)**2
         End Do
*
      EndIf
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
