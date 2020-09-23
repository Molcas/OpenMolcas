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
* Copyright (C) 1995,2001, Roland Lindh                                *
************************************************************************
      SubRoutine PCMgrd(
#define _CALLING_
#include "grd_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DCR                                                     *
*              XRysg1                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, May '95                                 *
*                                                                      *
*             Modified to PCM gradients September 2001, Lund, by       *
*             R. Lindh.                                                *
************************************************************************
      use PCM_arrays, only: PCM_SQ, PCMTess
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, XCff2D
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"

#include "grd_interface.fh"

*     Local variables
      Integer iDCRT(0:7)
      Real*8 Coori(3,4), CoorAC(3,2), C(3), TC(3)
      Logical NoLoop, JfGrad(3,4)
      Integer iAnga(4), iStb(0:7), JndGrd(3,4), lOp(4), iuvwx(4)
C     Character ChOper(0:7)*3
C     Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 151
      iPrint = nPrint(iRout)
      Call qEnter('PCMgrd')
*
      nRys = nHer
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(nOrdOp)
      If (nip-1.gt.nZeta*nArr) Then
         Write (6,*) 'nip-1.gt.nZeta*nArr'
         Call ErrTra
         Call Abend()
      End If
      nArray = nZeta*nArr - nip + 1
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = nOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      lOp(1) = kOp(1)
      lOp(2) = kOp(2)
*
      ipAOff = ipA
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
         ipAOff = ipAOff + nAlpha
      End Do
*
      ipBOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
      End Do
*
*---- Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
      Do iDAO = 1, nDAO
         Do iZeta = 1, nZeta
            Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
            DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
         End Do
      End Do
      If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*
*---- Generate stabilizor of C, i.e. a center of a tile.
*
      nStb=1
      iStb(0)=0
*
*     Loop over the tiles
*
      Do iTs = 1, nTs
         Q=PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
         NoLoop = Q.eq.Zero
         If (NoLoop) Go To 111
*------- Pick up the tile coordinates
         C(1:3)=PCMTess(1:3,iTs)

         If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
         Call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
         Fact = -Q*DBLE(nStabM) / DBLE(LmbdT)
*
         Call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
*
         iuvwx(3) = nStb
         iuvwx(4) = nStb
         Call ICopy(6,IndGrd,1,JndGrd,1)
         Do i = 1, 3
            Do j = 1, 2
               JfGrad(i,j) = IfGrad(i,j)
            End Do
         End Do
*
*------- No derivatives with respect to the third or fourth center.
*        The positions of the points in the external field are frozen.
*
         Call ICopy(3,[0],0,JndGrd(1,3),1)
         JfGrad(1,3) = .False.
         JfGrad(2,3) = .False.
         JfGrad(3,3) = .False.
         Call ICopy(3,[0],0,JndGrd(1,4),1)
         JfGrad(1,4) = .False.
         JfGrad(2,4) = .False.
         JfGrad(3,4) = .False.
         mGrad = 0
         Do iCar = 1, 3
            Do i = 1, 2
               If (JfGrad(iCar,i)) mGrad = mGrad + 1
            End Do
         End Do
         If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
         If (mGrad.eq.0) Go To 111
*
         Do lDCRT = 0, nDCRT-1
            lOp(3) = NrOpr(iDCRT(lDCRT))
            lOp(4) = lOp(3)
            Call OA(iDCRT(lDCRT),C,TC)
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
*
*
*           Compute integrals with the Rys quadrature.
*
            nDiff=1
            mRys=(la+lb+2+nDiff+nOrdOp)/2
            Eta=One
            EInv=One
            Call Rysg1(iAnga,mRys,nZeta,
     &                 Array(ipA),Array(ipB),[One],[One],
     &                 Zeta,ZInv,nZeta,
     &                 [Eta],[EInv],1,
     &                 P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                 Array(nip),nArray,
     &                 TNAI1,Fake,XCff2D,
     &                 Array(ipDAO),nDAO*nElem(nOrdOp),
     &                 Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)
*
*           Call RecPrt(' In PCMgrd:Grad',' ',Grad,nGrad,1)
         End Do  ! End loop over DCRs
*
111      Continue
      End Do     ! End loop over centers in the external field
*
      Call qExit('PCMgrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Final)
         Call Unused_integer(nRys)
         Call Unused_real_array(Ccoor)
         Call Unused_integer_array(lOper)
      End If
      End
