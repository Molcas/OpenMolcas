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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Sort_mck(A,B,iBas,jBas,kBas,lBas,iCmp,jCmp,kCmp,lCmp,
     &               iBasO,jBasO,kBasO,lBasO,
     &               iCmpO,jCmpO,kCmpO,lCmpO,
     &               nVec,nop,iAng,
     &               indgrd,indgrd2,ishll,C)
************************************************************************
*                                                                      *
*     This subroutine is a stupid solution on a easy problem, but it   *
*     should work and it doesnt take to much CPU time.                 *
*     eaw                                                              *
*                                                                      *
************************************************************************
      Use Basis_Info
      use Real_Spherical, only: iSphCr
      use Symmetry_Info, only: iOper, iChBas
      Implicit Real*8(a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
c#include "print.fh"
*
      Integer nop(4),iAng(4),indgrd(3,4,0:nirrep-1),
     &          indgrd2(3,4,0:nirrep-1),ishll(4)
      Real*8 A(iBas*jBas*kBas*lBas,
     &         iCmp,jCmp,kCmp,lCmp,nVec)
      Real*8 B(kBasO*kcmpO,lBasO,lcmpO,
     &         iBasO,iCmpO,jBasO,jCmpO*nvec)

      Real*8 prmt(0:7)
      Real*8 C(*)
*
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
*
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))
*
      rp=1.0d0
      Do iVec=1,nVec
        Do iC=1,iCmp
         ichbs=ichbas(ii+ic)
         If (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+ic))
         PrA= xPrmt(iOper(nOp(1)),iChBs)
         Do jC=1,jCmp
          jChBs = iChBas(jj+jc)
          If (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+jc))
          pRb = xPrmt(iOper(nOp(2)),jChBs)
          Do kC=1,kCmp
           kChBs = iChBas(kk+kc)
           If (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+kc))
           pTc = xPrmt(iOper(nOp(3)),kChBs)
           Do lC=1,lCmp
            lChBs = iChBas(ll+lC)
            If (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+lc))
            pTSd= xPrmt(iOper(nOp(4)),lChBs)
            qFctr=pTSd*pTc*pRb*pRa*rp
**EAW 970930
*
*        Some machines dont support more than 7 indexes, need to fool around
*        to make it work
*
*           Do iB=1,iBas
*            Do jB=1,jBas
*             Do kB=1,kBas
*              Do lB=1,lBas
*               B(kB,kC,lB,lC,iB,iC,jB,jC,iVec)=
*    &            qfctr*A(iB,jB,kB,lB,iC,jC,kC,lC,iVec)
*              End Do
*             End Do
*            End Do
*           End Do
            ijkl=0
            Do lB=1,lBas
             Do kB=1,kBas
              Do jB=1,jBas
               Do iB=1,iBas
                ijkl=ijkl+1
                B(kB+(kC-1)*kbaso,lB,lC,iB,iC,jB,jC+(iVec-1)*jcmpO)=
     &            qfctr*A(ijkl,iC,jC,kC,lC,iVec)
               End Do
              End Do
             End Do
            End Do

*EAW 970930
           End Do
          End Do
         End Do
        End Do
      End Do
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(indgrd)
         Call Unused_integer_array(indgrd2)
         Call Unused_real_array(C)
      End If
      End
