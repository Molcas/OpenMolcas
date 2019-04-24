************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine SODist(SOValue,mAO,nCoor,mBas,nCmp,
     &                  nDeg,MOValue,iShell,
     &                  nMOs,iAO,CMOs,nCMO,DoIt)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
      Real*8 SOValue(mAO*nCoor,mBas,nCmp*nDeg),
     &       MOValue(mAO*nCoor,nMOs),
     &       CMOs(nCMO)
      Integer DoIt(nMOs)
      Integer   iOff_MO(0:7), iOff_CMO(0:7), iTwoj(0:7)
      Character*80 Label
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      iRout=135
      iPrint=nPrint(iRout)
      Call QEnter('SODist')
      If (iPrint.ge.49) Then
         Write (6,*) 'SODist: MO-Coefficients'
         iOff=1
         Do iIrrep = 0, nIrrep-1
            If (nBas(iIrrep).gt.0) Then
               Write (6,*) ' Symmetry Block',iIrrep
               Call RecPrt(' ',' ',CMOs(iOff),nBas(iIrrep),nBas(iIrrep))
            End If
            iOff=iOff+nBas(iIrrep)**2
         End Do
      End If
*
*---- Compute some offsets
*
      itmp1=1
      itmp2=0
      Do iIrrep = 0, nIrrep-1
         iOff_MO(iIrrep)=itmp1
         iOff_CMO(iIrrep)=itmp2
         itmp1=itmp1+nBas(iIrrep)
         itmp2=itmp2+nBas(iIrrep)*nBas(iIrrep)
      End Do
*

      Do i1 = 1, nCmp
         iDeg=0
         Do iIrrep = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell)+i1),iTwoj(iIrrep)).eq.0)
     &         goto 100
            iDeg=iDeg+1
            iSO=iAOtSO(iAO+i1,iIrrep)
            iOff=(i1-1)*nDeg+iDeg

*
*---------- Distribute contribution to all MO's in this irrep
*
            iMO=iOff_MO(iIrrep)
            iCMO=iOff_CMO(iIrrep)+iSO
c          write(*,*) '> iCMO ',iCMO
            Call MyDGeMM(DoIt(iMO),
     &                 mAO*nCoor,nBas(iIrrep),mBas,
     &                 SOValue(1,1,iOff),mAO*nCoor,
     &                 CMOs(iCMO),nBas(iIrrep),
     &                 MOValue(1,iMO),mAO*nCoor)
 100     continue
          End Do
      End Do
*
      If (iPrint.ge.49) Then
         Write (Label,'(A)')'SODist: MOValue(mAO*nCoor,nMOs)'
         Call RecPrt(Label,' ',MOValue(1,1),mAO*nCoor,nMOs)
      End If
*
      Call GetMem('SODist ','CHEC','REAL',iDum,iDum)
      Call QExit('SODist')
      Return
      End

      SUBROUTINE MYDGEMM ( DoIt, M, N, K,
     $                     A, LDA, B, LDB,
     $                     C, LDC )
*     .. Scalar Arguments ..
      INTEGER            M, N, K, LDA, LDB, LDC
      Integer DoIt(*)
*     .. Array Arguments ..
      REAL*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*  Purpose
*  =======
*
*  DGEMM for a spesial case.
*
*
*     .. Local Scalars ..
      INTEGER            I, J, L
      REAL*8   TEMP
*     .. Parameters ..
      REAL*8   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     ..
*
*           Form  C := A*B + C.
*
            DO 90, J = 1, N
            if(DoIt(J).eq.1) then
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
            endif
   90       CONTINUE
*
      RETURN
      END
*
c      Subroutine SODist1(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
c     &                  iShell,nSOs,iAO,nCMO)
c      Implicit Real*8 (A-H,O-Z)
c#include "itmax.fh"
c#include "info.fh"
c*
c      Real*8  SOValue(mAO*nCoor,mBas,nCmp*nDeg),
c     &        SO(mAO*nCoor,nSOs),
c     &        TmpCMOs(nCMO)
c      Integer TmpDoIt(nSOs)
c*
c      Do k=1,nSOs
c         TmpDoIt(k) = 1
c      End Do
c*
c      call dcopy_(nCMO,0.0d0,0,TmpCMOs,1)
c*
c      iOff=0
c      Do i=0,nIrrep-1
c         iBas=nBas(i)
c         Do j=1,iBas
c            ii = iOff+(j-1)*iBas + j
c            TmpCMOs(ii) = 1.0d0
c         End Do
c         iOff=iOff+nBas(i)**2
c      End Do
c*
c      Call SODist(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
c     &            iShell,nSOs,iAO,TmpCMOs,nCMO,TmpDoIt)
c*
c      Return
c      End
c
*
      Subroutine SODist2(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
     &                  iShell,nSOs,iAO,TmpCMOs,nCMO,TmpDoIt)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*
      Real*8  SOValue(mAO*nCoor,mBas,nCmp*nDeg),
     &        SO(mAO*nCoor,nSOs),
     &        TmpCMOs(nCMO)
      Integer TmpDoIt(nSOs)
*
      Do k=1,nSOs
         TmpDoIt(k) = 1
      End Do
*
      call dcopy_(nCMO,[0.0d0],0,TmpCMOs,1)
*
      iOff=0
      Do i=0,nIrrep-1
         iBas=nBas(i)
         Do j=1,iBas
            ii = iOff+(j-1)*iBas + j
            TmpCMOs(ii) = 1.0d0
         End Do
         iOff=iOff+nBas(i)**2
      End Do
*
      Call SODist(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
     &            iShell,nSOs,iAO,TmpCMOs,nCMO,TmpDoIt)
*
      Return
      End
