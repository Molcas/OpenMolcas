!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine SODist(SOValue,mAO,nCoor,mBas,nCmp,nDeg,MOValue,
     &                  nMOs,iAO,CMOs,nCMO,DoIt)
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      use Constants
      Implicit None
      Integer mAO, nCoor, mBas, nCmp, nDeg, nCMO, nMOs
      Real*8 SOValue(mAO*nCoor,mBas,nCmp*nDeg),
     &       MOValue(mAO*nCoor,nMOs),CMOs(nCMO)
      Integer DoIt(nMOs)

      Integer   iOff_MO(0:7), iOff_CMO(0:7)
      Integer iIrrep, itmp1, itmp2, i1, iDeg, iSO, iOff, iMO, iCMO, iAO
#ifdef _DEBUGPRINT_
      Character(LEN=80) Label
#endif
!
#ifdef _DEBUGPRINT_
      Write (6,*) 'SODist: MO-Coefficients'
      iOff=1
      Do iIrrep = 0, nIrrep-1
         If (nBas(iIrrep).gt.0) Then
            Write (6,*) ' Symmetry Block',iIrrep
            Call RecPrt(' ',' ',CMOs(iOff),nBas(iIrrep),nBas(iIrrep))
         End If
         iOff=iOff+nBas(iIrrep)**2
      End Do
#endif
!
!---- Compute some offsets
!
      itmp1=1
      itmp2=0
      Do iIrrep = 0, nIrrep-1
         iOff_MO(iIrrep)=itmp1
         iOff_CMO(iIrrep)=itmp2
         itmp1=itmp1+nBas(iIrrep)
         itmp2=itmp2+nBas(iIrrep)*nBas(iIrrep)
      End Do
!

      Do i1 = 1, nCmp
         iDeg=0
         Do iIrrep = 0, nIrrep-1
            iSO=iAOtSO(iAO+i1,iIrrep)
            If (iSO<0) Cycle
            iDeg=iDeg+1
            iOff=(i1-1)*nDeg+iDeg

!
!---------- Distribute contribution to all MO's in this irrep
!
            iMO=iOff_MO(iIrrep)
            iCMO=iOff_CMO(iIrrep)+iSO
            Call MyDGeMM(DoIt(iMO),
     &                 mAO*nCoor,nBas(iIrrep),mBas,
     &                 SOValue(1,1,iOff),mAO*nCoor,
     &                 CMOs(iCMO),nBas(iIrrep),
     &                 MOValue(1,iMO),mAO*nCoor)
          End Do
      End Do
!
#ifdef _DEBUGPRINT_
      Write (Label,'(A)')'SODist: MOValue(mAO*nCoor,nMOs)'
      Call RecPrt(Label,' ',MOValue(1,1),mAO*nCoor,nMOs)
#endif
!
      End Subroutine SODist

      SUBROUTINE MYDGEMM ( DoIt, M, N, K,
     $                     A, LDA, B, LDB,
     $                     C, LDC )
      Use Constants, only: Zero
      Implicit None
!     .. Scalar Arguments ..
      INTEGER            M, N, K, LDA, LDB, LDC
      Integer DoIt(*)
!     .. Array Arguments ..
      REAL*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!  Purpose
!  =======
!
!  DGEMM for a spesial case.
!
!
!     .. Local Scalars ..
      INTEGER            J, L
!
!     Form  C := A*B + C.
!
      DO J = 1, N
         If (DoIt(J).ne.1) Cycle
         DO L = 1, K
            IF( B( L, J ).EQ.ZERO ) Cycle
            Call DAxPy_(M,B(L,J),A(:,L),1,C(:,J),1)
         End Do
      End Do
!
      END SUBROUTINE MYDGEMM
