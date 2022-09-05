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
      SubRoutine Cho_Reorder_RI(Vec,lVec,nVec,iSym)
      use ChoArr, only: iRS2F
      Implicit Real*8 (a-h,o-z)
      Real*8 Vec(lVec,nVec)
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Real*8, Allocatable :: Scr(:)
      Integer, Allocatable :: iF2RS(:)

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      If (nVec .lt. 1) Return
      If (lVec .lt. 1) Return
      If (lVec.ne.nnBstR(iSym,1) .or. nVec.gt.NumCho(iSym)) Then
         Call SysAbendMsg('Cho_Reorder_RI','Input argument error!',' ')
      End If
      If (nnShl .ne. nnShl_Tot) Then
         Call SysAbendMsg('Cho_Reorder_RI','Screening is not allowed!', &
     &                    '(nnShl.ne.nnShl_Tot)')
      End If

!     Set mapping from global address to reduced set.
!     -----------------------------------------------

      liF2RS = nBasT*(nBasT+1)/2
      Call mma_allocate(iF2RS,liF2RS,Label='iF2RS')
      iF2RS(:)=0
      Do iRS = 1,nnBstR(iSym,1)
         iRS_tot = iiBstR(iSym,1) + iRS
         na = iRS2F(1,iRS_tot)
         nb = iRS2F(2,iRS_tot)
         nab = iTri(na,nb)
         iF2RS(nab) = iRS
      End Do

!     Reorder.
!     --------

      lScr = lVec
      Call mma_allocate(Scr,lScr,Label='Scr')
      Do iVec = 1,nVec

         Scr(:)=Vec(:,iVec)
         kFrom = 0
         Do iSymb = 1,nSym

            iSyma = MulD2h(iSymb,iSym)

            If (iSyma .gt. iSymb) Then
               Do ib = 1,nBas(iSymb)
                  nb = iBas(iSymb) + ib
                  Do ia = 1,nBas(iSyma)
                     na = iBas(iSyma) + ia
                     nab = iTri(na,nb)
                     iRS = iF2RS(nab)
#if defined (_DEBUGPRINT_)
                     If (iRS.lt.1 .or. iRS.gt.nnBstR(iSym,1)) Then
                        Call SysAbendMsg('Cho_Reorder_RI',              &
     &                                   'Index out of bounds',' ')
                     End If
#endif
                     kFrom = kFrom + 1
                     Vec(iRS,iVec) = Scr(kFrom)
                  End Do
               End Do
            Else If (iSyma .eq. iSymb) Then
               Do ia = 1,nBas(iSyma)
                  na = iBas(iSyma) + ia
                  Do ib = 1,ia
                     nb = iBas(iSymb) + ib
                     nab = iTri(na,nb)
                     iRS = iF2RS(nab)
#if defined (_DEBUGPRINT_)
                     If (iRS.lt.1 .or. iRS.gt.nnBstR(iSym,1)) Then
                        Call SysAbendMsg('Cho_Reorder_RI',              &
     &                                   'Index out of bounds',' ')
                     End If
#endif
                     kFrom = kFrom + 1
                     Vec(iRS,iVec) = Scr(kFrom)
                  End Do
               End Do
            End If

         End Do

      End Do
      Call mma_deallocate(Scr)
      Call mma_deallocate(iF2RS)

      End
