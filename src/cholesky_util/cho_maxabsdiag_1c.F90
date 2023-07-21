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

subroutine Cho_MaxAbsDiag_1C(Diag,iLoc,DGMax)
!
! Specialization for 1-Center approximation: only find max for
! 1-center diagonals.

use ChoArr, only: iSP2F, iAtomShl
use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed

implicit real*8(a-h,o-z)
real*8 Diag(*)
#include "cholesky.fh"
character*17 SecNam
parameter(SecNam='Cho_MaxAbsDiag_1C')
logical LocDbg
#ifdef _DEBUGPRINT_
parameter(LocDbg=.true.)
#else
parameter(LocDbg=.false.)
#endif

if (iLoc == 1) then
  do iSym=1,nSym
    DiaMax(iSym) = 0.0d0
    do iShlAB=1,nnShl
      call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
      if (iAtomShl(iShlA) == iAtomShl(iShlB)) then
        i1 = iiBstR(iSym,1)+iiBstRSh(iSym,iShlAB,1)+1
        i2 = i1+nnBstRSh(iSym,iShlAB,1)-1
        do i=i1,i2
          DiaMax(iSym) = max(DiaMax(iSym),Diag(i))
        end do
      end if
    end do
    DiaMaxT(iSym) = DiaMax(iSym)
  end do
else if ((iLoc == 2) .or. (iLoc == 3)) then
  do iSym=1,nSym
    DiaMax(iSym) = 0.0d0
    do iShlAB=1,nnShl
      call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
      if (iAtomShl(iShlA) == iAtomShl(iShlB)) then
        i1 = iiBstR(iSym,iLoc)+iiBstRSh(iSym,iShlAB,iLoc)+1
        i2 = i1+nnBstRSh(iSym,iShlAB,iLoc)-1
        do i=i1,i2
          DiaMax(iSym) = max(DiaMax(iSym),Diag(IndRed(i,iLoc)))
        end do
      end if
    end do
    DiaMaxT(iSym) = 0.0d0
    do iShlAB=1,nnShl
      call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
      if (iAtomShl(iShlA) == iAtomShl(iShlB)) then
        i1 = iiBstR(iSym,1)+iiBstRSh(iSym,iShlAB,1)+1
        i2 = i1+nnBstRSh(iSym,iShlAB,1)-1
        do i=i1,i2
          DiaMaxT(iSym) = max(DiaMaxT(iSym),Diag(i))
        end do
      end if
    end do
  end do
else
  write(LuPri,*) SecNam,': unknown reduced set, iLoc = ',iLoc
  call Cho_Quit('Unknown reduced set in '//SecNam,104)
end if

DGMax = DiaMax(1)
do iSym=2,nSym
  DGMax = max(DGMax,DiaMax(iSym))
end do

if (LocDbg) then
  write(LuPri,*) SecNam,': in reduced set ',iLoc,':'
  write(LuPri,*) 'DiaMax  = ',(DiaMax(iSym),iSym=1,nSym)
  write(LuPri,*) 'DiaMaxT = ',(DiaMaxT(iSym),iSym=1,nSym)
  write(LuPri,*) 'DGMax   = ',DGMax
end if

end subroutine Cho_MaxAbsDiag_1C
