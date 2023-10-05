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

use Cholesky, only: DiaMax, DiaMaxT, iAtomShl, iiBstR, iiBstRSh, IndRed, iSP2F, LuPri, nnBstRSh, nnShl, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: iLoc
real(kind=wp), intent(out) :: DGMax
integer(kind=iwp) :: i, i1, i2, iShlA, iShlAB, iShlB, iSym
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_MaxAbsDiag_1C'

if (iLoc == 1) then
  do iSym=1,nSym
    DiaMax(iSym) = Zero
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
    DiaMax(iSym) = Zero
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
    DiaMaxT(iSym) = Zero
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
