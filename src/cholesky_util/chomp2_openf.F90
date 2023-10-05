!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_OpenF(iOpt,iTyp,iSym)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
!          delete (iOpt=3) Cholesky vector files for MP2 program
!          (full vectors).
!          For iOpt=0, the units are initialized (to -1).
!          iTyp=1: transformed Cholesky vectors.
!          iTyp=2: vectors from (ai|bj) decomposition.

use ChoMP2, only: DoDens, lUnit_F, nPQ_prod, nT1am, nTypF
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iOpt, iTyp, iSym
character(len=4) :: FullNm
character(len=3) :: BaseNm
character(len=*), parameter :: SecNam = 'ChoMP2_OpenF'

if ((iTyp < 1) .or. (iTyp > nTypF)) call SysAbendMsg(SecNam,'iTyp error',' ')

! Initialize units and return for iOpt=0.
! ---------------------------------------

if (iOpt == 0) then
  lUnit_F(iSym,iTyp) = -1
  return
end if

! Open or close files.
! --------------------

if (iOpt == 1) then
  if ((nT1am(iSym) > 0) .or. (DoDens .and. (nPQ_prod(iSym) > 0))) then
    if (lUnit_F(iSym,iTyp) < 1) then
      call ChoMP2_GetBaseNm(BaseNm,iTyp)
      write(FullNm,'(A3,I1)') BaseNm,iSym
      lUnit_F(iSym,iTyp) = 7
      call daName_MF_WA(lUnit_F(iSym,iTyp),FullNm)
    end if
  else
    lUnit_F(iSym,iTyp) = -1
  end if
else if (iOpt == 2) then
  if (lUnit_F(iSym,iTyp) > 0) then
    call daClos(lUnit_F(iSym,iTyp))
    lUnit_F(iSym,iTyp) = -1
  end if
else if (iOpt == 3) then
  if (lUnit_F(iSym,iTyp) > 0) then
    call daEras(lUnit_F(iSym,iTyp))
    lUnit_F(iSym,iTyp) = -1
  end if
else
  call SysAbendMsg(SecNam,'iOpt out of bounds',' ')
end if

end subroutine ChoMP2_OpenF
