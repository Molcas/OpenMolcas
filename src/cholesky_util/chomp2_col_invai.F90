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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Col_Invai(ai,iSymai,a,iSyma,i,iSymi)
!
! Thomas Bondo Pedersen, Dec. 2007.
!
! Purpose: calculate indices a and i (incl. symmetries)
!          from compound index ai of symmetry iSymai.

use Symmetry_Info, only: Mul
use Cholesky, only: nSym
use ChoMP2, only: iT1am, nOcc, nVir
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ai, iSymai
integer(kind=iwp), intent(out) :: a, iSyma, i, iSymi
integer(kind=iwp) :: iSym, i_, ai_1, ai_2
#ifdef _DEBUGPRINT_
character(len=*), parameter :: SecNam = 'ChoMP2_Col_Invai'
#endif

! Find iSyma and iSymi.
! ---------------------

iSymi = -999999
iSyma = -999999
iSym = nSym+1
do while (iSym > 1)
  iSym = iSym-1
  iSymi = iSym
  iSyma = Mul(iSymi,iSymai)
  if ((nOcc(iSymi) > 0) .and. (nVir(iSyma) > 0) .and. (ai > iT1Am(iSyma,iSymi))) iSym = 0 ! Found! -- break loop
end do

#ifdef _DEBUGPRINT_
if ((iSymi < 1) .or. (iSymi > nSym) .or. (iSyma < 1) .or. (iSyma > nSym)) call SysAbendMsg(SecNam,'bug detected','[1]')
if ((nOcc(iSymi) < 1) .or. (nVir(iSyma) < 1)) call SysAbendMsg(SecNam,'bug detected','[2]')
#endif

! Find a and i.
! -------------

i = -999999
a = -999999
i_ = 0
do while (i_ < nOcc(iSymi))
  i_ = i_+1
  ai_1 = iT1Am(iSyma,iSymi)+nVir(iSyma)*(i_-1)+1
  ai_2 = ai_1+nVir(iSyma)-1
  if ((ai >= ai_1) .and. (ai <= ai_2)) then
    i = i_
    a = ai-ai_1+1
    i_ = nOcc(iSymi)+1 ! Found! -- break loop
  end if
end do

#ifdef _DEBUGPRINT_
if ((i < 1) .or. (i > nOcc(iSymi)) .or. (a < 1) .or. (a > nVir(iSyma))) call SysAbendMsg(SecNam,'bug detected','[3]')
#endif

end subroutine ChoMP2_Col_Invai
