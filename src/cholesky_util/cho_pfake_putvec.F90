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

subroutine Cho_PFake_PutVec(Vec,InfV,nVec,iSym,iV1)

use Cholesky, only: Cho_AdrVec, LuCho_G
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: Vec(*)
integer(kind=iwp), intent(inout) :: InfV(2,*)
integer(kind=iwp), intent(in) :: nVec, iSym, iV1
integer(kind=iwp) :: iAdr, iOpt, iPos, iVec, lTot
character(len=*), parameter :: SecNam = 'Cho_PFake_PutVec'

if (nVec < 1) return

if (CHO_ADRVEC == 1) then
  iOpt = 1
  lTot = sum(InfV(1,iV1:iV1+nVec-1))
  iAdr = InfV(2,iV1)
  call dDAFile(LuCho_G(iSym),iOpt,Vec,lTot,iAdr)
  do iVec=iV1+1,iV1+nVec
    InfV(2,iVec) = InfV(2,iVec-1)+InfV(1,iVec-1)
  end do
else if (CHO_ADRVEC == 2) then
  iPos = 1
  do iVec=iV1,iV1+nVec-1
    iOpt = 1
    lTot = InfV(1,iVec)
    iAdr = InfV(2,iVec)
    call dDAFile(LuCho_G(iSym),iOpt,Vec(iPos),lTot,iAdr)
    iPos = iPos+InfV(1,iVec)
    InfV(2,iVec+1) = iAdr
  end do
else
  call Cho_Quit('Illegal CHO_ADRVEC in '//SecNam,102)
end if

end subroutine Cho_PFake_PutVec
