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

subroutine InsReaW3(aSGrp,beSGrp,bSGrp,length)
! this routine does:
! Check which W3 file corresponds to given combination of indices,
! - increase the parameter, checking the overall number of
!   required integrals if these block is conted first time
! - and set corresponding InqW3 to T

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpa, DimSGrpbe, InqW3, no
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aSGrp, beSGrp, bSGrp
integer(kind=iwp), intent(inout) :: length
integer(kind=iwp) :: abeSGrp, dima, dimb, dimbe

!0 def dimensions

dima = DimSGrpa(aSGrp)
dimb = DimSGrpa(bSGrp)
dimbe = DimSGrpbe(beSGrp)

if (aSGrp > beSGrp) then
  !1 case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
  abeSGrp = nTri_Elem(aSGrp-1)+beSGrp
  if (.not. InqW3(abeSGrp,bSGrp)) then
    InqW3(abeSGrp,bSGrp) = .true.
    length = length+dima*dimbe*dimb*no
  end if

else if (aSGrp == beSGrp) then
  !2 case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
  abeSGrp = nTri_Elem(aSGrp-1)+beSGrp
  if (.not. InqW3(abeSGrp,bSGrp)) then
    InqW3(abeSGrp,bSGrp) = .true.
    length = length+no*nTri_Elem(dima)*dimb
  end if

else
  !3 case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
  abeSGrp = nTri_Elem(beSGrp-1)+aSGrp
  if (.not. InqW3(abeSGrp,bSGrp)) then
    InqW3(abeSGrp,bSGrp) = .true.
    length = length+dima*dimbe*dimb*no
  end if
end if

return

end subroutine InsReaW3
