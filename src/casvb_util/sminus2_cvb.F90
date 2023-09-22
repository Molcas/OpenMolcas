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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine sminus2_cvb(bikfrom,bikto,nel,nalffrom,ndetfrom,nalfto,ndetto,nvec,xdetto,ioccfrom,ioccto)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nel, nalffrom, ndetfrom, nalfto, ndetto, nvec, xdetto(0:nel,0:nalfto), ioccfrom(nalffrom), ioccto(nalfto)
real(kind=wp) :: bikfrom(ndetfrom,nvec), bikto(ndetto,nvec)
integer(kind=iwp) :: iexc, indfrom, indto
integer(kind=iwp), external :: minind_cvb

call fzero(bikto,ndetto*nvec)

! Determinant (to) weight array:
call weightfl_cvb(xdetto,nalfto,nel)
if (ndetto /= xdetto(nel,nalfto)) then
  write(u6,*) ' Discrepancy in NDET:',ndetto,xdetto(nel,nalfto)
  call abend_cvb()
end if

call loopstr0_cvb(ioccfrom,indfrom,nalffrom,nel)
do
  call imove_cvb(ioccfrom(2),ioccto,nalfto)
  do iexc=1,nalffrom
    indto = minind_cvb(ioccto,nalfto,nel,xdetto)
    call daxpy_(nvec,One,bikfrom(indfrom,1),ndetfrom,bikto(indto,1),ndetto)
    if (iexc < nalffrom) ioccto(iexc) = ioccfrom(iexc)
  end do
  call loopstr_cvb(ioccfrom,indfrom,nalffrom,nel)
  if (indfrom == 1) exit
end do

return

end subroutine sminus2_cvb
