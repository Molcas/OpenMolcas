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

subroutine sumabdistt(n,idtot)
! this routine distributes work for n records among
! nprocab processors with the frequency, corresponding
! to ideffab values
!
! n     - # of records to be distributed (I)
! idtot - distribution vector (O)
!         (idtot(i) -  # of records to be realized by i-th node)

use ccsd_global, only: ideffab, nprocab
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(out) :: idtot(nprocab)
integer(kind=iwp) :: i, imax, max_, ntot
real(kind=wp) :: rsum

!1 distribute recordsc according to eff. coefficients

rsum = sum(ideffab(1:nprocab))

idtot(:) = int(((ideffab(1:nprocab)*n)/rsum)+Half)

!2 do corrections, if roundoff errors caused some differences

do
  ntot = sum(idtot(:))

  if (ntot > n) then
    ! ubrat treba (z najvacsieho dielu)
    max_ = idtot(1)
    imax = 1
    do i=1,nprocab
      if (max_ < idtot(i)) then
        max_ = idtot(i)
        imax = i
      end if
    end do
    idtot(imax) = idtot(imax)-1
  else if (ntot < n) then
    ! pridat treba (k najvacsiemu dielu)
    max_ = idtot(1)
    imax = 1
    do i=1,nprocab
      if (max_ < idtot(i)) then
        max_ = idtot(i)
        imax = i
      end if
    end do
    idtot(imax) = idtot(imax)+1
  else
    exit
  end if
end do

return

end subroutine sumabdistt
