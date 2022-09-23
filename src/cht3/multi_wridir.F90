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

subroutine multi_wridir(G,lg,ifile,ias,last)
! See multi_readir
!
! PV/LAOG, 22 may 2003.

use ChT3_global, only: IOPT, nblock
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lg, ifile, ias
real(kind=wp), intent(in) :: G(lg)
integer(kind=iwp), intent(out) :: last
integer(kind=iwp) :: iloc, irest, k, kas

iloc = 1
irest = lg
kas = ias

do while (irest > 0)
  k = min(irest,nblock)
  if (kas <= iopt(2)) then
    write(ifile,rec=kas) G(iloc:iloc+k-1)
  else
    write(ifile+1,rec=kas-iopt(2)) G(iloc:iloc+k-1)
  end if
  iloc = iloc+k
  irest = irest-k
  kas = kas+1
end do
last = kas-1

return

end subroutine multi_wridir
