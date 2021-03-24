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

subroutine PriTim(TimDat,nFld_par,nProcs)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFld_par, nProcs
real(kind=iwp), intent(in) :: TimDat(2*nFld_par,nProcs)
integer(kind=iwp) :: iFld, inode

write(u6,'(1X,A)') 'Timing statistics of individual servers:'
write(u6,'(1X,A5,6X,9A13)') ' node','   wait Dens.','     k2 stuff','     k4 stuff','  reduce Grad','    T O T A L', &
                            '      # tasks','  # shl quads'
do inode=1,nProcs
  write(u6,'(1X,I5,A6,5F13.2,2F13.0)') inode,' CPU  ',(TimDat(iFld,inode),iFld=1,nFld_par-1),TimDat(nFld_par,inode), &
                                       TimDat(2*nFld_par,inode)
  write(u6,'(1X,5X,A6,9F13.2)') ' Wall ',(TimDat(nFld_par+iFld,inode),iFld=1,nFld_par-1)
end do

return

end subroutine PriTim
