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

subroutine getAOs2(lhigh)
!bs get expansions of atomic orbitals in contracted functions

use AMFI_global, only: AOcoeffs, charge, Lmax, noccorb, occup
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lhigh
integer(kind=iwp) :: closedshells(0:Lmax), i, lrun, openshells(0:Lmax)

call getocc_ao(int(charge),closedshells,openshells)
AOcoeffs(:,:,0:lhigh) = Zero
!BS write(u6,*) 'Orbitals for mean-field'
do lrun=0,lhigh
  !BS write(u6,'(A3,I3)') 'L= ',lrun
  occup(1:closedshells(lrun),lrun) = Two
  do i=1,closedshells(lrun)
    AOcoeffs(i,i,lrun) = One
  end do
  noccorb(lrun) = closedshells(lrun)
  if (openshells(lrun) > 0) then
    i = closedshells(lrun)+1
    occup(i,lrun) = real(openshells(lrun),kind=wp)/real(lrun+lrun+1,kind=wp)
    AOcoeffs(i,i,lrun) = One
    noccorb(lrun) = i
  end if
  !BS if (noccorb(lrun) > 0) then
  !BS   write(u6,'(A,I3)') 'number of orbitals ',noccorb(lrun)
  !BS   do iorbital=1,noccorb(lrun)
  !BS     write(u6,'(A,8F8.4)') 'OCCUPATION: ',(occup(iorbital,lrun),iorbital=1,noccorb(lrun))
  !BS   end do
  !BS end if
end do

return

end subroutine getAOs2
