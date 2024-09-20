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

subroutine Mk_SO2cI(SO2cI,iSO2Shell,nSOs)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSOs, iSO2Shell(nSOs)
integer(kind=iwp), intent(out) :: SO2cI(2,nSOs)
integer(kind=iwp) :: Indx, iShell, iSO, nShell

!                                                                      *
!***********************************************************************
!                                                                      *
! Set up table SO to contiguous index over the shell

!write(u6,*) 'Enter SO2cI'
SO2cI(:,:) = 0
call Nr_Shells(nShell)
do iShell=1,nShell

  ! Generate contigues index for this shell

  Indx = 0
  do iSO=1,nSOs
    if (iSO2Shell(iSO) == iShell) then
      Indx = Indx+1
      SO2cI(1,iSO) = Indx
    end if
  end do

  ! Store dimension for this shell

  do iSO=1,nSOs
    if (iSO2Shell(iSO) == iShell) SO2cI(2,iSO) = Indx
  end do

end do ! iShell
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_SO2cI
