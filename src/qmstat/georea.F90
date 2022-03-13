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

subroutine GeoRea(nskipp,quantum)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nskipp
logical(kind=iwp) :: quantum
#include "maxi.fh"
#include "qminp.fh"
integer(kind=iwp) :: iDisk, j
real(kind=wp) :: Dum(1)

!----------------------------------------------------------------------*
! Enter.                                                               *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Read!                                                                *
!----------------------------------------------------------------------*
iDisk = 0
!If we are to skip something.
if ((nSkipp /= 0) .and. (iPrint >= 4)) write(u6,*) ' Reading from configuration ',nskipp,'.'
do j=1,nSkipp+1
  if ((j /= 1) .and. (iRead /= 9)) then
    call dDaFile(9,2,Dum,1,iDisk) !Etot
    call dDaFile(9,2,Dum,1,iDisk) !Ract
    call dDaFile(9,2,Dum,1,iDisk) !GamOld
    call dDaFile(9,2,Dum,1,iDisk) !Gam
    call dDaFile(9,2,Dum,1,iDisk) !ESub
  end if
  ! If this is a sampfile we do not care about
  ! the induced dipoles, so we just read them to get rid of them.
  !if (iRead == 9) then
  !  do i=1+nPol,IndMa
  !  end do
  !end if
end do

!----------------------------------------------------------------------*
! Exit.                                                                *
!----------------------------------------------------------------------*
return
! Avoid unused argument warnings
if (.false.) call Unused_logical(quantum)

end subroutine GeoRea
