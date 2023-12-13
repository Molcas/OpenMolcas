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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine page_in(KeyWord,nConf,Vector,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load any vector for further use by the Davidson diagonalization  *
!     which has been saved by the write through cache mechanism        *
!                                                                      *
!     calling arguments:                                               *
!     KeyWord : character(len=llab)                                    *
!               record identifier                                      *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     Vector  : array of real*8                                        *
!               any vector of length nConf                             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use davctl_mod, only: disk_address, llab, LblStk, memory_vectors, mxDiskStk, mxMemStk
use Definitions, only: wp, iwp, u6

implicit none
character(len=llab), intent(in) :: KeyWord
integer(kind=iwp), intent(in) :: nConf, LuDavid
real(kind=wp), intent(out) :: Vector(nConf)
integer(kind=iwp) :: iDisk, iStk, nStk
#include "rasdim.fh"

! check input arguments
if (nConf < 0) then
  write(u6,*) 'page_in: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if

! search for a matching record identifier
nStk = 0
do iStk=1,(mxMemStk+mxDiskStk)
  if (LblStk(iStk) == KeyWord) then
    nStk = iStk
    exit
  end if
end do
if (nStk == 0) then
  write(u6,*) 'page_in: nStk equal 0'
  write(u6,*) 'nStk = ',nStk
  call Abend()
end if

if (nStk <= mxMemStk) then
  Vector(:) = memory_vectors(:,nStk)
else
  iDisk = disk_address(nStk-mxMemStk)
  call DDaFile(LuDavid,2,Vector,nConf,iDisk)
end if

return

end subroutine page_in
