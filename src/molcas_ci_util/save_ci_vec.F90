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

subroutine Save_CI_vec(iRoot,nConf,CI_vec,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Save a CI vector                                                 *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     iRoot   : integer                                                *
!               root number                                            *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     CI_vec  : array of real                                          *
!               CI vector                                              *
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

use timers, only: TimePage
use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, nkeep, on_disk, save_mode
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iRoot, nConf, LuDavid
real(kind=wp), intent(_IN_) :: CI_vec(nConf)
integer(kind=iwp) :: CI_vec_PageNo, CI_vec_RecNo, iDisk
real(kind=wp) :: dum1, dum2, dum3, Time(2)
character(len=llab) :: KeyWord
integer(kind=iwp), external :: PageNo, RecNo
#include "rasdim.fh"

call Timing(Time(1),dum1,dum2,dum3)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Save_CI_vec: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(u6,*) 'Save_CI_vec: iRoot less than 0'
  write(u6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > nkeep) then
  write(u6,*) 'Save_CI_vec: iRoot greater than nkeep'
  write(u6,*) 'iRoot, nkeep = ',iRoot,nkeep
  call Abend()
end if

! the diagonalization can be run in core:
! copy the CI vector to new memory location
if (save_mode == in_core) then
  CI_vec_RecNo = RecNo(2,iRoot)
  memory_vectors(:,CI_vec_RecNo) = CI_vec(:)
end if

! the diagonalization must be run out of core:
! save the CI vector on disk
if (save_mode == on_disk) then
  CI_vec_RecNo = RecNo(2,iRoot)
  iDisk = disk_address(CI_vec_RecNo)
  call DDaFile(LuDavid,1,CI_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the CI vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  CI_vec_PageNo = PageNo(iRoot)
  write(KeyWord,'(A,I4.4)') 'CI_vec',CI_vec_PageNo
  call page_out(KeyWord,nConf,CI_vec,LuDavid)
end if

call Timing(Time(2),dum1,dum2,dum3)
TimePage = TimePage+Time(2)-Time(1)

return

end subroutine Save_CI_vec
