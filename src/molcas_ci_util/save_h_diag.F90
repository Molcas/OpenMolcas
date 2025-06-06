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

subroutine Save_H_diag(nConf,H_diag,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Save the diagonal approximation of the CI Hamiltonian for        *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     H_diag  : array of real                                          *
!               diagonal approximation of the CI Hamiltonian           *
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
use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, on_disk, save_mode
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nConf, LuDavid
real(kind=wp), intent(_IN_) :: H_diag(nConf)
integer(kind=iwp) :: H_diag_RecNo, iDisk
real(kind=wp) :: dum1, dum2, dum3, Time(2)
character(len=llab) :: KeyWord
integer(kind=iwp), external :: RecNo
#include "rasdim.fh"

call Timing(Time(1),dum1,dum2,dum3)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Save_H_diag: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if

! the diagonalization can be run in core:
! copy vector to new memory location
if (save_mode == in_core) then
  H_diag_RecNo = RecNo(1,1)
  memory_vectors(:,H_diag_RecNo) = H_diag(:)
end if

! the diagonalization must be run out of core:
! save H_diag on disk
if (save_mode == on_disk) then
  H_diag_RecNo = RecNo(1,1)
  iDisk = disk_address(H_diag_RecNo)
  call DDaFile(LuDavid,1,H_diag,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save H_diag
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  KeyWord = 'H_diag'
  call page_out(KeyWord,nConf,H_diag,LuDavid)
end if

call Timing(Time(2),dum1,dum2,dum3)
TimePage = TimePage+Time(2)-Time(1)

return

end subroutine Save_H_diag
