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

subroutine Load_H_diag(nConf,H_diag,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load the diagonal approximation of the CI Hamiltonian for        *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     H_diag  : array of real*8                                        *
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

use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, on_disk, save_mode
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nConf, LuDavid
real(kind=wp), intent(out) :: H_diag(nConf)
integer(kind=iwp) :: H_diag_RecNo, iDisk
real(kind=wp) :: dum1, dum2, dum3
character(len=llab) :: KeyWord
integer(kind=iwp), external :: RecNo
#include "rasdim.fh"
#include "timers.fh"

call Timing(WTC_1,dum1,dum2,dum3)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Load_H_diag: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if

! the diagonalization can be run in core:
! copy H_diag to new memory location
if (save_mode == in_core) then
  H_diag_RecNo = RecNo(1,1)
  H_diag(:) = memory_vectors(:,H_diag_RecNo)
end if

! the diagonalization must be run out of core:
! load H_diag from disk
if (save_mode == on_disk) then
  H_diag_RecNo = RecNo(1,1)
  iDisk = disk_address(H_diag_RecNo)
  call DDaFile(LuDavid,2,H_diag,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save H_diag
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  KeyWord = 'H_diag'
  call page_in(KeyWord,nConf,H_diag,LuDavid)
end if

call Timing(WTC_2,dum1,dum2,dum3)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Load_H_diag
