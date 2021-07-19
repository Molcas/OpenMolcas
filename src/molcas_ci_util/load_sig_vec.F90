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

subroutine Load_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load a sigma vector                                              *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     iRoot   : integer                                                *
!               root number                                            *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     Sig_vec : array of real*8                                        *
!               sigma vector                                           *
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

use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, nkeep, on_disk, save_mode
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iRoot, nConf, LuDavid
real(kind=wp), intent(out) :: Sig_vec(nConf)
integer(kind=iwp) :: iDisk, Sig_vec_PageNo, Sig_vec_RecNo
character(len=llab) :: KeyWord
integer(kind=iwp), external :: PageNo, RecNo
#include "rasdim.fh"
#include "timers.fh"

call Timing(WTC_1,Swatch,Swatch,Swatch)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Load_Sig_vec: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(u6,*) 'Load_Sig_vec: iRoot less than 0'
  write(u6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > nkeep) then
  write(u6,*) 'Load_Sig_vec: iRoot greater than nkeep'
  write(u6,*) 'iRoot, nkeep = ',iRoot,nkeep
  call Abend()
end if

! the diagonalization can be run in core:
! copy the sigma vector to new memory location
if (save_mode == in_core) then
  Sig_vec_RecNo = RecNo(3,iRoot)
  Sig_vec(:) = memory_vectors(:,Sig_vec_RecNo)
end if

! the diagonalization must be run out of core:
! load the sigma vector from disk
if (save_mode == on_disk) then
  Sig_vec_RecNo = RecNo(3,iRoot)
  iDisk = disk_address(Sig_vec_RecNo)
  call DDaFile(LuDavid,2,Sig_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the sigma vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  Sig_vec_PageNo = PageNo(iRoot)
  write(KeyWord,'(A,I4.4)') 'Sig_vec',Sig_vec_PageNo
  call page_in(KeyWord,nConf,Sig_vec,LuDavid)
end if

call Timing(WTC_2,Swatch,Swatch,Swatch)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Load_Sig_vec
