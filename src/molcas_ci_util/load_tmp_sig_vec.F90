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

subroutine Load_tmp_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load a temporary sigma vector                                    *
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

use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, n_Roots, on_disk, save_mode
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iRoot, nConf, LuDavid
real(kind=wp), intent(out) :: Sig_vec(nConf)
integer(kind=iwp) :: iDisk, tmp_Sig_vec_RecNo
real(kind=wp) :: dum1, dum2, dum3
character(len=llab) :: KeyWord
integer(kind=iwp), external :: RecNo
#include "rasdim.fh"
#include "timers.fh"

call Timing(WTC_1,dum1,dum2,dum3)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Load_tmp_Sig_vec: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(u6,*) 'Load_tmp_Sig_vec: iRoot less than 0'
  write(u6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > n_Roots) then
  write(u6,*) 'Load_tmp_Sig_vec: iRoot greater than nRoots'
  write(u6,*) 'iRoot = ',n_Roots
  call Abend()
end if

! the diagonalization can be run in core:
! copy the sigma vector to new memory location
if (save_mode == in_core) then
  tmp_Sig_vec_RecNo = RecNo(5,iRoot)
  Sig_vec(:) = memory_vectors(:,tmp_Sig_vec_RecNo)
end if

! the diagonalization must be run out of core:
! load the sigma vector from disk
if (save_mode == on_disk) then
  tmp_Sig_vec_RecNo = RecNo(5,iRoot)
  iDisk = disk_address(tmp_Sig_vec_RecNo)
  call DDaFile(LuDavid,2,Sig_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the sigma vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  write(KeyWord,'(A,I4.4)') 'tmp_Sig_vec',iRoot
  call page_in(KeyWord,nConf,Sig_vec,LuDavid)
end if

call Timing(WTC_2,dum1,dum2,dum3)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Load_tmp_Sig_vec
