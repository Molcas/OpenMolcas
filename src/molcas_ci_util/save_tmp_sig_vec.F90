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

subroutine Save_tmp_Sig_vec(iRoot,nConf,Sig_vec,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Save a temporary sigma vector                                    *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     iRoot   : integer                                                *
!               root number                                            *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     Sig_vec : array of real                                          *
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

use timers, only: TimePage
use davctl_mod, only: disk_address, in_core, llab, memory_vectors, mixed_mode_1, mixed_mode_2, n_Roots, on_disk, save_mode
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iRoot, nConf, LuDavid
real(kind=wp), intent(_IN_) :: Sig_vec(nConf)
integer(kind=iwp) :: iDisk, tmp_Sig_vec_RecNo
real(kind=wp) :: dum1, dum2, dum3, Time(2)
character(len=llab) :: KeyWord
integer(kind=iwp), external :: RecNo
#include "rasdim.fh"

call Timing(Time(1),dum1,dum2,dum3)

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Save_tmp_Sig_vec: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(u6,*) 'Save_tmp_Sig_vec: iRoot less than 0'
  write(u6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > n_Roots) then
  write(u6,*) 'Save_tmp_Sig_vec: iRoot greater than nRoots'
  write(u6,*) 'iRoot, nRoots = ',iRoot,n_Roots
  call Abend()
end if

! the diagonalization can be run in core:
! copy the sigma vector to new memory location
if (save_mode == in_core) then
  tmp_Sig_vec_RecNo = RecNo(5,iRoot)
  memory_vectors(:,tmp_Sig_vec_RecNo) = Sig_vec(:)
end if

! the diagonalization must be run out of core:
! save the sigma vector on disk
if (save_mode == on_disk) then
  tmp_Sig_vec_RecNo = RecNo(5,iRoot)
  iDisk = disk_address(tmp_Sig_vec_RecNo)
  call DDaFile(LuDavid,1,Sig_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the sigma vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  write(KeyWord,'(A,I4.4)') 'tmp_Sig_vec',iRoot
  call page_out(KeyWord,nConf,Sig_vec,LuDavid)
end if

call Timing(Time(2),dum1,dum2,dum3)
TimePage = TimePage+Time(2)-Time(1)

return

end subroutine Save_tmp_Sig_vec
