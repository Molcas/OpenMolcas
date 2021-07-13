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

implicit integer(A-Z)
real*8 Sig_vec(nConf)
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"
character*16 KeyWord

call Timing(WTC_1,Swatch,Swatch,Swatch)

! check input arguments
if (nConf < 0) then
  write(6,*) 'Save_tmp_Sig_vec: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(6,*) 'Save_tmp_Sig_vec: iRoot less than 0'
  write(6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > n_Roots) then
  write(6,*) 'Save_tmp_Sig_vec: iRoot greater than nRoots'
  write(6,*) 'iRoot, nRoots = ',iRoot,n_Roots
  call Abend()
end if

! the diagonalization can be run in core:
! copy the sigma vector to new memory location
if (save_mode == in_core) then
  tmp_Sig_vec_RecNo = RecNo((5),iRoot)
  iMem = memory_address(tmp_Sig_vec_RecNo)
  call dCopy_(nConf,Sig_vec,1,Work(iMem),1)
end if

! the diagonalization must be run out of core:
! save the sigma vector on disk
if (save_mode == on_disk) then
  tmp_Sig_vec_RecNo = RecNo((5),iRoot)
  iDisk = disk_address(tmp_Sig_vec_RecNo)
  call DDaFile(LuDavid,1,Sig_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the sigma vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  KeyWord = ''
  write(KeyWord,'(A,I4.4)') 'tmp_Sig_vec',iRoot
  call page_out(KeyWord,nConf,Sig_vec,LuDavid)
end if

call Timing(WTC_2,Swatch,Swatch,Swatch)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Save_tmp_Sig_vec
