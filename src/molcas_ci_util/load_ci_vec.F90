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

subroutine Load_CI_vec(iRoot,nConf,CI_vec,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load a CI vector                                                 *
!     further use by the Davidson diagonalization scheme               *
!                                                                      *
!     calling arguments:                                               *
!     iRoot   : integer                                                *
!               root number                                            *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     CI_vec  : array of real*8                                        *
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

implicit integer(A-Z)
real*8 CI_vec(nConf)
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"
character*16 KeyWord

call Timing(WTC_1,Swatch,Swatch,Swatch)

! check input arguments
if (nConf < 0) then
  write(6,*) 'Load_CI_vec: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if
if (iRoot < 0) then
  write(6,*) 'Load_CI_vec: iRoot less than 0'
  write(6,*) 'iRoot = ',iRoot
  call Abend()
end if
if (iRoot > nkeep) then
  write(6,*) 'Load_CI_vec: iRoot greater than nkeep'
  write(6,*) 'iRoot, nkeep = ',iRoot,nkeep
  call Abend()
end if

! the diagonalization can be run in core:
! copy the CI vector to new memory location
if (save_mode == in_core) then
  CI_vec_RecNo = RecNo((2),iRoot)
  iMem = memory_address(CI_vec_RecNo)
  call dCopy_(nConf,Work(iMem),1,CI_vec,1)
end if

! the diagonalization must be run out of core:
! load the CI vector from disk
if (save_mode == on_disk) then
  CI_vec_RecNo = RecNo((2),iRoot)
  iDisk = disk_address(CI_vec_RecNo)
  call DDaFile(LuDavid,2,CI_vec,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save
! the CI vector
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  CI_vec_PageNo = PageNo(iRoot)
  KeyWord = ''
  write(KeyWord,'(A,I4.4)') 'CI_vec',CI_vec_PageNo
  call page_in(KeyWord,nConf,CI_vec,LuDavid)
end if

call Timing(WTC_2,Swatch,Swatch,Swatch)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Load_CI_vec
