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

implicit integer(A-Z)
real*8 H_diag(nConf)
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"
character*16 KeyWord
call Timing(WTC_1,Swatch,Swatch,Swatch)

! check input arguments
if (nConf < 0) then
  write(6,*) 'Load_H_diag: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if

! the diagonalization can be run in core:
! copy H_diag to new memory location
if (save_mode == in_core) then
  H_diag_RecNo = RecNo((1),(1))
  iMem = memory_address(H_diag_RecNo)
  call dCopy_(nConf,Work(iMem),1,H_diag,1)
end if

! the diagonalization must be run out of core:
! load H_diag from disk
if (save_mode == on_disk) then
  H_diag_RecNo = RecNo((1),(1))
  iDisk = disk_address(H_diag_RecNo)
  call DDaFile(LuDavid,2,H_diag,nConf,iDisk)
end if

! the diagonalization may be run in mixed mode:
! use the write through cache mechanism to load and save H_diag
if ((save_mode == mixed_mode_1) .or. (save_mode == mixed_mode_2)) then
  KeyWord = ''
  write(KeyWord,'(A)') 'H_diag'
  call page_in(KeyWord,nConf,H_diag,LuDavid)
end if

call Timing(WTC_2,Swatch,Swatch,Swatch)
WTC_2 = WTC_2-WTC_1
WTC_3 = WTC_3+WTC_2

return

end subroutine Load_H_diag
