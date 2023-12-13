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

subroutine page_out(KeyWord,nConf,Vector,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Save any vector for further use by the Davidson diagonalization  *
!     Labels identifying the vectors are kept in a stack and to        *
!     minimize a write through cache strategy is applied               *
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

use davctl_mod, only: disk_address, LblStk, llab, memory_vectors, mixed_mode_1, mixed_mode_2, mxDiskStk, mxMemStk, nDiskStk, &
                      nMemStk, save_in_memory, save_mode
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
character(len=llab), intent(in) :: KeyWord
integer(kind=iwp), intent(in) :: nConf, LuDavid
real(kind=wp), intent(_IN_) :: Vector(nConf)
integer(kind=iwp) :: iDisk, iStk, nStk
#include "rasdim.fh"

! check input arguments
if (nConf < 0) then
  write(u6,*) 'page_out: nConf less than 0'
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

! there is a matching record identifier:
! overwrite the current record
if (nStk /= 0) then
  if (nStk <= mxMemStk) then
    memory_vectors(:,nStk) = Vector(:)
  else
    iDisk = disk_address(nStk-mxMemStk)
    call DDaFile(LuDavid,1,Vector,nConf,iDisk)
  end if
end if

! there is no matching record identifier:
! create a new record
if (nStk == 0) then
  if (save_mode == mixed_mode_1) then
    if (KeyWord(1:6) == 'CI_vec') then
      if (save_in_memory) then
        nMemStk = nMemStk+1
        memory_vectors(:,nMemStk) = Vector(:)
        LblStk(nMemStk) = KeyWord
        if (nMemStk == mxMemStk) save_in_memory = .false.
      else
        nMemStk = nMemStk+1
        if (nMemStk > mxMemStk) nMemStk = 1
        nDiskStk = nDiskStk+1
        if (nDiskStk > mxDiskStk) nDiskStk = 1
        iDisk = disk_address(nDiskStk)
        call DDaFile(LuDavid,1,memory_vectors(:,nMemStk),nConf,iDisk)
        memory_vectors(:,nMemStk) = Vector(:)
        LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
        LblStk(nMemStk) = KeyWord
      end if
    else
      nDiskStk = nDiskStk+1
      if (nDiskStk > mxDiskStk) nDiskStk = 1
      iDisk = disk_address(nDiskStk)
      call DDaFile(LuDavid,1,Vector,nConf,iDisk)
      LblStk(mxMemStk+nDiskStk) = KeyWord
    end if
  end if
  if (save_mode == mixed_mode_2) then
    if (save_in_memory) then
      nMemStk = nMemStk+1
      memory_vectors(:,nMemStk) = Vector(:)
      LblStk(nMemStk) = KeyWord
      if (nMemStk == mxMemStk) save_in_memory = .false.
    else
      nMemStk = nMemStk+1
      if (nMemStk > mxMemStk) nMemStk = 1
      nDiskStk = nDiskStk+1
      if (nDiskStk > mxDiskStk) nDiskStk = 1
      iDisk = disk_address(nDiskStk)
      call DDaFile(LuDavid,1,memory_vectors(:,nMemStk),nConf,iDisk)
      memory_vectors(:,nMemStk) = Vector(:)
      LblStk(mxMemStk+nDiskStk) = LblStk(nMemStk)
      LblStk(nMemStk) = KeyWord
    end if
  end if
end if

return

end subroutine page_out
