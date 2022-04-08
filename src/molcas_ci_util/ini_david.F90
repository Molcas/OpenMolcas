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

subroutine Ini_David(nRoots,nConf,nDet,nSel,n_keep,ntAsh,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Prepare address tables for further use by the Davidson           *
!     diagonalization scheme which is written in such a way that       *
!     the mechanism by which I/O is done is hidden in the subroutines  *
!     call load_xxx and save_xxx, where xxx stands may be one of       *
!     the following choices: H_diag, CI_vec, Sig_vec or Tmp_vec.       *
!     If possible (there is enough memory) a write through cache       *
!     mechanism is applied, that is to say all accessible memory is    *
!     used as a RAM-disk and dumped to physical disk in a FIFO mode.   *
!                                                                      *
!     calling arguments:                                               *
!     lRoots  : integer                                                *
!               number of roots to be optimized                        *
!     nConf   : integer                                                *
!               length of the CI vector in the CSF basis               *
!     nDet    : integer                                                *
!               length of the CI vector in the determinant basis       *
!     ntAsh   : integer                                                *
!               total number of active orbitals                        *
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

use davctl_mod, only: disk_address, in_core, istart, LblStk, memory_vectors, mixed_mode_1, mixed_mode_2, mxDiskStk, mxMemStk, &
                      n_Roots, nDiskStk, nkeep, nMemStk, nvec, on_disk, save_in_memory, save_mode
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRoots, nConf, nDet, nSel, n_keep, ntAsh, LuDavid
integer(kind=iwp) :: CI_vec_RecNo, H_diag_RecNo, iDisk, iRoot, lTmp1, lTmp2, lTmp3, Max_free_Mem, Max_used_Mem, Memory_Needed, &
                     nStk, Sig_vec_RecNo, tmp_CI_vec_RecNo, tmp_Sig_vec_RecNo
real(kind=wp) :: Dum(1)
integer(kind=iwp), external :: RecNo
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf_lucia.fh"

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Ini_David: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (nRoots < 0) then
  write(u6,*) 'Ini_David: nRoots less than zero'
  write(u6,*) 'nRoots = ',nRoots
  call Abend()
end if
if (nRoots > mxRoot) then
  write(u6,*) 'Ini_David: nRoots greater than mxRoot'
  write(u6,*) 'nRoots, mxRoot = ',nRoots,mxRoot
  call Abend()
end if
if (nDet < 0) then
  write(u6,*) 'Ini_David: nDet less than zero'
  write(u6,*) 'nDet = ',nDet
  call Abend()
end if
if (ntAsh < 0) then
  write(u6,*) 'Ini_David: ntAsh less than 0'
  write(u6,*) 'ntAsh = ',ntAsh
  call Abend()
end if
if (ntAsh > mxAct) then
  write(u6,*) 'Ini_David: ntAsh greater than mxAct'
  write(u6,*) 'ntAsh, mxAct = ',ntAsh,mxAct
  call Abend()
end if
n_Roots = nRoots
nkeep = n_keep
! If unitialized, determine a reasonable nkeep
if (nkeep == 0) then
  nkeep = (2*mxRoot)*nRoots
  nkeep = min(nkeep,400)
  nkeep = max(nkeep,3*nRoots)
  nkeep = min(nkeep,2*mxRoot)
end if

istart = 0
nvec = nkeep

! check the amount of available memory and decide which algorithm
! is to be used to save intermediate results
MxMemStk = 0
MxDiskStk = 0
call mma_maxDBLE(Max_free_Mem)
Max_free_Mem = Max_free_Mem-3*(nDet+4)
Max_free_Mem = Max_free_Mem-3*(nConf+4)
Max_free_Mem = Max_free_Mem-2*(ntAsh**3+4)
Max_free_Mem = Max_free_Mem-5*(ntAsh**2+4)
Max_used_Mem = (1+2*nKeep+2*nRoots)*(nConf+4)
! Calculate how much memory is needed in the rest of the Davidson
Memory_Needed = 0
if (ntAsh == 0) then
  Memory_Needed = 0
else if (nSel == nConf) then
  Memory_Needed = 2*nSel+nSel*nSel
else
  ! First: davctl
  Memory_Needed = 2*nSel+nSel*nSel
  ! Now: david5
  lTmp1 = nKeep
  lTmp2 = lTmp1*lTmp1
  lTmp3 = (lTmp2+lTmp1)/2
  Memory_Needed = Memory_Needed+5*nDet+lTmp1+3*lTmp2+2*lTmp3+3*nRoots*nSel
  ! Then: lucia_util
  Memory_Needed = Memory_Needed+Memory_Needed_Lucia
end if

if (Max_free_Mem < (nConf+4+Memory_Needed)) then
  MxMemStk = 0
  MxDiskStk = 1+2*nkeep+2*nRoots
  save_mode = on_disk
else if (Max_free_Mem >= Max_used_Mem+Memory_Needed) then
  MxMemStk = 1+2*nkeep+2*nRoots
  MxDiskStk = 0
  save_mode = in_core
else
  MxMemStk = Max_free_Mem/(nConf+4+Memory_Needed)
  MxDiskStk = 1+2*nkeep+2*nRoots-mxMemStk
  save_mode = mixed_mode_2
  if (mxMemStk < (nkeep+1)) save_mode = mixed_mode_1
end if
nMemStk = 0
nDiskStk = 0

call mma_allocate(disk_address,mxDiskStk,label='disk_address')
call mma_allocate(memory_vectors,nConf,mxMemStk,label='memory_vectors')

select case (save_mode)
  case (in_core)
    ! the diagonalization can be run in core,
    ! there's nothing to be done, everything is already allocated and indexed

  case (on_disk)
    ! the diagonalization must be run out of core:
    ! allocate disk space for all vectors that will be needed
    iDisk = 0
    H_diag_RecNo = RecNo(1,1)
    disk_address(H_diag_RecNo) = iDisk
    Dum(1) = Zero
    call DDafile(LuDavid,0,Dum,nConf,iDisk)
    do iRoot=1,nkeep
      CI_vec_RecNo = RecNo(2,iRoot)
      disk_address(CI_vec_RecNo) = iDisk
      call DDafile(LuDavid,0,Dum,nConf,iDisk)
    end do
    do iRoot=1,nKeep
      Sig_vec_RecNo = RecNo(3,iRoot)
      disk_address(Sig_vec_RecNo) = iDisk
      call DDaFile(LuDavid,0,Dum,nConf,iDisk)
    end do
    do iRoot=1,nRoots
      tmp_CI_vec_RecNo = RecNo(4,iRoot)
      disk_address(tmp_CI_vec_RecNo) = iDisk
      call DDaFile(LuDavid,0,Dum,nConf,iDisk)
    end do
    do iRoot=1,nRoots
      tmp_Sig_vec_RecNo = RecNo(5,iRoot)
      disk_address(tmp_Sig_vec_RecNo) = iDisk
      call DDaFile(LuDavid,0,Dum,nConf,iDisk)
    end do
  case (mixed_mode_1,mixed_mode_2)
    ! the diagonalization may be run in mixed mode:
    ! allocate memory and disk space for all vectors that will be needed
    iDisk = 0
    Dum(1) = Zero
    do nStk=1,mxDiskStk
      disk_address(nStk) = iDisk
      call DDaFile(LuDavid,0,Dum,nConf,iDisk)
    end do
    call mma_allocate(LblStk,mxDiskStk+mxMemStk,label='LblStk')
    LblStk(:) = ''
    save_in_memory = .true.
  case default
    call Abend()
end select

return

end subroutine Ini_David
