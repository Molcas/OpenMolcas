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

subroutine Term_David(ICICH,iter,lRoots,nConf,Vector,JOBIPH,LuDavid,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Terminate the Davidson diagonalization                           *
!                                                                      *
!     calling arguments:                                               *
!     ICICH   : integer                                                *
!               switch enabling root selection                         *
!     JOBIPH  : integer                                                *
!               logical unit number of the JOBIPH file                 *
!     iDisk   : integer                                                *
!               disk address of the first CI vector on JOBIPH          *
!     iter    : integer                                                *
!               iteration count of the final result                    *
!     nConf   : integer                                                *
!               length of the CI vector in the CSF basis               *
!     Vector  : array of real*8                                        *
!               temporary vector of length nConf                       *
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

use davctl_mod, only: disk_address, LblStk, memory_vectors
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ICICH, iter, lRoots, nConf, JOBIPH, LuDavid
integer(kind=iwp), intent(inout) :: iDisk
real(kind=wp), intent(out) :: Vector(nConf)
integer(kind=iwp) :: iRoot
real(kind=wp), allocatable :: Ovlp1(:,:), Ovlp2(:,:)
#include "rasdim.fh"

! check input arguments
if (nConf < 0) then
  write(u6,*) 'Term_David: nConf less than 0'
  write(u6,*) 'nConf = ',nConf
  call Abend()
end if
if (iter < 0) then
  write(u6,*) 'Term_David: iter less than 0'
  write(u6,*) 'iter = ',iter
  call Abend()
end if
if (iter > mxCiIt) then
  write(u6,*) 'Term_David: iter greater than mxCiIt'
  write(u6,*) 'iter, mxCiIt = ',iter,mxCiIt
  call Abend()
end if

! Restore the final CI vectors and save them for further use.
! If the root selectioning option has been enabled calculate
! also the overlap elemtents with the test vectors
if (ICICH == 1) then
  call mma_allocate(Ovlp1,lRoots,lRoots,label='CIovlp1')
  call mma_allocate(Ovlp2,lRoots,lRoots,label='CIovlp2')
  Ovlp1(:,:) = Zero
  Ovlp2(:,:) = Zero
end if
do iRoot=1,lRoots
  call Load_tmp_CI_vec(iRoot,nConf,Vector,LuDavid)
  call DDaFile(JOBIPH,1,Vector,nConf,iDisk)
  if (ICICH == 1) then
    call CIovlp(iRoot,Ovlp1,Ovlp2,Vector)
  end if
end do

! If the root selectioning option has been enabled
! make a new choice of the current roots
if (ICICH == 1) then
  call CIselect(Ovlp1,Ovlp2)
  call mma_deallocate(Ovlp1)
  call mma_deallocate(Ovlp2)
end if

call mma_deallocate(disk_address)
call mma_deallocate(memory_vectors)
if (allocated(LblStk)) call mma_deallocate(LblStk)

return

end subroutine Term_David
