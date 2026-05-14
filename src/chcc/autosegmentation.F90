!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine autoSegmentation(Nprocs,maxspace,Jal1,Jal2,NvGrp,NvSGrp,NchBlk,wrksize,maxdim)
! Issues:
!
! 1) o2v4: n'(n'+1)/2 + n' "poltaskov"
! 2) odhadnut overhead? on-the-fly vs precalculate?

use chcc_global, only: maxGrp, maxSGrp, printkey
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Nprocs, maxspace, NchBlk
integer(kind=iwp), intent(out) :: Jal1, Jal2, NvGrp, NvSGrp, wrksize, maxdim
real(kind=wp) :: eff
real(kind=wp), parameter :: eff_thrs = 80.0_wp
logical(kind=iwp), parameter :: requireEfficiency = .true.

write(u6,*)
write(u6,*) '==============================='
write(u6,*) 'Autogenerating segmentation'
write(u6,'(A,i8)') ' Nprocs: ',Nprocs
write(u6,*)

! reset small segmentation
NvSGrp = 1

! get rough estimated of N', e.g. Nprocs == nntasks
NvGrp = int(sqrt(Two*Nprocs))
if (printkey >= 2) write(u6,'(A,i4)') ' 1st Np estimate: ',NvGrp

! if no such Np, get the smallest Np for which nntasks > Nprocs
if ((NvGrp*NvGrp/2) < Nprocs) then
  NvGrp = NvGrp+1
  if (printkey >= 2) write(u6,'(A,i4)') ' Corrected 1st Np estimate: ',NvGrp
end if

! make correction for efficiency, if required
if (requireEfficiency) then
  call findNextEffSeg(NvGrp,eff,Nprocs,eff_thrs,maxGrp,printkey)
  if (printkey >= 2) write(u6,'(A,i4,A,f6.2)') ' (maybe) further correction for efficiency: ',NvGrp,', efficiency: ',eff*100
end if

write(u6,*)

! does the segmentation fit into memory (start with Npp = 1)?
call checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)

do while (wrksize > maxspace)

  if (printkey >= 10) write(u6,'(A,i13)') ' Not enough memory. Max: ',maxspace,', Current: ',wrksize

  ! increase small segmentation, if possible

  if ((NvSGrp < 8) .and. ((NvGrp*(NvSGrp+1)) <= maxSGrp)) then
    NvSGrp = NvSGrp+1
    if (printkey >= 10) write(u6,'(A,i4)') ' Npp increased: ',NvSGrp

  else
    ! if not, increase large segmentation and reset small

    ! reset small segmentation in any case
    NvSGrp = 1

    ! increment large segmentation
    NvGrp = NvGrp+1
    if (printkey >= 10) write(u6,'(A,i4)') ' Np increased: ',NvGrp

    ! make correction for efficiency, if required
    if (requireEfficiency) then
      call findNextEffSeg(NvGrp,eff,Nprocs,eff_thrs,maxGrp,printkey)
      if (printkey >= 10) write(u6,'(A,i4)') ' Np increased (corrected for efficiency): ',NvGrp
    end if

    ! check new large segmentation

    if (NvGrp > maxGrp) then
      write(u6,*)
      write(u6,*) ' No suitable segmentation found, quitting'
      call abend()
    end if

  end if ! increase small or large

  ! print final results and
  ! recompute memory requirements
  if (printkey >= 10) write(u6,'(2(A,i4))') ' Current Np: ',NvGrp,', Npp: ',NvSgrp

  call checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)

end do ! we fit to memory

write(u6,*)
write(u6,'(A,2(i4),i8)') ' Final segmentation (Large/Small/Cholesky): ',NvGrp,NvSGrp,NchBlk
write(u6,*) '==============================='
write(u6,*)

return

end subroutine autoSegmentation
