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

implicit none
#include "chcc1.fh"
integer Nprocs, NvGrp, NvSGrp, NchBlk
real*8 eff_thrs, eff
parameter(eff_thrs=80.0d0)
logical requireEfficiency
parameter(requireEfficiency=.true.)
integer Jal1, Jal2
integer wrksize, maxspace, maxdim

write(6,*)
write(6,*) '==============================='
write(6,*) 'Autogenerating segmentation'
write(6,'(A,i8)') ' Nprocs: ',Nprocs
write(6,*)

! reset small segmentation
NvSGrp = 1

! get rough estimated of N', e.g. Nprocs == nntasks
NvGrp = int(sqrt(2.0d0*Nprocs))
if (printkey >= 2) write(6,'(A,i4)') ' 1st Np estimate: ',NvGrp

! if no such Np, get the smallest Np for which nntasks > Nprocs
if ((NvGrp*NvGrp/2.0d0) < Nprocs) then
  NvGrp = NvGrp+1
  if (printkey >= 2) write(6,'(A,i4)') ' Corrected 1st Np estimate: ',NvGrp
end if

! make correction for efficiency, if required
if (requireEfficiency) then
  call findNextEffSeg(NvGrp,eff,Nprocs,eff_thrs,maxGrp,printkey)
  if (printkey >= 2) write(6,'(A,i4,A,f6.2)') ' (maybe) further correction for efficiency: ',NvGrp,', efficiency: ',eff*100
end if

write(6,*)

! does the segmentation fit into memory (start with Npp = 1)?
call checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)

12 continue
if (wrksize > maxspace) then

  if (printkey >= 10) write(6,'(A,i13)') ' Not enough memory. Max: ',maxspace,', Current: ',wrksize

  ! increase small segmentation, if possible

  if ((NvSGrp < 8) .and. ((NvGrp*(NvSGrp+1)) <= maxSGrp)) then
    NvSGrp = NvSGrp+1
    if (printkey >= 10) write(6,'(A,i4)') ' Npp increased: ',NvSGrp

  else
    ! if not, increase large segmentation and reset small

    ! reset small segmentation in any case
    NvSGrp = 1

    ! increment large segmentation
    NvGrp = NvGrp+1
    if (printkey >= 10) write(6,'(A,i4)') ' Np increased: ',NvGrp

    ! make correction for efficiency, if required
    if (requireEfficiency) then
      call findNextEffSeg(NvGrp,eff,Nprocs,eff_thrs,maxGrp,printkey)
      if (printkey >= 10) write(6,'(A,i4)') ' Np increased (corrected for efficiency): ',NvGrp
    end if

    ! check new large segmentation

    if (NvGrp > maxGrp) then
      write(6,*)
      write(6,*) ' No suitable segmentation found, quitting'
      call abend()
    end if

  end if ! increase small or large

  ! print final results and
  ! recompute memory requirements
  if (printkey >= 10) write(6,'(2(A,i4))') ' Current Np: ',NvGrp,', Npp: ',NvSgrp

  call checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)
  goto 12

end if ! we fit to memory

write(6,*)
write(6,'(A,2(i4),i8)') ' Final segmentation (Large/Small/Cholesky): ',NvGrp,NvSGrp,NchBlk
write(6,*) '==============================='
write(6,*)

return

end subroutine autoSegmentation
