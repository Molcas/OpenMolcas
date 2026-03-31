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
! Copyright (C) 1984,1989, Per Ake Malmqvist                           *
!***********************************************************************

subroutine PART(SXY,TRA1,TRA2)
! PURPOSE: SXY CONTAINS THE NONSECONDARY PART OF THE MO OVERLAP
! MATRIX. UPON RETURN, TRA1 AND TRA2 WILL CONTAIN THE COEFFICIENTS
! FOR SEQUENTIAL SINGLE-ORBITAL TRANSFORMATIONS (VI.2, MY IJQC ARTICLE)
! TO BIORTHONORMAL ORBITALS. SXY, TRA1 AND TRA2 ARE SYMMETRY-BLOCKED.
! ORIGINAL VERSION, MALMQUIST 84-04-04
! RASSCF VERSION,   MALMQUIST 89-11-15

use rasdef, only: NRS1, NRS2, NRS3
use Symmetry_Info, only: nIrrep
use rassi_data, only: NISH, NOSH, NSXY, NTRA
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: SXY(NSXY), TRA1(NTRA), TRA2(NTRA)
integer(kind=iwp) :: II, ISY, N, NBLOCK, NDIMEN, NOMAX, NSIZE(4)
integer(kind=iwp), allocatable :: ScrPiv(:)
real(kind=wp), allocatable :: ScrBuf(:), ScrMat(:)

NOMAX = 0
do ISY=1,nIrrep
  NOMAX = max(NOSH(ISY),NOMAX)
end do
call mma_allocate(SCRMAT,NOMAX*NOMAX,Label='ScrMat')
call mma_allocate(SCRPIV,2*NOMAX,Label='ScrPiv')
call mma_allocate(SCRBUF,NOMAX,Label='ScrBuf')
II = 1
do ISY=1,nIrrep
  NDIMEN = NOSH(ISY)
  if (NDIMEN == 0) cycle
  NBLOCK = 0
  N = NISH(ISY)
  if (N > 0) then
    NBLOCK = NBLOCK+1
    NSIZE(NBLOCK) = N
  end if
  N = NRS1(ISY)
  if (N > 0) then
    NBLOCK = NBLOCK+1
    NSIZE(NBLOCK) = N
  end if
  N = NRS2(ISY)
  if (N > 0) then
    NBLOCK = NBLOCK+1
    NSIZE(NBLOCK) = N
  end if
  N = NRS3(ISY)
  if (N > 0) then
    NBLOCK = NBLOCK+1
    NSIZE(NBLOCK) = N
  end if
  call PART1(NDIMEN,NBLOCK,NSIZE,SXY(II),TRA1(II),TRA2(II),SCRMAT,SCRPIV,SCRBUF)
  II = II+NDIMEN**2
end do
call mma_deallocate(SCRMAT)
call mma_deallocate(SCRPIV)
call mma_deallocate(SCRBUF)

end subroutine PART
