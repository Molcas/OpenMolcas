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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKSMAT()
! Set up S matrices for cases 1..13.

use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG, VERBOSE
use caspt2_global, only: DREF, LUSOLV, PREF
use caspt2_module, only: NASHT, NG3
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, byte

implicit none
integer(kind=iwp) :: iLUID
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: G3(:)

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,*) ' Construct S matrices'
end if

if (NASHT > 0) then
  !SVC: print header for debug info
  if (IPRGLB >= DEBUG) then
    write(u6,'("DEBUG> ",A)') 'CASE SYM S-MATRIX NORM'
    write(u6,'("DEBUG> ",A)') '==== === ============='
  end if
  ! For the cases A and C, begin by reading in the local storage
  !  part of the three-electron density matrix G3:
  call mma_allocate(G3,NG3,Label='G3')
  call PT2_GET(NG3,'GAMMA3',G3)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  call MKSA(DREF,size(DREF),PREF,size(PREF),NG3,G3,idxG3)
  call MKSC(DREF,size(DREF),PREF,size(PREF),NG3,G3,idxG3)

  call mma_deallocate(G3)
  call mma_deallocate(idxG3)

  !-SVC20100902: For the remaining cases that do not need G3, use replicate arrays
  call MKSB(DREF,size(DREF),PREF,size(PREF))
  call MKSD(DREF,size(DREF),PREF,size(PREF))
  call MKSE(DREF,size(DREF))
  call MKSF(PREF,size(PREF))
  call MKSG(DREF,size(DREF))
end if

call MKSH()

end subroutine MKSMAT
