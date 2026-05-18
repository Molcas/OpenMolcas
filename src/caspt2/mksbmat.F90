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

subroutine MKSBMAT()
! Set up S and B matrices for cases 1..13.

use definitions, only: iwp, wp, u6, Byte
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG, VERBOSE
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_global, only: DREF, PREF
use caspt2_global, only: LUSOLV
use caspt2_module, only: NASHT
use caspt2_module, only: NG1, NG2, NG3

implicit none
integer(kind=Byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: F1(:), F2(:), F3(:), FD(:), FP(:), G3(:)
integer(kind=iwp) iLUID
logical(kind=iwp) Single_set_of_PCO

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,*) ' Construct S and B matrices'
end if

! Single_set_of_PCO=.TRUE.
Single_set_of_PCO = .false.
if (Single_set_of_PCO) then
  call MKSMAT()
  call MKBMAT()
  return
end if

if (NASHT /= 0) then
  !SVC: print header for debug info
  if (IPRGLB >= DEBUG) then
    write(u6,'("DEBUG> ",A)') 'CASE SYM S/B-MATRIX NORM'
    write(u6,'("DEBUG> ",A)') '==== === ==============='
  end if
  ! For the cases A and C, begin by reading in the local storage
  !  part of the three-electron density matrix G3:

  !-SVC20100902: For the remaining cases that do not need G3, use replicate arrays

  call mma_allocate(F1,NG1,Label='F1')
  call PT2_GET(NG1,'DELTA1',F1)

  call mma_allocate(FD,size(DREF),Label='FD')
  call MKDREF_RPT2(NASHT,F1,FD,size(DREF))

  call mma_deallocate(F1)

  call mma_allocate(F2,NG2,Label='F2')
  call PT2_GET(NG2,'DELTA2',F2)

  call mma_allocate(FP,size(PREF),Label='FP')
  call MKPREF_RPT2(NASHT,F2,FP,size(PREF))

  call mma_deallocate(F2)

  call mma_allocate(G3,NG3,Label='G3')
  call PT2_GET(NG3,'GAMMA3',G3)

  call mma_allocate(F3,NG3,Label='F3')
  call PT2_GET(NG3,'DELTA3',F3)

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  call MKSA(DREF,size(DREF),PREF,size(PREF),NG3,G3,idxG3)
  call MKBA(DREF,size(DREF),PREF,size(PREF),FD,FP,NG3,F3,idxG3)

  call MKSC(DREF,size(DREF),PREF,size(PREF),NG3,G3,idxG3)
  call MKBC(DREF,size(DREF),PREF,size(PREF),FD,FP,NG3,F3,idxG3)

  call mma_deallocate(F3)
  call mma_deallocate(G3)
  call mma_deallocate(idxG3)

  call MKSB(DREF,size(DREF),PREF,size(PREF))
  call MKBB(DREF,size(DREF),PREF,size(PREF),FD,FP)

  call MKSD(DREF,size(DREF),PREF,size(PREF))
  call MKBD(DREF,size(DREF),PREF,size(PREF),FD,FP)

  call MKSE(DREF,size(DREF))
  call MKBE(DREF,size(DREF),FD)

  call MKSF(PREF,size(PREF))
  call MKBF(DREF,size(DREF),PREF,size(PREF),FP)

  call MKSG(DREF,size(DREF))
  call MKBG(DREF,size(DREF),FD)

  call mma_deallocate(FP)
  call mma_deallocate(FD)

end if

call MKSH()
call MKBH()

end subroutine MKSBMAT
