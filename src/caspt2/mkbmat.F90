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

subroutine MKBMAT()
! Set up B matrices for cases 1..13.

use PrintLevel, only: DEBUG, VERBOSE
use caspt2_global, only: DREF, iPrGlb, LUSOLV, PREF
use caspt2_module, only: NASHT, NG1, NG2, NG3
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, byte

implicit none
integer(kind=iwp) :: iLUID
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: F1(:), F2(:), F3(:), FD(:), FP(:)

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,*) ' Construct B matrices'
end if

if (NASHT /= 0) then

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

  call mma_allocate(F3,NG3,Label='F3')
  call PT2_GET(NG3,'DELTA3',F3)

  if (IPRGLB >= DEBUG) then
    write(u6,'("DEBUG> ",A)') 'CASE SYM B-MATRIX NORM'
    write(u6,'("DEBUG> ",A)') '==== === ============='
  end if

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  call MKBA(DREF,size(DREF),PREF,size(PREF),FD,FP,NG3,F3,idxG3)
  call MKBC(DREF,size(DREF),PREF,size(PREF),FD,FP,NG3,F3,idxG3)

  call mma_deallocate(F3)
  call mma_deallocate(idxG3)

  call MKBB(DREF,size(DREF),PREF,size(PREF),FD,FP)
  call MKBD(DREF,size(DREF),PREF,size(PREF),FD,FP)
  call MKBE(DREF,size(DREF),FD)
  call MKBF(DREF,size(DREF),PREF,size(PREF),FP)
  call MKBG(DREF,size(DREF),FD)

  call mma_deallocate(FP)
  call mma_deallocate(FD)

end if

call MKBH()

end subroutine MKBMAT
