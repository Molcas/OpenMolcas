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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Mar. 13, 2020, created this file.               *
! ****************************************************************

subroutine RotState()

use rasscf_global, only: ICMSP, ITER, IXMSP, LROOTS, IADR15, Ener
use PrintLevel, only: DEBUG, USUAL
use output_ras, only: LF, IPRLOC
use general_data, only: JOBIPH, NCONF
use Molcas, only: MxRoot
use RASDim, only: MxIter
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer NHrot     ! storing info in H0_Rotate.txt
integer NRState   ! storing info in Do_Rotate.txt
integer rcidisk
integer JRoot, IPRLEV
character(len=18) :: MatInfo
real*8, allocatable :: CIVEC(:,:), CIScr(:,:), HScr(:), State(:), HRot(:,:)
integer i, iad15

IPRLEV = IPRLOC(3)

if (IPRLEV >= USUAL) then
  write(LF,*)
  write(LF,*) repeat('=',71)
  write(LF,*)
  write(LF,'(11X,A)') 'Do_Rotate.txt is found in scratch directory.'
  if (IXMSP == 1) then
    write(LF,'(11X,A)') 'Following properties are for XMS intermediate states.'
  else if (ICMSP == 1) then
    write(LF,'(11X,A)') 'Following properties are for CMS intermediate states.'
  else
    write(LF,'(11X,A)') 'Following properties are for intermediate states'
    write(LF,'(11X,A)') ' obtained from the user-supplied rotation matrix'
  end if
end if

NRState = lRoots**2
NHRot = NRState

call mma_allocate(CIVec,nConf,lRoots,Label='CIVec')
call mma_allocate(CIScr,nConf,lRoots,Label='CIScr')
call mma_allocate(HScr,NHRot,Label='HScr')
call mma_allocate(State,NRState,Label='State')
call mma_allocate(HRot,lRoots,lRoots,Label='HRot')

!JB read rotation matrix in Do_Rotate.txt
call ReadMat2('ROT_VEC',MatInfo,State,lRoots,lRoots,7,18,'T')
if (IPRLEV >= DEBUG) then
  write(LF,*) 'rotation matrix'
  call RecPrt(' ',' ',State,lRoots,lRoots)
end if
HRot(:,:) = 0.0d0
NHRot = lRoots**2
do I=1,lRoots
  HRot(I,I) = ENER(I,ITER)
end do
call DGEMM_('t','n',lRoots,lRoots,lRoots,1.0d0,State,lRoots,HRot,lRoots,0.0d0,HScr,lRoots)
call DGEMM_('n','n',lRoots,lRoots,lRoots,1.0d0,HScr,lRoots,State,lRoots,0.0d0,HRot,lRoots)
call PrintMat2('ROT_HAM',MatInfo,HRot,lRoots,lRoots,7,18,'T')
if (IPRLEV >= DEBUG) then
  write(LF,'(6X,A)') 'Rotated Hamiltonian matrix '
  call RecPrt('HRot',' ',hRot,lRoots,lRoots)
end if

!JB read CI vector from jobiph
rcidisk = IADR15(4)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,CIScr(:,jRoot),nConf,rcidisk)
end do
call DGEMM_('n','n',NConf,lRoots,lRoots,1.0d0,CIScr,nConf,State,lRoots,0.0d0,CIVec,nConf)

! updating final energies as those for rotated states
rcidisk = IADR15(4)
do I=1,lRoots
  call DDafile(JOBIPH,1,CIVec(:,I),nConf,rcidisk)
  ENER(I,ITER) = HRot(I,I)
end do
IAD15 = IADR15(6)
call DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)

if (IPRLEV >= DEBUG) then
  write(LF,'(A)') 'Printing the coeff of the first CSF for each state'
  do I=1,lRoots
    write(LF,*) CIVec(1,I)
  end do
end if

call mma_deallocate(HScr)
call mma_deallocate(CIScr)
call mma_deallocate(State)
call mma_deallocate(CIVec)
call mma_deallocate(HRot)

if (IPRLEV >= USUAL) then
  write(LF,*)
  write(LF,*) repeat('=',71)
end if

end subroutine RotState
