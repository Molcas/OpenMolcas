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

subroutine CRYSTALFIELD(ESOJ,DIPSO,S_SO,nDIMcf,iDIM,nlanth,zmagn2,iopt,GRAD,iprint)

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: iprint, iDIM
integer, intent(in) :: nDIMcf, iopt, nlanth
real(kind=8), intent(in) :: ESOJ(nDIMcf)
real(kind=8), intent(in) :: ZMAGN2(3,3)
complex(kind=8), intent(in) :: DIPSO(3,nDIMcf,nDIMcf)
complex(kind=8), intent(in) :: S_SO(3,nDIMcf,nDIMcf)
logical, intent(in) :: GRAD
! local variables
integer :: info
real(kind=8), allocatable :: wtmp(:)
complex(kind=8), allocatable :: DIPJ(:,:,:)
complex(kind=8), allocatable :: SJ(:,:,:), ztmp(:,:)
integer :: i, j
real(kind=8), allocatable :: gtens(:), zmagn(:,:)

write(6,'(/)')
write(6,'(100A)') ('%',i=1,95)
if (mod(nDIMcf,2) == 1) then
  write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J = ',(nDIMcf-1)/2,'.'
else
  write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J = ',(nDIMcf-1),'/2.'
end if
write(6,'(100A)') ('%',i=1,95)
write(6,*)

call mma_allocate(gtens,3,'gtens')
call mma_allocate(zmagn,3,3,'zmagn')
call dcopy_(3,[0.0_wp],0,gtens,1)
call dcopy_(3*3,[0.0_wp],0,zmagn,1)
if (iopt == 1) then
  ! coordinate system for decomposition of the CF matrix identic to the coordinate system
  ! of the main magnetic axes of the ground multiplet (NDIM(1))
  call atens(DIPSO(1:3,1:idim,1:idim),idim,GTENS,ZMAGN,1)
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(iDIM,2) == 0) then
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseuDospin S = |',iDIM-1,'/2> multiplet.'
  else
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseuDospin S = |',(iDIM-1)/2,'> multiplet.'
  end if

else if (iopt == 2) then
  ! coordinate system for decomposition of the CF matrix identic to the coordinate system
  ! of the main magnetic axes of the ground multiplet (NDIM(1))
  call atens(DIPSO,nDIMCF,GTENS,ZMAGN,1)
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(nDIMCF,2) == 0) then
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic J = |',nDIMCF-1,'/2> multiplet'
  else
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic J = |',(nDIMCF-1)/2,'> multiplet'
  end if

else if (iopt == 3) then
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  write(6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
  call dcopy_(3*3,zmagn2,1,zmagn,1)
else
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the initial coordinate system.'
  do i=1,3
    ZMAGN(i,i) = 1.0_wp
  end do
end if ! axisoption

! rotate the momentum:
call mma_allocate(DIPJ,3,nDIMcf,nDIMcf,'DIPJ')
call mma_allocate(SJ,3,nDIMcf,nDIMcf,'SJ')
call zcopy_(3*nDIMcf*nDIMcf,[(0.0_wp,0.0_wp)],0,DIPJ,1)
call zcopy_(3*nDIMcf*nDIMcf,[(0.0_wp,0.0_wp)],0,SJ,1)
call rotmom2(DIPSO,nDIMCF,ZMAGN,DIPJ)
call rotmom2(S_SO,nDIMCF,ZMAGN,SJ)

write(6,'(a)') 'Rotation matrix from the initial coordinate system to the employed coordinate system is:'

if ((iopt == 1) .or. (iopt == 2)) then
  write(6,'(70a)') ('-',i=1,67),'|'
  write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
  write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
  write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(ZMAGN(j,2),j=1,3),'|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
  write(6,'(83a)') ('-',i=1,67),'|'
  write(6,'(A,I3)') 'Quantization axis is Zm.'

else if (iopt == 3) then

  write(6,'(70a)') ('-',i=1,67),'|'
  write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system defined in the input','|'
  write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
  write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(ZMAGN(j,2),j=1,3),'|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
  write(6,'(83a)') ('-',i=1,67),'|'
  write(6,'(A,I3)') 'Quantization axis is Zm.'

else

  write(6,'(A)') 'IDENTITY matrix.'
  write(6,'(A,I3)') 'Quantization axis is the initial z axis.'
end if

if (IPRINT > 2) then
  call prMom('CRYSTALFIELD::   DIPJ(l,i,j)',DIPJ,nDIMcf)
  call prMom('CRYSTALFIELD::     SJ(l,i,j)',SJ,nDIMcf)
end if

call mma_allocate(ztmp,nDIMcf,nDIMcf,'z')
call mma_allocate(wtmp,nDIMcf,'w')
wtmp(:) = 0.0_wp
ztmp(:,:) = (0.0_wp,0.0_wp)
info = 0
call diag_c2(DIPJ,nDIMcf,info,wtmp,ztmp)
do i=1,nDIMcf
  write(6,'(A,i2,A,4F20.15)') 'energy: ',i,' : ',wtmp(i),wtmp(i)+wtmp(nDIMcf-i+1)
end do
call mma_deallocate(ztmp)
call mma_deallocate(wtmp)

call CRYSTALFIELD_1(nDIMcf,nlanth,DIPJ,ESOJ,GRAD,iprint)

call mma_deallocate(DIPJ)
call mma_deallocate(SJ)
call mma_deallocate(gtens)
call mma_deallocate(zmagn)

return

end subroutine CRYSTALFIELD
