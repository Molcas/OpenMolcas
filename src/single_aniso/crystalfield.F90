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

subroutine CRYSTALFIELD(ESOJ,DIPSO,S_SO,nDIMcf,d,nlanth,zmagn2,iopt,GRAD,iprint)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDIMcf, d, nlanth, iopt, iprint
real(kind=wp), intent(in) :: ESOJ(nDIMcf), ZMAGN2(3,3)
complex(kind=wp), intent(in) :: DIPSO(3,nDIMcf,nDIMcf), S_SO(3,nDIMcf,nDIMcf)
logical(kind=iwp), intent(in) :: GRAD
integer(kind=iwp) :: i, info, j
real(kind=wp) :: gtens(3), zmagn(3,3)
real(kind=wp), allocatable :: wtmp(:)
complex(kind=wp), allocatable :: DIPJ(:,:,:), SJ(:,:,:), TMP(:,:,:), ztmp(:,:)

write(u6,'(/)')
write(u6,'(A)') repeat('%',95)
if (mod(nDIMcf,2) == 1) then
  write(u6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J = ',(nDIMcf-1)/2,'.'
else
  write(u6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J = ',(nDIMcf-1),'/2.'
end if
write(u6,'(A)') repeat('%',95)
write(u6,*)

if (iopt == 1) then
  ! coordinate system for decomposition of the CF matrix identic to the coordinate system
  ! of the main magnetic axes of the ground multiplet (NDIM(1))
  call mma_allocate(TMP,3,d,d,label='TMP')
  TMP(:,:,:) = DIPSO(:,1:d,1:d)
  call atens(TMP,d,GTENS,ZMAGN,1)
  call mma_deallocate(TMP)
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(d,2) == 0) then
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudospin S = |',d-1,'/2> multiplet.'
  else
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudospin S = |',(d-1)/2,'> multiplet.'
  end if

else if (iopt == 2) then
  ! coordinate system for decomposition of the CF matrix identic to the coordinate system
  ! of the main magnetic axes of the ground multiplet (NDIM(1))
  call atens(DIPSO,nDIMCF,GTENS,ZMAGN,1)
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(nDIMCF,2) == 0) then
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic J = |',nDIMCF-1,'/2> multiplet'
  else
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic J = |',(nDIMCF-1)/2,'> multiplet'
  end if

else if (iopt == 3) then
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  write(u6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
  zmagn(:,:) = zmagn2(:,:)
else
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the initial coordinate system.'
  call unitmat(zmagn,3)
end if ! axisoption

! rotate the momentum:
call mma_allocate(DIPJ,3,nDIMcf,nDIMcf,'DIPJ')
call mma_allocate(SJ,3,nDIMcf,nDIMcf,'SJ')
call rotmom2(DIPSO,nDIMCF,ZMAGN,DIPJ)
call rotmom2(S_SO,nDIMCF,ZMAGN,SJ)

write(u6,'(a)') 'Rotation matrix from the initial coordinate system to the employed coordinate system is:'

if ((iopt == 1) .or. (iopt == 2)) then
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(u6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
  write(u6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
  write(u6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(ZMAGN(j,2),j=1,3),'|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,I3)') 'Quantization axis is Zm.'

else if (iopt == 3) then

  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(u6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system defined in the input','|'
  write(u6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
  write(u6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(ZMAGN(j,2),j=1,3),'|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,I3)') 'Quantization axis is Zm.'

else

  write(u6,'(A)') 'IDENTITY matrix.'
  write(u6,'(A,I3)') 'Quantization axis is the initial z axis.'
end if

if (IPRINT > 2) then
  call prMom('CRYSTALFIELD::   DIPJ(l,i,j)',DIPJ,nDIMcf)
  call prMom('CRYSTALFIELD::     SJ(l,i,j)',SJ,nDIMcf)
end if

call mma_allocate(ztmp,nDIMcf,nDIMcf,'z')
call mma_allocate(wtmp,nDIMcf,'w')
call diag_c2(DIPJ,nDIMcf,info,wtmp,ztmp)
do i=1,nDIMcf
  write(u6,'(A,i2,A,4F20.15)') 'energy: ',i,' : ',wtmp(i),wtmp(i)+wtmp(nDIMcf-i+1)
end do
call mma_deallocate(ztmp)
call mma_deallocate(wtmp)

call CRYSTALFIELD_1(nDIMcf,nlanth,DIPJ,ESOJ,GRAD,iprint)

call mma_deallocate(DIPJ)
call mma_deallocate(SJ)

return

end subroutine CRYSTALFIELD
