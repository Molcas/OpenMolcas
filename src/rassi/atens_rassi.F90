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

subroutine ATENS_RASSI(moment,ndim,gtens,maxes,IPGLOB)
!------------------------------------------------
!  ndim   -- size of the matrices
!  moment -- matrix of size (3,ndim,ndim) of the moment (magnetic, spin or angular)
!  gtens  -- array of size (3) keeping the main values of the A tensor ( dsqrt(main_values) )
!  maxes  -- array of size (3,3) keeping the main axes of the A tensor writen in
!            the right coordinate system (Determinant = +1)
!  IPGLOB -- the print level of the subroutine
!----------------------------------------------

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Six, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ndim, IPGLOB
complex(kind=wp) :: moment(3,ndim,ndim)
real(kind=wp) :: gtens(3), maxes(3,3)
integer(kind=iwp) :: i, ic1, ic2, info, j, k
real(kind=wp) :: A_TENS_TERM(3,3), Det_gtens, diff12, diff23, MAIN(3), W(3), Z(3,3), ZR(3,3)
complex(kind=wp) :: AC_TENS(3,3)
complex(kind=wp), allocatable :: A_TEMP(:,:,:,:)

! initialization:

call mma_allocate(A_TEMP,3,3,ndim,ndim,Label='A_TEMP')

AC_TENS(:,:) = cZero
A_TENS_TERM(:,:) = Zero
A_TEMP(:,:,:,:) = cZero

do j=1,ndim
  do i=1,ndim
    do ic2=1,3
      do ic1=1,3
        A_temp(ic1,ic2,i,j) = A_temp(ic1,ic2,i,j)+sum(moment(ic1,i,:)*moment(ic2,:,j))
      end do
    end do
  end do
end do

if (IPGLOB >= 4) then
  write(u6,'(/)')
  write(u6,'(5X,A)') 'BPMOMENT(ic1,ic2):'
  write(u6,*)
  do ic1=1,3
    do i=1,ndim
      do j=1,ndim
        write(u6,*) moment(ic1,i,j)
      end do
    end do
  end do

  write(u6,'(/)')
  write(u6,'(5X,A)') 'A_TEMP(ic1,ic2):'
  write(u6,*)
  do ic1=1,3
    do ic2=1,3
      do i=1,ndim
        do j=1,ndim
          write(u6,*) A_temp(ic1,ic2,i,j)
        end do
      end do
    end do
  end do
end if

do ic1=1,3
  do ic2=1,3
    do k=1,ndim
      Ac_tens(ic1,ic2) = Ac_tens(ic1,ic2)+A_temp(ic1,ic2,k,k)
    end do
  end do
end do
call mma_deallocate(A_TEMP)
do ic1=1,3
  A_TENS_TERM(ic1,:) = Six*real(Ac_tens(ic1,:)+Ac_tens(:,ic1))/real(ndim**3-ndim,kind=wp)
end do

!if (IPGLOB > 2) then
write(u6,'(/)')
write(u6,'(5X,A)') 'A_TENS_TERM(ic1,ic2):'
write(u6,*)
do ic1=1,3
  !write(u6,'(5X,3(2F14.7,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
  write(u6,'(5X,3(2F21.14,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
end do
!endif

! Diagonalization of A_tens - g tensors

W(:) = Zero
Z(:,:) = Zero
info = 0

call DIAG_R2_RASSI(A_TENS_TERM,3,info,w,z)
if (INFO /= 0) return

if (all(w(:) < Zero)) then
  write(u6,'(2x,A)') 'ALL EIGENVALUES OF THE A-TENSOR ARE NEGATIVE'
  write(u6,'(2X,A)') 'THIS IS A VERY UNUSUAL SITUATION. PLEASE CHECK MANUALLY'
  write(u6,'(2x,A)') 'THE FOLLOWING PART OF THE PSEUDOSPIN SECTION'
  write(u6,'(2x,A)') 'MUST BE DISREGARDED. THE RESULTS ARE NOT TRUSTABLE.'
  return
end if

if (IPGLOB >= 4) then
  write(u6,*)
  write(u6,'(4x,A)') 'A_TENS_TERM TENSOR:'
  write(u6,'(2a)') repeat('-',56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(4a,3x,a)') repeat('-',19),'|',repeat('-',36),'|','Xm, Ym, Zm -- main magnetic axes'
  write(u6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z','|'
  write(u6,'(6a)') repeat('-',19),'|',repeat('-',4),'|',repeat('-',31),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',w(1),' | Xm |',(z(j,1),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',w(2),' | Ym |',(z(j,2),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',w(3),' | Zm |',(z(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',56),'|'
end if

do I=1,3
  if (W(I) < Zero) W(I) = 1.0e-15_wp
  MAIN(i) = sqrt(W(i))
end do

if (IPGLOB >= 4) write(u6,'(5x,a,3F9.5)') 'EIGenValues after DSPEV:',(W(I),I=1,3)

! Check the sign of the coordinate system. if CS is Left-handed,
! then change it to RIGHT-handed
Det_gtens = Zero
ZR(:,:) = Z(:,:)
Det_gtens = ZR(1,1)*(ZR(2,2)*ZR(3,3)-ZR(2,3)*ZR(3,2))-ZR(1,2)*(ZR(2,1)*ZR(3,3)-ZR(2,3)*ZR(3,1))+ &
            ZR(1,3)*(ZR(2,1)*ZR(3,2)-ZR(2,2)*ZR(3,1))
if (Det_gtens < Zero) then
  Z(:,1) = -Z(:,1)
  if (IPGLOB > 2) write(u6,'(a)') 'The original coordinate system was LEFT-handed. It has been changed to the RIGHT-handed'
end if
diff12 = MAIN(2)-MAIN(1)
diff23 = MAIN(3)-MAIN(2)
if (IPGLOB >= 4) then
  write(u6,'(5x,a,3F19.15)') 'diff12 = ',diff12
  write(u6,'(5x,a,3F19.15)') 'diff23 = ',diff23
end if

gtens(:) = Zero
maxes(:,:) = Zero
! set the main Z axis:
if (diff12 > diff23) then
  gtens(3) = MAIN(1)
  gtens(2) = MAIN(2)
  gtens(1) = MAIN(3)
  if (Z(3,1) >= Zero) then
    maxes(:,3) = Z(:,1)
    maxes(:,1) = Z(:,3)
  else if (Z(3,1) < Zero) then
    maxes(:,3) = -Z(:,1)
    maxes(:,1) = Z(:,3)
  end if
else if (diff23 > diff12) then
  gtens(3) = MAIN(3)
  gtens(2) = MAIN(2)
  gtens(1) = MAIN(1)
  if (Z(3,3) >= Zero) then
    maxes(:,3) = Z(:,3)
    maxes(:,1) = Z(:,1)
  else if (Z(3,3) < Zero) then
    maxes(:,3) = -Z(:,3)
    maxes(:,1) = Z(:,1)
  end if
end if
maxes(1,2) = maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
maxes(2,2) = maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
maxes(3,2) = maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

if (IPGLOB > 2) then
  write(u6,*)
  write(u6,'(20X,A)') 'A-TENSOR:'
  write(u6,*)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  xx    xy    xz  |',(A_TENS_TERM(1,ic2),ic2=1,3)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  yx    yy    yz  |',(A_TENS_TERM(2,ic2),ic2=1,3)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  zx    zy    zz  |',(A_TENS_TERM(3,ic2),ic2=1,3)
end if

if (IPGLOB > 2) then
  write(u6,*)
  write(u6,'(4x,A)') 'g TENSOR:'
  write(u6,'(2a)') repeat('-',56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(4a,3x,a)') repeat('-',19),'|',repeat('-',36),'|','Xm, Ym, Zm -- main magnetic axes'
  write(u6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z','|'
  write(u6,'(6a)') repeat('-',19),'|',repeat('-',4),'|',repeat('-',31),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',56),'|'
  !call Add_Info('GTENS_MAIN',gtens,3,5)
end if

end subroutine ATENS_RASSI
