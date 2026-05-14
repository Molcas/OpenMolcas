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

subroutine pa_pseudo(M,S,d,iopt,zfin,MF,SF,coord,iprint)
! d - dimension of the pseudospin
! M(3,d,d) -- magnetic moment in the initial basis
! S(3,d,d) -- spin moment in the initial basis
! iopt - option for choosing the quantization axis
!       = 1 : quantization axis is the main magnetic axis of the entire manifold (d)
!       = 2 : quantization axis is the main magnetic axis of the ground doublet (low-lying two states)
!       = 3 : quantization axis is provided by the user by means of coord
!       = 4 : quantization axis is the unit matrix, i.e. the original Z axis
!  coord(3,3): matrix specifying the rotation of initial coordinate system to the
!              coordinate system of used for determination of the quantization axis
!              by diagonalizing the Zeeman hamiltonian
! zfin - pseudospin eigenfunctions
! MF(3,d,d) -- magnetic moment in the pseudospin basis
! SF(3,d,d) -- spin moment in the pseudospin basis

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: d, iopt, iprint
complex(kind=wp), intent(in) :: M(3,d,d), S(3,d,d)
complex(kind=wp), intent(out) :: zfin(d,d), MF(3,d,d), SF(3,d,d)
real(kind=wp), intent(inout) :: coord(3,3)
integer(kind=iwp) :: i, i1, i2, info, j, k
real(kind=wp) :: coord2(3,3), det, gtens(3), maxes(3,3)
complex(kind=wp) :: MM(3,2,2)
real(kind=wp), allocatable :: w(:)
complex(kind=wp), allocatable :: dipso2(:,:,:), hzee(:,:), s_so2(:,:,:), z(:,:)
real(kind=wp), external :: FindDetR
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

call mma_allocate(w,d,label='w')
call mma_allocate(dipso2,3,d,d,label='dipso2')
call mma_allocate(hzee,d,d,label='hzee')
call mma_allocate(s_so2,3,d,d,label='s_so2')
call mma_allocate(z,d,d,label='z')

! set to zero important variables:
! choose the quantization axis:
if (iopt == 1) then
  call atens(M,d,gtens,maxes,2)

else if (iopt == 2) then
  MM(:,:,:) = M(:,1:2,1:2)
  call atens(MM,2,gtens,maxes,2)

else if (iopt == 3) then
  write(u6,'(A)') 'User provided the COORD to the PA_PSEUDo:'

else if (iopt == 4) then
  write(u6,'(A)') 'COORD is set to unity'
  call unitmat(coord,3)
else
  write(u6,'(A)') 'check the iopt parameter to PA_PSEUDo!!!'
  !call Abort()
end if
!ccccccccccccccccccccccccccccccccccc
write(u6,'(A)') 'input "coord" matrix to PA_PSEUDo'
do i=1,3
  write(u6,'(3F20.14)') (coord(j,i),j=1,3)
end do
! check If "coord" is empty:
coord2(:,:) = coord(:,:)
det = FindDetR(coord2,3)
write(u6,'(A, f20.13)') 'det = ',det
if (abs(det-One) < 1.0e-4_wp) then  ! 'coord is not empty'
  maxes(:,:) = coord(:,:)
else                                ! 'coord is empty'
  coord(:,:) = maxes(:,:)
end if
write(u6,'(A)') 'Employed  axes for diagonalization of the Zeeman Hamiltonian:'
do i=1,3
  write(u6,'(3F20.14)') (maxes(j,i),j=1,3)
end do

do i=1,d
  do j=1,d
    do k=1,3
      dipso2(:,i,j) = dipso2(:,i,j)+M(k,i,j)*(maxes(k,:)*cOne)
      s_so2(:,i,j) = s_so2(:,i,j)+S(k,i,j)*(maxes(k,:)*cOne)
    end do
  end do
end do

hzee(:,:) = hzee(:,:)-dipso2(3,:,:)
info = 0
w(:) = Zero
z(:,:) = cZero
call diag_c2(hzee,d,info,w,z)
if (info /= 0) then
  write(u6,'(5x,a)') 'diagonalization of the zeeman hamiltonian failed.'
else

  call spin_phase(dipso2,d,z,zfin)

  if (iprint > 2) then
    write(u6,'(5X,A)') 'MAIN VALUES OF THE ZEEMAN HAMILTONIAN:'
    write(u6,*)
    if (mod(d,2) == 1) then
      do I=1,d
        write(u6,'(3X,A,I3,A,F17.3)') '|',(d-1)/2+(1-I),'> = ',W(i)
      end do
    else
      do I=1,d
        write(u6,'(3X,A,I3,A,F17.3)') '|',(d-1)-2*(I-1),'/2 > = ',W(i)
      end do
    end if
    write(u6,*)
    write(u6,'(1X,A)') 'EIGENFUNCTIONS OF THE EFFECTIVE SPIN:'
    write(u6,*)
    if (mod(d,2) == 1) then
      do i=1,d
        write(u6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',(d-1)/2+(1-i),' > :',(zfin(j,i),j=1,d)
      end do
    else
      do i=1,d
        write(u6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',(d-1)-2*(i-1),'/2 > :',(zfin(j,i),j=1,d)
      end do
    end if
  end if ! printing of eigenfunctions with identical phase.

  MF(:,:,:) = Zero
  SF(:,:,:) = Zero
  do i=1,d
    do j=1,d
      do i1=1,d
        do i2=1,d
          MF(:,i,j) = MF(:,i,j)+dipso2(:,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
          SF(:,i,j) = SF(:,i,j)+s_so2(:,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
        end do
      end do
    end do
  end do

  if (IPRINT > 2) then
    call prMom('PseudoSpin basis:  MAGNETIC MOMENT: MF :',MF,d)
    call prMom('PseudoSpin basis:      SPIN MOMENT: SF :',SF,d)
  end if
end if

call mma_deallocate(w)
call mma_deallocate(dipso2)
call mma_deallocate(hzee)
call mma_deallocate(s_so2)
call mma_deallocate(z)

return

end subroutine pa_pseudo
