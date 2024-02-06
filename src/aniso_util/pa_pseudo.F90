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

subroutine pa_pseudo(M,S,dim,iopt,zfin,MF,SF,coord,iprint)
! dim - dimension of the pseuDospin
! M(3,dim,dim) -- magnetic moment in the initial basis
! S(3,dim,dim) -- spin moment in the initial basis
! iopt - option for choosing the quantization axis
!       = 1 : quantization axis is the main magnetic axis of the entire manIfold (dim)
!       = 2 : quantization axis is the main magnetic axis of the ground Doublet (low-lying two states)
!       = 3 : quantization axis is provided by the user by means of coord
!       = 4 : quantization axis is the unit matrix, i.e. the original Z axis
!  coord(3,3): matrix specIfying the rotation of initial coordinate system to the
!              coordinate system of used for determination of the quantization axis
!              by diagonalizing the Zeeman hamiltonian
! zfin - pseuDospin eigenfunctions
! MF(3,dim,dim) -- magnetic moment in the pseuDospin basis
! SF(3,dim,dim) -- spin moment in the pseuDospin basis

use Constants, only: Zero, One, cOne
use Definitions, only: wp, u6

implicit none
integer :: dim, info, i, j, k, l, i1, i2, iopt
integer :: iprint
real(kind=8) :: gtens(3), w(dim), maxes(3,3), det, FindDetR
real(kind=8) :: coord(3,3), coord2(3,3)
complex(kind=8) :: M(3,dim,dim), S(3,dim,dim), z(dim,dim)
complex(kind=8) :: MF(3,dim,dim), SF(3,dim,dim), zfin(dim,dim)
complex(kind=8) :: dipso2(3,dim,dim)
complex(kind=8) :: s_so2(3,dim,dim)
complex(kind=8) :: hzee(dim,dim)
external FindDetR
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! set to zero important variables:
! choose the quantization axis:
if (iopt == 1) then
  gtens(:) = Zero
  maxes(:,:) = Zero
  call atens(M,dim,gtens,maxes,2)

else if (iopt == 2) then
  gtens(:) = Zero
  maxes(:,:) = Zero
  call atens(M(1:3,1:2,1:2),2,gtens,maxes,2)

else if (iopt == 3) then
  write(u6,'(A)') 'User provided the COORD to the PA_PSEUDo:'

else if (iopt == 4) then
  write(u6,'(A)') 'COORD is set to unity'
  coord(:,:) = Zero
  do i=1,3
    coord(i,i) = One
  end do
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
call rZeroMatrix(coord2,3)
do i=1,3
  do j=1,3
    coord2(i,j) = coord(i,j)
  end do
end do
det = FindDetR(coord2,3)
write(u6,'(A, f20.13)') 'det = ',det
if (abs(det-One) < 1.0e-4_wp) then  ! 'coord is not empty'
  call rZeroMatrix(maxes,3)
  do i=1,3
    do j=1,3
      maxes(i,j) = coord(i,j)
    end do
  end do
else                                   ! 'coord is empty'
  call rZeroMatrix(coord,3)
  do i=1,3
    do j=1,3
      coord(i,j) = maxes(i,j)
    end do
  end do
end if
write(u6,'(A)') 'Employed  axes for diagonalization of the Zeeman Hamiltonian:'
do i=1,3
  write(u6,'(3F20.14)') (maxes(j,i),j=1,3)
end do

call cZeroMoment(s_so2,dim)
call cZeroMoment(dipso2,dim)
do l=1,3
  do i=1,dim
    do j=1,dim
      do k=1,3
        dipso2(l,i,j) = dipso2(l,i,j)+M(k,i,j)*(maxes(k,l)*cOne)
        s_so2(l,i,j) = s_so2(l,i,j)+S(k,i,j)*(maxes(k,l)*cOne)
      end do
    end do
  end do
end do

call cZeroMatrix(hzee,dim)
do i=1,dim
  do j=1,dim
    hzee(i,j) = hzee(i,j)-cOne*dipso2(3,i,j)
  end do
end do
info = 0
call rZeroVector(w,dim)
call cZeroMatrix(z,dim)
call diag_c2(hzee,dim,info,w,z)
if (info /= 0) then
  write(u6,'(5x,a)') 'diagonalization of the zeeman hamiltonian failed.'
  go to 199
end if

call cZeroMatrix(zfin,dim)
call spin_phase(dipso2,dim,z,zfin)

if (iprint > 2) then
  write(u6,'(5X,A)') 'MAIN VALUES OF THE ZEEMAN HAMILTONIAN:'
  write(u6,*)
  if (mod(dim,2) == 1) then
    do I=1,dim
      write(u6,'(3X,A,I3,A,F17.3)') '|',(dim-1)/2+(1-I),'> = ',W(i)
    end do
  else
    do I=1,dim
      write(u6,'(3X,A,I3,A,F17.3)') '|',(dim-1)-2*(I-1),'/2 > = ',W(i)
    end do
  end if
  write(u6,*)
  write(u6,'(1X,A)') 'EIGENFUNCTIONS OF THE EFFECTIVE SPIN:'
  write(u6,*)
  if (mod(dim,2) == 1) then
    do i=1,dim
      write(u6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',(dim-1)/2+(1-i),' > :',(zfin(j,i),j=1,dim)
    end do
  else
    do i=1,dim
      write(u6,'(10x,a,i2,a,10x,20(2f16.12,2x))') 'eigenvector of |',(dim-1)-2*(i-1),'/2 > :',(zfin(j,i),j=1,dim)
    end do
  end if
end if ! printing of eigenfunctions with identical phase.

call cZeroMoment(MF,dim)
call cZeroMoment(SF,dim)
do l=1,3
  do i=1,dim
    do j=1,dim
      do i1=1,dim
        do i2=1,dim
          MF(l,i,j) = MF(l,i,j)+dipso2(l,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
          SF(l,i,j) = SF(l,i,j)+s_so2(l,i1,i2)*conjg(zfin(i1,i))*zfin(i2,j)
        end do
      end do
    end do
  end do
end do

if (IPRINT > 2) then
  call prMom('PseudoSpin basis:  MAGNETIC MOMENT: MF :',MF,dim)
  call prMom('PseudoSpin basis:      SPIN MOMENT: SF :',SF,dim)
end if

199 continue
return

end subroutine pa_pseudo
