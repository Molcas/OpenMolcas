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

subroutine SPIN_PHASE_RASSI(IPGLOB,DIPSO2,GMAIN,nDIM,ZIN,ZOUT)
! The RASSI program gives a random phase to the spin-orbit functions.

! This routine performs a simple check with the obtained spin functions,
! in order to determine the phase of the spin functions.
! IF the phase is not the same, then the spin functions will be multiplied
! with the corresponding coefficient that sets the same phase to all spin
! eigenfunctions

use spin_constants, only: Setup_Spin_Moment_Matrix, Spin
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IPGLOB, nDIM
complex(kind=wp) :: DIPSO2(3,nDIM,nDIM), ZIN(nDIM,nDIM), ZOUT(nDIM,nDIM)
real(kind=wp) :: GMAIN(3)
integer(kind=iwp) :: i, i2, j, l, ms1, ms2, NPAR
real(kind=wp), allocatable :: ALFA(:)
complex(kind=wp), allocatable :: PHS(:,:,:), Spin2(:,:,:), PHSA(:,:), PHSA2(:,:)

call Setup_Spin_Moment_Matrix()
! Determine the Parity:
NPAR = mod(nDIM,2)

call mma_allocate(PHS,3,nDIM,nDIM,Label='PHS')
call mma_allocate(Spin2,3,nDIM,nDIM,Label='Spin2')
call mma_allocate(PHSA,nDIM,nDIM,Label='PHSA')
call mma_allocate(PHSA2,nDIM,nDIM,Label='PHSA2')
call mma_allocate(ALFA,nDIM,Label='ALFA')

! Change the basis of the magnetic moment matrices from the RASSI functions to the
! effective spin eigenfunctions-- eigenfunctions of the Mu_Z

PHS(:,:,:) = cZero
SPIN2(:,:,:) = cZero
PHSA(:,:) = cZero
PHSA2(:,:) = cZero

do l=1,3
  do i=1,nDIM
    do j=1,nDIM
      do i2=1,nDIM
        PHS(l,i,j) = PHS(l,i,j)+sum(DIPSO2(l,:,i2)*conjg(ZIN(:,i))*ZIN(i2,j))
      end do
    end do
  end do
end do

! Rewrite the Spin m.e. in a new basis:

i = 0
do ms1=(nDIM-NPAR)/2,-(nDIM-NPAR)/2,-1
  if ((ms1 == 0) .and. (NPAR == 0)) cycle
  i = i+1
  j = 0
  do ms2=(nDIM-NPAR)/2,-(nDIM-NPAR)/2,-1
    if ((ms2 == 0) .and. (NPAR == 0)) cycle
    j = j+1
    Spin2(1,i,j) = Spin(1,nDIM,ms1,ms2)
    Spin2(2,i,j) = Spin(2,nDIM,ms1,ms2)
    Spin2(3,i,j) = Spin(3,nDIM,ms1,ms2)
  end do
end do

do i=1,nDIM
  do j=1,nDIM
    if (Spin2(1,i,j) == cZero) cycle
    PHSA2(i,j) = PHS(1,i,j)/(GMAIN(1)*Spin2(1,i,j))
    PHSA(i,j) = -log(PHSA2(i,j))*Onei
  end do
end do

ALFA(1) = Zero
do I=2,nDIM
  ALFA(I) = ALFA(I-1)+real(PHSA(I-1,I))
end do

do i=1,nDIM
  ZOUT(:,I) = exp(-ALFA(i)*Onei)*ZIN(:,i)
end do

if (IPGLOB > 2) then
  write(u6,'(/)')
  write(u6,'( 5x,a)') 'MAGNETIC MOMENT MATRIX ELEMENTS IN THE BASIS OF SPIN EIGENFUNCTIONS -- PHS(ic,i,j)'
  write(u6,*)
  do l=1,3
    write(u6,'(5x,a,i2)') 'PROJECTION=',l
    write(u6,*)
    do i=1,nDIM
      write(u6,'(16(2F12.8,2x))') (PHS(l,i,j),j=1,nDIM)
    end do
  end do
  write(u6,'(5x,a)') 'Spin2(1,i,j)'
  do i=1,nDIM
    write(u6,'(3x,16(2F12.6,2x))') (Spin2(1,i,j),j=1,nDIM)
  end do
  write(u6,'(5x,a)') 'Spin2(2,i,j)'
  do i=1,nDIM
    write(u6,'(3x,16(2F12.6,2x))') (Spin2(2,i,j),j=1,nDIM)
  end do
  write(u6,'(5x,a)') 'PHSA2(i,j)'
  do i=1,nDIM
    write(u6,'(3x,16(2F12.6,2x))') (PHSA2(i,j),j=1,nDIM)
  end do
  write(u6,'(5x,a)') 'PHSA(i,j)'
  do i=1,nDIM
    write(u6,'(3x,16(2F12.6,2x))') (PHSA(i,j),j=1,nDIM)
  end do
  write(u6,'(5x,a)') 'ALFA'
  write(u6,'(3x,16(2F12.6,2x))') (ALFA(j),j=1,nDIM)
  write(u6,'(5x,a)') 'ZOUT'
  do j=1,nDIM
    write(u6,'(3x,16(2F12.6,2x))') (ZOUT(i,j),i=1,nDIM)
  end do
end if

call mma_deallocate(PHS)
call mma_deallocate(Spin2)
call mma_deallocate(PHSA)
call mma_deallocate(PHSA2)
call mma_deallocate(ALFA)

end subroutine SPIN_PHASE_RASSI
