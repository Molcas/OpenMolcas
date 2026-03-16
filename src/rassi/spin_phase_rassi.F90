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

subroutine SPIN_PHASE_RASSI(IPGLOB,DIPSO2,GMAIN,DIM,ZIN,ZOUT)
! The RASSI program gives a random phase to the spin-orbit functions.

! This routine performs a simple check with the obtained spin functions,
! in order to determine the phase of the spin functions.
! IF the phase is not the same, then the spin functions will be multiplied
! with the corresponding coefficient that sets the same phase to all spin
! eigenfunctions

use spin_constants, only: Setup_Spin_Moment_Matrix, Spin
use Constants, only: Zero, cZero, Onei
use Definitions, only: u6

implicit none
integer l, i, j, i1, i2, NPAR, ms1, ms2, DIM, IPGLOB
real*8 GMAIN(3), ALFA(DIM)
complex*16 PHS(3,DIM,DIM), ZIN(DIM,DIM), DIPSO2(3,DIM,DIM), Spin2(3,DIM,DIM), PHSA(DIM,DIM), PHSA2(DIM,DIM), ZOUT(DIM,DIM)

call Setup_Spin_Moment_Matrix()
! Determine the Parity:
NPAR = mod(DIM,2)

! Change the basis of the magnetic moment matrices from the RASSI functions to the
! effective spin eigenfunctions-- eigenfunctions of the Mu_Z

PHS(:,:,:) = cZero
SPIN2(:,:,:) = cZero
PHSA(:,:) = cZero
PHSA2(:,:) = cZero

do l=1,3
  do i=1,DIM
    do j=1,DIM
      do i1=1,DIM
        do i2=1,DIM
          PHS(l,i,j) = PHS(l,i,j)+DIPSO2(l,i1,i2)*conjg(ZIN(i1,i))*ZIN(i2,j)
        end do
      end do
    end do
  end do
end do

! Rewrite the Spin m.e. in a new basis:

i = 0
do ms1=(DIM-NPAR)/2,-(DIM-NPAR)/2,-1
  if ((ms1 == 0) .and. (NPAR == 0)) goto 18
  i = i+1
  j = 0
  do ms2=(DIM-NPAR)/2,-(DIM-NPAR)/2,-1
    if ((ms2 == 0) .and. (NPAR == 0)) goto 17
    j = j+1
    Spin2(1,i,j) = Spin(1,DIM,ms1,ms2)
    Spin2(2,i,j) = Spin(2,DIM,ms1,ms2)
    Spin2(3,i,j) = Spin(3,DIM,ms1,ms2)
17  continue
  end do
18 continue
end do

do i=1,DIM
  do j=1,DIM
    if (Spin2(1,i,j) == cZero) goto 20
    PHSA2(i,j) = PHS(1,i,j)/(GMAIN(1)*Spin2(1,i,j))
    PHSA(i,j) = -log(PHSA2(i,j))*Onei
20  continue
  end do
end do

do I=1,DIM
  ALFA(I) = Zero
end do

do I=2,DIM
  ALFA(I) = ALFA(I-1)+real(PHSA(I-1,I))
end do

do i=1,DIM
  do J=1,DIM
    ZOUT(J,I) = cZero
    ZOUT(J,I) = exp(-ALFA(i)*Onei)*ZIN(j,i)
  end do
end do

if (IPGLOB > 2) then
  write(u6,'(/)')
  write(u6,'( 5x,a)') 'MAGNETIC MOMENT MATRIX ELEMENTS IN THE BASIS OF SPIN EIGENFUNCTIONS -- PHS(ic,i,j)'
  write(u6,*)
  do l=1,3
    write(u6,'(5x,a,i2)') 'PROJECTION=',l
    write(u6,*)
    do i=1,DIM
      write(u6,'(16(2F12.8,2x))') (PHS(l,i,j),j=1,DIM)
    end do
  end do
  write(u6,'(5x,a)') 'Spin2(1,i,j)'
  do i=1,DIM
    write(u6,'(3x,16(2F12.6,2x))') (Spin2(1,i,j),j=1,DIM)
  end do
  write(u6,'(5x,a)') 'Spin2(2,i,j)'
  do i=1,DIM
    write(u6,'(3x,16(2F12.6,2x))') (Spin2(2,i,j),j=1,DIM)
  end do
  write(u6,'(5x,a)') 'PHSA2(i,j)'
  do i=1,DIM
    write(u6,'(3x,16(2F12.6,2x))') (PHSA2(i,j),j=1,DIM)
  end do
  write(u6,'(5x,a)') 'PHSA(i,j)'
  do i=1,DIM
    write(u6,'(3x,16(2F12.6,2x))') (PHSA(i,j),j=1,DIM)
  end do
  write(u6,'(5x,a)') 'ALFA'
  write(u6,'(3x,16(2F12.6,2x))') (ALFA(j),j=1,DIM)
  write(u6,'(5x,a)') 'ZOUT'
  do j=1,DIM
    write(u6,'(3x,16(2F12.6,2x))') (ZOUT(i,j),i=1,DIM)
  end do
end if

end subroutine SPIN_PHASE_RASSI
