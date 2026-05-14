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

subroutine TRPAD3(MAT,FACTOR,NDIM)
! MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
!
! With some considerations of effective cache use for large matrices

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(inout) :: MAT(NDIM,NDIM)
real(kind=wp), intent(in) :: FACTOR
integer(kind=iwp) :: I, IBLK, IEND, IOFF, J, JBLK, JEND, JOFF, LBLK, NBLK
real(kind=wp) :: FAC2
integer(kind=iwp), parameter :: IMET = 2

FAC2 = One-FACTOR**2

if (IMET == 1) then

  ! No blocking

  ! Lower half
  do J=1,NDIM
    MAT(J:NDIM,J) = MAT(J:NDIM,J)+FACTOR*MAT(J,J:NDIM)
  end do
  ! Upper half
  if (abs(FACTOR) /= One) then
    FAC2 = One-FACTOR**2
    do I=1,NDIM
      MAT(1:I-1,I) = FACTOR*MAT(I,1:I-1)+FAC2*MAT(1:I-1,I)
    end do
  else if (FACTOR == One) then
    do I=1,NDIM
      MAT(1:I-1,I) = MAT(I,1:I-1)
    end do
  else if (FACTOR == -One) then
    do I=1,NDIM
      MAT(1:I-1,I) = -MAT(I,1:I-1)
    end do
  end if

else if (IMET == 2) then

  ! Simple blocking of matrix

  LBLK = 40
  NBLK = NDIM/LBLK
  if (NBLK*LBLK < NDIM) NBLK = NBLK+1
  IOFF = 1-LBLK
  !write(u6,*) 'NBLK ',nblk
  do IBLK=1,NBLK
    if (IBLK == -1) write(u6,*) 'IBLK = ',IBLK
    IOFF = IOFF+LBLK
    IEND = min(IOFF+LBLK-1,NDIM)
    JOFF = 1-LBLK
    do JBLK=1,IBLK
      JOFF = JOFF+LBLK
      JEND = min(JOFF+LBLK-1,NDIM)
      ! Lower half
      do I=IOFF,IEND
        if (IBLK == JBLK) JEND = I
        MAT(I,JOFF:JEND) = MAT(I,JOFF:JEND)+FACTOR*MAT(JOFF:JEND,I)
      end do
      ! Upper half
      if (abs(FACTOR) /= One) then
        FAC2 = One-FACTOR**2
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I
          MAT(JOFF:JEND,I) = FACTOR*MAT(I,JOFF:JEND)+FAC2*MAT(JOFF:JEND,I)
        end do
      else if (FACTOR == One) then
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I-1
          MAT(JOFF:JEND,I) = MAT(I,JOFF:JEND)
        end do
      else if (FACTOR == -One) then
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I
          MAT(JOFF:JEND,I) = -MAT(I,JOFF:JEND)
        end do
      end if
      ! End of loop over blocks
    end do
  end do

end if
! End of IMET branching

end subroutine TRPAD3
