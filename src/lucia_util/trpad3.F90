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
use Definitions, only: u6

implicit real*8(A-H,O-Z)
real*8 MAT(NDIM,NDIM)

FAC2 = One-FACTOR**2

!IWAY = 1
IWAY = 2
if (IWAY == 1) then

  ! No blocking

  ! Lower half
  do J=1,NDIM
    do I=J,NDIM
      MAT(I,J) = MAT(I,J)+FACTOR*MAT(J,I)
    end do
  end do
  ! Upper half
  if (abs(FACTOR) /= One) then
    FAC2 = One-FACTOR**2
    do I=1,NDIM
      do J=1,I-1
        MAT(J,I) = FACTOR*MAT(I,J)+FAC2*MAT(J,I)
      end do
    end do
  else if (FACTOR == One) then
    do I=1,NDIM
      do J=1,I-1
        MAT(J,I) = MAT(I,J)
      end do
    end do
  else if (FACTOR == -One) then
    do I=1,NDIM
      do J=1,I-1
        MAT(J,I) = -MAT(I,J)
      end do
    end do
  end if
else if (IWAY == 2) then
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
        do J=JOFF,JEND
          MAT(I,J) = MAT(I,J)+FACTOR*MAT(J,I)
        end do
      end do
      ! Upper half
      if (abs(FACTOR) /= One) then
        FAC2 = One-FACTOR**2
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I
          do J=JOFF,JEND
            MAT(J,I) = FACTOR*MAT(I,J)+FAC2*MAT(J,I)
          end do
        end do
      else if (FACTOR == One) then
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I-1
          do J=JOFF,JEND
            MAT(J,I) = MAT(I,J)
          end do
        end do
      else if (FACTOR == -One) then
        do I=IOFF,IEND
          if (IBLK == JBLK) JEND = I
          do J=JOFF,JEND
            MAT(J,I) = -MAT(I,J)
          end do
        end do
      end if
      ! End of loop over blocks
    end do
  end do
  ! End of IWAY branching
end if

end subroutine TRPAD3
