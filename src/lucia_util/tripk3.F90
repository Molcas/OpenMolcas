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

subroutine TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SGN)
! REFORMATING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
!
! IWAY = 1 : FULL TO PACKED
!            LOWER HALF OF AUTPAK IS STORED IN APAK
! IWAY = 2 : PACKED TO FULL FORM
!            APAK STORED IN LOWER HALF
!             SGN * APAK TRANSPOSED IS STORED IN UPPPER PART
! NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
!
! Some considerations on cache minimization used for IMET = 2 Loop

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IWAY, MATDIM, NDIM
real(kind=wp), intent(inout) :: AUTPAK(MATDIM,MATDIM), APAK(*)
real(kind=wp), intent(in) :: SGN
integer(kind=iwp) :: I, IBLK, IEND, IJ, IJOFF, IMET, IOFF, IOFF2, J, JBLK, JEND, JOFF, LBLK, NBLK, NTEST

! To get rid of annoying and incorrect compiler warnings
IOFF = 0
JOFF = 0

! Packing : No problem with cache misses

if (IWAY == 1) then
  IJ = 0
  do J=1,NDIM
    APAK(IJ+J:IJ+NDIM) = AUTPAK(J:NDIM,J)
    IJ = IJ+NDIM-J
  end do
end if

! Unpacking : cache misses can occur so two routes

if (IWAY == 2) then
  ! Use blocked algorithm
  IMET = 2
  if (IMET == 1) then
    ! No blocking
    IJ = 0
    do J=1,NDIM
      do I=J,NDIM
        AUTPAK(J,I) = SGN*APAK(IJ+I)
        AUTPAK(I,J) = APAK(IJ+I)
      end do
      AUTPAK(J,J+1:NDIM) = SGN*APAK(IJ+J+1:IJ+NDIM)
      AUTPAK(J:NDIM,J) = APAK(IJ+J:IJ+NDIM)
      IJ = IJ+NDIM-J
    end do
  else if (IMET == 2) then
    ! Blocking
    LBLK = 40
    NBLK = MATDIM/LBLK
    if (LBLK*NBLK < MATDIM) NBLK = NBLK+1
    do JBLK=1,NBLK
      if (JBLK == 1) then
        JOFF = 1
      else
        JOFF = JOFF+LBLK
      end if
      JEND = min(JOFF+LBLK-1,MATDIM)
      do IBLK=JBLK,NBLK
        if (IBLK == JBLK) then
          IOFF = JOFF
        else
          IOFF = IOFF+LBLK
        end if
        IEND = min(IOFF+LBLK-1,MATDIM)
        do J=JOFF,JEND
          if (IBLK == JBLK) then
            IOFF2 = J
          else
            IOFF2 = IOFF
          end if
          IJOFF = (J-1)*MATDIM-J*(J-1)/2
          AUTPAK(J,IOFF2:IEND) = SGN*APAK(IJOFF+IOFF2:IJOFF+IEND)
          AUTPAK(IOFF2:IEND,J) = APAK(IJOFF+IOFF2:IJOFF+IEND)
        end do
        ! End of loop over I and J
      end do
    end do
    ! End of loop over blocks of I and J
  end if
end if

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' AUTPAK AND APAK FROM TRIPK3'
  call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
  call PRSM2(APAK,NDIM)
end if

end subroutine TRIPK3
