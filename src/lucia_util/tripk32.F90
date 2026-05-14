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

!#define _DEBUGPRINT_
subroutine TRIPK32(AUTPAK,APAK,MATDIM,NDIM,SGN)
! REFORMATING BETWEEN LOWER TRIANGULAR PACKING
! AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
!
! PACKED TO FULL FORM
! APAK STORED IN LOWER HALF
!  SGN * APAK TRANSPOSED IS STORED IN UPPPER PART
! NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
!
! Some considerations on cache minimization used for IMET = 2 Loop

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(in) :: APAK(*), SGN
integer(kind=iwp), intent(in) :: MATDIM, NDIM
real(kind=wp), intent(inout) :: AUTPAK(MATDIM,MATDIM)
integer(kind=iwp) :: I, IBLK, IEND, IJ, IJOFF, IMET, IOFF, IOFF2, J, JBLK, JEND, JOFF, LBLK, NBLK

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
        IJOFF = (J-1)*MATDIM-nTri_Elem(J-1)
        AUTPAK(J,IOFF2:IEND) = SGN*APAK(IJOFF+IOFF2:IJOFF+IEND)
        AUTPAK(IOFF2:IEND,J) = APAK(IJOFF+IOFF2:IJOFF+IEND)
      end do
      ! End of loop over I and J
    end do
  end do
  ! End of loop over blocks of I and J
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' AUTPAK AND APAK FROM TRIPK32'
call WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
call PRSM2(APAK,NDIM)
#endif

end subroutine TRIPK32
