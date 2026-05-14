!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine WRTTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,ISC)
! Print a batch of TTS blocks as given by IBLOCK
!
! ISC = 1 : In slater determinant form
! ISC = 2 : In Combination        form
!
! Jeppe Olsen, August 1995

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: BLOCKS(*)
integer(kind=iwp), intent(in) :: NBLOCK, IBLOCK(8,NBLOCK), NSMST, NSASO(NSMST,*), NSBSO(NSMST,*), ISC
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IOFF, IPACK, JBLOCK, NELMNT, NIA, NIB

write(u6,*) ' Batch of blocks'
write(u6,*) ' ==============='
write(u6,*)
write(u6,'(A,I4)') ' Number of blocks in batch ',NBLOCK

do JBLOCK=1,NBLOCK

  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  if (IBLOCK(1,JBLOCK) > 0) then

    if (ISC == 1) then
      IOFF = IBLOCK(5,JBLOCK)
    else
      IOFF = IBLOCK(6,JBLOCK)
    end if

    ! Is this block diagonal
    if ((ISC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP)) then
      IPACK = 1
    else
      IPACK = 0
    end if
    NIA = NSASO(IASM,IATP)
    NIB = NSBSO(IBSM,IBTP)
    !write(u6,*) ' iatp ibtp iasm ibsm nia nib ',iatp,ibtp,iasm,ibsm,nia,nib

    if (IPACK == 1) then
      NELMNT = nTri_Elem(NIA)
      if (NELMNT /= 0) then
        write(u6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
        write(u6,'(A)') '  ============================'
        call PRSM2(BLOCKS(IOFF),NIA)
      end if
    else
      NELMNT = NIA*NIB
      if (NELMNT /= 0) then
        write(u6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
        write(u6,'(A)') '  ============================'
        call WRTMAT(BLOCKS(IOFF),NIA,NIB,NIA,NIB)
      end if
    end if

  end if
end do

end subroutine WRTTTS
