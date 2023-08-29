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

subroutine CHO_RS2RS(IMAP,LMAP,IRS2,IRS3,IRED3,ISYM)
!
! Purpose: set up mapping between reduced sets stored at IRS2 and
!          IRS3 (IRED3 is the reduced set id of IRS3).
!
! WARNING: for IRED3 = 1, INDRED is reset here!!!!

use Cholesky, only: iiBstR, iiBstRSh, IndRed, nnBstR, nnBstRSh, nnShl
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LMAP, IRS2, IRS3, IRED3, ISYM
integer(kind=iwp), intent(out) :: IMAP(LMAP)
integer(kind=iwp) :: I, I1, I2, IAB, IAB1, IAB2, ISHLAB, JAB, KAB, LAB, LAB1, LAB2, LAST, MAB, N2, N3
character(len=*), parameter :: SECNAM = 'CHO_RS2RS'

! Check input.
! ------------

if ((IRS2 < 1) .or. (IRS2 > 3) .or. (IRS3 < 1) .or. (IRS3 > 3)) then
  call CHO_QUIT('Index error in '//SECNAM,104)
else if (LMAP < NNBSTR(ISYM,IRS2)) then
  call CHO_QUIT('Dimension error in '//SECNAM,104)
end if

! For IRED3 = 1, INDRED array addresses into shell pair. We hence
! need to reset it (as warned about above).
! ---------------------------------------------------------------

if (IRED3 == 1) then
  I1 = IIBSTR(ISYM,IRS3)+1
  I2 = I1+NNBSTR(ISYM,IRS3)-1
  do I=I1,I2
    IndRed(I,IRS3) = I
  end do
end if

! Set up mapping array.
! ---------------------

IMAP(1:NNBSTR(ISYM,IRS2)) = 0
do ISHLAB=1,NNSHL
  N2 = NNBSTRSH(ISYM,ISHLAB,IRS2)
  N3 = NNBSTRSH(ISYM,ISHLAB,IRS3)
  if ((N2 > 0) .and. (N3 > 0)) then
    if (N2 < N3) then
      IAB1 = IIBSTRSH(ISYM,ISHLAB,IRS2)+1
      IAB2 = IAB1+N2-1
      LAST = 0
      do IAB=IAB1,IAB2
        JAB = INDRED(IIBSTR(ISYM,IRS2)+IAB,IRS2)
        KAB = LAST
        do while (KAB < NNBSTRSH(ISYM,ISHLAB,IRS3))
          KAB = KAB+1
          LAB = IIBSTRSH(ISYM,ISHLAB,IRS3)+KAB
          MAB = INDRED(IIBSTR(ISYM,IRS3)+LAB,IRS3)
          if (MAB == JAB) then
            IMAP(IAB) = LAB
            LAST = KAB
            KAB = NNBSTRSH(ISYM,ISHLAB,IRS3)
          end if
        end do
      end do
    else
      LAB1 = IIBSTRSH(ISYM,ISHLAB,IRS3)+1
      LAB2 = LAB1+N3-1
      LAST = 0
      do LAB=LAB1,LAB2
        MAB = INDRED(IIBSTR(ISYM,IRS3)+LAB,IRS3)
        KAB = LAST
        do while (KAB < NNBSTRSH(ISYM,ISHLAB,IRS2))
          KAB = KAB+1
          IAB = IIBSTRSH(ISYM,ISHLAB,IRS2)+KAB
          JAB = INDRED(IIBSTR(ISYM,IRS2)+IAB,IRS2)
          if (JAB == MAB) then
            IMAP(IAB) = LAB
            LAST = KAB
            KAB = NNBSTRSH(ISYM,ISHLAB,IRS2)
          end if
        end do
      end do
    end if
  end if
end do

end subroutine CHO_RS2RS
