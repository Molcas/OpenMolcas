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

subroutine ADD_SKAIIB(SB,NI,NIA,SKAIIB,NKA,NIB,I,ISCA,SSCA)
! Update Transposed sigma block with contributions for given orbital index j
! from the matrix S(Ka,i,Ib)
!
! S(Ib,Isca(Ka)) =  S(Ib,Isca(Ka)) + Ssca(Ka)*S(Ka,I,Ib)
!
! For efficient processing of alpha-beta loop

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NI, NIA, NKA, NIB, I, ISCA(*)
real(kind=wp), intent(inout) :: SB(NIB,NIA)
real(kind=wp), intent(in) :: SKAIIB(*), SSCA(*)
integer(kind=iwp) :: IADR, IADR0, IB, ICBL, ICEND, ICOFF, ICONST, IROW, KA, LBLK, NBLK
real(kind=wp) :: S

! To get rid of annoying and incorrect compiler warnings
ICOFF = 0

!LBLK = 100
LBLK = 40
NBLK = NIB/LBLK
if (LBLK*NBLK < NIB) NBLK = NBLK+1
do ICBL=1,NBLK
  if (ICBL == 1) then
    ICOFF = 1
  else
    ICOFF = ICOFF+LBLK
  end if
  ICEND = min(ICOFF+LBLK-1,NIB)
  ICONST = NKA*NI
  IADR0 = (I-1)*NKA+(ICOFF-1-1)*NKA*NI
  if (ICEND > ICOFF) then
    ! Use form with Inner loop over IB
    do KA=1,NKA
      if (ISCA(KA) /= 0) then
        S = SSCA(KA)
        IROW = ISCA(KA)
        !IADR = KA+(I-1)*NKA+(ICOFF-1-1)*NKA*NI
        IADR = IADR0+KA
        do IB=ICOFF,ICEND
          ! Address of S(Ka,i,Ib)
          IADR = IADR+ICONST
          SB(IB,IROW) = SB(IB,IROW)+S*SKAIIB(IADR)
        end do
      end if
    end do
  else
    ! Form with no loop over IB
    do KA=1,NKA
      if (ISCA(KA) /= 0) then
        S = SSCA(KA)
        IROW = ISCA(KA)
        IADR = IADR0+KA+ICONST
        !do IB=ICOFF,ICEND
        ! Address of S(Ka,i,Ib)
        !IADR = IADR+ICONST
        SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
        !end do
      end if
    end do
  end if
  ! End of test of ICOFF=ICEND
end do

end subroutine ADD_SKAIIB
