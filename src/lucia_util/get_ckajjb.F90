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

subroutine GET_CKAJJB(CB,NJ,NJA,CKAJJB,NKA,NJB,J,ISCA,SSCA)
! Obtain for given orbital index j the gathered matrix
!
! C(Ka,j,Jb) = SSCA(Ka)C(Jb,ISCA(Ka))
!
! For efficient processing of alpha-beta loop

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NJ, NJA, NKA, NJB, J, ISCA(*)
real(kind=wp), intent(in) :: CB(NJB,NJA), SSCA(*)
real(kind=wp), intent(_OUT_) :: CKAJJB(*)
integer(kind=iwp) :: IADR, IADR0, ICBL, ICEND, ICOFF, ICONST, IROW, JB, KA, LBLK, NBLK
real(kind=wp) :: S

!write(u6,*) ' From GET_CKAJJB'
!LBLK = 100
LBLK = 40
NBLK = NJB/LBLK
if (LBLK*NBLK < NJB) NBLK = NBLK+1
ICOFF = 1
do ICBL=1,NBLK
  if (ICBL > 1) ICOFF = ICOFF+LBLK
  ICEND = min(ICOFF+LBLK-1,NJB)
  ICONST = NKA*NJ
  IADR0 = (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
  if (ICEND > ICOFF) then
    ! Inner loop over JB
    do KA=1,NKA
      if (ISCA(KA) /= 0) then
        S = SSCA(KA)
        IROW = ISCA(KA)
        !IADR = KA+(J-1)*NKA+(ICOFF-1-1)*NKA*NJ
        IADR = IADR0+KA
        do JB=ICOFF,ICEND
          ! Address of C(Ka,j,Jb)
          IADR = IADR+ICONST
          CKAJJB(IADR) = S*CB(JB,IROW)
        end do
      else
        IADR = IADR0+KA
        do JB=ICOFF,ICEND
          !IADR = KA+(J-1)*NKA+(JB-1)*NKA*NJ
          IADR = IADR+ICONST
          CKAJJB(IADR) = Zero
        end do
      end if
    end do
  else
    ! No inner loop over JB
    do KA=1,NKA
      if (ISCA(KA) /= 0) then
        S = SSCA(KA)
        IROW = ISCA(KA)
        !IADR = KA+(J-1)*NKA+(ICOFF-1-1)*NKA*NJ
        IADR = IADR0+KA
        !do JB=ICOFF,ICEND
        ! Address of C(Ka,j,Jb)
        IADR = IADR+ICONST
        CKAJJB(IADR) = S*CB(ICOFF,IROW)
        !end do
      else
        IADR = IADR0+KA
        !do JB=ICOFF,ICEND
        !  IADR = KA + (J-1)*NKA+(JB-1)*NKA*NJ
        IADR = IADR+ICONST
        CKAJJB(IADR) = Zero
        !end do
      end if
    end do
  end if
  ! End of test ICEND,ICOFF
end do

end subroutine GET_CKAJJB
