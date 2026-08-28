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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

subroutine SPNCOM_MCLR(iwork,NOPEN,MS2,NDET,IABDET,IABUPP,IFLAG,PSSIGN)
! Combinations of nopen unpaired electrons.Required
! spin projection MS2/2.
! JO 21-7-84
!    IFLAG = 1 : Only combinations (in IABDET)
!    IFLAG = 2 : combinations and upper dets
!    IFLAG = 3 : Only upper dets
! A few revisions october 1988
! Upper dets added feb 1989
! Changed to combinations June 1992
!
! If PSSIGN differs from 0, spin combinations are assumed.
! we select as the unique determinants those with first electron
! having alpha spin

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NOPEN, MS2, IFLAG
integer(kind=iwp), intent(out) :: IWORK(NOPEN+1), NDET
integer(kind=iwp), intent(inout) :: IABDET(NOPEN,*), IABUPP(NOPEN,*)
real(kind=wp) :: PSSIGN
integer(kind=iwp) :: ADD, I, IEL, J, lUPPER, MS2L, NALPHA, NUPPER

! LENGTH OF IWORK MUST BE AT LEAST NOPEN+1

NDET = 0
NUPPER = 0

! Determinants are considered as binary numbers,1=alpha,0=beta

IWORK(:) = 0
! Loop over all possible binary numbers
do I=1,2**NOPEN
  ! 1 : NEXT BINARY NUMBER
  ADD = 1
  J = 0
  do while (ADD == 1)
    J = J+1
    if (IWORK(J) == 1) then
      IWORK(J) = 0
    else
      IWORK(J) = 1
      ADD = 0
    end if
  end do

  ! 2 :  CORRECT SPIN PROJECTION ?
  NALPHA = sum(IWORK(1:NOPEN))

  if ((2*NALPHA-NOPEN == MS2) .and. ((PSSIGN == Zero) .or. (IWORK(1) /= 0))) then
    if (IFLAG < 3) then
      NDET = NDET+1
      IABDET(:,NDET) = IWORK(1:NOPEN)
    end if
    if (IFLAG > 1) then
      ! UPPER DET ?
      MS2L = 0
      LUPPER = 1

      do IEL=1,NOPEN
        if (IWORK(IEL) == 1) then
          MS2L = MS2L+1
        else
          MS2L = MS2L-1
        end if
        if (MS2L < 0) LUPPER = 0
      end do
      if (LUPPER == 1) then
        NUPPER = NUPPER+1
        IABUPP(:,NUPPER) = IWORK(1:NOPEN)
      end if
    end if
  end if

end do

!XMSD2 = real(MS2,kind=wp)*Half

end subroutine SPNCOM_MCLR
