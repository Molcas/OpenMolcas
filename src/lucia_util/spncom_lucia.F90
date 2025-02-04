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

subroutine SPNCOM_LUCIA(NOPEN,MS2,NDET,IABDET,IABUPP,IFLAG,PSSIGN,IPRCSF)
! Combinations of nopen unpaired electrons.Required
! spin projection MS2/2.

use lucia_data, only: MXPORB
use Constants, only: Zero, Half
use Definitions, only: wp, u6

implicit none
integer NOPEN, MS2, NDET, IFLAG, IPRCSF
real*8 PSSIGN
integer IABDET(NOPEN,*), IABUPP(NOPEN,*)
integer ADD
! Should have length of max number of open orbitals
integer IWORK(MXPORB)
! LENGTH OF IWORK MUST BE AT LEAST NOPEN
integer NTEST, NUPPER, MX, IFIRST, I, IZERO, J, NALPHA, MS2L, LUPPER, IEL, K
real*8 XMSD2

NTEST = 0
NTEST = max(NTEST,IPRCSF)
NDET = 0
NUPPER = 0

! Determinants are considered as binary numbers,1=alpha,0=beta

MX = 2**NOPEN
call ISETVC(IWORK,0,NOPEN)
IFIRST = 1
! Loop over all possible binary numbers
do I=1,MX
  if (IFIRST == 1) then
    ! Initial number
    IZERO = 0
    call ISETVC(IWORK,IZERO,NOPEN)
    IFIRST = 0
  else
    ! Next number
    ADD = 1
    J = 0
190 continue
    J = J+1
    if (IWORK(J) == 1) then
      IWORK(J) = 0
    else
      IWORK(J) = 1
      ADD = 0
    end if
    if (ADD == 1) goto 190
  end if

  ! 2 :  CORRECT SPIN PROJECTION ?
  NALPHA = 0
  do J=1,NOPEN
    NALPHA = NALPHA+IWORK(J)
  end do

  if ((2*NALPHA-NOPEN == MS2) .and. ((PSSIGN == Zero) .or. (IWORK(1) /= 0))) then
    if (IFLAG < 3) then
      NDET = NDET+1
      call ICOPVE(IWORK,IABDET(1,NDET),NOPEN)
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
        call ICOPVE(IWORK,IABUPP(1,NUPPER),NOPEN)
      end if
    end if
  end if

end do

XMSD2 = Half*real(MS2,kind=wp)

if ((NTEST >= 5) .and. (IFLAG /= 3)) then
  write(u6,1010) NOPEN,NDET,XMSD2
  write(u6,*)
  write(u6,'(A)') '  Combinations :'
  write(u6,'(A)') '  =============='
  write(u6,*)
  do J=1,NDET
    write(u6,1020) J,(IABDET(K,J),K=1,NOPEN)
  end do
end if

if ((IFLAG > 1) .and. (NTEST >= 5)) then
  write(u6,*)
  write(u6,'(A)') ' Upper determinants'
  write(u6,'(A)') ' =================='
  write(u6,*)
  do J=1,NUPPER
    write(u6,1020) J,(IABUPP(K,J),K=1,NOPEN)
  end do
end if

return
1010 format('0',2X,I3,' Unpaired electrons give ',I5,/,'           combinations with spin projection ',F12.7)
1020 format('0',I5,2X,30I2,/,(1X,7X,30I2))

end subroutine SPNCOM_LUCIA
