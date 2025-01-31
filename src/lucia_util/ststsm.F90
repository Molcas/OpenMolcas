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

subroutine STSTSM(STSTSX,STSTDX,NSMST)
! construct  STSTSX and STSTDX giving
! symmetry of sx (dx) connecting two given string symmetries

implicit real*8(A-H,O-Z)
integer STSTSX(NSMST,NSMST), STSTDX(NSMST,NSMST)

do ILSTSM=1,NSMST
  do IRSTSM=1,NSMST
    call SYMCOM(1,5,ISXSM,IRSTSM,ILSTSM)
    call SYMCOM(1,6,IDXSM,IRSTSM,ILSTSM)
    STSTSX(ILSTSM,IRSTSM) = ISXSM
    STSTDX(ILSTSM,IRSTSM) = IDXSM
  end do
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' STSTSM : STSTSX, STSTDX'
  call IWRTMA(STSTSX,NSMST,NSMST,NSMST,NSMST)
  call IWRTMA(STSTDX,NSMST,NSMST,NSMST,NSMST)
end if

end subroutine STSTSM
