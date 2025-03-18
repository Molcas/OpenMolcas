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

subroutine NRASDT(MNRS1,MXRS1,MNRS3,MXRS3,ITOTSM,NSMST,NOCTPA,NOCTPB,IEL1A,IEL1B,NSSOA,NSSOB,IEL3A,IEL3B,NCOMB,XNCOMB,MXSB,MXSOOB, &
                  IBLTP)
! Number of combimations with symmetry ITOTSM and
!       MNRS1 - MXRS1 elecs in RAS1
!       MNRS3 - MXRS3 elecs in RAS3
!
! In view of the limited range of I*4, the number of dets
! is returned as integer and  real*8
!
! MXSB is largest UNPACKED symmetry block
! MXSOOB is largest UNPACKED symmetry-type-type block
!
! Updated with IBLTP, Summer of 93

use Symmetry_Info, only: Mul

implicit none
integer MNRS1, MXRS1, MNRS3, MXRS3, ITOTSM, NSMST, NOCTPA, NOCTPB
integer IEL1A(*), IEL1B(*)
integer NSSOA(NOCTPA,*), NSSOB(NOCTPB,*)
integer IEL3A(*), IEL3B(*)
integer NCOMB
real*8 XNCOMB
integer MXSB, MXSOOB
integer IBLTP(*)
! local variables
integer IASM, LSB, IBSM, ISYM, IATP, MXBTP, IBTP, IEL1, IEL3, LTTSBL, LTTSUP, NTEST

MXSB = 0
MXSOOB = 0
NCOMB = 0
XNCOMB = 0.0d0
do IASM=1,NSMST
  if (IBLTP(IASM) == 0) goto 300
  IBSM = Mul(IASM,ITOTSM)
  LSB = 0
  if (IBSM /= 0) then
    if (IBLTP(IASM) == 2) then
      ISYM = 1
    else
      ISYM = 0
    end if
    do IATP=1,NOCTPA
      if (ISYM == 1) then
        MXBTP = IATP
      else
        MXBTP = NOCTPB
      end if
      do IBTP=1,MXBTP
        IEL1 = IEL1A(IATP)+IEL1B(IBTP)
        IEL3 = IEL3A(IATP)+IEL3B(IBTP)
        if ((IEL1 >= MNRS1) .and. (IEL1 <= MXRS1) .and. (IEL3 >= MNRS3) .and. (IEL3 <= MXRS3)) then
          ! Size of unpacked block
          LTTSUP = NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
          ! Size of packed block
          if ((ISYM == 0) .or. (IATP /= IBTP)) then
            LTTSBL = NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
          else
            LTTSBL = NSSOA(IATP,IASM)*(NSSOA(IATP,IASM)+1)/2
          end if
          NCOMB = NCOMB+LTTSBL
          LSB = LSB+LTTSUP
          MXSOOB = max(MXSOOB,LTTSUP)
          if ((ISYM == 0) .or. (IATP /= IBTP)) then
            XNCOMB = XNCOMB+dble(NSSOA(IATP,IASM))*dble(NSSOB(IBTP,IBSM))
          else
            XNCOMB = XNCOMB+dble(NSSOA(IATP,IASM))*(dble(NSSOB(IBTP,IBSM))+1.0d0)/2.0d0
          end if
        end if
      end do
    end do
    MXSB = max(MXSB,LSB)
  end if
300 continue
end do

NTEST = 0
if (NTEST /= 0) write(6,*) ' NCOMB and XNCOMB ',NCOMB,XNCOMB

end subroutine NRASDT
