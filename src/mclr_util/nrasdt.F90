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

!#define _DEBUGPRINT_
subroutine NRASDT(MNRS1,MXRS1,MNRS3,MXRS3,ITOTSM,NSM,NOCTPA,NOCTPB,IEL1A,IEL1B,NSSOA,NSSOB,IEL3A,IEL3B,NCOMB,XNCOMB,MXSB,MXSOOB, &
                  IBLTP)
! Number of combimations with symmetry ITOTSM and
!       MNRS1 - MXRS1 elecs in RAS1
!       MNRS3 - MXRS3 elecs in RAS3
!
! In view of the limited range of I*4, the number of dets
! is returned as integer and real(kind=wp)
!
! MXSB is largest UNPACKED symmetry block
! MXSOOB is largest UNPACKED symmetry-type-type block
!
! Updated with IBLTP, Summer of 93

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MNRS1, MXRS1, MNRS3, MXRS3, ITOTSM, NSM, NOCTPA, NOCTPB, IEL1A(*), IEL1B(*), NSSOA(NOCTPA,*), &
                                 NSSOB(NOCTPB,*), IEL3A(*), IEL3B(*), IBLTP(*)
integer(kind=iwp), intent(out) :: NCOMB, MXSB, MXSOOB
real(kind=wp), intent(out) :: XNCOMB
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IEL1, IEL3, ISYM, LSB, LTTSBL, LTTSUP, MXBTP

MXSB = 0
MXSOOB = 0
NCOMB = 0
XNCOMB = Zero
do IASM=1,NSM
  if (IBLTP(IASM) == 0) cycle
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
            LTTSBL = nTri_Elem(NSSOA(IATP,IASM))
          end if
          NCOMB = NCOMB+LTTSBL
          LSB = LSB+LTTSUP
          MXSOOB = max(MXSOOB,LTTSUP)
          if ((ISYM == 0) .or. (IATP /= IBTP)) then
            XNCOMB = XNCOMB+real(NSSOA(IATP,IASM),kind=wp)*real(NSSOB(IBTP,IBSM),kind=wp)
          else
            XNCOMB = XNCOMB+real(NSSOA(IATP,IASM),kind=wp)*real(NSSOB(IBTP,IBSM)+1,kind=wp)*Half
          end if
        end if
      end do
    end do
    MXSB = max(MXSB,LSB)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' NCOMB and XNCOMB ',NCOMB,XNCOMB
#endif

end subroutine NRASDT
