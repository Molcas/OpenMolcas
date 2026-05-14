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
subroutine NGASDT(ITOTSM,NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,NCOMB,XNCOMB,MXSOOB,IBLTP,NTTSBL,LCOL,IOCOC,MXSOOB_AS)
! Number of combinations with symmetry ITOTSM and
! occupation between IOCCMN and IOCCMX
!
! In view of the limited range of I*4, the number of dets
! is returned as integer and real*8
!
! IOCOC  : Allowed combinations of alpha and beta types
! NSSOA  : Number of strings per supergroup and symmetry
! IBLTP  : Block types
!
! MXSB is largest UNPACKED symmetry block
! MXSOOB is largest UNPACKED symmetry-type-type block
! NTTSBL is number of TTS blocks in vector
! LCOL is the sum of the number of columns in each block
!
! Winter 94/95
! May 1999 : Loops restructrured to sym,type,type (leftmost innerst)
!            MXSB not calculated

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ITOTSM, NSMST, NOCTPA, NOCTPB, NSSOA(NSMST,NOCTPA), NSSOB(NSMST,NOCTPB), IBLTP(NSMST), &
                                 IOCOC(NOCTPA,NOCTPB)
integer(kind=iwp), intent(out) :: NCOMB, MXSOOB, NTTSBL, LCOL, MXSOOB_AS
real(kind=wp) :: XNCOMB
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, ISYM, LASTR, LBSTR, LTTS_AS, LTTSBL, LTTSUP

#ifdef _DEBUGPRINT_
write(u6,*) ' NGASDT speaking'
write(u6,*) ' ==============='
write(u6,*) ' NOCTPA,NOCTPB ',NOCTPA,NOCTPB
write(u6,*) ' ITOTSM ',ITOTSM
write(u6,*) ' IOCOC matrix'
call IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
write(u6,*) ' Number of alpha and beta strings'
call IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
call IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
#endif

MXSOOB = 0
MXSOOB_AS = 0
NCOMB = 0
XNCOMB = Zero
NTTSBL = 0
LCOL = 0

do IATP=1,NOCTPA
  do IBTP=1,NOCTPB

    if (IOCOC(IATP,IBTP) == 1) then

      LTTS_AS = 0
      do IASM=1,NSMST
        if (IBLTP(IASM) == 0) cycle
        IBSM = Mul(IASM,ITOTSM)
        if (IBSM /= 0) then
          if (IBLTP(IASM) == 2) then
            ISYM = 1
          else
            ISYM = 0
          end if
          if ((ISYM == 1) .and. (IBTP > IATP)) cycle
          LASTR = NSSOA(IASM,IATP)
          LBSTR = NSSOB(IBSM,IBTP)
          ! Size of unpacked block
          LTTSUP = LASTR*LBSTR
          ! Size of packed block
          if ((ISYM == 0) .or. (IATP /= IBTP)) then
            LTTSBL = LASTR*LBSTR
            XNCOMB = XNCOMB+real(LASTR,kind=wp)*real(LBSTR,kind=wp)
          else
            LTTSBL = nTri_Elem(LASTR)
            XNCOMB = XNCOMB+Half*real(LASTR+1,kind=wp)*real(LASTR,kind=wp)
          end if
          LTTS_AS = LTTS_AS+LTTSUP
          NCOMB = NCOMB+LTTSBL
          MXSOOB = max(MXSOOB,LTTSUP)
          NTTSBL = NTTSBL+1
          LCOL = LCOL+NSSOB(IBSM,IBTP)
        end if
      end do
      MXSOOB_AS = max(MXSOOB_AS,LTTS_AS)
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' NGASDT : NCOMB XNCOMB, NTTSBL',NCOMB,XNCOMB,NTTSBL
#endif

end subroutine NGASDT
