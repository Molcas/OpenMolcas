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

subroutine NGASDT(IOCCMN,IOCCMX,NGAS,ITOTSM,NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,IAOCC,IBOCC,MXPNGAS,NCOMB,XNCOMB,MXSB,MXSOOB,IBLTP, &
                  NTTSBL,LCOL,IOCOC,MXSOOB_AS)
! Number of combinations with symmetry ITOTSM and
! occupation between IOCCMN and IOCCMX
!
! In view of the limited range of I*4, the number of dets
! is returned as integer and real*8
!
! MXSB is largest UNPACKED symmetry block
! MXSOOB is largest UNPACKED symmetry-type-type block
! NTTSBL is number of TTS blocks in vector
! LCOL is the sum of the number of columns in each block
!
! Winter 94/95
! May 1999 : Loops restructrured to sym,type,type (leftmost innerst)
!            MXSB not calculated

use Constants, only: Zero, Half
use Definitions, only: wp, u6

implicit real*8(A-H,O-Z)
! Allowed combinations of alpha and beta types
integer IOCOC(NOCTPA,NOCTPB)
! Occupation constraints
dimension IOCCMN(NGAS), IOCCMX(NGAS)
! Occupation of alpha and beta strings
dimension IAOCC(MXPNGAS,*), IBOCC(MXPNGAS,*)
! Number of strings per supergroup and symmetry
dimension NSSOA(NSMST,*), NSSOB(NSMST,*)
! block types
dimension IBLTP(*)

NTEST = 0
if (NTEST >= 5) then
  write(u6,*) ' NGASDT speaking'
  write(u6,*) ' ==============='
  write(u6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
  write(u6,*) ' ITOTSM ',ITOTSM
  write(u6,*) ' Upper and lower occupation constraints'
  call IWRTMA(IOCCMN,1,NGAS,1,NGAS)
  call IWRTMA(IOCCMX,1,NGAS,1,NGAS)
  write(u6,*) ' IOCOC matrix'
  call IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  write(u6,*) ' Number of alpha and beta strings'
  call IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
  call IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
end if

MXSB = 0
MXSOOB = 0
MXSOOB_AS = 0
NCOMB = 0
XNCOMB = Zero
NTTSBL = 0
LCOL = 0

do IATP=1,NOCTPA
  do IBTP=1,NOCTPB

    if (NTEST >= 10) then
      write(u6,*) ' Alpha super group and beta super group'
      call IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
      call IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
    end if

    if (IOCOC(IATP,IBTP) == 1) then

      LTTS_AS = 0
      do IASM=1,NSMST
        if (IBLTP(IASM) == 0) goto 300
        call SYMCOM(2,IASM,IBSM,ITOTSM)
        if (IBSM /= 0) then
          if (IBLTP(IASM) == 2) then
            ISYM = 1
          else
            ISYM = 0
          end if
          if ((ISYM == 1) .and. (IBTP > IATP)) goto 300
          LASTR = NSSOA(IASM,IATP)
          LBSTR = NSSOB(IBSM,IBTP)
          ! Size of unpacked block
          LTTSUP = LASTR*LBSTR
          ! Size of packed block
          if ((ISYM == 0) .or. (IATP /= IBTP)) then
            LTTSBL = LASTR*LBSTR
            XNCOMB = XNCOMB+real(LASTR,kind=wp)*real(LBSTR,kind=wp)
          else
            LTTSBL = LASTR*(LASTR+1)/2
            XNCOMB = XNCOMB+Half*real(LASTR+1,kind=wp)*real(LASTR,kind=wp)
          end if
          LTTS_AS = LTTS_AS+LTTSUP
          NCOMB = NCOMB+LTTSBL
          MXSOOB = max(MXSOOB,LTTSUP)
          NTTSBL = NTTSBL+1
          LCOL = LCOL+NSSOB(IBSM,IBTP)
        end if
300     continue
      end do
      MXSOOB_AS = max(MXSOOB_AS,LTTS_AS)
    end if
  end do
end do

if (NTEST >= 1) write(u6,*) ' NGASDT : NCOMB XNCOMB, NTTSBL',NCOMB,XNCOMB,NTTSBL

end subroutine NGASDT
