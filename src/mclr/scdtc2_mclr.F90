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
subroutine SCDTC2_MCLR(RASVEC,ISM,ICBLTP,NSM,NOCTPA,NOCTPB,NSASO,NSBSO,IOCOC,IDC,IWAY,IMMLST)
! Scale elements of a RAS vector to transfer between
! combinations and packed determinants
! IWAY = 1 : dets to combs
! IWAY = 2 : combs to dets
! Combination storage mode is defined BY IDC
!
! General symmetry version, Feb 1991

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(inout) :: RASVEC(*)
integer(kind=iwp), intent(in) :: ISM, ICBLTP(*), NSM, NOCTPA, NOCTPB, NSASO(NOCTPA,*), NSBSO(NOCTPB,*), IOCOC(NOCTPA,NOCTPB), IDC, &
                                 IWAY, IMMLST(*)
integer(kind=iwp) :: IASM, IATP, IBASE, IBSM, IBTP, IBTPMX, NELMNT, NIA, NIB
real(kind=wp) :: FACTOR
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQ2I = sqrt(Half)

#ifdef _DEBUGPRINT_
write(u6,*) ' Information from SCDTC2'
write(u6,*) ' ======================='
write(u6,*) ' Input vector'
call WRTRS2_MCLR(RASVEC,ISM,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSM)
#endif

IBASE = 1
do IASM=1,NSM

  IBSM = Mul(IASM,ISM)
  if ((IBSM == 0) .or. (ICBLTP(IASM) == 0)) cycle
  do IATP=1,NOCTPA
    if (ICBLTP(IASM) == 2) then
      IBTPMX = IATP
    else
      IBTPMX = NOCTPB
    end if
    NIA = NSASO(IATP,IASM)
    do IBTP=1,IBTPMX
      if (IOCOC(IATP,IBTP) == 0) cycle
      ! Number of elements in this block
      NIB = NSBSO(IBTP,IBSM)
      if ((ICBLTP(IASM) == 2) .and. (IATP == IBTP)) then
        NELMNT = nTri_Elem(NIA)
      else
        NELMNT = NIA*NIB
      end if

      if (IDC == 2) then
        if (IWAY == 1) then
          FACTOR = SQ2
        else
          FACTOR = SQ2I
        end if
        RASVEC(IBASE:IBASE+NELMNT-1) = FACTOR*RASVEC(IBASE:IBASE+NELMNT-1)
        if ((IASM == IBSM) .and. (IATP == IBTP)) then
          FACTOR = One/FACTOR
          call SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
        end if
      else if ((IDC == 3) .and. (IMMLST(IASM) /= IASM)) then
        if (IWAY == 1) then
          FACTOR = SQ2
        else
          FACTOR = SQ2I
        end if
        RASVEC(IBASE:IBASE+NELMNT-1) = FACTOR*RASVEC(IBASE:IBASE+NELMNT-1)
      else if (IDC == 4) then
        !Ml Ms combinations
        if (IWAY == 1) then
          if (IASM == IBSM) then
            FACTOR = SQ2
          else
            FACTOR = Two
          end if
        else !if (IWAY == 2) then
          if (IASM == IBSM) then
            FACTOR = SQ2I
          else
            FACTOR = Half
          end if
        end if
        RASVEC(IBASE:IBASE+NELMNT-1) = FACTOR*RASVEC(IBASE:IBASE+NELMNT-1)
        if (IATP == IBTP) then
          if (IWAY == 1) then
            FACTOR = SQ2I
          else !if (IWAY == 2) then
            FACTOR = SQ2
          end if
          call SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
        end if
      end if

      IBASE = IBASE+NELMNT
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Scaled vector'
call WRTRS2_MCLR(RASVEC,ISM,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSM)
#endif

return

end subroutine SCDTC2_MCLR
