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
subroutine SCDTC2_MCLR(RASVEC,ISMOST,ICBLTP,NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IOCOC,IDC,IWAY,IMMLST)
! Scale elements of a RAS vector to transfer between
! combinations and packed determinants
! IWAY = 1 : dets to combs
! IWAY = 2 : combs to dets
! Combination storage mode is defined BY IDC
!
! General symmetry version, Feb 1991

implicit real*8(A-H,O-Z)
dimension RASVEC(*), NSASO(NOCTPA,*), NSBSO(NOCTPB,*)
dimension IOCOC(NOCTPA,NOCTPB)
dimension ISMOST(*), ICBLTP(*), IMMLST(*)

#ifdef _DEBUGPRINT_
write(6,*) ' Information from SCDTC2'
write(6,*) ' ======================='
write(6,*) ' Input vector'
call WRTRS2_MCLR(RASVEC,ISMOST,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
#endif

SQ2 = sqrt(2.0d0)
SQ2I = 1.0d0/SQ2

IBASE = 1
do IASM=1,NSMST

  IBSM = ISMOST(IASM)
  if ((IBSM == 0) .or. (ICBLTP(IASM) == 0)) goto 200
  do IATP=1,NOCTPA
    if (ICBLTP(IASM) == 2) then
      IBTPMX = IATP
    else
      IBTPMX = NOCTPB
    end if
    NIA = NSASO(IATP,IASM)
    do IBTP=1,IBTPMX
      if (IOCOC(IATP,IBTP) == 0) goto 50
      ! Number of elements in this block
      call xflush(6)
      NIB = NSBSO(IBTP,IBSM)
      if ((ICBLTP(IASM) == 2) .and. (IATP == IBTP)) then
        NELMNT = NIA*(NIA+1)/2
      else
        NELMNT = NIA*NIB
      end if

      if (IDC == 2) then
        if (IWAY == 1) then
          FACTOR = SQ2
        else
          FACTOR = SQ2I
        end if
        call DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
        if ((IASM == IBSM) .and. (IATP == IBTP)) then
          FACTOR = 1.0d0/FACTOR
          call SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
        end if
      else if ((IDC == 3) .and. (IMMLST(IASM) /= IASM)) then
        if (IWAY == 1) then
          FACTOR = SQ2
        else
          FACTOR = SQ2I
        end if
        call DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
        !Ml Ms combinations
        call xflush(6)
      else if (IDC == 4) then
        if (IWAY == 1) then
          if (IASM == IBSM) then
            FACTOR = SQ2
          else
            FACTOR = 2.0d0
          end if
        else if (IWAY == 2) then
          if (IASM == IBSM) then
            FACTOR = SQ2I
          else
            FACTOR = 0.5d0
          end if
        end if
        call DSCAL_(NELMNT,FACTOR,RASVEC(IBASE),1)
        if (IATP == IBTP) then
          if (IWAY == 1) then
            FACTOR = SQ2I
          else if (IWAY == 2) then
            FACTOR = SQ2
          end if
          call SCLDIA(RASVEC(IBASE),FACTOR,NIA,1)
        end if
      end if

      IBASE = IBASE+NELMNT
50    continue
    end do
  end do
200 continue
end do

#ifdef _DEBUGPRINT_
write(6,*) ' Scaled vector'
call xflush(6)
call WRTRS2_MCLR(RASVEC,ISMOST,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
#endif

return

end subroutine SCDTC2_MCLR
