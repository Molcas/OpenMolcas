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

subroutine WRTRS2(VECTOR,ISMOST,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSMST)
! Write RAS vector . Storage form is defined by ICBLTP

implicit real*8(A-H,O-Z)

dimension VECTOR(*)
dimension IOCOC(NOCTPA,NOCTPB)
dimension NSASO(NSMST,*), NSBSO(NSMST,*)
dimension ICBLTP(*), ISMOST(*)

IBASE = 1
do IASM=1,NSMST
  IBSM = ISMOST(IASM)
  if ((IBSM == 0) .or. (ICBLTP(IASM) == 0)) goto 1000

  do IATP=1,NOCTPA
    if (ICBLTP(IASM) == 2) then
      IBTPMX = IATP
    else
      IBTPMX = NOCTPB
    end if
    NAST = NSASO(IASM,IATP)

    do IBTP=1,IBTPMX
      if (IOCOC(IATP,IBTP) == 0) goto 800
      NBST = NSBSO(IBSM,IBTP)
      if ((ICBLTP(IASM) == 2) .and. (IATP == IBTP)) then
        ! Diagonal block
        NELMNT = NAST*(NAST+1)/2
        if (NELMNT /= 0) then
          write(6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
          write(6,'(A)') '  ============================'
          call PRSM2(VECTOR(IBASE),NAST)
          IBASE = IBASE+NELMNT
        end if
      else
        NELMNT = NAST*NBST
        if (NELMNT /= 0) then
          write(6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
          write(6,'(A)') '  ============================'
          call WRTMAT(VECTOR(IBASE),NAST,NBST,NAST,NBST)
          IBASE = IBASE+NELMNT
        end if
      end if
800   continue
    end do
  end do
1000 continue
end do

end subroutine WRTRS2
