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

subroutine WRTRS2_MCLR(VECTOR,ISM,ICBLTP,IOCOC,NOCTPA,NOCTPB,NSASO,NSBSO,NSM)
! Write RAS vector. Storage form is defined by ICBLTP

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: VECTOR(*)
integer(kind=iwp), intent(in) :: ISM, ICBLTP(*), NOCTPA, NOCTPB, IOCOC(NOCTPA,NOCTPB), NSASO(NOCTPA,*), NSBSO(NOCTPB,*), NSM
integer(kind=iwp) :: IASM, IATP, IBASE, IBSM, IBTP, IBTPMX, NAST, NBST, NELMNT

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
    NAST = NSASO(IATP,IASM)

    do IBTP=1,IBTPMX
      if (IOCOC(IATP,IBTP) == 0) cycle
      NBST = NSBSO(IBTP,IBSM)
      if ((ICBLTP(IASM) == 2) .and. (IATP == IBTP)) then
        ! Diagonal block
        NELMNT = nTri_Elem(NAST)
        if (NELMNT /= 0) then
          write(u6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
          write(u6,'(A)') '  ============================'
          !do i=1,nTri_Elem(NAST)  !yma
          !  if (abs(VECTOR(IBASE+i-1)) < 1.0e-16_wp) VECTOR(IBASE+i-1) = Zero
          !  write(u6,*) 'vector',i,VECTOR(IBASE+i-1)
          !end do
          call PRSM2(VECTOR(IBASE),NAST)
          IBASE = IBASE+NELMNT
        end if
      else
        NELMNT = NAST*NBST
        if (NELMNT /= 0) then
          write(u6,'(A,3I3)') '  Iasm iatp ibtp : ',IASM,IATP,IBTP
          write(u6,'(A)') '  ============================'
          call WRTMAT(VECTOR(IBASE),NAST,NBST,NAST,NBST)
          IBASE = IBASE+NELMNT
        end if
      end if
    end do
  end do
end do

return

end subroutine WRTRS2_MCLR
