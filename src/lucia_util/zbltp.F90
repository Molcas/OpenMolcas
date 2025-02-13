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

subroutine ZBLTP(ISMOST,MAXSYM,IDC,ICBLTP,IMMLST)
! Generate vector ICBLTP giving type of each block
!
! ICBLTP gives type of symmetry block :
! = 0 : symmetry block is not included
! = 1 : symmetry block is included, all OO types
! = 2 : symmetry block is included, lower OO types

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MAXSYM, ISMOST(MAXSYM), IDC, IMMLST(MAXSYM)
integer(kind=iwp), intent(out) :: ICBLTP(MAXSYM)
integer(kind=iwp) :: IASYM, IBSYM, NTEST

! Changed to simplify structure for IDC <= 2
if (IDC <= 2) then
  ! No spatial degeneracy
  do IASYM=1,MAXSYM
    IBSYM = ISMOST(IASYM)
    if ((IDC == 2) .and. (IBSYM > IASYM)) then
      ! Symmetry block excluded
      ICBLTP(IASYM) = 0
    else if (((IDC == 2) .and. (IASYM > IBSYM)) .or. (IDC == 1)) then
      ! Complete symmetry block included
      ICBLTP(IASYM) = 1
    else
      ! Lower half  symmetry block included
      ICBLTP(IASYM) = 2
    end if
  end do
else
  ! Also spatial degeneracy
  do IASYM=1,MAXSYM

    IBSYM = ISMOST(IASYM)
    if (IBSYM == 0) cycle
    if ((((IDC == 2) .or. (IDC == 4)) .and. (IBSYM > IASYM)) .or. ((IDC == 3) .and. (IMMLST(IASYM) > IASYM))) then
      ! Symmetry block excluded
      ICBLTP(IASYM) = 0
    else if (((IDC == 2) .and. (IASYM > IBSYM)) .or. (IDC == 1) .or. ((IDC == 3) .and. (IASYM >= IMMLST(IASYM)))) then
      ! Complete symmetry block included
      ICBLTP(IASYM) = 1
    else
      ! Lower half  symmetry block included
      ICBLTP(IASYM) = 2
    end if
  end do
end if
! End of IDC switch

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' Block type of symmetry blocks'
  call IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
end if

end subroutine ZBLTP
