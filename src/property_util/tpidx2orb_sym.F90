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

subroutine tpidx2orb_sym(TYPEINDEX,NB,NF,NI,N1,N2,N3,NS,ND)
!SVC: read orbital partition info from a typeindex array
!     corresponding to a _specific_symmetry_ (so these variables are
!     scalars!!
!     A typeindex array consists of integers with 7 possible values
!     corresponding to the types 'fi123sd' -> '1234567'.

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: typeindex(*), NB
integer(kind=iwp), intent(out) :: NF, NI, N1, N2, N3, NS, ND
integer(kind=iwp) :: IB, ITYPE

NF = 0
NI = 0
N1 = 0
N2 = 0
N3 = 0
NS = 0
ND = 0
! switch typeindex:
do IB=1,NB
  ITYPE = TYPEINDEX(IB)
  if (ITYPE == 1) then
    NF = NF+1
  else if (ITYPE == 2) then
    NI = NI+1
  else if (ITYPE == 3) then
    N1 = N1+1
  else if (ITYPE == 4) then
    N2 = N2+1
  else if (ITYPE == 5) then
    N3 = N3+1
  else if (ITYPE == 6) then
    NS = NS+1
  else if (ITYPE == 7) then
    ND = ND+1
  else
    write(u6,*) 'TPIDX2ORB_SYM: unknown type index number'
    call AbEnd()
  end if
end do

end subroutine tpidx2orb_sym
