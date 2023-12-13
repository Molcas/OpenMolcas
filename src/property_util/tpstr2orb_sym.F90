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

subroutine tpstr2orb_sym(TYPESTRING,NB,NF,NI,N1,N2,N3,NS,ND)
!SVC: read orbital partition info from a typestring array
!     corresponding to a _specific_symmetry_ (so these variables are
!     scalars!!
!     A typestring array consists of characters of 'fi123sd'

use Definitions, only: iwp, u6

implicit none
character, intent(in) :: typestring(*)
integer(kind=iwp), intent(in) :: NB
integer(kind=iwp), intent(out) :: NF, NI, N1, N2, N3, NS, ND
integer(kind=iwp) :: IB
character :: CTYPE

NF = 0
NI = 0
N1 = 0
N2 = 0
N3 = 0
NS = 0
ND = 0
! switch typeindex:
do IB=1,NB
  CTYPE = TYPESTRING(IB)
  call UPCASE(CTYPE)
  if (CTYPE == 'F') then
    NF = NF+1
  else if (CTYPE == 'I') then
    NI = NI+1
  else if (CTYPE == '1') then
    N1 = N1+1
  else if (CTYPE == '2') then
    N2 = N2+1
  else if (CTYPE == '3') then
    N3 = N3+1
  else if (CTYPE == 'S') then
    NS = NS+1
  else if (CTYPE == 'D') then
    ND = ND+1
  else
    write(u6,*) 'TPSTR2ORB_SYM: unknown type index character '//CTYPE
    call AbEnd()
  end if
end do

end subroutine tpstr2orb_sym
