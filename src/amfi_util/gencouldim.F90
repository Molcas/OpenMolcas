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

subroutine gencoulDIM(l1,l2,l3,l4,makemean,icont4)
!bs SUBROUTINE to calculate the dimemsion of the radial integral
!bs arrays. BASICALLY GENCOUL WITHOUT EXPLICIT INTEGRAL CALCULATION
!bs integrals for the four angular momenta l1-l4

use AMFI_global, only: Lblocks, Lfirst, Llast, Lstarter, Lvalues, nblock, ncontrac
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: l1, l2, l3, l4
logical(kind=iwp), intent(in) :: makemean
integer(kind=iwp), intent(out) :: icont4
integer(kind=iwp) :: incl1, incl3, Lanf, Lend, nanz

!bs first of all, this routine determines, for which L
!bs values the radial integrals have to be solved
!bs initialize the number of blocks for the different
!bs l-combinations
!bs no (ss|ss) contributions
if ((l1 == 0) .and. (l2 == 0) .and. (l3 == 0) .and. (l4 == 0)) return
! no integrals for <ss|ss>
if (makemean) then
  nblock = 1  ! sp sp are the first, so the first block
  Lstarter(1) = 1
else
  call SysAbendMsg('gencoulDIM','only mean-field with this version',' ')
end if
!bs keep track of L-values for later purposes
Lvalues(1) = l1
Lvalues(2) = l2
Lvalues(3) = l3
Lvalues(4) = l4
!bs now nanz is given the new value
nanz = ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)

!bs prepare the powers needed for cfunctx
!
! There are seven different CASES of integrals following
!  (   A  --  C)
!
! The structure is the same for all cases, therefore comments can be found only on case A

!bs ####################################################################
!bs   the (+2) cases          CASE A
!bs ####################################################################
incl1 = 1  !  Those increments define the case
incl3 = 1
!bs determine the possible L-values for the integrals by checking for triangular equation

call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)

!bs returns first and last L-values (Lanf,Lend), for which
!bs radial integrals have to be calculated
if (Lend-Lanf >= 0) then
  !bs if there are blocks
  Lblocks(1) = (Lend-Lanf)/2+1 ! L increases in steps of 2, due to parity conservation
  Lfirst(1) = Lanf
  Llast(1) = Lend
else
  Lblocks(1) = 0
end if
!bs ####################################################################
!bs   the (0) cases         CASE  B
!bs ####################################################################
incl1 = 0
incl3 = 0
call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
if (Lend-Lanf >= 0) then
  Lblocks(2) = (Lend-Lanf)/2+1
  Lfirst(2) = Lanf
  Llast(2) = Lend
  Lblocks(3) = (Lend-Lanf)/2+1
  Lfirst(3) = Lanf
  Llast(3) = Lend
else
  Lblocks(2) = 0
  Lblocks(3) = 0
end if
Lstarter(2) = Lstarter(1)+nanz*Lblocks(1)
Lstarter(3) = Lstarter(2)+nanz*Lblocks(2)
!bs ####################################################################
!bs   the (-2) cases      CASE C
!bs ####################################################################
if ((l1 == 0) .or. (l3 == 0)) then
  Lblocks(4) = 0
else
  incl1 = -1
  incl3 = -1
  call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
  if (Lend-Lanf >= 0) then
    Lblocks(4) = (Lend-Lanf)/2+1
    Lfirst(4) = Lanf
    Llast(4) = Lend
  else
    Lblocks(4) = 0
  end if
end if
Lstarter(4) = Lstarter(3)+nanz*Lblocks(3)

!BS now the whole purpose of this routine

icont4 = Lstarter(4)+nanz*Lblocks(4)

return

end subroutine gencoulDIM
