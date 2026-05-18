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

subroutine GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,nCMO,OCC,nOCC,GDAO,nGDAO)
! General density matrix, in the sense of using a specified
! but arbitrary range of orbitals in each symmetry, and an
! occupation number.

! The input arguments:
!    NSYM=Nr of symmetry blocks,
!    NBAS(i), i=1..NSYM: Number of basis functions, equal
!       also to total number of orbitals, in each symmetry.
!    CNAT(i),i=1..sum(NBAS(i)**2,i=1..NSYM): The CMO
!       coefficients, stored as square symmetry blocks CN(a,q)
!       where a is basis function, q is orbital index.
!    ISTART(i), i=1..NSYM: Orbital number to start using
!        in each symmetry block, numbered 1,2,..NBAS(i)
!    NUSE(i), i=1..NSYM: How many orbitals to use,
!    OCC(i),i=1..sum(NBAS(i),i=1..NSYM): The natural
!        occupation number of each orbital.
! Output: Array GDAO.

! Computes the AO density matrix, using formula
!  D(a,b) = sum( CN(a,p)*xn(p)*CN(b,p), p=1..n)
! where a,b are basis function indices, p is an MO index,
!  CN() and xn() are CMO coefficients and occupation
! numbers of natural orbitals, and the indices and matrices
! are defined within symmetry blocks.
! The symmetry blocks of D are then stored triangularly
! after each other in the array GDAO.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSYM, NBAS(NSYM), ISTART(NSYM), NUSE(NSYM), nCMO, nOcc, nGDAO
real(kind=wp), intent(in) :: CNAT(nCMO), OCC(nOcc)
real(kind=wp), intent(out) :: GDAO(nGDAO)
integer(kind=iwp) :: IOEND, ICEND, IDAB, ISYM, NB, NW, IW1, IW2, IA, IB, IW
real(kind=wp) :: DAB

IOEND = 0
ICEND = 0
IDAB = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB < 1) cycle
  call DCOPY_((NB*(NB+1))/2,[Zero],0,GDAO(IDAB+1),1)
  NW = NUSE(ISYM)
  if (NW > 0) then
    IW1 = ISTART(ISYM)
    IW2 = IW1-1+NW
    do IA=1,NB
      do IB=1,IA
        IDAB = IDAB+1
        DAB = GDAO(IDAB)
        do IW=IW1,IW2
          DAB = DAB+OCC(IOEND+IW)*CNAT(ICEND+IA+NB*(IW-1))*CNAT(ICEND+IB+NB*(IW-1))
        end do
        GDAO(IDAB) = DAB
      end do
    end do
  else
    IDAB = IDAB+(NB*(NB+1))/2
  end if
  IOEND = IOEND+NB
  ICEND = ICEND+NB**2
end do

end subroutine GDMAT
