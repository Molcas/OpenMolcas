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

function nMemAm(nShBF,nIrrep,nS,jS,iOffA,Out_of_Core)

use Index_Functions, only: nTri_Elem
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nMemAm
integer(kind=iwp), intent(in) :: nIrrep, nS, nShBf(0:nIrrep-1,nS), jS
integer(kind=iwp), intent(out) :: iOffA(4,0:nIrrep-1)
logical(kind=iwp), intent(in) :: Out_of_Core
integer(kind=iwp) :: iIrrep, lS, ni, nj, nl, nn

!if (Out_Of_Core) then
!  Only a subblock of the upper triangular storage.
!else
!  The whole A matrix is in core! Upper triangular storage.
!endif

nMemAm = 0
do iIrrep=0,nIrrep-1
  ! Pick up the number of basis functions of shell jS in this irrep
  nj = nShBf(iIrrep,jS)
  ! Compute the number of basis functions up to and excluding this shell
  ni = 0
  do lS=1,jS-1
    ni = ni+nShBf(iIrrep,lS)
  end do
  nl = ni+nj
  ! Offset to where this block starts
  if (Out_Of_Core) then
    iOffA(1,iIrrep) = nMemAm
  else
    iOffA(1,iIrrep) = nMemAm+nTri_Elem(ni)
  end if
  ! # of basis functions for the jS shell
  iOffA(2,iIrrep) = nj
  ! # of basis functions up to and including jS shell
  iOffA(4,iIrrep) = nl
  ! Update nMemAm with the size of the whole block for this irrep.
  do lS=jS+1,nS
    nl = nl+nShBf(iIrrep,lS)
  end do
  if (Out_Of_Core) then
    nn = nl-nj
    nMemAm = nMemAm+nTri_Elem(nl)-nTri_Elem(nn)
  else
    nMemAm = nMemAm+nTri_Elem(nl)
  end if
end do

end function nMemAm
