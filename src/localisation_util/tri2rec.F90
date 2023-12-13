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

subroutine Tri2Rec(OvlTri,OvlRec,nBas,Debug)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas
real(kind=wp), intent(in) :: OvlTri(*)
real(kind=wp), intent(out) :: OvlRec(nBas,nBas)
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: iBas, ifr, ito, jBas, nlig

ifr = 1
do nlig=1,nBas
  ito = ifr+nlig-1
  OvlRec(1:nlig,nlig) = OvlTri(ifr:ito)
  ifr = ito+1
end do

do iBas=1,nBas
  do jBas=nBas,iBas,-1
    OvlRec(jBas,iBas) = OvlRec(iBas,jBas)
  end do
end do

if (Debug) call RecPrt('OvlRec ',' ',OvlRec,nBas,nBas)

return

end subroutine Tri2Rec
