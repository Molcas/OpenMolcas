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

subroutine GF_Mult(G,F,GF,nDoF)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDoF
real(kind=wp), intent(in) :: G(nDoF**2), F(nDoF**2)
real(kind=wp), intent(out) :: GF(nTri_Elem(nDoF))
integer(kind=iwp) :: ii, ij, ijT, iX, jj, jX
real(kind=wp) :: XMass_i, XMass_j

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the GF-matrix (actually G^(1/2)FG^(1/2))

do iX=1,nDoF
  ii = (iX-1)*nDoF+iX
  XMass_i = sqrt(G(ii))
  do jX=1,iX
    jj = (jX-1)*nDoF+jX
    XMass_j = sqrt(G(jj))
    ij = (jX-1)*nDoF+iX
    ijT = nTri_Elem(iX-1)+jX
    GF(ijT) = XMass_i*XMass_j*F(ij)
  end do
end do
#ifdef _DEBUGPRINT_
call TriPrt('G^(1/2)FG^(1/2)',' ',GF,nDoF)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GF_Mult
