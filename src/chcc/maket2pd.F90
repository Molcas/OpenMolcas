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

subroutine MakeT2pd(T2p,Tau,aGrp,aSGrp)
! this routine does:
! Make T2p((a>b)",i>=j )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
!                     from  Tau((a>=b)',i,j)
!
! parameter description:
! T2p   - T2+ array (O)
! Tau   - Tau array (I)
! xGrp  - Group of a (I)
! xSGrp - SubGroup of a (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpa, DimSGrpa, no
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, aSGrp
real(kind=wp), intent(out) :: T2p(DimSGrpa(aSGrp),nTri_Elem(no))
real(kind=wp), intent(in) :: Tau(nTri_Elem(DimGrpa(aGrp)),no,no)
integer(kind=iwp) :: dimabp, dimap, dimapp, dimi, dimij

dimi = no
dimij = nTri_Elem(no)

dimap = DimGrpa(aGrp)
dimabp = nTri_Elem(dimap)

dimapp = DimSGrpa(aSGrp)

call makeT2pdHlp(T2p,Tau,aGrp,aSGrp,dimi,dimij,dimapp,dimabp)

return

end subroutine MakeT2pd
