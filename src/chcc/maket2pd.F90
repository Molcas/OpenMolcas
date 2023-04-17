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

use chcc_global, only: DimGrpa, DimSGrpa, no
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: T2p(1), Tau(1)
integer(kind=iwp) :: aGrp, aSGrp
integer(kind=iwp) :: dimabp, dimap, dimapp, dimi, dimij

dimi = no
dimij = no*(no+1)/2

dimap = DimGrpa(aGrp)
dimabp = dimap*(dimap+1)/2

dimapp = DimSGrpa(aSGrp)

call makeT2pdHlp(T2p(1),Tau(1),aGrp,aSGrp,dimi,dimij,dimapp,dimabp)

return

end subroutine MakeT2pd
