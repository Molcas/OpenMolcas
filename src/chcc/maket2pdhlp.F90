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

subroutine makeT2pdHlp(T2p,Tau,aGrp,aSGrp,dimi,dimij,dimapp,dimabp)
! this routine does:
! define T2(+)((aa)",ij) = Tau((ab)',i,j) + Tau((ab)',j,i)
! for the case: aGrp=bGrp and aSGrp=bSGrp
!
! parameter description:
! T2p   - array for T2+ (O)
! Tau   - array for Tau (I)
! xGrp  - Groups of a',b' (I)
! xSGrp - SubGroups of a",b" (I)
! dimx  - Dimension of i,(i>=j),a",a',(a>=b)' (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpa, GrpaLow
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, aSGrp, dimi, dimij, dimapp, dimabp
real(kind=wp), intent(out) :: T2p(dimapp,dimij)
real(kind=wp), intent(in) :: Tau(dimabp,dimi,dimi)
integer(kind=iwp) :: abp, ap, app, appAdd, i, ij, j

!1 def appAdd

appAdd = 0
if (aSGrp /= Grpalow(aGrp)) then
  do i=Grpalow(aGrp),aSGrp-1
    appAdd = appAdd+DimSGrpa(i)
  end do
end if

!2 define T2+(aa",ij)

ij = 0
do i=1,dimi
  do j=1,i
    ij = ij+1
    do app=1,dimapp
      ap = appAdd+app
      abp = nTri_Elem(ap)
      T2p(app,ij) = Tau(abp,i,j)+Tau(abp,j,i)
    end do
  end do
end do

T2p(:,:) = Half*T2p(:,:)

return

end subroutine makeT2pdHlp
