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

subroutine makeT2ptHlp1(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,dimi,dimij,dimapp,dimabpp,dimabp)
! this routine does:
! define T2(+-)((ab)",ij) = Tau((ab)',i,j) +- Tau((ab)',j,i)
! for the case: aGrp=bGrp and aSGrp=bSGrp
!
! parameter description:
! T2p   - array for T2+ (O)
! Tau   - array for Tau (I)
! xGrp  - Groups of a',b' (I)
! xSGrp - SubGroups of a",b" (I)
! key   - 0 - T2+; 1 = T2- will be produced
! dimx  - Dimension of i,(i>=j),a",(a>(or>=)b)",a',(a>=b)' (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpa, GrpaLow
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, bGrp, aSGrp, bSGrp, key, dimi, dimij, dimapp, dimabpp, dimabp
real(kind=wp), intent(out) :: T2p(dimabpp,dimij)
real(kind=wp), intent(in) :: Tau(dimabp,dimi,dimi)
integer(kind=iwp) :: abp, abpp, ap, app, appAdd, bppAdd, i, ij

!1 def appAdd,bppAdd

appAdd = 0
if (aSGrp /= Grpalow(aGrp)) then
  do i=Grpalow(aGrp),aSGrp-1
    appAdd = appAdd+DimSGrpa(i)
  end do
end if

bppAdd = 0
if (bSGrp /= Grpalow(bGrp)) then
  do i=Grpalow(bGrp),bSGrp-1
    bppAdd = bppAdd+DimSGrpa(i)
  end do
end if

!2 define T2(+,-)

if (key == 0) then

  !2.1 define T2+
  ij = 0
  do i=1,dimi
    abpp = 0
    do app=2,dimapp
      ap = appAdd+app
      abp = nTri_Elem(ap-1)+bppAdd
      T2p(abpp+1:abpp+app-1,ij+1:ij+i) = Tau(abp+1:abp+app-1,i,1:i)+Tau(abp+1:abp+app-1,1:i,i)
      abpp = abpp+app-1
    end do
    ij = ij+i
  end do

else

  !2.2 define T2-
  ij = 0
  do i=2,dimi
    abpp = 0
    do app=2,dimapp
      ap = appAdd+app
      abp = nTri_Elem(ap-1)+bppAdd
      T2p(abpp+1:abpp+app-1,ij+1:ij+i-1) = Tau(abp+1:abp+app-1,i,1:i-1)-Tau(abp+1:abp+app-1,1:i-1,i)
      abpp = abpp+app-1
    end do
    ij = ij+i-1
  end do

end if

T2p(:,:) = Half*T2p(:,:)

return

end subroutine makeT2ptHlp1
