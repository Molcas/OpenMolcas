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

subroutine MkT_C78d(T2,Tp,Tm,dimbe,dimbepp,addbepp,no)
! this routine does:
! T2n(be',ga',u,v) <<-
! C7                    T2+(bega",uv)
! C8                    T2-(bega",uv)
! for beGrp=gaGrp, beSGrp=gaSGrp
! N.B. calc only contributions to be',ga' (not ga',be')

use Index_Functions, only: nTri_Elem
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dimbepp, addbepp, no
real(kind=wp), intent(inout) :: T2(dimbe,dimbe,no,no)
real(kind=wp), intent(in) :: Tp(nTri_Elem(dimbepp),nTri_Elem(no)), Tm(nTri_Elem(dimbepp-1),nTri_Elem(no-1))
integer(kind=iwp) :: be, bega, bep, u, uv, v
real(kind=wp) :: fact

!1 Distribute symmetric T2+ on proper positions

uv = 0
do u=1,no
  do v=1,u
    uv = uv+1
    if (u == v) then
      fact = Half
    else
      fact = One
    end if

    ! case be"/=ga"
    bep = addbepp+1
    do be=2,dimbepp
      bega = nTri_Elem(be-1)
      bep = bep+1

      T2(bep,addbepp+1:addbepp+be-1,u,v) = T2(bep,addbepp+1:addbepp+be-1,u,v)+Tp(bega+1:bega+be-1,uv)*fact
      T2(bep,addbepp+1:addbepp+be-1,v,u) = T2(bep,addbepp+1:addbepp+be-1,v,u)+Tp(bega+1:bega+be-1,uv)*fact
      !T2(addbepp+1:addbepp+be-1,bep,u,v) = T2(addbepp+1:addbepp+be-1,bep,u,v)+Tp(bega+1:bega+be-1,uv)*fact
      !T2(addbepp+1:addbepp+be-1,bep,v,u) = T2(addbepp+1:addbepp+be-1,bep,v,u)+Tp(bega+1:bega+be-1,uv)*fact

    end do

    ! case be=ga
    bep = addbepp
    do be=1,dimbepp
      bega = nTri_Elem(be)
      bep = bep+1

      T2(bep,bep,u,v) = T2(bep,bep,u,v)+Tp(bega,uv)*fact
      T2(bep,bep,v,u) = T2(bep,bep,v,u)+Tp(bega,uv)*fact

    end do

  end do
end do

!2 Distribute anti-symmetric T2- on proper positions

uv = 0
do u=2,no
  do v=1,u-1
    uv = uv+1

    bep = addbepp+1
    bega = 0
    do be=2,dimbepp
      bep = bep+1

      T2(bep,addbepp+1:addbepp+be-1,u,v) = T2(bep,addbepp+1:addbepp+be-1,u,v)+Tm(bega+1:bega+be-1,uv)
      T2(bep,addbepp+1:addbepp+be-1,v,u) = T2(bep,addbepp+1:addbepp+be-1,v,u)-Tm(bega+1:bega+be-1,uv)
      !T2(addbepp+1:addbepp+be-1,bep,u,v) = T2(addbepp+1:addbepp+be-1,bep,u,v)-Tm(bega+1:bega+be-1,uv)
      !T2(addbepp+1:addbepp+be-1,bep,v,u) = T2(addbepp+1:addbepp+be-1,bep,v,u)+Tm(bega+1:bega+be-1,uv)

      bega = bega+be-1
    end do

  end do
end do

return

end subroutine MkT_C78d
