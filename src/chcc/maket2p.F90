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

subroutine MakeT2p(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,keyT)
! this routine does:
! Make T2p(i>=j,(a>b)" )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
!                     from  Tau((a>=b)',i,j)
! or Transposed (T(ab",ij)
!
! parameter description:
! T2p   - T2+ array (O)
! Tau   - Tau array (I)
! xGrp  - Group of a,b (I)
! xSGrp - SubGroup of a,b (I)
! keyT  - 0 - make T(ij,ab")
!         1 - make T(ab",ij)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpa, DimSGrpa, no
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: T2p(*)
real(kind=wp), intent(in) :: Tau(*)
integer(kind=iwp), intent(in) :: aGrp, bGrp, aSGrp, bSGrp, keyT
integer(kind=iwp) :: dimabp, dimabpp, dimap, dimapp, dimbp, dimbpp, dimi, dimij

dimi = no
dimij = nTri_Elem(no)

dimap = DimGrpa(aGrp)
dimbp = DimGrpa(bGrp)
if (aGrp == bGrp) then
  dimabp = nTri_Elem(dimap)
else
  dimabp = dimap*dimbp
end if

dimapp = DimSGrpa(aSGrp)
dimbpp = DimSGrpa(bSGrp)
if (aSGrp == bSGrp) then
  dimabpp = nTri_Elem(dimapp-1)
else
  dimabpp = dimapp*dimbpp
end if

if (keyT == 0) then
  ! T+(ij,ab") case

  if (aGrp == bGrp) then
    if (aSGrp == bSGrp) then
      call makeT2pHlp1(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimabpp,dimabp)
    else
      call makeT2pHlp2(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimbpp,dimabp)
    end if
  else
    call makeT2pHlp3(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimbpp,dimap,dimbp)
  end if

else
  ! T+(ab",ij) case

  if (aGrp == bGrp) then
    if (aSGrp == bSGrp) then
      call makeT2ptHlp1(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimabpp,dimabp)
    else
      call makeT2ptHlp2(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimbpp,dimabp)
    end if
  else
    call makeT2ptHlp3(T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,0,dimi,dimij,dimapp,dimbpp,dimap,dimbp)
  end if

end if

return

end subroutine MakeT2p
