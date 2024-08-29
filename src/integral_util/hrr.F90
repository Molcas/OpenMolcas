!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine HRR(la,lb,A,B,tgt,nPrim,nTrgt,ipIn)
!***********************************************************************
!                                                                      *
! Object: to generate the contracted integrals (one or two-electron)   *
!         from (a+b,0| or (0,a+b| type of integrals. The approach      *
!         implemented here is the same as that one described by Rys,   *
!         Dupuis and King. In this they implemented the method for the *
!         2D-integrals. The method was later implemented to be applied *
!         on the contracted integrals (Head-Gordon & Pople).           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             February '90                                             *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden.                                         *
!             Modified to not use pointers June '91                    *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: la, lb, nPrim, nTrgt
real(kind=wp), intent(in) :: A(3), B(3)
real(kind=wp), intent(inout) :: tgt(nPrim,nTrgt)
integer(kind=iwp), intent(out) :: ipIn
integer(kind=iwp) :: ia, ia1, ia1b, iab, iab1, iaMax, iaMin, ib, ib1, ib1Max, Ind1, ipRslt
real(kind=wp) :: AB(3), ABSqrt
! Statement function for canonical indices
integer(kind=iwp) :: ix, ixyz, iz, nElem
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
Ind1(ixyz,ix,iz) = ixyz*(ixyz+1)*(ixyz+2)/6+(ixyz-ix)*(ixyz-ix+1)/2+iz+1

! Fast exit if HRR will not be applied.
if ((la == 0) .or. (lb == 0)) then
  ipIn = 1
  return
end if

AB(:) = A(:)-B(:)
if (lb > la) AB(:) = -AB(:)
ABSqrt = sqrt(AB(1)**2+AB(2)**2+AB(3)**2)
if (ABSqrt == Zero) then
  call OCHRR(tgt,nPrim,nTrgt,la,lb,ipIn)
else

  ib1Max = min(la,lb)
  iaMin = max(la,lb)
  ipRslt = 0
  do ib1=1,ib1Max
    ib = ib1-1
    iaMax = la+lb-ib1
    do ia=iaMax,iaMin,-1
      ia1 = ia+1
      if (mod(ib1,2) == 0) then
        ! Low loading of target integrals.
        iab1 = nElem(ib1)*(Ind1(ia,ia,0)-Ind1(iaMin,iaMin,0))
        ! High access of source integrals.
        iab = nTrgt-nElem(ib)*(Ind1(iaMax+1,0,iaMax+1)-Ind1(ia,ia,0)+1)
        ia1b = nTrgt-nElem(ib)*(Ind1(iaMax+1,0,iaMax+1)-Ind1(ia1,ia1,0)+1)
      else
        ! High loading of target integrals.
        iab1 = nTrgt-nElem(ib1)*(Ind1(iaMax,0,iaMax)-Ind1(ia,ia,0)+1)
        ! Low access of source integrals.
        iab = nElem(ib)*(Ind1(ia,ia,0)-Ind1(iaMin,iaMin,0))
        ia1b = nElem(ib)*(Ind1(ia1,ia1,0)-Ind1(iaMin,iaMin,0))
      end if
      ipRslt = iab1

      ! Generate this block of integrals with the HRR

      call HRR1(tgt(1,iab1+1),nElem(ia)*nElem(ib1),tgt(1,ia1b+1),nElem(ia1)*nElem(ib),AB,tgt(1,iab+1),nElem(ia)*nElem(ib),ia,ib, &
                ia1,ib1,nPrim,la,lb)

    end do
  end do
  ipIn = nPrim*ipRslt+1
end if

return

end subroutine HRR
