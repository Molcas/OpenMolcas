!********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 18, 2021, created this file.               *
! ****************************************************************
subroutine TranslateDens(Pi,dRho_dr,dPi,l_tanhr,nRho,mGrid,nPi,ndRho_dr,nEGrad,DoGrad)

use nq_Grid, only: GradRho, nGradRho, Rho
use nq_pdft, only: dZdR, fta, ftb, ftc, lft, lGGA, OneMZ, OnePZ, Pass1, Pass2, Pass3, RatioA, RhoAB, ThrsFT, ThrsNT, ThrsOMR, &
                   ThrsRho, ZetaA
use Constants, only: Zero, One, Two, Three, Four, Five, Six, Eight, Twelve, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nRho, mGrid, nPi, ndRho_dr, nEGrad
real(kind=wp) :: Pi(nPi,mGrid), dRho_dr(ndRho_dr,mGrid,nEGrad), dPi(nPi,nEGrad,mGrid)
logical(kind=iwp) :: l_tanhr, DoGrad
! Input: nRho mGrid nPi ndRho_dr nEGrad Pi dPi DoGrad l_tanhr
! Input & Output: dRho_dr
integer(kind=iwp) :: iEGrad, iGrid, iOff1, ngragri, nGRho
real(kind=wp) :: Diff1, dRatio(mGrid*nEGrad), dRhodx(mGrid), dRhody(mGrid), dRhodz(mGrid), dZeta(mGrid*nEGrad), ftx23(mGrid), & !IFG
                 fty23(mGrid), ftz23(mGrid), GraddZdR, GradRatio, GradRatioX, GradRatioY, GradRatioZ, GradRhoAB(mGrid*nEGrad), & !IFG
                 GradRhoX(mGrid*nEGrad), GradRhoY(mGrid*nEGrad), GradRhoZ(mGrid*nEGrad), GradZetaX, GradZetaY, GradZetaZ, & !IFG
                 RatioX(mGrid), RatioY(mGrid), RatioZ(mGrid), Rd2RdRho2, Rd2RdRhodPi, Rd2ZdR2, Rd2ZdRdZ, RdRdPi, RdRdRho, RRatio, & !IFG
                 tanhrx(mGrid), tanhry(mGrid), tanhrz(mGrid), TempR, XAdd, YAdd, ZAdd, ZetaX, ZetaY, ZetaZ !IFG
! PassX
! Pass1. Total density is greater than thresRho
! Pass2. Do translation
! Pass3. Do full translation
! if Pass1 is false, Pass2, 3, are both false.
! if Pass1 is true, Pass2 and 3 cannot both be true (can both be
! false).
!***********************************************************************

nGRho = nGradRho
!***********************************************************************
! calculating total density at each grid
!***********************************************************************
call DCopy_(mGrid,Rho(1,1),nRho,RhoAB,1)
call DAXPY_(mGrid,One,Rho(2,1),nRho,RhoAB,1)

!***********************************************************************
! calculating x, y, z components of density gradient
!***********************************************************************
if (lGGA) then
  call DCopy_(mGrid,GradRho(1,1),nGRho,dRhodx,1)
  call DAXPY_(mGrid,One,GradRho(4,1),nGRho,dRhodx,1)
  call DCopy_(mGrid,GradRho(2,1),nGRho,dRhody,1)
  call DAXPY_(mGrid,One,GradRho(5,1),nGRho,dRhody,1)
  call DCopy_(mGrid,GradRho(3,1),nGRho,dRhodz,1)
  call DAXPY_(mGrid,One,GradRho(6,1),nGRho,dRhodz,1)
end if

!***********************************************************************
! Ratio and Zeta at each grid point
!***********************************************************************
call FZero(ZetaA,mGrid)
call FZero(RatioA,mGrid)
call FZero(dZdR,mGrid)
do iGrid=1,mGrid
  Pass1(iGrid) = .false.
  Pass2(iGrid) = .false.
  Pass3(iGrid) = .false.
end do

if (.not. lft) then
  do iGrid=1,mGrid
    if (RhoAB(iGrid) >= ThrsRho) then
      Pass1(iGrid) = .true.
      RRatio = Four*Pi(1,iGrid)/(RhoAB(iGrid)**2)
      RatioA(iGrid) = Rratio
      if (l_tanhr) RRatio = tanh(RRatio)
      if ((One-Rratio) > ThrsOMR) then
        ZetaA(iGrid) = sqrt(One-Rratio)
        Pass2(iGrid) = .true.
        dZdR(iGrid) = -Half/ZetaA(iGrid)
      end if
    end if
  end do
else
  do iGrid=1,mGrid
    if (RhoAB(iGrid) >= ThrsRho) then
      Pass1(iGrid) = .true.
      RRatio = Four*Pi(1,iGrid)/(RhoAB(iGrid)**2)
      RatioA(iGrid) = Rratio
      if (RRatio < ThrsFT) then ! do t-translation
        ZetaA(iGrid) = sqrt(One-Rratio)
        Pass2(iGrid) = .true.
        dZdR(iGrid) = -Half/ZetaA(iGrid)
      else if (RRatio <= ThrsNT) then ! do ft-translation
        Diff1 = RRatio-ThrsNT
        ZetaA(iGrid) = (fta*Diff1**2+ftb*Diff1+ftc)*Diff1**3
        Pass3(iGrid) = .true.
        dZdR(iGrid) = (Five*fta*Diff1**2+Four*ftb*Diff1+Three*ftc)*Diff1**2
      end if
    end if
  end do
end if

!***********************************************************************
! (1 + zeta)/2 and (1 - zeta)/2
!***********************************************************************
call DCopy_(mGrid,[Half],0,OnePZ,1)
call DCopy_(mGrid,[Half],0,OneMZ,1)
call DAXPY_(mGrid,Half,ZetaA,1,OnePZ,1)
call DAXPY_(mGrid,-Half,ZetaA,1,OneMZ,1)

!***********************************************************************
! translating rho_a and rho_b
!***********************************************************************
do iGrid=1,mGrid
  if (Pass1(iGrid)) then
    Rho(1,iGrid) = OnePZ(iGrid)*RhoAB(iGrid)
    Rho(2,iGrid) = OneMZ(iGrid)*RhoAB(iGrid)
  end if
end do

!***********************************************************************
! translating gradient component of rho_a and rho_b
!***********************************************************************
if (lGGA) then
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      GradRho(1,iGrid) = OnePZ(iGrid)*dRhodX(iGrid)
      GradRho(2,iGrid) = OnePZ(iGrid)*dRhodY(iGrid)
      GradRho(3,iGrid) = OnePZ(iGrid)*dRhodZ(iGrid)
      GradRho(4,iGrid) = OneMZ(iGrid)*dRhodX(iGrid)
      GradRho(5,iGrid) = OneMZ(iGrid)*dRhodY(iGrid)
      GradRho(6,iGrid) = OneMZ(iGrid)*dRhodZ(iGrid)
    end if
  end do

  if (lft) then
    do iGrid=1,mGrid
      if (Pass1(iGrid)) then
        RatioX(iGrid) = (Four*Pi(2,iGrid)/RhoAB(iGrid)-Two*RatioA(iGrid)*dRhodX(iGrid))/RhoAB(iGrid)
        RatioY(iGrid) = (Four*Pi(3,iGrid)/RhoAB(iGrid)-Two*RatioA(iGrid)*dRhodY(iGrid))/RhoAB(iGrid)
        RatioZ(iGrid) = (Four*Pi(4,iGrid)/RhoAB(iGrid)-Two*RatioA(iGrid)*dRhodZ(iGrid))/RhoAB(iGrid)
      else
        RatioX(iGrid) = Zero
        RatioY(iGrid) = Zero
        RatioZ(iGrid) = Zero
      end if
      ftx23(iGrid) = Half*RhoAB(iGrid)*dZdR(iGrid)*RatioX(iGrid)
      fty23(iGrid) = Half*RhoAB(iGrid)*dZdR(iGrid)*RatioY(iGrid)
      ftz23(iGrid) = Half*RhoAB(iGrid)*dZdR(iGrid)*RatioZ(iGrid)
    end do
    call DaXpY_(mGrid,One,ftx23,1,GradRho(1,1),6)
    call DaXpY_(mGrid,One,fty23,1,GradRho(2,1),6)
    call DaXpY_(mGrid,One,ftz23,1,GradRho(3,1),6)
    call DaXpY_(mGrid,-One,ftx23,1,GradRho(4,1),6)
    call DaXpY_(mGrid,-One,fty23,1,GradRho(5,1),6)
    call DaXpY_(mGrid,-One,ftz23,1,GradRho(6,1),6)
  end if
end if

!********************************************************************
! Additional terms in the tanh translation
!********************************************************************
if (l_tanhr) then
  call FZero(tanhrx,mGrid)
  call FZero(tanhry,mGrid)
  call FZero(tanhrz,mGrid)
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      RRatio = RatioA(iGrid)
      TempR = Four*Pi(1,iGrid)/RhoAB(iGrid)
      TanhrX(iGrid) = (RRatio**2-One)*(Pi(2,iGrid)-(dRhodX(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
      TanhrY(iGrid) = (RRatio**2-One)*(Pi(3,iGrid)-(dRhodY(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
      TanhrZ(iGrid) = (RRatio**2-One)*(Pi(4,iGrid)-(dRhodZ(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
    end if
  end do
  call DAXPY_(mGrid,One,TanhrX,1,GradRho(1,1),nRho)
  call DAXPY_(mGrid,-One,TanhrX,1,GradRho(4,1),nRho)
  call DAXPY_(mGrid,One,TanhrY,1,GradRho(2,1),nRho)
  call DAXPY_(mGrid,-One,TanhrY,1,GradRho(5,1),nRho)
  call DAXPY_(mGrid,One,TanhrZ,1,GradRho(3,1),nRho)
  call DAXPY_(mGrid,-One,TanhrZ,1,GradRho(6,1),nRho)
end if

!***********************************************************************
! calculating terms needed in gradient calculation
!***********************************************************************
! if not doing gradient, code ends here
if (.not. DoGrad) return
!***********************************************************************
! calculating density gradient wrt geometrical changes
!***********************************************************************
ngragri = mGrid*nEGrad
call DCopy_(ngragri,dRho_dr(1,1,1),ndRho_dr,GradRhoAB,1)
call DAXPY_(ngragri,One,dRho_dr(2,1,1),ndRho_dr,GradRhoAB,1)

if (lGGA) then
  call DCopy_(ngragri,dRho_dr(3,1,1),ndRho_dr,GradRhoX,1)
  call DAXPY_(ngragri,One,dRho_dr(6,1,1),ndRho_dr,GradRhoX,1)
  call DCopy_(ngragri,dRho_dr(4,1,1),ndRho_dr,GradRhoY,1)
  call DAXPY_(ngragri,One,dRho_dr(7,1,1),ndRho_dr,GradRhoY,1)
  call DCopy_(ngragri,dRho_dr(5,1,1),ndRho_dr,GradRhoZ,1)
  call DAXPY_(ngragri,One,dRho_dr(8,1,1),ndRho_dr,GradRhoZ,1)
end if

!***********************************************************************
! dRatio and dZeta at each grid point
!***********************************************************************
! Calculate dRatio
call Fzero(dRatio,nGraGri)
do iGrid=1,mGrid
  if (Pass1(iGrid)) then
    do iEGrad=1,nEGrad
      IOff1 = (iEGrad-1)*mGrid
      dRatio(IOff1+iGrid) = Four*dPi(1,iEGrad,iGrid)/(RhoAB(iGrid)**2)-Eight*Pi(1,iGrid)*GradRhoAB(IOff1+iGrid)/(RhoAB(iGrid)**3)
    end do
  end if
end do
! Calculate dZeta
call Fzero(dZeta,nGraGri)
do iGrid=1,mGrid
  call DAxpy_(nEGrad,dZdR(iGrid),dRatio(iGrid),mGrid,dZeta(iGrid),mGrid)
end do

do iEGrad=1,nEGrad
  IOff1 = (iEGrad-1)*mGrid
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      dRho_dr(1,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoAB(IOff1+iGrid)+Half*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
      dRho_dr(2,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoAB(IOff1+iGrid)-Half*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
    end if
  end do
end do

if (lGGA) then
  do iEGrad=1,nEGrad
    IOff1 = (iEGrad-1)*mGrid
    do iGrid=1,mGrid
      if (Pass1(iGrid)) then
        dRho_dr(3,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoX(IOff1+iGrid)+Half*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
        dRho_dr(6,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoX(IOff1+iGrid)-Half*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
        dRho_dr(4,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoY(IOff1+iGrid)+Half*dZeta(IOFf1+iGrid)*dRhody(iGrid)
        dRho_dr(7,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoY(IOff1+iGrid)-Half*dZeta(IOFf1+iGrid)*dRhody(iGrid)
        dRho_dr(5,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoZ(IOff1+iGrid)+Half*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
        dRho_dr(8,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoZ(IOff1+iGrid)-Half*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
      end if
    end do
  end do
  if (lft) then
    do iGrid=1,mGrid
      if (.not. Pass1(iGrid)) cycle
      if (.not. (Pass2(iGrid) .or. Pass3(iGrid))) cycle
      ZetaX = dZdR(iGrid)*RatioX(iGrid)
      ZetaY = dZdR(iGrid)*RatioY(iGrid)
      ZetaZ = dZdR(iGrid)*RatioZ(iGrid)
      RdRdRho = -Two*RatioA(iGrid)/RhoAB(iGrid)
      RdRdPi = Four/RhoAB(iGrid)**2
      Rd2RdRho2 = -Three*RdRdRho/RhoAB(iGrid)
      Rd2RdRhodPi = -Two*RdRdPi/RhoAB(iGrid)
      Rd2ZdRdZ = Zero
      Rd2ZdR2 = Zero
      if (Pass2(iGrid)) Rd2ZdRdZ = Half/ZetaA(iGrid)**2
      if (Pass3(iGrid)) then
        Diff1 = RatioA(iGrid)-ThrsNT
        Rd2ZdR2 = (20.0_wp*fta*Diff1**2+Twelve*ftb*Diff1+Six*ftc)*Diff1
      end if
      do iEGrad=1,nEGrad
        IOff1 = (iEGrad-1)*mGrid
        GradRatio = dRatio(iOff1+iGrid)
        GradRatioX = (Rd2RdRho2*dRhodX(iGrid)+Rd2RdRhodPi*Pi(2,iGrid))*GradRhoAB(iOff1+iGrid)+ &
                     Rd2RdRhodPi*dRhodX(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoX(iOff1+iGrid)+RdRdPi*dPi(2,iEGrad,iGrid)

        GradRatioY = (Rd2RdRho2*dRhodY(iGrid)+Rd2RdRhodPi*Pi(3,iGrid))*GradRhoAB(iOff1+iGrid)+ &
                     Rd2RdRhodPi*dRhodY(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoY(iOff1+iGrid)+RdRdPi*dPi(3,iEGrad,iGrid)

        GradRatioZ = (Rd2RdRho2*dRhodZ(iGrid)+Rd2RdRhodPi*Pi(4,iGrid))*GradRhoAB(iOff1+iGrid)+ &
                     Rd2RdRhodPi*dRhodZ(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoZ(iOff1+iGrid)+RdRdPi*dPi(4,iEGrad,iGrid)
        GraddZdR = Zero
        if (Pass2(iGrid)) then
          GraddZdR = Rd2ZdRdZ*dZeta(iOff1+iGrid)
        else if (Pass3(iGrid)) then
          GraddZdR = Rd2ZdR2*GradRatio
        end if

        GradZetax = GraddZdR*RatioX(iGrid)+dZdR(iGrid)*GradRatioX
        GradZetaY = GraddZdR*RatioY(iGrid)+dZdR(iGrid)*GradRatioY
        GradZetaZ = GraddZdR*RatioZ(iGrid)+dZdR(iGrid)*GradRatioZ

        XAdd = Half*(RhoAB(iGrid)*GradZetaX+ZetaX*GradRhoAB(iOff1+iGrid))
        YAdd = Half*(RhoAB(iGrid)*GradZetaY+ZetaY*GradRhoAB(iOff1+iGrid))
        ZAdd = Half*(RhoAB(iGrid)*GradZetaZ+ZetaZ*GradRhoAB(iOff1+iGrid))

        dRho_dr(3,iGrid,iEGrad) = dRho_dr(3,iGrid,iEGrad)+XAdd
        dRho_dr(6,iGrid,iEGrad) = dRho_dr(6,iGrid,iEGrad)-XAdd
        dRho_dr(4,iGrid,iEGrad) = dRho_dr(4,iGrid,iEGrad)+YAdd
        dRho_dr(7,iGrid,iEGrad) = dRho_dr(7,iGrid,iEGrad)-YAdd
        dRho_dr(5,iGrid,iEGrad) = dRho_dr(5,iGrid,iEGrad)+ZAdd
        dRho_dr(8,iGrid,iEGrad) = dRho_dr(8,iGrid,iEGrad)-ZAdd
      end do
    end do
  end if
end if

return

end subroutine TranslateDens
