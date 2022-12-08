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
subroutine TranslateDens(Pi,dRho_dr,dPi,l_tanhr,mGrid,nPi,ndRho_dr,nEGrad,DoGrad)

use nq_Grid, only: GradRho, Rho
use nq_pdft, only: dZdR, fta, ftb, ftc, lft, lGGA, OneMZ, OnePZ, Pass1, Pass2, Pass3, RatioA, RhoAB, ThrsFT, ThrsNT, ThrsOMR, &
                   ThrsRho, ZetaA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Five, Six, Eight, Twelve, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid, nPi, ndRho_dr, nEGrad
real(kind=wp), intent(in) :: Pi(nPi,mGrid), dPi(nPi,nEGrad,mGrid)
real(kind=wp), intent(inout) :: dRho_dr(ndRho_dr,mGrid,nEGrad)
logical(kind=iwp), intent(in) :: l_tanhr, DoGrad
integer(kind=iwp) :: iEGrad, iGrid
real(kind=wp) :: Diff1, GraddZdR, GradRatio, GradRatioX, GradRatioY, GradRatioZ, GradZetaX, GradZetaY, GradZetaZ, Rd2RdRho2, &
                 Rd2RdRhodPi, Rd2ZdR2, Rd2ZdRdZ, RdRdPi, RdRdRho, RRatio, TempR, XAdd, YAdd, ZAdd, ZetaX, ZetaY, ZetaZ
real(kind=wp), allocatable :: dRatio(:,:), dRhodx(:), dRhody(:), dRhodz(:), dZeta(:,:), ftx23(:), fty23(:), ftz23(:), &
                              GradRhoAB(:,:), GradRhoX(:,:), GradRhoY(:,:), GradRhoZ(:,:), RatioX(:), RatioY(:), RatioZ(:), &
                              tanhrx(:), tanhry(:), tanhrz(:)
! PassX
! Pass1. Total density is greater than thresRho
! Pass2. Do translation
! Pass3. Do full translation
! if Pass1 is false, Pass2, 3, are both false.
! if Pass1 is true, Pass2 and 3 cannot both be true (can both be
! false).
!***********************************************************************

!***********************************************************************
! calculating total density at each grid
!***********************************************************************
RhoAB(:) = Rho(1,1:mGrid)+Rho(2,1:mGrid)

!***********************************************************************
! calculating x, y, z components of density gradient
!***********************************************************************
if (lGGA) then
  call mma_allocate(dRhodx,mGrid,label='dRhodx')
  call mma_allocate(dRhody,mGrid,label='dRhody')
  call mma_allocate(dRhodz,mGrid,label='dRhodz')
  dRhodx(:) = GradRho(1,1:mGrid)+GradRho(4,1:mGrid)
  dRhody(:) = GradRho(2,1:mGrid)+GradRho(5,1:mGrid)
  dRhodz(:) = GradRho(3,1:mGrid)+GradRho(6,1:mGrid)
end if

!***********************************************************************
! Ratio and Zeta at each grid point
!***********************************************************************
ZetaA(:) = Zero
RatioA(:) = Zero
dZdR(:) = Zero
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
OnePZ(:) = Half*(One+ZetaA)
OneMZ(:) = Half*(One-ZetaA)

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
    call mma_allocate(RatioX,mGrid,label='RatioX')
    call mma_allocate(RatioY,mGrid,label='RatioY')
    call mma_allocate(RatioZ,mGrid,label='RatioZ')
    call mma_allocate(ftx23,mGrid,label='ftx23')
    call mma_allocate(fty23,mGrid,label='fty23')
    call mma_allocate(ftz23,mGrid,label='ftz23')
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
    GradRho(1,1:mGrid) = GradRho(1,1:mGrid)+ftx23
    GradRho(2,1:mGrid) = GradRho(2,1:mGrid)+fty23
    GradRho(3,1:mGrid) = GradRho(3,1:mGrid)+ftz23
    GradRho(4,1:mGrid) = GradRho(4,1:mGrid)-ftx23
    GradRho(5,1:mGrid) = GradRho(5,1:mGrid)-fty23
    GradRho(6,1:mGrid) = GradRho(6,1:mGrid)-ftz23
    call mma_deallocate(ftx23)
    call mma_deallocate(fty23)
    call mma_deallocate(ftz23)
  end if
end if

!********************************************************************
! Additional terms in the tanh translation
!********************************************************************
if (l_tanhr) then
  call mma_allocate(tanhrx,mGrid,label='tanhrx')
  call mma_allocate(tanhry,mGrid,label='tanhry')
  call mma_allocate(tanhrz,mGrid,label='tanhrz')
  tanhrx(:) = Zero
  tanhry(:) = Zero
  tanhrz(:) = Zero
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      RRatio = RatioA(iGrid)
      TempR = Four*Pi(1,iGrid)/RhoAB(iGrid)
      TanhrX(iGrid) = (RRatio**2-One)*(Pi(2,iGrid)-(dRhodX(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
      TanhrY(iGrid) = (RRatio**2-One)*(Pi(3,iGrid)-(dRhodY(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
      TanhrZ(iGrid) = (RRatio**2-One)*(Pi(4,iGrid)-(dRhodZ(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
    end if
  end do
  GradRho(1,1:mGrid) = GradRho(1,1:mGrid)+TanhrX
  GradRho(2,1:mGrid) = GradRho(2,1:mGrid)+TanhrY
  GradRho(3,1:mGrid) = GradRho(3,1:mGrid)+TanhrZ
  GradRho(4,1:mGrid) = GradRho(4,1:mGrid)-TanhrX
  GradRho(5,1:mGrid) = GradRho(5,1:mGrid)-TanhrY
  GradRho(6,1:mGrid) = GradRho(6,1:mGrid)-TanhrZ
  call mma_deallocate(tanhrx)
  call mma_deallocate(tanhry)
  call mma_deallocate(tanhrz)
end if

!***********************************************************************
! calculating terms needed in gradient calculation
!***********************************************************************
! if not doing gradient, code ends here

if (DoGrad) then
  !*********************************************************************
  ! calculating density gradient wrt geometrical changes
  !*********************************************************************
  call mma_allocate(GradRhoAB,mGrid,nEGrad,label='GradRhoAB')
  GradRhoAB(:,:) = dRho_dr(1,:,:)+dRho_dr(2,:,:)

  if (lGGA) then
    call mma_allocate(GradRhoX,mGrid,nEGrad,label='GradRhoX')
    call mma_allocate(GradRhoY,mGrid,nEGrad,label='GradRhoY')
    call mma_allocate(GradRhoZ,mGrid,nEGrad,label='GradRhoZ')
    GradRhoX(:,:) = dRho_dr(3,:,:)+dRho_dr(6,:,:)
    GradRhoY(:,:) = dRho_dr(4,:,:)+dRho_dr(7,:,:)
    GradRhoZ(:,:) = dRho_dr(5,:,:)+dRho_dr(8,:,:)
  end if

  !***********************************************************************
  ! dRatio and dZeta at each grid point
  !***********************************************************************
  call mma_allocate(dRatio,mGrid,nEGrad,label='dRatio')
  call mma_allocate(dZeta,mGrid,nEGrad,label='dZeta')
  ! Calculate dRatio
  dRatio(:,:) = Zero
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      do iEGrad=1,nEGrad
        dRatio(iGrid,iEGrad) = Four*dPi(1,iEGrad,iGrid)/(RhoAB(iGrid)**2)- &
                               Eight*Pi(1,iGrid)*GradRhoAB(iGrid,iEGrad)/(RhoAB(iGrid)**3)
      end do
    end if
  end do
  ! Calculate dZeta
  do iGrid=1,mGrid
    dZeta(iGrid,:) = dZdR(iGrid)*dRatio(iGrid,:)
  end do

  do iEGrad=1,nEGrad
    do iGrid=1,mGrid
      if (Pass1(iGrid)) then
        dRho_dr(1,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoAB(iGrid,iEGrad)+Half*dZeta(iGrid,iEGrad)*RhoAB(iGrid)
        dRho_dr(2,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoAB(iGrid,iEGrad)-Half*dZeta(iGrid,iEGrad)*RhoAB(iGrid)
      end if
    end do
  end do

  if (lGGA) then
    do iEGrad=1,nEGrad
      do iGrid=1,mGrid
        if (Pass1(iGrid)) then
          dRho_dr(3,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoX(iGrid,iEGrad)+Half*dZeta(iGrid,iEGrad)*dRhodx(iGrid)
          dRho_dr(6,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoX(iGrid,iEGrad)-Half*dZeta(iGrid,iEGrad)*dRhodx(iGrid)
          dRho_dr(4,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoY(iGrid,iEGrad)+Half*dZeta(iGrid,iEGrad)*dRhody(iGrid)
          dRho_dr(7,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoY(iGrid,iEGrad)-Half*dZeta(iGrid,iEGrad)*dRhody(iGrid)
          dRho_dr(5,iGrid,iEGrad) = OnePZ(iGrid)*GradRhoZ(iGrid,iEGrad)+Half*dZeta(iGrid,iEGrad)*dRhodz(iGrid)
          dRho_dr(8,iGrid,iEGrad) = OneMZ(iGrid)*GradRhoZ(iGrid,iEGrad)-Half*dZeta(iGrid,iEGrad)*dRhodz(iGrid)
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
          GradRatio = dRatio(iGrid,iEGrad)
          GradRatioX = (Rd2RdRho2*dRhodX(iGrid)+Rd2RdRhodPi*Pi(2,iGrid))*GradRhoAB(iGrid,iEGrad)+ &
                       Rd2RdRhodPi*dRhodX(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoX(iGrid,iEGrad)+RdRdPi*dPi(2,iEGrad,iGrid)

          GradRatioY = (Rd2RdRho2*dRhodY(iGrid)+Rd2RdRhodPi*Pi(3,iGrid))*GradRhoAB(iGrid,iEGrad)+ &
                       Rd2RdRhodPi*dRhodY(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoY(iGrid,iEGrad)+RdRdPi*dPi(3,iEGrad,iGrid)

          GradRatioZ = (Rd2RdRho2*dRhodZ(iGrid)+Rd2RdRhodPi*Pi(4,iGrid))*GradRhoAB(iGrid,iEGrad)+ &
                       Rd2RdRhodPi*dRhodZ(iGrid)*dPi(1,iEGrad,iGrid)+RdRdRho*GradRhoZ(iGrid,iEGrad)+RdRdPi*dPi(4,iEGrad,iGrid)
          GraddZdR = Zero
          if (Pass2(iGrid)) then
            GraddZdR = Rd2ZdRdZ*dZeta(iGrid,iEGrad)
          else if (Pass3(iGrid)) then
            GraddZdR = Rd2ZdR2*GradRatio
          end if

          GradZetax = GraddZdR*RatioX(iGrid)+dZdR(iGrid)*GradRatioX
          GradZetaY = GraddZdR*RatioY(iGrid)+dZdR(iGrid)*GradRatioY
          GradZetaZ = GraddZdR*RatioZ(iGrid)+dZdR(iGrid)*GradRatioZ

          XAdd = Half*(RhoAB(iGrid)*GradZetaX+ZetaX*GradRhoAB(iGrid,iEGrad))
          YAdd = Half*(RhoAB(iGrid)*GradZetaY+ZetaY*GradRhoAB(iGrid,iEGrad))
          ZAdd = Half*(RhoAB(iGrid)*GradZetaZ+ZetaZ*GradRhoAB(iGrid,iEGrad))

          dRho_dr(3,iGrid,iEGrad) = dRho_dr(3,iGrid,iEGrad)+XAdd
          dRho_dr(6,iGrid,iEGrad) = dRho_dr(6,iGrid,iEGrad)-XAdd
          dRho_dr(4,iGrid,iEGrad) = dRho_dr(4,iGrid,iEGrad)+YAdd
          dRho_dr(7,iGrid,iEGrad) = dRho_dr(7,iGrid,iEGrad)-YAdd
          dRho_dr(5,iGrid,iEGrad) = dRho_dr(5,iGrid,iEGrad)+ZAdd
          dRho_dr(8,iGrid,iEGrad) = dRho_dr(8,iGrid,iEGrad)-ZAdd
        end do
      end do
    end if
    call mma_deallocate(GradRhoX)
    call mma_deallocate(GradRhoY)
    call mma_deallocate(GradRhoZ)
  end if
  call mma_deallocate(GradRhoAB)
  call mma_deallocate(dRatio)
  call mma_deallocate(dZeta)
end if

if (lGGA) then
  call mma_deallocate(dRhodx)
  call mma_deallocate(dRhody)
  call mma_deallocate(dRhodz)
  if (lft) then
    call mma_deallocate(RatioX)
    call mma_deallocate(RatioY)
    call mma_deallocate(RatioZ)
  end if
end if

return

end subroutine TranslateDens
