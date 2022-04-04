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

subroutine ClasClas(iCNum,nClas,FP,GP,DT,FI,Dist,DistIm,Elene,Edisp,Exrep,E2Die,ExDie)

use qmstat_global, only: CordIm, Cordst, Disp, Exdt1, ExdTal, nCent, nCha, nPart, nPol, Pol, QIm, QSta
use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCNum, nClas
real(kind=wp), intent(out) :: FP(3,nPol*nPart), GP(3,nPol*nPart), DT(3,nPol*nPart), FI(3,nPol*nPart), &
                              Dist(nCent,nCent,nTri_Elem(nClas-1)), DistIm(nCent,nClas,nCent,nClas), Elene, Edisp, Exrep, E2Die, &
                              ExDie
integer(kind=iwp) :: i, ii, ij, Ind, Ind1, indF, IndMa, indR, j, jj, k, l, nSize
real(kind=wp) :: Adisp, aLim, Dampfunk, Epoll, F, Q, Q1, Q2, r, r3, ri, Sum1, Sum2, Sum3, Sum4, X, Y, Z
real(kind=wp), parameter :: Const = 2.2677_wp, ExLim = Ten ! What is Const?
real(kind=wp), external :: ExNemo

!----------------------------------------------------------------------*
! Compute the distance matrices between the classical centers and the  *
! classical image centers.                                             *
!----------------------------------------------------------------------*
Adisp = Disp(1,2)
nSize = nTri_Elem(nClas-1)
Ind = 0
do ii=iCNum+2,nPart
  do jj=iCNum+1,ii-1
    Ind = Ind+1
    do ij=1,nCent
      i = (ii-1)*nCent+ij
      do k=1,nCent
        j = (jj-1)*nCent+k
        r = (Cordst(1,i)-Cordst(1,j))**2+(Cordst(2,i)-Cordst(2,j))**2+(Cordst(3,i)-Cordst(3,j))**2
        Dist(k,ij,Ind) = One/sqrt(r)
      end do
    end do
  end do
end do
do ii=iCNum+1,nPart
  do i=1,nCent
    k = i+(ii-1)*nCent
    do jj=iCNum+1,nPart
      do j=1,nCent
        l = j+(jj-1)*nCent
        r = (CordIm(1,k)-Cordst(1,l))**2+(CordIm(2,k)-Cordst(2,l))**2+(CordIm(3,k)-Cordst(3,l))**2
        DistIm(j,jj-iCNum,i,ii-iCNum) = One/sqrt(r)
      end do
    end do
  end do
end do
!----------------------------------------------------------------------*
! Compute the pairwise interaction between the solvent. Classical all  *
! the  way... early NEMO all the way.                                  *
!----------------------------------------------------------------------*
aLim = One/ExLim
! The electrostatic part

! This loop ONLY works for the early Nemo model of water.
! If the solvent model is changed this loop must be rewritten.
Sum1 = Zero
Sum2 = Zero
Sum3 = Zero
Sum4 = Zero
do i=1,nSize
  Sum1 = Sum1+Dist(2,2,i)*Qsta(1)*Qsta(1) !H-H
  Sum2 = Sum2+Dist(3,2,i)*Qsta(1)*Qsta(2) !H-H
  Sum3 = Sum3+Dist(4,2,i)*Qsta(1)*Qsta(3) !H-V
  Sum4 = Sum4+Dist(5,2,i)*Qsta(1)*Qsta(4) !H-V
  Sum1 = Sum1+Dist(2,3,i)*Qsta(2)*Qsta(1) !H-H
  Sum2 = Sum2+Dist(3,3,i)*Qsta(2)*Qsta(2) !H-H
  Sum3 = Sum3+Dist(4,3,i)*Qsta(2)*Qsta(3) !H-V
  Sum4 = Sum4+Dist(5,3,i)*Qsta(2)*Qsta(4) !H-V
  Sum1 = Sum1+Dist(4,4,i)*Qsta(3)*Qsta(3) !V-V
  Sum2 = Sum2+Dist(5,4,i)*Qsta(3)*Qsta(4) !V-V
  Sum3 = Sum3+Dist(2,4,i)*Qsta(3)*Qsta(1) !V-H
  Sum4 = Sum4+Dist(3,4,i)*Qsta(3)*Qsta(2) !V-H
  Sum1 = Sum1+Dist(4,5,i)*Qsta(4)*Qsta(3) !V-V
  Sum2 = Sum2+Dist(5,5,i)*Qsta(4)*Qsta(4) !V-V
  Sum3 = Sum3+Dist(2,5,i)*Qsta(4)*Qsta(1) !V-H
  Sum4 = Sum4+Dist(3,5,i)*Qsta(4)*Qsta(2) !V-H
end do
Elene = Sum1+Sum2+Sum3+Sum4

Sum1 = Zero
Sum2 = Zero
Sum3 = Zero
Sum4 = Zero
! The dispersion, now with damping.
do i=1,nSize
  DampFunk = One-exp(-One/(Dist(1,1,i)*Const))**4
  Sum1 = Sum1+Dist(1,1,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(2,1,i)*Const))**4
  Sum2 = Sum2+Dist(2,1,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(3,1,i)*Const))**4
  Sum2 = Sum2+Dist(3,1,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(1,2,i)*Const))**4
  Sum2 = Sum2+Dist(1,2,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(2,2,i)*Const))**4
  Sum3 = Sum3+Dist(2,2,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(3,2,i)*Const))**4
  Sum3 = Sum3+Dist(3,2,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(1,3,i)*Const))**4
  Sum2 = Sum2+Dist(1,3,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(2,3,i)*Const))**4
  Sum3 = Sum3+Dist(2,3,i)**6*DampFunk
  DampFunk = One-exp(-One/(Dist(3,3,i)*Const))**4
  Sum3 = Sum3+Dist(3,3,i)**6*DampFunk
end do
Edisp = Sum1*Disp(1,1)+Sum2*Disp(1,2)+Sum3*Disp(2,2)
!The exchange repulsion
Exrep = Zero
do i=1,nSize
  if (Dist(1,1,i) > aLim) Exrep = Exrep+ExNemo(1,1,Dist(1,1,i))
  if (Dist(2,1,i) > aLim) Exrep = Exrep+ExNemo(1,2,Dist(2,1,i))
  if (Dist(3,1,i) > aLim) Exrep = Exrep+ExNemo(1,2,Dist(3,1,i))
  if (Dist(1,2,i) > aLim) Exrep = Exrep+ExNemo(1,2,Dist(1,2,i))
  if (Dist(2,2,i) > aLim) Exrep = Exrep+ExNemo(2,2,Dist(2,2,i))
  if (Dist(3,2,i) > aLim) Exrep = Exrep+ExNemo(2,2,Dist(3,2,i))
  if (Dist(1,3,i) > aLim) Exrep = Exrep+ExNemo(1,2,Dist(1,3,i))
  if (Dist(2,3,i) > aLim) Exrep = Exrep+ExNemo(2,2,Dist(2,3,i))
  if (Dist(3,3,i) > aLim) Exrep = Exrep+ExNemo(2,2,Dist(3,3,i))
end do
!----------------------------------------------------------------------*
! Compute pair-wise interaction with image charges.                    *
!----------------------------------------------------------------------*
Sum1 = Zero
Sum2 = Zero
do i=iCNum+1,nPart
  do j=nCent-nCha+1,nCent !Only count over charged centers.
    Q1 = QIm((i-1)*nCent+j) !The image charge.
    do k=nCent-nCha+1,nCent
      Q2 = QSta(k-nCent+nCha)
      ! Here is the electrostatic interaction computed.
      ! Observe the difference with the real charges,
      ! since here interaction between ALL real-image charge pair is computed.
      do l=iCNum+1,nPart
        Sum1 = Sum1+Q1*Q2*DistIm(k,l-iCNum,j,i-iCNum)
        !                                                                      should this be Const=2.2677 ?
        Sum1 = Sum1-Adisp*DistIm(k,l-iCNum,j,i-iCNum)**6*(One-exp(-One/(DistIm(k,l-iCNum,j,i-iCNum)*2.9677_wp)**6))
      end do
    end do
  end do
  ! Include a repulsion with the boundary to prevent the waters to merge
  ! into the dielectric continuum. Its construction is such that the
  ! repulsion only is between the particle and the image of the particle,
  ! no other repulsion over the boundary.
  Sum2 = Sum2+ExNemo(1,1,DistIm(1,i-iCNum,1,i-iCNum))*Exdt1
  Sum2 = Sum2+ExNemo(1,2,DistIm(2,i-iCNum,2,i-iCNum))
  Sum2 = Sum2+ExNemo(1,2,DistIm(3,i-iCNum,3,i-iCNum))
end do
! The half is added since what we actually have
! computed is the interaction between charge and a part
! of its reaction field (recall:0.5*q*fi_q).
E2Die = sum1*Half
EXDie = sum2*Half*ExdTal
!----------------------------------------------------------------------*
! Compute the static electric field on the polarizabilities and obtain *
! initial guess of induced dipoles.                                    *
!----------------------------------------------------------------------*
IndMa = nPol*nPart
! Explanation: FP-field GP plus reaction field,
!              GP-field from real charges on polarizable centers,
!              FI-induced field.
FP(:,:) = Zero
GP(:,:) = Zero
DT(:,:) = Zero
FI(:,:) = Zero
! Real centers: The field at the polarizabilities - no reaction field.
Ind = 0
do ii=iCNum+2,nPart
  do jj=iCNum+1,ii-1
    Ind = Ind+1
    do ij=1,nCent
      i = (ii-1)*nCent+ij
      do l=1,nCent
        j = (jj-1)*nCent+l
        X = Cordst(1,i)-Cordst(1,j)
        Y = Cordst(2,i)-Cordst(2,j)
        Z = Cordst(3,i)-Cordst(3,j)
        ri = Dist(l,ij,Ind)**3
        if ((ij > nCent-nCha) .and. (l <= nPol)) then
          ! Given that ij is
          ! counting on centers with charges and
          ! l on polarizable centers, then compute
          ! the field from charge ij on center l.
          Q1 = Qsta(ij-nCent+nCha)
          Ind1 = (jj-1)*nPol+l
          GP(1,Ind1) = GP(1,Ind1)+x*Q1*ri
          GP(2,Ind1) = GP(2,Ind1)+y*Q1*ri
          GP(3,Ind1) = GP(3,Ind1)+z*Q1*ri
        end if
        if ((l > nCent-nCha) .and. (ij <= nPol)) then
          ! If ij is
          ! on center with polarizability and l is on
          ! center with charge, then compute the field
          ! from charge l on center ij.
          Q2 = Qsta(l-nCent+nCha)
          Ind1 = (ii-1)*nPol+ij
          GP(1,Ind1) = GP(1,Ind1)-x*Q2*ri
          GP(2,Ind1) = GP(2,Ind1)-y*Q2*ri
          GP(3,Ind1) = GP(3,Ind1)-z*Q2*ri
        end if
      end do
    end do
  end do
end do
Epoll = Zero
! Compute polarization energy.
! This is only for checking, and will
! not enter the energy expression.
do i=1+nPol*iCNum,IndMa
  k = i-((i-1)/nPol)*nPol
  FP(:,i) = GP(:,i)
  F = (FP(1,i)**2+FP(2,i)**2+FP(3,i)**2)*Pol(k)
  Epoll = Epoll+F
end do
Epoll = -Epoll*Half
! Image centers: The field at the polarizabilities - reaction field to the point charges added.
do jj=iCNum+1,nPart
  do ii=1,nCent
    i = ii+(jj-1)*nCent
    Q = Qim(i)
    do k=1,nPol
      indR = k+iCNum*nCent
      indF = k+iCNum*nPol
      do j=1,nClas
        x = CordIm(1,i)-Cordst(1,indR)
        y = CordIm(2,i)-Cordst(2,indR)
        z = CordIm(3,i)-Cordst(3,indR)
        r3 = DistIm(k,j,ii,jj-iCNum)**3
        FP(1,indF) = FP(1,indF)+x*Q*r3
        FP(2,indF) = FP(2,indF)+y*Q*r3
        FP(3,indF) = FP(3,indF)+z*Q*r3
        indR = indR+nCent
        indF = indF+nPol
      end do
    end do
  end do
end do
! We obtain an initial guess of the induced dipoles on the solvent.
do i=1+nPol*iCnum,IndMa
  k = i-((i-1)/nPol)*nPol
  DT(:,i) = FP(:,i)*Pol(k)
end do

return

end subroutine ClasClas
