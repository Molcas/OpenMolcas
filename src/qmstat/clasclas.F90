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

subroutine ClasClas(iCNum,iCStart,ncParm,iFP,iGP,iDT,iFI,iDist,iDistIm,Elene,Edisp,Exrep,E2Die,ExDie)

use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iCNum, iCStart, ncParm, iFP(3), iGP(3), iDT(3), iFI(3), iDist, iDistIm
real(kind=wp) :: Elene, Edisp, Exrep, E2Die, ExDie
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ii, ij, Inc, Inc2, Ind, Ind1, indF, IndMa, indR, indSep, j, jj, Jnd, k, l, nClas, nSize, nSizeIm
real(kind=wp) :: Adisp, aLim, Dampfunk, Epoll, F, Q, Q1, Q2, r, r3, ri, Sum1, Sum2, Sum3, Sum4, Sum5, X, Y, Z
character(len=20) :: MemLaaaa, Memlaaab, Memlaabe, Memlabel
character(len=2) :: ChCo
real(kind=wp), parameter :: Const = 2.2677_wp, ExLim = Ten ! What is Const?
real(kind=wp), external :: ExNemo

!----------------------------------------------------------------------*
! Compute the distance matrices between the classical centers and the  *
! classical image centers.                                             *
!----------------------------------------------------------------------*
nClas = nPart-iCNum
Adisp = Disp(1,2)
nSize = (nClas*(nClas-1)/2)*(nCent**2) !Get memory
call GetMem('DistMat','Allo','Real',iDist,nSize)
nSizeIm = (nClas*nCent)**2
call GetMem('DistMatIm','Allo','Real',iDistIm,nSizeIm)
Ind = 0
do ii=iCNum+2,nPart
  do jj=iCNum+1,ii-1
    do ij=1,nCent
      i = (ii-1)*nCent+ij
      do k=1,nCent
        Ind = Ind+1
        j = (jj-1)*nCent+k
        r = 0
        do l=1,3
          r = (Cordst(i,l)-Cordst(j,l))**2+r
        end do
        Work(iDist+Ind-1) = One/sqrt(r)
      end do
    end do
  end do
end do
Jnd = 0
do i=iCStart,nCent*nPart
  do j=iCStart,nCent*nPart
    Jnd = Jnd+1
    r = 0
    do k=1,3
      r = (CordIm(i,k)-Cordst(j,k))**2+r
    end do
    Work(iDistIm+Jnd-1) = One/sqrt(r)
  end do
end do
!----------------------------------------------------------------------*
! Compute the pairwise interaction between the solvent. Classical all  *
! the  way... early NEMO all the way.                                  *
!----------------------------------------------------------------------*
Elene = Zero
Edisp = Zero
Exrep = Zero
aLim = One/ExLim
Sum1 = Zero
Sum2 = Zero
Sum3 = Zero
Sum4 = Zero
! The electrostatic part

! This loop ONLY works for the early Nemo model of water.
! If the solvent model is changed this loop must be rewritten.
do i=1,nSize,nCent**2
  Sum1 = Sum1+Work(iDist+i-1+6)*Qsta(1)*Qsta(1) !H-H
  Sum2 = Sum2+Work(iDist+i-1+7)*Qsta(1)*Qsta(2) !H-H
  Sum3 = Sum3+Work(iDist+i-1+8)*Qsta(1)*Qsta(3) !H-V
  Sum4 = Sum4+Work(iDist+i-1+9)*Qsta(1)*Qsta(4) !H-V
  Sum1 = Sum1+Work(iDist+i-1+11)*Qsta(2)*Qsta(1) !H-H
  Sum2 = Sum2+Work(iDist+i-1+12)*Qsta(2)*Qsta(2) !H-H
  Sum3 = Sum3+Work(iDist+i-1+13)*Qsta(2)*Qsta(3) !H-V
  Sum4 = Sum4+Work(iDist+i-1+14)*Qsta(2)*Qsta(4) !H-V
  Sum1 = Sum1+Work(iDist+i-1+18)*Qsta(3)*Qsta(3) !V-V
  Sum2 = Sum2+Work(iDist+i-1+19)*Qsta(3)*Qsta(4) !V-V
  Sum3 = Sum3+Work(iDist+i-1+16)*Qsta(3)*Qsta(1) !V-H
  Sum4 = Sum4+Work(iDist+i-1+17)*Qsta(3)*Qsta(2) !V-H
  Sum1 = Sum1+Work(iDist+i-1+23)*Qsta(4)*Qsta(3) !V-V
  Sum2 = Sum2+Work(iDist+i-1+24)*Qsta(4)*Qsta(4) !V-V
  Sum3 = Sum3+Work(iDist+i-1+21)*Qsta(4)*Qsta(1) !V-H
  Sum4 = Sum4+Work(iDist+i-1+22)*Qsta(4)*Qsta(2) !V-H
end do
Elene = Sum1+Sum2+Sum3+Sum4

Sum1 = Zero
Sum2 = Zero
Sum3 = Zero
Sum4 = Zero
Sum5 = Zero
! The dispersion, now with damping.
do i=1,nSize,nCent**2
  DampFunk = One-exp(-One/(Work(iDist+i-1)*Const))**4
  Sum1 = Sum1+Work(iDist+i-1)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+1)*Const))**4
  Sum2 = Sum2+Work(iDist+i-1+1)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+2)*Const))**4
  Sum3 = Sum3+Work(iDist+i-1+2)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+11)*Const))**4
  Sum4 = Sum4+Work(iDist+i-1+11)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+7)*Const))**4
  Sum5 = Sum5+Work(iDist+i-1+7)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+5)*Const))**4
  Sum2 = Sum2+Work(iDist+i-1+5)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+10)*Const))**4
  Sum3 = Sum3+Work(iDist+i-1+10)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+6)*Const))**4
  Sum4 = Sum4+Work(iDist+i-1+6)**6*DampFunk
  DampFunk = One-exp(-One/(Work(iDist+i-1+12)*Const))**4
  Sum5 = Sum5+Work(iDist+i-1+12)**6*DampFunk
end do
Edisp = Sum1*Disp(1,1)+(Sum2+Sum3)*Disp(1,2)+(Sum4+Sum5)*Disp(2,2)
!The exchange repulsion
do i=1,nSize,nCent**2
  if (Work(iDist+i-1) > aLim) Exrep = Exrep+ExNemo(1,1,Work(iDist+i-1))
  if (Work(iDist+i-1+1) > aLim) Exrep = Exrep+ExNemo(1,2,Work(iDist+i-1+1))
  if (Work(iDist+i-1+2) > aLim) Exrep = Exrep+ExNemo(1,2,Work(iDist+i-1+2))
  if (Work(iDist+i-1+5) > aLim) Exrep = Exrep+ExNemo(1,2,Work(iDist+i-1+5))
  if (Work(iDist+i-1+6) > aLim) Exrep = Exrep+ExNemo(2,2,Work(iDist+i-1+6))
  if (Work(iDist+i-1+7) > aLim) Exrep = Exrep+ExNemo(2,2,Work(iDist+i-1+7))
  if (Work(iDist+i-1+10) > aLim) Exrep = Exrep+ExNemo(1,2,Work(iDist+i-1+10))
  if (Work(iDist+i-1+11) > aLim) Exrep = Exrep+ExNemo(2,2,Work(iDist+i-1+11))
  if (Work(iDist+i-1+12) > aLim) Exrep = Exrep+ExNemo(2,2,Work(iDist+i-1+12))
end do
!----------------------------------------------------------------------*
! Compute pair-wise interaction with image charges.                    *
!----------------------------------------------------------------------*
Sum1 = Zero
Sum2 = Zero
do i=iCNum+1,nPart
  do j=nCent-nCha+1,nCent !Only count over charged centers.
    Q1 = QIm((i-1)*nCent+j) !The image charge.
    Inc = ncParm*nCent*(i-(iCNum+1))+(j-1)*ncParm !Counting elements.
    do k=nCent-nCha+1,nCent
      Inc2 = Inc+k
      Q2 = QSta(k-nCent+nCha)
      ! Here is the electrostatic interaction computed.
      ! Observe the difference with the real charges,
      ! since here interaction between ALL real-image charge pair is computed.
      do l=iCNum+1,nPart
        Sum1 = Sum1+Q1*Q2*Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)
        Sum1 = Sum1-Adisp*Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)**6* &
               (One-exp(-One/(Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)*2.9677_wp)**6)) !should this be Const=2.2677 ?
      end do
    end do
  end do
  ! Include a repulsion with the boundary to prevent the waters to merge
  ! into the dielectric continuum. Its construction is such that the
  ! repulsion only is between the particle and the image of the particle,
  ! no other repulsion over the boundary.
  Sum2 = Sum2+ExNemo(1,2,Work(iDistIm-1+(i-(iCNum+1))*nClas*nCent**2+(nClas+i-(iCNum+1))*nCent+2))
  Sum2 = Sum2+ExNemo(1,2,Work(iDistIm-1+(i-(iCNum+1))*nClas*nCent**2+(2*nClas+i-(iCNum+1))*nCent+3))
  Sum2 = Sum2+ExNemo(1,1,Work(iDistIm-1+(i-(iCNum+1))*nClas*nCent**2+1+(i-(iCNum+1))*nCent))*Exdt1
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
do i=1,3   !Allocate memory
  write(ChCo,'(I2.2)') i
  write(MemLabel,*) 'FP'//ChCo
  write(MemLaabe,*) 'GP'//ChCo
  write(MemLaaab,*) 'DT'//ChCo
  write(MemLaaaa,*) 'FI'//ChCo
  ! Explanation: iFP-field iGP plus reaction field,
  !              iGP-field from real charges on polarizable centers,
  !              iFi-induced field.
  call GetMem(MemLabel,'Allo','Real',iFP(i),IndMa)
  call GetMem(MemLaabe,'Allo','Real',iGP(i),IndMa)
  call GetMem(MemLaaab,'Allo','Real',iDT(i),IndMa)
  call GetMem(MemLaaaa,'Allo','Real',iFi(i),IndMa)
end do
do j=1,3
  do i=0,Indma-1 !Set some zeros
    Work(iFI(j)+i) = Zero
    Work(iGP(j)+i) = Zero
    Work(iDT(j)+i) = Zero
    Work(iFP(j)+i) = Zero
  end do
end do
! Real centers: The field at the polarizabilities - no reaction field.
Ind = 0
do ii=iCNum+2,nPart
  do jj=iCNum+1,ii-1
    do ij=1,nCent
      i = (ii-1)*nCent+ij
      do l=1,nCent
        j = (jj-1)*nCent+l
        Ind = Ind+1
        X = Cordst(i,1)-Cordst(j,1)
        Y = Cordst(i,2)-Cordst(j,2)
        Z = Cordst(i,3)-Cordst(j,3)
        ri = Work(iDist+Ind-1)**3
        if ((ij > nCent-nCha) .and. (l <= nPol)) then
          ! Given that ij is
          ! counting on centers with charges and
          ! l on polarizable centers, then compute
          ! the field from charge ij on center l.
          Q1 = Qsta(ij-nCent+nCha)
          Ind1 = (jj-1)*nPol+l
          Work(iGP(1)+Ind1-1) = Work(iGP(1)+Ind1-1)+x*Q1*ri
          Work(iGP(2)+Ind1-1) = Work(iGP(2)+Ind1-1)+y*Q1*ri
          Work(iGP(3)+Ind1-1) = Work(iGP(3)+Ind1-1)+z*Q1*ri
        end if
        if ((l > nCent-nCha) .and. (ij <= nPol)) then
          ! If ij is
          ! on center with polarizability and l is on
          ! center with charge, then compute the field
          ! from charge l on center ij.
          Q2 = Qsta(l-nCent+nCha)
          Ind1 = (ii-1)*nPol+ij
          Work(iGP(1)+Ind1-1) = Work(iGP(1)+Ind1-1)-x*Q2*ri
          Work(iGP(2)+Ind1-1) = Work(iGP(2)+Ind1-1)-y*Q2*ri
          Work(iGP(3)+Ind1-1) = Work(iGP(3)+Ind1-1)-z*Q2*ri
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
  Work(iFP(1)+i-1) = Work(iGP(1)+i-1)
  Work(iFP(2)+i-1) = Work(iGP(2)+i-1)
  Work(iFP(3)+i-1) = Work(iGP(3)+i-1)
  F = (Work(iFP(1)+i-1)**2+Work(iFP(2)+i-1)**2+Work(iFP(3)+i-1)**2)*Pol(k)
  Epoll = Epoll+F
end do
Epoll = -Epoll*Half
! Image centers: The field at the polarizabilities - reaction field to the point charges added.
do i=iCStart,nCent*nPart
  Q = Qim(i)
  do k=1,nPol
    indSep = (i-iCStart)*nCent*(nPart-iCNum)+k
    indR = k+iCNum*nCent
    indF = k+iCNum*nPol
    do j=nPol*(iCNum+1),IndMa,nPol
      x = CordIm(i,1)-Cordst(indR,1)
      y = CordIm(i,2)-Cordst(indR,2)
      z = CordIm(i,3)-Cordst(indR,3)
      r3 = Work(iDistIm-1+indSep)**3
      Work(iFP(1)+indF-1) = Work(iFP(1)+indF-1)+x*Q*r3
      Work(iFP(2)+indF-1) = Work(iFP(2)+indF-1)+y*Q*r3
      Work(iFP(3)+indF-1) = Work(iFP(3)+indF-1)+z*Q*r3
      indR = indR+nCent
      indSep = indSep+nCent
      indF = indF+nPol
    end do
  end do
end do
! We obtain an initial guess of the induced dipoles on the solvent.
do i=1+nPol*iCnum,IndMa
  k = i-((i-1)/nPol)*nPol
  do l=1,3
    Work(iDt(l)+i-1) = Work(iFP(l)+i-1)*Pol(k)
  end do
end do

return

end subroutine ClasClas
