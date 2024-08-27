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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SymAdp(iAng,iCmp,jCmp,kCmp,lCmp,Shijij,iShll,iShell,iAO,kOp,ijkl,Aux,nAux,AOInt,SOInt,nSOInt,Done)
!***********************************************************************
!                                                                      *
!  Object: to transform the integrals in AO basis to symmetry adapted  *
!          orbitals , SO. This is done by accumulating the AO inte-    *
!          grals onto the SO integrals.                                *
!          Observe that one of the operator is the Unit operation      *
!          performed on center A. However, since we scramble the order *
!          we do not really know which center is which. However, the   *
!          Unit operator will always give a factor of one. Hence, this *
!          is a convenient way to do the symmetry transformation with- *
!          out having to know the order.                               *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          This code is never executed in the no symmetry case!!!      *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr

implicit none
integer, intent(In) :: ijkl, iCmp, jCmp, kCmp, lCmp, nSOInt, nAux
real*8, intent(In) :: AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
real*8, intent(InOut) :: SOInt(ijkl,nSOInt), Aux(nAux)
logical, intent(In) :: Shijij
integer, intent(In) :: iAng(4), iShell(4), iShll(4), kOp(4), iAO(4)
logical, intent(Out) :: Done
logical Shij, Shkl, Qij, Qkl, Qijij
integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
integer k12, k34, ii, jj, kk, ll, i1, j, ix, jCmpMx, i2, jChBs, iChBs, MemSO2, i12, i3, lCmpMx, kChBs, i4, iAux, j1, j2Max, j2, &
        j12, j3, j4, lChBs, i34
real*8 pEa, pRb, pTc, pTSd, Xb, Xg, Xa
! Statement Function
integer ixyz, iOff
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

Done = .false.
k12 = 0
k34 = 0
ii = iOff(iAng(1))
jj = iOff(iAng(2))
kk = iOff(iAng(3))
ll = iOff(iAng(4))
MemSO2 = 1
#ifdef _DEBUGPRINT_
call RecPrt(' In SymAdp: AOInt ',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

Shij = iShell(1) == iShell(2)
Shkl = iShell(3) == iShell(4)
do i1=1,iCmp
  do j=0,nIrrep-1
    ix = 0
    if (iAOtSO(iAO(1)+i1,j) > 0) ix = 2**j
    iSym(j) = ix
  end do
  jCmpMx = jCmp
  if (Shij) jCmpMx = i1
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pEa = Prmt(iOper(kOp(1)),iChBs)
  do i2=1,jCmpMx
    do j=0,nIrrep-1
      ix = 0
      if (iAOtSO(iAO(2)+i2,j) > 0) ix = 2**j
      jSym(j) = ix
    end do
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pRb = Prmt(iOper(kOp(2)),jChBs)*pEa
    Qij = i1 == i2
    if (iShell(2) > iShell(1)) then
      i12 = jCmp*(i1-1)+i2
    else
      i12 = iCmp*(i2-1)+i1
    end if
    do i3=1,kCmp
      do j=0,nIrrep-1
        ix = 0
        if (iAOtSO(iAO(3)+i3,j) > 0) ix = 2**j
        kSym(j) = ix
      end do
      lCmpMx = lCmp
      if (Shkl) lCmpMx = i3
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pTc = Prmt(iOper(kOp(3)),kChBs)*pRb
      do i4=1,lCmpMx
        do j=0,nIrrep-1
          ix = 0
          if (iAOtSO(iAO(4)+i4,j) > 0) ix = 2**j
          lSym(j) = ix
        end do
        Qkl = i3 == i4
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        pTSd = Prmt(iOper(kOp(4)),lChBs)*pTc
        if (iShell(4) > iShell(3)) then
          i34 = lCmp*(i3-1)+i4
        else
          i34 = kCmp*(i4-1)+i3
        end if
        if (Shijij .and. (i34 > i12)) Go To 400
        Qijij = Shijij .and. (i12 == i34)

        ! Loop over irreps which are spanned by the basis function.
        ! Again, the loop structure is restricted to ensure unique
        ! integrals.

        iAux = 0
        do j1=0,nIrrep-1
          if (iSym(j1) == 0) Go To 110
          Xa = dble(iChTbl(j1,kOp(1)))*pTSd
          j2Max = nIrrep-1
          if (Shij .and. Qij) j2Max = j1
          do j2=0,j2Max
            if (jSym(j2) == 0) Go To 210
            Xb = dble(iChTbl(j2,kOp(2)))*Xa
            j12 = ieor(j1,j2)
            if (Qijij) then
              if (Shij .and. Qij) then
                k12 = j1*(j1+1)/2+j2+1
              else if (Shij) then
                k12 = nIrrep*j1+j2+1
              else if (iShell(1) > iShell(2)) then
                k12 = nIrrep*j1+j2+1
              else
                k12 = nIrrep*j2+j1+1
              end if
            end if
            do j3=0,nIrrep-1
              if (kSym(j3) == 0) Go To 310
              j4 = ieor(j12,j3)
              if (lSym(j4) == 0) Go To 310
              if (Shkl .and. Qkl .and. (j4 > j3)) Go To 310
              if (Qijij) then
                if (Shkl .and. Qkl) then
                  k34 = j3*(j3+1)/2+j4+1
                else if (Shkl) then
                  k34 = nIrrep*j3+j4+1
                else if (iShell(3) > iShell(4)) then
                  k34 = nIrrep*j3+j4+1
                else
                  k34 = nIrrep*j4+j3+1
                end if
                if (Qijij .and. (k34 > k12)) Go To 310
              end if
              Xg = dble(iChTbl(j3,kOp(3)))*Xb
              iAux = iAux+1
              Aux(iAux) = dble(iChTbl(j4,kOp(4)))*Xg

310           continue
            end do
210         continue
          end do
110       continue
        end do

#       ifdef _DEBUGPRINT_
        call RecPrt(' Aux',' ',Aux,iAux,1)
#       endif
        if (iAux /= 0) then
          if (iAux /= 1) then
            call DNaXpY(iAux,ijkl,Aux,1,AOInt(1,i1,i2,i3,i4),1,0,SOInt(1,MemSO2),1,ijkl)
          else
            call DaXpY_(ijkl,Aux(1),AOInt(1,i1,i2,i3,i4),1,SOInt(1,MemSO2),1)
          end if
          MemSO2 = MemSO2+iAux
        end if

400     continue
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' On exit from SymAdp: SOInt ',' ',SOInt,ijkl,nSOInt)
#endif

end subroutine SymAdp
