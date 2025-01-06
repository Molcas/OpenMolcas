!*********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993,1998, Roland Lindh                                *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine FckAcc(iAng,iCmp_,Shijij,iShll,iShell,kOp,nijkl,AOInt,TwoHam,nDens,Scrt,nScrt,iAO,iAOst,iBas,jBas,kBas,lBas,Dij,ij1, &
                  ij2,ij3,ij4,Dkl,kl1,kl2,kl3,kl4,Dik,ik1,ik2,ik3,ik4,Dil,il1,il2,il3,il4,Djk,jk1,jk2,jk3,jk4,Djl,jl1,jl2,jl3,jl4, &
                  DoCoul,DoExch,ExFac)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmetry adapted Fock matrix.                        *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          In addition to this complication we have that the order of  *
!          indices in the integrals are not ordered canonically but    *
!          rather in an order such that the contraction step will be   *
!          optimal. Hence, special care has to be taken when tracing   *
!          the density with the integrals so that both entities have   *
!          the same order.                                             *
!                                                                      *
!          The Fock matrix is computed in lower triangular form.       *
!                                                                      *
!          The density matrix is not folded if the shell indices and   *
!          the angular indices are identical.                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. February '93                            *
!                                                                      *
!     Modified July '98 in Tokyo by R. Lindh                           *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Symmetry_Info, only: iChBas, iOper, nIrrep, Prmt
use Gateway_Info, only: CutInt, ThrInt
use k2_arrays, only: FT
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAng(4), iCmp_(4), iShll(4), iShell(4), kOp(4), nijkl, nDens, nScrt, iAO(4), iAOst(4), iBas, &
                                 jBas, kBas, lBas, ij1, ij2, ij3, ij4, kl1, kl2, kl3, kl4, ik1, ik2, ik3, ik4, il1, il2, il3, il4, &
                                 jk1, jk2, jk3, jk4, jl1, jl2, jl3, jl4
logical(kind=iwp), intent(in) :: Shijij, DoCoul, DoExch
real(kind=wp), intent(in) :: AOInt(nijkl,iCmp_(1),iCmp_(2),iCmp_(3),iCmp_(4)), ExFac
real(kind=wp), intent(inout) :: TwoHam(nDens)
real(kind=wp), target, intent(inout) :: Scrt(nScrt)
real(kind=wp), target, intent(in) :: Dij(ij1*ij2+1,ij3,ij4), Dkl(kl1*kl2+1,kl3,kl4), Dik(ik1*ik2+1,ik3,ik4), &
                                     Dil(il1*il2+1,il3,il4), Djk(jk1*jk2+1,jk3,jk4), Djl(jl1*jl2+1,jl3,jl4)
integer(kind=iwp) :: i1, i12, i2, i3, i34, i4, iChBs, iCmp, iCmpa(4), ii, iIrrep, ijkl, iOpt, ip, ipF, iSym(4), ix, j, jChBs, &
                     jCmp, jCmpMx, jj, kChBs, kCmp, kk, kOp2(4), lChBs, lCmp, lCmpMx, ll, mijkl, nF
real(kind=wp) :: D_ij, D_ik, D_il, D_jk, D_jl, D_kl, Fac_ij, Fac_ik, Fac_il, Fac_jk, Fac_jl, Fac_kl, Fact, pEa, pFctr, pRb, pTc, &
                 pTSd, Vij, Vijkl, Vik, Vil, Vjk, Vjl, Vkl
logical(kind=iwp) :: iQij, iQik, iQil, iQjk, iQjl, iQkl, iShij, iShik, iShil, iShjk, iShjl, iShkl, lFij, lFik, lFil, lFjk, lFjl, &
                     lFkl, Qijij
real(kind=wp), pointer :: Fij(:,:,:), Fik(:,:,:), Fil(:,:,:), Fjk(:,:,:), Fjl(:,:,:), Fkl(:,:,:), pDij(:), pDik(:), pDil(:), &
                          pDjk(:), pDjl(:), pDkl(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Since I sometimes use Scrt as an anchor to reach into the
! density matrix a check if I'm out of bounds does not work!

if ((.not. DoCoul) .and. (.not. DoExch)) return

iCmp = iCmp_(1)
jCmp = iCmp_(2)
kCmp = iCmp_(3)
lCmp = iCmp_(4)
#ifdef _DEBUGPRINT_
call RecPrt('FckAcc:AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
write(u6,*) 'Dij=',sum(Dij(1:ij1*ij2,:,:))
write(u6,*) 'Dkl=',sum(Dkl(1:kl1*kl2,:,:))
write(u6,*) 'Dik=',sum(Dik(1:ik1*ik2,:,:))
write(u6,*) 'Dil=',sum(Dil(1:il1*il2,:,:))
write(u6,*) 'Djk=',sum(Djk(1:jk1*jk2,:,:))
write(u6,*) 'Djl=',sum(Djl(1:jl1*jl2,:,:))
call RecPrt('AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
#endif

if (iBas*jBas*kBas*lBas > nScrt) then
  call WarningMessage(2,'FckAcc: nScrt too small!')
  call Abend()
end if
ii = nTri3_Elem(iAng(1))
jj = nTri3_Elem(iAng(2))
kk = nTri3_Elem(iAng(3))
ll = nTri3_Elem(iAng(4))

kOp2(1) = iOper(kOp(1))
kOp2(2) = iOper(kOp(2))
kOp2(3) = iOper(kOp(3))
kOp2(4) = iOper(kOp(4))
iCmpa(1) = iCmp
iCmpa(2) = jCmp
iCmpa(3) = kCmp
iCmpa(4) = lCmp
lFij = .false.
lFkl = .false.
lFik = .false.
lFjl = .false.
lFil = .false.
lFjk = .false.

ipF = 1
nF = iBas*jBas*iCmpa(1)*iCmpa(2)
Fij(1:iBas*jBas,1:iCmpa(1),1:iCmpa(2)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
nF = kBas*lBas*iCmpa(3)*iCmpa(4)
Fkl(1:kBas*lBas,1:iCmpa(3),1:iCmpa(4)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
nF = iBas*kBas*iCmpa(1)*iCmpa(3)
Fik(1:iBas*kBas,1:iCmpa(1),1:iCmpa(3)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
nF = jBas*lBas*iCmpa(2)*iCmpa(4)
Fjl(1:jBas*lBas,1:iCmpa(2),1:iCmpa(4)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
nF = iBas*lBas*iCmpa(1)*iCmpa(4)
Fil(1:iBas*lBas,1:iCmpa(1),1:iCmpa(4)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
nF = jBas*kBas*iCmpa(2)*iCmpa(3)
Fjk(1:jBas*kBas,1:iCmpa(2),1:iCmpa(3)) => FT(ipF:ipF+nF-1)

ipF = ipF+nF
FT(1:ipF-1) = Zero

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

iShij = (iShell(1) == iShell(2))
iShkl = (iShell(3) == iShell(4))
iShik = (iShell(1) == iShell(3))
iShil = (iShell(1) == iShell(4))
iShjk = (iShell(2) == iShell(3))
iShjl = (iShell(2) == iShell(4))
mijkl = iBas*jBas*kBas*lBas
do i1=1,iCmp
  ix = 0
  do j=0,nIrrep-1
    if (iAOtSO(iAO(1)+i1,j) > 0) ix = ibset(ix,j)
  end do
  iSym(1) = ix
  jCmpMx = jCmp
  if (iShij) jCmpMx = i1
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pEa = Prmt(iOper(kOp(1)),iChBs)
  do i2=1,jCmpMx
    ix = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) ix = ibset(ix,j)
    end do
    iSym(2) = ix
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pRb = Prmt(iOper(kOp(2)),jChBs)
    if (iShell(2) > iShell(1)) then
      i12 = jCmp*(i1-1)+i2
    else
      i12 = iCmp*(i2-1)+i1
    end if
    do i3=1,kCmp
      ix = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(3)+i3,j) > 0) ix = ibset(ix,j)
      end do
      iSym(3) = ix
      lCmpMx = lCmp
      if (iShkl) lCmpMx = i3
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pTc = Prmt(iOper(kOp(3)),kChBs)
      do i4=1,lCmpMx
        ix = 0
        do j=0,nIrrep-1
          if (iAOtSO(iAO(4)+i4,j) > 0) ix = ibset(ix,j)
        end do
        iSym(4) = ix
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        pTSd = Prmt(iOper(kOp(4)),lChBs)
        if (iShell(4) > iShell(3)) then
          i34 = lCmp*(i3-1)+i4
        else
          i34 = kCmp*(i4-1)+i3
        end if
        if (Shijij .and. (i34 > i12)) cycle
        vijkl = Zero
        do ijkl=1,mijkl
          vijkl = max(vijkl,abs(AOInt(ijkl,i1,i2,i3,i4)))
        end do
        !vijkl = DNrm2_(iBas*jBas*kBas*lBas,AOInt(1,i1,i2,i3,i4),1)
        if (vijkl < CutInt) cycle

        Qijij = (Shijij .and. (i12 == i34))
        iQij = (iShij .and. (i1 == i2))
        iQkl = (iShkl .and. (i3 == i4))
        iQik = (iShik .and. (i1 == i3))
        iQil = (iShil .and. (i1 == i4))
        iQjk = (iShjk .and. (i2 == i3))
        iQjl = (iShjl .and. (i2 == i4))
        pFctr = pEa*pRb*pTc*pTSd
        !***************************************************************

        Fac_ij = pFctr
        Fac_kl = pFctr
        Fac_ik = -Quart*pFctr
        Fac_jl = -Quart*pFctr
        Fac_il = -Quart*pFctr
        Fac_jk = -Quart*pFctr
        if (iQij) Fac_ij = Fac_ij*Half
        if (iQkl) Fac_kl = Fac_kl*Half
        if (iQjl .and. (.not. iQik)) Fac_ik = Fac_ik*Two
        if (iQik .and. (.not. iQjl)) Fac_jl = Fac_jl*Two
        if (iQjk .and. (.not. iQil)) Fac_il = Fac_il*Two
        if (iQil .and. (.not. iQjk)) Fac_jk = Fac_jk*Two

        D_ij = Two
        if (iQij) D_ij = One
        D_kl = Two
        if (iQkl) D_kl = One
        D_ik = Two
        if (iQik) D_ik = One
        D_jl = Two
        if (iQjl) D_jl = One
        D_il = Two
        if (iQil) D_il = One
        D_jk = Two
        if (iQjk) D_jk = One

        Fac_ij = Fac_ij*D_kl
        Fac_kl = Fac_kl*D_ij
        Fac_jl = Fac_jl*D_ik
        Fac_ik = Fac_ik*D_jl
        Fac_jk = Fac_jk*D_il
        Fac_il = Fac_il*D_jk

        if (Qijij) then
          Fac_kl = Zero
          Fac_jk = Zero
        end if
        if (iQij .and. iQkl) then
          Fac_jl = Zero
          Fac_il = Zero
          Fac_jk = Zero
        end if
        if (iQij .or. iQkl) then
          Fac_il = Zero
          Fac_jk = Zero
        end if
        !***************************************************************

        iOpt = 0
        ip = 1
        if (DoCoul .and. (iand(iSym(1),iSym(2)) /= 0) .and. (iand(iSym(3),iSym(4)) /= 0)) then
          iOpt = iOpt+1

          if (iShell(3) < iShell(4)) then
            vkl = Dkl(kl1*kl2+1,i4,i3)
            pDkl(1:kl1*kl2) => Scrt(ip:ip+kl1*kl2-1)
            call DGetMO(Dkl(1,i4,i3),kl1,kl1,kl2,pDkl,kl2)
            ip = ip+kl1*kl2
          else
            vkl = Dkl(kl1*kl2+1,i3,i4)
            pDkl(1:kl1*kl2) => Dkl(1:kl1*kl2,i3,i4)
          end if
          if (iShell(1) < iShell(2)) then
            vij = Dij(ij1*ij2+1,i2,i1)
            pDij(1:ij1*ij2) => Scrt(ip:ip+ij1*ij2-1)
            call DGeTMO(Dij(1,i2,i1),ij1,ij1,ij2,pDij,ij2)
            ip = ip+ij1*ij2
          else
            vij = Dij(ij1*ij2+1,i1,i2)
            pDij(1:ij1*ij2) => Dij(1:ij1*ij2,i1,i2)
          end if
          if ((vkl*vijkl*abs(Fac_ij) < ThrInt) .and. (vij*vijkl*abs(Fac_kl) < ThrInt)) then
            iOpt = iOpt-1
          else
            lFij = .true.
            lFkl = .true.
          end if
        end if
        if (DoExch .and. (iand(iSym(1),iSym(3)) /= 0) .and. (iand(iSym(2),iSym(4)) /= 0)) then
          iOpt = iOpt+2

          if (iShell(2) < iShell(4)) then
            vjl = Djl(jl1*jl2+1,i4,i2)
            pDjl(1:jl1*jl2) => Scrt(ip:ip+jl1*jl2-1)
            call DGeTMO(Djl(1,i4,i2),jl1,jl1,jl2,pDjl,jl2)
            ip = ip+jl1*jl2
          else
            vjl = Djl(jl1*jl2+1,i2,i4)
            pDjl(1:jl1*jl2) => Djl(1:jl1*jl2,i2,i4)
          end if
          if (iShell(1) < iShell(3)) then
            vik = Dik(ik1*ik2+1,i3,i1)
            pDik(1:ik1*ik2) => Scrt(ip:ip+ik1*ik2-1)
            call DGeTMO(Dik(1,i3,i1),ik1,ik1,ik2,pDik,ik2)
            ip = ip+ik1*ik2
          else
            vik = Dik(ik1*ik2+1,i1,i3)
            pDik(1:ik1*ik2) => Dik(1:ik1*ik2,i1,i3)
          end if
          if ((vik*vijkl*abs(Fac_jl) < ThrInt) .and. (vjl*vijkl*abs(Fac_ik) < ThrInt)) then
            iOpt = iOpt-2
          else
            lFik = .true.
            lFjl = .true.
          end if
        end if
        if (DoExch .and. (iand(iSym(1),iSym(4)) /= 0) .and. (iand(iSym(2),iSym(3)) /= 0)) then
          iOpt = iOpt+4

          if (iShell(2) < iShell(3)) then
            vjk = Djk(jk1*jk2+1,i3,i2)
            pDjk(1:jk1*jk2) => Scrt(ip:ip+jk1*jk2-1)
            call DGeTMO(Djk(1,i3,i2),jk1,jk1,jk2,pDjk,jk2)
            ip = ip+jk1*jk2
          else
            vjk = Djk(jk1*jk2+1,i2,i3)
            pDjk(1:jk1*jk2) => Djk(1:jk1*jk2,i2,i3)
          end if
          if (iShell(1) < iShell(4)) then
            vil = Dil(il1*il2+1,i4,i1)
            pDil(1:il1*il2) => Scrt(ip:ip+il1*il2-1)
            call DGeTMO(Dil(1,i4,i1),il1,il1,il2,pDil,il2)
            ip = ip+il1*il2
          else
            vil = Dil(il1*il2+1,i1,i4)
            pDil(1:il1*il2) => Dil(1:il1*il2,i1,i4)
          end if
          if ((vil*vijkl*abs(Fac_jk) < ThrInt) .and. (vjk*vijkl*abs(Fac_il) < ThrInt)) then
            iOpt = iOpt-4
          else
            lFil = .true.
            lFjk = .true.
          end if
        end if
        if (ip-1 > nScrt) then
          call WarningMessage(2,'FckAcc: nScrt too small!')
          call Abend()
        end if

        !*************************************************************
        select case (iOpt)
          case (1)
            call Fck1(AOInt(:,i1,i2,i3,i4),pDij,Fij(:,i1,i2),Fac_ij,pDkl,Fkl(:,i3,i4),Fac_kl)
          case (2)
            call Fck2(AOInt(:,i1,i2,i3,i4),pDik,Fik(:,i1,i3),Fac_ik,pDjl,Fjl(:,i2,i4),Fac_jl)
          case (3)
            call Fck3(AOInt(:,i1,i2,i3,i4),pDij,Fij(:,i1,i2),Fac_ij,pDkl,Fkl(:,i3,i4),Fac_kl,pDik,Fik(:,i1,i3),Fac_ik,pDjl, &
                      Fjl(:,i2,i4),Fac_jl)
          case (4)
            call Fck4(AOInt(:,i1,i2,i3,i4),pDil,Fil(:,i1,i4),Fac_il,pDjk,Fjk(:,i2,i3),Fac_jk)
          case (5)
            call Fck5(AOInt(:,i1,i2,i3,i4),pDij,Fij(:,i1,i2),Fac_ij,pDkl,Fkl(:,i3,i4),Fac_kl,pDil,Fil(:,i1,i4),Fac_il,pDjk, &
                      Fjk(:,i2,i3),Fac_jk)
          case (6)
            call Fck6(AOInt(:,i1,i2,i3,i4),pDik,Fik(:,i1,i3),Fac_ik,pDjl,Fjl(:,i2,i4),Fac_jl,pDil,Fil(:,i1,i4),Fac_il,pDjk, &
                      Fjk(:,i2,i3),Fac_jk)
          case (7)
            call Fck7(AOInt(:,i1,i2,i3,i4),pDij,Fij(:,i1,i2),Fac_ij,pDkl,Fkl(:,i3,i4),Fac_kl,pDik,Fik(:,i1,i3),Fac_ik,pDjl, &
                      Fjl(:,i2,i4),Fac_jl,pDil,Fil(:,i1,i4),Fac_il,pDjk,Fjk(:,i2,i3),Fac_jk)
          case default
        end select
        nullify(pDij,pDkl,pDik,pDil,pDjk,pDjl)
        !***************************************************************

      end do
    end do
  end do
end do

iIrrep = 0
Fact = One
if (lFij) &
  call FckDst(TwoHam,nDens,Fij,iBas,jBas,iCmpa(1),iCmpa(2),kOp2(1),kOp2(2),iIrrep,iShij,iAO(1),iAO(2),iAOst(1),iAOst(2),Fact)
if (lFkl) &
  call FckDst(TwoHam,nDens,Fkl,kBas,lBas,iCmpa(3),iCmpa(4),kOp2(3),kOp2(4),iIrrep,iShkl,iAO(3),iAO(4),iAOst(3),iAOst(4),Fact)
if (lFik) &
  call FckDst(TwoHam,nDens,Fik,iBas,kBas,iCmpa(1),iCmpa(3),kOp2(1),kOp2(3),iIrrep,iShik,iAO(1),iAO(3),iAOst(1),iAOst(3),Fact)
if (lFjl) &
  call FckDst(TwoHam,nDens,Fjl,jBas,lBas,iCmpa(2),iCmpa(4),kOp2(2),kOp2(4),iIrrep,iShjl,iAO(2),iAO(4),iAOst(2),iAOst(4),Fact)
if (lFil) &
  call FckDst(TwoHam,nDens,Fil,iBas,lBas,iCmpa(1),iCmpa(4),kOp2(1),kOp2(4),iIrrep,iShil,iAO(1),iAO(4),iAOst(1),iAOst(4),Fact)
if (lFjk) &
  call FckDst(TwoHam,nDens,Fjk,jBas,kBas,iCmpa(2),iCmpa(3),kOp2(2),kOp2(3),iIrrep,iShjk,iAO(2),iAO(3),iAOst(2),iAOst(3),Fact)

nullify(Fij,Fkl,Fik,Fil,Fjk,Fjl)

contains

subroutine Fck1(AOInt,Dij,Fij,Cij,Dkl,Fkl,Ckl)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dij(iBas,jBas), Cij, Dkl(kBas,lBas), Ckl
  real(kind=wp), intent(inout) :: Fij(iBas,jBas), Fkl(kBas,lBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_kl, F_kl, Vijkl

  do l=1,lBas
    do k=1,kBas
      F_kl = Zero
      D_kl = Dkl(k,l)*Cij
      do j=1,jBas
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)

          Fij(i,j) = Fij(i,j)+D_kl*Vijkl
          F_kl = F_kl+Dij(i,j)*Vijkl

        end do
      end do
      Fkl(k,l) = Fkl(k,l)+Ckl*F_kl
    end do
  end do

  return

end subroutine Fck1

subroutine Fck2(AOInt,Dik,Fik,Cik,Djl,Fjl,Cjl)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dik(iBas,kBas), Cik, Djl(jBas,lBas), Cjl
  real(kind=wp), intent(inout) :: Fik(iBas,kBas), Fjl(jBas,lBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jl, F_jl, Vijkl

  do l=1,lBas
    do k=1,kBas
      do j=1,jBas
        F_jl = Zero
        D_jl = Djl(j,l)*Cik
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)*ExFac

          Fik(i,k) = Fik(i,k)+D_jl*Vijkl
          F_jl = F_jl+Dik(i,k)*Vijkl

        end do
        Fjl(j,l) = Fjl(j,l)+Cjl*F_jl
      end do
    end do
  end do

  return

end subroutine Fck2

subroutine Fck3(AOInt,Dij,Fij,Cij,Dkl,Fkl,Ckl,Dik,Fik,Cik,Djl,Fjl,Cjl)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dij(iBas,jBas), Cij, Dkl(kBas,lBas), Ckl, Dik(iBas,kBas), Cik, &
                               Djl(jBas,lBas), Cjl
  real(kind=wp), intent(inout) :: Fij(iBas,jBas), Fkl(kBas,lBas), Fik(iBas,kBas), Fjl(jBas,lBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jl, D_kl, F_jl, F_kl, Vijkl

  do l=1,lBas
    do k=1,kBas
      F_kl = Zero
      D_kl = Dkl(k,l)*Cij
      do j=1,jBas
        F_jl = Zero
        D_jl = Djl(j,l)*Cik
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)

          Fij(i,j) = Fij(i,j)+D_kl*Vijkl
          F_kl = F_kl+Dij(i,j)*Vijkl

          Fik(i,k) = Fik(i,k)+D_jl*Vijkl*ExFac
          F_jl = F_jl+Dik(i,k)*Vijkl

        end do
        Fjl(j,l) = Fjl(j,l)+Cjl*F_jl*ExFac
      end do
      Fkl(k,l) = Fkl(k,l)+Ckl*F_kl
    end do
  end do

  return

end subroutine Fck3

subroutine Fck4(AOInt,Dil,Fil,Cil,Djk,Fjk,Cjk)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dil(iBas,lBas), Cil, Djk(jBas,kBas), Cjk
  real(kind=wp), intent(inout) :: Fil(iBas,lBas), Fjk(jBas,kBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jk, F_jk, Vijkl

  do l=1,lBas
    do k=1,kBas
      do j=1,jBas
        F_jk = Zero
        D_jk = Djk(j,k)*Cil
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)*ExFac

          Fil(i,l) = Fil(i,l)+D_jk*Vijkl
          F_jk = F_jk+Dil(i,l)*Vijkl

        end do
        Fjk(j,k) = Fjk(j,k)+Cjk*F_jk
      end do
    end do
  end do

  return

end subroutine Fck4

subroutine Fck5(AOInt,Dij,Fij,Cij,Dkl,Fkl,Ckl,Dil,Fil,Cil,Djk,Fjk,Cjk)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dij(iBas,jBas), Cij, Dkl(kBas,lBas), Ckl, Dil(iBas,lBas), Cil, &
                               Djk(jBas,kBas), Cjk
  real(kind=wp), intent(inout) :: Fij(iBas,jBas), Fkl(kBas,lBas), Fil(iBas,lBas), Fjk(jBas,kBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jk, D_kl, F_jk, F_kl, Vijkl

  do l=1,lBas
    do k=1,kBas
      F_kl = Zero
      D_kl = Dkl(k,l)*Cij
      do j=1,jBas
        F_jk = Zero
        D_jk = Djk(j,k)*Cil
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)

          Fij(i,j) = Fij(i,j)+D_kl*Vijkl
          F_kl = F_kl+Dij(i,j)*Vijkl

          Fil(i,l) = Fil(i,l)+D_jk*Vijkl*ExFac
          F_jk = F_jk+Dil(i,l)*Vijkl

        end do
        Fjk(j,k) = Fjk(j,k)+Cjk*F_jk*ExFac
      end do
      Fkl(k,l) = Fkl(k,l)+Ckl*F_kl
    end do
  end do

  return

end subroutine Fck5

subroutine Fck6(AOInt,Dik,Fik,Cik,Djl,Fjl,Cjl,Dil,Fil,Cil,Djk,Fjk,Cjk)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dik(iBas,kBas), Cik, Djl(jBas,lBas), Cjl, Dil(iBas,lBas), Cil, &
                               Djk(jBas,kBas), Cjk
  real(kind=wp), intent(inout) :: Fik(iBas,kBas), Fjl(jBas,lBas), Fil(iBas,lBas), Fjk(jBas,kBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jk, D_jl, F_jk, F_jl, Vijkl

  do l=1,lBas
    do k=1,kBas
      do j=1,jBas
        F_jl = Zero
        D_jl = Djl(j,l)*Cik
        F_jk = Zero
        D_jk = Djk(j,k)*Cil
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)

          Fik(i,k) = Fik(i,k)+D_jl*Vijkl*ExFac
          F_jl = F_jl+Dik(i,k)*Vijkl

          Fil(i,l) = Fil(i,l)+D_jk*Vijkl*ExFac
          F_jk = F_jk+Dil(i,l)*Vijkl

        end do
        Fjl(j,l) = Fjl(j,l)+Cjl*F_jl*ExFac
        Fjk(j,k) = Fjk(j,k)+Cjk*F_jk*ExFac
      end do
    end do
  end do

  return

end subroutine Fck6

subroutine Fck7(AOInt,Dij,Fij,Cij,Dkl,Fkl,Ckl,Dik,Fik,Cik,Djl,Fjl,Cjl,Dil,Fil,Cil,Djk,Fjk,Cjk)

  real(kind=wp), intent(in) :: AOInt(iBas,jBas,kBas,lBas), Dij(iBas,jBas), Cij, Dkl(kBas,lBas), Ckl, Dik(iBas,kBas), Cik, &
                               Djl(jBas,lBas), Cjl, Dil(iBas,lBas), Cil, Djk(jBas,kBas), Cjk
  real(kind=wp), intent(inout) :: Fij(iBas,jBas), Fkl(kBas,lBas), Fik(iBas,kBas), Fjl(jBas,lBas), Fil(iBas,lBas), Fjk(jBas,kBas)
  integer(kind=iwp) :: i, j, k, l
  real(kind=wp) :: D_jk, D_jl, D_kl, F_jk, F_jl, F_kl, Vijkl

  do l=1,lBas
    do k=1,kBas
      F_kl = Zero
      D_kl = Dkl(k,l)*Cij
      do j=1,jBas
        F_jl = Zero
        D_jl = Djl(j,l)*Cik
        F_jk = Zero
        D_jk = Djk(j,k)*Cil
        do i=1,iBas
          Vijkl = AOInt(i,j,k,l)

          Fij(i,j) = Fij(i,j)+D_kl*Vijkl
          F_kl = F_kl+Dij(i,j)*Vijkl

          Fik(i,k) = Fik(i,k)+D_jl*Vijkl*ExFac
          F_jl = F_jl+Dik(i,k)*Vijkl

          Fil(i,l) = Fil(i,l)+D_jk*Vijkl*ExFac
          F_jk = F_jk+Dil(i,l)*Vijkl

        end do
        Fjl(j,l) = Fjl(j,l)+Cjl*F_jl*ExFac
        Fjk(j,k) = Fjk(j,k)+Cjk*F_jk*ExFac
      end do
      Fkl(k,l) = Fkl(k,l)+Ckl*F_kl
    end do
  end do

  return

end subroutine Fck7

end subroutine FckAcc
