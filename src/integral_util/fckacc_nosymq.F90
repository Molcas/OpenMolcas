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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine FckAcc_NoSymq(iCmp,jCmp,kCmp,lCmp,Shijij,iShell,nijkl,AOInt,FMat,DMat,nDens,iAO,iAOst,iBas,jBas,kBas,lBas,DoCoul, &
                         DoExch,Dij,Dkl,Dik,Dil,Djk,Djl,ExFac)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmetry adapted Fock matrix.                        *
!                                                                      *
!          This version uses square density and fock matrices          *
!          Modifications by HJW, 28.12.98                              *
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
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Gateway_Info, only: ThrInt
use Constants, only: Zero, One, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp, jCmp, kCmp, lCmp, iShell(4), nijkl, nDens, iAO(4), iAOst(4), iBas, jBas, kBas, lBas
logical(kind=iwp), intent(in) :: Shijij, DoCoul, DoExch
real(kind=wp), intent(in) :: AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), DMat(nDens), Dij, Dkl, Dik, Dil, Djk, Djl, ExFac
real(kind=wp), intent(inout) :: FMat(nDens)
integer(kind=iwp), parameter :: nCBMax = 200
integer(kind=iwp) :: i, i1, i2, i3, i4, ic, icb, ii, ij, ij_, ijkl, ik, il, Indx(3,nCBMax,4), iSO, j, jcb, jj, jk, jl, jSO, k, &
                     kcb, kk, kl, kl_, kSO, l, lcb, ll, lSO, nBasx(4), ncb_Max, nCmpx(4), nij, ntg
real(kind=wp) :: AOijkl, DMax, Fac, Fac_C, Fac_E, Thr
logical(kind=iwp) :: Shij, Shkl
#ifdef _DEBUGPRINT_
real(kind=wp), external :: DDot_
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
if ((.not. DoExch) .and. (.not. DoCoul)) return
#ifdef _DEBUGPRINT_
!write(u6,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
write(u6,*) iCmp,jCmp,kCmp,lCmp
write(u6,*) 'iAO=',iAO
write(u6,*) 'iAOst=',iAOst
write(u6,*) 'iShell=',iShell
write(u6,*) DoCoul,DoExch,Shijij
write(u6,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,[One],0),DDot_(nDens,DMat,1,[One],0)
write(u6,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
#endif

DMax = max(Dij,Dkl,Dik,Dil,Djk,Djl)
if (DMax <= Zero) return
thr = ThrInt/DMax

Shij = (iShell(1) == iShell(2))
Shkl = (iShell(3) == iShell(4))
ntg = nbas(0)
Fac = One
if (Shij) Fac = Fac*Half
if (Shkl) Fac = Fac*Half
if (Shijij) Fac = Fac*Half
Fac_C = Four*Fac
Fac_E = -Fac*Exfac
!if (nijkl /= ibas*jbas*kbas*lbas) call SysHalt('fckacc_nosymq')

nCmpx(1) = iCmp
nCmpx(2) = jCmp
nCmpx(3) = kCmp
nCmpx(4) = lCmp
nBasx(1) = iBas
nBasx(2) = jBas
nBasx(3) = kBas
nBasx(4) = lBas
nCB_Max = max(iCmp*iBas,jCmp*jBas,kCmp*kBas,lCmp*lBas)
if (nCB_Max > nCBMax) then
  call WarningMessage(2,'FckAcc_NoSymq: nCB_Max > nCBMax')
  write(u6,*) 'nCB_Max=',nCB_Max
  call Abend()
end if
do ii=1,4
  iCB = 0
  do iC=1,nCmpx(ii)
    iSO = iAOtSO(iAO(ii)+iC,0)+iAOSt(ii)
    do i=1,nBasx(ii)
      iCB = iCB+1
      Indx(1,iCB,ii) = iC
      Indx(2,iCB,ii) = i
      Indx(3,iCB,ii) = iSO
      iSO = iSO+1
    end do
  end do
end do

nij = iBas*jBas

if (.not. DoExch) then
  do lCB=1,lCmp*lBas
    i4 = Indx(1,lCB,4)
    l = Indx(2,lCB,4)
    lSO = Indx(3,lCB,4)
    ll = (lSO-1)*ntg

    do kCB=1,kCmp*kBas
      i3 = Indx(1,kCB,3)
      k = Indx(2,kCB,3)
      kSO = Indx(3,kCB,3)
      kk = (kSO-1)*ntg
      kl = ll+kSO

      kl_ = (l-1)*kBas+k

      do jCB=1,jCmp*jBas
        i2 = Indx(1,jCB,2)
        j = Indx(2,jCB,2)
        jSO = Indx(3,jCB,2)
        jj = (jSO-1)*ntg

        do iCB=1,iCmp*iBas
          i1 = Indx(1,iCB,1)
          i = Indx(2,iCB,1)
          iSO = Indx(3,iCB,1)
          ij = jj+iSO

          ij_ = (j-1)*iBas+i
          ijkl = (kl_-1)*nij+ij_
          AOijkl = Fac_C*AOInt(ijkl,i1,i2,i3,i4)

          if (abs(AOijkl) > Thr) then

            FMat(ij) = FMat(ij)+AOijkl*DMat(kl)
            FMat(kl) = FMat(kl)+AOijkl*DMat(ij)

          end if
        end do
      end do
    end do
  end do
else if (.not. DoCoul) then
  do lCB=1,lCmp*lBas
    i4 = Indx(1,lCB,4)
    l = Indx(2,lCB,4)
    lSO = Indx(3,lCB,4)
    ll = (lSO-1)*ntg

    do kCB=1,kCmp*kBas
      i3 = Indx(1,kCB,3)
      k = Indx(2,kCB,3)
      kSO = Indx(3,kCB,3)
      kk = (kSO-1)*ntg

      kl_ = (l-1)*kBas+k

      do jCB=1,jCmp*jBas
        i2 = Indx(1,jCB,2)
        j = Indx(2,jCB,2)
        jSO = Indx(3,jCB,2)
        jk = kk+jSO
        jl = ll+jSO

        do iCB=1,iCmp*iBas
          i1 = Indx(1,iCB,1)
          i = Indx(2,iCB,1)
          iSO = Indx(3,iCB,1)
          ik = kk+iSO
          il = ll+iSO

          ij_ = (j-1)*iBas+i
          ijkl = (kl_-1)*nij+ij_
          AOijkl = Fac_E*AOInt(ijkl,i1,i2,i3,i4)

          if (abs(AOijkl) > Thr) then

            FMat(ik) = FMat(ik)+AOijkl*DMat(jl)
            FMat(il) = FMat(il)+AOijkl*DMat(jk)
            FMat(jl) = FMat(jl)+AOijkl*DMat(ik)
            FMat(jk) = FMat(jk)+AOijkl*DMat(il)

          end if
        end do
      end do
    end do
  end do
else
  do lCB=1,lCmp*lBas
    i4 = Indx(1,lCB,4)
    l = Indx(2,lCB,4)
    lSO = Indx(3,lCB,4)
    ll = (lSO-1)*ntg

    do kCB=1,kCmp*kBas
      i3 = Indx(1,kCB,3)
      k = Indx(2,kCB,3)
      kSO = Indx(3,kCB,3)
      kk = (kSO-1)*ntg
      kl = ll+kSO

      kl_ = (l-1)*kBas+k

      do jCB=1,jCmp*jBas
        i2 = Indx(1,jCB,2)
        j = Indx(2,jCB,2)
        jSO = Indx(3,jCB,2)
        jj = (jSO-1)*ntg
        jk = kk+jSO
        jl = ll+jSO

        do iCB=1,iCmp*iBas
          i1 = Indx(1,iCB,1)
          i = Indx(2,iCB,1)
          iSO = Indx(3,iCB,1)
          ik = kk+iSO
          il = ll+iSO
          ij = jj+iSO

          ij_ = (j-1)*iBas+i
          ijkl = (kl_-1)*nij+ij_
          AOijkl = AOInt(ijkl,i1,i2,i3,i4)
          if (abs(AOijkl) > Thr) then

            FMat(ij) = FMat(ij)+Fac_C*AOijkl*DMat(kl)
            FMat(kl) = FMat(kl)+Fac_C*AOijkl*DMat(ij)
            FMat(ik) = FMat(ik)+Fac_E*AOijkl*DMat(jl)
            FMat(il) = FMat(il)+Fac_E*AOijkl*DMat(jk)
            FMat(jl) = FMat(jl)+Fac_E*AOijkl*DMat(ik)
            FMat(jk) = FMat(jk)+Fac_E*AOijkl*DMat(il)

          end if
        end do
      end do
    end do
  end do
end if

!write(u6,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,One,0),DDot_(nDens,DMat,1,One,0)

return

end subroutine FckAcc_NoSymq
