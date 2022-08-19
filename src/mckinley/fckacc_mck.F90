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
! Copyright (C) 1993,1998, Roland Lindh                                *
!***********************************************************************

#define _USE_OLD_CODE_
#ifdef _USE_OLD_CODE_

subroutine FckAcc_mck(iAng,iCmp,jCmp,kCmp,lCmp,Shijij,iShll,iShell,kOp,nijkl,AOInt,TwoHam,nDens,Scrt,nScrt,iAO,iAOst,iBas,jBas, &
                      kBas,lBas,Dij,ij1,ij2,ij3,ij4,Dkl,kl1,kl2,kl3,kl4,Dik,ik1,ik2,ik3,ik4,Dil,il1,il2,il3,il4,Djk,jk1,jk2,jk3, &
                      jk4,Djl,jl1,jl2,jl3,jl4,FT,nFT,tfact,iCar,iCent,pert,indgrd,ipdisp)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmatry adapted Fock matrix.                        *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          In addition to this complication we have that the order of  *
!          indicies in the integrals are not ordered canonically but   *
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

use Basis_Info
use Symmetry_Info, only: nIrrep, iChTbl, iOper, iChBas
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Gateway_Info, only: CutInt

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "disp2.fh"
real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), TwoHam(nDens), Scrt(nScrt), FT(nFT), Dij(ij1*ij2+1,ij3,ij4), Dkl(kl1*kl2+1,kl3,kl4), &
       Dik(ik1*ik2+1,ik3,ik4), Dil(il1*il2+1,il3,il4), Djk(jk1*jk2+1,jk3,jk4), Djl(jl1*jl2+1,jl3,jl4)
!logical Qij, Qkl
logical Shij, Shkl, Shijij, Qijij, iShij, iShkl, iQij, iQkl, iQik, iShik, iQil, iShil, iQjk, iShjk, iQjl, iShjl, lFij, lFkl, lFik, &
        lFjl, lFil, lFjk
integer iAng(4), iShell(4), iShll(4), kOp(4), kOp2(4), iAO(4), iAOst(4), iCmpa(4)
logical Pert(0:nIrrep-1)
integer indgrd(3,4,0:nirrep-1), ipdisp(*)
! Local Arrays
integer iSym(4,0:7), iTwoj(0:7)
real*8 Prmt(0:7)
data iTwoj/1,2,4,8,16,32,64,128/
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Function
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
xPrmt(i,j) = Prmt(iand(i,j))

!iprint = 0

!write(6,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
!if (iPrint >= 49) then
!  if ((DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1) > Zero) .or. (DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0) > Zero)) then
!    call RecPrt('Dik','(5G20.10)',Dik,ik1*ik2+1,ik3*ik4)
!    write(6,'(A,2G20.10)') &
!      ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
!  end if
!end if
!if (iPrint >= 99) then
!  call RecPrt('FckAcc:AOInt','(5G20.10)',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
!  write(6,'(A,G20.10)') 'Dij=',XDot(Dij,ij1,ij2,ij3,ij4)
!  write(6,'(A,G20.10)') 'Dkl=',XDot(Dkl,kl1,kl2,kl3,kl4)
!  write(6,'(A,G20.10)') 'Dik=',XDot(Dik,ik1,ik2,ik3,ik4)
!  write(6,'(A,G20.10)') 'Dil=',XDot(Dil,il1,il2,il3,il4)
!  write(6,'(A,G20.10)') 'Djk=',XDot(Djk,jk1,jk2,jk3,jk4)
!  write(6,'(A,G20.10)') 'Djl=',XDot(Djl,jl1,jl2,jl3,jl4)
!end if

!write(6,'(A,8L1)') 'Pert=',Pert
if (iBas*jBas*kBas*lBas > nScrt) then
  write(6,*) 'FckAcc_McK: iBas*jBas*kBas*lBas > nScrt'
  write(6,*) 'iBas,jBas,kBas,lBas,nScrt=',iBas,jBas,kBas,lBas,nScrt
  call Abend()
end if
ii = iOff(iAng(1))
jj = iOff(iAng(2))
kk = iOff(iAng(3))
ll = iOff(iAng(4))
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

ipFij = 1
nFij = iBas*jBas*iCmpa(1)*iCmpa(2)

ipFkl = ipFij+nFij
nFkl = kBas*lBas*iCmpa(3)*iCmpa(4)

ipFik = ipFkl+nFkl
nFik = iBas*kBas*iCmpa(1)*iCmpa(3)

ipFjl = ipFik+nFik
nFjl = jBas*lBas*iCmpa(2)*iCmpa(4)

ipFil = ipFjl+nFjl
nFil = iBas*lBas*iCmpa(1)*iCmpa(4)

ipFjk = ipFil+nFil
nFjk = jBas*kBas*iCmpa(2)*iCmpa(3)

call dcopy_(nFij+nFkl+nFik+nFjl+nFil+nFjk,[Zero],0,FT(ipFij),1)

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

Shij = iShell(1) == iShell(2)
Shkl = iShell(3) == iShell(4)
iShij = iShell(1) == iShell(2)
iShkl = iShell(3) == iShell(4)
iShik = iShell(1) == iShell(3)
iShil = iShell(1) == iShell(4)
iShjk = iShell(2) == iShell(3)
iShjl = iShell(2) == iShell(4)
do i1=1,iCmp
  do j=0,nIrrep-1
    iSym(1,j) = 0
    if (iAOtSO(iAO(1)+i1,j) > 0) iSym(1,j) = iTwoj(j)
  end do
  jCmpMx = jCmp
  if (Shij) jCmpMx = i1
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pEa = xPrmt(iOper(kOp(1)),iChBs)
  do i2=1,jCmpMx
    do j=0,nIrrep-1
      iSym(2,j) = 0
      if (iAOtSO(iAO(2)+i2,j) > 0) iSym(2,j) = iTwoj(j)
    end do
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pRb = xPrmt(iOper(kOp(2)),jChBs)
    !Qij = i1 == i2
    if (iShell(2) > iShell(1)) then
      i12 = jCmp*(i1-1)+i2
    else
      i12 = iCmp*(i2-1)+i1
    end if
    do i3=1,kCmp
      do j=0,nIrrep-1
        iSym(3,j) = 0
        if (iAOtSO(iAO(3)+i3,j) > 0) iSym(3,j) = iTwoj(j)
      end do
      lCmpMx = lCmp
      if (Shkl) lCmpMx = i3
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pTc = xPrmt(iOper(kOp(3)),kChBs)
      do i4=1,lCmpMx
        do j=0,nIrrep-1
          iSym(4,j) = 0
          if (iAOtSO(iAO(4)+i4,j) > 0) iSym(4,j) = iTwoj(j)
        end do
        !Qkl = i3 == i4
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        pTSd = xPrmt(iOper(kOp(4)),lChBs)
        if (iShell(4) > iShell(3)) then
          i34 = lCmp*(i3-1)+i4
        else
          i34 = kCmp*(i4-1)+i3
        end if
        if (Shijij .and. (i34 > i12)) cycle
        Qijij = Shijij .and. (i12 == i34)
        iQij = iShij .and. (i1 == i2)
        iQkl = iShkl .and. (i3 == i4)
        iQik = iShik .and. (i1 == i3)
        iQil = iShil .and. (i1 == i4)
        iQjk = iShjk .and. (i2 == i3)
        iQjl = iShjl .and. (i2 == i4)
        pFctr = pEa*pRb*pTc*pTSd

        mFij = 0
        mFkl = 0
        mFik = 0
        mFjl = 0
        mFil = 0
        mFjk = 0
        do iIrrep=0,nIrrep-1
          if ((iSym(1,iIrrep) /= 0) .and. (iSym(2,iIrrep) /= 0)) mFkl = mFkl+1
          if ((iSym(3,iIrrep) /= 0) .and. (iSym(4,iIrrep) /= 0)) mFij = mFij+1
          if ((iSym(1,iIrrep) /= 0) .and. (iSym(3,iIrrep) /= 0)) mFjl = mFjl+1
          if ((iSym(2,iIrrep) /= 0) .and. (iSym(4,iIrrep) /= 0)) mFik = mFik+1
          if ((iSym(1,iIrrep) /= 0) .and. (iSym(4,iIrrep) /= 0)) mFjk = mFjk+1
          if ((iSym(2,iIrrep) /= 0) .and. (iSym(3,iIrrep) /= 0)) mFil = mFil+1
        end do
        if (mFij+mFkl+mFik+mFjl+mFil+mFjk == 0) cycle

        vijkl = DNrm2_(iBas*jBas*kBas*lBas,AOInt(1,i1,i2,i3,i4),1)
        if (vijkl < CutInt) cycle
        !***************************************************************
        !                                                              *
        ! Fij = hij + Sum(kl) Dkl {(ij|kl)-1/2(ik|jl)}                 *
        !                                                              *
        ! or                                                           *
        !                                                              *
        ! Fij = hij + Sum(k=>l) Dkl {(2-d(kl)} P(ij|kl)                *
        !                                                              *
        ! where P(ij|kl)=(ij|kl)-1/4(ik|jl)-1/4(il|jk)                 *
        !                                                              *
        ! or in the case of no sum restriction                         *
        !                                                              *
        ! P(ij|kl) = (ij|kl) - 1/2(ik|jl)                              *
        !                                                              *
        ! Coloumb Contributions                                        *
        !                                                              *
        ! Fij = Dkl * (ij|kl)                                          *
        !                                                              *
        ! Order density matrix in accordance with the integrals        *
        !                                                              *
        !***************************************************************
        if (mFij /= 0) then
          if (iShell(3) < iShell(4)) then
            vkl = Dkl(kBas*lBas+1,i4,i3)
          else
            vkl = Dkl(kBas*lBas+1,i3,i4)
          end if

          ! Pickup the right column of the density matrix and
          ! change order if not canonical.

          qFctr = One
          ipFij1 = ((i2-1)*iCmpa(1)+i1-1)*iBas*jBas+ipFij
          Fac = One
          if (iQij) Fac = Half
          D_kl = Two
          if (iQkl) D_kl = One
          Fac = Fac*D_kl
          if (vkl*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
            if (iShell(3) < iShell(4)) then
              call DGetMO(Dkl(1,i4,i3),kl1,kl1,kl2,Scrt,kl2)
              call dGeMV_('N',iBas*jBas,kBas*lBas,Fac*qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,Scrt,1,One,FT(ipFij1),1)
            else
              call dGeMV_('N',iBas*jBas,kBas*lBas,Fac*qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,Dkl(1,i3,i4),1,One,FT(ipFij1),1)
            end if
            !call RecPrt('Fij',' ',FT(ipFij1),iBas,jBas)
            lFij = .true.
          end if
        end if
        ! Fkl = Dij * (ij|kl)
        if ((.not. Qijij) .and. (mFkl /= 0)) then
          if (iShell(1) < iShell(2)) then
            vij = Dij(iBas*jBas+1,i2,i1)
          else
            vij = Dij(iBas*jBas+1,i1,i2)
          end if
          qFctr = One
          ipFkl1 = ((i4-1)*iCmpa(3)+i3-1)*kBas*lBas+ipFkl
          Fac = One
          if (iQkl) Fac = Half
          D_ij = Two
          if (iQij) D_ij = One
          Fac = Fac*D_ij
          if (vij*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
            if (iShell(1) < iShell(2)) then
              call DGeTMO(Dij(1,i2,i1),ij1,ij1,ij2,Scrt,ij2)
              call dGeMV_('T',iBas*jBas,kBas*lBas,Fac*qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,Scrt,1,One,FT(ipFkl1),1)
            else
              call dGeMV_('T',iBas*jBas,kBas*lBas,Fac*qFctr*pFctr,AOInt(1,i1,i2,i3,i4),iBas*jBas,Dij(1,i1,i2),1,One,FT(ipFkl1),1)
            end if
            !call RecPrt('Fkl',' ',FT(ipFkl1),kBas,lBas)
            lFkl = .true.
          end if
        end if

        ! Exchange contributions

        ! Change the order ijkl to ikjl. Make sure also that
        ! the index pairs are canonical.

        ipD = 1+iBas*jBas*kBas*lBas
        np = ipD-1+max(iBas*lBas,jBas*lBas,iBas*kBas,jBas*kBas)
        if (np > nScrt) then
          write(6,*) 'FckAcc_McK: np > nScrt'
          write(6,*) 'np,nScrt=',np,nScrt
          call Abend()
        end if
        if (mFik+mFjl /= 0) then
          if ((mFik /= 0) .and. (iShell(2) < iShell(4))) then
            vjl = Djl(jBas*lBas+1,i4,i2)
          else if (mFik /= 0) then
            vjl = Djl(jBas*lBas+1,i2,i4)
          else
            vjl = Zero
          end if
          if ((mFjl /= 0) .and. (iShell(1) < iShell(3))) then
            vik = Dik(iBas*kBas+1,i3,i1)
          else if (mFjl /= 0) then
            vik = Dik(iBas*kBas+1,i1,i3)
          else
            vik = Zero
          end if
          if ((vik*vijkl/Four >= CutInt) .or. (vjl*vijkl/Four >= CutInt)) then
            do j=1,jBas
              nij = (j-1)*iBas+1
              do k=1,kBas
                nik = (k-1)*iBas+1
                do l=1,lBas
                  nkl = (l-1)*kBas+k
                  njl = (l-1)*jBas+j
                  iOut = (nkl-1)*iBas*jBas+nij
                  iIn = (njl-1)*iBas*kBas+nik
                  call dcopy_(iBas,AOInt(iOut,i1,i2,i3,i4),1,Scrt(iIn),1)
                end do
              end do
            end do
            !***********************************************************
            !                                                          *
            ! Fik = - 1/4 * Djl * P(jl|ik)                             *
            !                                                          *
            ! P(jl|ik) = (jl|ik) - 1/4(ji|lk) - 1/4(jk|il)             *
            !                                                          *
            ! P(jl|ik) = (jl|ik) - 1/2(ji|lk)                          *
            !                                                          *
            ! Change factor if                                         *
            ! a) asymmetrical P matrix is implied                      *
            ! b) if the two exchange integrals in the symmetrical      *
            !    P matrix are identical.                               *
            !                                                          *
            !***********************************************************
            if (mFik /= 0) then
              qFctr = One
              ipFik1 = ((i3-1)*iCmpa(1)+i1-1)*iBas*kBas+ipFik
              Fac = -Quart
              if (iQjl .and. (.not. iQik)) Fac = -Half
              D_jl = Two
              if (iQjl) D_jl = One
              Fac = Fac*D_jl
              if (vjl*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
                if (iShell(2) < iShell(4)) then
                  call DGeTMO(Djl(1,i4,i2),jl1,jl1,jl2,Scrt(ipD),jl2)
                  call dGeMV_('N',iBas*kBas,jBas*lBas,Fac*qFctr*pFctr,Scrt,iBas*kBas,Scrt(ipD),1,One,FT(ipFik1),1)
                else
                  call dGeMV_('N',iBas*kBas,jBas*lBas,Fac*qFctr*pFctr,Scrt,iBas*kBas,Djl(1,i2,i4),1,One,FT(ipFik1),1)
                end if
                !call RecPrt('Fik',' ',FT(ipFik1),iBas,kBas)
                lFik = .true.
              end if
            end if
            if ((.not. iQij) .or. (.not. iQkl)) then
              !*********************************************************
              !                                                        *
              ! Fjl = - 1/4 * Dik * P(jl|ik)                           *
              !                                                        *
              ! P(jl|ik) = (jl|ik) - 1/4(ji|lk) - 1/4(jk|il)           *
              !                                                        *
              ! P(jl|ik) = (jl|ik) - 1/2(ji|lk)                        *
              !                                                        *
              !*********************************************************
              if (mFjl /= 0) then
                qFctr = One
                ipFjl1 = ((i4-1)*iCmpa(2)+i2-1)*jBas*lBas+ipFjl
                Fac = -Quart
                if (iQik .and. (.not. iQjl)) Fac = -Half
                D_ik = Two
                if (iQik) D_ik = One
                Fac = Fac*D_ik
                if (vik*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
                  if (iShell(1) < iShell(3)) then
                    call DGeTMO(Dik(1,i3,i1),ik1,ik1,ik2,Scrt(ipD),ik2)
                    call dGeMV_('T',iBas*kBas,jBas*lBas,Fac*qFctr*pFctr,Scrt,iBas*kBas,Scrt(ipD),1,One,FT(ipFjl1),1)
                  else
                    call dGeMV_('T',iBas*kBas,jBas*lBas,Fac*qFctr*pFctr,Scrt,iBas*kBas,Dik(1,i1,i3),1,One,FT(ipFjl1),1)
                  end if
                  !call RecPrt('Fjl',' ',FT(ipFjl1),jBas,lBas)
                  lFjl = .true.
                end if
              end if
            end if
          end if
        end if

        ! Change order ijkl to iljk

        if ((.not. iQij) .and. (.not. iQkl) .and. (mFil+mFjk /= 0)) then
          if ((mFil /= 0) .and. (iShell(2) < iShell(3))) then
            vjk = Djk(jBas*kBas+1,i3,i2)
          else if (mFil /= 0) then
            vjk = Djk(jBas*kBas+1,i2,i3)
          else
            vjk = Zero
          end if
          if ((mFjk /= 0) .and. (iShell(1) < iShell(4))) then
            vil = Dil(iBas*lBas+1,i4,i1)
          else if (mFjk /= 0) then
            vil = Dil(iBas*lBas+1,i1,i4)
          else
            vil = Zero
          end if
          if ((vil*vijkl/Four >= CutInt) .or. (vjk*vijkl/Four >= CutInt)) then
            i = 1
            do j=1,jBas
              nij = (j-1)*iBas+i
              do k=1,kBas
                njk = (k-1)*jBas+j
                do l=1,lBas
                  nkl = (l-1)*kBas+k
                  nil = (l-1)*iBas+i
                  ijkl = (nkl-1)*iBas*jBas+nij
                  iljk = (njk-1)*iBas*lBas+nil
                  call dcopy_(iBas,AOInt(ijkl,i1,i2,i3,i4),1,Scrt(iljk),1)
                end do
              end do
            end do
            ! Fil = - 1/4 * Djk * (ij|kl)
            if (mFil /= 0) then
              qFctr = One
              ipFil1 = ((i4-1)*iCmpa(1)+i1-1)*iBas*lBas+ipFil
              Fac = -Quart
              if (iQjk .and. (.not. iQil)) Fac = -Half
              D_jk = Two
              if (iQjk) D_jk = One
              Fac = Fac*D_jk
              if (vjk*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
                if (iShell(2) < iShell(3)) then
                  call DGeTMO(Djk(1,i3,i2),jk1,jk1,jk2,Scrt(ipD),jk2)
                  call dGeMV_('N',iBas*lBas,jBas*kBas,Fac*qFctr*pFctr,Scrt,iBas*lBas,Scrt(ipD),1,One,FT(ipFil1),1)
                else
                  call dGeMV_('N',iBas*lBas,jBas*kBas,Fac*qFctr*pFctr,Scrt,iBas*lBas,Djk(1,i2,i3),1,One,FT(ipFil1),1)
                end if
                !call RecPrt('Fil',' ',FT(ipFil1),iBas,lBas)
                lFil = .true.
              end if
            end if
            ! Fjk = - 1/4 * Dil * (ij|kl)
            if ((.not. Qijij) .and. (mFjk /= 0)) then
              qFctr = One
              ipFjk1 = ((i3-1)*iCmpa(2)+i2-1)*jBas*kBas+ipFjk
              Fac = -Quart
              if (iQil .and. (.not. iQjk)) Fac = -Half
              D_il = Two
              if (iQil) D_il = One
              Fac = Fac*D_il
              if (vil*vijkl*abs(Fac*qFctr*pFctr) >= CutInt) then
                if (iShell(1) < iShell(4)) then
                  call DGeTMO(Dil(1,i4,i1),il1,il1,il2,Scrt(ipD),il2)
                  call dGeMV_('T',iBas*lBas,jBas*kBas,Fac*qFctr*pFctr,Scrt,iBas*lBas,Scrt(ipD),1,One,FT(ipFjk1),1)
                else
                  call dGeMV_('T',iBas*lBas,jBas*kBas,Fac*qFctr*pFctr,Scrt,iBas*lBas,Dil(1,i1,i4),1,One,FT(ipFjk1),1)
                end if
                !call RecPrt('Fjk',' ',FT(ipFjk1),jBas,kBas)
                lFjk = .true.
              end if
            end if
          end if
        end if

      end do
    end do
  end do
end do

!write(6,*) ' Fij'
nnIrrep = nIrrep
if (sIrrep) nnIrrep = 1
do iIrrep=0,nnIrrep-1
  !write(6,'(I2,L1)') iIrrep,pert(iIrrep)
  if (pert(iIrrep)) then
    ip = ipDisp(abs(indgrd(iCar,iCent,iIrrep)))
    rCh = xPrmt(iOper(kOp(iCent)),iChBas(1+iCar))*dble(iChTbl(iIrrep,kOp(iCent)))
    Fact = tfact*rCh
    !write(6,*) 'Level ij'
    if (lFij) call FckDst(TwoHam(ip),ndens,FT(ipFij),iBas,jBas,iCmpa(1),iCmpa(2),kOp2(1),kOp2(2),iIrrep,iShij,iAO(1),iAO(2), &
                          iAOst(1),iAOst(2),fact)
    !write(6,*) 'Level kl'
    if (lFkl) call FckDst(TwoHam(ip),ndens,FT(ipFkl),kBas,lBas,iCmpa(3),iCmpa(4),kOp2(3),kOp2(4),iIrrep,iShkl,iAO(3),iAO(4), &
                          iAOst(3),iAOst(4),fact)
    !write(6,*) 'Level ik'
    if (lFik) call FckDst(TwoHam(ip),ndens,FT(ipFik),iBas,kBas,iCmpa(1),iCmpa(3),kOp2(1),kOp2(3),iIrrep,iShik,iAO(1),iAO(3), &
                          iAOst(1),iAOst(3),fact)
    !write(6,*) 'Level jl'
    if (lFjl) call FckDst(TwoHam(ip),ndens,FT(ipFjl),jBas,lBas,iCmpa(2),iCmpa(4),kOp2(2),kOp2(4),iIrrep,iShjl,iAO(2),iAO(4), &
                          iAOst(2),iAOst(4),fact)
    !write(6,*) 'Level il'
    if (lFil) call FckDst(TwoHam(ip),ndens,FT(ipFil),iBas,lBas,iCmpa(1),iCmpa(4),kOp2(1),kOp2(4),iIrrep,iShil,iAO(1),iAO(4), &
                          iAOst(1),iAOst(4),fact)
    !write(6,*) 'Level jk'
    if (lFjk) call FckDst(TwoHam(ip),ndens,FT(ipFjk),jBas,kBas,iCmpa(2),iCmpa(3),kOp2(2),kOp2(3),iIrrep,iShjk,iAO(2),iAO(3), &
                          iAOst(2),iAOst(3),fact)
    !if (DDot_(3468,TwoHam,1,TwoHam,1) > Zero) then
    !  if (Abs(DDot_(3468,TwoHam,1,One,0)) > 1.0D-16) then
    !    write(6,'(A,G20.6,G20.7)') 'TwoHam=',DDot_(3468,TwoHam,1,TwoHam,1),DDot_(3468,TwoHam,1,One,0)
    !  else
    !    write(6,'(A,G20.6)') 'TwoHam=',DDot_(3468,TwoHam,1,TwoHam,1)
    !  end if
    !  call FZero(TwoHam,3468)
    !end if
  end if
end do

return

end subroutine FckAcc_Mck

#else

subroutine FckAcc_Mck(iAng,iCmp,jCmp,kCmp,lCmp,Shijij,iShll,iShell,kOp,nijkl,AOInt,TwoHam,nDens,Scrt,nScrt,iAO,iAOst,iBas,jBas, &
                      kBas,lBas,Dij,ij1,ij2,ij3,ij4,Dkl,kl1,kl2,kl3,kl4,Dik,ik1,ik2,ik3,ik4,Dil,il1,il2,il3,il4,Djk,jk1,jk2,jk3, &
                      jk4,Djl,jl1,jl2,jl3,jl4,FT,nFT,tfact,iCar,iCent,pert,indgrd,ipdisp)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmatry adapted Fock matrix.                        *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          In addition to this complication we have that the order of  *
!          indicies in the integrals are not ordered canonically but   *
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

use Basis_Info
use Symmetry_Info, only: nIrrep, iChTbl, iChBas
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Gateway_Info, only: ThrInt, CutInt

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "disp2.fh"
#include "print.fh"
real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), TwoHam(nDens), Scrt(nScrt), FT(nFT), Dij(ij1*ij2+1,ij3,ij4), Dkl(kl1*kl2+1,kl3,kl4), &
       Dik(ik1*ik2+1,ik3,ik4), Dil(il1*il2+1,il3,il4), Djk(jk1*jk2+1,jk3,jk4), Djl(jl1*jl2+1,jl3,jl4)
logical Shijij, Qijij, iShij, iShkl, iQij, iQkl, iQik, iShik, iQil, iShil, iQjk, iShjk, iQjl, iShjl, lFij, lFkl, lFik, lFjl, lFil, &
        lFjk
integer iAng(4), iShell(4), iShll(4), kOp(4), kOp2(4), iAO(4), iAOst(4), iCmpa(4)
logical Pert(0:nIrrep-1)
integer indgrd(3,4,0:nirrep-1), ipdisp(*)
!Local Arrays
integer iSym(4)
real*8 Prmt(0:7)
!character*72 Label
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Function
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
xPrmt(i,j) = Prmt(iand(i,j))
!iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

!iRout = 38
!iPrint = nPrint(iRout)
!
!if (iPrint >= 49) then
!   write(6,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
!end if
!if (iPrint >= 99) then
!   call RecPrt('FckAcc:AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
!   write(6,*) 'Dij=',XDot(Dij,ij1,ij2,ij3,ij4)
!   write(6,*) 'Dkl=',XDot(Dkl,kl1,kl2,kl3,kl4)
!   write(6,*) 'Dik=',XDot(Dik,ik1,ik2,ik3,ik4)
!   write(6,*) 'Dil=',XDot(Dil,il1,il2,il3,il4)
!   write(6,*) 'Djk=',XDot(Djk,jk1,jk2,jk3,jk4)
!   write(6,*) 'Djl=',XDot(Djl,jl1,jl2,jl3,jl4)
!end if

ExFac = One
ThrInt = 0.0d0
if (iBas*jBas*kBas*lBas > nScrt) then
  write(6,*) 'FckAcc: nScrt too small!'
  call Abend()
end if
ii = iOff(iAng(1))
jj = iOff(iAng(2))
kk = iOff(iAng(3))
ll = iOff(iAng(4))

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

ipFij = 1
nFij = iBas*jBas*iCmpa(1)*iCmpa(2)

ipFkl = ipFij+nFij
nFkl = kBas*lBas*iCmpa(3)*iCmpa(4)

ipFik = ipFkl+nFkl
nFik = iBas*kBas*iCmpa(1)*iCmpa(3)

ipFjl = ipFik+nFik
nFjl = jBas*lBas*iCmpa(2)*iCmpa(4)

ipFil = ipFjl+nFjl
nFil = iBas*lBas*iCmpa(1)*iCmpa(4)

ipFjk = ipFil+nFil
nFjk = jBas*kBas*iCmpa(2)*iCmpa(3)

call dcopy_(nFij+nFkl+nFik+nFjl+nFil+nFjk,Zero,0,FT(ipFij),1)

ipDij = 1
ipDkl = 1
ipDik = 1
ipDil = 1
ipDjk = 1
ipDjl = 1
ipFij1 = 1
ipFkl1 = 1
ipFik1 = 1
ipFil1 = 1
ipFjk1 = 1
ipFjl1 = 1

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

iShij = iShell(1) == iShell(2)
iShkl = iShell(3) == iShell(4)
iShik = iShell(1) == iShell(3)
iShil = iShell(1) == iShell(4)
iShjk = iShell(2) == iShell(3)
iShjl = iShell(2) == iShell(4)
mijkl = iBas*jBas*kBas*lBas
do i1=1,iCmp
  ix = 0
  do j=0,nIrrep-1
    if (iAOtSO(iAO(1)+i1,j) > 0) ix = ior(ix,2**j)
  end do
  iSym(1) = ix
  jCmpMx = jCmp
  if (iShij) jCmpMx = i1
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pEa = xPrmt(iOper(kOp(1)),iChBs)
  do i2=1,jCmpMx
    ix = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) ix = ior(ix,2**j)
    end do
    iSym(2) = ix
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pRb = xPrmt(iOper(kOp(2)),jChBs)
    if (iShell(2) > iShell(1)) then
      i12 = jCmp*(i1-1)+i2
    else
      i12 = iCmp*(i2-1)+i1
    end if
    do i3=1,kCmp
      ix = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(3)+i3,j) > 0) ix = ior(ix,2**j)
      end do
      iSym(3) = ix
      lCmpMx = lCmp
      if (iShkl) lCmpMx = i3
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pTc = xPrmt(iOper(kOp(3)),kChBs)
      do i4=1,lCmpMx
        ix = 0
        do j=0,nIrrep-1
          if (iAOtSO(iAO(4)+i4,j) > 0) ix = ior(ix,2**j)
        end do
        iSym(4) = ix
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        pTSd = xPrmt(iOper(kOp(4)),lChBs)
        if (iShell(4) > iShell(3)) then
          i34 = lCmp*(i3-1)+i4
        else
          i34 = kCmp*(i4-1)+i3
        end if
        if (Shijij .and. (i34 > i12)) cycle
        vijkl = 0
        do ijkl=1,mijkl
          vijkl = max(vijkl,abs(AOInt(ijkl,i1,i2,i3,i4)))
        end do
        !vijkl = DNrm2_(iBas*jBas*kBas*lBas,AOInt(1,i1,i2,i3,i4),1)
        !if (vijkl < CutInt) cycle

        Qijij = Shijij .and. (i12 == i34)
        iQij = iShij .and. (i1 == i2)
        iQkl = iShkl .and. (i3 == i4)
        iQik = iShik .and. (i1 == i3)
        iQil = iShil .and. (i1 == i4)
        iQjk = iShjk .and. (i2 == i3)
        iQjl = iShjl .and. (i2 == i4)
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

        !write(6,*)
        !write(6,*) 'iShell(1),iShell(2),i1,i2=',iShell(1),iShell(2),i1,i2
        !write(6,*) 'Dij=',Dij(1,i1,i2)
        !write(6,*) 'Fac_ij,iQij=',Fac_ij,iQij
        !write(6,*)
        !write(6,*) 'iShell(3),iShell(4),i3,i4=',iShell(3),iShell(4),i3,i4
        !write(6,*) 'Dkl=',Dkl(1,i3,i4)
        !write(6,*) 'Fac_kl,iQkl=',Fac_kl,iQkl
        !write(6,*)
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
        if ((iand(iSym(1),iSym(2)) /= 0) .or. (iand(iSym(3),iSym(4)) /= 0)) then
          iOpt = iOpt+1

          if (iShell(3) < iShell(4)) then
            vkl = Dkl(kBas*lBas+1,i4,i3)
            ipDkl = ip
            ip = ip+kBas*lBas
            call DGetMO(Dkl(1,i4,i3),kl1,kl1,kl2,Scrt(ipDkl),kl2)
          else
            vkl = Dkl(kBas*lBas+1,i3,i4)
            loc1 = (ixLoc(Dkl(1,i3,i4))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDkl = 1+loc1/loc2
          end if
          ipFij1 = ((i2-1)*iCmpa(1)+i1-1)*iBas*jBas+ipFij
          if (iShell(1) < iShell(2)) then
            vij = Dij(iBas*jBas+1,i2,i1)
            ipDij = ip
            ip = ip+iBas*jBas
            call DGeTMO(Dij(1,i2,i1),ij1,ij1,ij2,Scrt(ipDij),ij2)
          else
            vij = Dij(iBas*jBas+1,i1,i2)
            loc1 = (ixLoc(Dij(1,i1,i2))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDij = 1+loc1/loc2
          end if
          ipFkl1 = ((i4-1)*iCmpa(3)+i3-1)*kBas*lBas+ipFkl
          if ((vkl*vijkl*abs(Fac_ij) < ThrInt) .and. (vij*vijkl*abs(Fac_kl) < ThrInt)) then
            iOpt = iOpt-1
          else
            lFij = .true.
            lFkl = .true.
          end if
        end if

        if ((iand(iSym(1),iSym(3)) /= 0) .or. (iand(iSym(2),iSym(4)) /= 0)) then
          iOpt = iOpt+2

          if (iShell(2) < iShell(4)) then
            vjl = Djl(jBas*lBas+1,i4,i2)
            ipDjl = ip
            ip = ip+jBas*lBas
            call DGeTMO(Djl(1,i4,i2),jl1,jl1,jl2,Scrt(ipDjl),jl2)
          else
            vjl = Djl(jBas*lBas+1,i2,i4)
            loc1 = (ixLoc(Djl(1,i2,i4))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDjl = 1+loc1/loc2
          end if
          ipFik1 = ((i3-1)*iCmpa(1)+i1-1)*iBas*kBas+ipFik
          if (iShell(1) < iShell(3)) then
            vik = Dik(iBas*kBas+1,i3,i1)
            ipDik = ip
            ip = ip+iBas*kBas
            call DGeTMO(Dik(1,i3,i1),ik1,ik1,ik2,Scrt(ipDik),ik2)
          else
            vik = Dik(iBas*kBas+1,i1,i3)
            loc1 = (ixLoc(Dik(1,i1,i3))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDik = 1+loc1/loc2
          end if
          ipFjl1 = ((i4-1)*iCmpa(2)+i2-1)*jBas*lBas+ipFjl
          if ((vik*vijkl*abs(Fac_jl) < ThrInt) .and. (vjl*vijkl*abs(Fac_ik) < ThrInt)) then
            iOpt = iOpt-2
          else
            lFik = .true.
            lFjl = .true.
          end if
        end if

        if ((iand(iSym(1),iSym(4)) /= 0) .or. (iand(iSym(2),iSym(3)) /= 0)) then
          iOpt = iOpt+4

          if (iShell(2) < iShell(3)) then
            vjk = Djk(jBas*kBas+1,i3,i2)
            ipDjk = ip
            ip = ip+jBas*kBas
            call DGeTMO(Djk(1,i3,i2),jk1,jk1,jk2,Scrt(ipDjk),jk2)
          else
            vjk = Djk(jBas*kBas+1,i2,i3)
            loc1 = (ixLoc(Djk(1,i2,i3))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDjk = 1+loc1/loc2
          end if
          ipFil1 = ((i4-1)*iCmpa(1)+i1-1)*iBas*lBas+ipFil
          if (iShell(1) < iShell(4)) then
            vil = Dil(iBas*lBas+1,i4,i1)
            ipDil = ip
            ip = ip+iBas*lBas
            call DGeTMO(Dil(1,i4,i1),il1,il1,il2,Scrt(ipDil),il2)
          else
            vil = Dil(iBas*lBas+1,i1,i4)
            loc1 = (ixLoc(Dil(1,i1,i4))-ixLoc(Scrt))
            loc2 = (ixLoc(Scrt(2))-ixLoc(Scrt(1)))
            ipDil = 1+loc1/loc2
          end if
          ipFjk1 = ((i3-1)*iCmpa(2)+i2-1)*jBas*kBas+ipFjk
          if ((vil*vijkl*abs(Fac_jk) < ThrInt) .and. (vjk*vijkl*abs(Fac_il) < ThrInt)) then
            iOpt = iOpt-4
          else
            lFil = .true.
            lFjk = .true.
          end if
        end if
        if (ip-1 > nScrt) then
          write(6,*) 'FckAcc: nScrt too small!'
          call Abend()
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        select case (iOpt)
          case (1)
            call Fck1(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDij),FT(ipFij1),Fac_ij,Scrt(ipDkl),FT(ipFkl1),Fac_kl,ExFac)
          case (2)
            call Fck2(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDik),FT(ipFik1),Fac_ik,Scrt(ipDjl),FT(ipFjl1),Fac_jl,ExFac)
          case (3)
            call Fck3(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDij),FT(ipFij1),Fac_ij,Scrt(ipDkl),FT(ipFkl1),Fac_kl, &
                      Scrt(ipDik),FT(ipFik1),Fac_ik,Scrt(ipDjl),FT(ipFjl1),Fac_jl,ExFac)
          case (4)
            call Fck4(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDil),FT(ipFil1),Fac_il,Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
          case (5)
            call Fck5(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDij),FT(ipFij1),Fac_ij,Scrt(ipDkl),FT(ipFkl1),Fac_kl, &
                      Scrt(ipDil),FT(ipFil1),Fac_il,Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
          case (6)
            call Fck6(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDik),FT(ipFik1),Fac_ik,Scrt(ipDjl),FT(ipFjl1),Fac_jl, &
                      Scrt(ipDil),FT(ipFil1),Fac_il,Scrt(ipDjk),FT(ipFjk1),Fac_jk,ExFac)
          case (7)
            call Fck7(AOInt(1,i1,i2,i3,i4),iBas,jBas,kBas,lBas,Scrt(ipDij),FT(ipFij1),Fac_ij,Scrt(ipDkl),FT(ipFkl1),Fac_kl, &
                      Scrt(ipDik),FT(ipFik1),Fac_ik,Scrt(ipDjl),FT(ipFjl1),Fac_jl,Scrt(ipDil),FT(ipFil1),Fac_il,Scrt(ipDjk), &
                      FT(ipFjk1),Fac_jk,ExFac)
          case default
            ! nothing
        end select
        !***************************************************************

      end do
    end do
  end do
end do

nnIrrep = nIrrep
if (sIrrep) nnIrrep = 1
do iIrrep=0,nnIrrep-1

  if (pert(iIrrep)) then
    ip = ipDisp(abs(indgrd(iCar,iCent,iIrrep)))
    rCh = xPrmt(iOper(kOp(iCent)),iChBas(1+iCar))*dble(iChTbl(iIrrep,kOp(iCent)))
    Fact = tfact*rCh
    !write(6,*) 'Level ij'
    if (lFij) call FckDst(TwoHam(ip),ndens,FT(ipFij),iBas,jBas,iCmpa(1),iCmpa(2),kOp2(1),kOp2(2),iIrrep,iShij,iAO(1),iAO(2), &
                          iAOst(1),iAOst(2),fact)
    !write(6,*) 'Level kl'
    if (lFkl) call FckDst(TwoHam(ip),ndens,FT(ipFkl),kBas,lBas,iCmpa(3),iCmpa(4),kOp2(3),kOp2(4),iIrrep,iShkl,iAO(3),iAO(4), &
                          iAOst(3),iAOst(4),fact)
    !write(6,*) 'Level ik'
    if (lFik) call FckDst(TwoHam(ip),ndens,FT(ipFik),iBas,kBas,iCmpa(1),iCmpa(3),kOp2(1),kOp2(3),iIrrep,iShik,iAO(1),iAO(3), &
                          iAOst(1),iAOst(3),fact)
    !write(6,*) 'Level jl'
    if (lFjl) call FckDst(TwoHam(ip),ndens,FT(ipFjl),jBas,lBas,iCmpa(2),iCmpa(4),kOp2(2),kOp2(4),iIrrep,iShjl,iAO(2),iAO(4), &
                          iAOst(2),iAOst(4),fact)
    !write(6,*) 'Level il'
    if (lFil) call FckDst(TwoHam(ip),ndens,FT(ipFil),iBas,lBas,iCmpa(1),iCmpa(4),kOp2(1),kOp2(4),iIrrep,iShil,iAO(1),iAO(4), &
                          iAOst(1),iAOst(4),fact)
    !write(6,*) 'Level jk'
    if (lFjk) call FckDst(TwoHam(ip),ndens,FT(ipFjk),jBas,kBas,iCmpa(2),iCmpa(3),kOp2(2),kOp2(3),iIrrep,iShjk,iAO(2),iAO(3), &
                          iAOst(2),iAOst(3),fact)
  end if
end do

return

end subroutine FckAcc_mck

#endif
