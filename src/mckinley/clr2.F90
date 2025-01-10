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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine Clr2(XrIn,rOut,nACO,nSD,iSDi,iSDj,nDisp,nTemp,Temp)

! nACO: number of active orbitals

use McKinley_global, only: ipDisp3, ipMO
use Index_Functions, only: iTri, nTri_Elem
use pso_stuff, only: G2
use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: iOper, nIrrep
use Basis_Info, only: nBas
use Disp, only: lDisp
use Etwas, only: nASh
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: XrIn(*)
real(kind=wp), intent(inout) :: rOut(*)
integer(kind=iwp), intent(in) :: nACO, nSD, iSDi(0:nSD), iSDj(0:nSD), nDisp, nTemp
real(kind=wp), intent(inout) :: Temp(nTemp)

call Clr2_inner(XrIn,Temp)

! This is to avoid an explicit interface
contains

subroutine Clr2_inner(XrIn,Temp)

  real(kind=wp), intent(in), target :: XrIn(*)
  real(kind=wp), intent(inout), target :: Temp(nTemp)
  integer(kind=iwp) :: i, ia, iaoi, iaoj, iAsh, iB, ibas, iC, iCmp, id, iDisp, iE, ih, iiii, iij, iIrr, ij1, ij12, ij2, ipF, &
                       ipFKL, ipi, ipj, ipM, ipm2, ipp(0:7), iS, iShell(2), iSO, j, ja, jAsh, jB, jbas, jC, jCmp, jh, jIrr, jis, &
                       js, k, kAsh, kIrr, kl, kls, klt, l, lAsh, lIrr, lMax, lsl, lSO, mIrr, na(0:7), ni, nj, nnA, nX, nXrIn
  real(kind=wp) :: fact, rd
  logical(kind=iwp) :: Process
  real(kind=wp), pointer :: rIn(:,:,:,:,:), Temp1(:,:,:), Temp2(:), Temp3(:,:,:), Temp4(:,:,:), Temp5(:,:,:), Temp6(:)
  integer(kind=iwp), external :: NrOpr

  ibas = iSDi(3)
  jbas = iSDj(3)
  iCmp = iSDi(2)
  jCmp = iSDj(2)
  iAOi = iSDi(7)
  iAOj = iSDj(7)

  nXrIn = iBas*iCmp*jBas*jCmp*nIrrep*nTri_Elem(nACO)*nDisp

  ni = iBas*iCmp
  nj = jBas*jCmp

  rIn(1:ni,1:nj,0:nIrrep-1,1:nTri_Elem(nACO),1:nDisp) => XrIn(1:nXrIn)

  iS = 1
  iE = iBas*iCmp*nACO
  Temp1(1:iBas,1:iCmp,1:nACO) => Temp(iS:iE)

  iS = iE+1
  iE = iE+nACO**2
  Temp2(1:nACO**2) => Temp(iS:iE)
  Temp2(:) = Zero

  iS = iE+1
  iE = iE+jBas*jCmp*nACO
  Temp3(1:jBas,1:jCmp,1:nACO) => Temp(iS:iE)
  iS = iE+1
  iE = iE+iBas*iCmp*nACO
  Temp4(1:iBas,1:iCmp,1:nACO) => Temp(iS:iE)
  Temp4(:,:,:) = Zero
  iS = iE+1
  iE = iE+jBas*jCmp*nACO
  Temp5(1:jBas,1:jCmp,1:nACO) => Temp(iS:iE)
  Temp5(:,:,:) = Zero
  iS = iE+1

  nX = 1
  iE = iE+nX
  Temp6(1:nX) => Temp(iS:iE) ! dummy set

  Process = .false.
  call PckMo2(temp6,icmp,iBas,jcmp,jBas,iaoi,iaoj,Process,nX)
  iE = iE-1+nX
  Temp6(1:nX) => Temp(iS:iE)

  iShell(1) = iSDi(11)
  iShell(2) = iSDj(11)

  nnA = 0
  do iS=0,nIrrep-1
    nA(iS) = nNA
    nnA = nnA+nAsh(is)
  end do

  ipi = 1

  ipj = ipi+nACO*ibas*icmp

  Process = .true.
  call PckMo2(temp6,icmp,iBas,jcmp,jBas,iaoi,iaoj,Process,nX)

  id = 0
  do mIrr=0,nIrrep-1
    iiii = 0
    do iS=0,nIrrep-1
      js = nrOpr(ieor(ioper(is),iOper(mIrr)))
      ipp(is) = iiii
      iiii = nbas(is)*nash(js)+iiii
    end do
    do iDisp=1,lDisp(mIrr)
      iD = id+1
      ia = 1
      do iIrr=0,nIrrep-1
        kl = 0
        k = 0
        do kIrr=0,nIrrep-1
          do kAsh=1,nAsh(kIrr)
            k = k+1
            l = 0
            do lIrr=0,kIrr
              kls = ieor(iOper(kIrr),iOper(lIrr))
              jIrr = nropr(ieor(ieor(iOper(iIrr),iOper(mIrr)),kls))
              ja = 1
              do j=0,jIrr-1
                ja = ja+nAsh(j)
              end do

              ! Symmetry of Q matrix

              jis = nropr(ieor(iOper(iIrr),ioper(mIrr)))

              lMax = nAsh(lIrr)
              if (lIrr == kirr) lmax = kash
              do lAsh=1,lMax
                l = lash+nA(lIrr)
                kl = itri(k,l)

                ! id,iirr,jirr,kA,lA

                if (nash(jirr) /= 0) &
                  call DGEMM_('N','N',ni,nAsh(jIrr),nj, &
                              One,rIn(:,:,iIrr,kl,id),ni, &
                                  Temp6(ipj+(ja-1)*jcmp*jBas:nX),nj, &
                              Zero,Temp1,ni)
                if (nash(iirr) /= 0) &
                  call DGEMM_('T','N',nash(iIrr),nAsh(jIrr),ni, &
                              One,Temp6(ipi+(ia-1)*icmp*ibas:nX),ni, &
                                  Temp1,ni, &
                              Zero,Temp2,nash(iirr))

                do iC=1,iCmp
                  do iB=1,iBas
                    do i=1,nAsh(jis)
                      ih = i+na(jis)
                      Temp4(iB,ic,i) = Zero
                      do iAsh=1,nAsh(jirr)
                        jh = iash+na(jirr)
                        fact = One
                        iij = itri(ih,jh)
                        if ((iij >= kl) .and. (k == l)) fact = Two
                        if ((iij < kl) .and. (ih == jh)) fact = Two
                        if (k /= l) FacT = fact*Two
                        rd = G2(itri(iij,kl),1)*fact
                        Temp4(iB,ic,i) = Temp4(ib,ic,i)+Temp1(ib,ic,iash)*rd
                      end do
                    end do
                  end do
                end do

                ipF = ipDisp3(id)-1+ipp(iirr)
                do jAsh=1,nAsh(jis)
                  do iC=1,iCmp
                    lSO = iAOtSO(iAOi+iC,iIrr)
                    if (lso > 0) then
                      do iB=1,iBas
                        lsl = lSO+ib-1
                        ipFKL = ipF+(jAsh-1)*nBas(iIrr)+lsl
                        rOut(ipFKL) = rOut(ipFKL)+Temp4(ib,ic,jash)
                      end do
                    end if
                  end do
                end do

                if (iShell(1) /= iShell(2)) then
                  if (nash(jirr) /= 0) &
                    call DGEMM_('T','N',nj,nAsh(jIrr),ni, &
                                One,rIn(:,:,jIrr,kl,id),ni, &
                                    Temp6(ipi+(ja-1)*icmp*ibas:nX),ni, &
                                Zero,Temp3,nj)
                  if (nash(iirr) /= 0) &
                    call DGEMM_('T','N',nAsh(iirr),nAsh(jirr),nj, &
                                One,Temp6(ipj+(ia-1)*jcmp*jBas:nX),nj, &
                                    Temp3,nj, &
                                One,Temp2,nAsh(iirr))

                  do jC=1,jCmp
                    do jB=1,jBas
                      do i=1,nAsh(jis)
                        ih = i+na(jis)
                        Temp5(jB,jc,i) = Zero
                        do iAsh=1,nAsh(jirr)
                          jh = iash+na(jirr)
                          fact = One
                          iij = itri(ih,jh)
                          if ((iij >= kl) .and. (k == l)) fact = Two
                          if ((iij < kl) .and. (ih == jh)) fact = Two
                          if (k /= l) FacT = fact*Two
                          rd = G2(itri(iij,kl),1)*fact
                          Temp5(jB,jc,i) = Temp5(jb,jc,i)+Temp3(jb,jc,iash)*rd
                        end do
                      end do
                    end do
                  end do

                  ipf = ipDisp3(id)-1+ipp(iirr)
                  do iAsh=1,nAsh(jis)
                    do jC=1,jCmp
                      iSO = iAOtSO(iAOj+jC,iIrr)
                      if (iso > 0) then
                        do jB=1,jBas
                          i = iSO+jb-1
                          ipFKL = ipF+(iAsh-1)*nBas(iIrr)+i
                          rOut(ipFKL) = rOut(ipFKL)+Temp5(jb,jc,iash)
                        end do
                      end if
                    end do
                  end do
                end if

                ! Distribute integrals

                if (iirr >= jirr) then
                  klt = itri(k,l)
                  do iAsh=1,nash(iirr)
                    ij1 = iAsh+na(iirr)
                    do jAsh=1,nash(jirr)
                      ij2 = jAsh+na(jirr)
                      if (ij1 >= ij2) then
                        ij12 = itri(ij1,ij2)
                        if (ij12 <= klt) then
                          ipM = ipMO(iD)+iTri(ij12,klt)-1
                          ipm2 = iash+(jash-1)*nash(iirr)
                          rOut(ipm) = rOut(ipm)+Temp2(ipm2)
                        end if
                      end if
                    end do
                  end do
                end if
              end do ! lash
            end do ! lirr
          end do ! kash
        end do ! kirr
        ia = ia+nAsh(iIrr)
      end do ! iirr
    end do ! ndisp
  end do ! msym

end subroutine Clr2_inner

end subroutine Clr2
