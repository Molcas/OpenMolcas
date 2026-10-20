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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************
!
! See cgs_mod.F90 for details

subroutine CGS_x(nDensC,nConf1,nRoots,iSym,jspin,iter,ipCI,ipDia,ipS1,ipS2,ipST,ipCIT,ipCId,ipPre2,reco,Fancy,Kappa,dKappa,Sigma, &
                 Temp4,Sc1,iPre,delta,resk,resci,deltak,deltac)

use ipPage, only: ipIn, opOut, W
use ISRotation, only: DMInvISR, InvSCF, ISR
use cgs_mod, only: CGSvec, PCGS
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: nDensC, nConf1, nRoots, jspin, iPre(*)
integer(kind=iwp), intent(inout) :: iSym, iter
integer(kind=iwp), intent(in) :: ipCI, ipDia, ipS1, ipS2, ipST, ipCIT, ipCId, ipPre2
real(kind=wp), intent(inout) :: reco, Fancy(*), Kappa(*), dKappa(*), Sigma(*), Temp4(*), Sc1(*), delta, resk, resci, deltak, deltac
real(kind=wp) :: ralpha, rAlphaC, rAlphaK, rbeta
real(kind=wp), external :: ddot_

select case (PCGS)
  case (4)
    !! precondition p (K-1*p)
    call ipIn(ipS2)
    call DMinvCI_SA(CGSvec%ipPvec,W(ipCId)%A,Fancy)
    call opOut(ipCI)
    call opOut(ipDia)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,CGSvec%Pvec,dKappa,Sc1,iSym,iter)
    call opOut(ipPre2)

    if (.not. InvSCF) call DMInvISR(ISR%Pvec,ISR%p)

    !! A*K-1*p
    call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)
  case (5)
    if (.not. InvSCF) ISR%p(:,:) = ISR%Pvec(:,:)
    !! A*p
    call TimesE2(CGSvec%Pvec,CGSvec%ipPvec,1,reco,jspin,ipS2,dKappa,ipCId)

    !! precondition A*p (K-1*A*p)
    call ipIn(ipS2)
    call DMinvCI_SA(ipCId,W(ipS1)%A,Fancy)
    call opOut(ipCI)
    call opOut(ipDia)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,dKappa,Temp4,Sc1,iSym,iter)
    call opOut(ipPre2)

    if (.not. InvSCF) then
      ISR%p(:,:) = ISR%Ap(:,:)
      call DMInvISR(ISR%p,ISR%Ap)
    end if
end select

!! compute alpha
rAlphaK = ddot_(nDensC,CGSvec%R0,1,Temp4,1)
rAlphaC = ddot_(nConf1*nroots,W(CGSvec%ipR0)%A,1,W(ipS1)%A,1)
if (.not. InvSCF) rAlphaC = rAlphaC+ddot_(nRoots**2,ISR%R0,1,ISR%Ap,1)
rAlpha = delta/(rAlphaK+rAlphaC)

! Algorithm 4: q_k = u_k - alpha*(A*K-1*p)
! Algorithm 5: q_k = u_k - alpha*(K-1*A*p)
CGSvec%Qvec(1:nDensC) = CGSvec%Uvec(1:nDensC)-ralpha*Temp4(1:nDenSC)
W(CGSvec%ipQvec)%A(1:nConf1*nRoots) = W(CGSvec%ipUvec)%A(1:nConf1*nRoots)-ralpha*W(ipS1)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%Qvec(:,:) = ISR%Uvec(:,:)-ralpha*ISR%Ap(:,:)

!! construct u = u + q
CGSvec%Uvec(1:nDensC) = CGSvec%Uvec(1:nDensC)+CGSvec%Qvec(1:nDensC)
W(CGSvec%ipUvec)%A(1:nConf1*nRoots) = W(CGSvec%ipUvec)%A(1:nConf1*nRoots)+W(CGSvec%ipQvec)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%Uvec(:,:) = ISR%Uvec(:,:)+ISR%Qvec(:,:)

select case (PCGS)
  case (4)
    !! precondition u + q
    call ipIn(ipS2)
    call DMinvCI_SA(CGSvec%ipUvec,W(ipCId)%A,Fancy)
    call opOut(ipCI)
    call opOut(ipDia)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,CGSvec%Uvec,dKappa,Sc1,iSym,iter)
    call opOut(ipPre2)

    if (.not. InvSCF) call DMInvISR(ISR%Uvec,ISR%p)

    !! x_k+1 = x_k + alpha*K-1*(u+q)
    Kappa(1:nDensC) = Kappa(1:nDensC)+ralpha*dKappa(1:nDensC)
    W(ipCIT)%A(1:nConf1*nRoots) = W(ipCIT)%A(1:nConf1*nRoots)+ralpha*W(ipCId)%A(1:nConf1*nRoots)
    if (.not. InvSCF) ISR%Xvec(:,:) = ISR%Xvec(:,:)+ralpha*ISR%p(:,:)

    !! A*K-1*(u+q)
    call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)
  case (5)
    !! x_k+1 = x_k + alpha*(u+q)
    Kappa(1:nDensC) = Kappa(1:nDensC)+ralpha*CGSvec%Uvec(1:nDensC)
    W(ipCIT)%A(1:nConf1*nRoots) = W(ipCIT)%A(1:nConf1*nRoots)+ralpha*W(CGSvec%ipUvec)%A(1:nConf1*nRoots)
    if (.not. InvSCF) ISR%Xvec(:,:) = ISR%Xvec(:,:)+ralpha*ISR%Uvec(:,:)

    !! A*(u+q)
    if (.not. InvSCF) ISR%p(:,:) = ISR%Uvec(:,:)
    call TimesE2(CGSvec%Uvec,CGSvec%ipUvec,1,reco,jspin,ipS2,dKappa,ipCId)

    !! precondition A*(u+q)
    !! rather than computing K-1*r_k+1, accumulate -K^1*A*(u+q)
    call ipIn(ipS2)
    call DMinvCI_SA(ipCId,W(ipS1)%A,Fancy)
    call opOut(ipCI)
    call opOut(ipDia)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,dKappa,Temp4,Sc1,iSym,iter)
    call opOut(ipPre2)

    if (.not. InvSCF) then
      ISR%p(:,:) = ISR%Ap(:,:)
      call DMInvISR(ISR%p,ISR%Ap)
    end if
end select

!! r_k+1 = r_k - alpha*A*(u+q)
Sigma(1:nDensC) = Sigma(1:nDensC)-ralpha*Temp4(1:nDensC)
resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
W(ipST)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)-ralpha*W(ipS1)%A(1:nConf1*nRoots)
resci = sqrt(ddot_(nconf1*nroots,W(ipST)%A,1,W(ipST)%A,1))
if (.not. InvSCF) then
  ISR%Rvec(:,:) = ISR%Rvec(:,:)-ralpha*ISR%Ap(:,:)
  resci = resci+sqrt(ddot_(nRoots**2,ISR%Rvec,1,ISR%Rvec,1))
end if

!! compute beta
deltaK = ddot_(nDensC,CGSvec%R0,1,Sigma,1)
deltaC = ddot_(nConf1*nroots,W(CGSvec%ipR0)%A,1,W(ipST)%A,1)
if (.not. InvSCF) deltaC = deltaC+ddot_(nRoots**2,ISR%R0,1,ISR%Rvec,1)
rbeta = (deltaC+deltaK)/delta
delta = deltac+deltaK

!! u_k = r_k + beta*q
CGSvec%Uvec(1:nDensC) = Sigma(1:nDensC)+rbeta*CGSvec%Qvec(1:nDensC)
W(CGSvec%ipUvec)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)+rbeta*W(CGSvec%ipQvec)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%Uvec(:,:) = ISR%Rvec(:,:)+rbeta*ISR%Qvec(:,:)

!! p = u + beta*q + beta*beta*p
CGSvec%Pvec(1:nDensC) = CGSvec%Uvec(1:nDensC)+rbeta*CGSvec%Qvec(1:nDensC)+(rbeta**2)*CGSvec%Pvec(1:nDensC)
W(CGSvec%ipPvec)%A(1:nConf1*nRoots) = W(CGSvec%ipUvec)%A(1:nConf1*nRoots)+rbeta*W(CGSvec%ipQvec)%A(1:nConf1*nRoots)+ &
                                      (rbeta**2)*W(CGSvec%ipPvec)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%Pvec(:,:) = ISR%Uvec(:,:)+rbeta*ISR%Qvec(:,:)+rbeta*rbeta*ISR%Pvec(:,:)

end subroutine CGS_x
