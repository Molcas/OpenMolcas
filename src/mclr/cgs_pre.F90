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
! Prepare conjugate gradient squared-related things

subroutine CGS_pre(nDensC,nConf1,nRoots,Kappa,dKappa,Sigma,Fancy,Sc1,ipCI,ipCId,ipDia,ipPre2,iPre,ipST,ipS1,delta,iSym,iter)

use ipPage, only: ipIn, opOut, W
use ISRotation, only: DMInvISR, InvSCF, ISR
use cgs_mod, only: CGSvec, PCGS
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: nDensC, nConf1, nRoots, ipCI, ipCId, ipDia, ipPre2, iPre(*), ipST, ipS1
real(kind=wp), intent(inout) :: Kappa(*), dKappa(*), Sigma(*), Fancy(*), Sc1(*)
real(kind=wp), intent(out) :: delta
integer(kind=iwp), intent(inout) :: iSym, iter
real(kind=wp) :: deltaC, deltaK
real(kind=wp), external :: ddot_

! R0
CGSvec%R0(1:nDensC) = Kappa(1:nDensC)
W(CGSvec%ipR0)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%R0(:,:) = ISR%Rvec(:,:)

if (PCGS == 5) then
  !! precondition p (K-1*p)
  call ipIn(ipST)
  call DMinvCI_SA(CGSvec%ipR0,W(ipST)%A,Fancy)
  call opOut(ipCI)
  call opOut(ipdia)

  call ipIn(ipPre2)
  call DMInvKap(W(ipPre2)%A,iPre,CGSvec%R0,Kappa,Sc1,iSym,iter)
  call opOut(ipPre2)

  if (.not. InvSCF) call DMInvISR(ISR%R0,ISR%Rvec)

  !! r0 <-- K-1*r
  CGSvec%R0(1:nDensC) = Kappa(1:nDensC)
  W(CGSvec%ipR0)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)
  if (.not. InvSCF) ISR%R0(:,:) = ISR%Rvec(:,:)
end if

! u_k
CGSvec%Uvec(1:nDensC) = Kappa(1:nDensC)
W(CGSvec%ipUvec)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)
if (.not. InvSCF) ISR%Uvec(:,:) = ISR%Rvec(:,:)

! p_k
dKappa(1:nDensC) = CGSvec%Uvec(1:nDensC)
CGSvec%Pvec(1:nDensC) = dKappa(1:nDensC)
W(ipCId)%A(1:nConf1*nRoots) = W(CGSvec%ipUvec)%A(1:nConf1*nRoots)
W(CGSvec%ipPvec)%A(1:nConf1*nRoots) = W(ipCId)%A(1:nConf1*nRoots)
if (.not. InvSCF) then
  ISR%Pvec(:,:) = ISR%Uvec(:,:)
  if (PCGS == 4) ISR%Rvec(:,:) = ISR%R0(:,:)
  ISR%Xvec(:,:) = Zero
end if

deltaC = ddot_(nConf1*nroots,W(CGSvec%ipR0)%A,1,W(CGSvec%ipR0)%A,1)
if (.not. InvSCF) deltaC = deltaC+ddot_(nRoots**2,ISR%R0,1,ISR%R0,1)
deltaK = ddot_(nDensC,CGSvec%R0,1,CGSvec%R0,1)
Kappa(1:nDensC) = Zero
delta = deltaC+deltaK

if (PCGS == 5) then
  Sigma(1:nDensC) = CGSvec%R0(1:nDensC)
  W(ipS1)%A(1:nConf1*nRoots) = W(CGSvec%ipR0)%A(1:nConf1*nRoots)
end if

end subroutine CGS_pre
