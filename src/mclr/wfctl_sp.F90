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

subroutine WfCtl_sp(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,iRHSCIDISP)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!***********************************************************************

use ipPage, only: ipclose, ipget, ipin, ipin1, ipnout, ipout, opout, W
use MCLR_Data, only: FIMO, G1m, G2mp, Int2, ipCI, ipDia, LuTemp, MS2P, n1Dens, nConf1, nDens, nDensC, nNA, RMS, SFock, XISPSM
use MCLR_procedures, only: CISigma
use input_mclr, only: Debug, Eps, ERASSCF, Fail, iBreak, iMethod, kPrint, nCSF, nDisp, nIter, PotNuc, rIn_Ene, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp), &
                                  iRHSCIDisp(nDisp)
integer(kind=iwp) :: i1, iDis, iDisp, iLen, ipCID, ipCIT, ipS1, ipS2, ipST, Iter, j1, Left, lLine, lPaper, nConf3
real(kind=wp) :: d_0, Delta, Delta0, DeltaC, DeltaK, Ec, rAlpha, rAlphaC, rAlphaK, rBeta, rCHC, rdum(1), Res, ResCI, ResK, rGrad
logical(kind=iwp) :: cnvrgd, lPrint
character(len=8) :: Fmt2
real(kind=wp), allocatable :: Dens(:), dKappa(:), Kappa(:), Pens(:), rmoaa(:), rmoaa2(:), Sc1(:), Sc2(:), Sc3(:), Sc4(:), &
                              Sigma(:), Temp1(:), Temp2(:), Temp4(:)
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

lPaper = 132
lLine = 120
left = (lPaper-lLine)/2
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

fail = .false.
idis = 0
lprint = .false.
nconf1 = 0

if (btest(kprint,1)) lprint = .true.
if (iMethod == 2) call InCSFSD(State_Sym,State_sym)

nconf1 = ncsf(State_Sym)
nconf3 = nint(xispsm(State_SYM,1))
call Setup_MCLR(1)

!                            [2]
! Calculate the diagonal of E    and store in core/disc

if (imethod > 0) then
  if (nconf1 > 1) call CIDia_MCLR(State_Sym,rCHC)
  call ipout(ipdia)

  ! Allocate disk/memory space

  ! This areas should be addressed through ipin
  ! ipout will page them out to disk and free the memory area
  ! if page mode is used

  ! opout will release the memory area without update the disk

  ips1 = ipget(nconf3)
  ips2 = ipget(nconf3)
  ipst = ipget(nconf3)
  ipcit = ipget(nconf1)
  ipcid = ipget(nconf1)

else
  call Untested('WfCtl_sp')
  ! These are uninitialized!
  ips1 = -1
  ips2 = -1
  ipst = -1
  ipcid = -1
end if

idisp = 1

! Allocate areas for scratch and state variables

call mma_allocate(Kappa,nDensC,Label='Kappa')
call mma_allocate(SFock,nDens,Label='SFock')
call mma_allocate(dKappa,nDensC,Label='dKappa')
call mma_allocate(Sigma,nDensC,Label='Sigma')
call mma_allocate(Temp1,nDensC,Label='Temp1')
call mma_allocate(Temp2,nDensC,Label='Temp2')
call mma_allocate(Temp4,nDens,Label='Temp4')
call mma_allocate(Sc1,nDens,Label='Sc1')
call mma_allocate(Sc2,nDens,Label='Sc2')
call mma_allocate(Sc3,nDens,Label='Sc3')
call mma_allocate(Sc4,nDensC,Label='Sc4')
Kappa(:) = Zero
dKappa(:) = Zero
Sigma(:) = Zero
if (iMethod == 2) then
  call mma_allocate(Dens,n1Dens,Label='Dens')
  call mma_allocate(Pens,nna**4,Label='Pens')
  call mma_allocate(rmoaa,nna**4,Label='rmoaa')
  call mma_allocate(rmoaa2,nna**4,Label='rmoaa2')
end if

!-----------------------------------------------------------------------
!
! Calculate RHS
!
!-----------------------------------------------------------------------

!call Pre_SP(Pre2,1)
call FockGen_sp(Zero,G1m,G2mp,SFock,Temp4,1)
call ipin(ipST)
W(ipST)%A(1:nconf1) = Zero

if (lprint) write(u6,*) '       Iteration         Delta     Res(kappa) Res(CI)'
iLen = nDensC
iRHSDisp(iDisp) = iDis
call Compress(Temp4,Sigma,1)
Sigma(:) = -sqrt(OneHalf)*dble(ms2p)*Sigma(:)
call UnCompress(Sigma,Temp4,1)
call dDaFile(LuTemp,1,Sigma,iLen,iDis)
if (iMethod == 2) then
  call ipin(ipCIT)
  W(ipCIT)%A(1:nConf1) = Zero
end if
call ipout(ipcit)
call ipin(ipST)
if (iMethod == 2) then
  ilen = nconf1
  iRHSCIDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
end if
W(ipST)%A(1:nConf1) = -W(ipST)%A(1:nConf1)
Sigma(:) = -Sigma(:)

call DMInvKap_sp(Sigma,dKappa)

call ipin(ipCId)
if (nconf1 > 1) then
  call DMinvCI(ipST,W(ipCId)%A,rCHC,1)
else
  call ipin(ipST)
  W(ipCid)%A(1:nConf1) = W(ipST)%A(1:nConf1)
end if

if ((iMethod == 2) .and. (nconf1 /= 0)) then
  call ipin(ipST)
  call ipin(ipCId)
  deltaC = ddot_(nConf1,W(ipST)%A,1,W(ipCId)%A,1)
  call ipout(ipcid)
else
  deltac = Zero
end if
deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
Kappa(:) = Zero
delta = deltac+deltaK
delta0 = delta
iter = 1
!-----------------------------------------------------------------------

!**********************************************************
!          I   T   E   R   A   T   I   O   N   S          *
!**********************************************************

cnvrgd = .true.
do
  !if (delta == Zero) exit
  read(u5,*) i1,j1
  if (i1 > 0) then
    dKappa(:) = Zero
    dKappa(i1) = One
  else
    call ipin(ipCID)
    W(ipCID)%A(1:nConf1) = Zero
    W(ipCID)%A(i1) = One
    call ipout(ipcid)
  end if
  !************************************************************

  call RInt_SP(dKappa,rmoaa,rmoaa2,Temp4,Sc4)

  if ((i1 > 0) .and. (j1 > 0)) write(u6,*) 'Kap_sig',Sc4(j1)
  if (nconf1 > 1) then
    call opout(-1)
    call CISigma(1,State_Sym,state_sym,Temp4,nDens,rmoaa,nna**4,rmoaa2,nna**4,ipCI,ipS1,.true.)
    call opout(-1)
    call ipin(ipCI)
    call ipin(ipS1)
    rGrad = ddot_(nconf1,W(ipCI)%A,1,W(ipS1)%A,1)
    W(ipS1)%A(1:nConf1) = rms*sqrt(OneHalf)*Two*(rgrad*W(ipCI)%A(1:nConf1)-W(ipS1)%A(1:nConf1))

    if ((i1 > 0) .and. (j1 < 0)) write(u6,*) 'CI_sig',W(ipS1)%A(j1)

    call opout(-1)

    call CISigma(0,State_Sym,state_sym,FIMO,size(FIMO),Int2,size(Int2),rdum,1,ipCId,ipS2,.true.)
    call opout(-1)
    EC = rin_ene+potnuc-ERASSCF(1)

    call ipin(ipCId)
    call ipin(ipS2)
    W(ipS2)%A(1:nConf1) = Two*(W(ipS2)%A(1:nConf1)+EC*W(ipCId)%A(1:nConf1))
    if ((i1 < 0) .and. (j1 < 0)) write(u6,*) 'CI_sig',W(ipS2)%A(j1)

    call ipin(ipCI)
    call ipin(ipCid)
    call SpinDens(W(ipCI)%A,W(ipCid)%A,State_Sym,State_sym,Pens,rdum,rdum,rdum,rdum,Dens,rdum,1)

    d_0 = ddot_(nconf1,W(ipCid)%A,1,W(ipci)%A,1)
    call FockGen_sp(d_0,Dens,Pens,Sc3,Sc1,1)
    Sc2(:) = -rms*sqrt(OneHalf)*Sc1(:)

    call Compress(Sc1,Sc4,1)
    if ((i1 < 0) .and. (j1 > 0)) write(u6,*) 'CI_sig',Sc4(j1)
    cycle

  end if

  !*********************************************************************
  !
  ! Sc1  kappa-> kappa
  ! Sc3  CI -> kappa
  ! S1   kappa -> CI
  ! S2   CI -> CI
  ! dKap present step
  ! Kap  kappaX
  ! CIT  CIX
  ! CId  present step
  !
  ! Add together
  !
  !*********************************************************************

  call Untested('WfCtl_sp')
  ! IFG: Given the sizes assumed later, I think this assignment is wrong
  !      and it should use the Compress call instead
  !Temp4(:) = Sc2(:)
  call Compress(Sc2,Temp1,1)
  Temp2(:) = dKappa(:)
  if (nconf1 > 1) then
    call ipin1(ipS1,nconf1)
    call ipin1(ipS2,nconf1)
    W(ipS1)%A(1:nConf1) = W(ipS1)%A(1:nConf1)+W(ipS2)%A(1:nConf1)
  else
    call ipin1(ipS1,nconf1)
    W(ipS1)%A(1:nconf1) = Zero
  end if

  !---------------------------------------------------------------------
  !
  !                ######   #####   #####
  !                #     # #     # #     #
  !                #     # #       #
  !                ######  #       #  ####
  !                #       #       #     #
  !                #       #     # #     #
  !                #        #####   #####
  !
  !---------------------------------------------------------------------
  !*********************************************************************
  !
  !
  !            delta
  ! rAlpha=------------
  !        dKappa:dSigma
  !
  !---------------------------------------------------------------------
  rAlphaC = Zero
  rAlphaK = ddot_(nDensC,Temp1,1,Temp2,1)
  if (nconf1 /= 0) then
    call ipin(ipS1)
    call ipin(ipCId)
    rAlphaC = ddot_(nConf1,W(ipS1)%A,1,W(ipCId)%A,1)
  end if
  rAlpha = delta/(rAlphaK+ralphaC)

  !--------------------------------------------------------------------*

  ! Kappa=Kappa+rAlpha*dKappa
  ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa

  Kappa(:) = Kappa(:)+ralpha*Temp2(:)
  Sigma(:) = Sigma(:)-ralpha*Temp1(:)
  resk = sqrt(ddot_(nDensC,Temp1,1,Temp1,1))
  resci = Zero
  if (nconf1 /= 0) then
    call ipin(ipCId)
    call ipin(ipCIT)
    W(ipCIT)%A(1:nConf1) = W(ipCIT)%A(1:nConf1)+ralpha*W(ipCId)%A(1:nConf1)
    call ipout(ipcit)
    call ipin1(ipST,nconf1)
    call ipin(ipS1)
    W(ipST)%A(1:nConf1) = W(ipST)%A(1:nConf1)-ralpha*W(ipS1)%A(1:nConf1)
    call opout(ipS1)
    call ipin(ipST)
    resci = sqrt(ddot_(nconf1,W(ipST)%A,1,W(ipST)%A,1))
  end if

  ! Precondition......
  !    -1
  ! S=M  Sigma

  call opout(ipcid)
  call ipin(ipS2)
  if (nconf1 > 1) then
    call DMinvCI(ipST,W(ipS2)%A,rCHC,1)
  else
    call ipin(ipST)
    W(ipS2)%A(1:nconf1) = W(ipST)%A(1:nconf1)
  end if

  call opout(ipci)
  call opout(ipdia)

  call DMInvKap_sp(Sigma,Sc4)

  !--------------------------------------------------------------------*
  !      s:Sigma
  ! Beta=-------
  !       delta
  !
  ! delta=s:sigma
  !
  ! dKappa=s+Beta*dKappa

  if ((iMethod == 2) .and. (nconf1 /= 0)) then
    call ipin(ipST)
    call ipin(ipS2)
    deltaC = ddot_(nConf1,W(ipST)%A,1,W(ipS2)%A,1)
    call ipout(ipST)
  else
    deltaC = Zero
  end if

  deltaK = ddot_(nDensC,Sigma,1,Sc4,1)
  if (imethod /= 2) then
    rBeta = deltaK/delta
    delta = deltaK
    Temp2(:) = rBeta*Temp2(:)+Sc4(:)
  else
    rbeta = (deltac+deltaK)/delta
    delta = deltac+deltaK
    call ipin(ipCID)
    call ipin(ipS2)
    W(ipCID)%A(1:nConf1) = rBeta*W(ipCID)%A(1:nConf1)+W(ipS2)%A(1:nConf1)
    Temp2(:) = rBeta*Temp2(:)+Sc4(:)
    call opout(ipS2)
    call ipout(ipCID)
  end if

  !    ######  #    #  #####        #####    ####    ####
  !    #       ##   #  #    #       #    #  #    #  #    #
  !    #####   # #  #  #    #       #    #  #       #
  !    #       #  # #  #    #       #####   #       #  ###
  !    #       #   ##  #    #       #       #    #  #    #
  !    ######  #    #  #####        #        ####    ####
  !
  !--------------------------------------------------------------------*

  dKappa(:) = Temp2(:)

  res = Zero ! dummy initialize
  if (iBreak == 1) then
    if (abs(delta) < abs(Eps**2*delta0)) exit
  else if (ibreak == 2) then
    res = sqrt(resk**2+resci**2)
    if (res < abs(Eps)) exit
  else
    if ((abs(delta) < abs(Eps**2*delta0)) .and. (res < abs(Eps))) exit
  end if
  if (iter >= niter) then
    cnvrgd = .false.
    exit
  end if
  if (lprint) write(u6,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)') '     ',iter,'       ',delta/delta0,resk,resci,deltac,deltak

  iter = iter+1

end do

!***********************************************************************

if (.not. cnvrgd) then
  write(u6,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
  fail = .true.
else
  write(u6,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
  call ipnout(-1)
end if
write(u6,*)
iLen = nDensC
iKapDisp(iDisp) = iDis
call dDaFile(LuTemp,1,Kappa,iLen,iDis)
iSigDisp(iDisp) = iDis
call dDaFile(LuTemp,1,Sigma,iLen,iDis)
if (iMethod == 2) then
  ilen = nconf1
  iCIDisp(iDisp) = iDis
  call ipin(ipCIT)
  call dDaFile(LuTemp,1,W(ipCIT)%A,iLen,iDis)
  iCISigDisp(iDisp) = iDis
  call ipin(ipST)
  call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
end if

call mma_deallocate(Temp4)
call mma_deallocate(Temp2)
call mma_deallocate(Temp1)
call mma_deallocate(dKappa)
call mma_deallocate(Sigma)
call mma_deallocate(Kappa)
call mma_deallocate(Sc4)
call mma_deallocate(Sc3)
call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
if (iMethod == 2) then
  call mma_deallocate(rmoaa2)
  call mma_deallocate(rmoaa)
  call mma_deallocate(Pens)
  call mma_deallocate(Dens)
end if

! Free all memory and remove from disk all data
! related to this symmetry

if (imethod == 2) call ipclose(ipci)

call Exp_Close()

if (debug) then
  write(u6,*) '********************************************************************************'
  write(u6,*)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine WfCtl_sp
