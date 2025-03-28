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

use Exp, only: Exp_Close
use Arrays, only: SFock, G1m, G2mp, Int2, FIMO
use ipPage, only: W
use MCLR_Data, only: nConf1, nDens2, nNA, nDensC, nDens, ipCI, n1Dens
use MCLR_Data, only: RMS, rAlpha
use MCLR_Data, only: ipDia
use MCLR_Data, only: LuTemp
use MCLR_Data, only: XISPSM
use MCLR_Data, only: MS2P
use input_mclr, only: nDisp, Fail, State_Sym, iMethod, rIn_Ene, PotNuc, iBreak, Eps, nIter, Debug, ERASSCF, kPrint, nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: u5, u6

implicit none
integer iKapDisp(nDisp), isigDisp(nDisp)
integer iCIDisp(nDisp), iCIsigDisp(nDisp)
integer iRHSDisp(nDisp), iRHSCIDisp(nDisp)
character(len=8) Fmt2
logical lPrint, cnvrgd
real*8 rdum(1)
real*8 d_0
real*8, allocatable :: Kappa(:), dKappa(:), Sigma(:), Temp1(:), Temp2(:), Temp3(:), Temp4(:), Sc1(:), Sc2(:), Sc3(:), Dens(:), &
                       Pens(:), rmoaa(:), rmoaa2(:), Pre2(:)
integer lPaper, lLine, Left, iDis, nConf3, ipS1, ipS2, ipST, ipCIT, ipCID, iDisp, iLen, Iter, i1, j1
real*8 DeltaC, DeltaK, Delta, Delta0, rGrad, Ec, rAlphaC, rAlphaK, ResK, ResCI, rBeta, Res, rCHC
real*8, external :: DDot_
integer, external :: ipGet
!----------------------------------------------------------------------*
interface
  subroutine CISigma(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    integer iispin, iCsym, iSSym
    integer nInt1, nInt2s, nInt2a
    real*8, target :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    integer ipCI1, ipCI2
    logical Have_2_el
  end subroutine CISigma
  subroutine FockGen_sp(d_0,rDens1,rdens2,Fock,fockout,idsym)
    real*8 d_0
    real*8 rDens1(*), rdens2(*), Fock(*), fockout(*)
    integer idsym
  end subroutine FockGen_sp
end interface

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

if (iand(kprint,2) == 2) lprint = .true.
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

end if

idisp = 1

! Allocate areas for scratch and state variables

call mma_allocate(Kappa,nDens2+6,Label='Kappa')
call mma_allocate(SFock,nDens2+6,Label='SFock')
call mma_allocate(dKappa,nDens2+6,Label='dKappa')
call mma_allocate(Sigma,nDens2+6,Label='Sigma')
call mma_allocate(Temp1,nDens2+6,Label='Temp1')
call mma_allocate(Temp2,nDens2+6,Label='Temp2')
call mma_allocate(Temp3,nDens2+6,Label='Temp3')
call mma_allocate(Temp4,nDens2+6,Label='Temp4')
call mma_allocate(Sc1,nDens2+6,Label='Sc1')
call mma_allocate(Sc2,nDens2+6,Label='Sc2')
call mma_allocate(Sc3,nDens2+6,Label='Sc3')
call mma_allocate(Pre2,nDensC,Label='Pre2')
Temp1(1:nDens2) = Zero
Kappa(1:nDens2) = Zero
dKappa(1:nDens2) = Zero
Sigma(1:nDens2) = Zero
if (iMethod == 2) then
  call mma_allocate(Dens,n1dens,Label='Dens')
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
call dcopy_(nconf1,[Zero],0,W(ipST)%Vec,1)

if (lprint) write(u6,*) '       Iteration         Delta     Res(kappa) Res(CI)'
iLen = nDensC
iRHSDisp(iDisp) = iDis
call Compress(Temp4,Sigma,1)
call DSCAL_(ndensc,-sqrt(OneHalf)*dble(ms2p),Sigma,1)
call UnCompress(Sigma,Temp4,1)
call dDaFile(LuTemp,1,Sigma,iLen,iDis)
if (iMethod == 2) then
  call ipin(ipCIT)
  call dcopy_(nConf1,[Zero],0,W(ipCIT)%Vec,1)
end if
call ipout(ipcit)
call ipin(ipST)
if (iMethod == 2) then
  ilen = nconf1
  iRHSCIDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
end if
call DSCAL_(nConf1,-One,W(ipST)%Vec,1)
call DSCAL_(nDensC,-One,Sigma,1)

call DMInvKap_sp(Sigma,dKappa,1)

call ipin(ipCId)
if (nconf1 > 1) then
  call DMinvCI(ipST,W(ipCId)%Vec,rCHC,1)
else
  call ipin(ipST)
  call dcopy_(nconf1,W(ipST)%Vec,1,W(ipCid)%Vec,1)
end if

if ((iMethod == 2) .and. (nconf1 /= 0)) then
  call ipin(ipST)
  call ipin(ipCId)
  deltaC = ddot_(nConf1,W(ipST)%Vec,1,W(ipCId)%Vec,1)
  call ipout(ipcid)
else
  deltac = Zero
end if
deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
Kappa(1:nDens) = Zero
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
    dKappa(1:nDens2) = Zero
    dKappa(i1) = One
  else
    call ipin(ipCID)
    W(ipCID)%Vec(1:nConf1) = Zero
    W(ipCID)%Vec(i1) = One
    call ipout(ipcid)
  end if
  !************************************************************

  call RInt_SP(dKappa,rmoaa,rmoaa2,Temp4,Sc2)

  if ((i1 > 0) .and. (j1 > 0)) write(u6,*) 'Kap_sig',Sc2(j1)
  if (nconf1 > 1) then
    call opout(-1)
    call CISigma(1,State_Sym,state_sym,Temp4,nDens2,rmoaa,size(rmoaa),rmoaa2,size(rmoaa2),ipCI,ipS1,.true.)
    call opout(-1)
    call ipin(ipCI)
    call ipin(ipS1)
    rGrad = ddot_(nconf1,W(ipCI)%Vec,1,W(ipS1)%Vec,1)
    call daxpy_(nConf1,-rgrad,W(ipCI)%Vec,1,W(ipS1)%Vec,1)
    call dscal_(nconf1,-rms*sqrt(OneHalf)*Two,W(ipS1)%Vec,1)

    if ((i1 > 0) .and. (j1 < 0)) write(u6,*) 'CI_sig',W(ipS1)%Vec(j1)

    call opout(-1)
    if (nconf1 > 1) then
      call CISigma(0,State_Sym,state_sym,FIMO,size(FIMO),Int2,size(Int2),rdum,1,ipCId,ipS2,.true.)
      call opout(-1)
      EC = rin_ene+potnuc-ERASSCF(1)

      call ipin(ipCId)
      call ipin(ipS2)
      call DaXpY_(nConf1,EC,W(ipCId)%Vec,1,W(ipS2)%Vec,1)
      call DSCAL_(nConf1,Two,W(ipS2)%Vec,1)
      if ((i1 < 0) .and. (j1 < 0)) write(u6,*) 'CI_sig',W(ipS2)%Vec(j1)

      call ipin(ipCI)
      call ipin(ipCid)
      call SpinDens(W(ipCI)%Vec,W(ipCid)%Vec,State_Sym,State_sym,Pens,rdum,rdum,rdum,rdum,Dens,rdum,1)

      d_0 = ddot_(nconf1,W(ipCid)%Vec,1,W(ipci)%Vec,1)
      call FockGen_sp(d_0,Dens,Pens,Sc3,Sc1,1)
      call DSCAL_(ndens2,-rms*sqrt(OneHalf),Sc1,1)

      call Compress(Sc1,Sc3,1)
      if ((i1 < 0) .and. (j1 > 0)) write(u6,*) 'CI_sig',Sc3(j1)
      cycle
    end if

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

  if (nconf1 > 1) then
    call DZaXpY(nDens,One,Sc2,1,Sc3,1,Temp4,1)
  else
    call dcopy_(nDens,Sc2,1,Temp4,1)
  end if
  call dcopy_(nDens,dKappa,1,Temp2,1)
  if (nconf1 > 1) then
    call ipin1(ipS1,nconf1)
    call ipin1(ipS2,nconf1)
    call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipS1)%Vec,1)
  else
    call ipin1(ipS1,nconf1)
    W(ipS1)%Vec(1:nconf1) = Zero
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
  rAlphaK = Zero
  rAlphaK = ddot_(nDensC,Temp4,1,Temp2,1)
  if (nconf1 /= 0) then
    call ipin(ipS1)
    call ipin(ipCId)
    rAlphaC = ddot_(nConf1,W(ipS1)%Vec,1,W(ipCId)%Vec,1)
  end if
  rAlpha = delta/(rAlphaK+ralphaC)

  !--------------------------------------------------------------------*

  ! Kappa=Kappa+rAlpha*dKappa
  ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa

  call DaxPy_(nDensC,ralpha,Temp2,1,Kappa,1)
  call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
  resk = sqrt(ddot_(nDensC,Temp4,1,Temp4,1))
  resci = Zero
  if (nconf1 /= 0) then
    call ipin(ipCId)
    call ipin(ipCIT)
    call DaXpY_(nConf1,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
    call ipout(ipcit)
    call ipin1(ipST,nconf1)
    call ipin(ipS1)
    call DaXpY_(nConf1,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
    call opout(ipS1)
    call ipin(ipST)
    resci = sqrt(ddot_(nconf1,W(ipST)%Vec,1,W(ipST)%Vec,1))
  end if

  ! Precondition......
  !    -1
  ! S=M  Sigma

  call opout(ipcid)
  call ipin(ipS2)
  if (nconf1 > 1) then
    call DMinvCI(ipST,W(ipS2)%Vec,rCHC,1)
  else
    call ipin(ipST)
    call dcopy_(nconf1,W(ipST)%Vec,1,W(ipS2)%Vec,1)
  end if

  call opout(ipci)
  call opout(ipdia)

  call DMInvKap_sp(Sigma,Sc2,1)

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
    deltaC = ddot_(nConf1,W(ipST)%Vec,1,W(ipS2)%Vec,1)
    call ipout(ipST)
  else
    deltaC = Zero
  end if

  deltaK = ddot_(nDensC,Sigma,1,Sc2,1)
  if (imethod /= 2) then
    rBeta = deltaK/delta
    delta = deltaK
    call DScal_(nDensC,rBeta,Temp2,1)
    call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
  else
    rbeta = (deltac+deltaK)/delta
    delta = deltac+deltaK
    call ipin(ipCID)
    call DScal_(nConf1,rBeta,W(ipCID)%Vec,1)
    call DScal_(nDensC,rBeta,Temp2,1)
    call ipin(ipS2)
    call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
    call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
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

  call dcopy_(ndensc,Temp2,1,dKappa,1)

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
iLen = ndensC
iKapDisp(iDisp) = iDis
call dDaFile(LuTemp,1,Kappa,iLen,iDis)
iSigDisp(iDisp) = iDis
call dDaFile(LuTemp,1,Sigma,iLen,iDis)
if (iMethod == 2) then
  ilen = nconf1
  iCIDisp(iDisp) = iDis
  call ipin(ipCIT)
  call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)
  iCISigDisp(iDisp) = iDis
  call ipin(ipST)
  call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
end if

call mma_deallocate(Temp4)
call mma_deallocate(Temp3)
call mma_deallocate(Temp2)
call mma_deallocate(Temp1)
call mma_deallocate(dKappa)
call mma_deallocate(Sigma)
call mma_deallocate(Kappa)
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
