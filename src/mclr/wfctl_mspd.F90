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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

subroutine WfCtl_MSPD(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,converged,iPL)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!***********************************************************************

use Exp, only: Exp_Close
use ipPage, only: W
use cmslag, only: ResQaaLag2
use MCLR_Data, only: nConf1, nDens2, nDensC, nDens, ipCI
use MCLR_Data, only: ipDia
use MCLR_Data, only: ISNAC, OVERRIDE, IRLXROOT, ISMECIMSPD, NACSTATES
use MCLR_Data, only: LuTemp, LuQDat
use MCLR_Data, only: XISPSM
use input_mclr, only: nDisp, Fail, lSave, State_Sym, iMethod, iBreak, Eps, nIter, Debug, kPrint, nCSF, nRoots, TwoStep, StepType, &
                      iAddressQDat, nAsh, nRS2
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer iKapDisp(nDisp), isigDisp(nDisp)
integer iCIDisp(nDisp), iCIsigDisp(nDisp)
integer iRHSDisp(nDisp)
logical converged(8)
integer iPL
#include "rasdim.fh"
#include "warnings.h"
logical CI
character(len=8) Fmt2
integer opOut
logical lPrint, cnvrgd
real*8 rchc(mxroot)
real*8, allocatable :: Kappa(:), dKappa(:), Sigma(:), Temp3(:), Temp4(:), Sc1(:), Sc2(:), Fancy(:)
integer LURot, IsFreeUnit, JRoot
external IsFreeUnit
character(len=16) :: VecName
real*8 R1, R2, DeltaC, DeltaK, Delta, Delta0, ReCo, rAlphaC, rAlphaK, rAlpha, rEsk, rEsci, rBeta, Res
real*8, external :: DDot_
integer lPaper, lLine, Left, iDis, Lu_50, iDisp, iSym, nConf3, iRC, ipS1, ipS2, ipST, ipCIT, ipCID, nPre2, iLen, Iter, ipPre2, &
        jSpin, i
integer, external :: ipClose, ipGet, ipIn, ipOut, ipNOut
integer, external :: nPre

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
call StatusLine('MCLR: ','Computing Lagrangian multipliers for MS-PDFT')

lPaper = 132
lLine = 120
left = (lPaper-lLine)/2
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*

iDis = 0

fail = .false.
do i=1,8
  Converged(i) = .true.
end do
!MGD I think this is nice when printed...
lprint = .true.
reco = -One
Lu_50 = 50
if (lSAVE) then
  call DANAME(Lu_50,'RESIDUALS')
  write(u6,*) 'WfCtl_MSPD: SAVE option not implemented'
  call Abend()
end if
if (iand(kprint,2) == 2) lprint = .true.
isym = 1
nconf1 = ncsf(State_Sym)

CI = .false.
if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.

! Initiate CSF <-> SD
call InCSFSD(ieor(iSym-1,State_Sym-1)+1,State_sym,.false.)

! Calculate length of the density, Fock and Kappa matrix etc
! notice that this matrices are not necessarily symmetric.
! Store pointers.
!
! Input:
!        iSym: Symmetry of perturbation
!
! Output: Commonblocks (Pointers.fh)

nConf3 = nint(max(xispsm(State_SYM,1),xispsm(State_SYM,1)))

call Setup_MCLR(iSym)

! MSPDFT version ot lines 433-441 in mclr/rdinp_mclr.f
! MSPDFT version of reading NAC command in the mclr/rdinp_mclr.f
call Get_lScalar('isCMSNAC',isNAC)
if (isNAC) call Get_iArray('cmsNACstates',NACstates,2)
if (isNAC) call Get_lScalar('isMECIMSPD',isMECIMSPD)
if (isNAC) override = .true.

! Determine if we should page CI vectors
!                            [2]
! Calculate the diagonal of E    and store in core/disc

call mma_allocate(FANCY,nroots**3,Label='FANCY')
call CIDia_SA(State_Sym,rCHC,Fancy)

irc = ipOut(ipdia)

! Allocate disk/memory space

! This areas should be addressed through ipIn
! ipOut will page them out to disk and free the memory area
! if page mode is used

! opOut will release the memory area without update the disk

ipS1 = ipGet(nconf3*nroots)
ipS2 = ipGet(nconf3*nroots)
ipST = ipGet(nconf3*nroots)
ipCIT = ipGet(nconf1*nroots)
ipCId = ipGet(nconf1*nroots)

npre2 = npre(isym)
ipPre2 = ipGet(npre2)
irc = ipIn(ipPre2)
if (TwoStep .and. (StepType == 'RUN2')) then
  ! fetch data from LuQDAT and skip the call to "Prec"
  call ddafile(LuQDAT,2,W(ipPre2)%Vec,npre2,iaddressQDAT)
else
  call Prec(W(ipPre2)%Vec,isym)
  irc = ipOut(ippre2)
end if
if (TwoStep .and. (StepType == 'RUN1')) then
  ! save the computed data in "Prec" to LuQDAT and skip the
  ! following part of this function
  irc = ipIn(ipPre2)
  call ddafile(LuQDAT,1,W(ipPre2)%Vec,npre2,iaddressQDAT)
else

  ! OK START WORKING
  !
  !idisp = 1
  jspin = 0

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kappa,nDens2+6,Label='Kappa')
  call mma_allocate(dKappa,nDens2+6,Label='dKappa')
  call mma_allocate(Sigma,nDens2+6,Label='Sigma')
  call mma_allocate(Temp3,nDens2+6,Label='Temp3')
  call mma_allocate(Temp4,nDens2+6,Label='Temp4')
  call mma_allocate(Sc1,nDens2+6,Label='Sc1')
  call mma_allocate(Sc2,nDens2+6,Label='Sc2')

  cnvrgd = .true.
  do iDisp=1,nDisp
    Kappa(1:nDens2) = Zero
    dKappa(1:nDens2) = Zero
    Sigma(1:nDens2) = Zero

    !-------------------------------------------------------------------
    !
    ! Calculate RHS for the perturbation
    !
    !-------------------------------------------------------------------

    ! (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)

    if (debug) then
      if (isNAC) then
        write(u6,*) 'States: ',NACstates(1),NACstates(2)
      else
        write(u6,*) 'State: ',irlxroot
      end if
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (isNAC) then
      LURot = 233
      LURot = IsFreeUnit(LURot)
      call Molcas_Open(LURot,'ROT_VEC')
      do JRoot=1,nRoots
        read(LURot,*) VecName
      end do
      read(LURot,*) VecName
      close(LURot)
      if (VecName == 'CMS-PDFT') then
        call RHS_CMS_NAC(Temp4,W(ipST)%Vec)
        call DMinvCI_SA(ipST,W(ipS2)%Vec,Fancy)
      else
        write(u6,'(6X,A)') 'Error: Lagrangian Not Implemented for MS-PDFT'
        write(u6,'(6X,A)') '       Other Than CMS-PDFT'
        call xQuit(_RC_INPUT_ERROR_)
      end if
    else
      LURot = 233
      LURot = IsFreeUnit(LURot)
      call Molcas_Open(LURot,'ROT_VEC')
      do JRoot=1,nRoots
        read(LURot,*) VecName
      end do
      read(LURot,*) VecName
      close(LURot)
      if (VecName == 'CMS-PDFT') then
        call RHS_CMS(Temp4,W(ipST)%Vec)
        call DMinvCI_SA(ipST,W(ipS2)%Vec,Fancy)
      else
        write(u6,'(6X,A)') 'Error: Lagrangian Not Implemented for MS-PDFT'
        write(u6,'(6X,A)') '       Other Than CMS-PDFT'
        call xQuit(_RC_INPUT_ERROR_)
      end if
    end if

    irc = opOut(ipci)

    if (lprint) write(u6,*) '       Iteration       Delta           Res(kappa)       Res(CI)          DeltaK           DeltaC'
    iLen = nDensC
    iRHSDisp(iDisp) = iDis
    call Compress(Temp4,Sigma,iSym)
    r1 = ddot_(nDensc,Sigma,1,Sigma,1)
    if (debug) write(u6,*) 'Hi how about r1',r1
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)

    irc = ipIn(ipCIT)
    !irc = ipIn(ipST)
    irc = ipIn(ipCID)
    call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
    ! already initialized in rhs_mspdft
    !call dcopy_(nConf1*nroots,[Zero],0,W(ipST)%Vec,1)
    !call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
    !call dcopy_(nConf1*nroots,W(ipST)%Vec,1,W(ipCID)%Vec,1)
    irc = ipOut(ipCIT)
    call DSCAL_(nDensC,-One,Sigma,1)

    irc = ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,Kappa,nDens2+6,Temp3,nDens2+6,isym,iter)

    irc = opOut(ippre2)
    r2 = ddot_(ndensc,Kappa,1,Kappa,1)
    if (debug) write(u6,*) 'In that case I think that r2 should be:',r2
    if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

    call dcopy_(ndensC,Kappa,1,dKappa,1)

    ! In MS-PDFT deltaC is no longer zero initially
    !deltaC = Zero
    ! Use following two lines instead
    call DMinvCI_SA(ipST,W(ipS2)%Vec,Fancy)
    call dcopy_(nConf1*nroots,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
    deltaC = ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
    irc = ipOut(ipcid)
    deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
    Kappa(1:nDens) = Zero
    delta = deltac+deltaK
    delta0 = delta
    iter = 1
    !-------------------------------------------------------------------

    do
      if (delta == Zero) exit
      call TimesE2MSPDFT(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

      !-----------------------------------------------------------------
      !
      !            delta
      ! rAlpha=------------
      !        dKappa:dSigma
      !
      !-----------------------------------------------------------------

      rAlphaK = Zero
      rAlphaK = ddot_(nDensC,Temp4,1,dKappa,1)
      rAlphaC = Zero
      irc = ipIn(ipS1)
      irc = ipIn(ipCId)
      rAlphaC = ddot_(nConf1*nroots,W(ipS1)%Vec,1,W(ipCId)%Vec,1)
      rAlpha = delta/(rAlphaK+rAlphaC)

      !----------------------------------------------------------------*

      ! Kappa=Kappa+rAlpha*dKappa
      call DaxPy_(nDensC,ralpha,dKappa,1,Kappa,1)
      ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
      call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
      resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
      resci = Zero
      irc = ipIn(ipCIT)
      call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
      irc = ipOut(ipcit)
      ! ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
      irc = ipIn(ipS1)
      irc = ipIn(ipST)
      call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
      irc = opOut(ipS1)
      irc = ipIn(ipST)
      resci = sqrt(ddot_(nconf1*nroots,W(ipST)%Vec,1,W(ipST)%Vec,1))

      !----------------------------------------------------------------*

      ! Precondition......
      !    -1
      ! S=M  Sigma

      irc = opOut(ipcid)

      irc = ipIn(ipS2)
      call DMinvCI_SA(ipST,W(ipS2)%Vec,Fancy)
      irc = opOut(ipci)
      irc = opOut(ipdia)

      irc = ipIn(ipPre2)
      call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,Sc2,nDens2+6,Sc1,nDens2+6,iSym,iter)
      irc = opOut(ippre2)

      !----------------------------------------------------------------*
      !      s:Sigma (k+1)     s:Sigma (k+1)
      ! Beta=-------        =  -------------
      !       delta  (k)        s:Sigma (k)
      !
      ! delta=s:sigma
      !
      ! dKappa=s+Beta*dKappa

      deltaC = ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
      irc = ipOut(ipST)

      deltaK = ddot_(nDensC,Sigma,1,Sc2,1)
      if (.not. CI) then
        rBeta = deltaK/delta
        delta = deltaK
        call DScal_(nDensC,rBeta,dKappa,1)
        call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
      else
        rbeta = (deltac+deltaK)/delta
        delta = deltac+deltaK
        irc = ipIn(ipCID)
        call DScal_(nConf1*nroots,rBeta,W(ipCID)%Vec,1)
        call DScal_(nDensC,rBeta,dKappa,1)
        call DaXpY_(nConf1*nroots,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
        call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
        irc = opOut(ipS2)
        irc = ipOut(ipCID)
      end if

      !  ######  #    #  #####        #####    ####    ####
      !  #       ##   #  #    #       #    #  #    #  #    #
      !  #####   # #  #  #    #       #    #  #       #
      !  #       #  # #  #    #       #####   #       #  ###
      !  #       #   ##  #    #       #       #    #  #    #
      !  ######  #    #  #####        #        ####    ####
      !
      !----------------------------------------------------------------*

      res = Zero ! dummy initialize
      if (iBreak == 1) then
        if (abs(delta) < abs(Eps**2*delta0)) exit
      else if (iBreak == 2) then
        res = sqrt(resk**2+resci**2+ResQaaLag2)
        if (doDMRG) res = sqrt(resk**2)
        if (res < abs(Eps)) exit
      else
        if ((abs(delta) < abs(Eps**2*delta0)) .and. (res < abs(Eps))) exit
      end if
      if (iter >= niter) then
        cnvrgd = .false.
        exit
      end if
      if (lprint) write(u6,Fmt2//'I7,4X,ES17.9,ES17.9,ES17.9,ES17.9,ES17.9)') iter,delta/delta0,resk,resci,deltac,deltak
      iter = iter+1

    end do

    !*******************************************************************

    if (.not. cnvrgd) then
      write(u6,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
      converged(isym) = .false.
      fail = .true.
    else
      if (iPL >= 2) write(u6,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
      irc = ipnout(-1)
    end if

    if (iPL >= 2) write(u6,*)
    iLen = ndensC
    iKapDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Kappa,iLen,iDis)
    iSigDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)
    ilen = nconf1*nroots
    iCIDisp(iDisp) = iDis

    irc = ipin(ipCIT)
    call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)

    ! MGD This last call seems unused, so I comment it

    !call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,Temp4,ipS2)
    iCISigDisp(iDisp) = iDis
    irc = ipin(ipST)
    call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
  end do ! iDisp

  call mma_deallocate(Temp4)
  call mma_deallocate(Temp3)
  call mma_deallocate(Sigma)
  call mma_deallocate(dKappa)
  call mma_deallocate(Kappa)
  call mma_deallocate(Sc2)
  call mma_deallocate(Sc1)
end if

! Free all memory and remove from disk all data
! related to this symmetry

call mma_deallocate(Fancy)

irc = ipclose(ipdia)
if (.not. CI) irc = ipclose(ipPre2)

call Exp_Close()

if (debug) then
  write(u6,*) '********************************************************************************'
  write(u6,*)
end if
if (doDMRG) then  ! yma
  call dmrg_spc_change_mclr(RGras2(1:8),nash)
  call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer(irc)
#endif

end subroutine WfCtl_MSPD
