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

subroutine WfCtl_SA(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,converged,iPL)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!***********************************************************************

use ipPage, only: ipclose, ipget, ipin, ipnout, ipout, opout, W
use gugx, only: SGS, CIS, EXS
use MCLR_Data, only: nConf1, nDens2, nDensC, nDens, ipCI
use MCLR_Data, only: ipDia
use MCLR_Data, only: ISNAC, IRLXROOT, NACSTATES
use MCLR_Data, only: LuTemp, LuQDat
use MCLR_Data, only: XISPSM
use input_mclr, only: nDisp, Fail, lSave, nSym, PT2, State_Sym, iMethod, iBreak, Eps, nIter, iSpin, Debug, kPrint, nCSF, &
                      iAddressQDat, NROOTS, TWOSTEP, STEPTYPE, nConf, nActEl, nAsh, nElec3, nHole1, nRS1, nRS2, nRS3
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
logical CI
character(len=8) Fmt2
logical lPrint, cnvrgd
real*8 rchc(mxroot)
real*8, allocatable :: Kappa(:), dKappa(:), Sigma(:), Temp3(:), Temp4(:), Sc1(:), Sc2(:), Fancy(:), SLag(:), wrk(:)
real*8 R1, R2, DeltaC, DeltaK, Delta, Delta0, ReCo, rAlphaC, rAlphaK, rAlpha, rEsk, rEsci, rBeta, Res
real*8, external :: DDot_
integer lPaper, lLine, Left, iDis, Lu_50, iDisp, iSym, nConf3, ipS1, ipS2, ipST, ipCIT, ipCID, nPre2, iLen, Iter, ipPre2, jSpin, &
        i, iR
integer, external :: nPre

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
call StatusLine('MCLR: ','Computing Lagrangian multipliers for SA-CASSCF')

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
if (lSAVE) call DANAME(Lu_50,'RESIDUALS')
if (lSAVE) then
  write(u6,*) 'WfCtl_SA: SAVE option not implemented'
  call Abend()
end if
if (iand(kprint,2) == 2) lprint = .true.
isym = 1
nconf1 = ncsf(State_Sym)

CI = .false.
if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.

! Initiate CSF <-> SD
call InCSFSD(ieor(iSym-1,State_Sym-1)+1,State_sym)

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

! Determine if we should page CI vectors
!                            [2]
! Calculate the diagonal of E    and store in core/disc

call mma_allocate(FANCY,nroots**3,Label='FANCY')
call CIDia_SA(State_Sym,rCHC,Fancy)

call ipOut(ipdia)

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
call ipIn(ipPre2)
if (TwoStep .and. (StepType == 'RUN2')) then
  ! fetch data from LuQDAT and skip the call to "Prec"
  call ddafile(LuQDAT,2,W(ipPre2)%A,npre2,iaddressQDAT)
else
  call Prec(W(ipPre2)%A,isym)
  call ipOut(ippre2)
end if
if (TwoStep .and. (StepType == 'RUN1')) then
  ! save the computed data in "Prec" to LuQDAT and skip the
  ! following part of this function
  call ipIn(ipPre2)
  call ddafile(LuQDAT,1,W(ipPre2)%A,npre2,iaddressQDAT)
else

  ! OK START WORKING

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
    call mma_allocate(SLag,nRoots**2,Label='SLag')
    SLag(1:nRoots**2) = Zero
    if (PT2) call RHS_PT2(Kappa,W(ipST)%A,Slag)

    if (isNAC) then
      call RHS_NAC(Temp4,SLag)
    else
      call RHS_SA(Temp4,SLag)
    end if

    if (PT2) then
      call DaXpY_(nDens2,One,Kappa,1,Temp4,1)
      Kappa(1:nDens2) = Zero
    end if
    call mma_deallocate(SLag)
    call opOut(ipci)

    if (lprint) write(u6,*) '       Iteration       Delta           Res(kappa)       Res(CI)          DeltaK           DeltaC'
    iLen = nDensC
    iRHSDisp(iDisp) = iDis
    call Compress(Temp4,Sigma,iSym)
    r1 = ddot_(nDensc,Sigma,1,Sigma,1)
    if (PT2) R1 = R1+DDot_(nConf1*nRoots,W(ipST)%A,1,W(ipST)%A,1)
    if (debug) write(u6,*) 'Hi how about r1',r1
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)

    if (PT2) then
      call DSCAL_(nConf1*nRoots,-One,W(ipST)%A,1)
      if (CI) then
        ! The order of CSF coefficients in CASPT2 and MCLR is somehow
        ! different, so the CI lagrangian computed in CASPT2 must be
        ! reordered so that it can be used here.
        call mma_allocate(wrk,nConf1,Label='wrk')
        do iR=1,nRoots
          call DCopy_(nConf1,W(ipST)%A(1+nConf1*(iR-1):nConf1*iR),1,wrk,1)
          call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,wrk,1,State_Sym,State_Sym)
          NCSF(1:nSym) = CIS%NCSF(1:nSym)
          NCONF = CIS%NCSF(State_Sym)
          call mkGuga_Free(SGS,CIS,EXS)

          call DCopy_(nConf1,wrk,1,W(ipST)%A(1+nConf1*(iR-1):nConf1*iR),1)
        end do
        call mma_deallocate(wrk)

        ! precondition (z0 = M^{-1}*r0)
        call DMinvCI_sa(ipST,W(ipS2)%A,fancy)
        call opOut(ipci)
        call opOut(ipdia)
        ! z0 <= p0
        call DCopy_(nConf1*nRoots,W(ipS2)%A,1,W(ipCId)%A,1)
      end if
    else
      call ipIn(ipCIT)
      call ipIn(ipST)
      call ipIn(ipCID)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%A,1)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipST)%A,1)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%A,1)
    end if
    call ipOut(ipCIT)
    call DSCAL_(nDensC,-One,Sigma,1)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,Sigma,nDens2+6,Kappa,nDens2+6,Temp3,nDens2+6,isym,iter)

    call opOut(ippre2)
    r2 = ddot_(ndensc,Kappa,1,Kappa,1)
    if (PT2) R2 = R2+DDot_(nConf1*nRoots,W(ipS2)%A,1,W(ipS2)%A,1)
    if (debug) write(u6,*) 'In that case I think that r2 should be:',r2
    if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

    call dcopy_(ndensC,Kappa,1,dKappa,1)

    deltaC = Zero
    if (PT2) deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipS2)%A,1)
    call ipOut(ipcid)
    deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
    Kappa(1:nDens) = Zero
    delta = deltac+deltaK
    delta0 = delta
    iter = 1
    !-------------------------------------------------------------------

    cnvrgd = .true.
    do
      if (delta == Zero) exit

      call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

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
      call ipIn(ipS1)
      call ipIn(ipCId)
      rAlphaC = ddot_(nConf1*nroots,W(ipS1)%A,1,W(ipCId)%A,1)
      rAlpha = delta/(rAlphaK+rAlphaC)

      !----------------------------------------------------------------*

      ! Kappa=Kappa+rAlpha*dKappa
      call DaxPy_(nDensC,ralpha,dKappa,1,Kappa,1)
      ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
      call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
      resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
      resci = Zero
      call ipIn(ipCIT)
      call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%A,1,W(ipCIT)%A,1)
      call ipOut(ipcit)
      ! ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
      call ipIn(ipS1)
      call ipIn(ipST)
      call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%A,1,W(ipST)%A,1)
      call opOut(ipS1)
      call ipIn(ipST)
      resci = sqrt(ddot_(nconf1*nroots,W(ipST)%A,1,W(ipST)%A,1))

      !----------------------------------------------------------------*
      ! Precondition......
      !    -1
      ! S=M  Sigma

      call opOut(ipcid)

      call ipIn(ipS2)
      call DMinvCI_SA(ipST,W(ipS2)%A,Fancy)
      call opOut(ipci)
      call opOut(ipdia)

      call ipIn(ipPre2)
      call DMInvKap(W(ipPre2)%A,Sigma,nDens2+6,Sc2,nDens2+6,Sc1,nDens2+6,iSym,iter)
      call opOut(ippre2)

      !----------------------------------------------------------------*
      !      s:Sigma (k+1)     s:Sigma (k+1)
      ! Beta=-------        =  -------------
      !       delta  (k)        s:Sigma (k)
      !
      ! delta=s:sigma
      !
      ! dKappa=s+Beta*dKappa

      deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipS2)%A,1)
      call ipOut(ipST)

      deltaK = ddot_(nDensC,Sigma,1,Sc2,1)
      if (.not. CI) then
        rBeta = deltaK/delta
        delta = deltaK
        call DScal_(nDensC,rBeta,dKappa,1)
        call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
      else
        rbeta = (deltac+deltaK)/delta
        delta = deltac+deltaK
        call ipIn(ipCID)
        call DScal_(nConf1*nroots,rBeta,W(ipCID)%A,1)
        call DScal_(nDensC,rBeta,dKappa,1)
        call DaXpY_(nConf1*nroots,One,W(ipS2)%A,1,W(ipCID)%A,1)
        call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
        call opOut(ipS2)
        call ipOut(ipCID)
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
        res = sqrt(resk**2+resci**2)
        if (doDMRG) res = sqrt(resk**2)
        if (res < abs(Eps)) exit
      else
        if ((abs(delta) < abs(Eps**2*delta0)) .and. (res < abs(Eps))) exit
      end if
      if (iter >= niter) then
        cnvrgd = .false.
        exit
      end if
      if (lprint) write(u6,Fmt2//'I7,4X,ES17.9,ES17.9,ES17.9,ES17.9,ES17.9)') iter,delta/delta0,resk,resci,deltak,deltac
      iter = iter+1

    end do

    !*******************************************************************

    if (.not. cnvrgd) then
      write(u6,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
      converged(isym) = .false.
      fail = .true.
    else
      if (iPL >= 2) write(u6,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
      call ipnout(-1)
    end if

    if (iPL >= 2) write(u6,*)
    iLen = ndensC
    iKapDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Kappa,iLen,iDis)
    iSigDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)
    ilen = nconf1*nroots
    iCIDisp(iDisp) = iDis

    call ipin(ipCIT)
    call dDaFile(LuTemp,1,W(ipCIT)%A,iLen,iDis)

    !MGD This last call seems unused, so I comment it

    !call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,Temp4,ipS2)
    iCISigDisp(iDisp) = iDis
    call ipin(ipST)
    call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
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

call ipclose(ipdia)
if (.not. CI) call ipclose(ipPre2)

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

end subroutine WfCtl_SA
