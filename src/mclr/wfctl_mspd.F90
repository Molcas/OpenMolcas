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

use Symmetry_Info, only: Mul
use ipPage, only: ipclose, ipget, ipin, ipnout, ipout, opout, W
use MCLR_Data, only: ipCI, ipDia, IRLXROOT, ISMECIMSPD, ISNAC, LuQDat, LuTemp, NACSTATES, nConf1, nDens, nDensC, OVERRIDE, &
                     ResQaaLag2, XISPSM
use input_mclr, only: Debug, Eps, Fail, iAddressQDat, iBreak, iMethod, kPrint, lSave, nAsh, nCSF, nDisp, nIter, nRoots, nRS2, &
                      State_Sym, StepType, TwoStep
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp)
logical(kind=iwp), intent(out) :: converged(8)
integer(kind=iwp), intent(in) :: iPL
#include "warnings.h"
integer(kind=iwp) :: iDis, iDisp, iLen, ipCID, ipCIT, ipPre2, ipS1, ipS2, ipST, iSym, Iter, JRoot, jSpin, Left, lLine, lPaper, &
                     Lu_50, LURot, nConf3, niPre2, nPre2
real(kind=wp) :: Delta, Delta0, DeltaC, DeltaK, R1, R2, rAlpha, rAlphaC, rAlphaK, rBeta, ReCo, Res, rEsci, rEsk
logical(kind=iwp) :: CI, cnvrgd, lPrint
character(len=16) :: VecName
character(len=8) :: Fmt2
integer(kind=iwp), allocatable :: iPre(:)
real(kind=wp), allocatable :: dKappa(:), Fancy(:), Kappa(:), rCHC(:), Sc1(:), Sc2(:), Sigma(:), Temp3(:), Temp4(:)
integer(kind=iwp), external :: IsFreeUnit, niPre, nPre
real(kind=wp), external :: DDot_

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
Converged(:) = .true.
!MGD I think this is nice when printed...
lprint = .true.
reco = -One
Lu_50 = 50
if (lSAVE) then
  call DANAME(Lu_50,'RESIDUALS')
  write(u6,*) 'WfCtl_MSPD: SAVE option not implemented'
  call Abend()
end if
if (btest(kprint,1)) lprint = .true.
isym = 1
nconf1 = ncsf(State_Sym)

CI = .false.
if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.

! Initiate CSF <-> SD
call InCSFSD(Mul(iSym,State_Sym),State_sym)

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
call mma_allocate(rCHC,nRoots,Label='rCHC')
call CIDia_SA(State_Sym,rCHC,Fancy)
call mma_deallocate(rCHC)

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

nipre2 = nipre(isym)
npre2 = npre(isym)
ipPre2 = ipGet(npre2)
call ipIn(ipPre2)
call mma_allocate(iPre,nipre2,Label='iPre')
if (TwoStep .and. (StepType == 'RUN2')) then
  ! fetch data from LuQDAT and skip the call to "Prec"
  call ddafile(LuQDAT,2,W(ipPre2)%A,npre2,iaddressQDAT)
  call idafile(LuQDAT,2,iPre,nipre2,iaddressQDAT)
else
  call Prec(W(ipPre2)%A,iPre,isym)
  call ipOut(ippre2)
end if
if (TwoStep .and. (StepType == 'RUN1')) then
  ! save the computed data in "Prec" to LuQDAT and skip the
  ! following part of this function
  call ipIn(ipPre2)
  call ddafile(LuQDAT,1,W(ipPre2)%A,npre2,iaddressQDAT)
  call idafile(LuQDAT,1,iPre,nipre2,iaddressQDAT)
else

  ! OK START WORKING
  !
  !idisp = 1
  jspin = 0

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kappa,nDensC,Label='Kappa')
  call mma_allocate(dKappa,nDensC,Label='dKappa')
  call mma_allocate(Sigma,nDensC,Label='Sigma')
  call mma_allocate(Temp3,nDensC,Label='Temp3')
  call mma_allocate(Temp4,nDens,Label='Temp4')
  call mma_allocate(Sc1,nDens,Label='Sc1')
  call mma_allocate(Sc2,nDensC,Label='Sc2')

  cnvrgd = .true.
  do iDisp=1,nDisp
    Kappa(:) = Zero
    dKappa(:) = Zero
    Sigma(:) = Zero

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
        call RHS_CMS_NAC(Temp4,W(ipST)%A)
        call DMinvCI_SA(ipST,W(ipS2)%A,Fancy)
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
        call RHS_CMS(Temp4,W(ipST)%A)
        call DMinvCI_SA(ipST,W(ipS2)%A,Fancy)
      else
        write(u6,'(6X,A)') 'Error: Lagrangian Not Implemented for MS-PDFT'
        write(u6,'(6X,A)') '       Other Than CMS-PDFT'
        call xQuit(_RC_INPUT_ERROR_)
      end if
    end if

    call opOut(ipci)

    if (lprint) write(u6,*) '       Iteration       Delta           Res(kappa)       Res(CI)          DeltaK           DeltaC'
    iLen = nDensC
    iRHSDisp(iDisp) = iDis
    call Compress(Temp4,Sigma,iSym)
    r1 = ddot_(nDensc,Sigma,1,Sigma,1)
    if (debug) write(u6,*) 'Hi how about r1',r1
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)

    call ipIn(ipCIT)
    !call ipIn(ipST)
    call ipIn(ipCID)
    W(ipCIT)%A(1:nConf1*nroots) = Zero
    ! already initialized in rhs_mspdft
    !W(ipST)%A(1:nConf1*nroots) = Zero
    !W(ipCID)%A(1:nConf1*nroots) = Zero
    !W(ipCID)%A(1:nConf1*nroots) = W(ipST)%A(1:nConf1*nroots)
    call ipOut(ipCIT)
    Sigma(:) = -Sigma(:)

    call ipIn(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,Sigma,Kappa,Temp4,isym,iter)

    call opOut(ippre2)
    r2 = ddot_(nDensC,Kappa,1,Kappa,1)
    if (debug) write(u6,*) 'In that case I think that r2 should be:',r2
    if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

    dKappa(:) = Kappa(:)

    ! In MS-PDFT deltaC is no longer zero initially
    !deltaC = Zero
    ! Use following two lines instead
    call DMinvCI_SA(ipST,W(ipS2)%A,Fancy)
    W(ipCID)%A(1:nConf1*nroots) = W(ipS2)%A(1:nConf1*nroots)
    deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipS2)%A,1)
    call ipOut(ipcid)
    deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
    Kappa(:) = Zero
    delta = deltac+deltaK
    delta0 = delta
    iter = 1
    !-------------------------------------------------------------------

    do
      if (delta == Zero) exit
      call TimesE2MSPDFT(dKappa,ipCId,1,reco,jspin,ipS2,Temp3,ipS1)

      !-----------------------------------------------------------------
      !
      !            delta
      ! rAlpha=------------
      !        dKappa:dSigma
      !
      !-----------------------------------------------------------------

      rAlphaK = Zero
      rAlphaK = ddot_(nDensC,Temp3,1,dKappa,1)
      rAlphaC = Zero
      call ipIn(ipS1)
      call ipIn(ipCId)
      rAlphaC = ddot_(nConf1*nroots,W(ipS1)%A,1,W(ipCId)%A,1)
      rAlpha = delta/(rAlphaK+rAlphaC)

      !----------------------------------------------------------------*

      ! Kappa=Kappa+rAlpha*dKappa
      Kappa(:) = Kappa(:)+ralpha*dKappa(:)
      ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
      Sigma(:) = Sigma(:)-ralpha*Temp3(:)
      resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
      resci = Zero
      call ipIn(ipCIT)
      W(ipCIT)%A(1:nConf1*nroots) = W(ipCIT)%A(1:nConf1*nroots)+ralpha*W(ipCId)%A(1:nConf1*nroots)
      call ipOut(ipcit)
      ! ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
      call ipIn(ipS1)
      call ipIn(ipST)
      W(ipST)%A(1:nConf1*nroots) = W(ipST)%A(1:nConf1*nroots)-ralpha*W(ipS1)%A(1:nConf1*nroots)
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
      call DMInvKap(W(ipPre2)%A,iPre,Sigma,Sc2,Sc1,iSym,iter)
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
        dKappa(:) = rBeta*dKappa(:)+Sc2(:)
      else
        rbeta = (deltac+deltaK)/delta
        delta = deltac+deltaK
        call ipIn(ipCID)
        W(ipCID)%A(1:nConf1*nroots) = rBeta*W(ipCID)%A(1:nConf1*nroots)+W(ipS2)%A(1:nConf1*nroots)
        dKappa(:) = rBeta*dKappa(:)+Sc2(:)
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
      call ipnout(-1)
    end if

    if (iPL >= 2) write(u6,*)
    iLen = nDensC
    iKapDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Kappa,iLen,iDis)
    iSigDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)
    ilen = nconf1*nroots
    iCIDisp(iDisp) = iDis

    call ipin(ipCIT)
    call dDaFile(LuTemp,1,W(ipCIT)%A,iLen,iDis)

    ! MGD This last call seems unused, so I comment it

    !call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,Temp3,ipS2)
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
call mma_deallocate(iPre)

call ipclose(ipdia)
if (.not. CI) call ipclose(ipPre2)

if (debug) then
  write(u6,*) '********************************************************************************'
  write(u6,*)
end if
if (doDMRG) then  ! yma
  nash(:) = RGras2(:)
  nrs2(:) = RGras2(:)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine WfCtl_MSPD
