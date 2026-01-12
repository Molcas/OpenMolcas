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

use Symmetry_Info, only: Mul
use ipPage, only: ipclose, ipget, ipin, ipnout, ipout, opout, W
use gugx, only: CIS, EXS, SGS
use MCLR_Data, only: ipCI, ipDia, IRLXROOT, ISNAC, LuQDat, LuTemp, NACSTATES, nConf1, nDens, nDensC, XISPSM
use input_mclr, only: Debug, Eps, Fail, iAddressQDat, iBreak, iMethod, iSpin, kPrint, lSave, nActEl, nAsh, nConf, nCSF, nDisp, &
                      nElec3, nHole1, nIter, NROOTS, nRS1, nRS2, nRS3, nSym, PT2, State_Sym, STEPTYPE, TWOSTEP
use PCM_grad, only: def_solv, do_RF, iStpPCM, PCM_grad_CLag, PCM_grad_PT2
use ISRotation, only: DMInvISR, InvSCF, ISR, ISR_final, ISR_init, ISR_projection, ISR_RHS
use cgs_mod, only: CGS, CGS_init, CGS_final
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp)
logical(kind=iwp), intent(out) :: converged(8)
integer(kind=iwp), intent(in) :: iPL
integer(kind=iwp) :: iDis, iDisp, iLen, ipCID, ipCIT, ipPre2, ipS1, ipS2, ipST, iR, iSym, Iter, jSpin, Left, lLine, lPaper, Lu_50, &
                     nConf3, niPre2, nPre2
real(kind=wp) :: Delta, Delta0, DeltaC, DeltaK, R1, R2, rAlpha, rAlphaC, rAlphaK, rBeta, ReCo, Res, rEsci, rEsk, rEss
logical(kind=iwp) :: CI, cnvrgd, lPrint
character(len=8) :: Fmt2
integer(kind=iwp), allocatable :: iPre(:)
real(kind=wp), allocatable :: dKappa(:), Fancy(:), Kappa(:), rCHC(:), Sc1(:), Sc2(:), Sigma(:), SLag(:,:), Temp3(:), Temp4(:), &
                              wrk(:)
integer(kind=iwp), external :: niPre, nPre
real(kind=wp), external :: DDot_

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
Converged(:) = .true.
!MGD I think this is nice when printed...
lprint = .true.
reco = -One
Lu_50 = 50
if (lSAVE) call DANAME(Lu_50,'RESIDUALS')
if (lSAVE) then
  write(u6,*) 'WfCtl_SA: SAVE option not implemented'
  call Abend()
end if
if (btest(kprint,1)) lprint = .true.
isym = 1
nconf1 = ncsf(State_Sym)
!! iStpPCM has been set somewhere; just a reminder
iStpPCM = 1

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

! Determine if we should page CI vectors
!                            [2]
! Calculate the diagonal of E    and store in core/disc

call mma_allocate(FANCY,nroots**3,Label='FANCY')
call mma_allocate(rCHC,nRoots,Label='rCHC')
call CIDia_SA(State_Sym,rCHC,Fancy)
call mma_deallocate(rCHC)

call ipOut(ipdia)

! initialize CGS and some for InvSCF/InvEne
call ISR_Init(iPL,do_RF,def_solv)

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

  !idisp = 1
  jspin = 0

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kappa,nDensC,Label='Kappa')
  call mma_allocate(dKappa,nDensC,Label='dKappa')
  call mma_allocate(Sigma,nDensC,Label='Sigma')
  call mma_allocate(Temp3,nDens,Label='Temp3')
  call mma_allocate(Temp4,nDensC,Label='Temp4')
  call mma_allocate(Sc1,nDens,Label='Sc1')
  call mma_allocate(Sc2,nDensC,Label='Sc2')
  if (CGS) call CGS_init(nconf1,nRoots,nDens)

  do iDisp=1,nDisp
    Kappa(:) = Zero
    dKappa(:) = Zero
    Sigma(:) = Zero
    Sc1(:) = Zero

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
    call mma_allocate(SLag,nRoots,nRoots,Label='SLag')
    SLag(:,:) = Zero
    if (PT2) then
      call RHS_PT2(Sc1,W(ipST)%A,SLag)
      ! The order of CSF coefficients in CASPT2 and MCLR is somehow
      ! different, so the CI lagrangian computed in CASPT2 must be
      ! reordered so that it can be used here.
      call mma_allocate(wrk,nConf1,Label='wrk')
      do iR=1,nRoots
        wrk(:) = W(ipST)%A(nConf1*(iR-1)+1:nConf1*iR)
        call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,wrk,1,State_Sym,State_Sym)
        NCSF(1:nSym) = CIS%NCSF(1:nSym)
        NCONF = CIS%NCSF(State_Sym)
        call mkGuga_Free(SGS,CIS,EXS)
        W(ipST)%A(nConf1*(iR-1)+1:nConf1*iR) = wrk(:)
      end do
      call mma_deallocate(wrk)
      ! SLag given by the CASPT2 routine should be treated as the
      ! initial rotation of the reference states.
      ! If we need to treat parameters as internal rotations, they
      ! will be computed later using configuration parameters.
      if (do_RF) call PCM_grad_PT2()
    else
      W(ipST)%A(1:nconf3*nRoots) = Zero
    end if

    if (do_RF .or. PT2) then
      !! SA-CASSCF/PCM may not be stationary wrt CI coefficients
      !! Use ipS2 for the moment (some intermediate vectors are SD but output vectors are CSF)
      if (do_RF) then
        call PCM_grad_CLag(1,ipCI,ipS2)
        W(ipS2)%A(1:nconf1*nRoots) = W(ipS2)%A(1:nconf1*nRoots)+W(ipST)%A(1:nconf1*nRoots)
      else if (PT2) then
        W(ipS2)%A(1:nconf1*nRoots) = W(ipST)%A(1:nconf1*nRoots)
      end if
      if (InvSCF) then
        call ISR_RHS(W(ipCI)%A,W(ipS2)%A)
        call ISR_Projection(W(ipCI)%A,W(ipS2)%A)
        SLag(:,:) = SLag(:,:)+ISR%Rvec(:,:)
        ISR%Rvec(:,:) = Zero
      end if
      W(ipST)%A(1:nconf1*nRoots) = W(ipS2)%A(1:nconf1*nRoots)
      W(ipS2)%A(1:nconf3*nRoots) = Zero
      ! If SLag is nonzero, the contribution comes from the rotation of the model state (i.e., MS-type CASPT2)
      ! but should not be taken as the initial residue of the internal rotation, so SLag should not be added to Rvec.
      !if (.not. InvSCF) ISR%Rvec(:,:) = ISR%Rvec(:,:)+SLag
    end if

    if (isNAC) then
      call RHS_NAC(Temp3,SLag,ipS2)
    else
      call RHS_SA(Temp3,SLag,ipS2)
    end if
    call mma_deallocate(SLag)

    ! Add implicit CI derivative contributions (if evaluated in rhs_sa/nac)
    ! due to the non-iterative internal state rotations
    if (do_RF) then
      W(ipST)%A(1:nconf1*nRoots) = W(ipST)%A(1:nconf1*nRoots)+W(ipS2)%A(1:nconf1*nRoots)
      call ISR_RHS(W(ipCI)%A,W(ipST)%A)
      call ISR_Projection(W(ipCI)%A,W(ipST)%A)
    end if

    ! Temp3: SA-CASSCF contributions
    ! Sc1  : CASPT2 contributions
    if (PT2) then
      Temp3(:) = Temp3(:)+Sc1(:)
      Sc1(:) = Zero
    end if

    call opOut(ipci)

    if (lprint) write(u6,*) '       Iteration       Delta           Res(kappa)       Res(CI)          DeltaK           DeltaC'
    iLen = nDensC
    iRHSDisp(iDisp) = iDis
    call Compress(Temp3,Sigma,iSym)
    r1 = ddot_(nDensC,Sigma,1,Sigma,1)
    if (PT2 .or. do_RF) R1 = R1+DDot_(nConf1*nRoots,W(ipST)%A,1,W(ipST)%A,1)
    if (debug) write(u6,*) 'Hi how about r1',r1
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)

    if (PT2 .or. do_RF) then !! .or. unequal_SA) then
      W(ipST)%A(1:nConf1*nRoots) = -W(ipST)%A(1:nConf1*nRoots)
      if (CI) then
        ! precondition (z0 = M^{-1}*r0)
        if (CGS) then
          W(ipS2)%A(1:nConf1*nRoots) = W(ipST)%A(1:nConf1*nRoots)
        else
          call DMinvCI_sa(ipST,W(ipS2)%A,fancy)
        end if
        call opOut(ipci)
        call opOut(ipdia)
        ! z0 <= p0
        W(ipCId)%A(1:nConf1*nRoots) = W(ipS2)%A(1:nConf1*nRoots)
      end if
    else
      call ipIn(ipCIT)
      call ipIn(ipST)
      call ipIn(ipCID)
      W(ipCIT)%A(1:nConf1*nroots) = Zero
      W(ipST)%A(1:nConf1*nroots) = Zero
      W(ipCID)%A(1:nConf1*nroots) = Zero
    end if
    call ipOut(ipCIT)
    Sigma(:) = -Sigma(:)

    iter = 1
    call ipIn(ipPre2)
    if (CGS) then
      Kappa(1:nDensC) = Sigma(1:nDensC)
    else
      call DMInvKap(W(ipPre2)%A,iPre,Sigma,Kappa,Temp3,isym,iter)
    end if

    call opOut(ippre2)
    r2 = ddot_(nDensC,Kappa,1,Kappa,1)
    if (PT2) R2 = R2+DDot_(nConf1*nRoots,W(ipS2)%A,1,W(ipS2)%A,1)
    if (debug) write(u6,*) 'In that case I think that r2 should be:',r2
    if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

    dKappa(:) = Kappa(:)

    if (.not. InvSCF) then
      ISR%Rvec(:,:) = -ISR%Rvec(:,:)
      if (.not. CGS) then
        call DMInvISR(ISR%Rvec,ISR%prec)
        ISR%p(:,:) = ISR%prec(:,:)
      end if
    end if

    if (CGS) then
      call cgs_pre(nDensC,nConf1,nRoots,Kappa,dKappa,Sigma,Fancy,Sc1,ipCI,ipCId,ipDia,ipPre2,iPre,ipST,ipS1,delta,iSym,iter)
    else
      deltaC = Zero
      if (PT2) deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipS2)%A,1)
      call ipOut(ipcid)
      deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
      Kappa(:) = Zero
      delta = deltac+deltaK
    end if
    delta0 = delta
    iter = 1
    iStpPCM = 2
    !-------------------------------------------------------------------

    cnvrgd = .true.
    do
      if (delta == Zero) exit
      if (CGS) then
        !! (Preconditioned) Conjugate Gradient Squared (CGS) method, instead of the conventional PCG
        call cgs_x(nDensC,nConf1,nRoots,iSym,jspin,iter,ipCI,ipDia,ipS1,ipS2,ipST,ipCIT,ipCId,ipPre2,reco,Fancy,Kappa,dKappa, &
                   Sigma,Temp4,Sc1,iPre,delta,resk,resci,deltak,deltac)
        !! end of CGS
      else
        !! standard PCG starts here

        call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

        !---------------------------------------------------------------
        !
        !            delta
        ! rAlpha=------------
        !        dKappa:dSigma
        !
        !---------------------------------------------------------------

        rAlphaK = Zero
        rAlphaK = ddot_(nDensC,Temp4,1,dKappa,1)
        rAlphaC = Zero
        call ipIn(ipS1)
        call ipIn(ipCId)
        rAlphaC = ddot_(nConf1*nroots,W(ipS1)%A,1,W(ipCId)%A,1)
        if (.not. InvSCF) rAlphaC = rAlphaC+ddot_(nRoots**2,ISR%Ap,1,ISR%p,1)
        rAlpha = delta/(rAlphaK+rAlphaC)

        !--------------------------------------------------------------*

        ! Kappa=Kappa+rAlpha*dKappa
        Kappa(:) = Kappa(:)+ralpha*dKappa(:)
        ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
        Sigma(:) = Sigma(:)-ralpha*Temp4(:)
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
        if (.not. InvSCF) then
          ISR%Xvec(:,:) = ISR%Xvec(:,:)+ralpha*ISR%p(:,:)
          ISR%Rvec(:,:) = ISR%Rvec(:,:)-ralpha*ISR%Ap(:,:)
          ress = sqrt(ddot_(nRoots**2,ISR%Rvec,1,ISR%Rvec,1))
          resci = resci+ress
        end if

        !--------------------------------------------------------------*
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

        if (.not. InvSCF) call DMInvISR(ISR%Rvec,ISR%prec)

        !--------------------------------------------------------------*
        !      s:Sigma (k+1)     s:Sigma (k+1)
        ! Beta=-------        =  -------------
        !       delta  (k)        s:Sigma (k)
        !
        ! delta=s:sigma
        !
        ! dKappa=s+Beta*dKappa

        deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipS2)%A,1)
        call ipOut(ipST)
        if (.not. InvSCF) deltaC = deltaC+ddot_(nRoots**2,ISR%Rvec,1,ISR%prec,1)

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
          if (.not. InvSCF) ISR%p(:,:) = ISR%prec(:,:)+rbeta*ISR%p(:,:)
        end if
      end if ! CGS vs PCG

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
    iLen = nDensC
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

    !! Save the total internal rotations so that CIDens_sa can use
    !! the parameters to generate correct densities
    if (.not. InvSCF) ISR%p(:,:) = ISR%Xvec(:,:)
  end do ! iDisp

  call mma_deallocate(Temp4)
  call mma_deallocate(Temp3)
  call mma_deallocate(Sigma)
  call mma_deallocate(dKappa)
  call mma_deallocate(Kappa)
  call mma_deallocate(Sc2)
  call mma_deallocate(Sc1)
end if

if (CGS) call CGS_final()
call ISR_final()

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
iStpPCM = 3

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine WfCtl_SA
