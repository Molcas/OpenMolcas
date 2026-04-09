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

subroutine WfCtl_td(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,iRHSCIDISP,converged)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use ipPage, only: ipclose, ipget, ipin, ipin1, ipnout, ipout, opout, W
use MCLR_Data, only: CMO, FIMO, Int2, ipCI, ipDia, lDisp, LuTemp, n1Dens, n2Dens, nConf1, nDens, nDensC, XISPSM
use MCLR_procedures, only: CISigma_td
use input_mclr, only: Debug, Eps, ERASSCF, Fail, iBreak, iMethod, kPrint, lCalc, lSave, nCSF, nDisp, nIter, nSym, nTPert, Omega, &
                      PotNuc, PT2, rIn_Ene, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp), &
                                  iRHSCIDisp(nDisp)
logical(kind=iwp), intent(out) :: converged(8)
integer(kind=iwp) :: iDEnd, iDis, iDisp, iLen, ipCID, ipCIT, ipPre2, ipS1, ipS2, ipST, iSym, Iter, jDisp, jSpin, kkkSym, kkSym, &
                     Left, lLine, lPaper, Lu_50, nConf3, niPre2, nPre2, pstate_sym
real(kind=wp) :: Clock(4), D_0, D_1, D_2, Delta, Delta0, DeltaC, DeltaK, EC, R1, R2, rAlpha, rAlphaC, rAlphaK, rBeta, rDum(1), &
                 ReCo, Res, rEsci, rEsk, rGrad, Tim2, Tim3, Tim4
logical(kind=iwp) :: CI, cnvrgd, lPrint, Orb
character(len=8) :: Fmt2
integer(kind=iwp), allocatable :: iPre(:)
real(kind=wp), allocatable :: Dens(:), DigPrec(:), dKappa(:), Kappa(:), Pens(:), rmoaa(:), Sc1(:), Sc2(:), Sc3(:), Sc4(:), &
                              Sigma(:), Temp1(:), Temp2(:), Temp3(:), Temp4(:), TempTD(:)
integer(kind=iwp), parameter :: iTimeCC = 1, iTimeKK = 2, iTimeKC = 3, iTimeCK = 4
integer(kind=iwp), external :: niPre, nPre
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
TIM2 = Zero
TIM3 = Zero
TIM4 = Zero

CLOCK(:) = Zero
lPaper = 132
lLine = 120
left = (lPaper-lLine)/2
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
iDis = 0

fail = .false.
Converged(:) = .true.
lprint = .false.
LU_50 = 50
if (lSAVE) call DANAME(LU_50,'RESIDUALS')
if (btest(kprint,1)) lprint = .true.
iDisp = 0
kksym = 1
kkksym = nsym
if (PT2) kkkSym = 1

! Starting loop over all symmetries/PT

do iSym=kksym,kkksym
  PState_SYM = Mul(State_Sym,iSym)
  nconf1 = ncsf(PState_Sym)
  CI = .false.
  if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.

  if (CI .and. (nconf1 == 1) .and. (isym == 1)) CI = .false.
  ! Initiate CSF <-> SD
  if (CI) call InCSFSD(Mul(iSym,State_Sym),State_sym)

  ! Calculate length of the density, Fock and Kappa matrix etc
  ! notice that this matrixes not necessary are symmetric.
  ! Store pointers.
  !
  ! Input:
  !        iSym : Symmetry of perturbation
  !
  ! Output: Commonblocks (Pointers.fh)

  PState_SYM = Mul(State_Sym,iSym)
  !nConf2 = nint(xispsm(PState_SYM,1))
  !nConf2 = ndtasm(PState_SYM)
  nconf3 = nint(max(xispsm(PState_SYM,1),xispsm(State_SYM,1)))
  !nconf3 = max(ndtasm(PState_SYM),ndtasm(State_SYM))

  ! Setup is realted to symmetry treatment

  call Setup_MCLR(iSym)

  ! Determine if we should page CI vectors

  !                            [2]
  ! Calculate the diagonal of E    and store in core/disc

  if (CI) then

    ! CIDia_td calculates the <i|H|i> elements (diagonal)? From the CI-CI part of E?

    call CIDia_td(PState_Sym)
    call ipout(ipdia)

    ! Allocate disk/memory space

    ! This areas should be addressed through ipin
    ! ipout will page them out to disk and free the memory area
    ! if page mode is used

    ! opout will release the memory area without update the disk

    ! ipget is like getmem but handles weather the vector is on disk or not.

    ips1 = ipget(2*nconf3)
    ips2 = ipget(2*nconf3)
    ipst = ipget(2*nconf3)
    ipcit = ipget(2*nconf1)
    ipcid = ipget(2*nconf1)

  end if

  npre2 = npre(isym)
  nipre2 = nipre(isym)
  call mma_allocate(iPre,nipre2,Label='iPre')

  !OBS npre2 not def

  ipPre2 = ipget(npre2)

  call ipin(ipPre2)
  call Prec_dig(W(ipPre2)%A,iPre,isym)

  call mma_allocate(DigPrec,nDensC,Label='DigPrec')
  DigPrec(:) = Zero
  call Prec_td(W(ipPre2)%A,DigPrec,isym)

  call ipout(ippre2)

  ! OK START WORKING

  iDEND = lDisp(iSym)
  !if ((SewLab == 'NONE') .and. (.not. mckinley)) iDEND = 1

  ! Loop over all PT of sym isym

  do jDisp=1,iDEnd
    iDisp = iDisp+1
    jspin = 0
    if (btest(nTPert(idisp),0)) jSpin = 1
    if (jspin == 0) then
      nconf1 = ncsf(PState_Sym)
    else
      nConf1 = nint(xispsm(Pstate_Sym,1))
    end if
    if (.not. lCalc(iDisp)) then
      converged(isym) = .false.
      cycle
    end if

    ! Allocate areas for scratch and state variables

    call mma_allocate(Kappa,nDensC,Label='Kappa')
    call mma_allocate(dKappa,nDens,Label='dKappa')
    call mma_allocate(Sigma,nDensC,Label='Sigma')
    call mma_allocate(Temp1,nDens,Label='Temp1')
    call mma_allocate(Temp2,nDensC,Label='Temp2')
    call mma_allocate(Temp3,nDens,Label='Temp3')
    call mma_allocate(Temp4,nDens,Label='Temp4')
    call mma_allocate(Sc1,nDens,Label='Sc1')
    call mma_allocate(Sc2,nDens,Label='Sc2')
    call mma_allocate(Sc3,nDens,Label='Sc3')
    call mma_allocate(Sc4,nDensC,Label='Sc4')
    Temp1(:) = Zero
    Kappa(:) = Zero
    dKappa(:) = Zero
    Sc1(:) = Zero
    if (CI) then
      call mma_allocate(Dens,n1Dens,Label='Dens')
      call mma_allocate(Pens,n2Dens,Label='Pens')
      call mma_allocate(rmoaa,n2Dens,Label='Pens')
      Dens(:) = Zero
      Pens(:) = Zero
      rmoaa(:) = Zero
    end if

    !-------------------------------------------------------------------
    !
    ! Calculate RHS for the perturbation. The b^x vector!
    ! sigma is the orbital part and Temp4 is the cipart
    !
    !-------------------------------------------------------------------

    ! (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)

    call RHS_td(Sc1,Temp1,Temp3,Sc2,dKappa,Sc3,Temp4,ipST,iDisp,iSym-1,CMO,jdisp,CI)

    Temp4(:) = -Temp4(:)

    ! Make RHS twice as long and change sign on second part!

    if (CI) then
      call ipin(ipST)
      W(ipST)%A(nConf1+1:2*nConf1) = -W(ipST)%A(1:nConf1)
    end if

    call opout(ipci)

    if (lprint) write(u6,*) '       Iteration         Delta     Res(kappa) Res(CI)'
    iLen = nDensC
    iRHSDisp(iDisp) = iDis
    call Compress(Temp4,Sigma,iSym)
    r1 = Half*ddot_(nDensC,Sigma,1,Sigma,1)

    call UnCompress(Sigma,Temp4,iSym)
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)
    if (CI) then
      call ipin(ipCIT)
      W(ipCIT)%A(1:2*nConf1) = Zero
    end if
    call ipout(ipcit)
    if (CI) then
      ilen = 2*nconf1
      iRHSCIDisp(iDisp) = iDis
      call ipin(ipST)
      call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
      W(ipST)%A(1:2*nConf1) = -W(ipST)%A(1:2*nConf1)
    end if

    call DMInvKap_td(DigPrec,Sigma,Kappa)

    call opout(ippre2)
    ! kap:kap
    r2 = Half*ddot_(nDensC,Kappa,1,Kappa,1)
    if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

    call mma_deallocate(Temp1)
    call mma_allocate(Temp1,nDensC,Label='Temp1')

    ! dkap=Kap in matrix form
    call UnCompress(Kappa,dKappa,iSym)

    ! DMinvCI_td is related to cidia_td and inverts the precond CI-CI
    ! part of E
    ! Has to be modified: <i|H|i> --> <i|H|i>+w

    if (CI) then
      call ipin(ipST)
      call ipin(ipCid)
      call DMinvCI_td(W(ipST)%A,W(ipCid)%A,-omega,isym)
      call DMinvCI_td(W(ipST)%A(1+nConf1),W(ipCId)%A(1+nconf1),omega,isym)

      deltaC = Half*ddot_(2*nConf1,W(ipST)%A,1,W(ipCId)%A,1)
      call ipout(ipcid)
    else
      deltac = Zero
    end if
    deltaK = Half*ddot_(nDensC,Kappa,1,Sigma,1)
    Kappa(:) = Zero
    delta = deltac+deltaK

    delta0 = delta
    Orb = .true.
    ReCo = -One
    iter = 1
    !-------------------------------------------------------------------

    !**********************************************************
    !          I   T   E   R   A   T   I   O   N   S          *
    !**********************************************************

    cnvrgd = .true.
    do
      if (delta == Zero) exit

      !*****************************************************************
      !
      !       O R B I T A L    P A R T of the trial vector
      !
      !*****************************************************************

      if (orb) then

        !---------------------------------------------------------------
        !
        !           ~    ~
        ! Construct F,(ij|kl)
        !
        !---------------------------------------------------------------
        !
        ! j2 specifies which part of E I want to look at
        ! j2=0 --> K-K, j2=-1 --> CI-CI, These are antisym within themself
        ! j2>0 --> CI-K and K-CI, These parts are antisym between eachother

        call ipnout(-1)
        call RInt_ns(dKappa,rmoaa,Sc2,Temp4,isym,reco,jspin)

        call RInt_td(Sc2,dKappa,isym)

        Clock(iTimeKK) = Clock(iTimeKK)+Tim2

        !---------------------------------------------------------------
        !
        ! kappa->CI
        !
        ! H(kappa)|0>
        !
        !     [2]
        ! S1=E   k (kappa TO CI) <i|H|0>-<i|0>*Energy
        !
        !---------------------------------------------------------------

        ! This cisigma_td call gives <j|H(k)|0> and <j|H(k)t|0>

        if (CI) then
          ! Adjusted to timedep
          call CISigma_td(jspin,State_Sym,pstate_sym,Temp4,nDens,rmoaa,size(rmoaa),rdum,1,ipCI,ipS1,'T',.true.)
          Clock(iTimeKC) = Clock(iTimeKC)+Tim3

          ! This will give us a better
          ! convergence in the PCG. Notice that
          !
          ! ~Inactive     ~
          ! E        + <0|H|0> = 0
          !
          ! when the wavefunction is converged.

          ! These terms are to be able to handle less converged CASSCF wave func

          call ipin(ipS1)
          if (isym == 1) then
            call ipin(ipCI)
            rGrad = ddot_(nconf1,W(ipCI)%A,1,W(ipS1)%A,1)
            W(ipS1)%A(1:nConf1) = W(ipS1)%A(1:nConf1)-rGrad*W(ipCI)%A(1:nConf1)
            rGrad = ddot_(nconf1,W(ipCI)%A,1,W(ipS1)%A(1+nconf1),1)
            W(ipS1)%A(nConf1+1:2*nConf1) = W(ipS1)%A(nConf1+1:2*nConf1)-rGrad*W(ipCI)%A(1:nConf1)
          end if
          W(ipS1)%A(1:nConf1) = -Two*W(ipS1)%A(1:nConf1)
          W(ipS1)%A(nConf1+1:2*nConf1) = Two*W(ipS1)%A(nConf1+1:2*nConf1)

          call opout(ipCI)
          !*************************************************************

        end if  ! If ci

      end if
      !
      !*****************************************************************
      !
      !    C I    P A R T of the trial vector
      !
      !*****************************************************************
      !
      if (CI) then

        !***************************************************************
        !
        !       [2]
        ! S2 = E   CID  (CI TO CI) = <i|H|d> -E<i|d>
        !
        !***************************************************************

        call ipnout(-1)
        if (CI) call CISigma_td(0,PState_Sym,Pstate_sym,FIMO,size(FIMO),Int2,size(Int2),rdum,1,ipCId,ipS2,'S',.true.)

        ! I want the RASSCF energy of the ACTIVE electrons !!!!
        ! EC=-E[act]           E[RASSCF]=E[inact]+E[act]+E[nuc]

        EC = rin_ene+potnuc-ERASSCF(1)

        call ipin(ipCId)
        call ipin(ipS2)
        W(ipS2)%A(1:2*nConf1) = Two*(W(ipS2)%A(1:2*nConf1)+EC*W(ipCId)%A(1:2*nConf1))

        ! Add the wS contribution
        ! The (-) sign in both daxpys assumes that the two parts of ipcid are def with diff sign.
        ! This is not true for the debug option!! ipcid = 1 regardless of part which part
        ! The S-contribution will make E-wS loose its symmetry because E is sym and S
        ! is antisym.

        W(ipS2)%A(1:nConf1) = W(ipS2)%A(1:nConf1)-Two*omega*W(ipCId)%A(1:nConf1)
        W(ipS2)%A(nConf1+1:2*nConf1) = W(ipS2)%A(nConf1+1:2*nConf1)+Two*omega*W(ipCId)%A(nConf1+1:2*nConf1)
        Clock(iTimeCC) = Clock(iTimeCC)+Tim4

        call ipout(ips2)
        call opout(ipcid)
        call opout(ipci)

        !---------------------------------------------------------------
        !
        ! CI -> Kappa
        !
        !      [2]
        ! SC3=E   CID   (SC3=F(<d|E|0>+<0|E|d>)
        !
        !---------------------------------------------------------------

        call ipnout(-1)

        call CIDens_TD(ipCid,PState_Sym,Pens,Dens)     ! Jeppe's

        ! density for inactive= 2(<d|0>+<0|d>)

        d_0 = Zero

        ! This is just for debugging purpose.
        ! When we use it for actual calculations d_0 == 0

        ! This if statement is just for better convergence! Grad term
        ! Leave this for later!

        if (isym == 1) then
          call ipin(ipCid)
          call ipin(ipci)
          d_1 = ddot_(nconf1,W(ipCid)%A,1,W(ipci)%A,1)
          d_2 = ddot_(nconf1,W(ipCid)%A(1+nConf1),1,W(ipci)%A,1)
          d_0 = d_1+d_2
        end if

        ! Fockgen gives the Fock matrix, MO integrals and one index transformed
        ! MO integrals

        call Fockgen_td(d_0,Dens,Pens,Sc3,isym)

        !---------------------------------------------------------------

      end if

      !*****************************************************************
      !
      ! Sc1  kappa-> kappa
      ! Sc3  CI -> kappa
      ! S1   kappa -> CI
      ! S2   CI -> CI
      ! dKap present step
      ! Kap  kappaX
      ! CIT  CIX
      ! CId  present step
      !*****************************************************************
      !
      ! Add together
      !
      !*****************************************************************

      call ipnout(-1)
      if (CI) then   ! if (.false.) then
        Sc1(:) = Sc2(:)+Sc3(:)
      else
        Sc1(:) = Sc2(:)
      end if
      call Compress(Sc1,Temp1,isym)   ! ds
      call Compress(dKappa,Temp2,isym) ! DX

      ! S1 + S2 --> S1

      if (CI) then  !If (.false.) then
        call ipin1(ipS1,2*nconf1)
        call ipin1(ipS2,2*nconf1)
        W(ipS1)%A(1:2*nConf1) = W(ipS1)%A(1:2*nConf1)+W(ipS2)%A(1:2*nConf1)
        call opout(ips2)
      end if

      !-----------------------------------------------------------------
      !
      !            ######   #####   #####
      !            #     # #     # #     #
      !            #     # #       #
      !            ######  #       #  ####
      !            #       #       #     #
      !            #       #     # #     #
      !            #        #####   #####
      !
      !-----------------------------------------------------------------
      !*****************************************************************
      !
      !            delta
      ! rAlpha=------------
      !        dKappa:dSigma
      !
      !-----------------------------------------------------------------
      rAlphaC = Zero
      rAlphaK = Zero
      if (orb) rAlphaK = Half*ddot_(nDensC,Temp1,1,Temp2,1)
      if (CI) then
        call ipin(ipS1)
        call ipin(ipCId)
        rAlphaC = Half*ddot_(2*nConf1,W(ipS1)%A,1,W(ipCId)%A,1)
      end if
      rAlpha = delta/(rAlphaK+ralphaC)

      !----------------------------------------------------------------*

      ! Kappa=Kappa+rAlpha*dKappa
      ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa

      if (orb) then
        Kappa(:) = Kappa(:)+ralpha*Temp2(:)
        Sigma(:) = Sigma(:)-ralpha*Temp1(:)
        resk = sqrt(Half*ddot_(nDensC,Sigma,1,Sigma,1))
      end if
      resci = Zero

      if (CI) then
        call ipin(ipCId)
        call ipin(ipCIT)
        W(ipCIT)%A(1:2*nConf1) = W(ipCIT)%A(1:2*nConf1)+ralpha*W(ipCId)%A(1:2*nConf1)
        call ipout(ipcit)
        call ipin1(ipST,2*nconf1)
        call ipin(ipS1)
        W(ipST)%A(1:2*nConf1) = W(ipST)%A(1:2*nConf1)-ralpha*W(ipS1)%A(1:2*nConf1)
        call opout(ipS1)
        resci = sqrt(Half*ddot_(2*nconf1,W(ipST)%A,1,W(ipST)%A,1))
      end if

      !----------------------------------------------------------------*
      ! Precondition......
      !    -1
      ! S=M  Sigma

      call opout(ipcid)
      if (CI) then
        call ipin(ipST)
        call ipin(ipS2)
        call DMinvCI_td(W(ipST)%A,W(ipS2)%A,-omega,isym)
        call DMinvCI_td(W(ipST)%A(1+nConf1),W(ipS2)%A(1+nconf1),omega,isym)

      end if
      call opout(ipci)
      call opout(ipdia)

      call DMInvKap_td(DigPrec,Sigma,Sc4)

      call opout(ippre2)

      !----------------------------------------------------------------*
      !      s:Sigma
      ! Beta=-------
      !       delta
      !
      ! delta=s:sigma
      !
      ! dKappa=s+Beta*dKappa

      if (CI) then
        call ipin(ipST)
        call ipin(ipS2)
        deltaC = Half*ddot_(2*nConf1,W(ipST)%A,1,W(ipS2)%A,1)
        call ipout(ipST)
      else
        deltaC = Zero
      end if

      deltaK = Half*ddot_(nDensC,Sigma,1,Sc4,1)
      if (.not. CI) then
        rBeta = deltaK/delta
        delta = deltaK
        Temp2(:) = rBeta*Temp2(:)+Sc4(:)
      else
        rbeta = (deltac+deltaK)/delta
        delta = deltac+deltaK
        call ipin(ipCID)
        call ipin(ipS2)
        W(ipCID)%A(1:2*nConf1) = rBeta*W(ipCID)%A(1:2*nConf1)+W(ipS2)%A(1:2*nConf1)
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
      !----------------------------------------------------------------*

      call UnCompress(Temp2,dKappa,isym)

      ! iBreak is defined via include!

      res = Zero ! dummy initialize
      if (iBreak == 1) then
        ! This is the actual breaking!
        if (abs(delta) < abs(Eps**2*delta0)) exit
      else if (ibreak == 2) then
        res = sqrt(resk**2+resci**2)
        if (res < abs(Eps)) exit
      else
        if ((abs(delta) < abs(Eps**2*delta0)) .and. (res < abs(Eps))) exit
      end if

      ! This breaks the PCG iterations

      if (iter >= niter) then
        cnvrgd = .false.
        exit
      end if
      if (lprint) &
        write(u6,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)') '     ',iter,'       ',delta/delta0,resk,resci,deltac,deltak

      iter = iter+1

    end do

    !*******************************************************************

    if (.not. cnvrgd) then
      write(u6,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
      converged(isym) = .false.
      fail = .true.
    else
      write(u6,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
      call ipnout(-1)
      !stop 10
    end if
    call mma_allocate(TempTD,nDens,Label='TempTD')
    call Uncompress(Kappa,TempTD,isym)
    call mma_deallocate(TempTD)

    write(u6,*)
    iLen = nDensC
    iKapDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Kappa,iLen,iDis)
    iSigDisp(iDisp) = iDis
    call dDaFile(LuTemp,1,Sigma,iLen,iDis)
    if (CI) then
      ilen = 2*nconf1
      iCIDisp(iDisp) = iDis
      call ipin(ipCIT)
      call dDaFile(LuTemp,1,W(ipCIT)%A,iLen,iDis)
      iCISigDisp(iDisp) = iDis
      call ipin(ipST)
      call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
    end if

    call mma_deallocate(Temp4)
    call mma_deallocate(Temp3)
    call mma_deallocate(Temp2)
    call mma_deallocate(Temp1)
    call mma_deallocate(Sigma)
    call mma_deallocate(dKappa)
    call mma_deallocate(Kappa)
    call mma_deallocate(Sc4)
    call mma_deallocate(Sc3)
    call mma_deallocate(Sc2)
    call mma_deallocate(Sc1)
    if (CI) then
      call mma_deallocate(rmoaa)
      call mma_deallocate(Pens)
      call mma_deallocate(Dens)
    end if
  end do

  ! Free all memory and remove from disk all data
  ! related to this symmetry

  call mma_deallocate(iPre)

  if (CI) then
    call ipclose(ipdia)
  else
    call ipclose(ipPre2)
  end if

  call mma_deallocate(DigPrec)
  call Exp_Close()

end do
if (debug) then
  write(u6,*) '********************************************************************************'
  write(u6,*)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine WfCtl_td
