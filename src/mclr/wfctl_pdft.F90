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
! Copyright (C) 2017, Andrew M. Sand                                   *
!***********************************************************************

subroutine WfCtl_pdft(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,converged,iPL)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use ipPage, only: ipclose, ipget, ipin, ipnout, ipout, opout, W
use MCLR_Data, only: Do_Hybrid, ipCI, ipDia, ipMat, IRLXROOT, ISNAC, LuTemp, nAcPar, nAcPr2, NACSTATES, nConf1, nDens, nDensC, &
                     nNA, PDFT_Ratio, WF_Ratio, XISPSM
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: Debug, Eps, ERASSCF, Fail, iBreak, iMethod, kPrint, lSave, nAsh, nBas, nCSF, nDisp, nIter, nRoots, nRs2, &
                      nSym, ntAsh, State_Sym, Weight
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp)
logical(kind=iwp), intent(out) :: converged(8)
integer(kind=iwp), intent(in) :: iPL
integer(kind=iwp) :: i, iDis, iDisp, ij, iLen, iOff, ipCID, ipCIT, ipPre2, ipS1, ipS2, ipST, iS, iSym, Iter, j, ji, jS, jSpin, &
                     kSym, Left, lLine, lPaper, Lu_50, nConf3, nG1, nG2, nOrbAct, nPre2, nTri
real(kind=wp) :: Delta, Delta0, DeltaC, DeltaK, Diff, R1, R2, rAlpha, rAlphaC, rAlphaK, rBeta, rDum(1), rE, ReCo, Res, rEsci, &
                 rEsk, TRoot, WScale
logical(kind=iwp) :: CI, cnvrgd, lPrint
character(len=8) :: Fmt2
integer(kind=iwp), allocatable :: iPre(:)
real(kind=wp), allocatable :: dKappa(:), Fancy(:), FMO1(:), FMO1t(:), FMO2t(:), FOSq(:), FOTr(:), FT99(:), Kap_New(:), &
                              Kap_New_Temp(:), Kappa(:), lmroots(:), lmroots_new(:), P2PDFT(:), P2WF(:), rCHC(:), Sc1(:), Sc2(:), &
                              Sigma(:), Temp4(:), Temp5(:), WFOrb(:)
real(kind=wp), external :: DDot_
integer(kind=iwp), external :: niPre, nPre

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
call StatusLine('MCLR: ','Computing Lagrangian multipliers for MC-PDFT')

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
debug = .false.
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

call mma_allocate(Fancy,nRoots**3,Label='Fancy')
call mma_allocate(rCHC,nRoots,Label='rCHC')
!____________________
!AMS - what to do here?
! What should rCHC be?  is it computed with E(mcscf) or E(pdft)?
!____________________
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
ipCID = ipGet(nconf1*nroots)

call mma_allocate(iPre,nipre(isym),Label='iPre')
npre2 = npre(isym)
ipPre2 = ipGet(npre2)

call ipIn(ipPre2)
call Prec(W(ipPre2)%A,iPre,isym)
call ipOut(ippre2)

! OK START WORKING

!idisp = 1
jspin = 0

! Allocate areas for scratch and state variables

call mma_allocate(Kappa,nDensC,Label='Kappa')
call mma_allocate(dKappa,nDensC,Label='dKappa')
call mma_allocate(Sigma,nDensC,Label='Sigma')
call mma_allocate(Temp4,nDensC,Label='Temp4')
call mma_allocate(Temp5,nDens,Label='Temp5')
call mma_allocate(Sc1,nDens,Label='Sc1')
call mma_allocate(Sc2,nDensC,Label='Sc2')

! I think the lagrange multiplers are independent of the
! displacement, no?
nDisp = 1
do iDisp=1,nDisp
  Kappa(:) = Zero
  dKappa(:) = Zero
  Sigma(:) = Zero

  !---------------------------------------------------------------------
  !
  ! Calculate RHS for the perturbation
  !
  !---------------------------------------------------------------------

  ! (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)

  if (debug) then
    if (isNAC) then
      write(u6,*) 'States: ',NACstates(1),NACstates(2)
    else
      write(u6,*) 'State: ',irlxroot
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !AMS - I Think I can skip all of this RHS stuff - I'll read it in below.
  !if (isNAC) then
  !  call RHS_NAC(Temp5,rDum)
  !else
  !  call RHS_SA(Temp5,rDum)
  !end if

  !AMS _____________________________________________________
  ! Read in the Fock operator for the calculation of the CI part of the RHS ipF1 and ipF2.
  nOrbAct = sum(nAsh(1:nSym))
  nTri = 0
  do ksym=1,nsym
    nTri = nTri+nTri_Elem(nBas(ksym))
  end do
  nacpar = nTri_Elem(nOrbAct)
  call mma_allocate(FMO1t,nTri,Label='FMO1t')
  call mma_allocate(FMO1,nDens,Label='FMO1')
  nacpar = nTri_Elem(nnA)
  nacpr2 = nTri_Elem(nacpar)
  call mma_allocate(FMO2t,nacpr2,Label='FMO2t')

  call get_darray('F1_PDFT',FMO1t,nTri)

  ioff = 0
  do iS=1,nSym
    jS = iS
    if (nBas(is)*nBas(jS) /= 0) then
      if (iS == jS) then
        do i=1,nBas(iS)
          do j=1,i
            ioff = ioff+1
            ji = ipMat(is,js)-1+(i-1)*nbas(iS)+j
            FMO1(ji) = FMO1t(ioff)
            if (i /= j) then
              ij = ipMat(is,js)-1+(j-1)*nbas(iS)+i
              FMO1(ij) = FMO1t(ioff)
            end if
          end do
        end do
      else
        do i=1,nBas(iS)
          do j=1,nBas(jS)
            ioff = ioff+1
            ji = ipMat(is,js)-1+(i-1)*nbas(iSym)+j
            FMO1(ji) = FMO1t(ioff)
          end do
        end do
      end if
    end if
  end do

  call get_darray('F2_PDFT',FMO2t,nacpr2)

  call CISigma_sa(0,State_sym,State_sym,FMO1,nDens,FMO2t,size(FMO2t),rdum,1,ipci,ipST,.true.)
  call mma_deallocate(FMO2t)

  troot = (irlxroot-1)
  call ipin(ipST)
  call ipin(ipCI)
  do i=0,nroots-1
    if (i == troot) then
      W(ipST)%A(i*nconf1+1:(i+1)*nconf1) = W(ipST)%A(i*nconf1+1:(i+1)*nconf1)/weight(i+1)
      rE = ddot_(nconf1,W(ipST)%A(1+i*nconf1),1,W(ipCI)%A(1+i*nconf1),1)
      W(ipST)%A(i*nconf1+1:(i+1)*nconf1) = W(ipST)%A(i*nconf1+1:(i+1)*nconf1)-rE*W(ipCI)%A(i*nconf1+1:(i+1)*nconf1)

    else
      W(ipst)%A(i*nconf1+1:(i+1)*nconf1) = Zero
    end if
  end do

  W(ipST)%A(1:nconf1*nroots) = -Two*W(ipST)%A(1:nconf1*nroots)

  ! scaling the CI resp. for PDFT part in HMC-PDFT
  if (Do_Hybrid) W(ipST)%A(1:nconf1*nroots) = PDFT_Ratio*W(ipST)%A(1:nconf1*nroots)

  if (debug) then
    write(u6,*) 'RHS CI part:'
    do iS=1,nconf1*nroots
      write(u6,*) W(ipST)%A(iS)
    end do
  end if

  call mma_deallocate(FMO1t)
  call mma_deallocate(FMO1)

  ! Get the fock matrix needed for the determination of the orbital part of the RHS.

  call mma_allocate(FT99,nDens,Label='FT99')
  FT99(:) = Zero
  Temp5(:) = Zero
  call get_dArray('Fock_PDFT',FT99,nDens)
  do iS=1,nSym
    jS = Mul(iS,1)
    if (nBas(is)*nBas(jS) /= 0) &
      call DGeSub(FT99(ipMat(iS,jS)),nBas(iS),'N',FT99(ipMat(jS,iS)),nBas(jS),'T',Temp5(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
  end do
  Temp5(:) = -Two*Temp5(:)

  call mma_deallocate(FT99)
  if (Do_Hybrid) then
    ! scaling the orb resp. for PDFT part in HMC-PDFT
    Temp5(:) = PDFT_Ratio*Temp5(:)
    ! calculating the orb resp. for WF part in HMC-PDFT
    call mma_allocate(WForb,nDens,Label='WForb')
    ! saving Fock matrix for PDFT part in HMC-PDFT
    call mma_allocate(FOTr,nTri,Label='FOTr')
    call Get_dArray_chk('FockOcc',FOTr,nTri)
    ! note that the Fock matrix will be overwritten with the wf one
    ! ini rhs_sa
    call rhs_sa(WForb,rDum,ipS2)
    Temp5(:) = Temp5(:)+WF_Ratio*WForb(:)
    call mma_deallocate(WForb)
  end if

  if (debug) then
    write(u6,*) 'RHS orb part:'
    do iS=1,nDens
      write(u6,*) Temp5(iS)
    end do
  end if

  ! Also, along with this RHS stuff, the Fock_occ array already stored on
  ! the runfile needs to be replaced - switch triangular storage to square
  ! storage:

  if (Do_Hybrid) then
    ng1 = nTri_Elem(ntash)
    ng2 = nTri_Elem(ng1)
    call mma_allocate(FOSq,nDens,Label='FOSq')
    call Get_dArray_chk('FockOcc',FOsq,nDens)

    ! scaling fock for wf part
    ! adding fock for pdft part
    FOsq(1:nTri) = WF_Ratio*FOsq(1:nTri)+PDFT_Ratio*FOtr(:)

    call mma_allocate(P2PDFT,nG2,Label='P2PDFT')
    call mma_allocate(P2WF,nG2,Label='P2WF')

    call Get_dArray_chk('P2MOt',P2PDFT,nG2)
    call Get_dArray_chk('P2mo',P2WF,nG2)
    ! scaling P2 for pdft part'
    ! adding P2 for wf part'
    P2PDFT(:) = PDFT_Ratio*P2PDFT(:)+WF_Ratio*P2WF(:)

    call Put_dArray('P2MOt',P2PDFT,nG2)

    call Put_dArray('FockOcc',FOSq,nDens)

    call mma_deallocate(FOSq)
    call mma_deallocate(FOTr)
    call mma_deallocate(P2PDFT)
    call mma_deallocate(P2WF)

  else
    call mma_allocate(FOSq,nDens,Label='FOSq')
    FOSq(:) = Zero
    call Get_dArray_chk('FockOcc',FOSq,nTri)
    call Put_dArray('FockOcc',FOSq,nDens)

    call mma_deallocate(FOSq)
  end if

  ! This seems to calculate the RHS, at least for the orbital part.
  ! Now, my sigma_0 should be given by
  ! (RHS) - A*Kappa, where Kappa is my initial guess at the solution, x_0.
  ! So, should I be running a "TimesE2"-like subroutine here, to do the
  ! A*Kappa multiplication, before we go on to multiply things by the
  ! preconditioner inverse?
  !___________________________________________________________

  call opOut(ipci)

  if (lprint) write(u6,*) '       Iteration       Delta       Res(kappa)  Res(CI)     DeltaK      DeltaC'
  iLen = nDensC
  iRHSDisp(iDisp) = iDis
  call Compress(Temp5,Sigma,iSym)
  r1 = ddot_(nDensC,Sigma,1,Sigma,1)
  if (debug) write(u6,*) 'Hi how about r1',r1
  call dDaFile(LuTemp,1,Sigma,iLen,iDis)

  call ipIn(ipCIT)
  W(ipCIT)%A(1:nConf1*nroots) = Zero
  call ipIn(ipCID)
  W(ipCID)%A(1:nConf1*nroots) = Zero
  call ipOut(ipCIT)
  Sigma(:) = -Sigma(:)

  deltaC = Zero
  !AMS _________________________________________________________
  ! I need to read in the CI portion of the RHS here.
  if (CI) then
    call ipIn(ipS2)
    call DMinvCI_sa(ipST,W(ipS2)%A,Fancy)
  end if
  call ipin(ipST)
  call ipin(ipCId)
  W(ipCId)%A(1:nconf1*nroots) = W(ipST)%A(1:nconf1*nroots)
  !*******************
  !TRS
  call mma_allocate(lmroots,nroots,Label='lmroots')
  call mma_allocate(lmroots_new,nroots,Label='lmroots_new')
  call mma_allocate(kap_new,nDensC,Label='kap_new')
  call mma_allocate(kap_new_temp,nDensC,Label='kap_new_temp')

  Kap_New(:) = Zero
  Kap_New_Temp(:) = Zero

  call ipin(ipCI)
  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%A,nconf1,W(ipCId)%A(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)
  ! SA-SA rotations w/in SA space in eigen state basis
  if (debug) call recprt('lmroots',' ',lmroots,1,nroots)
  ! SA-SA rotations w/in SA space in CSF basis
  call dgemv_('N',nconf1,nroots,One,W(ipCI)%A,nconf1,lmroots,1,Zero,W(ipCId)%A(1+(irlxroot-1)*nconf1),1)
  ! SA-SA rotations w/in SA space for new lagrange multipliers
  do i=1,nroots
    if (i == irlxroot) then
      lmroots_new(i) = Zero
    else
      diff = (ERASSCF(i)-ERASSCF(irlxroot))
      wscale = (One/(Two*diff))*(One/weight(i))
      if (debug) then
        write(u6,*) 'diff',diff
        write(u6,*) 'wscale',wscale
        write(u6,*) 'weight',weight(i)
      end if
      lmroots_new(i) = wscale*lmroots(i)
    end if
  end do

  if (debug) call recprt('lmroots_new',' ',lmroots_new,1,nroots)
  ! SA-SA rotations w/in SA space for new lagrange multipliers in csf basis
  call dgemv_('N',nconf1,nroots,One,W(ipCI)%A,nconf1,lmroots_new,1,Zero,W(ipcid)%A(1+(irlxroot-1)*nconf1),1)

  ! First iter of PCG
  call TimesE2_(kap_new,ipCId,1,reco,jspin,ipS2,kap_new_temp,ipS1)

  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%A,nconf1,W(ipST)%A(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ipst this should be 1lmroots'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if

  call ipin(ipS1)
  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%A,nconf1,W(ipS1)%A(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ips1 this should be -lmroots'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if
  ! Initializing some of the elements of the PCG
  ! Modifying the response
  call ipIn(ipS1)
  call ipIn(ipST)
  W(ipST)%A(1:nConf1*nroots) = W(ipST)%A(1:nConf1*nroots)-W(ipS1)%A(1:nConf1*nroots)

  ! Kap part put into  sigma
  Sigma(:) = Sigma(:)-kap_new_temp(:)
  call ipIn(ipCId)
  call ipIn(ipCIT)
  W(ipCIT)%A(1:nConf1*nroots) = W(ipCIT)%A(1:nConf1*nroots)+W(ipCId)%A(1:nConf1*nroots)

  W(ipCId)%A(1:nconf1*nroots) = W(ipST)%A(1:nconf1*nroots)

  call opOut(ipci)
  call opOut(ipdia)

  call ipIn(ipPre2)
  call DMInvKap(W(ipPre2)%A,iPre,Sigma,dKappa,Sc1,iSym,iter)
  call opOut(ippre2)
  r2 = ddot_(nDensC,dKappa,1,dKappa,1)
  if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

  call mma_deallocate(kap_new)
  call mma_deallocate(kap_new_temp)
  call mma_deallocate(lmroots_new)
  !TRS
  !*********************
  call ipin(ipCI)
  call ipin(ipST)
  call DgeMV_('T',nconf1,nroots,One,W(ipci)%A,nconf1,W(ipST)%A(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ipst this should be zero'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if
  call mma_deallocate(lmroots)

  if (CI) then
    call ipin(ipCId)
    deltaC = ddot_(nConf1*nroots,W(ipST)%A,1,W(ipCId)%A,1)
    call ipout(ipcid)
  else
    deltaC = Zero
  end if
  !AMS_______________________________________________

  call ipOut(ipcid)
  deltaK = ddot_(nDensC,dKappa,1,Sigma,1)
  delta = deltac+deltaK
  !write(u6,*) 'deltac and deltak',deltac,deltak
  delta0 = delta
  iter = 1
  ! Naming System:
  ! Kappa: accumulates Lagrange multiplier orbital parts (dKappa * ralpha)
  ! dKappa: orbital input of Hessian matrix-vector;
  ! Temp4: orbital output of Hessian matrix-vector
  ! Sigma: accumulates error vector orbital part
  ! ipCIT: accumulates Lagrange multiplier CI parts (ipCId * ralpha)
  ! ipCId: CI input of Hessian matrix-vector;
  ! ipS1: CI output of Hessian matrix-vector
  ! ipST: accumulates error vector CI part

  !---------------------------------------------------------------------

  cnvrgd = .true.
  do
    if (delta == Zero) exit

    call TimesE2_(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

    !-------------------------------------------------------------------
    !
    !            delta
    ! rAlpha=------------
    !        dKappa:dSigma
    !
    !-------------------------------------------------------------------

    rAlphaK = Zero
    rAlphaK = ddot_(nDensC,Temp4,1,dKappa,1)
    rAlphaC = Zero
    call ipIn(ipS1)
    call ipIn(ipCId)
    rAlphaC = ddot_(nConf1*nroots,W(ipS1)%A,1,W(ipCId)%A,1)

    rAlpha = delta/(rAlphaK+rAlphaC)

    !------------------------------------------------------------------*

    ! Kappa=Kappa+rAlpha*dKappa
    Kappa(:) = Kappa(:)+ralpha*dKappa(:)
    ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
    Sigma(:) = Sigma(:)-ralpha*Temp4(:)
    resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))

    resci = Zero
    call ipIn(ipCIT)
    W(ipCIT)%A(1:nConf1*nroots) = W(ipCIT)%A(1:nConf1*nroots)+ralpha*W(ipCId)%A(1:nConf1*nroots)
    call ipOut(ipCIT)
    ! ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
    call ipIn(ipS1)
    call ipIn(ipST)
    W(ipST)%A(1:nConf1*nroots) = W(ipST)%A(1:nConf1*nroots)-ralpha*W(ipS1)%A(1:nConf1*nroots)
    call opOut(ipS1)
    resci = sqrt(ddot_(nconf1*nroots,W(ipST)%A,1,W(ipST)%A,1))

    !------------------------------------------------------------------*
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

    !------------------------------------------------------------------*
    !      s:Sigma (k+1)     s:Sigma (k+1)
    ! Beta=-------        =  -------------
    !       delta  (k)        s:Sigma (k)
    !
    ! delta=s:sigma
    !
    ! dKappa=s+Beta*dKappa

    call ipIn(ipST)
    call ipIn(ipS2)
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
      call ipIn(ipS2)
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
    !------------------------------------------------------------------*

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
    if (lprint) write(u6,Fmt2//'I7,7X,F12.7,F12.7,F12.7,F12.7,F12.7)') iter,delta/delta0,resk,resci,deltac,deltak
    iter = iter+1

  end do

  !*********************************************************************

  if (.not. cnvrgd) then
    write(u6,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
    converged(isym) = .false.
    fail = .true.
  else
    if (iPL >= 2) then
      write(u6,Fmt2//'I7,7X,F12.7,F12.7,F12.7,F12.7,F12.7)') iter,delta/delta0,resk,resci,deltac,deltak
      write(u6,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
    end if
    call ipnout(-1)
  end if

  if (iPL >= 2) write(u6,*)
  if (debug) then
    write(u6,*) 'outputs'
    write(u6,*) 'kappa'
    do i=1,nDensC
      write(u6,*) Kappa(i)
    end do
    call ipin(ipCIT)
    !W(ipCIT)%A(1:nconf1*nroots) = Zero
    write(u6,*) 'cit'
    do i=1,nconf1*nroots
      write(u6,*) W(ipCIT)%A(i)
    end do
  end if

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
end do

call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
call mma_deallocate(Temp5)
call mma_deallocate(Temp4)
call mma_deallocate(Sigma)
call mma_deallocate(dKappa)
call mma_deallocate(Kappa)
call mma_deallocate(Fancy)

! Free all memory and remove from disk all data
! related to this symmetry

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

end subroutine WfCtl_pdft
