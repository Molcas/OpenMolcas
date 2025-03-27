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

use Exp, only: Exp_Close
use ipPage, only: W
use PDFT_Util, only: Do_Hybrid, WF_Ratio, PDFT_Ratio
use MCLR_Data, only: nConf1, nDens2, nDensC, nDens, ipCI, nAcPar, nNA, nAcPr2, ipMat
use MCLR_Data, only: ipDia
use MCLR_Data, only: ISNAC, IRLXROOT, NACSTATES
use MCLR_Data, only: LuTemp
use MCLR_Data, only: XISPSM
use input_mclr, only: nDisp, Fail, lSave, nSym, State_Sym, iMethod, iBreak, Eps, nIter, Weight, Debug, ERASSCF, kPrint, nCSF, &
                      nRoots, ntAsh, nAsh, nBas, nRs2
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
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
integer opOut
logical lPrint, cnvrgd
real*8 rchc(mxroot)
real*8 rDum(1)
external IsFreeUnit
real*8, allocatable :: FOSq(:), FOTr(:)
real*8, allocatable :: Kappa(:), dKappa(:), Sigma(:), Temp4(:), Sc1(:), Sc2(:), Fancy(:), FMO1t(:), FMO1(:), FMO2t(:), FT99(:), &
                       Temp5(:)
real*8, allocatable :: lmroots(:), lmroots_new(:), Kap_New(:), Kap_New_Temp(:)
real*8, allocatable :: WFOrb(:), P2WF(:), P2PDFT(:)
real*8 R1, R2, DeltaC, DeltaK, Delta, Delta0, ReCo, rAlphaC, rAlphaK, rAlpha, rEsk, rEsci, rBeta, Res, TRoot, rE, Diff, WScale
real*8, external :: DDot_
integer lPaper, lLine, Left, iDis, Lu_50, iDisp, iSym, nConf3, iRC, ipS1, ipS2, ipST, ipCIT, ipCID, nPre2, iLen, Iter, ipPre2, &
        jSpin, i, nTri, nOrbAct, kSym, iOff, iS, jS, j, ji, ij, nG1, nG2
integer, external :: ipClose, ipGet, ipIn, ipOut, ipNOut
integer, external :: nPre
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    integer iispin, iCsym, iSSym
    integer nInt1, nInt2s, nInt2a
    real*8, target :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    integer ipCI1, ipCI2
    logical Have_2_el
  end subroutine CISigma_sa
  subroutine rhs_sa(Fock,SLag_pt2)
    real*8 Fock(*)
    real*8, optional :: SLag_pt2(*)
  end subroutine
end interface

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

! Determine if we should page CI vectors
!                            [2]
! Calculate the diagonal of E    and store in core/disc

call mma_allocate(Fancy,nRoots**3,Label='Fancy')
!____________________
!AMS - what to do here?
! What should rCHC be?  is it computed with E(mcscf) or E(pdft)?
!____________________
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
ipCID = ipGet(nconf1*nroots)

npre2 = npre(isym)
ipPre2 = ipGet(npre2)

irc = ipIn(ipPre2)
call Prec(W(ipPre2)%Vec,isym)
irc = ipOut(ippre2)

! OK START WORKING

!idisp = 1
jspin = 0

! Allocate areas for scratch and state variables

call mma_allocate(Kappa,nDens2+6,Label='Kappa')
call mma_allocate(dKappa,nDens2+6,Label='dKappa')
call mma_allocate(Sigma,nDens2+6,Label='Sigma')
call mma_allocate(Temp4,nDens2+6,Label='Temp4')
call mma_allocate(Sc1,nDens2+6,Label='Sc1')
call mma_allocate(Sc2,nDens2+6,Label='Sc2')

! I think the lagrange multiplers are independent of the
! displacement, no?
nDisp = 1
do iDisp=1,nDisp
  Kappa(1:nDens2) = Zero
  dKappa(1:nDens2) = Zero
  Sigma(1:nDens2) = Zero

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
  !  call RHS_NAC(Temp4)
  !else
  !  call RHS_SA(Temp4)
  !end if

  !AMS _____________________________________________________
  ! Read in the Fock operator for the calculation of the CI part of the RHS ipF1 and ipF2.
  nTri = 0
  nOrbAct = 0
  do ksym=1,nsym
    nTri = nTri+nBas(ksym)*(nBas(ksym)+1)/2
    nOrbAct = nOrbAct+nAsh(ksym)
  end do
  nacpar = nOrbAct*(nOrbAct+1)/2
  call mma_allocate(FMO1t,nTri,Label='FMO1t')
  call mma_allocate(FMO1,nDens2,Label='FMO1')
  nacpar = (nnA+1)*nnA/2
  nacpr2 = (nacpar+1)*nacpar/2
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

  call CISigma_sa(0,State_sym,State_sym,FMO1,nDens2,FMO2t,size(FMO2t),rdum,1,ipci,ipST,.true.)
  call mma_deallocate(FMO2t)

  troot = (irlxroot-1)
  irc = ipin(ipST)
  irc = ipin(ipCI)
  do i=0,nroots-1
    if (i == troot) then
      call Dscal_(nconf1,(1/weight(i+1)),W(ipST)%Vec(1+i*nconf1),1)
      rE = ddot_(nconf1,W(ipST)%Vec(1+i*nconf1),1,W(ipCI)%Vec(1+i*nconf1),1)
      call Daxpy_(nconf1,-rE,W(ipCI)%Vec(1+i*nconf1),1,W(ipST)%Vec(1+i*nconf1),1)

    else
      call dcopy_(nConf1,[Zero],0,W(ipst)%Vec(1+i*nconf1),1)
    end if
  end do

  call DSCAL_(nconf1*nroots,-Two,W(ipST)%Vec,1)

  ! scaling the CI resp. for PDFT part in HMC-PDFT
  if (Do_Hybrid) call DScal_(nconf1*nroots,PDFT_Ratio,W(ipST)%Vec,1)

  if (debug) then
    write(u6,*) 'RHS CI part:'
    do iS=1,nconf1*nroots
      write(u6,*) W(ipST)%Vec(iS)
    end do
  end if

  call mma_deallocate(FMO1t)
  call mma_deallocate(FMO1)

  ! Get the fock matrix needed for the determination of the orbital part of the RHS.

  call mma_allocate(FT99,nDens2,Label='FT99')
  call mma_allocate(Temp5,nDens2+6,Label='Temp5')
  FT99(:) = Zero
  Temp5(:) = Zero
  call get_dArray('Fock_PDFT',FT99,nDens2)
  do iS=1,nSym
    jS = ieor(iS-1,0)+1
    if (nBas(is)*nBas(jS) /= 0) &
      call DGeSub(FT99(ipMat(iS,jS)),nBas(iS),'N',FT99(ipMat(jS,iS)),nBas(jS),'T',Temp5(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
  end do
  call dcopy_(nDens2+6,Temp5,1,Temp4,1)
  call DSCAL_(ndens2+6,-Two,Temp4,1)

  call mma_deallocate(FT99)
  call mma_deallocate(Temp5)
  if (Do_Hybrid) then
    ! scaling the orb resp. for PDFT part in HMC-PDFT
    call DSCAL_(ndens2,PDFT_Ratio,Temp4,1)
    ! calculating the orb resp. for WF part in HMC-PDFT
    call mma_allocate(WForb,nDens2+6,Label='WForb')
    ! saving Fock matrix for PDFT part in HMC-PDFT
    call mma_allocate(FOTr,nTri,Label='FOTr')
    call Get_dArray_chk('FockOcc',FOTr,nTri)
    ! note that the Fock matrix will be overwritten with the wf one
    ! ini rhs_sa
    call rhs_sa(WForb)
    call dAXpY_(nDens2,WF_Ratio,WForb,1,Temp4,1)
    call mma_deallocate(WForb)
  end if

  if (debug) then
    write(u6,*) 'RHS orb part:'
    do iS=1,nDens2
      write(u6,*) Temp4(iS)
    end do
  end if

  ! Also, along with this RHS stuff, the Fock_occ array already stored on
  ! the runfile needs to be replaced - switch triangular storage to square
  ! storage:

  if (Do_Hybrid) then
    ng1 = (ntash+1)*ntash/2
    ng2 = (ng1+1)*ng1/2
    call mma_allocate(FOSq,nDens2,Label='FOSq')
    call Get_dArray_chk('FockOcc',FOsq,nDens2)

    ! scaling fock for wf part
    call DScal_(nTri,WF_Ratio,FOsq,1)

    ! adding fock for pdft part
    call daxpy_(ntri,pdft_ratio,fotr,1,fosq,1)

    call mma_allocate(P2PDFT,nG2,Label='P2PDFT')
    call mma_allocate(P2WF,nG2,Label='P2WF')

    call Get_dArray_chk('P2MOt',P2PDFT,nG2)
    ! scaling P2 for pdft part'
    call DScal_(nG2,PDFT_Ratio,P2PDFT,1)

    call Get_dArray_chk('P2mo',P2WF,nG2)
    ! adding P2 for wf part'
    call daxpy_(ng2,wf_ratio,P2WF,1,P2PDFT,1)

    call Put_dArray('P2MOt',P2PDFT,nG2)

    call Put_dArray('FockOcc',FOSq,ndens2)

    call mma_deallocate(FOSq)
    call mma_deallocate(FOTr)
    call mma_deallocate(P2PDFT)
    call mma_deallocate(P2WF)

  else
    call mma_allocate(FOSq,nDens2,Label='FOSq')
    call mma_allocate(FOTr,nTri,Label='FOTr')
    FOSq(:) = Zero
    call Get_dArray_chk('FockOcc',FOTr,nTri)
    call dcopy_(nTri,FOtr,1,FOSq,1)
    call Put_dArray('FockOcc',FOSq,ndens2)

    call mma_deallocate(FOSq)
    call mma_deallocate(FOTr)
  end if

  ! This seems to calculate the RHS, at least for the orbital part.
  ! Now, my sigma_0 should be given by
  ! (RHS) - A*Kappa, where Kappa is my initial guess at the solution, x_0.
  ! So, should I be running a "TimesE2"-like subroutine here, to do the
  ! A*Kappa multiplication, before we go on to multiply things by the
  ! preconditioner inverse?
  !___________________________________________________________

  irc = opOut(ipci)

  if (lprint) write(u6,*) '       Iteration       Delta       Res(kappa)  Res(CI)     DeltaK      DeltaC'
  iLen = nDensC
  iRHSDisp(iDisp) = iDis
  do iS=1,nDens2
  end do
  call Compress(Temp4,Sigma,iSym)
  r1 = ddot_(nDensc,Sigma,1,Sigma,1)
  if (debug) write(u6,*) 'Hi how about r1',r1
  call dDaFile(LuTemp,1,Sigma,iLen,iDis)

  irc = ipIn(ipCIT)
  call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
  irc = ipIn(ipCID)
  call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
  irc = ipOut(ipCIT)
  call DSCAL_(nDensC,-One,Sigma,1)

  deltaC = Zero
  !AMS _________________________________________________________
  ! I need to read in the CI portion of the RHS here.
  if (CI) then
    irc = ipIn(ipS2)
    call DMinvCI_sa(ipST,W(ipS2)%Vec,Fancy)
  end if
  irc = ipin(ipST)
  irc = ipin(ipCId)
  call dcopy_(nconf1*nroots,W(ipST)%Vec,1,W(ipCId)%Vec,1)
  !*******************
  !TRS
  call mma_allocate(lmroots,nroots,Label='lmroots')
  call mma_allocate(lmroots_new,nroots,Label='lmroots_new')
  call mma_allocate(kap_new,ndensc,Label='kap_new')
  call mma_allocate(kap_new_temp,ndens,Label='kap_new_temp')

  Kap_New(:) = Zero
  Kap_New_Temp(:) = Zero

  irc = ipin(ipCI)
  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%Vec,nconf1,W(ipCId)%Vec(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)
  ! SA-SA rotations w/in SA space in eigen state basis
  if (debug) call recprt('lmroots',' ',lmroots,1,nroots)
  ! SA-SA rotations w/in SA space in CSF basis
  call dgemv_('N',nconf1,nroots,One,W(ipCI)%Vec,nconf1,lmroots,1,Zero,W(ipCId)%Vec(1+(irlxroot-1)*nconf1),1)
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
  call dgemv_('N',nconf1,nroots,One,W(ipCI)%Vec,nconf1,lmroots_new,1,Zero,W(ipcid)%Vec(1+(irlxroot-1)*nconf1),1)

  ! First iter of PCG
  call TimesE2_(kap_new,ipCId,1,reco,jspin,ipS2,kap_new_temp,ipS1)

  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%Vec,nconf1,W(ipST)%Vec(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ipst this should be 1lmroots'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if

  irc = ipin(ipS1)
  call DgeMV_('T',nconf1,nroots,One,W(ipCI)%Vec,nconf1,W(ipS1)%Vec(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ips1 this should be -lmroots'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if
  ! Initializing some of the elements of the PCG
  ! Modifying the response
  irc = ipIn(ipS1)
  irc = ipIn(ipST)
  call DaXpY_(nConf1*nroots,-One,W(ipS1)%Vec,1,W(ipST)%Vec,1)

  ! Kap part put into  sigma
  call DaxPy_(nDensC,-One,kap_new_temp,1,Sigma,1)
  irc = ipIn(ipCId)
  irc = ipIn(ipCIT)
  call DaXpY_(nConf1*nroots,One,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)

  call dcopy_(nconf1*nroots,W(ipST)%Vec,1,W(ipCId)%Vec,1)

  irc = opOut(ipci)
  irc = opOut(ipdia)

  irc = ipIn(ipPre2)
  call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,dKappa,nDens2+6,Sc1,nDens2+6,iSym,iter)
  irc = opOut(ippre2)
  r2 = ddot_(ndensc,dKappa,1,dKappa,1)
  if (r2 > r1) write(u6,*) 'Warning perturbation number ',idisp,' might diverge'

  call mma_deallocate(kap_new)
  call mma_deallocate(kap_new_temp)
  call mma_deallocate(lmroots_new)
  !TRS
  !*********************
  irc = ipin(ipCI)
  irc = ipin(ipST)
  call DgeMV_('T',nconf1,nroots,One,W(ipci)%Vec,nconf1,W(ipST)%Vec(1+(irlxroot-1)*nconf1),1,Zero,lmroots,1)

  if (debug) then
    write(u6,*) 'lmroots_ipst this should be zero'
    call recprt('lmroots',' ',lmroots,1,nroots)
  end if
  call mma_deallocate(lmroots)

  if (CI) then
    irc = ipin(ipCId)
    deltaC = ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipCId)%Vec,1)
    irc = ipout(ipcid)
  else
    deltaC = Zero
  end if
  !AMS_______________________________________________

  irc = ipOut(ipcid)
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
    irc = ipIn(ipS1)
    irc = ipIn(ipCId)
    rAlphaC = ddot_(nConf1*nroots,W(ipS1)%Vec,1,W(ipCId)%Vec,1)

    rAlpha = delta/(rAlphaK+rAlphaC)

    !------------------------------------------------------------------*

    ! Kappa=Kappa+rAlpha*dKappa
    call DaxPy_(nDensC,ralpha,dKappa,1,Kappa,1)
    ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
    call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
    resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))

    resci = Zero
    irc = ipIn(ipCIT)
    call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
    irc = ipOut(ipCIT)
    ! ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
    irc = ipIn(ipS1)
    irc = ipIn(ipST)
    call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
    irc = opOut(ipS1)
    resci = sqrt(ddot_(nconf1*nroots,W(ipST)%Vec,1,W(ipST)%Vec,1))

    !------------------------------------------------------------------*
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

    !------------------------------------------------------------------*
    !      s:Sigma (k+1)     s:Sigma (k+1)
    ! Beta=-------        =  -------------
    !       delta  (k)        s:Sigma (k)
    !
    ! delta=s:sigma
    !
    ! dKappa=s+Beta*dKappa

    irc = ipIn(ipST)
    irc = ipIn(ipS2)
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
      irc = ipIn(ipS2)
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
    irc = ipnout(-1)
  end if

  if (iPL >= 2) write(u6,*)
  if (debug) then
    write(u6,*) 'outputs'
    write(u6,*) 'kappa'
    do i=1,ndens2
      write(u6,*) Kappa(i)
    end do
    irc = ipin(ipCIT)
    !call dcopy_(nconf1*nroots,Zero,0,W(ipCIT)%Vec,1)
    write(u6,*) 'cit'
    do i=1,nconf1*nroots
      write(u6,*) W(ipCIT)%Vec(i)
    end do
  end if

  iLen = ndensC

  iKapDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,Kappa,iLen,iDis)
  iSigDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,Sigma,iLen,iDis)
  ilen = nconf1*nroots
  iCIDisp(iDisp) = iDis

  irc = ipin(ipCIT)
  call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)

  !MGD This last call seems unused, so I comment it

  !call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,Temp4,ipS2)
  iCISigDisp(iDisp) = iDis
  irc = ipin(ipST)
  call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
end do

call mma_deallocate(Sc2)
call mma_deallocate(Sc1)
call mma_deallocate(Temp4)
call mma_deallocate(Sigma)
call mma_deallocate(dKappa)
call mma_deallocate(Kappa)
call mma_deallocate(Fancy)

! Free all memory and remove from disk all data
! related to this symmetry

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

end subroutine WfCtl_pdft
