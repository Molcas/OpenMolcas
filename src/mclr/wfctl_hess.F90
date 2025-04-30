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
!               2002, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine WfCtl_Hess(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,iRHSDisp,iRHSCIDISP,converged)
!***********************************************************************
!                                                                      *
!     called from: MCLR                                                *
!                                                                      *
!----------------------------------------------------------------------*
!     Parallelization of perturbations, RL 2002.                       *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use ipPage, only: ipclose, ipget, ipin, ipin1, ipnout, ipout, opout, W
use Para_Info, only: myRank, nProcs
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Spool, only: LuWr
use MCLR_Data, only: CMO, FIMO, Int2, ipCI, ipDia, lDisp, LuTemp, n1Dens, n2Dens, nConf1, nDens, nDensC, XISPSM
use MCLR_procedures, only: CISigma
use input_mclr, only: Eps, ERASSCF, Fail, iBreak, iMethod, kPrint, lSave, nAsh, nCSF, nDisp, nIter, nRs2, nSym, nTPert, PotNuc, &
                      PT2, rIn_Ene, State_Sym, TimeDep
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use MCLR_Data, only: nCMO
#endif

implicit none
integer(kind=iwp), intent(out) :: iKapDisp(nDisp), isigDisp(nDisp), iCIDisp(nDisp), iCIsigDisp(nDisp), iRHSDisp(nDisp), &
                                  iRHSCIDisp(nDisp)
logical(kind=iwp), intent(out) :: converged(8)
integer(kind=iwp) :: iD, iDis, iDisp, iLen, ipCID, ipCIT, ipPre2, ipS1, ipS2, ipST, iRank, istatus, iSym, iSym_Old, Iter, jDisp, &
                     jSpin, kkkSym, kkSym, Left, lLine, lPaper, LuWR_Save, nConf3, nPre2, pstate_sym
real(kind=wp) :: Clock(4), D_0, Delta, Delta0, DeltaC, DeltaK, EC, R1, R2, rAlpha, rAlphaC, rAlphaK, rBeta, rCHC, rdum(1), ReCo, &
                 Res, res_tmp, rEsci, rEsk, rGrad, Tim2, Tim3, Tim4
logical(kind=iwp) :: CI, cnvrgd, lPrint, Orb, Response
character(len=132) :: Line
character(len=72) :: SLine
character(len=8) :: Fmt2
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp) :: iglfail
real(kind=wp) :: dfail
#endif
integer(kind=iwp), allocatable :: iPre(:), List(:,:)
real(kind=wp), allocatable :: Dens(:), dKappa(:), Kappa(:), Pens(:), rmoaa(:), Sc1(:), Sc2(:), Sc3(:), Sc4(:), Sigma(:), Temp1(:), &
                              Temp2(:), Temp3(:), Temp4(:)
integer(kind=iwp), parameter :: iTimeCC = 1, iTimeKK = 2, iTimeKC = 3, iTimeCK = 4
integer(kind=iwp), external :: IsFreeUnit, niPre, nPre
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: Rsv_Tsk

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
SLine = 'Solving CP(CAS)HF equations'
call StatusLine('MCLR: ',SLine)

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
!if (lSAVE) call DANAME(50,'RESIDUALS')
if (lSAVE) then
  write(LuWr,*) 'WfCtl: SAVE option not implemented'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Start loop over the symmetry of the perturbations

if (btest(kprint,1)) lprint = .true.
#ifdef _DEBUGPRINT_
lprint = .true.
#endif
kksym = 1
kkksym = nsym
if (PT2) kkkSym = 1

! Set up parallelization over the loop over perturbations

call mma_allocate(List,2,nDisp,Label='List')
iDisp = 0
do iSym=kksym,kkksym
  List(1,iDisp+1:iDisp+lDisp(iSym)) = iSym
  List(2,iDisp+1:iDisp+lDisp(iSym)) = [(jDisp,jDisp=1,lDisp(iSym))]
  iDisp = iDisp+lDisp(iSym)
end do

! Change output unit

LuWr_save = LuWr
if (MyRank /= 0) then
  LuWr = 55
  LuWr = isFreeUnit(LuWr)
  call molcas_open(luwr,'Temp_OutPut')
end if

call Init_Tsk(id,nDisp)
#ifdef _MOLCAS_MPP_
! iglfail is global "array", a flag used to communicate to all processes
! if any of them failed, the communication need not be synchronous,
if (Is_Real_Par()) then
  if (.not. GA_Create(MT_DBL,1,1,'GlFail',0,0,iglfail)) call SysAbendMsg('wfctl_hess','failed to create global failed flag',' ')

  call GA_Zero(iglfail)
  dfail = Zero
end if
#endif

ipdia = 0
ipPre2 = 0
iSym_Old = 0
cnvrgd = .true.
do
  if (.not. Rsv_Tsk(id,iDisp)) exit
# ifdef _MOLCAS_MPP_
  ! Check if some other process failed
  if (Is_Real_Par()) then
    call GA_Get(iglfail,1,1,1,1,dfail,1)
    if (dfail > Zero) then
      fail = .true.
      exit
    end if
  end if
# endif
  iSym = List(1,iDisp)
  jDisp = List(2,iDisp)

  write(SLine,'(A,I3,A)') 'Solving CP(CAS)HF equations for perturbation ',iDisp,'.'
  call StatusLine('MCLR: ',SLine)

  !do iSym=kksym,kkksym

  ! Execute setup for symmetry block if a new one!

  if (iSym /= iSym_Old) then

    if (iSym_Old /= 0) then

      ! If a previous symmetry block was processed
      ! free all memory and remove from disk all data
      ! related to this symmetry

      call mma_deallocate(iPre)
      if (CI) then
        call ipclose(ipdia)
      else
        call ipclose(ipPre2)
      end if

      call Exp_Close()

    end if
    iSym_Old = iSym

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

    if (doDMRG) then  ! yma
      nash(:) = RGras2(:)
      nrs2(:) = RGras2(:)
    end if

    call Setup_MCLR(iSym)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !     Determine if we should page CI vectors

    !                                  [2]
    !       Calculate the diagonal of E    and store in core/disc

    if (CI) then
      call CIDia_MCLR(PState_Sym,rCHC)
      call ipout(ipdia)

      ! Allocate disk/memory space
      !
      ! This areas should be addressed through ipin
      ! ipout will page them out to disk and free the memory area
      ! if page mode is used
      !
      ! opout will release the memory area without update the disk

      ips1 = ipget(nconf3)
      ips2 = ipget(nconf3)
      ipst = ipget(nconf3)
      ipcit = ipget(nconf1)
      ipcid = ipget(nconf1)

    else
      ips1 = 0
      ips2 = 0
      ipst = 0
      ipcit = 0
      ipcid = 0
    end if

    call mma_allocate(iPre,nipre(isym),Label='iPre')
    npre2 = max(npre(isym),1)
    ! "max" just there to make sure that
    ! W(ipPre2) is allocated even if
    ! npre2(isym) is zero.
    ipPre2 = ipget(npre2)
    call ipin(ipPre2)
    call Prec(W(ipPre2)%A,iPre,isym)
    call ipout(ippre2)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! OK START WORKING, looping over the perturbations of the
    ! current symmetry

    !iDEND = lDisp(iSym)
    !if ((SewLab == 'NONE') .and. (.not. mckinley)) iDEND = 1

  end if

  !do jDisp=1,iDEnd
  !  iDisp = iDisp+1
  jspin = 0
  if (btest(nTPert(idisp),0)) jSpin = 1
  if (jspin == 0) then
    nconf1 = ncsf(PState_Sym)
  else
    nConf1 = nint(xispsm(Pstate_Sym,1))
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

    Dens(:) = Zero
    Pens(:) = Zero
  end if
  if (iMethod == 2) then
    call mma_allocate(rmoaa,n2Dens,Label='rmoaa')
  else
    call mma_allocate(rmoaa,1,Label='rmoaa')
  end if
  rmoaa(:) = Zero
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Calculate RHS for the perturbation
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)

  !if (PT2) then
  !  call RHS_PT2(Temp4,Temp4,Temp4)
  !else
  call RHS(Sc1,Temp1,Temp3,Sc2,dKappa,Sc3,Temp4,ipST,iDisp,iSym-1,CMO,jdisp,CI)
# ifdef _DEBUGPRINT_
  write(LuWr,*) 'After RHS'
  write(LuWr,*) 'Sc1=',DDot_(nDens,Sc1,1,Sc1,1)
  write(LuWr,*) 'Kappa=',DDot_(nDensC,Kappa,1,Kappa,1)
  write(LuWr,*) 'Sc2=',DDot_(nDens,Sc2,1,Sc2,1)
  write(LuWr,*) 'dKap=',DDot_(nDens,dKappa,1,dKappa,1)
  write(LuWr,*) 'CMO=',DDot_(nCMO,CMO,1,CMO,1)
# endif
  !end if
  call opout(ipci)
  call mma_deallocate(Temp1)
  call mma_allocate(Temp1,nDensC,Label='Temp1')

  write(LuWr,*) 'Process perturbation number ',iDisp
  if (lprint) write(LuWr,*) '       Iteration         Delta     Res(kappa) Res(CI)'
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write RHS to disk
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iLen = nDensC
  iRHSDisp(iDisp) = iDis
  call Compress(Temp4,Sigma,iSym)
  r1 = ddot_(nDensC,Sigma,1,Sigma,1)
  call UnCompress(Sigma,Temp4,iSym)
  call dDaFile(LuTemp,1,Sigma,iLen,iDis)
  if (CI) then
    call ipin(ipCIT)
    W(ipCIT)%A(1:nConf1) = Zero
  end if
  call ipout(ipCIT)
  if (CI) then
    ilen = nconf1
    iRHSCIDisp(iDisp) = iDis
    call ipin(ipST)
    call dDaFile(LuTemp,1,W(ipST)%A,iLen,iDis)
    W(ipST)%A(1:nConf1) = -W(ipST)%A(1:nConf1)
  end if

  Sigma(:) = -Sigma(:)

# ifdef _DEBUGPRINT_
  call ipin(ipST)
  write(LuWr,*) 'ST=',DDot_(nConf1,W(ipST)%A,1,W(ipST)%A,1)
  write(LuWr,*) 'Sigma=',DDot_(nDensC,Sigma,1,Sigma,1)
  write(LuWr,*) 'Kappa=',DDot_(nDensC,Kappa,1,Kappa,1)
  write(LuWr,*) 'End RHS!'
# endif

  iter = 1
  call ipin(ipPre2)
  call DMInvKap(W(ipPre2)%A,iPre,Sigma,Kappa,Temp3,isym,iter)
  call opout(ippre2)
  r2 = ddot_(nDensC,Kappa,1,Kappa,1)
# ifdef _DEBUGPRINT_
  write(LuWr,*) 'DMinvKap'
  write(LuWr,*) 'Kap=',DDot_(nDensC,Kappa,1,Kappa,1)
  write(LuWr,*) 'End DMinvKap'
# endif
  if (r2 > r1) write(LuWr,Fmt2//'A,I3,A)') 'Warning  perturbation number ',idisp,' might diverge!'
  call UnCompress(Kappa,dKappa,iSym)

  if (CI) then
    call ipin(ipCId)
    call DMinvCI(ipST,W(ipCId)%A,rCHC,isym)
    call ipin(ipST)
    deltaC = ddot_(nConf1,W(ipST)%A,1,W(ipCId)%A,1)
    call ipout(ipcid)
  else
    deltaC = Zero
  end if
  deltaK = ddot_(nDensC,Kappa,1,Sigma,1)
  Kappa(:) = Zero
  delta = deltac+deltaK
# ifdef _DEBUGPRINT_
  if (abs(DeltaC) < 1.0e-12_wp) DeltaC = Zero
  write(LuWr,*) 'DeltaK, DeltaC, Delta=',DeltaK,DeltaC,Delta
# endif

  if (delta /= Zero) then
    delta0 = delta
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'Delta0=',Delta0
    write(LuWr,*) 'Start ITERATIONS'
    write(LuWr,*) 'Orb,CI=',ORB,CI
#   endif
    Orb = .true.
    TimeDep = .false.
    ReCo = -One
  end if
  !                                                                    *
  !*********************************************************************
  !          I   T   E   R   A   T   I   O   N   S                     *
  !*********************************************************************
  !                                                                    *
  do
    if (delta == Zero) exit

    if (doDMRG) then  ! It can be deleted
      nash(:) = RGras2(:)
      nrs2(:) = RGras2(:)
    end if

    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    !       O R B I T A L    P A R T
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    if (orb) then
      !                                                                *
      !*****************************************************************
      !                                                                *
      !           ~    ~
      ! Construct F,(ij|kl)
      !                                                                *
      !*****************************************************************
      !                                                                *

      call ipnout(-1)
      call RInt_generic(dKappa,rmoaa,rdum,Sc2,Temp3,Temp4,Sc3,isym,reco,jspin)
      Clock(iTimeKK) = Clock(iTimeKK)+Tim2

      !                                                                *
      !*****************************************************************
      !                                                                *
      ! kappa->CI
      !
      ! H(kappa)|0>
      !
      !     [2]
      ! S1=E   k (kappa TO CI)
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (CI) then
        call CISigma(jspin,State_Sym,pstate_sym,Temp4,nDens,rmoaa,size(rmoaa),rdum,1,ipCI,ipS1,.true.)
        Clock(iTimeKC) = Clock(iTimeKC)+Tim3
#       ifdef _DEBUGPRINT_
        call ipin(ipCI)
        call ipin(ipS1)
        write(LuWr,*) 'CISigma'
        call RecPrt('CI','(3F10.4)',W(ipCI)%A,1,nConf1)
        call RecPrt('S1','(3F10.4)',W(ipS1)%A,1,nConf1)
        write(LuWr,*) 'CI=',DDot_(nConf1,W(ipCI)%A,1,W(ipCI)%A,1)
        write(LuWr,*) 'S1=',DDot_(nConf1,W(ipS1)%A,1,W(ipS1)%A,1)
        write(LuWr,*) 'End CISigma'
#       endif

        ! This will give us a better
        ! convergence in the PCG. Notice that
        !
        ! ~Inactive     ~
        ! E        + <0|H|0> = 0
        !
        ! when the wavefunction is converged.

        call ipin(ipS1)
        if (isym == 1) then
          call ipin(ipCI)
          rGrad = DDot_(nconf1,W(ipCI)%A,1,W(ipS1)%A,1)
          W(ipS1)%A(1:nConf1) = W(ipS1)%A(1:nConf1)-rgrad*W(ipCI)%A(1:nConf1)
        end if
        W(ipS1)%A(1:nConf1) = Two*W(ipS1)%A(1:nConf1)

        call opout(ipCI)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    end if
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    !       C I    P A R T
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    if (CI) then
      !                                                                *
      !*****************************************************************
      !                                                                *
      !       [2]
      ! S2 = E   CID  (CI TO CI) = <i|H|d> -E<i|d>
      !                                                                *
      !*****************************************************************
      !                                                                *
      call ipnout(-1)
      call CISigma(0,PState_Sym,Pstate_sym,FIMO,size(FIMO),Int2,size(Int2),rdum,1,ipCId,ipS2,.true.)
      EC = rin_ene+potnuc-ERASSCF(1)

      call ipin(ipCId)
      call ipin(ipS2)
      W(ipS2)%A(1:nConf1) = Two*(W(ipS2)%A(1:nConf1)+EC*W(ipCId)%A(1:nConf1))
      Clock(iTimeCC) = Clock(iTimeCC)+Tim4
      call ipout(ipS2)
      call opout(ipCId)
      call opout(ipCI)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! CI -> Kappa
      !
      !      [2]
      ! SC3=E   CID   (SC3=F(<d|E|0>+<0|E|d>)
      !                                                                *
      !*****************************************************************
      !                                                                *
      Response = .true.
      call ipnout(-1)
#     ifdef _DEBUGPRINT_
      call ipin(ipCI)
      call ipin(ipCId)
      call RecPrt('CI','(3F10.4)',W(ipCI)%A,1,nConf1)
      call RecPrt('CId','(3F10.4)',W(ipCId)%A,1,nConf1)
#     endif
      call CIDens(Response,ipCI,ipCId,State_sym,PState_Sym,Pens,Dens)     ! Jeppe's

#     ifdef _DEBUGPRINT_
      write(LuWr,*) 'After CIDens'
      write(LuWr,*) 'P=',DDot_(n2Dens,Pens,1,Pens,1)
      write(LuWr,*) 'De=',DDot_(n1Dens,Dens,1,Dens,1)
#     endif

      ! density for inactive= 2(<d|0>+<0|d>)

      d_0 = Zero

      ! This is just for debugging purpose.
      ! When we use it for actual calculations d_0 == 0

      if (isym == 1) then
        call ipin(ipCI)
        call ipin(ipCid)
        d_0 = ddot_(nconf1,W(ipCid)%A,1,W(ipCI)%A,1)
      end if
      if (Response) d_0 = d_0*Two

      call FockGen(d_0,Dens,Pens,Sc1,Sc3,isym)    ! Made
      !                                                                *
      !*****************************************************************
      !                                                                *
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Sc1  kappa-> kappa
    ! Sc3  CI -> kappa
    ! S1   kappa -> CI
    ! S2   CI -> CI
    ! dKap present step
    ! Kap  kappaX
    ! CIT  CIX
    ! CId  present step
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Add together
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'Add together!'
    write(LuWr,*) 'dKap=',DDot_(nDens,dKappa,1,dKappa,1)
    write(LuWr,*) 'Sc2=',DDot_(nDens,Sc2,1,Sc2,1)
    if (CI) then
      write(LuWr,*) 'Sc3=',DDot_(nDens,Sc3,1,Sc3,1)
      call ipin1(ipS2,nconf1)
      write(LuWr,*) 'S2=',DDot_(nConf1,W(ipS2)%A,1,W(ipS2)%A,1)
      call ipin1(ipS1,nconf1)
      write(LuWr,*) 'S1=',DDot_(nConf1,W(ipS1)%A,1,W(ipS1)%A,1)
    end if
    write(LuWr,*)
#   endif
    call ipnout(-1)
    if (CI) then
      Sc1(:) = Sc2(:)+Sc3(:)
    else
      Sc1(:) = Sc2(:)
    end if
    call Compress(Sc1,Temp1,isym)    ! ds
    call Compress(dKappa,Temp2,isym) ! DX
    if (CI) then
      call ipin1(ipS1,nconf1)
      call ipin1(ipS2,nconf1)
      W(ipS1)%A(1:nConf1) = W(ipS1)%A(1:nConf1)+W(ipS2)%A(1:nConf1)
      call opout(ipS2)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !                ######   #####   #####
    !                #     # #     # #     #
    !                #     # #       #
    !                ######  #       #  ####
    !                #       #       #     #
    !                #       #     # #     #
    !                #        #####   #####
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !            delta
    ! rAlpha=------------
    !        dKappa:dSigma

    rAlphaC = Zero
    rAlphaK = Zero
    call ipin(ipS1)
    call ipin(ipCId)
    if (orb) rAlphaK = DDot_(nDensC,Temp1,1,Temp2,1)
    if (CI) rAlphaC = DDot_(nConf1,W(ipS1)%A,1,W(ipCId)%A,1)
    rAlpha = delta/(rAlphaK+rAlphaC)
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'At PCG'
    write(LuWr,*) 'S1=',DDot_(nConf1,W(ipS1)%A,1,W(ipS1)%A,1)
    write(LuWr,*) 'CId=',DDot_(nConf1,W(ipCId)%A,1,W(ipCId)%A,1)
    write(LuWr,*) 'rAlphaK, rAlphaC, rAlpha=',rAlphaK,rAlphaC,rAlpha
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Kappa=Kappa+rAlpha*dKappa
    ! Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa

    if (orb) then
      Kappa(:) = Kappa(:)+ralpha*Temp2(:)
      Sigma(:) = Sigma(:)-ralpha*Temp1(:)
      resk = sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
    end if
    resci = Zero
    if (CI) then
      call ipin(ipCIT)
      W(ipCIT)%A(1:nConf1) = W(ipCIT)%A(1:nConf1)+ralpha*W(ipCId)%A(1:nConf1)
      call ipout(ipCIT)
      call ipin1(ipST,nconf1)
      call ipin(ipS1)
      W(ipST)%A(1:nConf1) = W(ipST)%A(1:nConf1)-ralpha*W(ipS1)%A(1:nConf1)
      call opout(ipS1)
      call ipin(ipST)
      resci = sqrt(ddot_(nconf1,W(ipST)%A,1,W(ipST)%A,1))
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Precondition......
    !    -1
    ! S=M  Sigma

    call opout(ipCID)
    call ipin(ipS2)
    if (CI) call DMinvCI(ipST,W(ipS2)%A,rCHC,isym)
    call opout(ipCI)
    call opout(ipdia)

    call ipin(ipPre2)
    call DMInvKap(W(ipPre2)%A,iPre,Sigma,Sc4,Sc1,iSym,iter)
    call opout(ippre2)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !      s:Sigma
    ! Beta=-------
    !       delta
    !
    ! delta=s:sigma
    !
    ! dKappa=s+Beta*dKappa

    if (CI) then
      call ipin(ipS2)
      call ipin(ipST)
      deltaC = ddot_(nConf1,W(ipST)%A,1,W(ipS2)%A,1)
      call ipout(ipST)
    else
      deltaC = Zero
    end if

    deltaK = ddot_(nDensC,Sigma,1,Sc4,1)
    if (.not. CI) then
      rBeta = deltaK/delta
      delta = deltaK
      Temp2(:) = rBeta*Temp2(:)+Sc4(:)
    else
      rbeta = (deltac+deltaK)/delta
      delta = deltac+deltaK
      call ipin(ipCID)
      W(ipCID)%A(1:nConf1) = rBeta*W(ipCID)%A(1:nConf1)
      Temp2(:) = rBeta*Temp2(:)
      call ipin(ipS2)
      W(ipCID)%A(1:nConf1) = W(ipCID)%A(1:nConf1)+W(ipS2)%A(1:nConf1)
      Temp2(:) = Temp2(:)+Sc4(:)
      call opout(ipS2)
      call ipout(ipCID)
    end if
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'rBeta, DeltaK, DeltaC=',rBeta,DeltaK,DeltaC
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !    ######  #    #  #####        #####    ####    ####
    !    #       ##   #  #    #       #    #  #    #  #    #
    !    #####   # #  #  #    #       #    #  #       #
    !    #       #  # #  #    #       #####   #       #  ###
    !    #       #   ##  #    #       #       #    #  #    #
    !    ######  #    #  #####        #        ####    ####
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call UnCompress(Temp2,dKappa,isym)

    res = Zero ! dummy initialize
    res_tmp = -One
    if (iBreak == 1) then
      if (abs(delta) < abs(Eps**2*delta0)) exit
    else if (ibreak == 2) then
      res = sqrt(resk**2+resci**2)
      if (doDMRG) then ! yma
        !write(u6,*) 'resk**2, resci**2 ',resk**2,resci**2
        res = sqrt(resk**2+resci**2)
        ! And a bit loose in DMRG case
        if (res < abs(Eps)) exit
        if (sqrt(resk**2) < abs(Eps)) then
          if (abs(res_tmp-sqrt(resci**2)) < 1.0e-6_wp) exit
        end if
        res_tmp = sqrt(resci**2)
      else
        if (res < abs(Eps)) exit
      end if
    else
      if ((abs(delta) < abs(Eps**2*delta0)) .and. (res < abs(Eps))) exit
    end if
    if (iter >= niter) then
      cnvrgd = .false.
      exit
    end if
    if (lprint) &
      write(LuWr,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)') '     ',iter,'       ',delta/delta0,resk,resci,deltac,deltak

    iter = iter+1

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (.not. cnvrgd) then
    write(LuWr,Fmt2//'A,I4,A)') 'No convergence for perturbation no: ',idisp,'. Increase Iter.'
    converged(isym) = .false.
    fail = .true.
#   ifdef _MOLCAS_MPP_
    ! Set the global flag to signal a process failed
    if (Is_Real_Par()) then
      dfail = One
      call GA_Acc(iglfail,1,1,1,1,dfail,1,One)
    end if
#   endif

    exit
  else
    write(LuWr,Fmt2//'A,I4,A,I4,A)') 'Perturbation no: ',idisp,' converged in ',iter-1,' steps.'
    call ipnout(-1)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !         Write response to disk                                     *
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  !write(LuWr,Fmt2//'A)') 'Writing response to one-file.'
  write(LuWr,*)
  iLen = nDensC
  iKapDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,Kappa,iLen,iDis)
  iSigDisp(iDisp) = iDis
  call dDaFile(LuTemp,1,Sigma,iLen,iDis)
  if (CI) then
    ilen = nconf1
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
  call mma_deallocate(rmoaa,safe='*')
  if (CI) then
    call mma_deallocate(Pens)
    call mma_deallocate(Dens)
  end if
end do
if (nProcs >= 2) then
  write(LuWr,*)
  write(LuWr,*) ' Perturbations were printed only by master node'
  write(LuWr,*)
end if

call Free_Tsk(id)
call mma_deallocate(List)
!                                                                      *
!***********************************************************************
!                                                                      *
! Free all memory and remove from disk all data
! related to the last symmetry

call mma_deallocate(iPre)
if (CI) then
  call ipclose(ipdia)
else
  call ipclose(ipPre2)
end if

call Exp_Close()
!                                                                      *
!***********************************************************************
!                                                                      *
! Flush the output from the other nodes.

do iRank=1,nProcs-1
  call GASync()
  if (iRank == MyRank) then
    rewind(LuWr)
    do
      read(LuWr,'(A)',iostat=istatus) Line
      if (istatus < 0) exit
      write(LuWr_Save,*) Line
    end do
  end if
  call GASync()
end do
if (MyRank /= 0) then
  close(LuWr)
  LuWr = LuWr_save
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(LuWr,*) '********************************************************************************'
write(LuWr,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Final synchronization of the fail flag
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GAdGOp_Scal(dfail,'max')
  Fail = Fail .or. (dfail > Zero)
end if
#endif
if (Fail) call Quit_OnConvError()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine WfCtl_Hess
