************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Anders Bernhardsson                                    *
*               2002, Roland Lindh                                     *
************************************************************************
      SubRoutine WfCtl_Hess(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                      iRHSDisp,iRHSCIDISP,converged)
************************************************************************
*                                                                      *
*                                                                      *
*     called from: MCLR                                                *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*     Parallelization of perturbations, RL 2002.                       *
*                                                                      *
************************************************************************
      use Exp, only: Exp_Close
      use Arrays, only: CMO, Int2, FIMO
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)
      External Rsv_Tsk
*
#include "standard_iounits.fh"
#include "real.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "machine.fh"
#include "stdalloc.fh"
#include "dmrginfo_mclr.fh"
*
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
#endif
      Logical Orb,CI,Response
      Parameter (iTimeCC = 1 )
      Parameter (iTimeKK = 2 )
      Parameter (iTimeKC = 3 )
      Parameter (iTimeCK = 4 )
#include "crun_mclr.fh"
      Character*8   Fmt2
      Character*132 Line
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp),iRHSCIDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      Integer pstate_sym,opout
      Logical lPrint,converged(8), Rsv_Tsk
      Real*8 Clock(4)
      Character*72 SLine
      Real*8 res_tmp
      Real*8 rdum(1)
      Real*8, Allocatable:: Kappa(:), dKappa(:), Sigma(:),
     &                      Temp1(:), Temp2(:), Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Sc3(:),
     &                      Dens(:), Pens(:), rmoaa(:)
      Integer, Allocatable:: List(:,:)

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*define _DEBUGPRINT_
*----------------------------------------------------------------------*
      SLine=' Solving CP(CAS)HF equations'
      Call StatusLine(' MCLR:',SLine)
*
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      TIM2=Zero
      TIM3=Zero
      TIM4=Zero
*
      CLOCK(:)=Zero
      lPaper=132
      lLine =120
      left=(lPaper-lLine)/2
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*

      iDis=0
*
      fail=.false.
      Converged(:)=.true.
      lprint=.false.
      idasave=0
*     If (SAVE) CALL DANAME(50,'RESIDUALS')
      If (SAVE) Then
         Write (LuWr,*) 'WfCtl: SAVE option not implemented'
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Start loop over the symmetry of the perturbations
*
      If (iAnd(kprint,2).eq.2) lprint=.true.
#ifdef _DEBUGPRINT_
      lprint=.true.
#endif
      kksym=1
      kkksym=nsym
      If (PT2) kkkSym=1
*
*     Set up parallelization over the loop over perturbations
*
      Call mma_allocate(List,2,nDisp,Label='List')
      iDisp=0
      Do iSym=kksym,kkksym
        iDEnd=lDisp(iSym)
        Do jDisp=1,iDEnd
           iDisp=iDisp+1
           List(1,iDisp)=iSym
           List(2,iDisp)=jDisp
        End Do
      End Do
*
*     Change output unit
*
      Call Get_MyRank(MyRank)
      LuWr_save=LuWr
      If (MyRank.ne.0) Then
         LuWr=55
         LuWr=isFreeUnit(LuWr)
         call molcas_open(luwr,'Temp_OutPut')
      End If
*
      Call Init_Tsk(id,nDisp)
#ifdef _MOLCAS_MPP_
*     iglfail is global "array", a flag used to communicate to all processes
*     if any of them failed, the communication need not be synchronous,
      If (Is_Real_Par()) Then
         If (.Not.GA_Create(MT_DBL,1,1,'GlFail',0,0,iglfail)) Then
           Call SysAbendMsg ('wfctl_hess',
     &                       'failed to create global failed flag',' ')

         End If
         Call GA_Zero(iglfail)
         dfail=Zero
      End If
#endif
*
      ipdia=0
      ipPre2=0
      iSym_Old=0
 888  If (.Not.Rsv_Tsk(id,iDisp)) Go To 999
#ifdef _MOLCAS_MPP_
*       Check if some other process failed
        If (Is_Real_Par()) Then
           Call GA_Get(iglfail,1,1,1,1,dfail,1)
           If (dfail.gt.Zero) Then
              fail=.true.
              Go To 999
           End If
        End If
#endif
        iSym =List(1,iDisp)
        jDisp=List(2,iDisp)
*
        Write (SLine,'(A,I3,A)')
     &        ' Solving CP(CAS)HF equations for perturbation ',
     &        iDisp,'.'
        Call StatusLine(' MCLR:',SLine)
*
C     Do iSym=kksym,kkksym
*
*       Execute setup for symmetry block if a new one!
*
        If (iSym.ne.iSym_Old) Then
*
           If (iSym_Old.ne.0) Then
*
*             If a previous symmetry block as processed
*             free all memory and remove from disk all data
*             related to this symmetry
*
              If (CI) Then
                 irc=ipclose(ipdia)
              Else
                 irc=ipclose(ipPre2)
              End If
*
              Call Exp_Close()
*
           End If
           iSym_Old=iSym
*
        PState_SYM=iEor(State_Sym-1,iSym-1)+1
        nconf1=ncsf(PState_Sym)
        CI=.false.
        If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.

        If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.
*          Initiate CSF <-> SD
        If (CI) Call InCSFSD(iEor(iSym-1,State_Sym-1)+1,
     &                       State_sym,.false.)
*
*
*          Calculate length of the density, Fock and Kappa matrix etc
*          notice that this matrixes not necessary are symmetric.
*          Store pointers.
*
*          Input:
*                 iSym : Symmetry of perturbation
*
*          Output: Commonblocks (Pointers.fh)
*
        PState_SYM=iEor(State_Sym-1,iSym-1)+1
        nConf2=nint(xispsm(PState_SYM,1))
*       nConf2=ndtasm(PState_SYM)
        nconf3=nint(Max(xispsm(PState_SYM,1),xispsm(State_SYM,1)))
*       nconf3=Max(ndtasm(PState_SYM),ndtasm(State_SYM))

        if(doDMRG)then  ! yma
          call dmrg_spc_change_mclr(RGras2(1:8),nash)
          call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
        end if

        Call Setup_MCLR(iSym)
*                                                                      *
************************************************************************
*                                                                      *
*       Determine if we should page CI vectors
*
*                                    [2]
*         Calculate the diagonal of E    and store in core/disc
*
        If (CI) Then
            Call CIDia(PState_Sym,rCHC)
            irc=ipout(ipdia)
*
*       Allocate disk/memory space
*
*
*       This areas should be addressed through ipin
*       ipout will page them out to disk and free the memory area
*       if page mode is used
*
*       opout will release the memory area without update the disk
*
           ips1 =ipget(nconf3)
           ips2 =ipget(nconf3)
           ipst =ipget(nconf3)
           ipcit=ipget(nconf1)
           ipcid=ipget(nconf1)
*
        Else
           ips1=0
           ips2=0
           ipst=0
           ipcit=0
           ipcid=0
        End If
*
        npre2=Max(npre(isym),1) ! "Max" just there to make sure that
!                                 W(ipPre2) is allocated even if
!                                 npre2(isym) is zero.
        ipPre2=ipget(npre2)
        irc=ipin(ipPre2)
        call Prec(W(ipPre2)%Vec,isym)
        irc=ipout(ippre2)
*                                                                      *
************************************************************************
*                                                                      *
*      OK START WORKING, looping over the perturbations of the
*      current symmetry
*
        iDEND=lDisp(iSym)
*       If (SewLab.eq.'NONE    '.and.(.not.mckinley)) iDEND=1
*
        End If
*
C       Do jDisp=1,iDEnd
C         iDisp=iDisp+1
          jspin=0
          If (iAnd(nTPert(idisp),1).eq.1) jSpin=1
          If (jspin.eq.0) Then
              nconf1=ncsf(PState_Sym)
          Else
              nConf1=nint(xispsm(Pstate_Sym,1))
          End If
*
*    Allocate areas for scratch and state variables
*
          Call mma_allocate(Kappa,nDens2+6,Label='Kappa')
          Call mma_allocate(dKappa,nDens2+6,Label='dKappa')
          Call mma_allocate(Sigma,nDens2+6,Label='Sigma')
          Call mma_allocate(Temp1,nDens2+6,Label='Temp1')
          Call mma_allocate(Temp2,nDens2+6,Label='Temp2')
          Call mma_allocate(Temp3,nDens2+6,Label='Temp3')
          Call mma_allocate(Temp4,nDens2+6,Label='Temp4')
          Call mma_allocate(Sc1,nDens2+6,Label='Sc1')
          Call mma_allocate(Sc2,nDens2+6,Label='Sc2')
          Call mma_allocate(Sc3,nDens2+6,Label='Sc3')
          Temp1(1:nDens2)=Zero
          Kappa(1:nDens2)=Zero
          dKappa(1:nDens2)=Zero
          Sigma(1:nDens2)=Zero
          If (CI) Then
             Call mma_allocate(Dens,n1dens,Label='Dens')
             Call mma_allocate(Pens,n2dens,Label='Pens')
*
             Dens(:)=Zero
             Pens(:)=Zero
          End If
          If (iMethod.eq.2) Then
             Call mma_allocate(rmoaa,n2dens,Label='rmoaa')
          Else
             Call mma_allocate(rmoaa,1,Label='rmoaa')
          End If
          rmoaa(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*         Calculate RHS for the perturbation
*                                                                      *
************************************************************************
*                                                                      *
*         (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)
*
          If (PT2) then
             Call RHS_PT2(Temp4,ipST)
          Else
             Call RHS(Sigma,Kappa,Temp1,
     &                Temp3,Sc2,dKappa,
     &                Sc3,
     &                Temp4,ipST,
     &                iDisp,iSym-1,CMO,jdisp,jspin,CI)
#ifdef _DEBUGPRINT_
             Write (LuWr,*) 'After RHS'
             Write (LuWr,*) 'Sigma=',DDot_(nDens,Sigma,1,Sigma,1)
             Write (LuWr,*) 'Kappa=',DDot_(nDens,Kappa,1,Kappa,1)
             Write (LuWr,*) 'Sc2=',DDot_(nDens,Sc2,1,Sc2,1)
             Write (LuWr,*) 'dKap=',DDot_(nDens,dKappa,1,dKappa,1)
             Write (LuWr,*) 'CMO=',DDot_(nCMO,CMO,1,CMO,1)
#endif
          End If
          irc=opout(ipci)
*
          Write (LuWr,*) 'Process perturbation number ',iDisp
          If (lprint) Write(LuWr,*)
     &      '       Iteration         Delta     Res(kappa) Res(CI)'
*                                                                      *
************************************************************************
*                                                                      *
*         Write RHS to disk                                            *
*                                                                      *
************************************************************************
*                                                                      *
          iLen=nDensC
          iRHSDisp(iDisp)=iDis
          Call Compress(Temp4,Sigma,iSym)
          r1=ddot_(ndensc,sigma,1,sigma,1)
          Call UnCompress(Sigma,Temp4,iSym)
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (CI) Then
             irc=ipin(ipCIT)
             W(ipCIT)%Vec(1:nConf1)=Zero
          End If
          irc=ipout(ipCIT)
          If (CI) Then
             ilen=nconf1
             iRHSCIDisp(iDisp)=iDis
             irc=ipin(ipST)
             Call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
             Call DSCAL_(nConf1,-One,W(ipST)%Vec,1)
          End If
*
          Call DSCAL_(nDensC,-One,Sigma,1)
*
#ifdef _DEBUGPRINT_
          irc=ipin(ipST)
          Write (LuWr,*) 'ST=',DDot_(nConf1,W(ipST)%Vec,1,W(ipST)%Vec,1)
          Write (LuWr,*) 'Sigma=',DDot_(nDensC,Sigma,1,Sigma,1)
          Write (LuWr,*) 'Kappa=',DDot_(nDens2,Kappa,1,Kappa,1)
          Write (LuWr,*) 'End RHS!'
#endif
*
          iter=1
          irc=ipin(ipPre2)
          Call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,
     &                  Kappa,nDens2+6,Temp3,nDens2+6,
     &                  isym,iter)
          irc=opout(ippre2)
          r2=ddot_(ndensc,Kappa,1,Kappa,1)
#ifdef _DEBUGPRINT_
          Write (LuWr,*) 'DMinvKap'
          Write (LuWr,*) 'Kap=',DDot_(nDens2,Kappa,1,Kappa,1)
          Write (LuWr,*) 'End DMinvKap'
#endif
          If (r2.gt.r1) Write(LuWr,Fmt2//'A,I3,A)')
     &         'Warning  perturbation number ',idisp,' might diverge!'
          Call UnCompress(Kappa,dKappa,iSym)
*
          If (CI) Then
             irc=ipin(ipCId)
             Call DMinvCI(ipST,W(ipCId)%Vec, rCHC,isym)
             irc=ipin(ipST)
             deltaC=ddot_(nConf1,W(ipST)%Vec,1,W(ipCId)%Vec,1)
             irc=ipout(ipcid)
          Else
             deltaC=Zero
          End If
          deltaK=ddot_(nDensC,Kappa,1,Sigma,1)
          Kappa(1:nDens)=Zero
          delta=deltac+deltaK
#ifdef _DEBUGPRINT_
          If (Abs(DeltaC).lt.1.0D-12) DeltaC=Zero
          Write (LuWr,*)'DeltaK, DeltaC, Delta=',DeltaK, DeltaC, Delta
#endif

          If (delta.eq.Zero) Goto 300
          delta0=delta
#ifdef _DEBUGPRINT_
          Write (LuWr,*) 'Delta0=',Delta0
          Write (LuWr,*) 'Start ITERATIONS'
          Write (LuWr,*) 'Orb,CI=',ORB,CI
#endif
          Orb=.true.
          TimeDep=.false.
          ReCo=-One
*                                                                      *
************************************************************************
*          I   T   E   R   A   T   I   O   N   S                       *
************************************************************************
*                                                                      *
200       Continue

          if(doDMRG)then  ! It can be deleted
            call dmrg_spc_change_mclr(RGras2(1:8),nash)
            call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
          end if

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*       O R B I T A L    P A R T
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
           If (orb) Then
*                                                                      *
************************************************************************
*                                                                      *
*                      ~    ~
*            Construct F,(ij|kl)
*                                                                      *
************************************************************************
*                                                                      *

             irc=ipnout(-1)
             Call RInt_generic(dKappa,rmoaa,rdum,
     &                         Sc2,Temp3,Temp4,
     &                         Sc3,isym,reco,jspin)
             Clock(iTimeKK)=Clock(iTimeKK)+Tim2

*                                                                      *
************************************************************************
*                                                                      *
*            kappa->CI
*
*            H(kappa)|0>
*
*                [2]
*            S1=E   k ( kappa TO CI )
*                                                                      *
************************************************************************
*                                                                      *
             If (CI) Then
                Call CISigma(jspin,State_Sym,pstate_sym,
     &                       Temp4,nDens2,rmoaa,SIZE(rmoaa),rdum,1,
     &                       ipCI,ipS1,.True.)
                Clock(iTimeKC)=Clock(iTimeKC)+Tim3
#ifdef _DEBUGPRINT_
                irc=ipin(ipCI )
                irc=ipin(ipS1 )
                Write (LuWr,*) 'CISigma'
                Call RecPrt('CI ','(3F10.4)',W(ipCI )%Vec,1,nConf1)
                Call RecPrt('S1 ','(3F10.4)',W(ipS1 )%Vec,1,nConf1)
                Write (LuWr,*) 'CI=',DDot_(nConf1,W(ipCI)%Vec,1,
     &                                            W(ipCI)%Vec,1)
                Write (LuWr,*) 'S1=',DDot_(nConf1,W(ipS1)%Vec,1,
     &                                            W(ipS1)%Vec,1)
                Write (LuWr,*) 'End CISigma'
#endif
*
*               This will give us a better
*               convergence in the PCG. Notice that
*
*               ~Inactive     ~
*               E        + <0|H|0> = 0
*
*               when the wavefunction is converged.
*
                irc=ipin(ipS1)
                If (isym.eq.1) Then
                   irc=ipin(ipCI)
                   rGrad=DDot_(nconf1,W(ipCI)%Vec,1,W(ipS1)%Vec,1)
                   call daxpy_(nConf1,-rgrad,W(ipCI)%Vec,1,
     &                                       W(ipS1)%Vec,1)
                End IF
                call dscal_(nconf1,Two,W(ipS1)%Vec,1)
*
                irc=opout(ipCI)
             End If
*                                                                      *
************************************************************************
*                                                                      *
           End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*       C I    P A R T
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
           If (CI) Then
*                                                                      *
************************************************************************
*                                                                      *
*                  [2]
*            S2 = E   CID  ( CI TO CI ) = <i|H|d> -E<i|d>
*                                                                      *
************************************************************************
*                                                                      *
              irc=ipnout(-1)
              Call CISigma(0,PState_Sym,Pstate_sym,
     &                     FIMO,SIZE(FIMO),Int2,SIZE(Int2),rdum,1,ipCId,
     &                     ipS2,.True.)
              EC=rin_ene+potnuc-ERASSCF(1)

              irc=ipin(ipCId)
              irc=ipin(ipS2)
              Call DaXpY_(nConf1,EC,W(ipCId)%Vec,1,W(ipS2)%Vec,1)
              Call DSCAL_(nConf1,Two,W(ipS2)%Vec,1)
              Clock(iTimeCC)=Clock(iTimeCC)+Tim4
              irc=ipout(ipS2)
              irc=opout(ipCId)
              irc=opout(ipCI)
*                                                                      *
************************************************************************
*                                                                      *
*             CI -> Kappa
*
*                  [2]
*             SC3=E   CID   (SC3=F(<d|E|0>+<0|E|d>)
*                                                                      *
************************************************************************
*                                                                      *
              Response=.true.
              irc=ipnout(-1)
#ifdef _DEBUGPRINT_
              irc=ipin(ipCI )
              irc=ipin(ipCId)
              Call RecPrt('CI ','(3F10.4)',W(ipCI )%Vec,1,nConf1)
              Call RecPrt('CId','(3F10.4)',W(ipCId)%Vec,1,nConf1)
#endif
              Call CIDens(Response,ipCI,ipCId,
     &                    State_sym,
     &                    PState_Sym,jspin,
     &                    Pens,Dens)     ! Jeppes

#ifdef _DEBUGPRINT_
              Write (LuWr,*) 'After CIDens'
              Write (LuWr,*) 'P=',DDot_(n2dens,Pens,1,Pens,1)
              Write (LuWr,*) 'De=',DDot_(n1dens,Dens,1,Dens,1)
#endif
*
*             density for inactive= 2(<d|0>+<0|d>)
*
              d_0=Zero
*
*             This is just for debugging purpose.
*             When we use it for actual calculations d_0 == 0
*
              If (isym.eq.1) Then
                 irc=ipin(ipCI)
                 irc=ipin(ipCid)
                 d_0=ddot_(nconf1,W(ipCid)%Vec,1,W(ipCI)%Vec,1)
              End If
              If (Response) d_0=d_0*Two
*

              Call FockGen(d_0,Dens,Pens,Sc1,Sc3,isym)    ! Made
*                                                                      *
************************************************************************
*                                                                      *
           End If
*                                                                      *
************************************************************************
*                                                                      *
*        Sc1  kappa-> kappa
*        Sc3  CI -> kappa
*        S1   kappa -> CI
*        S2   CI -> CI
*        dKap present step
*        Kap  kappaX
*        CIT  CIX
*        CId  present step
*                                                                      *
************************************************************************
*                                                                      *
*        Add together
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'Add together!'
           Write (LuWr,*) 'dKap=',DDot_(nDens,dKappa,1,dKappa,1)
           Write (LuWr,*) 'Sc2=',DDot_(nDens,Sc2,1,Sc2,1)
           If (CI) Then
              Write (LuWr,*) 'Sc3=',DDot_(nDens,Sc3,1,Sc3,1)
              irc=pin1(ipS2,nconf1)
              Write(LuWr,*)'S2=',DDot_(nConf1,W(ipS2)%Vec,1,
     &                                        W(ipS2)%Vec,1)
              irc=pin1(ipS1,nconf1)
              Write(LuWr,*)'S1=',DDot_(nConf1,W(ipS1)%Vec,1,
     &                                        W(ipS1)%Vec,1)
           End If
           Write (LuWr,*)
#endif
           irc=ipnout(-1)
           if (CI) then
              Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)
           Else
              call dcopy_(nDens,Sc2,1,Sc1,1)
           End If
           Call Compress(Sc1,Temp4,isym)   ! ds
           Call Compress(dKappa,Temp2,isym)  ! DX
           If (CI) Then
              irc=ipin1(ipS1,nconf1)
              irc=ipin1(ipS2,nconf1)
              Call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipS1)%Vec,1)
              irc=opout(ipS2)
           End If
*                                                                      *
************************************************************************
*                                                                      *
*                ######   #####   #####
*                #     # #     # #     #
*                #     # #       #
*                ######  #       #  ####
*                #       #       #     #
*                #       #     # #     #
*                #        #####   #####
*                                                                      *
************************************************************************
*                                                                      *
*                     delta
*          rAlpha=------------
*                 dKappa:dSigma
*
           rAlphaC=Zero
           rAlphaK=Zero
           irc=ipin(ipS1)
           irc=ipin(ipCId)
           If (orb) rAlphaK=DDot_(nDensC,Temp4,1,Temp2,1)
           If (CI)  rAlphaC=DDot_(nConf1,W(ipS1)%Vec,1,
     &                                   W(ipCId)%Vec,1)
           rAlpha=delta/(rAlphaK+rAlphaC)
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'At PCG'
           Write (LuWr,*) 'S1=',
     &                      DDot_(nConf1,W(ipS1 )%Vec,1,
     &                                   W(ipS1 )%Vec,1)
           Write (LuWr,*) 'CId=',
     &                      DDot_(nConf1,W(ipCId)%Vec,1,
     &                                   W(ipCId)%Vec,1)
           Write (LuWr,*) 'rAlphaK, rAlphaC, rAlpha=',
     &                     rAlphaK, rAlphaC, rAlpha
#endif
*                                                                      *
************************************************************************
*                                                                      *
*          Kappa=Kappa+rAlpha*dKappa
*          Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
*
           If (orb) Then
             Call DaxPy_(nDensC,ralpha,Temp2,1,Kappa,1)
             Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
              resk=sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
           End If
           resci=Zero
           If (CI) Then
              irc=ipin(ipCIT)
              Call DaXpY_(nConf1,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
              irc=ipout(ipCIT)
              irc=ipin1(ipST,nconf1)
              irc=ipin(ipS1)
              Call DaXpY_(nConf1,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
              irc=opout(ipS1)
              irc=ipin(ipST)
              resci=sqrt(ddot_(nconf1,W(ipST)%Vec,1,W(ipST)%Vec,1))
           End If
*                                                                      *
************************************************************************
*                                                                      *
*          Precondition......
*             -1
*          S=M  Sigma
*
          irc=opout(ipCID)
          irc=ipin(ipS2)
          If (CI) Call DMinvCI(ipST,W(ipS2)%Vec,rCHC,isym)
          irc=opout(ipCI)
          irc=opout(ipdia)
*
          irc=ipin(ipPre2)
          Call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,
     &                  Sc2,nDens2+6,Sc1,nDens2+6,iSym,iter)
          irc=opout(ippre2)
*                                                                      *
************************************************************************
*                                                                      *
*               s:Sigma
*          Beta=-------
*                delta
*
*          delta=s:sigma
*
*          dKappa=s+Beta*dKappa
*
           If (CI) Then
              irc=ipin(ipS2)
              irc=ipin(ipST)
              deltaC=ddot_(nConf1,W(ipST)%Vec,1,W(ipS2)%Vec,1)
              irc=ipout(ipST)
           Else
              deltaC=Zero
           End If
*
           deltaK=ddot_(nDensC,Sigma,1,Sc2,1)
           If (.not.CI) Then
              rBeta=deltaK/delta
              delta=deltaK
              Call DScal_(nDensC,rBeta,Temp2,1)
              Call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
           Else
              rbeta=(deltac+deltaK)/delta
              delta=deltac+deltaK
              irc=ipin(ipCID)
              Call DScal_(nConf1,rBeta,W(ipCID)%Vec,1)
              Call DScal_(nDensC,rBeta,Temp2,1)
              irc=ipin(ipS2)
              Call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
              Call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
              irc=opout(ipS2)
              irc=ipout(ipCID)
           End If
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'rBeta, DeltaK, DeltaC=',
     &                     rBeta, DeltaK, DeltaC
#endif
*                                                                      *
************************************************************************
*                                                                      *
*    ######  #    #  #####        #####    ####    ####
*    #       ##   #  #    #       #    #  #    #  #    #
*    #####   # #  #  #    #       #    #  #       #
*    #       #  # #  #    #       #####   #       #  ###
*    #       #   ##  #    #       #       #    #  #    #
*    ######  #    #  #####        #        ####    ####
*                                                                      *
************************************************************************
*                                                                      *
           Call UnCompress(Temp2,dKappa,isym)
*
           res=Zero ! dummy initialize
           res_tmp=-One
           If (iBreak.eq.1) Then
              If (abs(delta).lt.abs(Epsilon**2*delta0)) Goto 300
           Else If (ibreak.eq.2) Then
              res=sqrt(resk**2+resci**2)
              if (doDMRG) then ! yma
!                write(*,*)"resk**2, resci**2",resk**2,resci**2
                res=sqrt(resk**2+resci**2)
                ! And a bit loose in DMRG case
                If (res.lt.abs(epsilon)) Goto 300
                if (sqrt(resk**2).lt.abs(epsilon))then
                  if (abs(res_tmp-sqrt(resci**2)).lt.1.0e-06)then
                    goto 300
                  end if
                end if
                res_tmp=sqrt(resci**2)
              else
                If (res.lt.abs(epsilon)) Goto 300
              end if
           Else
              If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &            res.lt.abs(epsilon))  Goto 300
           End If
           If (iter.ge.niter) goto 210
           If (lprint)
     &     Write(LuWr,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)')
     &            '     ',
     &            iter,'       ',delta/delta0,resk,resci,deltac,deltak

           iter=iter+1

          Goto 200
*                                                                      *
************************************************************************
*                                                                      *
 210      Continue
          Write(LuWr,Fmt2//'A,I4,A)')
     &    'No convergence for perturbation no: ',
     &                      idisp,'. Increase Iter.'
          converged(isym)=.false.
          fail=.true.
#ifdef _MOLCAS_MPP_
*         Set the global flag to signal a process failed
          If (Is_Real_Par()) Then
             dfail=One
             Call GA_Acc(iglfail,1,1,1,1,dfail,1,One)
          End If
#endif
*
C         Goto 310
          Go To 999
*
 300      Write(LuWr,Fmt2//'A,I4,A,I4,A)')
     &           'Perturbation no: ',idisp,' converged in ',
     &             iter-1,' steps.'
          irc=ipnout(-1)
*310      Continue
*                                                                      *
************************************************************************
*                                                                      *
*         Write response to disk                                       *
*                                                                      *
************************************************************************
*                                                                      *
*
C         Write(LuWr,Fmt2//'A)')'Writing response to one-file.'
          Write(LuWr,*)
          iLen=ndensC
          iKapDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Kappa,iLen,iDis)
          iSigDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (CI) Then
             ilen=nconf1
             iCIDisp(iDisp)=iDis
             irc=ipin(ipCIT)
             Call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)
             iCISigDisp(iDisp)=iDis
             irc=ipin(ipST)
             Call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
          End If
*
          If (CI) Then
          End If

          Call mma_deallocate(Temp4)
          Call mma_deallocate(Temp3)
          Call mma_deallocate(Temp2)
          Call mma_deallocate(Temp1)
          Call mma_deallocate(Sigma)
          Call mma_deallocate(dKappa)
          Call mma_deallocate(Kappa)
          Call mma_deallocate(Sc3)
          Call mma_deallocate(Sc2)
          Call mma_deallocate(Sc1)
          If (allocated(rmoaa)) Call mma_deallocate(rmoaa)
          If (CI) Then
             Call mma_deallocate(Pens)
             Call mma_deallocate(Dens)
          End If
      Go To 888
 999  Continue
      Call Get_nProcs(nProcs)
      If(nProcs.ge.2) Then
       Write (LuWr,*)
       Write (LuWr,*) ' Perturbations were printed only by master node'
       Write (LuWr,*)
      End If

      Call Free_Tsk(id)
      Call mma_deallocate(List)
*                                                                      *
************************************************************************
*                                                                      *
*      Free all memory and remove from disk all data
*      related to the last symmetry
*
      If (CI) Then
         irc=ipclose(ipdia)
      Else
         irc=ipclose(ipPre2)
      End If
*
      Call Exp_Close()
*                                                                      *
************************************************************************
*                                                                      *
*     Flush the output from the other nodes.
*
      Do iRank = 1, nProcs-1
         Call GASync()
         If (iRank.eq.MyRank) Then
            ReWind(LuWr)
 777        Read(LuWr,'(A)',END=778) Line
            Write (LuWr_Save,*) Line
            Go To 777
 778        Continue
         End If
         Call GASync()
      End Do
      If (MyRank.ne.0) Then
         Close(LuWr)
         LuWr=LuWr_save
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write(LuWr,*)  '****************************************',
     &               '****************************************'
      Write(LuWr,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Final synchronization of the fail flag
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
         Call GAdGOp_Scal(dfail,'max')
         Fail=Fail.or.(dfail.gt.0.0d0)
      End If
#endif
      If (Fail) Call Quit_OnConvError()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
