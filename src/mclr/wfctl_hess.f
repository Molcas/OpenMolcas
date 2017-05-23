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
      Implicit Real*8 (a-h,o-z)
      External Rsv_Tsk
*
#include "standard_iounits.fh"
#include "WrkSpc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "csfbas_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "machine.fh"
#include "stdalloc.fh"
#include "dmrginfo_mclr.fh"
*
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
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Call QEnter('WfCtl')
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*define _DEBUG_
*----------------------------------------------------------------------*
      one=1.0d0
      SLine=' Solving CP(CAS)HF equations'
      Call StatusLine(' MCLR:',SLine)
*
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      TIM2=0.0D0
      TIM3=0.0D0
      TIM4=0.0D0
*
      call dcopy_(4,0.0d0,0,CLOCK,1)
      lPaper=132
      lLine =120
      left=(lPaper-lLine)/2
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*

      iDis=0
*
      fail=.false.
      Do i=1,8
       Converged(i)=.true.
      end do
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
#ifdef _DEBUG_
      lprint=.true.
#endif
      kksym=1
      kkksym=nsym
      If (PT2) kkkSym=1
*
*     Set up parallelization over the loop over perturbations
*
      Call Allocate_iWork(ipList,2*nDisp)
      iDisp=0
      Do iSym=kksym,kkksym
        iDEnd=lDisp(iSym)
        Do jDisp=1,iDEnd
           iDisp=iDisp+1
           iWork(ipList+(iDisp-1)*2  )=iSym
           iWork(ipList+(iDisp-1)*2+1)=jDisp
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
*
      ipdia=0
      ipPre2=0
      iSym_Old=0
 888  If (.Not.Rsv_Tsk(id,iDisp)) Go To 999
        iSym =iWork(ipList+(iDisp-1)*2  )
        jDisp=iWork(ipList+(iDisp-1)*2+1)
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
              If (iphx.ne.0) Then
                 Call Getmem('EXPHS','FREE','REAL',iphx,idum)
                 Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
                 Call Getmem('EXPLS','FREE','INTE',iplst,idum)
              End If
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
        iphx=0
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
        npre2=npre(isym)
        ipPre2=ipget(npre2)
        call Prec(Work(ipin(ipPre2)),isym)
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
          Call GetMem('kappa ','Allo','Real',ipKap  ,nDens2+6)
          Call GetMem('dkappa','Allo','Real',ipdKap ,nDens2+6)
          Call GetMem('sigma ','Allo','Real',ipSigma,nDens2+6)
          Call GetMem('Temp1 ','Allo','Real',ipTemp1,nDens2+6)
          Call GetMem('Temp2 ','Allo','Real',ipTemp2,nDens2+6)
          Call GetMem('Temp3 ','ALLO','Real',ipTemp3,nDens2+6)
          Call GetMem('Temp4 ','Allo','Real',ipTemp4,nDens2+6)
          Call Getmem('Scr1  ','ALLO','Real',ipSc1  ,nDens2+6)
          Call Getmem('Scr3  ','ALLO','Real',ipSc3  ,nDens2+6)
          Call Getmem('Scr2  ','ALLO','Real',ipSc2  ,nDens2+6)
          call dcopy_(nDens2,0.0d0,0,Work(ipTemp1),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipKap),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipsigma),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipdKap),1)
          If (CI) Then
             Call GetMem('1Dens','ALLO','Real',ipDe,n1dens)
             Call GetMem('2Dens','ALLO','Real',ipP,n2dens)
*
             call dcopy_(n1dens,0.0d0,0,Work(ipDe),1)
             call dcopy_(n2dens,0.0d0,0,Work(ipP),1)
          End If
          If (iMethod.eq.2) Then
             Call GetMem('RMOAA','ALLO','Real',iprmoaa,n2dens)
             Call FZero(Work(iprmoaa),n2dens)
          Else
             iprmoaa = ip_iDummy
          End If
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
             Call RHS_PT2(Work(ipTemp4),ipST)
          Else
             kain=ipkap

             Call RHS(Work(ipSigma),Work(ipKap),Work(ipTemp1),
     &                Work(ipTemp3),Work(ipSc2),Work(ipdKap),
     &                Work(ipSc3),
     &                Work(ipTemp4),ipST,
     &                iDisp,iSym-1,Work(ipCMO),jdisp,jspin,CI)
#ifdef _DEBUG_
             Write (LuWr,*) 'After RHS'
             Write (LuWr,*) 'Sigma=',DDot_(nDens,Work(ipSigma),1,
     &                                          Work(ipSigma),1)
             Write (LuWr,*) 'Kap=',DDot_(nDens,Work(ipKap),1,
     &                                        Work(ipKap),1)
             Write (LuWr,*) 'Sc2=',DDot_(nDens,Work(ipSc2),1,
     &                                        Work(ipSc2),1)
             Write (LuWr,*) 'dKap=',DDot_(nDens,Work(ipdKap),1,
     &                                         Work(ipdKap),1)
             Write (LuWr,*) 'CMO=',DDot_(nCMO,Work(ipCMO),1,
     &                                       Work(ipCMO),1)
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
          Call Compress(Work(ipTemp4),Work(ipSigma),iSym)
          r1=ddot_(ndensc,Work(ipsigma),1,Work(ipsigma),1)
          Call UnCompress(Work(ipSigma),Work(ipTemp4),iSym)
          Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
          If (CI) call dcopy_(nConf1,0.0d0,0,Work(ipin(ipCIT)),1)
          irc=ipout(ipCIT)
          If (CI) Then
             ilen=nconf1
             iRHSCIDisp(iDisp)=iDis
             Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
             Call DSCAL_(nConf1,-1.0d0,Work(ipin(ipST)),1)
          End If
*
          Call DSCAL_(nDensC,-1.0d0,Work(ipSigma),1)
*
#ifdef _DEBUG_
          Write (LuWr,*) 'ST=',DDot_(nConf1,Work(ipin(ipST)),1,
     &                                     Work(ipin(ipST)),1)
          Write (LuWr,*) 'Sigma=',DDot_(nDensC,Work(ipSigma),1,
     &                                        Work(ipSigma),1)
          Write (LuWr,*) 'Kap=',DDot_(nDens2,Work(ipKap),1,
     &                                      Work(ipKap),1)
          Write (LuWr,*) 'End RHS!'
#endif
*
          iter=1
          Call DMInvKap(Work(ipin(ipPre2)),Work(ipSigma),nDens2+6,
     &                  Work(ipKap),nDens2+6,work(ipTemp3),nDens2+6,
     &                  isym,iter)
          irc=opout(ippre2)
          r2=ddot_(ndensc,Work(ipKap),1,Work(ipKap),1)
#ifdef _DEBUG_
          Write (LuWr,*) 'DMinvKap'
          Write (LuWr,*) 'Kap=',DDot_(nDens2,Work(ipKap),1,
     &                                      Work(ipKap),1)
          Write (LuWr,*) 'End DMinvKap'
#endif
          If (r2.gt.r1) Write(LuWr,Fmt2//'A,I3,A)')
     &         'Warning  perturbation number ',idisp,' might diverge!'
          Call UnCompress(Work(ipKap),Work(ipdKap),iSym)
*
          If (CI) Call DMinvCI(ipST,Work(ipin(ipCId)), rCHC,isym)
*
          If (CI) Then
             deltaC=ddot_(nConf1,Work(ipin(ipST)),1,Work(ipin(ipCId)),1)
             irc=ipout(ipcid)
          Else
             deltaC=0.0d0
          End If
          deltaK=ddot_(nDensC,Work(ipKap),1,Work(ipSigma),1)
          call dcopy_(nDens,0.0d0,0,Work(ipKap),1)
          delta=deltac+deltaK
#ifdef _DEBUG_
          If (Abs(DeltaC).lt.1.0D-12) DeltaC=0.0D0
          Write (LuWr,*)'DeltaK, DeltaC, Delta=',DeltaK, DeltaC, Delta
#endif

          If (delta.eq.0.0D0) Goto 300
          delta0=delta
#ifdef _DEBUG_
          Write (LuWr,*) 'Delta0=',Delta0
          Write (LuWr,*) 'Start ITERATIONS'
          Write (LuWr,*) 'Orb,CI=',ORB,CI
#endif
          Orb=.true.
          TimeDep=.false.
          ReCo=-1.0d0
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
             Call RInt_generic(Work(ipdKap),Work(iprmoaa),rdum,
     &                         Work(ipSc2),Work(ipTemp3),Work(ipTemp4),
     &                         Work(ipSc3),isym,reco,jspin)
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
     &                       ipTemp4,iprmoaa,idum,
     &                       ipCI,ipS1,'N')
                Clock(iTimeKC)=Clock(iTimeKC)+Tim3
#ifdef _DEBUG_
                Write (LuWr,*) 'CISigma'
                Call RecPrt('CI ','(3F10.4)',Work(ipin(ipCI )),1,nConf1)
                Call RecPrt('S1 ','(3F10.4)',Work(ipin(ipS1 )),1,nConf1)
                Write (LuWr,*) 'CI=',DDot_(nConf1,Work(ipin(ipCI)),1,
     &                                           Work(ipin(ipCI)),1)
                Write (LuWr,*) 'S1=',DDot_(nConf1,Work(ipin(ipS1)),1,
     &                                           Work(ipin(ipS1)),1)
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
                If (isym.eq.1) Then
                   rGrad=DDot_(nconf1,Work(ipin(ipCI)),1,
     &                        Work(ipin(ipS1)),1)
                   call daxpy_(nConf1,-rgrad,
     &                        Work(ipin(ipCI)),1,Work(ipin(ipS1)),1)
                End IF
                call dscal_(nconf1,2.0d0,Work(ipin(ips1)),1)
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
     &                     ipFIMO,k2int,idum,
     &                     ipCId,ipS2,'N')
              EC=rin_ene+potnuc-ERASSCF(1)

              Call DaXpY_(nConf1,EC,Work(ipin(ipCId)),1,
     &                   Work(ipin(ipS2)),1)
              Call DSCAL_(nConf1,2.0d0,Work(ipin(ipS2)),1)
              Clock(iTimeCC)=Clock(iTimeCC)+Tim4
              irc=ipout(ips2)
              irc=opout(ipcid)
              irc=opout(ipci)
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
#ifdef _DEBUG_
              Call RecPrt('CI ','(3F10.4)',Work(ipin(ipCI )),1,nConf1)
              Call RecPrt('CId','(3F10.4)',Work(ipin(ipCId)),1,nConf1)
#endif
              Call CIDens(Response,ipCI,ipCId,
     &                    State_sym,
     &                    PState_Sym,jspin,
     &                    Work(ipP),Work(ipDe))     ! Jeppes

#ifdef _DEBUG_
              Write (LuWr,*) 'After CIDens'
              Write (LuWr,*) 'P=',DDot_(n2dens,Work(ipP),1,
     &                                        Work(ipP),1)
              Write (LuWr,*) 'De=',DDot_(n1dens,Work(ipDe),1,
     &                                         Work(ipDe),1)
#endif
*
*             density for inactive= 2(<d|0>+<0|d>)
*
              d_0=0.0d0
*
*             This is just for debugging purpose.
*             When we use it for actual calculations d_0 == 0
*
              If (isym.eq.1)
     &        d_0=ddot_(nconf1,Work(ipin(ipCid)),1,Work(ipin(ipCi)),1)
              If (Response) d_0=d_0*2.0d0
*

              Call FockGen(d_0,Work(ipDe),Work(ipP),
     &                     Work(ipSc1),Work(ipSc3),isym)    ! Made
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
#ifdef _DEBUG_
           Write (LuWr,*) 'Add together!'
           Write (LuWr,*) 'dKap=',DDot_(nDens,Work(ipdKap),1,
     &                                    Work(ipdKap),1)
           Write (LuWr,*) 'Sc2=',DDot_(nDens,Work(ipSc2),1,
     &                                   Work(ipSc2),1)
           If (CI) Then
              Write (LuWr,*) 'Sc3=',DDot_(nDens,Work(ipSc3),1,
     &                                      Work(ipSc3),1)
             Write(LuWr,*)'S2=',DDot_(nConf1,Work(ipin1(ipS2,nconf1)),1,
     &                                       Work(ipin1(ipS2,nconf1)),1)
             Write(LuWr,*)'S1=',DDot_(nConf1,Work(ipin1(ipS1,nconf1)),1,
     &                                       Work(ipin1(ipS1,nconf1)),1)
           End If
           Write (LuWr,*)
#endif
           irc=ipnout(-1)
           if (CI) then
              Call DZaXpY(nDens,One,Work(ipSc2),1,
     &                    Work(ipSc3),1,Work(ipSc1),1)
           Else
              call dcopy_(nDens,Work(ipSc2),1,Work(ipSc1),1)
           End If
           Call Compress(Work(ipSc1),Work(ipTemp4),isym)   ! ds
           Call Compress(Work(ipdKap),Work(ipTemp2),isym)  ! DX
           If (CI) Then
              Call DaXpY_(nConf1,1.0d0,Work(ipin1(ipS2,nconf1)),1,
     &                                Work(ipin1(ipS1,nconf1)),1)
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
           rAlphaC=0.0d0
           rAlphaK=0.0d0
           If (orb) rAlphaK=DDot_(nDensC,Work(ipTemp4),1,
     &                                  Work(ipTemp2),1)
           If (CI)  rAlphaC=DDot_(nConf1,Work(ipin(ipS1)),1,
     &                                  Work(ipin(ipCId)),1)
           rAlpha=delta/(rAlphaK+rAlphaC)
#ifdef _DEBUG_
           Write (LuWr,*) 'At PCG'
           Write (LuWr,*) 'S1=',
     &                      DDot_(nConf1,Work(ipin(ipS1 )),1,
     &                                  Work(ipin(ipS1 )),1)
           Write (LuWr,*) 'CId=',
     &                      DDot_(nConf1,Work(ipin(ipCId)),1,
     &                                  Work(ipin(ipCId)),1)
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
             Call DaxPy_(nDensC,ralpha,Work(ipTemp2),1,Work(ipKap),1)
             Call DaxPy_(nDensC,-ralpha,Work(ipTemp4),1,Work(ipSigma),1)
              resk=sqrt(ddot_(nDensC,Work(ipSigma),1,Work(ipSigma),1))
           End If
           resci=0.0d0
           If (CI) Then
              Call DaXpY_(nConf1,ralpha,Work(ipin(ipCId)),1,
     &                   Work(ipin(ipCIT)),1)
              irc=ipout(ipcit)
              Call DaXpY_(nConf1,-ralpha,Work(ipin(ipS1)),1,
     &                   Work(ipin1(ipST,nconf1)),1)
              irc=opout(ipS1)
              ip=ipin(ipst)
              resci=sqrt(ddot_(nconf1,Work(ip),1,
     &                        Work(ip),1))
           End If
*                                                                      *
************************************************************************
*                                                                      *
*          Precondition......
*             -1
*          S=M  Sigma
*
          irc=opout(ipcid)
          If (CI) Call DMinvCI(ipST,Work(ipin(ipS2)),rCHC,isym)
          irc=opout(ipci)
          irc=opout(ipdia)
*
           Call DMInvKap(Work(ipin(ipPre2)),work(ipSigma),nDens2+6,
     &                   work(ipSC2),nDens2+6,Work(ipSc1),nDens2+6,
     &                   iSym,iter)
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
              deltaC=ddot_(nConf1,Work(ipin(ipST)),1,Work(ipin(ipS2)),1)
              irc=ipout(ipST)
           Else
              deltaC=0.0d0
           End If
*
           deltaK=ddot_(nDensC,Work(ipSigma),1,Work(ipSc2),1)
           If (.not.CI) Then
              rBeta=deltaK/delta
              delta=deltaK
              Call DScal_(nDensC,rBeta,Work(ipTemp2),1)
              Call DaXpY_(nDensC,1.0d0,Work(ipsc2),1,Work(ipTemp2),1)
           Else
              rbeta=(deltac+deltaK)/delta
              delta=deltac+deltaK
              Call DScal_(nConf1,rBeta,Work(ipin(ipCID)),1)
              Call DScal_(nDensC,rBeta,Work(ipTemp2),1)
              Call DaXpY_(nConf1,1.0d0,Work(ipin(ipS2)),1,
     &                                Work(ipin(ipCID)),1)
              Call DaXpY_(nDensC,1.0d0,Work(ipsc2),1,Work(ipTemp2),1)
              irc=opout(ipS2)
              irc=ipout(ipCID)
           End If
#ifdef _DEBUG_
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
           Call UnCompress(Work(ipTemp2),Work(ipdKap),isym)
*
           res=0.0D0 ! dummy initialize
           res_tmp=-1.0D0
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
*          Call GetMem(' LIST ','LIST','REAL',iDum,iDum)

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
          Call dDaFile(LuTemp,1,Work(ipKap),iLen,iDis)
          iSigDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
          If (CI) Then
             ilen=nconf1
             iCIDisp(iDisp)=iDis
             Call dDaFile(LuTemp,1,Work(ipin(ipCIT)),iLen,iDis)
             iCISigDisp(iDisp)=iDis
             Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
          End If
*
          Call GetMem('Temp4  ','FREE','Real',ipTemp4,nDens)
          Call GetMem('Temp3  ','FREE','Real',ipTemp3,ndens)
          Call GetMem('Temp2  ','FREE','Real',ipTemp2,nDens)
          Call GetMem('Temp1  ','FREE','Real',ipTemp1,nDens)
          Call GetMem('sigma  ','FREE','Real',ipSigma,nDens)
          Call GetMem('dkappa ','FREE','Real',ipdKap ,nDens)
          Call GetMem('kappa  ','FREE','Real',ipKap  ,nDens)
          Call Getmem('Scr1   ','FREE','Real',ipSc1  ,nDens2)
          Call Getmem('Scr2   ','FREE','Real',ipSc2  ,nDens2)
          Call Getmem('Scr3   ','FREE','Real',ipSc3  ,nDens2)
          If (iMethod.eq.2) Then
             Call GetMem('RMOAA','FREE','Real',iprmoaa,n2dens)
          End If
          If (CI) Then
*            Call GetMem('RMOAA','FREE','Real',iprmoaa,n2dens)
             Call GetMem('2Dens','FREE','Real',ipP,n1dens)
             Call GetMem('1Dens','FREE','Real',ipDe,n2dens)
          End If
C        End Do ! jDisp
*
*        Free all memory and remove from disk all data
*        related to this symmetry
*
C        If (CI) Then
C           irc=ipclose(ipdia)
C        Else
C           irc=ipclose(ipPre2)
C        End If
*
C        If (iphx.ne.0) Then
C           Call Getmem('EXPHS','FREE','REAL',iphx,idum)
C           Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
C           Call Getmem('EXPLS','FREE','INTE',iplst,idum)
C        End If
*
C     End Do ! iSym
      Go To 888
 999  Continue
      Call Get_nProcs(nProcs)
      If(nProcs.ge.2) Then
       Write (LuWr,*)
       Write (LuWr,*) ' Perturbations were printed only by master node'
       Write (LuWr,*)
      End If

      Call Free_Tsk(id)
      Call Free_iWork(ipList)
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
      If (iphx.ne.0) Then
         Call Getmem('EXPHS','FREE','REAL',iphx,idum)
         Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
         Call Getmem('EXPLS','FREE','INTE',iplst,idum)
      End If
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
#ifdef _DEBUG_
      Write(LuWr,*)  '****************************************',
     &               '****************************************'
      Write(LuWr,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      If (Fail) Call Quit_OnConvError()
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('WfCtl')
      Return
      End
