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
************************************************************************
      SubRoutine WfCtl_td(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                    iRHSDisp,iRHSCIDISP,converged)
************************************************************************
*                                                                      *
*                                                                      *
*     called from: MCLR                                                *
*                                                                      *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "csfbas_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"

      Logical Orb,CI,Response
      Parameter (iTimeCC = 1 )
      Parameter (iTimeKK = 2 )
      Parameter (iTimeKC = 3 )
      Parameter (iTimeCK = 4 )
#include "crun_mclr.fh"
      Character(LEN=8)   Fmt2
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp),iRHSCIDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      Integer pstate_sym,opout
      Logical lPrint,converged(8)
      Real*8 Clock(4)
      Real*8, Allocatable:: DigPrec(:), Kappa(:), dKappa(:), Sigma(:),
     &                      Temp1(:), Temp2(:), Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Sc3(:), TempTD(:),
     &                      Dens(:), Pens(:), rmoaa(:)
*
*----------------------------------------------------------------------*
*     Start                                                            *
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
*
      fail=.false.
      Converged(:)=.true.
      lprint=.false.
      idasave=0
      LU_50 = 50
      If (SAVE) CALL DANAME(LU_50,'RESIDUALS')
      If (iAnd(kprint,2).eq.2) lprint=.true.
      iDisp=0
      kksym=1
      kkksym=nsym
      If (PT2) kkkSym=1

c
c Starting loop over all symmetries/PT
c
      Do 100 iSym=kksym,kkksym
        PState_SYM=iEor(State_Sym-1,iSym-1)+1
        nconf1=ncsf(PState_Sym)
        CI=.false.
        If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.

        If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.
*          Initiate CSF <-> SD
        If (CI)
     &     Call InCSFSD(iEor(iSym-1,State_Sym-1)+1,
     &                  State_sym,.false.)
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
c
c Setup is realted to symmtry treatment
c
*
        Call Setup_MCLR(iSym)
*
*
*       Determine if we should page CI vectors
*
*                                    [2]
*         Calculate the diagonal of E    and store in core/disc
*
        iphx=0
        If (CI) Then
c
c CIDia calculates the <i|H|i> elements (diagonal)? From the CI-CI part of E?
c
            Call CIDia_td(PState_Sym)
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
c      ipget is like getmem but handles weather the vector is on disk or not.
c
        ips1 =ipget(2*nconf3)
        ips2 =ipget(2*nconf3)
        ipst =ipget(2*nconf3)
        ipcit=ipget(2*nconf1)
        ipcid=ipget(2*nconf1)
C
        End If
*
*
       npre2=npre(isym)
c
*#ifndef DEBUG
C
C OBS npre2 not def
C
       ipPre2=ipget(npre2)
*
*
        call Prec_dig(Work(ipin(ipPre2)),isym)
*
*
        Call mma_allocate(DigPrec,nDensC,Label='DigPrec')
        DigPrec(:)=Zero
        Call Prec_td(Work(ipin(ipPre2)),DigPrec,isym)
*
       irc=ipout(ippre2)
*#endif
*
*      OK START WORKING
*
        iDEND=lDisp(iSym)
*       If (SewLab.eq.'NONE    '.and.(.not.mckinley)) iDEND=1
c
c Loop over all PT of sym isym
c
        Do 110 jDisp=1,iDEnd
          iDisp=iDisp+1
          jspin=0
          If (iAnd(nTPert(idisp),1).eq.1) jSpin=1
          if (jspin.eq.0) Then
           nconf1=ncsf(PState_Sym)
          else
           nConf1=nint(xispsm(Pstate_Sym,1))
          end if
          If (.not.lCalc(iDisp)) Then
             converged(isym)=.false.
             Goto 110
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
             Call mma_allocate(rmoaa,n2dens,Label='Pens')
             Dens(:)=Zero
             Pens(:)=Zero
             rmoaa(:)=Zero
          End If
*
*-----------------------------------------------------------------------------
*
*    Calculate RHS for the perturbation. The b^x vector!
*    sigma is the orbital part and Temp4 is the cipart
*
*-----------------------------------------------------------------------------
*
*                 (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)
*
           Call RHS_td(Sigma,Kappa,Temp1,Temp3,Sc2,dKappa,
     &                 Sc3,Temp4,ipST,
     &                 iDisp,iSym-1,Work(ipCMO),jdisp,jspin,CI)
*
*
           call dscal_(nDens2,-1.0d0,Temp4,1)
c
c Make RHS twice as long and change sign on second part!
c
           If (CI) Then
               call dcopy_(nConf1,Work(ipin(ipst)),1,
     &               Work(ipin(ipst)+nConf1),1)
               call dscal_(nConf1,-1.0d0,Work(ipin(ipst)+nConf1),1)
           End If
C
          irc=opout(ipci)
*
          If (lprint) Write(6,*)
     &      '       Iteration         Delta     Res(kappa) Res(CI)'
          iLen=nDensC
          iRHSDisp(iDisp)=iDis
          Call Compress(Temp4,Sigma,iSym)
          r1 = 0.50d0*ddot_(ndensc,Sigma,1,Sigma,1)
*
          Call UnCompress(Sigma,Temp4,iSym)
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (CI)
     &    call dcopy_(2*nConf1,[0.0d0],0,Work(ipin(ipCIT)),1)
          irc=ipout(ipcit)
          If (CI) Then
            ilen=2*nconf1
            iRHSCIDisp(iDisp)=iDis
            Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
            Call DSCAL_(2*nConf1,-1.0d0,Work(ipin(ipST)),1)
          End If
*
         Call DMInvKap_td(DigPrec,Sigma,Kappa)
C
          irc=opout(ippre2)
*         kap:kap
          r2 = 0.50d0*ddot_(ndensc,Kappa,1,Kappa,1)
          If (r2.gt.r1) Write(6,*) 'Warning ',
     &     ' perturbation number ',idisp,
     &     ' might diverge'
*
*         dkap=Kap in matrix form
          Call UnCompress(Kappa,dKappa,iSym)
c
c DMinvCI is related to cidia and inverts the precond CI-CI
c part of E
c Has to be modified: <i|H|i> --> <i|H|i>+w
c

          If (CI) Then
                Call DMinvCI_td(Work(ipin(ipST)),
     &                          Work(ipin(ipCid)), -omega,isym)
                Call DMinvCI_td(work(ipin(ipST)+nConf1),
     &                          work(ipin(ipCId)+nconf1),
     &                          omega,isym)
C
          End if

*
          If (CI) Then
             deltaC= 0.50d0*ddot_(2*nConf1,Work(ipin(ipST)),1,
     &                   Work(ipin(ipCId)),1)
             irc=ipout(ipcid)
          Else
            deltac=0.0d0
          End If
          deltaK= 0.50d0*ddot_(nDensC,Kappa,1,Sigma,1)
          Kappa(1:nDens)=Zero
          delta=deltac+deltaK

          If (delta.eq.0.0D0) Goto 300

          delta0=delta
          Orb=.true.
          ReCo=-1.0d0
          iter=1
*-----------------------------------------------------------------------------
*
*
***********************************************************
*          I   T   E   R   A   T   I   O   N   S          *
***********************************************************
*
*
200       Continue
*
*
**********************************************************************
*
*       O R B I T A L    P A R T of the trail vector
*
**********************************************************************
*



           If (orb) Then

*-----------------------------------------------------------------------------
*
*                      ~    ~
*            Construct F,(ij|kl)
*
*-----------------------------------------------------------------------------
c
c j2 specifies which part of E I want to look at
c j2=0 --> K-K, j2=-1 --> CI-CI, These are antisym within themself
c j2>0 --> CI-K and K-CI, These parts are antisym between eachother
c

             irc=ipnout(-1)
             Call RInt_ns(dKappa,
     &                 rmoaa,    ! OIT 2-el-int (active indexes)
     &                 Sc2,      ! SC2 contains the E*kappa
     &                 Temp4,    ! Contains OIT FI
     &                 isym,reco,jspin,rInEne) ! OIT integrals are used
*
             Call RInttd(Sc2,dKappa,isym)
c
             Clock(iTimeKK)=Clock(iTimeKK)+Tim2
*
*
*-----------------------------------------------------------------------------
*
*            kappa->CI
*
*            H(kappa)|0>
*
*                [2]
*            S1=E   k ( kappa TO CI ) <i|H|0>-<i|0>*Energy
*
*-----------------------------------------------------------------------------
*
c
c This cisigma call gives <j|H(k)|0> and <j|H(k)t|0>
c
             If (CI) Then
*              Adjusted to timedep
               Call CISigma_td(jspin,State_Sym,pstate_sym,
     &                      Temp4,rmoaa,idum,
     &                      ipCI,ipS1,'T')
               Clock(iTimeKC)=Clock(iTimeKC)+Tim3
*
*              This will give us a better
*              convergence in the PCG. Notice that
*
*              ~Inactive     ~
*              E        + <0|H|0> = 0
*
*              when the wavefunction is converged.
*
c
c These terms are to be able to handel less converged CASSCF wave func
c
               If (isym.eq.1) Then
                 rGrad=ddot_(nconf1,Work(ipin(ipCI)),1,
     &                   Work(ipin(ips1)),1)
                 call daxpy_(nConf1,-rGrad,
     &                    Work(ipin(ipCI)),1,Work(ipin(ipS1)),1)
                 rGrad=ddot_(nconf1,Work(ipin(ipCI)),1,
     &                   Work(ipin(ips1)+nconf1),1)
                 call daxpy_(nConf1,-rGrad,
     &                     Work(ipin(ipCI)),1,Work(ipin(ipS1)+nconf1),1)
               End if
                call dscal_(nconf1,-1.0d0,Work(ipin(ips1)),1)
                call dscal_(2*nconf1,2.0d0,Work(ipin(ips1)),1)
C
C

               irc=opout(ipCI)
************************************************************************
*
*
************************************************************************
*
             End If  ! If ci


           End If
*
************************************************************************
*
*        C I    P A R T of the trail vector
*
************************************************************************
*
           If (CI) Then


*-----------------------------------------------------------------------------
*
*                  [2]
*            S2 = E   CID  ( CI TO CI ) = <i|H|d> -E<i|d>
*
*-----------------------------------------------------------------------------
*
             irc=ipnout(-1)
             if (CI) Call CISigma_td(0,PState_Sym,Pstate_sym,
     &                    Work(ipFIMO),Work(k2int),idum,
     &                    ipCId,ipS2,'S')
*
c I want the RASSCF energy of the ACTIVE electrons !!!!
c EC=-E[act]           E[RASSCF]=E[inact]+E[act]+E[nuc]
c
             EC=rin_ene+potnuc-ERASSCF(1)
*
                Call DaXpY_(nConf1,EC,Work(ipin(ipCId)),1,
     &                  Work(ipin(ipS2)),1)
                Call DaXpY_(nConf1,EC,Work(ipin(ipCId)+nConf1),1,
     &                  Work(ipin(ipS2)+nConf1),1)
                call dscal_(2*nConf1,2.0d0,Work(ipin(ipS2)),1)
c
c Add the wS contribution
c The (-) sign in both daxpys assumes that the two parts of ipcid are def with diff sign.
c This is not true for the debug option!! ipcid = 1 regardless of part which part
c The S-contribution will make E-wS loose its symmetry because E is sym and S
c is antisym.
c
             Call DaXpY_(nConf1,-2.0d0*omega,Work(ipin(ipCId)),1,
     &                  Work(ipin(ipS2)),1)
             Call DaXpY_(nConf1,2.0d0*omega,
     &            Work(ipin(ipCId)+nConf1),1,Work(ipin(ipS2)+nConf1),1)
             Clock(iTimeCC)=Clock(iTimeCC)+Tim4
*
             irc=ipout(ips2)
             irc=opout(ipcid)
             irc=opout(ipci)
*
*-----------------------------------------------------------------------------
*
*            CI -> Kappa
*
*
*                 [2]
*            SC3=E   CID   (SC3=F(<d|E|0>+<0|E|d>)
*
*-----------------------------------------------------------------------------
*
             Response=.true.
             irc=ipnout(-1)
c
             Call CIDens_TD(ipCid,PState_Sym,Pens,Dens)     ! Jeppes
*
*            density for inactive= 2(<d|0>+<0|d>)
*
             d_0=0.0d0
*
*            This is just for debugging purpose.
*            When we use it for actual calculations d_0 == 0
*
c This if statement is just for better convergence! Grad term
c Leave this for later!
c
             If (isym.eq.1) Then
                d_1=ddot_(nconf1,Work(ipin(ipCid)),1,
     &                  Work(ipin(ipci)),1)
                d_2=ddot_(nconf1,Work(ipin(ipCid)+nConf1),1,
     &                  Work(ipin(ipci)),1)
                d_0 = d_1 + d_2
             End If
c
c Fockgen gives the Fock matrix, MO integrals and one index transformed
c MO integrals
c
                 Call Fockgen_td(d_0,Dens,Pens,Sc3,isym)
*
*-----------------------------------------------------------------------------
*
           End If
*
**********************************************************************
*
*        Sc1  kappa-> kappa
*        Sc3  CI -> kappa
*        S1   kappa -> CI
*        S2   CI -> CI
*        dKap present step
*        Kap  kappaX
*        CIT  CIX
*        CId  present step
**********************************************************************
*
*        Add together
*
*
**********************************************************************
C
           irc=ipnout(-1)
           if (CI) then   ! if (.false.) then
              Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)
           Else
              call dcopy_(nDens,Sc2,1,Sc1,1)
           End If
           Call Compress(Sc1,Temp4,isym)   ! ds
           Call Compress(dKappa,Temp2,isym) ! DX
c
c S1 + S2 --> S1
c
           If (CI) Then  !If (.false.) then
                     Call DaXpY_(2*nConf1,1.0d0,
     &               Work(ipin1(ipS2,2*nconf1)),1,
     &               Work(ipin1(ipS1,2*nconf1)),1)
                     irc=opout(ips2)
           End If
*
*-----------------------------------------------------------------------------
*
*                ######   #####   #####
*                #     # #     # #     #
*                #     # #       #
*                ######  #       #  ####
*                #       #       #     #
*                #       #     # #     #
*                #        #####   #####
*
*-----------------------------------------------------------------------------
**********************************************************************
*
*
*                     delta
*          rAlpha=------------
*                 dKappa:dSigma
*
*-----------------------------------------------------------------------------
           rAlphaC=0.0d0
           rAlphaK=0.0d0
           If (orb)
     &      rAlphaK=0.5d0*ddot_(nDensC,Temp4,1,Temp2,1)
           If (CI) Then
              rAlphaC=0.5d0*ddot_(2*nConf1,Work(ipin(ipS1)),1,
     &                   Work(ipin(ipCId)),1)
           End If
           rAlpha=delta/(rAlphaK+ralphaC)
*
*-------------------------------------------------------------------*
*
*          Kappa=Kappa+rAlpha*dKappa
*          Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
*
           If (orb) Then
            Call DaxPy_(nDensC,ralpha,Temp2,1,Kappa,1)
            Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
            resk=sqrt(0.5d0*ddot_(nDensC,Sigma,1,Sigma,1))
           End If
           resci=0.0d0
C
           If (CI) Then
              Call DaXpY_(2*nConf1,ralpha,Work(ipin(ipCId)),1,
     &                  Work(ipin(ipCIT)),1)
              irc=ipout(ipcit)
              Call DaXpY_(2*nConf1,-ralpha,Work(ipin(ipS1)),1,
     &                   Work(ipin1(ipST,2*nconf1)),1)
              irc=opout(ipS1)
              ip=ipin(ipst)
              resci=sqrt(0.5d0*ddot_(2*nconf1,Work(ip),1,
     &                        Work(ip),1))
           End If
*
*-------------------------------------------------------------------*
*
*          Precondition......
*             -1
*          S=M  Sigma
*
          irc=opout(ipcid)
          If (CI) Then
             Call DMinvCI_td(work(ipin(ipST)),
     &                          Work(ipin(ipS2)), -omega,isym)
             Call DMinvCI_td(work(ipin(ipST)+nConf1),
     &                          Work(ipin(ipS2)+nconf1),omega,isym)
C
          End if
          irc=opout(ipci)
          irc=opout(ipdia)
*
          Call DMInvKap_td(DigPrec,Sigma,Sc2)
C
           irc=opout(ippre2)
*
*-------------------------------------------------------------------*
*               s:Sigma
*          Beta=-------
*                delta
*
*          delta=s:sigma
*
*          dKappa=s+Beta*dKappa
*
           If (CI) Then
               deltaC=0.50d0*ddot_(2*nConf1,Work(ipin(ipST)),1,
     &                     Work(ipin(ipS2)),1)
           irc=ipout(ipST)
           else
           deltaC=0.0d0
           end if
*
           deltaK=0.50d0*ddot_(nDensC,Sigma,1,Sc2,1)
           If (.not.CI) Then
             rBeta=deltaK/delta
             delta=deltaK
             Call DScal_(nDensC,rBeta,Temp2,1)
             Call DaXpY_(nDensC,1.0d0,sc2,1,Temp2,1)
           Else
             rbeta=(deltac+deltaK)/delta
             delta=deltac+deltaK
             Call DScal_(2*nConf1,rBeta,Work(ipin(ipCID)),1)
             Call DScal_(nDensC,rBeta,Temp2,1)
             Call DaXpY_(2*nConf1,1.0d0,Work(ipin(ipS2)),1,
     &                             Work(ipin(ipCID)),1)
             Call DaXpY_(nDensC,1.0d0,sc2,1,Temp2,1)
             irc=opout(ipS2)
             irc=ipout(ipCID)
           End If
*

*    ######  #    #  #####        #####    ####    ####
*    #       ##   #  #    #       #    #  #    #  #    #
*    #####   # #  #  #    #       #    #  #       #
*    #       #  # #  #    #       #####   #       #  ###
*    #       #   ##  #    #       #       #    #  #    #
*    ######  #    #  #####        #        ####    ####
*


*
*-------------------------------------------------------------------*
*
           Call UnCompress(Temp2,dKappa,isym)
*
c
c iBreak is defined via include!
c
           res=0.0D0 ! dummy initialize
           If (iBreak.eq.1) Then
*             This is the actual breaking!
              If (abs(delta).lt.abs(Epsilon**2*delta0)) Goto 300
           Else If (ibreak.eq.2) Then
              res=sqrt(resk**2+resci**2)
              If (res.lt.abs(epsilon)) Goto 300
           Else
              If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &            res.lt.abs(epsilon))  Goto 300
           End If
c
c This breaks the PCG iterations by going to 210
c
           If (iter.ge.niter) goto 210
           If (lprint)
     &     Write(6,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)')
     &            '     ',
     &            iter,'       ',delta/delta0,resk,resci,deltac,deltak

           iter=iter+1

          Goto 200
*
**********************************************************************
*
 210      Continue
          Write(6,Fmt2//'A,I4,A)')
     &    'No convergence for perturbation no: ',
     &                      idisp,'. Increase Iter.'
          converged(isym)=.false.
          fail=.true.
          Goto 310
 300      Write(6,Fmt2//'A,I4,A,I4,A)')
     &           'Perturbation no: ',idisp,' converged in ',
     &             iter-1,' steps.'
          irc=ipnout(-1)
*          stop 10
 310      Continue
          Call mma_allocate(TempTD,nDens2,Label='TempTD')
          Call Uncompress(Kappa,TempTD,isym)
          Call mma_deallocate(TempTD)
C
          Write(6,*)
          iLen=ndensC
          iKapDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Kappa,iLen,iDis)
          iSigDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (CI) Then
            ilen=2*nconf1
            iCIDisp(iDisp)=iDis

            Call dDaFile(LuTemp,1,Work(ipin(ipCIT)),iLen,iDis)
            iCISigDisp(iDisp)=iDis
            Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
          End If
*
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
          If (CI) Then
             Call mma_deallocate(rmoaa)
             Call mma_deallocate(Pens)
             Call mma_deallocate(Dens)
          End If
 110     Continue
*
*        Free all memory and remove from disk all data
*        related to this symmetry
*
         If (CI) Then
            irc=ipclose(ipdia)
         Else
            irc=ipclose(ipPre2)
         End If
*
         Call mma_deallocate(DigPrec)
         If (iphx.ne.0) Then
            Call Getmem('EXPHS','FREE','REAL',iphx,idum)
            Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
            Call Getmem('EXPLS','FREE','INTE',iplst,idum)
         End If

 100  Continue
      If (debug) Then
      Write(6,*)  '****************************************',
     &            '****************************************'
      Write(6,*)
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
