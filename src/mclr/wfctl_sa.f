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
      SubRoutine WfCtl_SA(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                    iRHSDisp,converged,iPL)
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

#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "csfbas_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"

      Logical CI
#include "crun_mclr.fh"
      Character*8   Fmt2
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      Integer opOut
      Logical lPrint,converged(8)
      Real*8 rchc(mxroot)
*
      Call QEnter('WfCtl_SA')
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call StatusLine(' MCLR:',
     &                ' Computing Lagrangian multipliers for SA-CASSCF')
*

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
*MGD I think this is nice when printed...
      lprint=.true.
      idasave=0
      reco=-One
      Lu_50=50
      If (SAVE) CALL DANAME(Lu_50,'RESIDUALS')
      If (SAVE) Then
         Write (6,*) 'WfCtl_SA: SAVE option not implemented'
         Call Abend()
      End If
      If (iAnd(kprint,2).eq.2) lprint=.true.
      isym=1
      nconf1=ncsf(State_Sym)

      CI=.false.
      If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.

*          Initiate CSF <-> SD
           Call InCSFSD(iEor(iSym-1,State_Sym-1)+1,
     &                  State_sym,.false.)
*
*
*          Calculate length of the density, Fock and Kappa matrix etc
*          notice that this matrices are not necessarily symmetric.
*          Store pointers.
*
*          Input:
*                 iSym: Symmetry of perturbation
*
*          Output: Commonblocks (Pointers.fh)
*
      nConf2=nint(xispsm(State_SYM,1))
      nConf3=nint(Max(xispsm(State_SYM,1),xispsm(State_SYM,1)))

      Call Setup_MCLR(iSym)

*
*     Determine if we should page CI vectors
*                                [2]
*     Calculate the diagonal of E    and store in core/disc
*
      iphx=0
      Call Getmem('FANCY','ALLO','REAL',ips,nroots**3)
      Call CIDia_SA(State_Sym,rCHC,Work(ipS))

      irc=ipOut(ipdia)
*
*     Allocate disk/memory space
*
*
*     This areas should be addressed through ipIn
*     ipOut will page them out to disk and free the memory area
*     if page mode is used
*
*     opOut will release the memory area without update the disk
*
      ipS1 =ipGet(nconf3*nroots)
      ipS2 =ipGet(nconf3*nroots)
      ipST =ipGet(nconf3*nroots)
      ipCIT=ipGet(nconf1*nroots)
      ipCId=ipGet(nconf1*nroots)
*
      npre2=npre(isym)
      ipPre2=ipGet(npre2)

      If (TwoStep.and.(StepType.eq.'RUN2')) Then
         ! fetch data from LuQDAT and skip the call to "Prec"
         Call ddafile(LuQDAT,2,Work(ipIn(ipPre2)),npre2,iaddressQDAT)
      Else
         Call Prec(Work(ipIn(ipPre2)),isym)
         irc=ipOut(ippre2)
         If (TwoStep.and.(StepType.eq.'RUN1')) Then
            ! save the computed data in "Prec" to LuQDAT and skip the
            ! following part of this function
            Call ddafile(LuQDAT,1,Work(ipIn(ipPre2)),npre2,iaddressQDAT)
            Go To 193
         End If
      End If
*
*
*     OK START WORKING
*
*     idisp=1
      jspin=0
*
*     Allocate areas for scratch and state variables
*
      Call GetMem('kappa ','Allo','Real',ipKap  ,nDens2+6)
      Call GetMem('dkappa','Allo','Real',ipdKap ,nDens2+6)
      Call GetMem('sigma ','Allo','Real',ipSigma,nDens2+6)
      Call GetMem('Temp3 ','ALLO','Real',ipTemp3,nDens2+6)
      Call GetMem('Temp4 ','Allo','Real',ipTemp4,nDens2+6)
      Call Getmem('Scr1  ','ALLO','Real',ipSc1  ,nDens2+6)
      Call Getmem('Scr2  ','ALLO','Real',ipSc2  ,nDens2+6)
*
      do iDisp=1,nDisp
      call dcopy_(nDens2,[Zero],0,Work(ipKap),1)
      call dcopy_(nDens2,[Zero],0,Work(ipSigma),1)
      call dcopy_(nDens2,[Zero],0,Work(ipdKap),1)
*
*-----------------------------------------------------------------------------
*
*     Calculate RHS for the perturbation
*
*-----------------------------------------------------------------------------
*
*     (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)
*
      If (debug) Then
         If (isNAC) Then
            Write(6,*)'States: ',NACstates(1),NACstates(2)
         Else
            Write(6,*)'State: ',irlxroot
         EndIf
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If (isNAC) Then
        Call RHS_NAC(Work(ipTemp4))
      Else
        Call RHS_SA(Work(ipTemp4))
      End If
*
      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta           Res(kappa)     '//
     &      '  Res(CI)          DeltaK           DeltaC'
      iLen=nDensC
      iRHSDisp(iDisp)=iDis
      Call Compress(Work(ipTemp4),Work(ipSigma),iSym)
*     Call RecPrt('RHS',' ',Work(ipSigma),nDensc,1)
      r1=ddot_(nDensc,Work(ipSigma),1,Work(ipSigma),1)
      If(debug)Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
*
      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCIT)),1)
      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipST)),1)
      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCID)),1)
      irc=ipOut(ipCIT)
      Call DSCAL_(nDensC,-One,Work(ipSigma),1)
*
      Call DMInvKap(Work(ipIn(ipPre2)),Work(ipSigma),nDens2+6,
     &              Work(ipKap),nDens2+6,work(ipTemp3),nDens2+6,
     &              isym,iter)

      irc=opOut(ippre2)
      r2=ddot_(ndensc,Work(ipKap),1,Work(ipKap),1)
      If(debug)Write(6,*) 'In that case I think that r2 should be:',r2
      If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
      call dcopy_(ndensC,Work(ipKap),1,Work(ipdKap),1)
*
      deltaC=Zero
      irc=ipOut(ipcid)
      deltaK=ddot_(nDensC,Work(ipKap),1,Work(ipSigma),1)
      call dcopy_(nDens,[Zero],0,Work(ipKap),1)
      delta=deltac+deltaK
      delta0=delta
      iter=1
      If (delta.eq.Zero) Goto 300
*-----------------------------------------------------------------------------
*
200   Continue
*
         Call TimesE2(Work(ipdKap),ipCId,1,reco,jspin,ipS2,
     &                Work(ipTemp4),ipS1)
*
*-----------------------------------------------------------------------------
*
*                   delta
*        rAlpha=------------
*               dKappa:dSigma
*
*-----------------------------------------------------------------------------
*
         rAlphaK=Zero
         rAlphaK=ddot_(nDensC,Work(ipTemp4),1,Work(ipdKap),1)
         rAlphaC=Zero
         rAlphaC=ddot_(nConf1*nroots,Work(ipIn(ipS1)),1,
     &                              Work(ipIn(ipCId)),1)
         rAlpha=delta/(rAlphaK+rAlphaC)
*
*-------------------------------------------------------------------*
*
*        Kappa=Kappa+rAlpha*dKappa
         Call DaxPy_(nDensC,ralpha,Work(ipdKap),1,Work(ipKap),1)
*        Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
         Call DaxPy_(nDensC,-ralpha,Work(ipTemp4),1,Work(ipSigma),1)
         resk=sqrt(ddot_(nDensC,Work(ipSigma),1,Work(ipSigma),1))
         resci=Zero
         Call DaXpY_(nConf1*nroots,ralpha,Work(ipIn(ipCId)),1,
     &                                   Work(ipIn(ipCIT)),1)
         irc=ipOut(ipcit)
*        ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
         Call DaXpY_(nConf1*nroots,-ralpha,Work(ipIn(ipS1)),1,
     &                                    Work(ipIn(ipST)),1)
         irc=opOut(ipS1)
         ip=ipIn(ipst)
         resci=sqrt(ddot_(nconf1*nroots,Work(ip),1,
     &                                 Work(ip),1))

*-------------------------------------------------------------------*
*
*        Precondition......
*           -1
*        S=M  Sigma
*
         irc=opOut(ipcid)

         Call DMinvCI_SA(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
         irc=opOut(ipci)
         irc=opOut(ipdia)

         Call DMInvKap(Work(ipIn(ipPre2)),Work(ipSigma),nDens2+6,
     &                 Work(ipSc2),nDens2+6,Work(ipSc1),nDens2+6,
     &                 iSym,iter)
         irc=opOut(ippre2)
*
*-------------------------------------------------------------------*
*             s:Sigma (k+1)     s:Sigma (k+1)
*        Beta=-------        =  -------------
*              delta  (k)        s:Sigma (k)
*
*        delta=s:sigma
*
*        dKappa=s+Beta*dKappa
*
         deltaC=ddot_(nConf1*nroots,Work(ipIn(ipST)),1,
     &                             Work(ipIn(ipS2)),1)
         irc=ipOut(ipST)
*
         deltaK=ddot_(nDensC,Work(ipSigma),1,Work(ipSc2),1)
         If (.not.CI) Then
            rBeta=deltaK/delta
            delta=deltaK
            Call DScal_(nDensC,rBeta,Work(ipdKap),1)
            Call DaXpY_(nDensC,One,Work(ipSc2),1,Work(ipdKap),1)
         Else
            rbeta=(deltac+deltaK)/delta
            delta=deltac+deltaK
            Call DScal_(nConf1*nroots,rBeta,Work(ipIn(ipCID)),1)
            Call DScal_(nDensC,rBeta,Work(ipdKap),1)
            Call DaXpY_(nConf1*nroots,One,Work(ipIn(ipS2)),1,
     &                                     Work(ipIn(ipCID)),1)
            Call DaXpY_(nDensC,One,Work(ipSc2),1,Work(ipdKap),1)
            irc=opOut(ipS2)
            irc=ipOut(ipCID)
         End If

*    ######  #    #  #####        #####    ####    ####
*    #       ##   #  #    #       #    #  #    #  #    #
*    #####   # #  #  #    #       #    #  #       #
*    #       #  # #  #    #       #####   #       #  ###
*    #       #   ##  #    #       #       #    #  #    #
*    ######  #    #  #####        #        ####    ####
*
*-------------------------------------------------------------------*
*
*
         res=Zero ! dummy initialize
         If (iBreak.eq.1) Then
            If (abs(delta).lt.abs(Epsilon**2*delta0)) Goto 300
         Else If (iBreak.eq.2) Then
            res=sqrt(resk**2+resci**2)
            if (doDMRG) res=sqrt(resk**2)
            If (res.lt.abs(epsilon)) Goto 300
         Else
            If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &          res.lt.abs(epsilon))  Goto 300
         End If
         If (iter.ge.niter) goto 210
         If (lprint)
     &   Write(6,Fmt2//'I7,4X,ES17.9,ES17.9,ES17.9,ES17.9,ES17.9)')
     &          iter,delta/delta0,resk,resci,deltac,deltak
         iter=iter+1
*
         Goto 200
*
**********************************************************************
*
 210     Continue
         Write(6,Fmt2//'A,I4,A)')
     &         'No convergence for perturbation no: ',
     &          idisp,'. Increase Iter.'
         converged(isym)=.false.
         fail=.true.
         Goto 310
 300  Continue
      If (iPL.ge.2) Then
          Write(6,Fmt2//'A,I4,A,I4,A)')
     &          'Perturbation no: ',idisp,' converged in ',
     &          iter-1,' steps.'
      End If
      irc=ipnout(-1)
*
 310  Continue
      If (iPL.ge.2) Write(6,*)
      iLen=ndensC
      iKapDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Work(ipKap),iLen,iDis)
      iSigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
      ilen=nconf1*nroots
      iCIDisp(iDisp)=iDis
*
      Call dDaFile(LuTemp,1,Work(ipin(ipCIT)),iLen,iDis)
*
**MGD This last call seems unused, so I comment it
*
*      Call TimesE2(Work(ipKap),ipCIT,1,reco,jspin,ipS2,
*     &             Work(ipTemp4),ipS2)
      iCISigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
      end do ! iDisp
*
      Call Getmem('Scr2   ','FREE','Real',ipSc2  ,nDens2)
      Call Getmem('Scr1   ','FREE','Real',ipSc1  ,nDens2)
      Call GetMem('Temp4  ','FREE','Real',ipTemp4,nDens)
      Call GetMem('Temp3  ','FREE','Real',ipTemp3,ndens)
      Call GetMem('sigma  ','FREE','Real',ipSigma,nDens)
      Call GetMem('dkappa ','FREE','Real',ipdKap ,nDens)
      Call GetMem('kappa  ','FREE','Real',ipKap  ,nDens)
*
*     Free all memory and remove from disk all data
*     related to this symmetry
*
193   Continue
      Call Getmem('FANCY',  'FREE','REAL',ips,nroots**3)

      irc=ipclose(ipdia)
      If (.not.CI) irc=ipclose(ipPre2)
*
*     Call GetMem('PREC','FREE','Real',ipPRE,nDens2)
      If (iphx.ne.0) Then
         Call Getmem('EXPHS','FREE','REAL',iphx,idum)
         Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
         Call Getmem('EXPLS','FREE','INTE',iplst,idum)
      End If

      If (debug) Then
      Write(6,*)  '****************************************'//
     &            '****************************************'
      Write(6,*)
      End If
      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Call QExit('WfCtl_SA')
      Return
      End

      Subroutine TimesE2(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)

      Implicit Real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "dmrginfo_mclr.fh"
#include "real.fh"
#include "Input.fh"
      Integer opOut
      Real*8 Kap(*),KapOut(*)
      Dimension rdum(1)
*
      Call GetMem('RMOAA','ALLO','REAL',iprmoaa,n2dens)
      Call GetMem('SCR2','ALLO','REAL',ipSc2,ndens2)
      Call GetMem('SCR1','ALLO','REAL',ipSc1,ndens2)
      Call GetMem('SCR3','ALLO','REAL',ipSc3,ndens2)
      Call GetMem('SCR4','ALLO','REAL',ipTemp4,ndens2)
      Call GetMem('SCR3','ALLO','REAL',ipTemp3,ndens2)
*
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      Call Uncompress(Kap,Work(ipSC1),isym)

! Integral derivative !yma
      Call RInt_generic(Work(ipSC1),Work(iprmoaa),rdum,
     &                 Work(ipSc2),
     &                 Work(ipTemp3),Work(ipTemp4),Work(ipSc3),
     &                 isym,reco,jspin)

      Call Kap_CI(ipTemp4,iprmoaa,ipCIOUT)
      Call Ci_Ci(ipcid,ipS2)
      Call CI_KAP(ipCid,Work(ipSc1),Work(ipSc3),isym)

      Call DZaXpY(nDens,One,Work(ipSc2),1,
     &            Work(ipSc3),1,Work(ipSc1),1)
*
      Call Compress(Work(ipSc1),KapOut,isym)   ! ds
*     Call RecPrt('Ex',' ',KapOut,ndensC,1)
*
      Call DaXpY_(nConf1*nroots,One,
     &               Work(ipin(ipS2)),1,
     &               Work(ipin(ipCIOUT)),1)
      irc=opOut(ipCId)

*
      Call GetMem('RMOAA','FREE','REAL',iprmoaa,n2dens)
      Call GetMem('SCR2','FREE','REAL',ipSc2,ndens2)
      Call GetMem('SCR1','FREE','REAL',ipSc1,ndens2)
      Call GetMem('SCR3','FREE','REAL',ipSc3,ndens2)
      Call GetMem('SCR4','FREE','REAL',ipTemp4,ndens2)
      Call GetMem('SCR5','FREE','REAL',ipTemp3,ndens2)

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*
      Return
      End
