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
      SubRoutine WfCtl_sp(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                    iRHSDisp,iRHSCIDISP)
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
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "spin.fh"
#include "crun_mclr.fh"
      Character*8   Fmt2
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp),iRHSCIDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      integer opout
      Logical lPrint
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      one=1.0d0
*
*
      lPaper=132
      lLine =120
      left=(lPaper-lLine)/2
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

      fail=.false.
      idis=0
      lprint=.false.
      nconf1=0

      If (iAnd(kprint,2).eq.2) lprint=.true.
      If (iMethod.eq.2)
     &     Call InCSFSD(State_Sym,
     &                  State_sym,.false.)
*
        nconf1=ncsf(State_Sym)
        nConf2=nint(xispsm(State_SYM,1))
        nconf3=nint(xispsm(State_SYM,1))
        Call Setup_MCLR(1)
*
*                                    [2]
*         Calculate the diagonal of E    and store in core/disc
*
        If (imethod.gt.0) Then
        If (nconf1.gt.1)
     &       Call CIDia(State_Sym,rCHC)
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

        End If
*
*
        idisp=1
*
*      OK START WORKING
*
*
*    Allocate areas for scratch and state variables
*
          Call GetMem('kappa ','Allo','Real',ipKap  ,nDens2)
          Call GetMem('SpinFock','Allo','Real',ipFS ,nDens2)
          Call GetMem('dkappa','Allo','Real',ipdKap ,nDens2)
          Call GetMem('sigma ','Allo','Real',ipSigma,nDens2)
          Call GetMem('Temp1 ','Allo','Real',ipTemp1,nDens2)
          Call GetMem('Temp2 ','Allo','Real',ipTemp2,nDens2)
          Call GetMem('Temp3 ','ALLO','Real',ipTemp3,nDens2)
          Call GetMem('Temp4 ','Allo','Real',ipTemp4,nDens2)
          Call Getmem('Scr1  ','ALLO','Real',ipSc1  ,nDens2)
          Call Getmem('Scr3  ','ALLO','Real',ipSc3  ,nDens2)
          Call Getmem('Scr2  ','ALLO','Real',ipSc2  ,nDens2)
          Call Getmem('PRE  ','ALLO','Real',ipPre2  ,nDensC)
          call dcopy_(nDens2,0.0d0,0,Work(ipTemp1),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipKap),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipsigma),1)
          call dcopy_(nDens2,0.0d0,0,Work(ipdKap),1)
          If (iMethod.eq.2) Then
             Call GetMem('1Dens','ALLO','Real',ipDe,n1dens)
             Call GetMem('2Dens','ALLO','Real',ipP,nna**4)
             Call GetMem('2Dens','ALLO','Real',iprmoaa,nna**4)
             Call GetMem('2Dens','ALLO','Real',iprmoaa2,nna**4)
          End If
*
*-----------------------------------------------------------------------------
*
*    Calculate RHS
*
*-----------------------------------------------------------------------------
*
          Call Pre_SP(work(ippre2),1)
          Call FockGen_sp(0.0d0,Work(ipG1m),Work(ipG2mp),
     &                   Work(ipFS),Work(ipTemp4),1)
          call dcopy_(nconf1,0.0d0,0,Work(ipin(ipST)),1)
*
          If (lprint) Write(6,*)
     &      '       Iteration         Delta     Res(kappa) Res(CI)'
          iLen=nDensC
          iRHSDisp(iDisp)=iDis
          Call Compress(Work(ipTemp4),Work(ipSigma),1)
          Call DSCAL_(ndensc,-sqrt(1.5d0)*DBLE(ms2p),
     &               Work(ipSigma),1)
          Call UnCompress(Work(ipSigma),Work(ipTemp4),1)
          Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
          If (iMethod.eq.2)
     &    call dcopy_(nConf1,0.0d0,0,Work(ipin(ipCIT)),1)
          irc=ipout(ipcit)
          If (iMethod.eq.2) Then
            ilen=nconf1
            iRHSCIDisp(iDisp)=iDis
            Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
          End If
          Call DSCAL_(nConf1,-1.0d0,Work(ipin(ipST)),1)
          Call DSCAL_(nDensC,-1.0d0,Work(ipSigma),1)
*
          Call DMInvKap_sp(Work(ipPre2), Work(ipSigma),
     &                    Work(ipdKap),1)


          If (nconf1.gt.1) then
          Call DMinvCI(ipST,
     &                 Work(ipin(ipCId)),
     &                 rCHC,1)
          Else
            call dcopy_(nconf1,Work(ipin(ipST)),1,
     &              Work(ipin(ipCid)),1)
          End if

*
          If (iMethod.eq.2.and.nconf1.ne.0) Then
            deltaC=ddot_(nConf1,Work(ipin(ipST)),1,Work(ipin(ipCId)),1)
            irc=ipout(ipcid)
          Else
            deltac=0.0d0
          end if
          deltaK=ddot_(nDensC,Work(ipKap),1,Work(ipSigma),1)
          call dcopy_(nDens,0.0d0,0,Work(ipKap),1)
          delta=deltac+deltaK
*         If (delta.eq.0) Goto 300
          delta0=delta
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
             Read(5,*)  i1,j1
             If (i1.gt.0) Then
             call dcopy_(nDens2,0.0d0,0,Work(ipdkap),1)
             Work(ipdKap+i1-1)
     &       =1.0d0
             Else
              call dcopy_(nCONF1,0.0d0,0,Work(ipin(ipCID)),1)
              Work(ipin(ipCid)-i1-1)=1.0d0
              irc=ipout(ipcid)
             End If
*************************************************************
*

            Call RInt_SP(Work(ipDkap),
     &                   Work(iprmoaa),Work(iprmoaa2),
     &                   Work(ipTemp4),Work(ipSc2))

*
             If (i1.gt.0.and.j1.gt.0) Then
               Write(6,*) 'Kap_sig',Work(ipSc2-1+j1)
             End If
             If (nconf1.gt.1) Then
               irc=opout(-1)
               Call CISigma(1,State_Sym,state_sym,
     &                      ipTemp4,iprmoaa,iprmoaa2,
     &                      ipCI,ipS1,'N')
               irc=opout(-1)
               rGrad=ddot_(nconf1,Work(ipin(ipCI)),1,
     &                   Work(ipin(ips1)),1)
               call daxpy_(nConf1,-rgrad,
     &                    Work(ipin(ipCI)),1,Work(ipin(ipS1)),1)
               call dscal_(nconf1,-rms*sqrt(1.5d0)*2.0d0,
     &                    Work(ipin(ips1)),1)
*
             If (i1.gt.0.and.j1.lt.0) Then
               Write(6,*) 'CI_sig',Work(ipin(ips1)-1-j1)
             End If


               irc=opout(-1)
               If (nconf1.gt.1) Then
               Call CISigma(0,State_Sym,state_sym,
     &                    ipFIMO,k2int,idum,
     &                    ipCId,ipS2,'N')
               irc=opout(-1)
               EC=rin_ene+potnuc-ERASSCF(1)

               Call DaXpY_(nConf1,EC,Work(ipin(ipCId)),1,
     &                  Work(ipin(ipS2)),1)
               Call DSCAL_(nConf1,2.0d0,Work(ipin(ipS2)),1)
             If (i1.lt.0.and.j1.lt.0) Then
               Write(6,*) 'CI_sig',Work(ipin(ips2)-1-j1)
             End If
*

*
               Call SpinDens(Work(ipin(ipCI)),Work(ipin(ipCid)),
     &                   State_Sym,State_sym,
     &                   Work(ipp),rdum,rdum,rdum,rdum,
     &                   Work(ipDe),rdum,1)

              d_0=ddot_(nconf1,Work(ipin(ipCid)),1,Work(ipin(ipci)),1)
               Call FockGen_sp(d_0,Work(ipDe),Work(ipP),
     &                    Work(ipSc3),Work(ipSc1),1)
               Call DSCAL_(ndens2,-rms*sqrt(1.5d0),Work(ipSc1),1)

               Call Compress(Work(ipSc1),Work(ipSc3),1)
               If (i1.lt.0.and.j1.gt.0) Then
                Write(6,*) 'CI_sig',Work(ipsc3-1+j1)
               End If
               goto 200
           end if
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
*
*        Add together
*
*
**********************************************************************
*
*          Call GetMem('AKKK','CHECK','REAL',idum,idum)
           if (nconf1.gt.1) then
            Call DZaXpY(nDens,One,Work(ipSc2),1,
     &            Work(ipSc3),1,Work(ipTemp4),1)
           Else
            call dcopy_(nDens,Work(ipSc2),1,Work(ipTemp4),1)
           End If
           call dcopy_(nDens,Work(ipdKap),1,Work(ipTemp2),1)
           If (nconf1.gt.1) Then
             Call DaXpY_(nConf1,1.0d0,
     &            Work(ipin1(ipS2,nconf1)),1,
     &           Work(ipin1(ipS1,nconf1)),1)
           Else
             call dcopy_(nconf1,0.0d0,0,Work(ipin1(ipS1,nconf1)),1)
           End If

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
           rAlphaK=ddot_(nDensC,Work(ipTemp4),1,Work(ipTemp2),1)
           If (nconf1.ne.0)
     &      rAlphaC=ddot_(nConf1,Work(ipin(ipS1)),1,
     &                   Work(ipin(ipCId)),1)
           rAlpha=delta/(rAlphaK+ralphaC)
*
*-------------------------------------------------------------------*
*
*          Kappa=Kappa+rAlpha*dKappa
*          Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
*
           Call DaxPy_(nDensC,ralpha,Work(ipTemp2),1,Work(ipKap),1)
           Call DaxPy_(nDensC,-ralpha,Work(ipTemp4),1,Work(ipSigma),1)
           resk=sqrt(ddot_(nDensC,Work(ipTemp4),1,Work(ipTemp4),1))
           resci=0.0d0
           If (nconf1.ne.0) Then
             Call DaXpY_(nConf1,ralpha,Work(ipin(ipCId)),1,
     &                  Work(ipin(ipCIT)),1)
             irc=ipout(ipcit)
             Call DaXpY_(nConf1,-ralpha,Work(ipin(ipS1)),1,
     &                   Work(ipin1(ipST,nconf1)),1)
             irc=opout(ipS1)
             ip=ipin(ipst)
             resci=sqrt(ddot_(nconf1,Work(ip),1,
     &                        Work(ip),1))
           End If
*
*          Precondition......
*             -1
*          S=M  Sigma
*
           irc=opout(ipcid)
           If (nconf1.gt.1) Then
           Call DMinvCI(ipST,
     &                  Work(ipin(ipS2)),
     &                  rCHC,1)
           Else
            call dcopy_(nconf1,Work(ipin(ipST)),1,
     &                  Work(ipin(ipS2)),1)
           end if


           irc=opout(ipci)
           irc=opout(ipdia)

           Call DMInvKap_sp(Work(ipPre2),work(ipSigma),
     &                      work(ipSC2),1)
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
           If (iMethod.eq.2.and.nconf1.ne.0) Then
           deltaC=ddot_(nConf1,Work(ipin(ipST)),1,Work(ipin(ipS2)),1)
           irc=ipout(ipST)
           else
           deltaC=0.0d0
           end if
*
           deltaK=ddot_(nDensC,Work(ipSigma),1,Work(ipSc2),1)
           If (imethod.ne.2) Then
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
     &                             Work(ipin(ipCID)),1)
             Call DaXpY_(nDensC,1.0d0,Work(ipsc2),1,Work(ipTemp2),1)
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
           call dcopy_(ndensc,Work(ipTemp2),1,Work(ipdKap),1)
*
           res=0.0D0 ! dummu initialize
           If (iBreak.eq.1) Then
              If (abs(delta).lt.abs(Epsilon**2*delta0)) Goto 300
           Else If (ibreak.eq.2) Then
              res=sqrt(resk**2+resci**2)
              If (res.lt.abs(epsilon)) Goto 300
           Else
              If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &            res.lt.abs(epsilon))  Goto 300
           End If
           If (iter.ge.niter) goto 210
           If (lprint)
     &     Write(6,Fmt2//'A,i2,A,F12.7,F12.7,F12.7,F12.7,F12.7)')
     &            '     ',
     &            iter,'       ',delta/delta0,resk,resci,deltac,deltak

           iter=iter+1
*          Call GetMem(' LIST ','LIST','REAL',iDum,iDum)

          Goto 200
*
**********************************************************************
*
 210      Continue
          Write(6,Fmt2//'A,I4,A)')
     &    'No convergence for perturbation no: ',
     &                      idisp,'. Increase Iter.'
          fail=.true.
          Goto 310
 300      Write(6,Fmt2//'A,I4,A,I4,A)')
     &           'Perturbation no: ',idisp,' converged in ',
     &             iter-1,' steps.'
          irc=ipnout(-1)
 310      Continue
          Write(6,*)
          iLen=ndensC
          iKapDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Work(ipKap),iLen,iDis)
          iSigDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
          If (iMethod.eq.2) Then
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
             Call GetMem('2Dens','FREE','Real',iprmoaa,n2dens)
             Call GetMem('2Dens','FREE','Real',ipP,n1dens)
             Call GetMem('1Dens','FREE','Real',ipDe,n2dens)
          End If
*
*        Free all memory and remove from disk all data
*        related to this symmetry
*
         If (imethod.eq.2) irc=ipclose(ipci)
*
         If (iphx.ne.0) Then
          Call Getmem('EXPHS','FREE','REAL',iphx,idum)
          Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
          Call Getmem('EXPLS','FREE','INTE',iplst,idum)
         End If

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
