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
      use Exp, only: Exp_Close
      Implicit Real*8 (a-h,o-z)
*
#include "WrkSpc.fh"
#include "stdalloc.fh"

#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "spin.fh"
#include "real.fh"
#include "crun_mclr.fh"
      Character*8   Fmt2
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp),iRHSCIDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      integer opout
      Logical lPrint
      Real*8 rdum(1)
      Real*8, Allocatable:: Kappa(:), dKappa(:), Sigma(:),
     &                      Temp1(:), Temp2(:), Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Sc3(:),
     &                      Dens(:), Pens(:), rmoaa(:), rmoaa2(:),
     &                      SFock(:), Pre2(:)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
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
           If (nconf1.gt.1) Call CIDia(State_Sym,rCHC)
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
        idisp=1
*
*    Allocate areas for scratch and state variables
*
          Call mma_allocate(Kappa,nDens2+6,Label='Kappa')
          Call mma_allocate(SFock,nDens2+6,Label='SFock')
          ipFS = ip_of_Work(SFock)
          Call mma_allocate(dKappa,nDens2+6,Label='dKappa')
          Call mma_allocate(Sigma,nDens2+6,Label='Sigma')
          Call mma_allocate(Temp1,nDens2+6,Label='Temp1')
          Call mma_allocate(Temp2,nDens2+6,Label='Temp2')
          Call mma_allocate(Temp3,nDens2+6,Label='Temp3')
          Call mma_allocate(Temp4,nDens2+6,Label='Temp4')
          Call mma_allocate(Sc1,nDens2+6,Label='Sc1')
          Call mma_allocate(Sc2,nDens2+6,Label='Sc2')
          Call mma_allocate(Sc3,nDens2+6,Label='Sc3')
          Call mma_allocate(Pre2,nDensC,Label='Pre2')
          Temp1(1:nDens2)=Zero
          Kappa(1:nDens2)=Zero
          dKappa(1:nDens2)=Zero
          Sigma(1:nDens2)=Zero
          If (iMethod.eq.2) Then
             Call mma_allocate(Dens,n1dens,Label='Dens')
             Call mma_allocate(Pens,nna**4,Label='Pens')
             Call mma_allocate(rmoaa,nna**4,Label='rmoaa')
             Call mma_allocate(rmoaa2,nna**4,Label='rmoaa2')
          End If
*
*-----------------------------------------------------------------------
*
*    Calculate RHS
*
*-----------------------------------------------------------------------
*
          Call Pre_SP(Pre2,1)
          Call FockGen_sp(Zero,Work(ipG1m),Work(ipG2mp),
     &                   SFock,Temp4,1)
          call dcopy_(nconf1,[Zero],0,Work(ipin(ipST)),1)
*
          If (lprint) Write(6,*)
     &      '       Iteration         Delta     Res(kappa) Res(CI)'
          iLen=nDensC
          iRHSDisp(iDisp)=iDis
          Call Compress(Temp4,Sigma,1)
          Call DSCAL_(ndensc,-sqrt(1.5d0)*DBLE(ms2p),Sigma,1)
          Call UnCompress(Sigma,Temp4,1)
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (iMethod.eq.2)
     &    call dcopy_(nConf1,[Zero],0,Work(ipin(ipCIT)),1)
          irc=ipout(ipcit)
          If (iMethod.eq.2) Then
            ilen=nconf1
            iRHSCIDisp(iDisp)=iDis
            Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
          End If
          Call DSCAL_(nConf1,-One,Work(ipin(ipST)),1)
          Call DSCAL_(nDensC,-One,Sigma,1)
*
          Call DMInvKap_sp(Pre2,Sigma,dKappa,1)

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
            deltac=Zero
          end if
          deltaK=ddot_(nDensC,Kappa,1,Sigma,1)
          Kappa(1:nDens)=Zero
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
                dKappa(1:nDens2)=Zero
                dKappa(i1) =One
             Else
              call dcopy_(nCONF1,[Zero],0,Work(ipin(ipCID)),1)
              Work(ipin(ipCid)-i1-1)=One
              irc=ipout(ipcid)
             End If
*************************************************************
*

            Call RInt_SP(dKappa,rmoaa,rmoaa2,Temp4,Sc2)

*
             If (i1.gt.0.and.j1.gt.0) Then
               Write(6,*) 'Kap_sig',Sc2(j1)
             End If
             If (nconf1.gt.1) Then
               irc=opout(-1)
               Call CISigma(1,State_Sym,state_sym,
     &                      Temp4,rmoaa,rmoaa2,
     &                      ipCI,ipS1,'N')
               irc=opout(-1)
               rGrad=ddot_(nconf1,Work(ipin(ipCI)),1,
     &                   Work(ipin(ips1)),1)
               call daxpy_(nConf1,-rgrad,
     &                    Work(ipin(ipCI)),1,Work(ipin(ipS1)),1)
               call dscal_(nconf1,-rms*sqrt(1.5d0)*Two,
     &                    Work(ipin(ips1)),1)
*
             If (i1.gt.0.and.j1.lt.0) Then
               Write(6,*) 'CI_sig',Work(ipin(ips1)-1-j1)
             End If


               irc=opout(-1)
               If (nconf1.gt.1) Then
               Call CISigma(0,State_Sym,state_sym,
     &                    Work(ipFIMO),Work(k2int),rdum,
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
     &                       State_Sym,State_sym,
     &                       Pens,rdum,rdum,rdum,rdum,
     &                       Dens,rdum,1)

              d_0=ddot_(nconf1,Work(ipin(ipCid)),1,Work(ipin(ipci)),1)
               Call FockGen_sp(d_0,Dens,Pens,Sc3,Sc1,1)
               Call DSCAL_(ndens2,-rms*sqrt(1.5d0),Sc1,1)

               Call Compress(Sc1,Sc3,1)
               If (i1.lt.0.and.j1.gt.0) Then
                Write(6,*) 'CI_sig',Sc3(j1)
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
           if (nconf1.gt.1) then
            Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Temp4,1)
           Else
            call dcopy_(nDens,Sc2,1,Temp4,1)
           End If
           call dcopy_(nDens,dKappa,1,Temp2,1)
           If (nconf1.gt.1) Then
             Call DaXpY_(nConf1,One,
     &            Work(ipin1(ipS2,nconf1)),1,
     &           Work(ipin1(ipS1,nconf1)),1)
           Else
             call dcopy_(nconf1,[Zero],0,Work(ipin1(ipS1,nconf1)),1)
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
           rAlphaC=Zero
           rAlphaK=Zero
           rAlphaK=ddot_(nDensC,Temp4,1,Temp2,1)
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
           Call DaxPy_(nDensC,ralpha,Temp2,1,Kappa,1)
           Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
           resk=sqrt(ddot_(nDensC,Temp4,1,Temp4,1))
           resci=Zero
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

           Call DMInvKap_sp(Pre2,Sigma,Sc2,1)
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
           deltaC=Zero
           end if
*
           deltaK=ddot_(nDensC,Sigma,1,Sc2,1)
           If (imethod.ne.2) Then
             rBeta=deltaK/delta
             delta=deltaK
             Call DScal_(nDensC,rBeta,Temp2,1)
             Call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
           Else
             rbeta=(deltac+deltaK)/delta
             delta=deltac+deltaK
             Call DScal_(nConf1,rBeta,Work(ipin(ipCID)),1)
             Call DScal_(nDensC,rBeta,Temp2,1)
             Call DaXpY_(nConf1,One,Work(ipin(ipS2)),1,
     &                             Work(ipin(ipCID)),1)
             Call DaXpY_(nDensC,One,Sc2,1,Temp2,1)
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
           call dcopy_(ndensc,Temp2,1,dKappa,1)
*
           res=Zero ! dummu initialize
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
          Call dDaFile(LuTemp,1,Kappa,iLen,iDis)
          iSigDisp(iDisp)=iDis
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (iMethod.eq.2) Then
            ilen=nconf1
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
          Call mma_deallocate(dKappa)
          Call mma_deallocate(Sigma)
          Call mma_deallocate(Kappa)
          Call mma_deallocate(Sc3)
          Call mma_deallocate(Sc2)
          Call mma_deallocate(Sc1)
          If (iMethod.eq.2) Then
             Call mma_deallocate(rmoaa2)
             Call mma_deallocate(rmoaa)
             Call mma_deallocate(Pens)
             Call mma_deallocate(Dens)
          End If
*
*        Free all memory and remove from disk all data
*        related to this symmetry
*
         If (imethod.eq.2) irc=ipclose(ipci)
*
         Call Exp_Close()

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
