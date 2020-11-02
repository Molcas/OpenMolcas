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
      use Arrays, only: SFock, G1m, G2mp, Int2, FIMO
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)
*
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
      Real*8 d_0
      Real*8, Allocatable:: Kappa(:), dKappa(:), Sigma(:),
     &                      Temp1(:), Temp2(:), Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Sc3(:),
     &                      Dens(:), Pens(:), rmoaa(:), rmoaa2(:),
     &                      Pre2(:)
*
*----------------------------------------------------------------------*
*
       Interface
         SubRoutine FockGen_sp(d_0,rDens1,rdens2,Fock,fockout,idsym)
         Real*8 d_0
         Real*8 rDens1(*), rdens2(*), Fock(*), fockout(*)
         Integer idsym
         End SubRoutine FockGen_sp
       End Interface
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
          Call FockGen_sp(Zero,G1m,G2mp,SFock,Temp4,1)
          irc=ipin(ipST)
          call dcopy_(nconf1,[Zero],0,W(ipST)%Vec,1)
*
          If (lprint) Write(6,*)
     &      '       Iteration         Delta     Res(kappa) Res(CI)'
          iLen=nDensC
          iRHSDisp(iDisp)=iDis
          Call Compress(Temp4,Sigma,1)
          Call DSCAL_(ndensc,-sqrt(1.5d0)*DBLE(ms2p),Sigma,1)
          Call UnCompress(Sigma,Temp4,1)
          Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
          If (iMethod.eq.2) Then
             irc=ipin(ipCIT)
             call dcopy_(nConf1,[Zero],0,W(ipCIT)%Vec,1)
          End If
          irc=ipout(ipcit)
          irc=ipin(ipST)
          If (iMethod.eq.2) Then
            ilen=nconf1
            iRHSCIDisp(iDisp)=iDis
            Call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
          End If
          Call DSCAL_(nConf1,-One,W(ipST)%Vec,1)
          Call DSCAL_(nDensC,-One,Sigma,1)
*
          Call DMInvKap_sp(Pre2,Sigma,dKappa,1)

          irc=ipin(ipCId)
          If (nconf1.gt.1) then
             Call DMinvCI(ipST,W(ipCId)%Vec,rCHC,1)
          Else
            irc=ipin(ipST)
            call dcopy_(nconf1,W(ipST)%Vec,1,W(ipCid)%Vec,1)
          End if

*
          If (iMethod.eq.2.and.nconf1.ne.0) Then
            irc=ipin(ipST)
            irc=ipin(ipCId)
            deltaC=ddot_(nConf1,W(ipST)%Vec,1,W(ipCId)%Vec,1)
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
              irc=ipin(ipCID)
              W(ipCID)%Vec(1:nConf1)=Zero
              W(ipCID)%Vec(i1)=One
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
     &                      ipCI,ipS1)
               irc=opout(-1)
               irc=ipin(ipCI)
               irc=ipin(ipS1)
               rGrad=ddot_(nconf1,W(ipCI)%Vec,1,W(ipS1)%Vec,1)
               call daxpy_(nConf1,-rgrad,W(ipCI)%Vec,1,W(ipS1)%Vec,1)
               call dscal_(nconf1,-rms*sqrt(1.5d0)*Two,W(ipS1)%Vec,1)
*
             If (i1.gt.0.and.j1.lt.0) Then
               Write(6,*) 'CI_sig',W(ipS1)%Vec(j1)
             End If


               irc=opout(-1)
               If (nconf1.gt.1) Then
               Call CISigma(0,State_Sym,state_sym,
     &                    FIMO,Int2,rdum,
     &                    ipCId,ipS2)
               irc=opout(-1)
               EC=rin_ene+potnuc-ERASSCF(1)

               irc=ipin(ipCId)
               irc=ipin(ipS2)
               Call DaXpY_(nConf1,EC,W(ipCId)%Vec,1,W(ipS2)%Vec,1)
               Call DSCAL_(nConf1,2.0d0,W(ipS2)%Vec,1)
             If (i1.lt.0.and.j1.lt.0) Then
               Write(6,*) 'CI_sig',W(ipS2)%Vec(j1)
             End If
*

*
               irc=ipin(ipCI)
               irc=ipin(ipCid)
               Call SpinDens(W(ipCI)%Vec,W(ipCid)%Vec,
     &                       State_Sym,State_sym,
     &                       Pens,rdum,rdum,rdum,rdum,
     &                       Dens,rdum,1)

               d_0=ddot_(nconf1,W(ipCid)%Vec,1,W(ipci)%Vec,1)
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
              irc=ipin1(ipS1,nconf1)
              irc=ipin1(ipS2,nconf1)
              Call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipS1)%Vec,1)
           Else
              irc=ipin1(ipS1,nconf1)
              W(ipS1)%Vec(1:nconf1)=Zero
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
           If (nconf1.ne.0) Then
              irc=ipin(ipS1)
              irc=ipin(ipCId)
              rAlphaC=ddot_(nConf1,W(ipS1)%Vec,1,W(ipCId)%Vec,1)
           End If
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
             irc=ipin(ipCId)
             irc=ipin(ipCIT)
             Call DaXpY_(nConf1,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
             irc=ipout(ipcit)
             irc=ipin1(ipST,nconf1)
             irc=ipin(ipS1)
             Call DaXpY_(nConf1,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
             irc=opout(ipS1)
             irc=ipin(ipST)
             resci=sqrt(ddot_(nconf1,W(ipST)%Vec,1,W(ipST)%Vec,1))
           End If
*
*          Precondition......
*             -1
*          S=M  Sigma
*
           irc=opout(ipcid)
           irc=ipin(ipS2)
           If (nconf1.gt.1) Then
              Call DMinvCI(ipST,W(ipS2)%Vec,rCHC,1)
           Else
              irc=ipin(ipST)
              call dcopy_(nconf1,W(ipST)%Vec,1,W(ipS2)%Vec,1)
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
              irc=ipin(ipST)
              irc=ipin(ipS2)
              deltaC=ddot_(nConf1,W(ipST)%Vec,1,W(ipS2)%Vec,1)
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
             irc=ipin(ipCID)
             Call DScal_(nConf1,rBeta,W(ipCID)%Vec,1)
             Call DScal_(nDensC,rBeta,Temp2,1)
             irc=ipin(ipS2)
             Call DaXpY_(nConf1,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
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
            irc=ipin(ipCIT)
            Call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)
            iCISigDisp(iDisp)=iDis
            irc=ipin(ipST)
            Call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
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
