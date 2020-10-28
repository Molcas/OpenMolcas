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
* Copyright (C) 2017, Andrew M. Sand                                   *
*                                                                      *
************************************************************************
      SubRoutine WfCtl_pdft(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                    iRHSDisp,converged,iPL)
************************************************************************
*                                                                      *
*                                                                      *
*     called from: MCLR                                                *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
************************************************************************
      use Exp, Only: Exp_Close
      use ipPage, only: W
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
      Integer nacpar,nacpr2
      Logical lPrint,converged(8)
      Real*8 rchc(mxroot)
      Real*8 rDum(1)

      External IsFreeUnit
      Real*8, Allocatable:: FOSq(:), FOTr(:)
      Real*8, Allocatable:: Kappa(:), dKappa(:), Sigma(:),
     &                      Temp4(:), Sc1(:), Sc2(:), Fancy(:),
     &                      FMO1t(:), FMO1(:), FMO2t(:),
     &                      FT99(:), Temp5(:)
      Real*8, Allocatable:: lmroots(:), lmroots_new(:), Kap_New(:),
     &                      Kap_New_Temp(:)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call StatusLine(' MCLR:',
     &                ' Computing Lagrangian multipliers for MC-PDFT')
*

      lPaper=132
      lLine =120
      left=(lPaper-lLine)/2
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*

      iDis=0

      fail=.false.
      Converged(:)=.true.
*MGD I think this is nice when printed...
      lprint=.true.
      debug=.false.
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
      Call mma_allocate(Fancy,nRoots**3,Label='Fancy')
!____________________
!AMS - what to do here?
! What should rCHC be?  is it computed with E(mcscf) or E(pdft)?

!____________________
      Call CIDia_SA(State_Sym,rCHC,Fancy)
      irc=ipOut(ipdia)


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
      ipCID=ipGet(nconf1*nroots)

*

      npre2=npre(isym)
      ipPre2=ipGet(npre2)

      irc=ipIn(ipPre2)
      Call Prec(W(ipPre2)%Vec,isym)
      irc=ipOut(ippre2)

*
*     OK START WORKING
*
*     idisp=1
      jspin=0
*
*     Allocate areas for scratch and state variables
*
      Call mma_allocate(Kappa,nDens2+6,Label='Kappa')
      Call mma_allocate(dKappa,nDens2+6,Label='dKappa')
      Call mma_allocate(Sigma,nDens2+6,Label='Sigma')
      Call mma_allocate(Temp4,nDens2+6,Label='Temp4')
      Call mma_allocate(Sc1,nDens2+6,Label='Sc1')
      Call mma_allocate(Sc2,nDens2+6,Label='Sc2')

*
      !I think the lagrange multiplers are independent of the
      !displacement, no?
      nDisp = 1
      do iDisp=1,nDisp
      Kappa(1:nDens2)=Zero
      dKappa(1:nDens2)=Zero
      Sigma(1:nDens2)=Zero
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
!AMS - I Think I can skip all of this RHS stuff - I'll read it in
!below.
!      If (isNAC) Then
!        Call RHS_NAC(Temp4)
!      Else
!        Call RHS_SA(Temp4)
!        goto 538
!      End If
*
!AMS _____________________________________________________
!Read in the Fock operator for the calculation of the CI part of the RHS
!ipF1 and ipF2.
      LUTMP=87
      Call Molcas_Open(LUTMP,'TmpFock')
      nTri = 0
      nTri2 = 0
      nOrbAct = 0
      do ksym=1,nsym
        nTri = nTri + nBas(ksym)*(nBas(ksym)+1)/2
        nOrbAct = nOrbAct + nAsh(ksym)
      end do
      nacpar = nOrbAct*(nOrbAct+1)/2
      nTri2 = nacpar*(nacpar+1)/2
      Call mma_allocate(FMO1t,nTri,Label='FMO1t')
      Call mma_allocate(FMO1,nDens2,Label='FMO1')
      nacpar=(nnA+1)*nnA/2
      nacpr2=(nacpar+1)*nacpar/2
      Call mma_allocate(FMO2t,nacpr2,Label='FMO2t')
      do i=1,nTri
        read(LUTMP,*) FMO1t(i)
      end do

      ioff = 0
      do iS=1,nSym
        jS=iS
        If (nBas(is)*nBas(jS).ne.0) then
          if(iS.eq.jS) then
            do i=1,nBas(iS)
              do j=1,i
               ioff = ioff+1
               ji= ipMat(is,js)-1 +(i-1)*nbas(iS)+j
               FMO1(ji) = FMO1t(ioff)
               if (i.ne.j) then
                  ij= ipMat(is,js)-1 +(j-1)*nbas(iS)+i
                  FMO1(ij) = FMO1t(ioff)
               end if
              end do
            end do
          else
            do i=1,nBas(iS)
              do j=1,nBas(jS)
              ioff = ioff+1
              ji= ipMat(is,js)-1 +(i-1)*nbas(iSym)+j
              FMO1(ji) = FMO1t(ioff)
              end do
            end do
          end if
        End if
      end do



      do i=1,nacpr2
        read(LUTMP,*) FMO2t(i)
      end do
      Close(LUTMP)

      iprci = ipget(nconf3)
      Call CISigma_sa(0,State_sym,State_sym,FMO1,FMO2t,
     &                rdum,ipci,ipST,'N')
      Call mma_deallocate(FMO2t)

      troot = (irlxroot - 1)
      irc=ipin(ipST)
      irc=ipin(ipCI)
      Do i=0,nroots-1
        if (i.eq.troot) then
          Call Dscal_(nconf1,(1/weight(i+1)),W(ipST)%Vec(1+i*nconf1),1)
          rE=ddot_(nconf1,W(ipST)%Vec(1+i*nconf1),1,
     &                    W(ipCI)%Vec(1+i*nconf1),1)
          Call Daxpy_(nconf1,-rE,W(ipCI)%Vec(1+i*nconf1),1,
     &                           W(ipST)%Vec(1+i*nconf1),1)

        else
        call dcopy_(nConf1,[Zero],0,W(ipst)%Vec(1+i*nconf1),1)
        end if
      Enddo

      Call DSCAL_(nconf1*nroots,-2.0d0,W(ipST)%Vec,1)

      if (debug) then
      write(6,*) 'RHS CI part:'
      do iS=1,nconf1*nroots
        write(6,*) W(ipST)%Vec(iS)
      end do
      end if

      Call mma_deallocate(FMO1t)
      Call mma_deallocate(FMO1)

!Get the fock matrix needed for the determination of the orbital part of
!the RHS.

      Call mma_allocate(FT99,nDens2,Label='FT99')
      Call mma_allocate(Temp5,nDens2+6,Label='Temp5')
      FT99(:)=Zero
      Temp5(:)=Zero
      Call get_dArray('Fock_PDFT',FT99,nDens2)
      Do iS=1,nSym
         jS=iEOR(iS-1,0)+1
         If (nBas(is)*nBas(jS).ne.0) then
           Call DGeSub(FT99(ipMat(iS,jS)),nBas(iS),'N',
     &                 FT99(ipMat(jS,iS)),nBas(jS),'T',
     &                 Temp5(ipMat(iS,jS)),nBas(iS),
     &                 nBas(iS),nBas(jS))
         End If
      End Do
      Call dcopy_(nDens2+6,Temp5,1,Temp4,1)
      Call DSCAL_(ndens2+6,-Two,Temp4,1)

      Call mma_deallocate(FT99)
      Call mma_deallocate(Temp5)
      if (debug) then
      write(6,*) 'RHS orb part:'
      do iS=1,nDens2
        write(6,*) Temp4(iS)
      end do
      end if
!Also, along with this RHS stuff, the Fock_occ array already stored on
!the runfile needs to be replaced - switch triangular storage to square
!storage:
!
      Call mma_allocate(FOSq,nDens2,Label='FOSq')
      Call mma_allocate(FOTr,nTri  ,Label='FOTr')
      FOSq(:)=Zero
      Call Get_Fock_Occ(FOTr,nTri)
      Call dcopy_(nTri,FOtr,1,FOSq,1)
      Call Put_Fock_Occ(FOSq,ndens2)

      Call mma_deallocate(FOSq)
      Call mma_deallocate(FOTr)


!This seems to calculate the RHS, at least for the orbital part.
!Now, my sigma_0 should be given by
!(RHS) - A*Kappa, where Kappa is my initial guess at the solution, x_0.
!So, should I be running a "TimesE2"-like subroutine here, to do the
!A*Kappa multiplication, before we go on to multiply things by the
!preconditioner inverse?
!___________________________________________________________

      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta       Res(kappa)  Res(CI)'
     &    //'     DeltaK      DeltaC'
      iLen=nDensC
      iRHSDisp(iDisp)=iDis
      do iS=1,nDens2
      end do
      Call Compress(Temp4,Sigma,iSym)
      r1=ddot_(nDensc,Sigma,1,Sigma,1)
      If(debug)Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Sigma,iLen,iDis)

      irc=ipIn(ipCIT)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
      irc=ipIn(ipCID)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
      irc=ipOut(ipCIT)
      Call DSCAL_(nDensC,-One,Sigma,1)
*
*
*
*
      deltaC=Zero
!AMS _________________________________________________________
!I need to read in the CI portion of the RHS here.
      If (CI) Then
         irc=ipIn(ipS2)
         Call DMinvCI_sa(ipST,W(ipS2)%Vec,rdum(1),isym,Fancy)
      End If
      irc=ipin(ipST)
      irc=ipin(ipCId)
      Call dcopy_(nconf1*nroots,W(ipST)%Vec,1,W(ipCId)%Vec,1)
********************
*TRS
       Call mma_allocate(lmroots,nroots,Label='lmroots')
       Call mma_allocate(lmroots_new,nroots,Label='lmroots_new')
       Call mma_allocate(kap_new,ndensc,Label='kap_new')
       Call mma_allocate(kap_new_temp,ndens,Label='kap_new_temp')
*
      Kap_New(:)=Zero
      Kap_New_Temp(:)=Zero
*
      irc=ipin(ipCI)
      Call DgeMV_('T', nconf1, nroots, One, W(ipCI)%Vec,
     &            nconf1,W(ipCId)%Vec(1+(irlxroot-1)*nconf1),
     &            1,Zero,lmroots,1)
*SA-SA rotations w/in SA space in eigen state basis
      if(debug) Call recprt('lmroots',' ',lmroots,1,nroots)
*SA-SA rotations w/in SA space in CSF basis
        call dgemv_('N', nconf1, nroots, One, W(ipCI)%Vec,
     &              nconf1, lmroots,1,Zero,
     &              W(ipCId)%Vec(1+(irlxroot-1)*nconf1),1)
*SA-SA rotations w/in SA space for new lagrange multipliers
       do i=1, nroots
          if (i.eq.irlxroot) then
             lmroots_new(i)=Zero
          else
             diff=(ERASSCF(i)-ERASSCF(irlxroot))
             if(debug) write(6,*) 'diff',diff
             wscale = (One/(Two*diff))*(One/weight(i))
             if(debug) write(6,*)'wscale',wscale
             if(debug) write(6,*) 'weight', weight(i)
             lmroots_new(i)= wscale*lmroots(i)
          end if
       end do
*
      if(debug) Call recprt('lmroots_new',' ',lmroots_new,1,nroots)
*SA-SA rotations w/in SA space for new lagrange multipliers in csf basis
        call dgemv_('N', nconf1, nroots, One, W(ipCI)%Vec,
     &              nconf1, lmroots_new,1,Zero,
     &              W(ipcid)%Vec(1+(irlxroot-1)*nconf1),1)
*
*First iter of PCG
          Call TimesE2_(kap_new,ipCId,1,reco,jspin,ipS2,
     &                  kap_new_temp,ipS1)
*
        Call DgeMV_('T', nconf1, nroots, One, W(ipCI)%Vec,
     &              nconf1,W(ipST)%Vec(1+(irlxroot-1)*nconf1),
     &              1,Zero,lmroots,1)
*
       if (debug) then
          write(6,*) 'lmroots_ipst this should be 1lmroots'
          Call recprt('lmroots',' ',lmroots,1,nroots)
       end if
*
        irc=ipin(ipS1)
        Call DgeMV_('T', nconf1, nroots, One, W(ipCI)%Vec,
     &              nconf1,W(ipS1)%Vec(1+(irlxroot-1)*nconf1),
     &              1,Zero,lmroots,1)
*
       if (debug) then
          write(6,*) 'lmroots_ips1 thisshould be -lmroots'
          Call recprt('lmroots',' ',lmroots,1,nroots)
       end if
* Initializing some of the elements of the PCG
* Modifying the response
       Call DaXpY_(nConf1*nroots,-One,Work(ipIn(ipS1)),1,
     &                                Work(ipIn(ipST)),1)
*
*Kap part put into  sigma
      Call DaxPy_(nDensC,-One,kap_new_temp,1,Sigma,1)
      Call DaXpY_(nConf1*nroots,1.0d0,Work(ipIn(ipCId)),1,
     &                                   Work(ipIn(ipCIT)),1)
*
       Call dcopy_(nconf1*nroots,Work(ipin(ipst)),1,
     &                          Work(ipin(ipcid)),1)
*
         irc=opOut(ipci)
         irc=opOut(ipdia)
*
         Call DMInvKap(Work(ipIn(ipPre2)),Sigma,nDens2+6,
     &                 dKappa,nDens2+6,Sc1,nDens2+6,iSym,iter)
         irc=opOut(ippre2)
         r2=ddot_(ndensc,dKappa,1,dKappa,1)
         If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
*
       Call mma_deallocate(kap_new)
       Call mma_deallocate(kap_new_temp)
       Call mma_deallocate(lmroots_new)
*TRS
**********************
        Call DgeMV_('T', nconf1, nroots, One, work(ipin(ipci)),
     &              nconf1,work(ipin(ipst)+(irlxroot-1)*nconf1),
     &              1,Zero,lmroots,1)
*
      if(debug) then
         write(6,*) 'lmroots_ipst this should be zero'
         Call recprt('lmroots',' ',lmroots,1,nroots)
       end if
       Call mma_deallocate(lmroots)
*
*
      If (CI) Then
        deltaC=ddot_(nConf1*nroots,Work(ipin(ipST)),
     &   1,Work(ipin(ipCId)),1)
        irc=ipout(ipcid)
      Else
        deltaC=0.0d0
      End If
!AMS_______________________________________________

      irc=ipOut(ipcid)
      deltaK=ddot_(nDensC,dKappa,1,Sigma,1)
      delta=deltac+deltaK
*         write(6,*)'deltac and deltak', deltac,deltak
      delta0=delta
      iter=1
      If (delta.eq.Zero) Goto 300
* Naming System:
* Kappa: accumulates Lagrange multiplier orbital parts (dKappa * ralpha)
* dKappa: orbital input of Hessian matrix-vector;
* Temp4: orbital output of Hessian matrix-vector
* Sigma: accumulates error vector orbital part
* ipCIT: accumulates Lagrange multiplier CI parts (ipCId * ralpha)
* ipCId: CI input of Hessian matrix-vector;
* ipS1: CI output of Hessian matrix-vector
* ipST: accumulates error vector CI part


*-----------------------------------------------------------------------------
*
200   Continue
*
         Call TimesE2_(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)
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
         rAlphaK=ddot_(nDensC,Temp4,1,dKappa,1)
         rAlphaC=Zero
         rAlphaC=ddot_(nConf1*nroots,Work(ipIn(ipS1)),1,
     &                              Work(ipIn(ipCId)),1)
*
         rAlpha=delta/(rAlphaK+rAlphaC)
*
*-------------------------------------------------------------------*
*
*        Kappa=Kappa+rAlpha*dKappa
         Call DaxPy_(nDensC,ralpha,dKappa,1,Kappa,1)
*        Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
         Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
         resk=sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
*
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
*
*-------------------------------------------------------------------*
*
*        Precondition......
*           -1
*        S=M  Sigma
*
         irc=opOut(ipcid)

         Call DMinvCI_SA(ipST,Work(ipIn(ipS2)),rdum(1),isym,Fancy)

         irc=opOut(ipci)
         irc=opOut(ipdia)
*
         Call DMInvKap(Work(ipIn(ipPre2)),Sigma,nDens2+6,
     &                 Sc2,nDens2+6,Sc1,nDens2+6,iSym,iter)
         irc=opOut(ippre2)
*
*
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
*
*
         irc=ipOut(ipST)
*
         deltaK=ddot_(nDensC,Sigma,1,Sc2,1)
         If (.not.CI) Then
            rBeta=deltaK/delta
            delta=deltaK
            Call DScal_(nDensC,rBeta,dKappa,1)
            Call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
         Else
            rbeta=(deltac+deltaK)/delta
            delta=deltac+deltaK

            Call DScal_(nConf1*nroots,rBeta,Work(ipIn(ipCID)),1)
            Call DScal_(nDensC,rBeta,dKappa,1)
            Call DaXpY_(nConf1*nroots,One,Work(ipIn(ipS2)),1,
     &                                     Work(ipIn(ipCID)),1)
            Call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
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
     &   Write(6,Fmt2//'I7,7X,F12.7,F12.7,F12.7,F12.7,F12.7)')
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
        Write(6,Fmt2//'I7,7X,F12.7,F12.7,F12.7,F12.7,F12.7)')
     &          iter,delta/delta0,resk,resci,deltac,deltak
          Write(6,Fmt2//'A,I4,A,I4,A)')
     &          'Perturbation no: ',idisp,' converged in ',
     &          iter-1,' steps.'
      End If
      irc=ipnout(-1)
*
 310  Continue
      If (iPL.ge.2) Write(6,*)
      if (debug) then
       write(6,*) 'outputs'
       write(6,*) 'kappa'
       do i=1,ndens2
         write(6,*) Kappa(i)
       end do
*         Call dcopy_(nconf1*nroots,0.0d0,0,Work(ipin(ipCIT)),1)
       write(6,*) 'cit'
       do i=1,nconf1*nroots
         write(6,*) Work(ipin(ipCIT)-1+i)
       end do
      end if
*
      iLen=ndensC
*
      iKapDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Kappa,iLen,iDis)
      iSigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
      ilen=nconf1*nroots
      iCIDisp(iDisp)=iDis
*
      Call dDaFile(LuTemp,1,Work(ipin(ipCIT)),iLen,iDis)
*
**MGD This last call seems unused, so I comment it
*
*      Call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,
*     &             Temp4,ipS2)
      iCISigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Work(ipin(ipST)),iLen,iDis)
      end do
*
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Sigma)
      Call mma_deallocate(dKappa)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(Fancy)
*
*     Free all memory and remove from disk all data
*     related to this symmetry
*
      irc=ipclose(ipdia)
      If (.not.CI) irc=ipclose(ipPre2)
*
      Call Exp_Close()
*
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
      Return
      End

      Subroutine TimesE2_(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)
*
      Implicit Real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "Pointers.fh"
#include "dmrginfo_mclr.fh"
#include "real.fh"
#include "Input.fh"
      Integer opOut
      Real*8 Kap(*),KapOut(*)
      Real*8 rdum(1)
      Real*8, Allocatable:: RMOAA(:), Sc1(:), Sc2(:), Sc3(:),
     &                      Temp4(:), Temp3(:)
*
      Call mma_allocate(RMOAA,n2dens,Label='RMOAA')
      Call mma_allocate(Sc1,ndens2,Label='Sc1')
      Call mma_allocate(Sc2,ndens2,Label='Sc2')
      Call mma_allocate(Sc3,ndens2,Label='Sc3')
      Call mma_allocate(Temp3,ndens2,Label='Temp3')
      Call mma_allocate(Temp4,ndens2,Label='Temp4')
*
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      Call Uncompress(Kap,Sc1,isym)
      Call RInt_generic(SC1,rmoaa,rdum,Sc2,Temp3,Temp4,Sc3,
     &                 isym,reco,jspin)

      Call Kap_CI(Temp4,rmoaa,ipCIOUT)
      Call Ci_Ci(ipcid,ipS2)
      Call CI_KAP(ipCid,Sc1,Sc3,isym)

      Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)
*
      Call Compress(Sc1,KapOut,isym)   ! ds
*     Call RecPrt('Ex',' ',KapOut,ndensC,1)
*
      Call DaXpY_(nConf1*nroots,One,
     &               Work(ipin(ipS2)),1,
     &               Work(ipin(ipCIOUT)),1)
      irc=opOut(ipCId)

*
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sc3)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(RMOAA)

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*
      Return
      End

