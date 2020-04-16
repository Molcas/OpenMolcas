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
      Integer nacpar,nacpr2
      Logical lPrint,converged(8)
      Real*8 rchc(mxroot)

      Integer ipFT99,iptemp5
*
      Call QEnter('WfCtl_SA')
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
      Do i=1,8
       Converged(i)=.true.
      end do
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
      iphx=0
      Call Getmem('FANCY','ALLO','REAL',ips,nroots**3)
!____________________
!AMS - what to do here?
! What should rCHC be?  is it computed with E(mcscf) or E(pdft)?

!____________________
      Call CIDia_SA(State_Sym,rCHC,Work(ipS))
      Call xflush(6)
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

      Call Prec(Work(ipIn(ipPre2)),isym)
      irc=ipOut(ippre2)

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
      !I think the lagrange multiplers are independent of the
      !displacement, no?
      nDisp = 1
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
!AMS - I Think I can skip all of this RHS stuff - I'll read it in
!below.
!      If (isNAC) Then
!        Call RHS_NAC(Work(ipTemp4))
!      Else
!        Call RHS_SA(Work(ipTemp4))
!        goto 538
!      End If
*
!AMS _____________________________________________________
!Read in the Fock operator for the calculation of the CI part of the RHS
!ipF1 and ipF2.
      Call xflush(6)
      Open(unit=87,file='TmpFock',action='read',iostat=ios)
      if (ios.ne.0) write(6,*) "error opening file TmpFock"
      nTri = 0
      nTri2 = 0
      nOrbAct = 0
      do ksym=1,nsym
        nTri = nTri + nBas(ksym)*(nBas(ksym)+1)/2
        nOrbAct = nOrbAct + nAsh(ksym)
!      nTri = nTri + nAsh(ksym)*(nAsh(ksym)+1)/2
      end do
      !nTri2 = nTri*(nTri+1)/2
      nacpar = nOrbAct*(nOrbAct+1)/2
      nTri2 = nacpar*(nacpar+1)/2
      Call GetMem('FockOt ','Allo','Real',ipFMO1t,nTri)
      Call GetMem('FockO ','Allo','Real',ipFMO1,ndens2)
      Call GetMem('FockT ','Allo','Real',ipFMO2,n2dens)
      nacpar=(nnA+1)*nnA/2
      nacpr2=(nacpar+1)*nacpar/2
      Call GetMem('FockTt','Allo','Real',ipFMO2t,nacpr2)
      do i=1,nTri
      Call xflush(6)
        read(87,*) Work(ipFMO1t-1+i)
      end do

      ioff = 0
      do iS=1,nSym
        jS=iEOR(iS-1,state_sym-1)+1
        js=is
        If (nBas(is)*nBas(jS).ne.0) then
          if(is.eq.js) then
            do i=1,nBas(iS)
              do j=1,i
            Work(ipFMO1-2+ipMat(is,js) +(i-1)*nbas(iS)+j) =
     &                Work(ipFMO1t+ioff)
               if (i.ne.j) then
            Work(ipFMO1-2+ipMat(is,js) +(j-1)*nbas(iS)+i) =
     &                Work(ipFMO1t+ioff)
               end if
      Call xflush(6)
               ioff = ioff+1
              end do
            end do
          else
            do i=1,nBas(iS)
              do j=1,nBas(jS)
            Work(ipFMO1-1+ipMat(is,js) +(i-1)*nbas(iSym)+j) =
     &                Work(ipFMO1t+ioff)
              ioff = ioff+1
              end do
            end do
          end if
        End if
      end do



      do i=1,nacpr2
        read(87,*) Work(ipFMO2t-1+i)
      end do
      Close(87)

      Call xflush(6)
      iprci = ipget(nconf3)
      Call CISigma_sa(0,State_sym,State_sym,ipFMO1,ipFMO2t,0,
     & ipci,ipST,'N')
      Call GetMem('FockTt','Free','Real',ipFMO2t,nacpr2)

      troot = (irlxroot - 1)
      Do i=0,nroots-1
        if (i.eq.troot) then
          Call Dscal_(nconf1,(1/weight(i+1)),
     &            Work(ipin(ipST)+i*nconf1),1)
          rE=ddot_(nconf1,Work(ipin(ipST)+i*nconf1),1,
     &         Work(ipin(ipci)+i*nconf1),1)
          Call Daxpy_(nconf1,-rE,Work(ipin(ipCI)+i*nconf1),1,
     &              Work(ipin(ipST)+i*nconf1),1)

        else 
        call dcopy_(nConf1,[Zero],0,Work(ipIn(ipst)+i*nconf1),1)
        end if
      Enddo

*     all GetMem('FockTt','Free','Real',ipFMO2t,nacpr2)

      Call DSCAL_(nconf1*nroots,-2.0d0,Work(ipin(ipST)),1)

      if (debug) then
      write(6,*) 'RHS CI part:'
      do iS=1,nconf1*nroots
        write(6,*) Work(ipin(ipST)-1+iS)
      end do
      end if

      Call GetMem('FockOt ','Free','Real',ipFMO1t,nTri)
      Call GetMem('FockO ','Free','Real',ipFMO1,ndens2)
      Call GetMem('FockT ','Free','Real',ipFMO2,n2dens)

!Get the fock matrix needed for the determination of the orbital part of
!the RHS.

      Call GetMem('FockT ','Allo','Real',ipFT99,nDens2)
      Call GetMem('Temp5 ','Allo','Real',ipTemp5,nDens2+6)
      call dcopy_(nDens2,[Zero],0,Work(ipFT99),1)
      call dcopy_(nDens2+6,[Zero],0,Work(ipTemp5),1)
      Call get_dArray('Fock_PDFT',Work(ipFT99),nDens2)
      Do iS=1,nSym
         jS=iEOR(iS-1,0)+1
         If (nBas(is)*nBas(jS).ne.0) then
           Call DGeSub(Work(ipFT99-1+ipMat(iS,jS)),nBas(iS),'N',
     &                  Work(ipFT99-1+ipMat(jS,iS)),nBas(jS),'T',
     &                  Work(ipTemp5-1+ipMat(iS,jS)),nBas(iS),
     &                  nBas(iS),nBas(jS))
         end if
      End Do
      Call dcopy_(nDens2+6,Work(ipTemp5),1,Work(ipTemp4),1)
      Call DSCAL_(ndens2+6,-2.0d0,Work(ipTemp4),1)

      Call GetMem('FockT ','Free','Real',ipFT99,nDens2)
      Call GetMem('Temp5 ','Free','Real',ipTemp5,nDens2+6)
      if (debug) then
      write(6,*) 'RHS orb part:'
      do iS=1,nDens2
        write(6,*) Work(ipTemp4-1+iS)
      end do
      end if
!Also, along with this RHS stuff, the Fock_occ array already stored on
!the runfile needs to be replaced - switch triangular storage to square
!storage:
!
      Call GetMem('FockOSq ','Allo','Real',ipFOSq,nDens2)
      Call dcopy_(nDens2,[Zero],0,Work(ipFOSq),1)
      Call Get_Fock_Occ(ipFOtr,Length)
      Call dcopy_(Length,Work(ipFOtr),1,Work(ipFOSq),1)
      Call Put_Fock_Occ(Work(ipFOSq),ndens2)

!TEMP TEST
!      Call GetMem('FockTri2','ALLO','Real',ipFOtmp,Length)
!      Call dcopy_(nDens2,[Zero],0,Work(ipFOSq),1)
!!      Call Get_Darray('fock_tempo',Work(ipFOtmp),Length)
!      Call dcopy_(Length,Work(ipFOtmp),1,Work(ipFOSq),1)
!      Call Put_Darray('fock_tempS',Work(ipFOSq),ndens2)
!      Call GetMem('FockTri2','Free','Real',ipFOtmp,Length)

!END TEMP TEST

      Call GetMem('FockOSq ','Free','Real',ipFOSq,nDens2)
      Call GetMem('Dens','Free','Real',ipFOTr,Length)


!This seems to calculate the RHS, at least for the orbital part.
!Now, my sigma_0 should be given by
!(RHS) - A*Kappa, where Kappa is my initial guess at the solution, x_0.
!So, should I be running a "TimesE2"-like subroutine here, to do the
!A*Kappa multiplication, before we go on to multiply things by the
!preconditioner inverse?
!___________________________________________________________

      Call xFlush(6)

      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta       Res(kappa)  Res(CI)'
     &    //'     DeltaK      DeltaC'
      iLen=nDensC
      iRHSDisp(iDisp)=iDis
      do iS=1,nDens2
      end do
      Call Compress(Work(ipTemp4),Work(ipSigma),iSym)
      r1=ddot_(nDensc,Work(ipSigma),1,Work(ipSigma),1)
      If(debug)Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)

      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCIT)),1)
      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCID)),1)
      irc=ipOut(ipCIT)
      Call DSCAL_(nDensC,-One,Work(ipSigma),1)
*
*
*
*
      deltaC=Zero
!AMS _________________________________________________________
!I need to read in the CI portion of the RHS here.
      If (CI) Call DMinvCI_sa(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
      Call dcopy_(nconf1*nroots,Work(ipin(ipst)),1,
     &   Work(ipin(ipCId)),1)
********************
*TRS
       Call GetMem('lmroots','Allo','Real',lmroots,nroots)
       Call GetMem('lmroots_new','Allo','Real',lmroots_new,nroots)
       Call GetMem('kap_new','Allo','Real',kap_new,ndensc)
       Call GetMem('kap_new_temp','Allo','Real',kap_new_temp,ndens)
*
      call Dcopy(ndens,Zero,
     &     0,work(kap_new_temp),1)
*
      call Dcopy(ndensc,Zero,
     &     0,work(kap_new),1)
*
        Call DgeMV_('T', nconf1, nroots, 1.0d0, work(ipin(ipci)),
     & nconf1,work(ipin(ipcid)+(irlxroot-1)*nconf1),
     & 1,0.0d0,work(lmroots),1)
*SA-SA rotations w/in SA space in eigen state basis
      if(debug) then
      write(*,*) 'lmroots'
       do i=1,nroots
        write(*,*) work(lmroots-1+i)
      end do       
      end if
*SA-SA rotations w/in SA space in CSF basis
        call dgemv_('N', nconf1, nroots, 1.0d0, work(ipin(ipci)), 
     & nconf1, work(lmroots),1,0.0d0,
     & work(ipin(ipcid)+(irlxroot-1)*nconf1),1) 
*SA-SA rotations w/in SA space for new lagrange multipliers
       do i=0, nroots-1
       if (i.eq.(irlxroot-1)) then
       work(lmroots_new+i)=0.0d0
       else
         diff=(ERASSCF(i+1)-ERASSCF(irlxroot))
        if(debug) write(*,*) 'diff',diff
         wscale = (1.0d0/(2.0d0*diff))*(1.0d0/weight(i+1))
        if(debug) write(*,*)'wscale',wscale
        if(debug) write(*,*) 'weight', weight(i+1)        
        work(lmroots_new+i)= wscale*work(lmroots+i)
       end if
       end do
*
      if(debug) then
      write(*,*) 'lmroots_new'
       do i=1,nroots
        write(*,*) work(lmroots_new-1+i)
      end do
      end if
*SA-SA rotations w/in SA space for new lagrange multipliers in csf basis
        call dgemv_('N', nconf1, nroots, 1.0d0, work(ipin(ipci)),
     &  nconf1, work(lmroots_new),1,0.0d0,
     &  work(ipin(ipcid)+(irlxroot-1)*nconf1),1)
*
*First iter of PCG
          Call TimesE2_(work(kap_new),ipCId,1,reco,jspin,ipS2,
     &                work(kap_new_temp),ipS1)
*
        Call DgeMV_('T', nconf1, nroots, 1.0d0, work(ipin(ipci)),
     & nconf1,work(ipin(ipst)+(irlxroot-1)*nconf1),
     & 1,0.0d0,work(lmroots),1)
*
       if(debug) then
       write(*,*) 'lmroots_ipst this should be 1lmroots'
       do i=1,nroots
        write(*,*) work(lmroots-1+i)
       end do
       end if
*
        Call DgeMV_('T', nconf1, nroots, 1.0d0, work(ipin(ipci)),
     & nconf1,work(ipin(ips1)+(irlxroot-1)*nconf1),
     & 1,0.0d0,work(lmroots),1)
*
       if(debug) then
       write(*,*) 'lmroots_ips1 thisshould be -lmroots'
       do i=1,nroots
        write(*,*) work(lmroots-1+i)
       end do
       end if
* Initializing some of the elements of the PCG
* Modifying the response
       Call DaXpY_(nConf1*nroots,-1.0d0,Work(ipIn(ipS1)),1,
     &                                    Work(ipIn(ipST)),1)
*
*Kap part put into  ipsigma
      Call DaxPy_(nDensC,-1.0d0,Work(kap_new_temp),1,Work(ipSigma),1)
      Call DaXpY_(nConf1*nroots,1.0d0,Work(ipIn(ipCId)),1,
     &                                   Work(ipIn(ipCIT)),1)
*
       Call dcopy_(nconf1*nroots,Work(ipin(ipst)),
     &   1,Work(ipin(ipcid)),1)
*
         irc=opOut(ipci)
         irc=opOut(ipdia)
*
         Call DMInvKap(Work(ipIn(ipPre2)),Work(ipSigma),nDens2+6,
     &                 Work(ipdKap),nDens2+6,Work(ipSc1),nDens2+6,
     &                 iSym,iter)
         irc=opOut(ippre2)
         r2=ddot_(ndensc,Work(ipdKap),1,Work(ipdKap),1)
         If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
*
       Call GetMem('kap_new','Free','Real',kap_new,ndensc)
       Call GetMem('kap_new_temp','Free','Real',kap_new_temp,ndens)
       Call GetMem('lmroots_new','Free','Real',lmroots_new,nroots)
*TRS
**********************
        Call DgeMV_('T', nconf1, nroots, 1.0d0, work(ipin(ipci)),
     & nconf1,work(ipin(ipst)+(irlxroot-1)*nconf1),
     & 1,0.0d0,work(lmroots),1)
*
      if(debug) then
       write(*,*) 'lmroots_ipst this should be zero'
       do i=1,nroots
        write(*,*) work(lmroots-1+i)
       end do
       end if
       Call GetMem('lmroots','Free','Real',lmroots,nroots)
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
      deltaK=ddot_(nDensC,Work(ipdKap),1,Work(ipSigma),1)
      delta=deltac+deltaK
*         write(*,*)'deltac and deltak', deltac,deltak
      delta0=delta
      iter=1
      If (delta.eq.Zero) Goto 300
* Naming System:
* ipKap: accumulates Lagrange multiplier orbital parts (ipdKap * ralpha)
* ipdKap: orbital input of Hessian matrix-vector; 
* ipTemp4: orbital output of Hessian matrix-vector
* ipSigma: accumulates error vector orbital part 
* ipCIT: accumulates Lagrange multiplier CI parts (ipCId * ralpha)
* ipCId: CI input of Hessian matrix-vector; 
* ipS1: CI output of Hessian matrix-vector
* ipST: accumulates error vector CI part 


*-----------------------------------------------------------------------------
*
200   Continue
*
         Call TimesE2_(Work(ipdKap),ipCId,1,reco,jspin,ipS2,
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
*
         rAlpha=delta/(rAlphaK+rAlphaC)
*
*-------------------------------------------------------------------*
*
*        Kappa=Kappa+rAlpha*dKappa
         Call DaxPy_(nDensC,ralpha,Work(ipdKap),1,Work(ipKap),1)
*        Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
         Call DaxPy_(nDensC,-ralpha,Work(ipTemp4),1,Work(ipSigma),1)
         resk=sqrt(ddot_(nDensC,Work(ipSigma),1,Work(ipSigma),1))
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

         Call DMinvCI_SA(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))

         irc=opOut(ipci)
         irc=opOut(ipdia)
*
         Call DMInvKap(Work(ipIn(ipPre2)),Work(ipSigma),nDens2+6,
     &                 Work(ipSc2),nDens2+6,Work(ipSc1),nDens2+6,
     &                 iSym,iter)
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
         write(6,*) Work(ipKap-1+i)
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
      end do
*
      Call Getmem('Scr2   ','FREE','Real',ipSc2  ,nDens2)
      Call Getmem('Scr1   ','FREE','Real',ipSc1  ,nDens2)
      Call GetMem('Temp4  ','FREE','Real',ipTemp4,nDens)
      Call GetMem('Temp3  ','FREE','Real',ipTemp3,ndens)
      Call GetMem('sigma  ','FREE','Real',ipSigma,nDens)
      Call GetMem('dkappa ','FREE','Real',ipdKap ,nDens)
      Call GetMem('kappa  ','FREE','Real',ipKap  ,nDens)
      Call Getmem('FANCY',  'FREE','REAL',ips,nroots**3)
*
*     Free all memory and remove from disk all data
*     related to this symmetry
*
      irc=ipclose(ipdia)
      If (.not.CI) irc=ipclose(ipPre2)
*
*     Call GetMem('PREC','FREE','Real',ipPRE,nDens2)
      If (iphx.ne.0) Then
         Call Getmem('EXPHS','FREE','REAL',iphx,idum)
         Call Getmem('EXPHF','FREE','INTE',ipvt,idum)
         Call Getmem('EXPLS','FREE','INTE',iplst,idum)
      End If
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
      Call QExit('WfCtl_SA')
      Return
      End

      Subroutine TimesE2_(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)
*
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

