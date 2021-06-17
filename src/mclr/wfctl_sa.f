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
      ! nconf1: CSF
      ! nconf3: determinant
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
      If (PT2) Then
        Call Getmem('SLag  ','ALLO','Real',ipSLag ,nRoots*(nRoots-1)/2)
        Call DCopy_(nRoots*(nRoots-1)/2,[0.0d+00],0,Work(ipSLag),1)
        Call RHS_PT2(Work(ipKap),ipST,Work(ipIn(ipST)),Work(ipSLag))
        ! Call Dcopy_(ndens2,[Zero],0,Work(ipKap),1)
        ! Call Dcopy_(nRoots*(nRoots-1)/2,[Zero],0,Work(ipSLag),1)
      End If
      If (isNAC) Then
        Call RHS_NAC(Work(ipTemp4))
      Else
        Call RHS_SA(Work(ipTemp4),Work(ipSLag))
      End If
C     write(6,*) "casscf lag"
C     call sqprt(work(iptemp4),norb(1))
      If (PT2) Then
C       write(6,*) "pt2 lag"
C       call sqprt(work(ipKap),norb(1))
        Call Getmem('SLag  ','FREE','Real',ipSLag ,nRoots*(nRoots-1)/2)
        Call DaXpY_(nDens2,1.0D+00,Work(ipKap),1,Work(ipTemp4),1)
        Call DCopy_(nDens2,[0.0D+00],0,Work(ipKap),1)
      End If
      ! Call Dcopy_(ndens2,[Zero],0,Work(ipTemp4),1)
C     write(6,*) "casscf+pt2 lag"
C     call sqprt(work(iptemp4),norb(1))
*
      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta           Res(kappa)     '//
     &      '  Res(CI)          DeltaK           DeltaC'
      iLen=nDensC !! number of independent rotations
      iRHSDisp(iDisp)=iDis
      !! orbital rotation: maybe square to indepdent
      Call Compress(Work(ipTemp4),Work(ipSigma),iSym)
C     write(6,*) "compressed orbital Lagrangian"
C     do i = 1, ndensc
C       write(6,'(i3,f20.10)') i,work(ipsigma+i-1)
C     end do
C     write(6,*) "Work(ipIn(ipST))"
C     do i = 1, nconf3
C       write(6,'(i3,f20.10)') i,work(ipin(ipST)+i-1)
C     end do
C     Call RecPrt('RHS',' ',Work(ipSigma),nDensc,1)
      r1=ddot_(nDensc,Work(ipSigma),1,Work(ipSigma),1)
      If (PT2) R1 = R1 + DDot_(nConf1*nRoots,Work(ipIn(ipST)),1,
     *                                       Work(ipIn(ipST)),1)
C      write(6,*) "r1"
C     write(6,*)ddot_(nDensc,Work(ipSigma),1,Work(ipSigma),1)
C     write(*,*) DDot_(nConf1*nRoots,Work(ipIn(ipST)),1,
C    *                                       Work(ipIn(ipST)),1)
      If(debug)Write(6,*) 'Hi how about r1',r1
C     Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
*
      call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCIT)),1)
      If (PT2) then
        Call DSCAL_(nConf1*nRoots,-One,Work(ipIn(ipST)),1)
* Naming System:
* ipKap: accumulates Lagrange multiplier orbital parts (ipdKap * ralpha)
* ipdKap: orbital input of Hessian matrix-vector;
* ipTemp4: orbital output of Hessian matrix-vector
* ipSigma: accumulates error vector orbital part
* ipCIT: accumulates Lagrange multiplier CI parts (ipCId * ralpha)
* ipCId: CI input of Hessian matrix-vector;
* ipS1: CI output of Hessian matrix-vector
* ipST: accumulates error vector CI part
        !! Some post-process are required
        !! the CI Lagrangian has to be projected out
        !! (although not required)
        !! Also, state rotations are to be included
C       call slag(Work(ipIn(ipST)),slag)
        !! initial residue is r = b - Ax, where we have b.
        !! The Ax part has to be computed.
        If (CI) Then
       !  write(6,*) "Work(ipIn(ipST)) in CSF?"
       !  write(6,'(5f20.10)') (work(ipin(ipST)+i-1),i=1,nconf1)
       !  Call DCopy_(nConf3*nRoots,0.0D+00,0,Work(ipIn(ipS2)),1)
       !  Call CSF2SD(Work(ipIn(ipST)),Work(ipIn(ipS2)),iSym)
C     write(6,*) "Work(ipIn(ipS2)) in det."
C     do i = 1, nconf3
C       write(6,'(i3,f20.10)') i,work(ipin(ipS2)+i-1)
C     end do
       !  write(6,*) "Work(ipIn(ipS2)) in determinant?"
       !  write(6,'(5f20.10)') (work(ipin(ipS2)+i-1),i=1,nconf3)
       !  Call DCopy_(nConf3*nRoots,Work(ipIn(ipS2)),1,
     * !                            Work(ipIn(ipST)),1)
       !  !! input : ipST
       !  !! output: ipS2
       !  Call DMinvCI_sa(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
       !  Call DCopy_(nConf3*nRoots,Work(ipIn(ipS2)),1,
     * !                            Work(ipIn(ipST)),1)
C         Call DMinvCI_sa(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
C         Call DCopy_(nConf3*nRoots,Work(ipIn(ipS2)),1,
C    *                              Work(ipIn(ipST)),1)
C         Call CSF2SD(Work(ipIn(ipST)),Work(ipIn(ipS2)),iSym)
C         Call DCopy_(nConf3*nRoots,Work(ipIn(ipS2)),1,
C    *                              Work(ipIn(ipST)),1)
C         Call DCopy_(nConf3*nRoots,Work(ipIn(ipS2)),1,
C    *                              Work(ipIn(ipCID)),1)
          !! The order of CSF coefficients in CASPT2 and MCLR is somehow
          !! different, so the CI lagrangian computed in CASPT2 must be
          !! reordered so that it can be used here.
          Call Getmem('WRK   ','ALLO','Real',ipWRK  ,nConf1)
          Do iR = 1, nRoots
C         write(6,*) "Work(ipIn(ipST)) in CSF, Root =",ir
C         do i = 1, nconf1
C         write(6,'(i3,f20.10)') i,work(ipin(ipST)+i-1+nConf1*(iR-1))
C         end do
            Call DCopy_(nConf1,Work(ipIn(ipST)+nConf1*(iR-1)),1,
     *                  Work(ipWRK),1)
            Call GugaCtl_MCLR(ipWRK,1)
            Call DCopy_(nConf1,Work(ipWRK),1,
     *                  Work(ipIn(ipST)+nConf1*(iR-1)),1)
C         write(6,*) "Work(ipIn(ipST)) in CSF, Root =",ir
C         do i = 1, nconf1
C         write(6,'(i3,f20.10)') i,work(ipin(ipST)+i-1+nConf1*(iR-1))
C         end do
          End Do
          Call Getmem('WRK   ','FREE','Real',ipWRK  ,nConf1)
          if (.false.) then
          Call DaXpY_(nConf1*nRoots,-1.0D+00,
     *                Work(ipIn(ipS1)),1,Work(ipIn(ipST)),1)
          end if
C
C         call dcopy_(nconf1*nroots,0.0d+00,0,work(ipin(ipst)),1)
C         Call DScal_(nConf1*nRoots,-1.0D+00,Work(ipIn(ipST)),1)
          !! precondition (z0 = M^{-1}*r0)
          Call DMinvCI_sa(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
         irc=opOut(ipci)
         irc=opOut(ipdia)
       !!!Call DMinvCI(ipST,Work(ipIn(ipS2)),rdum,isym)
C         Call SD2CSF(ipS2,Work(ipIn(ipST)),iSym)
          !   ovl = ddot_(nconf1,work(ipin(ips2)),1,work(ipin(ipci)),1)
          !   call daxpy_(nconf1,-ovl,work(ipci),1,work(ipin(ips2)),1)
          !! z0 <= p0
          Call DCopy_(nConf1*nRoots,Work(ipIn(ipS2)),1,
     *                              Work(ipIn(ipCId)),1)
C         Call DCopy_(nConf1*nRoots,Work(ipIn(ipST)),1,
C    *                              Work(ipIn(ipCId)),1)
C         Call DCopy_(nConf1*nRoots,Work(ipIn(ipST)),1,
C    *                              Work(ipIn(ipS2 )),1)
        End If
C       write(6,*) "after preconditioning?"
C       write(6,*) "Work(ipIn(ipST))"
C       write(6,'(5f20.10)') (work(ipin(ipST)+i-1),i=1,nconf1)
C       write(6,*) "Work(ipIn(ipCIT))"
C       write(6,'(5f20.10)') (work(ipin(ipCIT)+i-1),i=1,nconf1)
C       write(6,*) "Work(ipIn(ipCId))"
C       write(6,'(5f20.10)') (work(ipin(ipCId)+i-1),i=1,nconf1)
C       write(6,*) "Work(ipIn(ipS2))"
C       write(6,'(5f20.10)') (work(ipin(ipS2)+i-1),i=1,nconf1)
C       Call TimesE2(Work(ipdKap),ipCId,1,reco,jspin,ipS2,
C    &              Work(ipTemp4),ipS1,0)
C       irc=ipout(ipcid)
C       Call DaXpY_(nConf1*nroots,-1.0D+00,
C    *              Work(ipIn(ipS1)),1,
C    *              Work(ipIn(ipST)),1)
C       Call DaXpY_(nConf1*nRoots,1.0D+00,
C    *              Work(ipIn(ipCId)),1,
C    *              Work(ipIn(ipCIT)),1)
C
C       Call DCopy_(nConf1*nRoots,Work(ipIn(ipST)),1,
C    *                            Work(ipIn(ipCID)),1)
      Else
        call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipST)),1)
        call dcopy_(nConf1*nroots,[Zero],0,Work(ipIn(ipCID)),1)
      End If
      irc=ipOut(ipCIT)
C     weight = 0.99331749869259100549d+00
      Call DSCAL_(nDensC,-One,Work(ipSigma),1)
*
      !! set x0 = 0
      !! r0 = b = RHS    // ipST  // ipSigma
      !! z0 = M^{-1}*r0  // ipS2  // ipKap
      !! p0 = z0         // ipCId // ipdKap
      Call DMInvKap(Work(ipIn(ipPre2)),Work(ipSigma),nDens2+6,
     &              Work(ipKap),nDens2+6,work(ipTemp3),nDens2+6,
     &              isym,iter)
C      write(6,*) "M^{-1}*r0"
C      do i = 1, ndensc
C       write(6,'(i3,f20.10)') i,work(ipkap+i-1)
C      end do
C
      IF (PT2) THEN
      END IF

      !! z0^T*z0
      irc=opOut(ippre2)
C     write(6,*) "r2"
C     write(6,*) ddot_(ndensc,Work(ipKap),1,Work(ipKap),1)
C     write(*,*) DDot_(nConf1*nRoots,Work(ipIn(ipS2)),1,
C    *                                       Work(ipIn(ipS2)),1)
      r2=ddot_(ndensc,Work(ipKap),1,Work(ipKap),1)
      If (PT2) R2 = R2 + DDot_(nConf1*nRoots,Work(ipIn(ipS2)),1,
     *                                       Work(ipIn(ipS2)),1)
      If(debug)Write(6,*) 'In that case I think that r2 should be:',r2
C     Write(6,*) 'In that case I think that r2 should be:',r2
      If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
      !! define p0 <- z0
      call dcopy_(ndensC,Work(ipKap),1,Work(ipdKap),1)
C       write(6,*) "CI, r0,z0,p0"
C       do i = 1, nconf1
C         write(6,'(i3,3f20.10)')
C    *    i,Work(ipIn(ipST)+i-1),Work(ipIn(ipS2)+i-1),
C    *    Work(ipIn(ipCId)+i-1)
C       end do
C       write(6,*) "ORB, r0,z0,p0"
C       do i = 1, ndensc
C         write(6,'(i3,3f20.10)')
C    *    i,Work(ipSigma+i-1),Work(ipKap+i-1),
C    *    Work(ipdKap+i-1)
C       end do
*
      ! r^T dot z
      ! r (residue) = ipST (det?)
      ! z (prec. r) = ipS2 (det?)  // ipSc2
      ! p (...)     = ipCId (CSF)  // ipdKap
      ! x (solution)= ipCIT (CSF)  // ipKap
      ! Ap          = ipS1 (det?)  // ipTemp4
      ! r_{k}z_{k}  = ipST*ipS2 = deltaC
      deltaC=Zero
      If (PT2) deltaC=ddot_(nConf1*nroots,Work(ipin(ipST)),1,
     &                                    Work(ipin(ipS2)),1)
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
         ! Compute Ap
         Call TimesE2(Work(ipdKap),ipCId,1,reco,jspin,ipS2,
     &                Work(ipTemp4),ipS1,0)
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
         !! new x  // x_{k+1} = x_{k} + alpha_{k}*p_{k}
         !! new r // r_{k+1} = r_{k} - alpha_{k}*A*p_{k}
*        Kappa=Kappa+rAlpha*dKappa
         !  new x of orbital
         Call DaxPy_(nDensC,ralpha,Work(ipdKap),1,Work(ipKap),1)
*        Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
         !  new r of orbital
         Call DaxPy_(nDensC,-ralpha,Work(ipTemp4),1,Work(ipSigma),1)
         resk=sqrt(ddot_(nDensC,Work(ipSigma),1,Work(ipSigma),1))
C        resk = resk - 7.921473959d-04
         !!
         resci=Zero
         !  new x of CI
         Call DaXpY_(nConf1*nroots,ralpha,Work(ipIn(ipCId)),1,
     &                                   Work(ipIn(ipCIT)),1)
         irc=ipOut(ipcit)
*        ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
         !  new r of CI
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

         !! new z for orbital // z_{k+1} = M^{-1}*r_{k+1}
         !! ipS2 = M^{-1}*ipST
         Call DMinvCI_SA(ipST,Work(ipIn(ipS2)),rdum,isym,work(ipS))
         irc=opOut(ipci)
         irc=opOut(ipdia)

C          write(6,*) "calling dminvkap in iter",iter
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
         !! beta_{k} = r_{k+1}^T*z_{k+1} / r_{k}^T*z_{k}
         !  where r_{k+1} = ipST*ipS2/
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
            !! p_{k+1} = z_{k+1} + beta_{k}*p_{k}
            !  ipCId   = ipS2 + beta*ipCId
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
     &          iter,delta/delta0,resk,resci,deltak,deltac
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
C       write(6,*) "Work(ipKap)"
C       do i = 1, ilen
C         write(6,'(i3,f20.10)') i,work(ipkap+i-1)
C       end do
      iSigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Work(ipSigma),iLen,iDis)
        !! Work(ipSigma) is the orbital residue
C       write(6,*) "Work(ipSigma)"
C       do i = 1, ilen
C         write(6,'(i3,f20.10)') i,work(ipSigma+i-1)
C       end do
      ilen=nconf1*nroots
      iCIDisp(iDisp)=iDis
*
      Call dDaFile(LuTemp,1,Work(ipin(ipCIT)),iLen,iDis)
C       write(6,*) "Work(ipin(ipCIT))"
C       do i = 1, ilen
C         write(6,'(i3,f20.10)') i,work(ipin(ipCIT)+i-1)
C       end do
C       ovl = ddot_(nconf1,work(ipin(ipcit)),1,work(ipin(ipci)),1)
C       write(6,*) "overlap",ovl
      !!------------
C       write(6,*) "Check A*Z-L"
       ! Call TimesE2(Work(ipKap),ipCIT,1,reco,jspin,ipS2,
     & !              Work(ipTemp4),ipS1,0)
         Call TimesE2(Work(ipKap),ipCIT,1,reco,jspin,ipS2,
     &                Work(ipTemp4),ipS1,2)
C       write(6,*) "A*Z for orbital"
C       do i = 1, ndensc
C         write(6,'(i3,f20.10)') i,work(iptemp4+i-1)
C       end do
C       write(6,*) "A*Z for CI"
C       do i = 1, nconf1*nroots
C         write(6,'(i3,f20.10)') i,work(ipIn(ipS1)+i-1)
C       end do
      !!------------
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

      Subroutine TimesE2(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut,
     *                   mode)

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
      If (ActRot) Then
      Call GetMem('SCR5','ALLO','REAL',ipTemp5,ntash*ntash)
      End If
*
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      If (Mode.eq.0) Then
        Call Uncompress(Kap,Work(ipSC1),isym)
      Else If (Mode.eq.1) Then
        Call DCopy_(ndens2,Kap,1,Work(ipSC1),1)
      End If

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
      If (Mode.eq.0) Then
        Call Compress(Work(ipSc1),KapOut,isym)   ! ds
      Else If (Mode.eq.1) Then
        Call DaXpY_(ndens2,1.0D+00,Work(ipSc1),1,KapOut,1)
      End If
*     Call RecPrt('Ex',' ',KapOut,ndensC,1)
*
      Call DaXpY_(nConf1*nroots,One,
     &               Work(ipin(ipS2)),1,
     &               Work(ipin(ipCIOUT)),1)
      irc=opOut(ipCId)
       do iR = 1, nroots
         do jR = 1, nroots
           ovl = ddot_(nconf1,work(ipin(ipciout)+(iR-1)*nconf1),1,
     *                        work(ipin(ipci)+(jR-1)*nconf1),1)
           call daxpy_(nconf1,-ovl,work(ipin(ipci)+(jR-1)*nconf1),1,
     *                             work(ipin(ipciout)+(iR-1)*nconf1),1)
         end do
        end do

*
      Call GetMem('RMOAA','FREE','REAL',iprmoaa,n2dens)
      Call GetMem('SCR2','FREE','REAL',ipSc2,ndens2)
      Call GetMem('SCR1','FREE','REAL',ipSc1,ndens2)
      Call GetMem('SCR3','FREE','REAL',ipSc3,ndens2)
      Call GetMem('SCR4','FREE','REAL',ipTemp4,ndens2)
      Call GetMem('SCR5','FREE','REAL',ipTemp3,ndens2)
      If (ActRot) Then
      Call GetMem('SCR5','FREE','REAL',ipTemp5,ntash*ntash)
      End If

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*
      Return
      End
      Subroutine SD2CSF(SD,CSF,is)
*
*  Transforms a CSF vector to slater determinants
*
      implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "csfbas_mclr.fh"
#include "WrkSpc.fh"
#include "cicisp_mclr.fh"
#include "Input.fh"
#include "spinfo_mclr.fh"
#include "ippage.fh"
*
      Real*8 CSF(*),SD(*)
*

      iiCOPY=0
      iprdia=0
C     nConf=Max(ncsf(is),ndtasm(iS))
      nConf=nint(Max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
      write(6,*) "nconf in sd2csf = ", nconf
      isym=iEor(is-1,State_Sym-1)+1
      i=2
      If (isym.eq.1) i=1
      write(6,*) "diskbased = ", diskbased
      If (diskbased) Then
         CALL CSDTVC_MCLR(SD,CSF,1,wORK(KDTOC),
     &                    iWORK(KICTS(i)),
     &                    IS,iiCOPY,IPRDIA)
      Else
         Call GetMem('CITEMP','ALLO','REAL',ipCTM,nConf)
         Call FZero(Work(ipCTM),nConf)
C        call dcopy_(ncsf(is),SD,1,Work(ipCTM),1)
         call dcopy_(nconf,SD,1,Work(ipCTM),1)
         CALL CSDTVC_MCLR(Work(ipCTM),CSF,2,WORK(KDTOC),
     &                    iWORK(KICTS(i)),
     &                    IS,iiCOPY,IPRDIA)
         Call GetMem('CITEMP','FREE','REAL',ipCTM,nConf)
      End If
*
      Return
      End
