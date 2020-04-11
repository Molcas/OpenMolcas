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
* Copyright (C) 2008, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE FNO_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                          vfrac,IFQCAN,DoMP2,EMP2,CMO,NCMO)
*****************************************************************************
*                                                                           *
*     Purpose:  setup of Frozen Natural Orbitals CASPT2 (FNO-CASPT2).       *
*                                                                           *
*     Author:   F. Aquilante  (Geneva, May  2008)                           *
*                                                                           *
*****************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nAsh(nSym),nSsh(nSym),
     &        nDel(nSym)
      Integer irc,IFQCAN
      Real*8  EMP2, vfrac
      Logical DoMP2
      Real*8  DeMP2
      Logical MP2_small
      Common / ChFNOPT/ DeMP2, MP2_small
*
      Integer ns_V(8), nAct(8)
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Real*8  TrDP(8), TrDF(8)
*
      Real*8 CMO(*)
*
      irc=0
      MP2_small=.false.
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
*----------------------------------------------------------------------*
*
      nBasT=0
      ntri=0
      nSQ=0
      nBmx=0
      mAsh=0
      nOrb=0
      nVV=0
      Do i=1,nSym
        nAct(i)=0
        ns_V(i)=0
        nBasT=nBasT+nBas(i)
        nOrb=nOrb+nFro(i)+nIsh(i)+nAsh(i)+nSsh(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nVV=nVV+nSsh(i)**2
        nBmx=Max(nBmx,nBas(i))
        mAsh=Max(mAsh,nAsh(i))
      End Do
      IF(nBasT.GT.mxBas) then
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
*
*----------------------------------------------------------------------*
*     Read the molecular orbitals from JobIph                          *
*----------------------------------------------------------------------*
      IF (IFQCAN.EQ.0) Then
         Write(6,'(/6X,A)')
     &   'Need pseudocanonical RASSCF orbitals. See OUTOrbitals keyword'
     & //' in RASSCF input section.'
         Call Abend
      EndIf
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,2*NCMO)
      iCMO=LCMO+NCMO
* This is not the best solution, but I wanted to avoid having to rewrite
* the indexing code below just to use the CMO array directly
      call dcopy_(NCMO,CMO,1,WORK(LCMO),1)
      Call GetMem('Eorb','Allo','Real',ipOrbE,4*nOrb)
      Call Get_darray('RASSCF OrbE',Work(ipOrbE),nOrb)
      iAoff=0
      Do iSym=1,nSym
         ipOrbE_=ipOrbE+iAoff+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (Work(ipOrbE_+k).lt.0.0d0) nAct(iSym)=nAct(iSym)+1
         End Do
         iAoff=iAoff+nBas(iSym)
      End Do
*
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnOrb(iSym)=nBas(iSym)
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)+nAct(iSym)
         nOA=nOA+lnOcc(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnDel(iSym)=nDel(iSym)
      End Do
*
      ip_ZZ=ipOrbE
      ipEorb=ipOrbE+nOrb
      ip_Z=ipEorb
      kEOcc=ipEorb+nOrb
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=ipOrbE+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(lnOcc(iSym),Work(ifr),1,Work(ito),1)
         ifr=ipOrbE+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),Work(ifr),1,Work(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do
      Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
      Call FZero(Work(ip_X),nVV+nOA)
      ip_Y = ip_X + nVV
*
      Call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                      ip_X,ip_Y)
      Call FZero(Work(iCMO),NCMO)
      iOff=0
      Do iSym=1,nSym
         kfr=LCMO+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),Work(kfr),1,Work(kto),1)
         kfr=LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),Work(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
*
*     Diagonalize the pseudodensity to get natural virtual orbitals
*     -------------------------------------------------------------
      iOff=0
      jOff=0
      Do iSym=1,nSym
         if (nSsh(iSym).gt.0) then
           jD=ip_X+iOff
*     Eigenvectors will be in increasing order of eigenvalues
           Call Eigen_Molcas(nSsh(iSym),Work(jD),Work(ip_Z),Work(ip_ZZ))
*     Reorder to get relevant eigenpairs first
           Do j=1,nSsh(iSym)/2
              Do i=1,nSsh(iSym)
                 lij=jD-1+nSsh(iSym)*(j-1)+i
                 kij=jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
                 tmp=Work(lij)
                 Work(lij)=Work(kij)
                 Work(kij)=tmp
              End Do
              tmp=Work(ip_Z-1+j)
              Work(ip_Z-1+j)=Work(ip_Z+nSsh(iSym)-j)
              Work(ip_Z+nSsh(iSym)-j)=tmp
           End Do
*
*     Compute new MO coeff. : X=C*U
           kfr=iCMO+jOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
           kto=LCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
           Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                        1.0d0,Work(kfr),nBas(iSym),
     &                              Work(jD),nSsh(iSym),
     &                        0.0d0,Work(kto),nBas(iSym))
           iOff=iOff+nSsh(iSym)**2
           TrDF(iSym)=ddot_(nSsh(iSym),Work(ip_Z),1,[1.0d0],0)
           ns_V(iSym)=int(vfrac*dble(nSsh(iSym)))
           TrDP(iSym)=ddot_(ns_V(iSym),Work(ip_Z),1,[1.0d0],0)
         endif
         jOff=jOff+nBas(iSym)**2
      End Do
      write(6,*)'------------------------------------------------------'
      write(6,*)'   Symm.     Trace     (Full Dmat)     (Partial Dmat) '
      write(6,*)'------------------------------------------------------'
      STrDF=0.0d0
      STrDP=0.0d0
      Do iSym=1,nSym
        write(6,'(4X,I4,15X,G13.6,4X,G13.6)') iSym,TrDF(iSym),TrDP(iSym)
        STrDF=STrDF+TrDF(iSym)
        STrDP=STrDP+TrDP(iSym)
      End Do
      write(6,*)'------------------------------------------------------'
      write(6,'(A,G13.6,4X,G13.6)')'   Sum :               ',STrDF,STrDP
      write(6,*)'------------------------------------------------------'
*
*     Update the nSsh, nDel for FNO-CASPT2
      Do iSym=1,nSym
         nDel(iSym)=nDel(iSym)+nSsh(iSym)-ns_V(iSym)
         nSsh(iSym)=ns_V(iSym)
      End Do
*
*     Write the resorted MOs back to JobIph
*
      IFQCAN=0 ! MOs need to be recanonicalized on exit
      call dcopy_(NCMO,WORK(LCMO),1,CMO,1)
*
*     Reset MP2_small for this new call to ChoMP2_Drv
      Call Check_Amp(nSym,lnOcc,nSsh,iSkip)
      MP2_small = DoMP2 .and. iSkip.gt.0
      If (MP2_small) Then
*
         Call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh,
     &                         ip_X,ip_Y)
*
         Call GetMem('iD_orb','Allo','Inte',ip_iD,nOrb)
         Do k=1,nOrb
            iWork(ip_iD-1+k) = k
         End Do
         lOff=0
         kOff=0
         jOff=0
         iOff=0
         Do iSym=1,nSym  ! canonical orb. in the reduced virtual space
            jD=ip_X+iOff
            Call Get_Can_Lorb(Work(kEVir+lOff),Work(ipOrbE+jOff),
     &                        nSsh(iSym),lnVir(iSym),
     &                        iWork(ip_iD),Work(jD),iSym)

            kfr=LCMO+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
            kto=iCMO+kOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
            nBx=Max(1,nBas(iSym))
            nSx=Max(1,nSsh(iSym))
            Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                         1.0d0,Work(kfr),nBx,
     &                               Work(jD),nSx,
     &                         0.0d0,Work(kto),nBx)

            lOff=lOff+lnVir(iSym)
            kOff=kOff+nBas(iSym)**2
            jOff=jOff+nSsh(iSym)
            iOff=iOff+lnVir(iSym)**2
         End Do
         Call GetMem('iD_orb','Free','Inte',ip_iD,nOrb)
         kEVir=ipOrbE
*
         EMP2=DeMP2
         DeMP2=0.0d0
         Call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 in truncated virtual space failed !'
           Call Abend
         Endif
         EMP2 = -1.0d0*(EMP2-DeMP2)
      EndIf
      Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
      Call GetMem('Eorb','Free','Real',ipOrbE,4*nOrb)
*
      CALL GETMEM('LCMO','FREE','REAL',LCMO,2*NCMO)
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine FnoCASPT2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            ip_X,ip_Y)
C
C     Purpose: put info in MP2 common blocks.
C
#include "implicit.fh"
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Integer ip_X, ip_Y
C
#include "corbinf.fh"
#include "chomp2_cfg.fh"
C
C
      nSym = mSym
C
      Do iSym = 1,nSym
         nOrb(iSym) = lnOrb(iSym)
         nOcc(iSym) = lnOcc(iSym)
         nFro(iSym) = lnFro(iSym)
         nDel(iSym) = lnDel(iSym)
         nExt(iSym) = lnVir(iSym)
      End Do
C
      ChoAlg=2
      DecoMP2=Decom_Def
      ThrMP2=-9.9D9
      SpanMP2=Span_Def
      MxQualMP2=MxQual_Def
      ChkDecoMP2=.false.
      ForceBatch=.false.
      Verbose=.false.
      SOS_mp2=.false.
      set_cd_thr=.true.
      OED_Thr=1.0d-8
      C_os=1.3d0
      EOSMP2=0.0d0
C
      DoFNO=.true.
      ip_Dab=ip_X
      ip_Dii=ip_Y
      l_Dab=nExt(1)
      l_Dii=nOcc(1)
      Do iSym=2,nSym
         l_Dab=l_Dab+nExt(iSym)**2
         l_Dii=l_Dii+nOcc(iSym)
      End Do
C
      Return
      End

************************************************************************
*                                                                      *
************************************************************************
      Subroutine Check_Amp(nSym,nOcc,nVir,iSkip)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nOcc(nSym), nVir(nSym), iSkip
      Integer nT1amTot, nT1am(8)

      MulD2h(i,j)=iEor(i-1,j-1) + 1

      iSkip=0
      nT1amTot=0
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            nT1am(iSym) = nT1am(iSym)
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
         nT1amTot = nT1amTot + nT1am(iSym)
      End Do

      If (nT1amTot .gt. 0) iSkip=1
      Return
      End
