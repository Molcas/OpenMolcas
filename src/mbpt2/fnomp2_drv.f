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
      SUBROUTINE FNOMP2_Drv(irc,EMP2,CMOI,EOcc,EVir)

#include "implicit.fh"
      Real*8 EMP2, CMOI(*), EOcc(*), EVir(*)
C
      Logical DoDens_
      Integer ChoAlg_
#include "orbinf2.fh"
#include "corbinf.fh"
#include "chomp2_cfg.fh"

      DoDens_= DoDens
      DoDens = .false.
      ChoAlg_= ChoAlg
      ChoAlg = 2

      CALL FNO_MP2(irc,nSym,nBas,nFro,nOcc,nExt,nDel,
     &                      CMOI,EOcc,EVir,vkept,DoMP2,XEMP2)
      If (irc .ne. 0) Then
         Write(6,*) 'FNO_MP2 returned ',irc
         Call SysAbendMsg('FNO_MP2',
     &                    'Non-zero return code from FNO_MP2',
     &                    ' ')
      EndIf

      ChoAlg = ChoAlg_
      DoDens = DoDens_
      DoFNO = .false.
      Call ChoMP2_Drv(irc,EMP2,CMOI,EOcc,EVir)
      EMP2=EMP2+XEMP2

      Return
      End
*****************************************************************************
*                                                                           *
*****************************************************************************

      SUBROUTINE FNO_MP2(irc,nSym,nBas,nFro,nIsh,nSsh,nDel,
     &                       CMOI,EOcc,EVir,vfrac,DoMP2,EMP2)
*****************************************************************************
*                                                                           *
*     Purpose:  setup of Frozen Natural Orbitals MP2 (FNO-MP2).             *
*                                                                           *
*     Author:   F. Aquilante  (Geneva, Nov  2008)                           *
*                                                                           *
*****************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nSsh(nSym),
     &        nDel(nSym),nAuxO(8)
      Real*8  CMOI(*), EOcc(*), EVir(*), EMP2, vfrac
      Logical DoMP2
      Real*8  DeMP2
      Logical MP2_small
      Common / ChFNOPT/ DeMP2, MP2_small
*
      Integer ns_V(8)
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Real*8  TrDF(8), TrDP(8)
*
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
      nOrb=0
      nVV=0
      Do i=1,nSym
        ns_V(i)=0
        nBasT=nBasT+nBas(i)
        nOrb=nOrb+nFro(i)+nIsh(i)+nSsh(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nVV=nVV+nSsh(i)**2
        nBmx=Max(nBmx,nBas(i))
      End Do
      IF(nBasT.GT.mxBas) then
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
      NCMO=nSQ
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,2*NCMO)
      iCMO=LCMO+NCMO
      CALL DCOPY_(NCMO,CMOI,1,WORK(LCMO),1)
*
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)
         nOA=nOA+lnOcc(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnOrb(iSym)=lnOcc(iSym)+lnVir(iSym)
         lnDel(iSym)=nDel(iSym)
      End Do
*
      Call GetMem('Eorb','Allo','Real',ipOrbE,4*nOrb)
      jOff=0
      kOff=0
      lOff=0
      Do iSym=1,nSym
         jp=ipOrbE+lOff+nFro(iSym)
         jOcc=jOff+1
         call dcopy_(nIsh(iSym),EOcc(jOcc),1,Work(jp),1)
         jVir=kOff+1
         jp=jp+nIsh(iSym)
         call dcopy_(nSsh(iSym),EVir(jVir),1,Work(jp),1)
         jOff=jOff+nIsh(iSym)
         kOff=kOff+nSsh(iSym)
         lOff=lOff+nBas(iSym)
      End Do
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
         call dcopy_(nIsh(iSym),Work(ifr),1,Work(ito),1)
         ifr=ipOrbE+ioff+nFro(iSym)+nIsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),Work(ifr),1,Work(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+nIsh(iSym)
         koff=koff+nSsh(iSym)
      End Do
      Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
      ip_Y=ip_X+nVV
      Call FZero(Work(ip_X),nVV+nOA)
*
      Call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y)
      Call FZero(Work(iCMO),NCMO)
      iOff=0
      Do iSym=1,nSym
         kfr=LCMO+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),Work(kfr),1,Work(kto),1)
         kfr=LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),Work(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
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
           kfr=iCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           kto=LCMO+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                        1.0d0,Work(kfr),nBas(iSym),
     &                              Work(jD),nSsh(iSym),
     &                        0.0d0,Work(kto),nBas(iSym))
           iOff=iOff+nSsh(iSym)**2
           TrDF(iSym)=ddot_(nSsh(iSym),Work(ip_Z),1,[1.0d0],0)
           ns_V(iSym)=int(vfrac*nSsh(iSym))
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
        write(6,'(4X,I4,14X,G13.6,5X,G13.6)') iSym,TrDF(iSym),TrDP(iSym)
        STrDF=STrDF+TrDF(iSym)
        STrDP=STrDP+TrDP(iSym)
      End Do
      write(6,*)'------------------------------------------------------'
      write(6,'(A,G13.6,5X,G13.6)')'   Sum :              ',STrDF,STrDP
      write(6,*)'------------------------------------------------------'
*
*     Update the nSsh, nDel for FNO-MP2
      nSsh_t=0
      Do iSym=1,nSym
         lnOrb(iSym)=lnOrb(iSym)-nSsh(iSym)+ns_V(iSym)
         nDel(iSym)=nDel(iSym)+nSsh(iSym)-ns_V(iSym)
         nSsh(iSym)=ns_V(iSym)
         nSsh_t=nSsh_t+nSsh(iSym)
         nAuxO(iSym)=nBas(iSym)-nDel(iSym)
      End Do
      Call Put_iArray('nDelPT',nDel,nSym)
      Call Put_iArray('nOrb',nAuxO,nSym)
*
*     Reset MP2_small for this new call to ChoMP2_Drv
      Call Check_Amp2(nSym,lnOcc,nSsh,iSkip)
      MP2_small = iSkip .gt. 0
      If (MP2_small) Then
*
         Call FnoMP2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh,ip_X,ip_Y)
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

            kfr=LCMO+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
            kto=1+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
            nBx=Max(1,nBas(iSym))
            nSx=Max(1,nSsh(iSym))
            Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                         1.0d0,Work(kfr),nBx,
     &                               Work(jD),nSx,
     &                         0.0d0,CMOI(kto),nBx)

            lOff=lOff+lnVir(iSym)
            kOff=kOff+nBas(iSym)**2
            jOff=jOff+nSsh(iSym)
            iOff=iOff+lnVir(iSym)**2
         End Do
         Call GetMem('iD_orb','Free','Inte',ip_iD,nOrb)
         kEVir=ipOrbE
*
*   Copy the new Evir to output array
         call dcopy_(nSsh_t,Work(kEVir),1,EVir,1)
*
         Write(6,*)
         Write(6,'(A,8I4)')
     & ' Secondary orbitals after selection:',(nSsh(i),i=1,nSym)
         Write(6,'(A,8I4)')
     & ' Deleted orbitals after selection:  ',(nDel(i),i=1,nSym)
         Write(6,*)
         Write(6,*) 'Energies of the active virtual orbitals '
         ii=0
         Do iSym=1,nSym
            If ( nSsh(iSym).ne.0 ) then
               Write(6,*)
               Write(6,'(A,I2,(T40,5F14.6))')
     &            ' symmetry species',iSym,(EVir(ii+k),k=1,nSsh(iSym))
               ii=ii+nSsh(iSym)
            End If
         End Do
         Write(6,*)
*
         EMP2=DeMP2
         DeMP2=0.0d0
         If (DoMP2) Call ChoMP2_Drv(irc,Dummy,CMOI,Work(kEOcc),
     &                                  Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 in truncated virtual space failed !'
           Call Abend
         Endif
         EMP2 = -1.0d0*(EMP2-DeMP2)
         If (DoMP2) Then
            DeMP2=EMP2
c            write(6,*)
c            write(6,'(1x,a,f18.10,a)')'FNO correction:       ',EMP2,
c     &                                '   (estimate)   '
c            write(6,*)
         Else
            EMP2=0.0d0
         EndIf
      Else
         write(6,*)
         write(6,*)'We found  ZERO amplitudes T(ai,bj) with the final '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      EndIf
*
*     Update runfile for subsequent calcs (e.g., CHCC)
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff
         ito=ipOrbE+koff+nFro(iSym)
         call dcopy_(nIsh(iSym),EOcc(ifr),1,Work(ito),1)
         ifr=1+joff
         ito=ito+nIsh(iSym)
         call dcopy_(nSsh(iSym),EVir(ifr),1,Work(ito),1)
         ioff=ioff+nIsh(iSym)
         joff=joff+nSsh(iSym)
         koff=koff+nBas(iSym)
      End Do
      Call Put_dArray('OrbE',Work(ipOrbE),nOrb)
      Call Put_dArray('Last orbitals',CMOI,NCMO)
*
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
      SubRoutine FnoMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                         ip_X,ip_Y)
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
