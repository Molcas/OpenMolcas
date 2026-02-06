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
      use definitions, only: iwp, wp
      use InputData, only: Input
      use ChoMP2, only: DeMP2, MP2_small, shf
      use stdalloc, only: mma_allocate, mma_deallocate
      use Molcas, only: MxBas
      Implicit Real*8 (A-H,O-Z)
*
      Integer(kind=iwp), intent(out):: irc
      Integer(kind=iwp), intent(in):: nSym
      Integer(kind=iwp), intent(in):: nBas(nSym),nFro(nSym),nIsh(nSym),
     &                                nAsh(nSym)
      Integer(kind=iwp), intent(inout):: nSsh(nSym),nDel(nSym)
      real(kind=wp), intent(in)::  vfrac
      Integer(kind=iwp), intent(inout):: IFQCAN
      real(kind=wp), intent(out)::  EMP2
      Logical(kind=iwp), intent(in):: DoMP2
      Integer(kind=iwp), intent(in):: NCMO
      real(kind=wp), intent(inout):: CMO(NCMO)
*
      Integer(kind=iwp) ns_V(8), nAct(8)
      Integer(kind=iwp) lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      real(kind=wp)  TrDP(8), TrDF(8)
      real(kind=wp), allocatable:: CMOX(:), DMAT(:), OrbE(:)
      Integer(kind=iwp), allocatable:: ID(:)
*
*
      irc=0
      MP2_small=.false.
      shf=Input%RegFNO
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
     &   'No pseudocanonical RASSCF orbitals found! '
     & //'I will proceed with FDIAG values.'
      EndIf
      CALL mma_allocate(CMOX,2*NCMO,Label='CMOX')
      iCMO=1+NCMO
* This is not the best solution, but I wanted to avoid having to rewrite
* the indexing code below just to use the CMO array directly
      call dcopy_(NCMO,CMO,1,CMOX,1)
      Call mma_allocate(OrbE,4*nOrb,Label='OrbE')
      ipOrbE=1
      Call Get_darray('RASSCF OrbE',OrbE(ipOrbE),nOrb)
      iAoff=0
      Do iSym=1,nSym
         ipOrbE_=ipOrbE+iAoff+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (OrbE(ipOrbE_+k).lt.0.0d0) nAct(iSym)=nAct(iSym)+1
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
         call dcopy_(lnOcc(iSym),OrbE(ifr),1,OrbE(ito),1)
         ifr=ipOrbE+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),OrbE(ifr),1,OrbE(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do
      Call mma_allocate(Dmat,nVV+nOA,LABEL='DMAT')
      DMAT(:)=0.0D0
      ip_X = 1
      ip_Y = ip_X + nVV
*
      Call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
      Call FZero(CMOX(iCMO),NCMO)
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMOX(kfr),1,CMOX(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMOX(kfr),1,CMOX(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,CMOX(iCMO),OrbE(kEOcc),OrbE(kEVir),
     &                   DMAT(ip_X),DMAT(ip_Y))
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend()
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend()
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
           Call Eigen_Molcas(nSsh(iSym),DMAT(jD),OrbE(ip_Z),OrbE(ip_ZZ))
*     Reorder to get relevant eigenpairs first
           Do j=1,nSsh(iSym)/2
              Do i=1,nSsh(iSym)
                 lij=jD-1+nSsh(iSym)*(j-1)+i
                 kij=jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
                 tmp=DMAT(lij)
                 DMAT(lij)=DMAT(kij)
                 DMAT(kij)=tmp
              End Do
              tmp=OrbE(ip_Z-1+j)
              OrbE(ip_Z-1+j)=OrbE(ip_Z+nSsh(iSym)-j)
              OrbE(ip_Z+nSsh(iSym)-j)=tmp
           End Do
*
*     Compute new MO coeff. : X=C*U
           kfr=iCMO+jOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
           kto=1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
           Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                        1.0d0,CMOX(kfr),nBas(iSym),
     &                              DMAT(jD),nSsh(iSym),
     &                        0.0d0,CMOX(kto),nBas(iSym))
           iOff=iOff+nSsh(iSym)**2
           TrDF(iSym)=ddot_(nSsh(iSym),OrbE(ip_Z),1,[1.0d0],0)
           If (vfrac.ge.0.0d0) Then
              ns_V(iSym)=int(vfrac*dble(nSsh(iSym)))
              TrDP(iSym)=ddot_(ns_V(iSym),OrbE(ip_Z),1,[1.0d0],0)
           Else
              ns_V(iSym)=nSsh(iSym)-1
              TrDP(iSym)=ddot_(ns_V(iSym),OrbE(ip_Z),1,[1.0d0],0)
              Delta_TrD=TrDP(iSym)-TrDF(iSym) ! this is negative
              Delta_TrD=Delta_TrD/TrDF(iSym)
              Do While (Delta_TrD.gt.vfrac)
                 ns_V(iSym)=ns_V(iSym)-1
                 TrDP(iSym)=ddot_(ns_V(iSym),OrbE(ip_Z),1,[1.0d0],0)
                 Delta_TrD=(TrDP(iSym)-TrDF(iSym))/TrDF(iSym)
              End Do
           End If
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
      call dcopy_(NCMO,CMOX,1,CMO,1)
*
*     Reset MP2_small for this new call to ChoMP2_Drv
      Call Check_Amp(nSym,lnOcc,nSsh,iSkip)
      MP2_small = DoMP2 .and. iSkip.gt.0
      If (MP2_small) Then
*
         Call FnoCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,nDel,nSsh)
*
         Call mma_allocate(iD,nOrb,Label='iD')
         Do k=1,nOrb
            iD(k) = k
         End Do
         lOff=0
         kOff=0
         jOff=0
         iOff=0
         Do iSym=1,nSym  ! canonical orb. in the reduced virtual space
            jD=ip_X+iOff
            Call Get_Can_Lorb(OrbE(kEVir+lOff),OrbE(ipOrbE+jOff),
     &                        nSsh(iSym),lnVir(iSym),
     &                        iD,DMAT(jD))

            kfr=1+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
            kto=iCMO+kOff+nBas(iSym)*(nFro(iSym)+lnOcc(iSym))
            nBx=Max(1,nBas(iSym))
            nSx=Max(1,nSsh(iSym))
            Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),
     &                         1.0d0,CMOX(kfr),nBx,
     &                               DMAT(jD),nSx,
     &                         0.0d0,CMOX(kto),nBx)

            lOff=lOff+lnVir(iSym)
            kOff=kOff+nBas(iSym)**2
            jOff=jOff+nSsh(iSym)
            iOff=iOff+lnVir(iSym)**2
         End Do
         Call mma_deallocate(iD)
         kEVir=ipOrbE
*
         EMP2=DeMP2
         DeMP2=0.0d0
         Call ChoMP2_Drv(irc,Dummy,CMOX(iCMO),OrbE(kEOcc),OrbE(kEVir),
     &                   DMAT(ip_X),DMAT(ip_Y))
         If(irc.ne.0) then
           write(6,*) 'MP2 in truncated virtual space failed !'
           Call Abend
         Endif
         EMP2 = -1.0d0*(EMP2-DeMP2)
      EndIf
      Call mma_deallocate(Dmat)
      Call mma_deallocate(OrbE)
*
      CALL mma_deallocate(CMOX)
*
      End SUBROUTINE FNO_CASPT2
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine FnoCASPT2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
C
C     Purpose: put info in MP2 common blocks.
C
      use definitions, only: iwp
      use constants, only: Zero
      Use ChoMP2, only: C_os, ChkDecoMP2, ChoAlg, Decom_Def, DecoMP2,
     &                  DoFNO, EOSMP2, ForceBatch, l_Dii, MxQual_Def,
     &                  MxQualMP2, OED_Thr, set_cd_thr, SOS_mp2,
     &                  Span_Def, SpanMP2, ThrMP2, Verbose
      use cOrbInf, only: nSym, nOrb, nOcc, nFro, nDel, nExt

      Implicit None
      Integer(kind=iwp), intent(in):: mSym
      Integer(kind=iwp), intent(in)::  lnOrb(8), lnOcc(8), lnFro(8),
     &                                 lnDel(8), lnVir(8)

      Integer(kind=iwp) :: iSym
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
      EOSMP2=Zero
C
      DoFNO=.true.
      l_Dii=nOcc(1)
      Do iSym=2,nSym
         l_Dii=l_Dii+nOcc(iSym)
      End Do
C
      End SubRoutine FnoCASPT2_putInf

************************************************************************
*                                                                      *
************************************************************************
      Subroutine Check_Amp(nSym,nOcc,nVir,iSkip)
      use definitions, only: iwp
      use SYmmetry_Info, only: Mul

      Implicit None
      integer(kind=iwp), intent(in)::  nSym, nOcc(nSym), nVir(nSym)
      Integer(kind=iwp), intent(out):: iSkip

      integer(kind=iwp) iSym, iSymi, iSyma
      Integer(kind=iwp) nT1amTot, nT1am(8)

      iSkip=0
      nT1amTot=0
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = Mul(iSymi,iSym)
            nT1am(iSym) = nT1am(iSym)
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
         nT1amTot = nT1amTot + nT1am(iSym)
      End Do

      If (nT1amTot .gt. 0) iSkip=1
      End Subroutine Check_Amp
