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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine NewOrb_SCF(Fock,nFock,CMO,nCMO,FOVMax,EOrb,nEOrb,Ovlp,
     &                      nFO,AllowFlip,Scram,nD)
************************************************************************
*                                                                      *
*     purpose: Diagonalize Fock matrix to get new orbitals             *
*                                                                      *
*     input:                                                           *
*       Fock    : Fock matrix of length nFock                          *
*       CMO     : orthonormal vectors from previous iteration of       *
*                 length nCMO                                          *
*                                                                      *
*     output:                                                          *
*       CMO     : orthonormal vectors in current iteration             *
*       FOVMax  : Max Element of occ/virt block in Fock Matrix,        *
*                 transformed into MO basis                            *
*       EOrb    : orbital energies of length nEOrb                     *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: ModFck, PickUp, SortEig, ChkOrt                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*     modified by M.Schuetz                                            *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Integer nFock,nCMO,nEOrb
      Real*8 Fock(nFock,nD),CMO(nCMO,nD),EOrb(nEOrb,nD),FOVMax
      Real*8, Dimension(:), Allocatable:: eConstr, FckM, FckS, HlfF,
     &                                    TraF, Scratch, Temp
      Integer, Dimension(:), Allocatable:: iFerm
      Logical AllowFlip
      Logical Scram, em_On
*
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
      Real*8 Ovlp(nFO)
      Save iSeed
      Data iSeed/13/
*
      Real*8 Fia
      Integer iCMO,iiBT,jEOr,iptr,nOrbmF,nOccmF,nVrt,ia,ij
*                                                                      *
************************************************************************
*                                                                      *
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _SPECIAL_DEBUG_
#ifdef _SPECIAL_DEBUG_
      Call DebugCMOx(CMO,nCMO,nD,nBas,nOrb,nSym,'NewOrb: CMO old')
#endif

*define _DEBUG_
#ifdef _DEBUG_
      Call NrmClc(Fock,nFock*nD,'NewOrb','Fock')
#endif
*
      nSdg=1
      If (Do_SpinAV) nSdg=2
      If (MxConstr.gt.0) Call mma_allocate(eConstr,nSdg*MxConstr,
     &                                     Label='eConstr')
*
      If(.not.AllowFlip) Then
         If(.not.DoHLgap) Then
            DoHLgap=.true.
            HLgap=FlipThr
         End If
      End If
*---- Allocate memory for orbital homeing
      If (.Not.FckAuf) Then
         Call mma_allocate(Scratch,MaxBas**2,Label='Scratch')
         Call mma_allocate(Temp,MaxBas**2,Label='TempX')
      End If
*---- Allocate memory for modified Fock matrix
      Call mma_allocate(FckM,nBT,Label='FckM')
*---- Allocate memory for squared Fock matrix
      Call mma_allocate(FckS,MaxBas**2,Label='FckSX')
*---- Allocate memory for half-transformed Fock matrix
      Call mma_allocate(HlfF,MaxBOF,Label='HlfF')
*---- Allocate memory for transformed Fock matrix (triangular)
      Call mma_allocate(TraF,MaxOrF*(MaxOrF + 1)/2,Label='TraF')
*---- Allocate memory for fermi index array
      Call mma_allocate(iFerm,nBB,Label='iFerm')
      Call Get_iArray('Fermion IDs',iFerm,nnB)
      iChk=0
      Do iiB = 1, nnB
         iChk=iChk+iFerm(iiB)
      End Do
      em_On=iChk.ne.0.and.iChk.ne.nnB
*
      FOVMax = Zero
      WarnCfg=.False.
      Do iD = 1, nD
*
*---- Modify Fock matrix
      call dcopy_(nBT,Fock(1,iD),1,FckM,1)
      If (nnFr.gt.0)
     &   Call ModFck(FckM,Ovlp,nBT,CMO(1,iD),nBO,nOcc(1,iD))
*---- Prediagonalize Fock matrix
      iAddGap = 0
      GapAdd  = 0.0d0
      If(.not.Aufb.and.DoHLgap) Then
         ij   = 1
         iCMO = 1
         jEOr = 1
         Ehomo = -1.0d6
         Elumo =  1.0d6
         Do iSym = 1, nSym
           iiBT=nBas(iSym)*(nBas(iSym)+1)/2
           nOrbmF=nOrb(iSym)-nFro(iSym)
           nOccmF=nOcc(iSym,iD)-nFro(iSym)
           nVrt=nOrb(iSym)-nOcc(iSym,iD)
           iCMO = iCMO + nBas(iSym)*nFro(iSym)
           jEOr = jEOr + nFro(iSym)
           Call Square(FckM(ij),FckS,1,nBas(iSym),nBas(iSym))
           If (nOccmF.gt.0) Then
              Call DGEMM_('N','N',
     &                    nBas(iSym),nOccmF,nBas(iSym),
     &                    1.0d0,FckS,nBas(iSym),
     &                    CMO(iCMO,iD),nBas(iSym),
     &                    0.0d0,HlfF,nBas(iSym))
              Call MxMt(CMO(iCMO,iD),   nBas(iSym),1,
     &                  HlfF,1,nBas(iSym),
     &                  TraF,
     &                  nOccmF,nBas(iSym))
#ifdef _DEBUG_
*             Call Triprt('Occupied Fock matrix in MO basis',
*    &                    '(20F10.4)',TraF,nOccmF)
#endif
              nOccmF=nOccmF-nConstr(iSym)
              Call NIdiag(TraF,CMO(iCMO,iD),nOccmF,nBas(iSym),0)
              nOccmF=nOccmF+nConstr(iSym)
#ifdef _DEBUG_
*             Call Triprt('Occupied Fock matrix in MO basis',
*    &                    '(20F10.4)',TraF,nOccmF)
#endif
              Do iBas = 1,nOccmF
                 ind=iBas*(iBas+1)/2
                 Ehomo=Max(Ehomo,TraF(ind))
              End Do
           End If
*
           iCMO = iCMO + nOccmF*nBas(iSym)
           jEOr = jEOr + nOccmF
           If (Do_SpinAV) Then
              nVrt=nVrt-nConstr(iSym)
              iCMO = iCMO + nConstr(iSym)*nBas(iSym)
              jEOr = jEOr + nConstr(iSym)
           EndIf
           If(nVrt.gt.0) Then
              Call DGEMM_('N','N',
     &                    nBas(iSym),nVrt,nBas(iSym),
     &                    1.0d0,FckS,nBas(iSym),
     &                    CMO(iCMO,iD),nBas(iSym),
     &                    0.0d0,HlfF,nBas(iSym))
              Call MxMt(CMO(iCMO,iD),   nBas(iSym),1,
     &                  HlfF,1,nBas(iSym),
     &                  TraF,
     &                  nVrt,nBas(iSym))
#ifdef _DEBUG_
*             Call Triprt('Virtual Fock matrix in MO basis',
*    &                    '(20F10.4)',TraF,nVrt)
#endif
              Call NIdiag(TraF,CMO(iCMO,iD),nVrt,nBas(iSym),0)
#ifdef _DEBUG_
*             Call Triprt('Virtual Fock matrix in MO basis',
*    &                    '(20F10.4)',TraF,nVrt)
#endif
              Do iBas = 1,nVrt
                 ind=iBas*(iBas+1)/2
                 Elumo=Min(Elumo,TraF(ind))
              End Do
           End If
           If (Do_SpinAV) Then
              nVrt=nVrt+nConstr(iSym)
              iCMO = iCMO - nConstr(iSym)*nBas(iSym)
              jEOr = jEOr - nConstr(iSym)
           EndIf
           iCMO = iCMO + nVrt*nBas(iSym)
           jEOr = jEOr + nVrt
           ij   = ij   + iiBT
         End Do
C        Write(6,'(a,F12.6)') 'E(homo)   ',Ehomo
C        Write(6,'(a,F12.6)') 'E(lumo)   ',Elumo
C        Write(6,'(a,F12.6)') 'E(gap)    ',Elumo-Ehomo
         WarnCfg=Elumo-Ehomo.lt.0.0d0.or.WarnCfg
         If(Elumo-Ehomo.lt.HLgap) Then
            iAddGap=1
            GapAdd=HLgap-Elumo+Ehomo
C           Write(6,'(a,F12.6)') 'E(add)    ',GapAdd
         End If
      End If
*---- Diagonalize Fock matrix in non-frozen molecular basis
      ij   = 1
      iCMO = 1
      jEOr = 1
      iOvlpOff = 1
      Do iSym = 1, nSym
        iiBT=nBas(iSym)*(nBas(iSym)+1)/2
        nOrbmF=nOrb(iSym)-nFro(iSym)
        nOccmF=nOcc(iSym,iD)-nFro(iSym)
        nVrt=nOrb(iSym)-nOcc(iSym,iD)
*------ Find the proper pointers to CMO and EOr
        iCMO = iCMO + nBas(iSym)*nFro(iSym)
        jEOr = jEOr + nFro(iSym)
        If (nOrbmF.gt.0) Then
           Call Square(FckM(ij),FckS,1,nBas(iSym),nBas(iSym))
*--------- Transform Fock matrix to the basis from previous iteration
           Call DGEMM_('N','N',
     &                 nBas(iSym),nOrbmF,nBas(iSym),
     &                 1.0d0,FckS,nBas(iSym),
     &                 CMO(iCMO,iD),nBas(iSym),
     &                 0.0d0,HlfF,nBas(iSym))
           Call MxMt(CMO(iCMO,iD),   nBas(iSym),1,
     &               HlfF,1,nBas(iSym),
     &               TraF,
     &               nOrbmF,nBas(iSym))
c             Call Triprt('Case3 Fock matrix in MO basis',
c    &                    '(20F10.4)',TraF,nOrbmF)
*
*--------- Constrained SCF section begins --------------
           kConstr=1
           Do iConstr=nConstr(iSym),1,-1
              nj=nOccmF-iConstr
              ifc=1+nj*(nj+1)/2
              eConstr(kConstr)=TraF(ifc+nj)
              Call FZero(TraF(ifc),nj)
              Do j=nj+1,nOrbmF-1
                 jj=1+j*(j+1)/2+nj
                 Traf(jj)=0.0d0
              End Do
              Traf(ifc+nj)=-0.666d6*dble(iConstr) ! for sorting
              kConstr=kConstr+1
           End Do
           If (Do_SpinAV) Then
              Do iConstr=0,nConstr(iSym)-1
                 nj=nOccmF+iConstr
                 ifc=1+nj*(nj+1)/2
                 eConstr(kConstr)=TraF(ifc+nj)
                 Call FZero(TraF(ifc),nj)
                 Do j=nj+1,nOrbmF-1
                    jj=1+j*(j+1)/2+nj
                    TraF(jj)=0.0d0
                 End Do
                 TraF(ifc+nj)=0.666d6*dble(kConstr) ! for sorting
                 kConstr=kConstr+1
              End Do
           EndIf
*--------- Constrained SCF section ends ----------------
*
*
*          get max element of Fock matrix
c          Do 400 iBas = 2, nOrbmF
c             Do 401 jBas = 1, iBas-1
c                ijBas = iBas*(iBas-1)/2 + jBas -1 + 1
c                FOVMax=Max(Abs(TraF(ijBas)),FOVMax)
c401          Continue
c400       Continue
*          get max element of occ/virt Block of Fock matrix
*
           If (Teee) Then
              Do iBas = 2, nOrbmF
                 Do jBas = 1, iBas-1
                    ijBas = iBas*(iBas-1)/2 + jBas
                    FOVMax=Max(Abs(TraF(ijBas)),FOVMax)
                 End Do
              End Do
           Else If ((nOccmF.gt.0).AND.(nVrt.gt.0)) Then
              iptr=1+nOccmF*(nOccmF+1)/2
              Do ia=1,nVrt
                 Fia=abs(TraF(iptr+IDAMAX_(nOccmF,TraF(iptr),1)-1))
                 FOVMax=Max(Fia,FOVMax)
                 iptr=iptr+nOccmF+ia
              End Do
           End If
*
*--------- Modify Fock matrix to enhance convergence
c          Call Triprt('Fock matrix in MO basis before modification',
c    &                 '(20F10.4)',TraF,nOrbmF)
*
*--- Add to homo lumo gap
*
           If(iAddGap.eq.1) Then
              Do iOrb=nOccmF+1,nOrbmF
                 ind=iOrb*(iOrb+1)/2
                 TraF(ind)=TraF(ind)+GapAdd
              End Do
           End If
           ind=1
           Do iOrb=1,nOrbmF
              Do jOrb=1,iOrb
*--- Scale OV elements of fock matrix
                 If(iOrb.gt.nOcc(iSym,iD) .and.
     &              jOrb.le.nOcc(iSym,iD)) Then
                    TraF(ind)=RotFac*TraF(ind)
                 End If
*--- Levelshift virtual diagonal matrix elements
                 If(iOrb.gt.nOcc(iSym,iD) .and. iOrb.eq.jOrb) Then
                    TraF(ind)=TraF(ind)+RotLev
                 End If
                 ind=ind+1
              End Do
           End Do
*--- Add scrambling
           If(Scram) Then
              ind=1
              Do iOrb=1,nOrbmF
                 Do jOrb=1,iOrb
                    If(iOrb.ne.jOrb) Then
                       q=ScrFac*(2.0d0*Random_Molcas(iSeed)-1.0d0)
                       TraF(ind)=TraF(ind)+q
                    End If
                    ind=ind+1
                 End Do
              End Do
           End If

*--- Scale OV elements if too big.
           indii=1
           indij=1
           Do iOrb=1,nOrbmF
              indjj=1
              Do jOrb=1,iOrb
                 If(iOrb.gt.nOcc(iSym,iD) .and.
     &              jOrb.le.nOcc(iSym,iD)) Then
                    tmp1=Max(Abs(TraF(indii)-TraF(indjj)),1.0d-3)
                    tmp2=Abs(TraF(indij)/tmp1)
                    If(tmp2.gt.RotMax .and.
     &                 Abs(TraF(indij)).gt.0.001d0) Then
                       TraF(indij)=TraF(indij)*RotMax/tmp2
                    End If
                 End If
                 indij=indij+1
                 indjj=indjj+jOrb+1
              End Do
              indii=indii+iOrb+1
           End Do
*
*--------- Constrained SCF section begins --------------
           kConstr=1
           Do iConstr=nConstr(iSym),1,-1
              nj=nOccmF-iConstr
              ifc=1+nj*(nj+1)/2
              Call FZero(TraF(ifc),nj)
              Do j=nj+1,nOrbmF-1
                 jj=1+j*(j+1)/2+nj
                 TraF(jj)=0.0d0
              End Do
              TraF(ifc+nj)=-0.666d6*dble(iConstr) ! for sorting
              kConstr=kConstr+1
           End Do
           If (Do_SpinAV) Then
              Do iConstr=0,nConstr(iSym)-1
                 nj=nOccmF+iConstr
                 ifc=1+nj*(nj+1)/2
                 Call FZero(TraF(ifc),nj)
                 Do j=nj+1,nOrbmF-1
                    jj=1+j*(j+1)/2+nj
                    TraF(jj)=0.0d0
                 End Do
                 TraF(ifc+nj)=0.666d6*dble(kConstr) ! for sorting
                 kConstr=kConstr+1
              End Do
           EndIf
*--------- Constrained SCF section ends ----------------
*
*
#ifdef _DEBUG_
*          Call Triprt('Fock matrix in MO basis after modification',
*    &                 '(10F10.4)',TraF,nOrbmF)
#endif

*
*--------- Diagonalize and form orbital energies
*
            Dummy=0.0D0
            iDum=0
*
*           Store the original CMOs for root following.
*
            Call DCopy_(nBas(iSym)**2,CMO(iCMO,iD),1,FckS,1)
*
            Call Diag_Driver('V','A','L',nOrbmF,TraF,
     &                       TraF,nOrbmF,Dummy,Dummy,iDum,
     &                       iDum,EOrb(jEOr,iD),CMO(iCMO,iD),
     &                       nBas(iSym),0,-1,
     &                       'J',nFound,iErr)
*
*           Fix standard phase pf the orbitals
*
            Do i = 1, nBas(iSym)
               tmp = OrbPhase(CMO(iCMO+(i-1)*nBas(iSym),iD),nBas(iSym))
            End Do
*define _DEBUG_
#ifdef _DEBUG_
            Call NrmClc(Fcks,nbas(iSym)*nOrb(iSym),'NewOrb','Old CMOs')
            Call NrmClc(CMO(iCMO,iD),nbas(iSym)*nOrb(iSym),
     &                  'NewOrb','New CMOs')
*           Call RecPrt('Old CMOs',' ',FckS,nBas(iSym),nOrb(iSym))
*           Call RecPrt('New CMOs',' ',CMO(iCMO,iD),nBas(iSym),
*    &                                              nOrb(iSym))
#endif
*                                                                      *
************************************************************************
*                                                                      *
*           Reorder the orbitals to preserve e/m partition.
*
            If (.Not.em_On) Go To 110 ! skip if all orbitals are of
                                      ! the same type.
*
            Do iOrb = 1, nOrb(iSym)-1    ! Loop over the old orbitals
*
*              Compute a check sum which is not zero if the orbital
*              is a muonic orbital.
*
               tmp=0.0D0
               Do kBas = 0, nBas(iSym)-1
                  tmp = tmp
     &                + DBLE(iFerm(jEOr+kBas))
     &                * ABS(FckS((iOrb-1)*nBas(iSym)+kBas+1))
               End Do
               Muon_i=0                  ! electronic
               If (tmp.ne.0.0D0) Muon_i=1! muonic
*              Write (6,*) 'iOrb,Muon_i,tmp=',iOrb,Muon_i,tmp
*
*              Loop over the new orbitals and test if it is of the
*              same type. i.e. fermionic or electronic.
*
               Do jOrb = iOrb, nOrb(iSym)
                  tmp=0.0D0
                  Do kBas = 0, nBas(iSym)-1
                     tmp = tmp
     &                   + DBLE(iFerm(jEOr+kBas))
     &                   * ABS(CMO(iCMO+(jOrb-1)*nBas(iSym)+kBas,iD))
                  End Do
*
                  Muon_j=0                  ! electronic
                  If (tmp.ne.0.0D0) Muon_j=1! muonic
*
*                 If orbital and fermion index are identical fine.
*
                  If (iOrb.eq.jOrb .and. Muon_i.eq.Muon_j) Go To 678
*
*                 If fermion index the same swap orbital and the
*                 corresponding orbital energy in the new list.
*
                  If (Muon_i.eq.Muon_j) Then
                     Tmp=EOrb(jEOr-1+iOrb,iD)
                     EOrb(jEOr-1+iOrb,iD)=EOrb(jEOr-1+jOrb,iD)
                     EOrb(jEOr-1+jOrb,iD)=Tmp
                     Call DSwap_(nBas(iSym),
     &                          CMO(iCMO+(iOrb-1)*nBas(iSym),iD),1,
     &                          CMO(iCMO+(jOrb-1)*nBas(iSym),iD),1)
                     Go To 678
                  End If
*
               End Do   !  jOrb
*
*              Arrive at this point when the iOrb'th orbital in the
*              new list is of the same type as in the old list.
*
 678           Continue
*
            End Do      !  iOrb
 110        Continue
*                                                                      *
************************************************************************
*                                                                      *
*           Order the orbitals in the same way as previously, that is,
*           do not populate according to the aufbau principle.
*
            If (FckAuf) Go To 120
*           Write (6,*) 'Follow the orbitals'
*
*           Form  C^+ S  C for the old orbitals
*
            Call FZero(Scratch,nOrb(iSym)*nBas(iSym))
            Call Square(Ovlp(iOvlpOff),Temp,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('T','N',
     &                 nOrb(iSym),nBas(iSym),nBas(iSym),
     &                 1.0D0,FckS,nBas(iSym),
     &                       Temp,nBas(iSym),
     &                 0.0D0,Scratch,nOrb(iSym))
*
            Do iOrb = 1, nOrb(iSym)-1  ! Loop over the old orbitals
*
*              Don't apply non-aufbau principle for the electrons!
*
               iOff = (iOrb-1)*nBas(iSym) + iCMO
*
*              Idenify orbital type.
*
               tmp=0.0D0
               Do kBas = 0, nBas(iSym)-1
                  tmp = tmp
     &                + DBLE(iFerm(jEOr+kBas))
     &                * ABS(FckS((iOrb-1)*nBas(iSym)+kBas+1))
               End Do
               Muon_i=0                  ! electronic
               If (tmp.ne.0.0D0) Muon_i=1! muonic
               If (Muon_i.eq.0) Cycle
*              Write (6,*) 'iOrb,Muon_i=',iOrb,Muon_i
*
               kOrb=0
               Tmp0= 0.0D0
               Do jOrb = 1, nOrb(iSym)   ! Loop over the new orbitals
                  jOff = (jOrb-1)*nBas(iSym) + iCMO
                  Tmp1=Abs(DDot_(nBas(iSym),Scratch(iOrb),nOrb(iSym),
     &                                      CMO(jOff,iD),1))
                  If (Tmp1.gt.Tmp0) Then
                     Tmp0=Tmp1
                     kOrb=jOrb
                  End If
               End Do
*              Write (6,*) 'kOrb,Tmp0=',kOrb,Tmp0
*
               If (iOrb.ne.kOrb) Then
                  ii = iOrb + jEOr - 1
                  kk = kOrb + jEOr - 1
                  tmp = EOrb(ii,iD)
                  EOrb(ii,iD)=EOrb(kk,iD)
                  EOrb(kk,iD)=tmp
                  kOff = (kOrb-1)*nBas(iSym) + iCMO
                  Call DSwap_(nBas(iSym),CMO(iOff,iD),1,
     &                                   CMO(kOff,iD),1)
               End If
            End Do
#ifdef _DEBUG_
            Call NrmClc(CMO(iCMO,iD),nbas(iSym)*nOrb(iSym),
     &                  'NewOrb','New CMOs')
*           Call RecPrt('New CMOs',' ',CMO(iCMO,iD),nBas(iSym),
*    &                                              nOrb(iSym))
#endif
 120        Continue
*                                                                      *
************************************************************************
*                                                                      *
*
*--------- Constrained SCF section begins --------------
*
            If (nConstr(iSym).gt.0) Then
               Do kConstr=1,nConstr(iSym)
                  iConstr=jEOr-1+kConstr
                  EOrb(iConstr,iD)=0.666d6*dble(kConstr)
               End Do
               Call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),nOccmF,
     &                      nBas(iSym))
               iConstr=1
               Do kConstr=nConstr(iSym),1,-1
                  lConstr=jEOr+nOccmF-kConstr
                  EOrb(lConstr,iD)=eConstr(iConstr)
                  iConstr=iConstr+1
               End Do
               If (Do_SpinAV) Then
                  Do kConstr=nConstr(iSym),1,-1
                     iConstr=jEOr+nOrbmF-kConstr
                     EOrb(iConstr,iD)=-0.666d6*dble(kConstr)
                  End Do
                  kCMO = iCMO + nOccmF*nBas(iSym)
                  kEOr = jEOr + nOccmF
                  Call SortEig(EOrb(kEOr,iD),CMO(kCMO,iD),nVrt,
     &                     nBas(iSym))
                  iConstr=nConstr(iSym)+1
                  Do kConstr=1,nConstr(iSym)
                     lConstr=kEOr+kConstr
                     EOrb(lConstr,iD)=eConstr(iConstr)
                     iConstr=iConstr+1
                  End Do
               EndIf
            EndIf
*--------- Constrained SCF section ends ----------------
        End If
*------ Update pointers
        iCMO = iCMO + nOrbmF*nBas(iSym)
        jEOr = jEOr + nOrbmF
        ij   = ij   + iiBT
        iOvlpOff = iOvlpOff + nBas(iSym)*(nBas(iSym)+1)/2
      End Do
*
*---- Check orthogonality
      Call ChkOrt(CMO(1,iD),nBO,Ovlp,nBT,Whatever)
*
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(iFerm)
      Call mma_deallocate(TraF)
      Call mma_deallocate(HlfF)
      Call mma_deallocate(FckS)
      Call mma_deallocate(FckM)
      If (.Not.FckAuf) Then
         Call mma_deallocate(Scratch)
         Call mma_deallocate(Temp)
      End If
      If (MxConstr.gt.0) Call mma_deallocate(eConstr)
*
#ifdef _SPECIAL_DEBUG_
      Call DebugCMOx(CMO,nCMO,nD,nBas,nOrb,nSym,'NewOrb: CMO new')
#endif
*define _DEBUG_
#ifdef _DEBUG_
*     Do iD = 1, nD
*        iOff=1
*        jOff=1
*        Do iSym = 1, nSym
*           Call RecPrt('CMO',' ',CMO(jOff,iD),
*    &                  nBas(iSym),nOrb(iSym))
*           iOff = iOff + nOrb(iSym)
*           jOff = jOff + nBas(iSym)*nOrb(iSym)
*        End Do
*     End Do
#endif

      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 8) = TimFld( 8) + (Cpu2 - Cpu1)
      Return
      End
