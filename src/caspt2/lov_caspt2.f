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
      SUBROUTINE Lov_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,NAME,
     &       nUniqAt,Thrs,IFQCAN,DoMP2,DoEnv,all_Vir,EMP2,CMO,NCMO)
************************************************************************
*                                                                      *
* Purpose:  setup of Localized occupied-virtual CASPT2 (LovCASPT2).    *
*           The CASPT2 correction to the energy will later be computed *
*           only for the "active region" of the molecule.              *
*           The MP2 correction due to the remaining frozen region      *
*           is computed here if DoMP2=.true.                           *
*           If DoEnv=.true. we compute the energy of the environment   *
*           as the total MP2 energy minus the MP2 energy of the        *
*           "active region" of the molecule.                           *
*                                                                      *
* Author:   F. Aquilante  (Geneva, Feb. 2008)                          *
*                                                                      *
************************************************************************
      use definitions, only: iwp, wp, u6
      use OneDat, only: sNoNuc, sNoOri
      use Molcas, only: LenIn, LenIn8, MxAtom, MxBas
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      integer(kind=iwp), intent(out):: irc
      integer(kind=iwp), intent(in)::  nSym
      integer(kind=iwp), intent(inout):: nFro(nSym),nIsh(nSym),
     &                                   nSsh(nSym),nDel(nSym)
      integer(kind=iwp), intent(in):: nBas(nSym), nAsh(nSym)
      Character(Len=LENIN8), intent(in):: NAME(*)
      integer(kind=iwp), intent(in):: nUniqAt
      real(kind=wp), intent(in)::  Thrs
      integer(kind=iwp), intent(inout):: IFQCAN
      Logical(kind=iwp), intent(inout):: DoMP2
      Logical(kind=iwp), intent(in):: DoEnv, all_Vir
      real(kind=wp), intent(out)::  EMP2
      integer(kind=iwp), intent(in):: NCMO
      real(kind=wp), intent(inout):: CMO(nCMO)

      Character(Len=LENIN) blank, NamAct(mxAtom)
      character(len=8) :: Label
      Logical(kind=iwp) ortho
      real(kind=wp)  TrA(8), TrF(8), TrX(8)
      integer(kind=iwp) ns_O(8), ns_V(8)
      integer(kind=iwp) lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      integer(kind=iwp), allocatable:: nBas_per_Atom(:), nBas_Start(:),
     &                                 D_A(:), D_Vir(:)
      real(kind=wp), allocatable:: SQ(:), SLT(:), CMOX(:), Q(:), Z(:)
      real(kind=wp), allocatable:: Saa(:), XMO(:), DMat(:), OrbE(:)
      real(kind=wp) Dumm, E2_ab, E2_Aonly, STrA, STrF, STrX, Thrd
      real(kind=wp), external:: DDot_
      integer(kind=iwp) i, iAt, iBat, iCMO, iComp, iDo, ie, ik, iloc,
     &                  iOff, iopt, ip_X, ip_Y, ipAsh, ipCMO, ipEorb,
     &                  ipQa, iQ, iQa, iSkip, iSQ, iSym, iSymLbl, iV,
     &                  jAt, jBas, jBat, jCMO, jDo, jjCMO, jjZ, jOff,
     &                  jQ, jZ, kBas, kEOcc, kEVir, kfr, kOff, kto,
     &                  l_nBas_per_Atom, l_nBas_Start, lBas, lOff, lsq,
     &                  ltri, mAsh, mOff, nActa, nAk, nBasT, nBat, nBk,
     &                  nBmx, nOA, nOrb, nSQ, nTri, nVV, ipOrbE, nBx
*
      irc=0
      EMP2=Zero
      blank='   '
      iDo=0
      jDo=0
      If (DoEnv .and. DoMP2) Then
         Call WarningMessage(1,'Both DoEnv and DoMP2 selected.')
         Write (u6,'(/,A)') ' DoMP2 will be ignored.'
         DoMP2=.false.
      EndIf
      If (all_Vir .and. DoMP2) Then
         Call WarningMessage(1,'Both VirAll and DoMP2 selected.')
         Write (u6,'(/,A)') ' DoMP2 will be ignored.'
         DoMP2=.false.
      EndIf
      Do iSym=1,nSym
         TrA(iSym)=0
         TrF(iSym)=0
         TrX(iSym)=0
      End Do
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
      Do i=1,nSym
        nBasT=nBasT+nBas(i)
        nOrb=nOrb+nFro(i)+nIsh(i)+nAsh(i)+nSsh(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nBmx=Max(nBmx,nBas(i))
        mAsh=Max(mAsh,nAsh(i))
      End Do
      IF(nBasT.GT.mxBas) then
       Write(u6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
*     nUniqAt = # of symm. unique atoms. Initialize NamAct to blanks.
*     ---------------------------------------------------------------

      If (nUniqAt.lt.1 .or. nUniqAt.gt.MxAtom) Then
         Write(u6,'(A,I9)') 'nUniqAt =',nUniqAt
         Call Abend()
      End If
      Do iAt=1,nUniqAt
         NamAct(iAt)=blank
      End Do

C     Allocate and get index arrays for basis functions per atom.
C     -----------------------------------------------------------

      l_nBas_per_Atom = nUniqAt
      l_nBas_Start    = nUniqAt
      Call mma_allocate(nBas_per_Atom,l_nBas_per_Atom,Label='nB/A')
      Call mma_allocate(nBas_Start,l_nBas_Start,Label='nBStart')
*
*----------------------------------------------------------------------*
*     Read the overlap matrix                                          *
*----------------------------------------------------------------------*
      CALL mma_allocate(SQ,nSQ,Label='SQ')
      CALL mma_allocate(SLT,nTri,Label='SLT')
      isymlbl=1
      iopt=ibset(ibset(0,sNoOri),sNoNuc)
      Label='Mltpl  0'
      iComp=1
      Call RdOne(irc,iopt,Label,iComp,SLT,isymlbl)
      If(irc.ne.0) return
      ltri=1
      lsq=1
      Do iSym=1,nSym
         Call Square(SLT(ltri),SQ(lsq),1,nBas(iSym),
     &                                               nBas(iSym))
         ltri=ltri+nBas(iSym)*(nBAs(iSym)+1)/2
         lsq=lsq+nBas(iSym)**2
      End Do
      CALL mma_deallocate(SLT)
*
      CALL mma_allocate(CMOX,2*NCMO,Label='CMOX')
      ipCMO=1+NCMO
* This is not the best solution, but I wanted to avoid having to rewrite
* the indexing code below just to use the CMO array directly
      call dcopy_(NCMO,CMO,1,CMOX,1)
      call dcopy_(NCMO,CMOX,1,CMOX(ipCMO),1)

*----------------------------------------------------------------------*
*     Compute Mulliken atomic charges of each active orbital           *
*             on each center to define the Active Site                 *
*----------------------------------------------------------------------*
      Call mma_allocate(Q,nUniqAt*(mAsh+1),Label='Q')
      ipQa=1+nUniqAt*mAsh
      Call Fzero(Q(ipQa),nUniqAt)
      Call mma_allocate(Z,nBmx*mAsh,Label='Z')
      lBas=0
      iOff=0
      Do iSym=1,nSym
         iSQ=1+iOff
         ipAsh=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         nBx=Max(1,nBas(iSym))
         Call DGEMM_('N','N',nBas(iSym),nAsh(iSym),nBas(iSym),
     &                      One,SQ(iSQ),nBx,
     &                          CMOX(ipAsh),nBx,
     &                      Zero,Z,nBx)
         jBas=lBas+1
         kBas=lBas+nBas(iSym)
         Call BasFun_Atom_Sym(nBas_per_Atom,nBas_Start,
     &                        Name,jBas,kBas,nUniqAt,.false.)
         Do ik=0,nAsh(iSym)-1
            nAk=nUniqAt*ik
            nBk=nBas(iSym)*ik
            jCMO=ipAsh+nBk-1
            jZ=nBk
            Do iAt=0,nUniqAt-1
               iBat=nBas_Start(1+iAt)
               jjCMO=jCMO+iBat
               jjZ=jZ+iBat
               nBat=nBas_per_Atom(1+iAt)
               iQ=1+nAk+iAt
               Q(iQ)=ddot_(nBat,CMOX(jjCMO),1,Z(jjZ),1)
            End Do
         End Do
         Do iAt=0,nUniqAt-1
            jQ=1+iAt
            iQa=ipQa+iAt
            Q(iQa) = Q(iQa)
     &             + ddot_(nAsh(iSym),Q(jQ),nUniqAt,Q(jQ),nUniqAt)
            If (sqrt(Q(iQa)).ge.Thrs) Then
               jBat=nBas_Start(1+iAt)+lBas
               NamAct(iAt+1)=Name(jBat)(1:LENIN)
            EndIf
         End Do
         lBas=lBas+nBas(iSym)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call mma_deallocate(Z)
      Call mma_deallocate(Q)

*     We have now completed the definition of the active site
*----------------------------------------------------------------------*
      Call mma_allocate(D_A,nUniqAt,Label='D_A')
      nActa=0
      Do iAt=1,nUniqAt
         If (NamAct(iAt).ne.blank) Then
            nActa=nActa+1
            D_A(nActa)=iAt
         EndIf
      End Do
      Do iAt=1,nActa
         jAt=D_A(iAt)
         NamAct(iAt)=NamAct(jAt)
      End Do
      Do iAt=nActa+1,nUniqAt
         NamAct(iAt)=blank
      End Do
      Write(u6,*)
      Write(u6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
      Write(u6,*)
      If (nActa.ne.0) Then
         Write(u6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
         Write(u6,*)
         Write(u6,*) (NamAct(i),i=1,nActa)
         Write(u6,*)
      ElseIf (.not.DoMP2 .and. .not.DoEnv) Then
         Write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
         Go To 2000
      Else
         Write(u6,'(A,18A4)') ' Selected atoms: *** None *** '
      EndIf

      Call mma_deallocate(D_A)
*----------------------------------------------------------------------*

      Call mma_allocate(OrbE,4*nOrb,Label='OrbE')
      ipOrbE=1
      Call Get_darray('RASSCF OrbE',OrbE(ipOrbE),nOrb)
      Call Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                    CMOX(ipCMO),OrbE(ipOrbE),TrX)
*
*---  MP2 calculation on the whole system (incompatible with DoMP2)
      If (DoEnv) Then
         Call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                           CMOX(ipCMO:),OrbE(ipOrbE:),E2_ab)
      EndIf
*----------------------------------------------------------------------*
*     Localize the inactive and virtual orbitals                       *
*                                                                      *
*        1) inactive orbitals ---> cholesky orbitals (orthonormal)     *
*        2) virtual orbitals ---> lin. indep. PAOs (non-orthonormal)   *
*                                                                      *
*----------------------------------------------------------------------*
      Thrd=1.d-06
      Call mma_allocate(D_vir,nBasT,Label='D_Vir')
      Call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nIsh,
     &                    nAsh,nSsh,CMOX(ipCMO),SQ,
     &                    D_vir)

      If(irc.ne.0) then
       write(u6,*) 'Localization failed in LovCASPT2'
       Call Abend()
      Endif

      ipEorb=ipOrbE+nOrb
      kEOcc=ipEorb+nOrb
      kEVir=kEOcc+nOrb
      Call mma_allocate(Xmo,2*NCMO,Label='XMO')
      iCMO=1+NCMO
      Call mma_allocate(Saa,nOrb,Label='Saa')
      Saa(:)=One


*     Inactive orbital selection                                       *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      lOff=0
      mOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*nFro(iSym)
         call dcopy_(nBas(iSym)*nIsh(iSym),CMOX(ipCMO+jOff),1,
     &                                    XMO(1+kOff),1)
         call dcopy_(nBas(iSym)*nIsh(iSym),CMOX(1+jOff),1,
     &                                    XMO(iCMO+kOff),1)
         jOff=lOff+nFro(iSym)
         call dcopy_(nIsh(iSym),OrbE(ipOrbE+jOff),1,
     &                         OrbE(ipEorb+mOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nIsh(iSym)
         lOff=lOff+nBas(iSym)
         mOff=mOff+nIsh(iSym)
      End Do
      ortho=.true.
*
      Call get_Orb_select(irc,XMO(iCMO),XMO,OrbE(ipEorb),
     &                        SQ,Saa,Name,NamAct,
     &                        nSym,nActa,nIsh,nBas,ortho,Thrs,ns_O)
      If(irc.ne.0) Return
      iOff=0
      kOff=0
      Do iSym=1,nSym
         lOff=iOff+nBas(iSym)*nFro(iSym)
         Do ik=nIsh(iSym),1,-1
            jOff=kOff+nBas(iSym)*(ik-1)
            call dcopy_(nBas(iSym),XMO(iCMO+jOff),1,
     &                            CMOX(1+lOff),1)
            lOff=lOff+nBas(iSym)
         End Do
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nIsh(iSym)
      End Do
      iloc=0
      loff=0
      Do iSym=1,nSym
         Do ik=nIsh(iSym),ns_O(iSym)+1,-1
            ie=ipEorb+loff+ik-1
            OrbE(kEOcc+iloc)=OrbE(ie)
            iloc=iloc+1
         End Do
         loff=loff+nIsh(iSym)
      End Do
      joff=0
      loff=0
      Do iSym=1,nSym
         koff=joff+nFro(iSym)+nIsh(iSym)-ns_O(iSym)
         Do ik=0,ns_O(iSym)-1
            ie=ipEorb+loff+ik
            OrbE(ipOrbE+koff+ik)=OrbE(ie)
         End Do
         loff=loff+nIsh(iSym)
         joff=joff+nBas(iSym)
      End Do

      If (all_Vir) Then

        Do iSym=1,nSym
         ns_V(iSym)=nSsh(iSym)
        End Do

      Else

*     Virtual orbital selection                                        *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      lOff=0
      mOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         call dcopy_(nBas(iSym)*nSsh(iSym),CMOX(ipCMO+jOff),1,
     &                                    XMO(1+kOff),1)
         call dcopy_(nBas(iSym)*nSsh(iSym),CMOX(1+jOff),1,
     &                                    XMO(iCMO+kOff),1)
         jOff=lOff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         call dcopy_(nSsh(iSym),OrbE(ipOrbE+jOff),1,
     &                         OrbE(ipEorb+mOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nSsh(iSym)
         lOff=lOff+nBas(iSym)
         mOff=mOff+nSsh(iSym)
      End Do
      ortho=.false.
      Call get_Saa(nSym,nBas,nSsh,SQ,XMO,Saa)
*
      Call get_Vir_select(irc,XMO(iCMO),XMO,OrbE(ipEorb),
     &                        SQ,Name,NamAct,D_vir,
     &                        nSym,nActa,nSsh,nBas,ortho,ns_V)
      If(irc.ne.0) Return
      Call mma_deallocate(D_vir)
      iOff=0
      kOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         call dcopy_(nBas(iSym)*nSsh(iSym),XMO(iCMO+kOff),1,
     &                                    CMOX(1+jOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nSsh(iSym)
      End Do
      iloc=0
      loff=0
      Do iSym=1,nSym
         Do ik=ns_V(iSym)+1,nSsh(iSym)
            ie=ipEorb+loff+ik-1
            OrbE(kEVir+iloc)=OrbE(ie)
            iloc=iloc+1
         End Do
         loff=loff+nSsh(iSym)
      End Do
      joff=0
      loff=0
      Do iSym=1,nSym
         koff=joff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         Do ik=0,ns_V(iSym)-1
            ie=ipEorb+loff+ik
            OrbE(ipOrbE+koff+ik)=OrbE(ie)
         End Do
         joff=joff+nBas(iSym)
         loff=loff+nSsh(iSym)
      End Do

      EndIf


*     MP2 calculation on the Frozen region                             *
*----------------------------------------------------------------------*
      If (DoMP2) Then

         iDo=0
         jDo=0
         nVV=0
         nOA=0
         Do iSym=1,nSym  ! setup info
            lnOrb(iSym)=nBas(iSym)
            lnOcc(iSym)=nIsh(iSym)-ns_O(iSym)
            lnFro(iSym)=nFro(iSym)+ns_O(iSym)
            lnDel(iSym)=nDel(iSym)+ns_V(iSym)
            lnVir(iSym)=nSsh(iSym)-ns_V(iSym)
            iDo=Max(iDo,lnOcc(iSym))
            jDo=Max(jDo,lnVir(iSym))
            nVV=nVV+lnVir(iSym)**2
            nOA=nOA+lnOcc(iSym)
         End Do
         If (Min(iDo,jDo).eq.0) goto 1000
*
         Call mma_allocate(Dmat,nVV+nOA,Label='DMat')
         ip_X=1
         ip_Y=ip_X+nVV
         DMat(:)=0.0D0
         Call FZero(XMO(iCMO),NCMO)
         iOff=0
         Do iSym=1,nSym
            kfr=1   +iOff+nBas(iSym)*nFro(iSym)
            kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
            call dcopy_(nBas(iSym)*lnOcc(iSym),CMOX(kfr),1,XMO(kto),1)
            kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym)
     &                                                     +ns_V(iSym))
            kto=kto+nBas(iSym)*lnOcc(iSym)
            call dcopy_(nBas(iSym)*lnVir(iSym),CMOX(kfr),1,XMO(kto),1)
            iOff=iOff+nBas(iSym)**2
         End Do
         Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
         If (iSkip.gt.0) Then
            Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            .true.)
            Call ChoMP2_Drv(irc,Dumm,XMO(iCMO),OrbE(kEOcc),OrbE(kEVir),
     &                      DMAT(ip_X),DMAT(ip_Y))
            Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            .false.)
            Call ChoMP2_Drv(irc,EMP2,XMO(iCMO),OrbE(kEOcc),OrbE(kEVir),
     &                      DMAT(ip_X),DMAT(ip_Y))
            If(irc.ne.0) then
              write(u6,*) 'Frozen region MP2 failed'
              Call Abend
            Endif
            iV=ip_X
            Do iSym=1,nSym
             TrF(iSym)=ddot_(lnVir(iSym),DMAT(iV),1+lnVir(iSym),[One],0)
              iV=iV+lnVir(iSym)**2
            End Do
         EndIf
         Call mma_deallocate(Dmat)
1000     Write(u6,*)

         If (nActa.eq.0) Then
            Write(u6,'(A,F18.10)')' Frozen region MP2 correction: ',EMP2
            Write(u6,*)
         EndIf
      EndIf
*                                                                      *
*----------------------------------------------------------------------*

      Call mma_deallocate(Saa)
      Call mma_deallocate(XMO)
*
*     Update the nFro, nIsh, nSsh, nDel for the Active site CASPT2
      Do iSym=1,nSym
         nFro(iSym)=nFro(iSym)+nIsh(iSym)-ns_O(iSym)
         nIsh(iSym)=ns_O(iSym)
         nDel(iSym)=nDel(iSym)+nSsh(iSym)-ns_V(iSym)
         nSsh(iSym)=ns_V(iSym)
         iDo=Max(iDo,nIsh(iSym))
         jDo=Max(jDo,nSsh(iSym))
      End Do

      Call Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                    CMOX,OrbE(ipOrbE),TrA)

      write(u6,*)
     &    '------------------------------------------------------'
      write(u6,*)
     &    ' Symm.  Tr(D):  Active        Frozen        Full      '
      write(u6,*)
     &    '------------------------------------------------------'
      STrA=Zero
      STrF=Zero
      STrX=Zero
      Do iSym=1,nSym
        If (DoEnv) TrF(iSym)=TrX(iSym) ! just a convention
        write(u6,'(2X,I4,10X,G11.4,3X,G11.4,3X,G11.4)') iSym,TrA(iSym),
     &       TrF(iSym),TrX(iSym)
        STrA=STrA+TrA(iSym)
        STrF=STrF+TrF(iSym)
        STrX=STrX+TrX(iSym)
      End Do
      write(u6,*)
     &    '------------------------------------------------------'
      write(u6,'(A,G11.4,3X,G11.4,3X,G11.4)')'          Sum:  ',
     & STrA,STrF,STrX
      write(u6,*)
     &    '------------------------------------------------------'
      write(u6,*)
*
      If (DoEnv) Then
         Call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                           CMOX,OrbE(ipOrbE),E2_Aonly)
         EMP2 = E2_ab - E2_Aonly
c         Write(6,'(A,F18.10)')' MP2 correction (environment): ',EMP2
c         Write(6,*)
      EndIf

      Call mma_deallocate(OrbE)

2000  Continue
      If (Min(iDo,jDo).eq.0) Then
         Write(u6,*)
         Write(u6,*)
     &       ' None of the inactive or virtual orbitals has been '
         Write(u6,*)
     &       ' assigned to the Active region of the molecule.    '
         Write(u6,*)
     &       ' This is presumably NOT what you want !!!          '
         Write(u6,*)' CASPT2 will Stop here. Bye Bye !! '
         Write(u6,*)
         Call Abend
      EndIf
*
      IF (IFQCAN.NE.0) IFQCAN=0 ! MOs to be recanonicalized on exit
      call dcopy_(NCMO,CMOX,1,CMO,1)

      Call mma_deallocate(CMOX)
      Call mma_deallocate(SQ)
      Call mma_deallocate(nBas_per_Atom)
      Call mma_deallocate(nBas_Start)

      End SUBROUTINE Lov_CASPT2
************************************************************************
*                                                                      *
************************************************************************
      Subroutine get_Saa(nSym,nBas,nOrb,Smn,Xmo,Saa)
      use definitions, only: iwp, wp
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      integer(kind=iwp), intent(in)::  nSym, nBas(nSym), nOrb(nSym)
      real(kind=wp), intent(in)::  Smn(*), Xmo(*)
      real(kind=wp), intent(out)::  Saa(*)

      real(kind=wp), Allocatable:: Z(:)
      integer(kind=iwp) mOb, iSym, iX, kX, lX, nBX, j, jK, jX, jZ, lk
      real(kind=wp), external:: DDot_
*
*
      mOb=nBas(1)*nOrb(1)
      Do iSym=2,nSym
         mOb=Max(mOb,nBas(iSym)*nOrb(iSym))
      End Do
      Call mma_allocate(Z,mOb,Label='Z')

      iX=1
      kX=1
      lX=1
      Do iSym=1,nSym
         nBx=Max(1,nBas(iSym))
         Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),
     &                      1.0d0,Smn(iX),nBx,
     &                            Xmo(kX),nBx,
     &                      0.0d0,Z,nBx)
         Do j=0,nOrb(iSym)-1
            jK=nBas(iSym)*j
            lk=kX+jK
            jZ=1+jK
            jX=lX+j
            Saa(jX)=ddot_(nBas(iSym),Xmo(lk),1,Z(jZ),1)
         End Do
         iX=iX+nBas(iSym)**2
         kX=kX+nBas(iSym)*nOrb(iSym)
         lX=lX+nOrb(iSym)
      End Do
      Call mma_deallocate(Z)
*
      End Subroutine get_Saa
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine LovCASPT2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            isFNO)
C
C     Purpose: put info in MP2 common blocks.
C
      use definitions, only: iwp
      use constants, only: Zero
      Use ChoMP2, only: C_os, ChkDecoMP2, ChoAlg, Decom_Def, DecoMP2,
     &                  DoFNO, EOSMP2, ForceBatch, l_Dii, MxQual_Def,
     &                  MxQualMP2, OED_Thr, set_cd_thr, shf, SOS_mp2,
     &                  Span_Def, SpanMP2, ThrMP2, Verbose
      use cOrbInf, only: nSym, nOrb, nOcc, nFro, nDel, nExt
      Implicit none
      Integer(kind=iwp), intent(in) :: mSym
      Integer(kind=iwp), intent(in) :: lnOrb(8), lnOcc(8), lnFro(8),
     &                                 lnDel(8), lnVir(8)
      Logical(kind=iwp), intent(in) ::  isFNO

      Integer(kind=iwp) iSym
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
      shf=Zero
C
      DoFNO=isFNO
      l_Dii=nOcc(1)
      Do iSym=2,nSym
         l_Dii=l_Dii+nOcc(iSym)
      End Do
C
      End SubRoutine LovCASPT2_putInf

      Subroutine Energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                         CMO,OrbE,E2_ab)

      use definitions, only: iwp, wp, u6
      use constants, only: Zero
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      integer(kind=iwp), intent(in):: nSym
      integer(kind=iwp), intent(in):: nBas(nSym), nFro(nSym), nIsh(nSym)
      integer(kind=iwp), intent(in):: nAsh(nSym), nSsh(nSym), nDel(nSym)
      real(kind=wp), intent(in)::  CMO(*), OrbE(*)
      real(kind=wp), intent(out)::  E2_ab

      integer(kind=iwp) nAct(8), lnOrb(8), lnOcc(8), lnFro(8), lnDel(8),
     &                  lnVir(8)
      real(kind=wp) Dummy(1)
      real(kind=wp), Allocatable:: Eorb(:), CMOX(:)
      integer(kind=iwp) iE, ifr, ioff, irc, iSkip, iSym, ito, joff, k,
     &                  kEOcc, kEVir, kfr, koff, kto, nBB, nOA, nOrb,
     &                  nVV
*
*
      Call Izero(nAct,nSym)
      nVV=0
      nOrb=0
      Do iSym=1,nSym
         iE=1+nOrb+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (OrbE(iE+k).lt.Zero) nAct(iSym)=nAct(iSym)+1
         End Do
         nVV=nVV+nSsh(iSym)**2
         nOrb=nOrb+nBas(iSym)
      End Do
*
      nBB=0
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnOrb(iSym)=nBas(iSym)
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)+nAct(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnDel(iSym)=nDel(iSym)
         nBB=nBB+nBas(iSym)**2
         nOA=nOA+lnOcc(iSym)
      End Do
*
      Call mma_allocate(Eorb,2*nOrb,Label='Eorb')
      kEOcc=1
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(lnOcc(iSym),OrbE(ifr),1,Eorb(ito),1)
         ifr=1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),OrbE(ifr),1,Eorb(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do
*
      Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.false.)
      Call mma_allocate(CMOX,nBB,Label='CMOX')
      CMOX(:)=Zero
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=1+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,CMOX(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,CMOX(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
*
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,E2_ab,CMOX,Eorb(kEOcc),Eorb(kEVir),
     &                   Dummy,Dummy)
         If(irc.ne.0) then
           write(u6,*) 'MP2 calculation failed in energy_AplusB !'
           Call Abend
         Endif
      Else
         write(u6,*)
         write(u6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(u6,*)'combinations of inactive and virtual orbitals !! '
         write(u6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend()
      Endif
      Call mma_deallocate(CMOX)
*
      Call mma_deallocate(Eorb)
*
      End Subroutine Energy_AplusB
