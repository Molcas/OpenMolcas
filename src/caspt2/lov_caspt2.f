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
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
      COMMON /LovCAS3/ STrA, STrF, STrX
*
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nAsh(nSym),nSsh(nSym),
     &        nDel(nSym)
      Integer irc,nUniqAt,IFQCAN
      Real*8  Thrs, EMP2
      Logical DoMP2, DoEnv, all_Vir
      Character(LENIN8) NAME(*)
      Character(LENIN) blank, NamAct(mxAtom)
      Logical ortho
      Real*8  TrA(8), TrF(8), TrX(8)
      Integer ns_O(8), ns_V(8)
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Real*8 CMO(*)
*
*
      irc=0
      EMP2=Zero
      blank='   '
      iDo=0
      jDo=0
      If (DoEnv .and. DoMP2) Then
         Call WarningMessage(1,'Both DoEnv and DoMP2 selected.')
         Write (6,'(/,A)') ' DoMP2 will be ignored.'
         DoMP2=.false.
      EndIf
      If (all_Vir .and. DoMP2) Then
         Call WarningMessage(1,'Both VirAll and DoMP2 selected.')
         Write (6,'(/,A)') ' DoMP2 will be ignored.'
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
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
*     nUniqAt = # of symm. unique atoms. Initialize NamAct to blanks.
*     ---------------------------------------------------------------

      If (nUniqAt.lt.1 .or. nUniqAt.gt.MxAtom) Then
         Write(6,'(A,I9)') 'nUniqAt =',nUniqAt
         Call Abend()
      End If
      Do iAt=1,nUniqAt
         NamAct(iAt)=blank
      End Do

C     Allocate and get index arrays for basis functions per atom.
C     -----------------------------------------------------------

      l_nBas_per_Atom = nUniqAt
      l_nBas_Start    = nUniqAt
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
*
*----------------------------------------------------------------------*
*     Read the overlap matrix                                          *
*----------------------------------------------------------------------*
      CALL GetMem('SMAT','ALLO','REAL',ipSQ,nSQ)
      CALL GetMem('SLT','ALLO','REAL',ipS,nTri)
      isymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
      If(irc.ne.0) return
      ltri=0
      lsq=0
      Do iSym=1,nSym
         Call Square(Work(ipS+ltri),Work(ipSQ+lsq),1,nBas(iSym),
     &                                               nBas(iSym))
         ltri=ltri+nBas(iSym)*(nBAs(iSym)+1)/2
         lsq=lsq+nBas(iSym)**2
      End Do
      CALL GetMem('SLT','FREE','REAL',ipS,nTri)
*
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,2*NCMO)
      ipCMO=LCMO+NCMO
* This is not the best solution, but I wanted to avoid having to rewrite
* the indexing code below just to use the CMO array directly
      call dcopy_(NCMO,CMO,1,WORK(LCMO),1)
      call dcopy_(NCMO,WORK(LCMO),1,WORK(ipCMO),1)

*----------------------------------------------------------------------*
*     Compute Mulliken atomic charges of each active orbital           *
*             on each center to define the Active Site                 *
*----------------------------------------------------------------------*
      Call GetMem('Qai','Allo','Real',ipQ,nUniqAt*(mAsh+1))
      ipQa=ipQ+nUniqAt*mAsh
      Call Fzero(Work(ipQa),nUniqAt)
      Call GetMem('Zm','Allo','Real',ipZ,nBmx*mAsh)
      lBas=0
      iOff=0
      Do iSym=1,nSym
         iSQ=ipSQ+iOff
         ipAsh=LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         nBx=Max(1,nBas(iSym))
         Call DGEMM_('N','N',nBas(iSym),nAsh(iSym),nBas(iSym),
     &                      One,Work(iSQ),nBx,
     &                          Work(ipAsh),nBx,
     &                      Zero,Work(ipZ),nBx)
         jBas=lBas+1
         kBas=lBas+nBas(iSym)
         Call BasFun_Atom_(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                     Name,jBas,kBas,nUniqAt,.false.)
         Do ik=0,nAsh(iSym)-1
            nAk=nUniqAt*ik
            nBk=nBas(iSym)*ik
            jCMO=ipAsh+nBk-1
            jZ=ipZ+nBk-1
            Do iAt=0,nUniqAt-1
               iBat=iWork(ip_nBas_Start+iAt)
               jjCMO=jCMO+iBat
               jjZ=jZ+iBat
               nBat=iWork(ip_nBas_per_Atom+iAt)
               iQ=ipQ+nAk+iAt
               Work(iQ)=ddot_(nBat,Work(jjCMO),1,Work(jjZ),1)
            End Do
         End Do
         Do iAt=0,nUniqAt-1
            jQ=ipQ+iAt
            iQa=ipQa+iAt
            Work(iQa) = Work(iQa)
     &                + ddot_(nAsh(iSym),Work(jQ),nUniqAt,
     &                                  Work(jQ),nUniqAt)
            If (sqrt(Work(iQa)).ge.Thrs) Then
               jBat=iWork(ip_nBas_Start+iAt)+lBas
               NamAct(iAt+1)=Name(jBat)(1:LENIN)
            EndIf
         End Do
         lBas=lBas+nBas(iSym)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call GetMem('Zm','Free','Real',ipZ,nBmx*mAsh)
      Call GetMem('Qai','Free','Real',ipQ,nUniqAt*(mAsh+1))

*     We have now completed the definition of the active site
*----------------------------------------------------------------------*
      Call GetMem('ID_A','Allo','Inte',iD,nUniqAt)
      nActa=0
      Do iAt=1,nUniqAt
         If (NamAct(iAt).ne.blank) Then
            iWork(iD+nActa)=iAt
            nActa=nActa+1
         EndIf
      End Do
      Do iAt=1,nActa
         jAt=iWork(iD+iAt-1)
         NamAct(iAt)=NamAct(jAt)
      End Do
      Do iAt=nActa+1,nUniqAt
         NamAct(iAt)=blank
      End Do
      Write(6,*)
      Write(6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
      Write(6,*)
      If (nActa.ne.0) Then
         Write(6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
         Write(6,*)
         Write(6,*) (NamAct(i),i=1,nActa)
         Write(6,*)
      ElseIf (.not.DoMP2 .and. .not.DoEnv) Then
         Write(6,'(A,18A4)') ' Selected atoms: *** None *** '
         Go To 2000
      Else
         Write(6,'(A,18A4)') ' Selected atoms: *** None *** '
      EndIf

      Call GetMem('ID_A','Free','Inte',iD,nUniqAt)
*----------------------------------------------------------------------*

      Call GetMem('Eorb','Allo','Real',ipOrbE,4*nOrb)
      Call Get_darray('RASSCF OrbE',Work(ipOrbE),nOrb)
      Call Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                    Work(ipCMO),Work(ipOrbE),TrX)
*
*---  MP2 calculation on the whole system (incompatible with DoMP2)
      If (DoEnv) Then
         Call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                           Work(ipCMO),Work(ipOrbE),E2_ab)
      EndIf
*----------------------------------------------------------------------*
*     Localize the inactive and virtual orbitals                       *
*                                                                      *
*        1) inactive orbitals ---> cholesky orbitals (orthonormal)     *
*        2) virtual orbitals ---> lin. indep. PAOs (non-orthonormal)   *
*                                                                      *
*----------------------------------------------------------------------*
      Thrd=1.d-06
      Call GetMem('ID_vir','Allo','Inte',iD_vir,nBasT)
      Call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nIsh,
     &                    nAsh,nSsh,Work(ipCMO),Work(ipSQ),
     &                    iWork(iD_vir))

      If(irc.ne.0) then
       write(6,*) 'Localization failed in LovCASPT2'
       Call Abend
      Endif

      ipEorb=ipOrbE+nOrb
      kEOcc=ipEorb+nOrb
      kEVir=kEOcc+nOrb
      Call GetMem('XMO','Allo','Real',ipXmo,2*NCMO)
      iCMO=ipXmo+NCMO
      Call GetMem('Saa','Allo','Real',ipSaa,nOrb)
      call dcopy_(nOrb,[One],0,Work(ipSaa),1)


*     Inactive orbital selection                                       *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      lOff=0
      mOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*nFro(iSym)
         call dcopy_(nBas(iSym)*nIsh(iSym),Work(ipCMO+jOff),1,
     &                                    Work(ipXMO+kOff),1)
         call dcopy_(nBas(iSym)*nIsh(iSym),Work(LCMO+jOff),1,
     &                                    Work(iCMO+kOff),1)
         jOff=lOff+nFro(iSym)
         call dcopy_(nIsh(iSym),Work(ipOrbE+jOff),1,
     &                         Work(ipEorb+mOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nIsh(iSym)
         lOff=lOff+nBas(iSym)
         mOff=mOff+nIsh(iSym)
      End Do
      ortho=.true.
*
      Call get_Orb_select(irc,Work(iCMO),Work(ipXMO),Work(ipEorb),
     &                        Work(ipSQ),Work(ipSaa),Name,NamAct,
     &                        nSym,nActa,nIsh,nBas,ortho,Thrs,ns_O)
      If(irc.ne.0) Return
      iOff=0
      kOff=0
      Do iSym=1,nSym
         lOff=iOff+nBas(iSym)*nFro(iSym)
         Do ik=nIsh(iSym),1,-1
            jOff=kOff+nBas(iSym)*(ik-1)
            call dcopy_(nBas(iSym),Work(iCMO+jOff),1,
     &                            Work(LCMO+lOff),1)
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
            Work(kEOcc+iloc)=Work(ie)
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
            Work(ipOrbE+koff+ik)=Work(ie)
         End Do
         loff=loff+nIsh(iSym)
         joff=joff+nBas(iSym)
      End Do

      If (all_Vir) Then
        Do iSym=1,nSym
         ns_V(iSym)=nSsh(iSym)
        End Do
        goto 999
      EndIf

*     Virtual orbital selection                                        *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      lOff=0
      mOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         call dcopy_(nBas(iSym)*nSsh(iSym),Work(ipCMO+jOff),1,
     &                                    Work(ipXMO+kOff),1)
         call dcopy_(nBas(iSym)*nSsh(iSym),Work(LCMO+jOff),1,
     &                                    Work(iCMO+kOff),1)
         jOff=lOff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         call dcopy_(nSsh(iSym),Work(ipOrbE+jOff),1,
     &                         Work(ipEorb+mOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nSsh(iSym)
         lOff=lOff+nBas(iSym)
         mOff=mOff+nSsh(iSym)
      End Do
      ortho=.false.
      Call get_Saa(nSym,nBas,nSsh,Work(ipSQ),Work(ipXMO),Work(ipSaa))
*
      Call get_Vir_select(irc,Work(iCMO),Work(ipXMO),Work(ipEorb),
     &                        Work(ipSQ),Name,NamAct,iWork(iD_vir),
     &                        nSym,nActa,nSsh,nBas,ortho,ns_V)
      If(irc.ne.0) Return
      Call GetMem('ID_vir','Free','Inte',iD_vir,nBasT)
      iOff=0
      kOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         call dcopy_(nBas(iSym)*nSsh(iSym),Work(iCMO+kOff),1,
     &                                    Work(LCMO+jOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nSsh(iSym)
      End Do
      iloc=0
      loff=0
      Do iSym=1,nSym
         Do ik=ns_V(iSym)+1,nSsh(iSym)
            ie=ipEorb+loff+ik-1
            Work(kEVir+iloc)=Work(ie)
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
            Work(ipOrbE+koff+ik)=Work(ie)
         End Do
         joff=joff+nBas(iSym)
         loff=loff+nSsh(iSym)
      End Do

999   Continue

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
         Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
         ip_Y=ip_X+nVV
         Call FZero(Work(ip_X),nVV+nOA)
         Call FZero(Work(iCMO),NCMO)
         iOff=0
         Do iSym=1,nSym
            kfr=LCMO+iOff+nBas(iSym)*nFro(iSym)
            kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
            call dcopy_(nBas(iSym)*lnOcc(iSym),Work(kfr),1,Work(kto),1)
            kfr=LCMO+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym)
     &                                                     +ns_V(iSym))
            kto=kto+nBas(iSym)*lnOcc(iSym)
            call dcopy_(nBas(iSym)*lnVir(iSym),Work(kfr),1,Work(kto),1)
            iOff=iOff+nBas(iSym)**2
         End Do
         Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
         If (iSkip.gt.0) Then
            Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            ip_X,ip_Y,.true.)
            Call ChoMP2_Drv(irc,Dumm,Work(iCMO),Work(kEOcc),Work(kEVir))
            Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            ip_X,ip_Y,.false.)
            Call ChoMP2_Drv(irc,EMP2,Work(iCMO),Work(kEOcc),Work(kEVir))
            If(irc.ne.0) then
              write(6,*) 'Frozen region MP2 failed'
              Call Abend
            Endif
            iV=ip_X
            Do iSym=1,nSym
             TrF(iSym)=ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),[One],0)
              iV=iV+lnVir(iSym)**2
            End Do
         EndIf
         Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
1000     Write(6,*)

         If (nActa.eq.0) Then
            Write(6,'(A,F18.10)')' Frozen region MP2 correction: ',EMP2
            Write(6,*)
         EndIf
      EndIf
*                                                                      *
*----------------------------------------------------------------------*

      Call GetMem('Saa','Free','Real',ipSaa,nOrb)
      Call GetMem('XMO','Free','Real',ipXmo,2*NCMO)
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
     &                    Work(LCMO),Work(ipOrbE),TrA)

      write(6,*)'------------------------------------------------------'
      write(6,*)' Symm.  Tr(D):  Active        Frozen        Full      '
      write(6,*)'------------------------------------------------------'
      STrA=Zero
      STrF=Zero
      STrX=Zero
      Do iSym=1,nSym
        If (DoEnv) TrF(iSym)=TrX(iSym) ! just a convention
        write(6,'(2X,I4,10X,G11.4,3X,G11.4,3X,G11.4)') iSym,TrA(iSym),
     &       TrF(iSym),TrX(iSym)
        STrA=STrA+TrA(iSym)
        STrF=STrF+TrF(iSym)
        STrX=STrX+TrX(iSym)
      End Do
      write(6,*)'------------------------------------------------------'
      write(6,'(A,G11.4,3X,G11.4,3X,G11.4)')'          Sum:  ',
     & STrA,STrF,STrX
      write(6,*)'------------------------------------------------------'
      write(6,*)
*
      If (DoEnv) Then
         Call energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                           Work(LCMO),Work(ipOrbE),E2_Aonly)
         EMP2 = E2_ab - E2_Aonly
c         Write(6,'(A,F18.10)')' MP2 correction (environment): ',EMP2
c         Write(6,*)
      EndIf

      Call GetMem('Eorb','Free','Real',ipOrbE,4*nOrb)

2000  Continue
      If (Min(iDo,jDo).eq.0) Then
         Write(6,*)
         Write(6,*)' None of the inactive or virtual orbitals has been '
         Write(6,*)' assigned to the Active region of the molecule.    '
         Write(6,*)' This is presumably NOT what you want !!!          '
         Write(6,*)' CASPT2 will Stop here. Bye Bye !! '
         Write(6,*)
         Call Abend
      EndIf
*
      IF (IFQCAN.NE.0) IFQCAN=0 ! MOs to be recanonicalized on exit
      call dcopy_(NCMO,WORK(LCMO),1,CMO,1)

      CALL GETMEM('LCMO','FREE','REAL',LCMO,2*NCMO)
      CALL GetMem('SMAT','FREE','REAL',ipSQ,nSQ)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine get_Saa(nSym,nBas,nOrb,Smn,Xmo,Saa)
      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nOrb(nSym)
      Real*8  Smn(*), Xmo(*), Saa(*)
#include "WrkSpc.fh"
*
*
      mOb=nBas(1)*nOrb(1)
      Do iSym=2,nSym
         mOb=Max(mOb,nBas(iSym)*nOrb(iSym))
      End Do
      Call GetMem('Z','Allo','Real',ipZ,mOb)

      iX=1
      kX=1
      lX=1
      Do iSym=1,nSym
         nBx=Max(1,nBas(iSym))
         Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),
     &                      1.0d0,Smn(iX),nBx,
     &                            Xmo(kX),nBx,
     &                      0.0d0,Work(ipZ),nBx)
         Do j=0,nOrb(iSym)-1
            jK=nBas(iSym)*j
            lk=kX+jK
            jZ=ipZ+jK
            jX=lX+j
            Saa(jX)=ddot_(nBas(iSym),Xmo(lk),1,Work(jZ),1)
         End Do
         iX=iX+nBas(iSym)**2
         kX=kX+nBas(iSym)*nOrb(iSym)
         lX=lX+nOrb(iSym)
      End Do
      Call GetMem('Z','Free','Real',ipZ,mOb)
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine LovCASPT2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                            ip_X,ip_Y,isFNO)
C
C     Purpose: put info in MP2 common blocks.
C
#include "implicit.fh"
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Integer ip_X, ip_Y
      Logical isFNO
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
      DoFNO=isFNO
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

      Subroutine Energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                              CMO,OrbE,E2_ab)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nIsh(nSym)
      Integer nAsh(nSym), nSsh(nSym), nDel(nSym)
      Real*8  CMO(*), OrbE(*)
#include "WrkSpc.fh"
      Integer nAct(8), lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
*
*
      Call Izero(nAct,nSym)
      nVV=0
      nOrb=0
      Do iSym=1,nSym
         iE=1+nOrb+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (OrbE(iE+k).lt.0.0d0) nAct(iSym)=nAct(iSym)+1
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
      Call GetMem('EOV','Allo','Real',ipEorb,2*nOrb)
      kEOcc=ipEorb
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(lnOcc(iSym),OrbE(ifr),1,Work(ito),1)
         ifr=1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),OrbE(ifr),1,Work(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do
*
      Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,iDummy,
     &                           jDummy,.false.)
      Call GetMem('CMON','Allo','Real',iCMO,nBB)
      Call FZero(Work(iCMO),nBB)
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
*
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,E2_ab,Work(iCMO),Work(kEOcc),Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 calculation failed in energy_AplusB !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
      Call GetMem('CMON','Free','Real',iCMO,nBB)
*
      Call GetMem('EOV ','Free','Real',ipEorb,2*nOrb)
*
      Return
      End
