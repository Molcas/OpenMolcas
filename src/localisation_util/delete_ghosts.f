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
* Copyright (C) 2010, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE Delete_Ghosts(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                         NAME,nUniqAt,ThrS,isCASPT2,CMO,EOrb)
************************************************************************
*                                                                      *
* Purpose:  Eliminates MOs of ghost atoms from PT2 treatment           *
*                                                                      *
* Author:   F. Aquilante  (Geneva, July 2010)                          *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
*
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nAsh(nSym),nSsh(nSym),
     &        nDel(nSym)
      Integer irc,nUniqAt
      Real*8  ThrS, CMO(*), EOrb(*)
      Logical isCASPT2
      Character*(LENIN4) NAME(*)
      Character*(LENIN) blank, NamAct(mxAtom), tmp
      Integer n_OK(8)
************************************************************************
      jD(i) = iWork(ip_iD-1+i)
************************************************************************
*
*
      irc=0
      blank='   '
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
*----------------------------------------------------------------------*
*
      nBasT=0
      ntri=0
      nSQ=0
      nBmx=0
      nSmx=0
      mAsh=0
      Do i=1,nSym
        nBasT=nBasT+nBas(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nBmx=Max(nBmx,nBas(i))
        nSmx=Max(nSmx,nSsh(i))
        mAsh=Max(mAsh,nIsh(i)+nAsh(i))
      End Do
      NCMO=nSQ
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
         ltri=ltri+nBas(iSym)*(nBas(iSym)+1)/2
         lsq=lsq+nBas(iSym)**2
      End Do
      CALL GetMem('SLT','FREE','REAL',ipS,nTri)
*
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      call dcopy_(NCMO,CMO,1,WORK(LCMO),1)

*----------------------------------------------------------------------*
*     Compute Mulliken atomic charges of each occupied orbital         *
*             on each center to define the Active Site                 *
*----------------------------------------------------------------------*
      Call GetMem('Qai','Allo','Real',ipQ,nUniqAt*(mAsh+1))
      ipQa=ipQ+nUniqAt*mAsh
      Call Fzero(Work(ipQa),nUniqAt)
      Call GetMem('Zm','Allo','Real',ipZ,nBmx*mAsh)
      lBas=0
      iOff=0
      Do iSym=1,nSym
         nOkk=nIsh(iSym)+nAsh(iSym)
         iSQ=ipSQ+iOff
         ipAsh=LCMO+iOff+nBas(iSym)*nFro(iSym)
         nBx=Max(1,nBas(iSym))
         Call DGEMM_('N','N',nBas(iSym),nOkk,nBas(iSym),
     &                      1.0d0,Work(iSQ),nBx,
     &                            Work(ipAsh),nBx,
     &                      0.0d0,Work(ipZ),nBx)
         jBas=lBas+1
         kBas=lBas+nBas(iSym)
         Call BasFun_Atom_(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                     Name,jBas,kBas,nUniqAt,.false.)
         Do ik=0,nOkk-1
            nAk=nUniqAt*ik
            nBk=nBas(iSym)*ik
            jCMO=ipAsh+nBk-1
            jZ=ipZ+nBk-1
            Do iAt=0,nUniqAt-1
               iBat=iWork(ip_nBas_Start+iAt)
               jjCMO=jCMO+iBat
               jjZ=jZ+iBat
               nBat=iWork(ip_nBas_per_Atom+iAt)
               iQQ=ipQ+nAk+iAt
               Work(iQQ)=ddot_(nBat,Work(jjCMO),1,Work(jjZ),1)
            End Do
         End Do
         Do iAt=0,nUniqAt-1
            jQ=ipQ+iAt
            iQa=ipQa+iAt
            Work(iQa) = Work(iQa)
     &                + ddot_(nOkk,Work(jQ),nUniqAt,
     &                            Work(jQ),nUniqAt)
            If (sqrt(Work(iQa)).ge.ThrS) Then
               jBat=iWork(ip_nBas_Start+iAt)+lBas
               If (iWork(ip_nBas_per_Atom+iAt) .gt. 0)
     &            NamAct(iAt+1)=Name(jBat)(1:LENIN)
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
         If (NamAct(iAt)(1:4).ne.blank) Then
            iWork(iD+nActa)=iAt
            nActa=nActa+1
         EndIf
      End Do
      Do iAt=1,nActa
         jAt=iWork(iD+iAt-1)
         NamAct(iAt)=NamAct(jAt)
      End Do
      Do iAt=nActa+1,nUniqAt
         NamAct(iAt)(1:4)=Trim(blank)
      End Do
      Write(6,*)
      Write(6,'(A,F6.3)') ' Threshold for atom selection: ',ThrS
      Write(6,*)
      If (nActa.ne.0) Then
         Write(6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
         Write(6,*)
         Write(6,*) (NamAct(i),i=1,nActa)
         Write(6,*)
      Else
         Write(6,*)' None of the occupied non-frozen orbitals has been '
         Write(6,*)' assigned to the Active region of the molecule.    '
         Write(6,*)' This is presumably NOT what you want !!!          '
         Write(6,*)' I will Stop here. Bye Bye !! '
         Write(6,*)
         Call Abend
      EndIf

      Call GetMem('ID_A','Free','Inte',iD,nUniqAt)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)

*----------------------------------------------------------------------*
*     Virtual orbital selection                                        *
*----------------------------------------------------------------------*
      Call GetMem('ID_vir','Allo','Inte',ip_iD,nBmx+2*nSmx)
      nSmall=nBmx**2+nSmx+3*nBmx*nSmx
      Call GetMem('Small','Allo','Real',iS,nSmall)
      iQ=iS+nBmx**2
      iC=iQ+nSmx
      iZ=iC+nBmx*nSmx
      iX=iZ+nBmx*nSmx
*
      iOff=0
      kOff=0
      Do iSym=1,nSym

         nBa=0
         Do ia=1,nBas(iSym)
            ja=ia+iOff
            tmp=Name(ja)(1:LENIN)
            Do j=1,nActa
               If(NamAct(j).eq.tmp) Then
                  iWork(ip_iD+nBa)=ia
                  nBa=nBa+1
               EndIf
            End Do
         End Do

         iCMO=LCMO+kOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         Do ia=1,nBa
            ifr=iCMO+jD(ia)-1
            ito=iC+ia-1
            call dcopy_(nSsh(iSym),Work(ifr),nBas(iSym),Work(ito),nBa)
         End Do

         iSQ=ipSQ+kOff
         Do ia=1,nBa
            jb=jD(ia)
            jfr=iSQ+nBas(iSym)*(jb-1)
            jto=iS+nBas(iSym)*(ia-1)
            call dcopy_(nBas(iSym),Work(jfr),1,Work(jto),1)
         End Do

         nBx=Max(1,nBas(iSym))
         nBax=Max(1,nBa)
         Call DGEMM_('T','N',nBa,nSsh(iSym),nBas(iSym),
     &                    1.0d0,Work(iS),nBx,
     &                          Work(iCMO),nBx,
     &                    0.0d0,Work(iZ),nBax)
         Do i=0,nSsh(iSym)-1
            jQ=iQ+i
            jC=iC+nBa*i
            jZ=iZ+nBa*i
            Work(jQ)=ddot_(nBa,Work(jC),1,Work(jZ),1)**2
         End Do
         n_OK(iSym)=0
         n_KO=0
         Do i=1,nSsh(iSym)
            jQ=iQ+i-1
            jfr=iCMO+nBas(iSym)*(i-1)
            If (sqrt(Work(jQ)).ge.ThrS) Then
               jX=iX+nBas(iSym)*n_OK(iSym)
               call dcopy_(nBas(iSym),Work(jfr),1,Work(jX),1)
               iWork(ip_iD+nBmx+n_OK(iSym))=i
               n_OK(iSym)=n_OK(iSym)+1
            Else
               jZ=iZ+nBas(iSym)*n_KO
               call dcopy_(nBas(iSym),Work(jfr),1,Work(jZ),1)
               iWork(ip_iD+nBmx+nSmx+n_KO)=i
               n_KO=n_KO+1
            EndIf
         End Do
*
         call dcopy_(nBas(iSym)*n_OK(iSym),Work(iX),1,Work(iCMO),1)
         kCMO=iCMO+nBas(iSym)*n_OK(iSym)
         call dcopy_(nBas(iSym)*n_KO,Work(iZ),1,Work(kCMO),1)
         If (.not.isCASPT2) Then
            jZ=iZ
            jOff=iOff+nFro(iSym)+nOkk
            Do i=nBmx+1,nBmx+n_OK(iSym)
               iv=jD(i)+jOff
               Work(jZ)=EOrb(iv)
               jZ=jZ+1
            End Do
            Do i=nBmx+nSmx+1,nBmx+nSmx+n_KO
               iv=jD(i)+jOff
               Work(jZ)=EOrb(iv)
            End Do
            call dcopy_(nSsh(iSym),Work(iZ),1,EOrb(1+jOff),1)
         EndIf
*
         iOff=iOff+nBas(iSym)
         kOff=kOff+nBas(iSym)**2
      End Do
      Call GetMem('Small','Free','Real',iS,nSmall)
      Call GetMem('ID_vir','Free','Inte',ip_iD,nBmx+2*nSmx)
*                                                                      *
*----------------------------------------------------------------------*
*
*     Update nSsh, nDel for the Active site PT2
      Do iSym=1,nSym
         nDel(iSym)=nDel(iSym)+nSsh(iSym)-n_OK(iSym)
         nSsh(iSym)=n_OK(iSym)
      End Do
*
      call dcopy_(NCMO,WORK(LCMO),1,CMO,1)
*
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
      CALL GetMem('SMAT','FREE','REAL',ipSQ,nSQ)
*
      Return
      End
