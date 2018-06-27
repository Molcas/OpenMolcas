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
      SubRoutine LovMP2_Drv(irc,EMP2,CMO,EOcc,EVir,NamAct,nActa,
     &                          Thrs,DoMP2,all_Vir)

************************************************************************
*                                                                      *
* Purpose:  setup of Localized occupied-virtual MP2  (LovMP2)          *
*           The MP2 correction to the energy will later be computed    *
*           only for the "active region" of the molecule.              *
*           The MP2 correction due to the remaining frozen region      *
*           is computed here if DoMP2=.true.                           *
*                                                                      *
* Author:   F. Aquilante  (Geneva, Jun. 2008)                          *
*                                                                      *
************************************************************************
#include "implicit.fh"
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "WrkSpc.fh"
      Real*8 OED_Thr, EOSMP2, C_os, XEMP2, Wref
      Common / ChSOS2 / OED_Thr, EOSMP2, C_os
      Common / ChMP24 / XEMP2, Wref

      Dimension CMO(*), EOcc(*), EVir(*)
      Real*8  Thrs
      Logical DoMP2, all_Vir
      Character*(LENIN8) Name(mxBas)
      Character*(LENIN) NamAct(*)
      Logical ortho
      Real*8  TrA(8), TrF(8), TrX(8)
      Integer ns_O(8), ns_V(8), nZero(8)
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Integer lnOcc2(8), lnFro2(8), lnDel2(8), lnVir2(8), nAuxO(8)

      Character*3  ThisNm
      Character*10 SecNam
      Parameter (SecNam = 'LovMP2_Drv', ThisNm = 'Dry')


      Call qEnter(ThisNm)
      irc=0
      EMP2=0.0d0
      EFRO=0.0d0
      EOSF=0.0d0
      iDo=0
      jDo=0
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
*----------------------------------------------------------------------*
*
      nxBasT=0
      ntri=0
      nSQ=0
      nBmx=0
      nxOrb=0
      Do i=1,nSym
        TrA(i)=0.0d0
        TrF(i)=0.0d0
        TrX(i)=0.0d0
        nZero(i)=0
        nxBasT=nxBasT+nBas(i)
        nxOrb=nxOrb+nFro(i)+nOcc(i)+nExt(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nBmx=Max(nBmx,nBas(i))
      End Do
      IF(nxBasT.GT.mxBas) then
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nxBasT)
*
*----------------------------------------------------------------------*
*     Read the overlap matrix                                          *
*----------------------------------------------------------------------*
      CALL GetMem('SMAT','ALLO','REAL',ipSQ,nSQ)
      CALL GetMem('SLT','ALLO','REAL',ipS,nTri)
      isymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
      If(irc.ne.0) Then
        Call qExit(ThisNm)
        return
      End If
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
      Write(6,'(A,F15.6)') ' Threshold for atom selection: ',Thrs
      Write(6,*)
      If (nActa.ne.0) Then
         Write(6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
         Write(6,*)
         Write(6,*) (NamAct(i),i=1,nActa)
         Write(6,*)
      ElseIf (.not.DoMP2) Then
         Write(6,'(A,18A4)') ' Selected atoms: *** None *** '
         Go To 2000
      Else
         Write(6,'(A,18A4)') ' Selected atoms: *** None *** '
      EndIf

*----------------------------------------------------------------------*
      Call Get_Tr_Dab(nSym,nBas,nFro,nOcc,nExt,nDel,
     &                CMO,EOcc,EVir,TrX)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Localize the inactive and virtual orbitals                       *
*                                                                      *
*        1) inactive orbitals ---> cholesky orbitals (orthonormal)     *
*        2) virtual orbitals ---> lin. indep. PAOs (non-orthonormal)   *
*                                                                      *
*----------------------------------------------------------------------*
      CALL GetMem('LCMO','ALLO','REAL',iCMO,3*nSQ)
      ipXMO=iCMO+nSQ
      iXMO=ipXMO+nSQ
      call dcopy_(nSQ,CMO,1,Work(ipXMO),1)
      Thrd=1.0d-06
      Call GetMem('ID_vir','Allo','Inte',iD_vir,nxBasT)
      Call Cho_ov_Loc(irc,Thrd,nSym,nBas,nFro,nOcc,
     &                    nZero,nExt,Work(ipXMO),Work(ipSQ),
     &                    iWork(iD_vir))

      If(irc.ne.0) then
       write(6,*) 'Localization failed in LovMP2'
       Call Abend
      Endif

      Call GetMem('Eorb','Allo','Real',kEOcc,2*nxOrb)
      kEVir=kEOcc+nxOrb
      call dcopy_(nxOrb,EOcc,1,Work(kEOcc),1)
      call dcopy_(nxOrb,EVir,1,Work(kEVir),1)

      Call GetMem('Saa','Allo','Real',ipSaa,nxOrb)
      call dcopy_(nxOrb,1.0d0,0,Work(ipSaa),1)


*     Occupied orbital selection                                       *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*nFro(iSym)
         call dcopy_(nBas(iSym)*nOcc(iSym),Work(ipXMO+jOff),1,
     &                                    Work(iXMO+kOff),1)
         call dcopy_(nBas(iSym)*nOcc(iSym),CMO(1+jOff),1,
     &                                    Work(iCMO+kOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nOcc(iSym)
      End Do
      ortho=.true.
*
      Call get_Orb_select(irc,Work(iCMO),Work(iXMO),Work(kEOcc),
     &                        Work(ipSQ),Work(ipSaa),Name,NamAct,
     &                        nSym,nActa,nOcc,nBas,ortho,Thrs,ns_O)
      If(irc.ne.0) Then
        Call qExit(ThisNm)
        Return
      End If
      iOff=0
      kOff=0
      Do iSym=1,nSym
         lOff=iOff+nBas(iSym)*nFro(iSym)
         Do ik=nOcc(iSym),1,-1
            jOff=kOff+nBas(iSym)*(ik-1)
            call dcopy_(nBas(iSym),Work(iCMO+jOff),1,
     &                            CMO(1+lOff),1)
            lOff=lOff+nBas(iSym)
         End Do
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nOcc(iSym)
      End Do
      iloc=0
      loff=0
      Do iSym=1,nSym
         Do ik=nOcc(iSym),1,-1
            ie=kEOcc+loff+ik-1
            iloc=iloc+1
            EOcc(iloc)=Work(ie)
         End Do
         loff=loff+nOcc(iSym)
      End Do
*
      If (all_Vir) Then
         Do iSym=1,nSym
            ns_V(iSym)=nExt(iSym)
         End Do
         goto 999
      EndIf
*
*     Virtual orbital selection                                        *
*----------------------------------------------------------------------*
      iOff=0
      kOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
         call dcopy_(nBas(iSym)*nExt(iSym),Work(ipXMO+jOff),1,
     &                                    Work(iXMO+kOff),1)
         call dcopy_(nBas(iSym)*nExt(iSym),CMO(1+jOff),1,
     &                                    Work(iCMO+kOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nExt(iSym)
      End Do
      ortho=.false.
*
      Call get_Vir_select(irc,Work(iCMO),Work(iXMO),Work(kEVir),
     &                        Work(ipSQ),Name,NamAct,iWork(iD_vir),
     &                        nSym,nActa,nExt,nBas,ortho,ns_V)
      If(irc.ne.0) Then
        Call qExit(ThisNm)
        Return
      End If
      Call GetMem('ID_vir','Free','Inte',iD_vir,nxBasT)
      iOff=0
      kOff=0
      Do iSym=1,nSym
         jOff=iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym))
         call dcopy_(nBas(iSym)*nExt(iSym),Work(iCMO+kOff),1,
     &                                    CMO(1+jOff),1)
         iOff=iOff+nBas(iSym)**2
         kOff=kOff+nBas(iSym)*nExt(iSym)
      End Do
      iloc=0
      loff=0
      Do iSym=1,nSym
         Do ik=1,nExt(iSym)
            ie=kEVir+loff+ik-1
            iloc=iloc+1
            EVir(iloc)=Work(ie)
         End Do
         loff=loff+nExt(iSym)
      End Do

999   Continue
      iDo=0
      jDo=0
      Do iSym=1,nSym  ! setup info
         lnOcc2(iSym)=nOcc(iSym)
         lnFro2(iSym)=nFro(iSym)
         lnDel2(iSym)=nDel(iSym)
         lnVir2(iSym)=nExt(iSym)
         lnOrb(iSym)=nBas(iSym)
         lnOcc(iSym)=nOcc(iSym)-ns_O(iSym)
         lnFro(iSym)=nFro(iSym)+ns_O(iSym)
         lnDel(iSym)=nDel(iSym)+ns_V(iSym)
         lnVir(iSym)=nExt(iSym)-ns_V(iSym)
         iDo=Max(iDo,lnOcc(iSym))
         jDo=Max(jDo,lnVir(iSym))
      End Do
      If (Min(iDo,jDo).eq.0) goto 1000

*     MP2 calculation on the Frozen region                             *
*----------------------------------------------------------------------*
      If (DoMP2) Then
*
         iloc=0
         jloc=0
         loff=0
         joff=0
         Do iSym=1,nSym
            Do ik=ns_V(iSym)+1,nExt(iSym)
               ie=loff+ik
               Work(kEVir+iloc)=EVir(ie)
               iloc=iloc+1
            End Do
            loff=loff+nExt(iSym)
            Do ik=1,lnOcc(iSym)
               ie=joff+ik
               Work(kEOcc+jloc)=EOcc(ie)
               jloc=jloc+1
            End Do
            joff=joff+nOcc(iSym)
         End Do
         iOff=0
         nVV=0
         nOA=0
         Do iSym=1,nSym
            kfr=1+iOff+nBas(iSym)*nFro(iSym)
            kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
            call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
            kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym)+ns_V(iSym))
            kto=kto+nBas(iSym)*lnOcc(iSym)
            call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
            iOff=iOff+nBas(iSym)**2
            nVV=nVV+lnVir(iSym)**2
            nOA=nOA+lnOcc(iSym)
         End Do
         Call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
         If (iSkip.gt.0) Then
          Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
          ip_Y=ip_X+nVV
          Call FZero(Work(ip_X),nVV+nOA)
          Call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,
     &                       ip_Y,.true.)
          Call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
          Call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,
     &                       ip_Y,.false.) ! compute energy and not Dab
          Call ChoMP2_Drv(irc,EFRO,Work(iCMO),Work(kEOcc),Work(kEVir))
          If(irc.ne.0) then
            write(6,*) 'Frozen region MP2 failed'
            Call Abend
          Else
            Write (6,'(A,F20.10,A)')
     &    ' Frozen region E2 contrib. = ',EFRO,' a.u.'
            EOSF=EOSMP2
            Write (6,'(A,F20.10,A)')
     &    ' (Opposite-Spin contrib.   = ',-EOSF,' )'
            Write (6,*)
          Endif
          iV=ip_X
          Do iSym=1,nSym
            TrF(iSym)=ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),1.0d0,0)
            iV=iV+lnVir(iSym)**2
          End Do
          Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
          Do iSym=1,nSym
             nOcc(iSym)=lnOcc2(iSym)
             nFro(iSym)=lnFro2(iSym)
             nDel(iSym)=lnDel2(iSym)
             nExt(iSym)=lnVir2(iSym)
          End Do
         EndIf

      EndIf
1000  Continue
*                                                                      *
*----------------------------------------------------------------------*
*
*     Update the nFro, nOcc, nExt, nDel for the Active site MP2
      Do iSym=1,nSym
         nFro(iSym)=nFro(iSym)+nOcc(iSym)-ns_O(iSym)
         nOcc(iSym)=ns_O(iSym)
         nDel(iSym)=nDel(iSym)+nExt(iSym)-ns_V(iSym)
         nExt(iSym)=ns_V(iSym)
         iDo=Max(iDo,nOcc(iSym))
         jDo=Max(jDo,nExt(iSym))
         nAuxO(iSym)=nBas(iSym)-nDel(iSym)
      End Do
      Call Put_iArray('nFroPT',nFro,nSym)
      Call Put_iArray('nDelPT',nDel,nSym)
      Call Put_iArray('nOrb',nAuxO,nSym)
*
      Call Check_Amp2(nSym,nOcc,nExt,iSkip)
      If (iSkip.gt.0) Then
         iOff=0
         lOff=0
         jk=1
         ja=1
         nVV=0
         nOA=0
         Do iSym=1,nSym
            Do ik=1,nOcc(iSym)
               kk=ik+iOff+lnOcc(iSym)
               EOcc(jk)=EOcc(kk)
               jk=jk+1
            End Do
            iOff=iOff+lnOcc(iSym)+nOcc(iSym)
            Do ia=1,nExt(iSym)
               ka=ia+lOff
               EVir(ja)=EVir(ka)
               ja=ja+1
            End Do
            lOff=lOff+lnVir(iSym)+nExt(iSym)
            nVV=nVV+nExt(iSym)**2
            nOA=nOA+nOcc(iSym)
         End Do
*
         Write(6,*)
         Write(6,'(A,8I4)')
     & ' Frozen orbitals after selection :  ',(nFro(i),i=1,nSym)
         Write(6,'(A,8I4)')
     & ' Occupied orbitals after selection: ',(nOcc(i),i=1,nSym)
         Write(6,*)
         Write(6,*) 'Energies of the active occupied orbitals '
         ii=0
         Do iSym=1,nSym
            If ( nOcc(iSym).ne.0 ) then
               Write(6,*)
               Write(6,'(A,I2,(T40,5F14.6))')
     &            ' symmetry species',iSym,(EOcc(ii+k),k=1,nOcc(iSym))
               ii=ii+nOcc(iSym)
            End If
         End Do
         Write(6,*)
*
         Write(6,*)
         Write(6,'(A,8I4)')
     & ' Secondary orbitals after selection:',(nExt(i),i=1,nSym)
         Write(6,'(A,8I4)')
     & ' Deleted orbitals after selection:  ',(nDel(i),i=1,nSym)
         Write(6,*)
         Write(6,*) 'Energies of the active virtual orbitals '
         ii=0
         Do iSym=1,nSym
            If ( nExt(iSym).ne.0 ) then
               Write(6,*)
               Write(6,'(A,I2,(T40,5F14.6))')
     &            ' symmetry species',iSym,(EVir(ii+k),k=1,nExt(iSym))
               ii=ii+nExt(iSym)
            End If
         End Do
         Write(6,*)
*
         Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
         ip_Y=ip_X+nVV
         Call FZero(Work(ip_X),nVV+nOA)
         Call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,ip_X,ip_Y,
     &                      .true.)
         Call ChoMP2_Drv(irc,Dummy,CMO,EOcc,EVir)
         Call LovMP2_putInf(nSym,lnOrb,nOcc,nFro,nDel,nExt,ip_X,ip_Y,
     &                      .false.)
         Wref=0.0d0
         Call ChoMP2_Drv(irc,EMP2,CMO,EOcc,EVir)
         If(irc.ne.0) then
           write(6,*) 'LovMP2 failed'
           Call Abend
         Endif
         iV=ip_X
         Do iSym=1,nSym
           TrA(iSym)=ddot_(nExt(iSym),Work(iV),1+nExt(iSym),1.0d0,0)
           iV=iV+nExt(iSym)**2
         End Do
         Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
      EndIf
      write(6,*)'------------------------------------------------------'
      write(6,*)' Symm.   Tr(D):  Active        Frozen         Full    '
      write(6,*)'------------------------------------------------------'
      STrA=0.0d0
      STrF=0.0d0
      STrX=0.0d0
      Do iSym=1,nSym
        write(6,'(2X,I4,10X,G10.4,4X,G10.4,4X,G10.4)') iSym,TrA(iSym),
     &       TrF(iSym),TrX(iSym)
        STrA=STrA+TrA(iSym)
        STrF=STrF+TrF(iSym)
        STrX=STrX+TrX(iSym)
      End Do
      write(6,*)'------------------------------------------------------'
      write(6,'(A,G10.4,4X,G10.4,4X,G10.4)')'          Sum:  ',
     & STrA,STrF,STrX
      write(6,*)'------------------------------------------------------'
      write(6,*)
*
      Write (6,'(A,F20.10,A)')
     &    ' Active region E2 contrib. = ',EMP2,' a.u.'
            Write (6,'(A,F20.10,A)')
     &    ' (Opposite-Spin contrib.   = ',-EOSMP2,' )'
      Write (6,*)
*
      XEMP2= EFRO
      EMP2 = EMP2 + EFRO
      EOSMP2 = EOSMP2 + EOSF
*
*     Update runfile for subsequent calcs (e.g., CHCC)
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff
         ito=kEOcc+koff+nFro(iSym)
         call dcopy_(nOcc(iSym),EOcc(ifr),1,Work(ito),1)
         ifr=1+joff
         ito=ito+nOcc(iSym)
         call dcopy_(nExt(iSym),EVir(ifr),1,Work(ito),1)
         ioff=ioff+nOcc(iSym)
         joff=joff+nExt(iSym)
         koff=koff+nBas(iSym)
      End Do
      Call Put_dArray('OrbE',Work(kEOcc),nOrb)
      Call Put_dArray('Last orbitals',CMO,nSQ)
*
      Call GetMem('LCMO','Free','Real',iCMO,3*nSQ)
      Call GetMem('Saa','Free','Real',ipSaa,nxOrb)
      Call GetMem('Eorb','Free','Real',kEOcc,2*nxOrb)

2000  Continue
      If (Min(iDo,jDo).eq.0) Then
         Write(6,*)
         Write(6,*)' None of the occupied or virtual orbitals has been '
         Write(6,*)' assigned to the Active region of the molecule.    '
         Write(6,*)' This is presumably NOT what you want !!!          '
         Write(6,*)' MP2 will Stop here. Bye Bye !! '
         Write(6,*)
         Call Abend
      EndIf

      CALL GetMem('SMAT','FREE','REAL',ipSQ,nSQ)
      Call qExit(ThisNm)
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Check_Amp2(nSym,nOcc,nVir,iSkip)

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
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine LovMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,
     &                         ip_X,ip_Y,isFNO)
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
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Get_Tr_Dab(nSym,nBas,nFro,nIsh,nSsh,nDel,
     &                      CMO,EOcc,EVir,TrD)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nIsh(nSym)
      Integer nSsh(nSym), nDel(nSym)
      Real*8  CMO(*), EOcc(*), EVir(*), TrD(nSym)
#include "WrkSpc.fh"
      Integer lnOrb(8), lnFro(8), lnOcc(8), lnVir(8), lnDel(8)
*
*
*
      nVV=0
      nBB=0
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnOrb(iSym)=nBas(iSym)
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnDel(iSym)=nDel(iSym)
         nVV=nVV+lnVir(iSym)**2
         nBB=nBB+nBas(iSym)**2
         nOA=nOA+lnOcc(iSym)
      End Do
*
      Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
      ip_Y=ip_X+nVV
      Call FZero(Work(ip_X),nVV+nOA)
*
      Call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y,
     &                   .true.)
      Call GetMem('CMON','Allo','Real',iCMO,nBB)
      Call FZero(Work(iCMO),nBB)
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
*
      Call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,Work(iCMO),EOcc,EVir)
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of occupied and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
      Call GetMem('CMON','Free','Real',iCMO,nBB)
*
      iV=ip_X
      Do iSym=1,nSym
        TrD(iSym)=ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),1.0d0,0)
        iV=iV+lnVir(iSym)**2
      End Do
      Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
*
      Return
      End
