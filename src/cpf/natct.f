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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE NATCT(H,LIC0)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION H(LIC0)
      Character*72 Header
*
      NSUM =0
      N2SUM=0
      n2Tri = 0
      nbMax = 0
      DO 9 ISYM=1,NSYM
         nbMax = Max(nbMax,nBas(iSym))
         NSUM =NSUM +NBAS(ISYM)
         N2SUM=N2SUM+NBAS(ISYM)**2
         n2Tri = n2Tri + nBas(iSym)*(nBas(iSym)+1)/2
9     CONTINUE

*     Read MO coefficients
      IDISK=ITOC17(1)
      CALL dDAFILE(Lu_TraOne,2,H(LW(87)),N2SUM,IDISK)
      IF (LW(87)+N2SUM-1.ge.LW(88)) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)')'*** ERROR IN SUBROUTINE NATCT ***'
         WRITE(6,'(6X,A)')'NO SPACE LEFT TO GENERATE FINAL ORBITALS'
         WRITE(6,*)
      CALL XFLUSH(6)
         Call Abend
      ENDIF
*
*     Loop over irreps and compute natural orbitals
*
      IOCC=LW(90)
      ICMO=LW(87)
      DO 10 M=1,NSYM
*        set occupation number of orbitals prefrozen in MOTRA
         CALL DCOPY_(NBAS(M),0.0D0,0,H(IOCC),1)
*        skip orbitals prefrozen in MOTRA
         CALL DCOPY_(NPFRO(M),2.0D0,0,H(IOCC),1)
         CALL NATORB(H(LW(62)),H(ICMO+NBAS(M)*NPFRO(M)),H(LW(88)),
     &               H(LW(89)),H(LW(89)),H(IOCC+NPFRO(M)),M)
         CALL DCOPY_(NORB(M)*NBAS(M),H(LW(89)),1,
     &              H(ICMO+NBAS(M)*NPFRO(M)),1)
         ICMO=ICMO+NBAS(M)**2
         IOCC=IOCC+NBAS(M)
10    CONTINUE

      LW91A = LW(91)
      LW91B = LW91A + n2Sum
      If (LW91B+n2Tri-1.gt.Lic) Then
         WRITE(6,*) ' Not enough core in NATCT'
      CALL XFLUSH(6)
         Call ErrTra
         Call Abend
      End If
      Call RelEne(ErelMV,ErelDC,nSym,nBas,H(LW(87)),
     *            H(LW(90)),H(LW91A),H(LW91B))

      EREL=ERELMV+ERELDC
      WRITE(6,'(/,5X,A)') 'FIRST ORDER RELATIVISTIC CORRECTIONS'
      WRITE(6,'(5X,A,F17.8)')
     *            'MASS-VELOCITY        ', ErelMV
      WRITE(6,'(5X,A,F17.8)')
     *            '1-EL DARWIN CONTACT  ', ErelDC
      WRITE(6,'(5X,A,F17.8)')
     *            'TOTAL REL. CORRECTION', Erel
      IF (ISDCI.EQ.1) THEN
        WRITE(6,'(5X,A,F17.8)')
     *            'REL. CI ENERGY       ', ETOT+Erel
        WRITE(6,'(5X,A,F17.8)')
     *            'REL. CI+Q ENERGY     ',DETOT+Erel
      ELSE
        WRITE(6,'(5X,A,F17.8)')
     *            'TOTAL REL. ENERGY    ', ETOT+Erel
      END IF
      CALL XFLUSH(6)
      CALL PRWF(H(LW(1)),H(LW(2)),H(LW(3)),H(LW(26)))
      If (iCPF.eq.1) Then
         Header=' CPF natural orbitals'
      Else If (iSDCI.eq.1) Then
         Header=' SDCI natural orbitals'
      Else If (iNCPF.eq.1) Then
         Header=' ACPF natural orbitals'
      Else
         Header=' MCPF natural orbitals'
      End If
      Call Primo(Header,.True.,.False.,1.0D-4,dum,nSym,nBas,nBas,
     *           Name,dum,H(LW(90)),H(LW(87)),-1)
*
*     Read the overlap matrix in ao basis
      iiRC=-1
      iOpt = 6
      Call RdOne(iiRC,iOpt,'MLTPL  0',1,H(LW(91)),iDum)
      If (iiRC.ne.0) Then
         Write (6,*) 'Natct: Error reading overlap matrix!'
         Call QTrace
         Call Abend
      End If
      Call Charge(nSym,nBas,Name,H(LW(87)),H(LW(90)),H(LW(91)),2,.True.,
     &            .True.)
      Call Prpt_old(nSym,nBas,nSum,n2Sum,H(LW(87)),H(LW(90)),
     *          Lic-LW(92)+1,H(LW(92)))
*
      If (iCPF.eq.1) Then
         Header='* CPF NO COEFS'
      Else If (iSDCI.eq.1) Then
         Header='* SDCI NO COEFS'
      Else If (iNCPF.eq.1) Then
         Header='* ACPF NO COEFS'
      Else
         Header='* MCPF NO COEFS'
      End If
      Call WrVec('CPFORB',Lu_CPFORB,'CO',nSym,nBas,nBas,
     & H(LW(87)), H(LW(90)), Dummy, iDummy, Header)
*
      RETURN
      END
