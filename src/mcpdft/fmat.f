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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine Fmat_m(CMO,PUVX,D,D1A,FI,FA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Update the Fock matrix for the active orbitals and transform     *
*     it to MO basis as well as the matrix FI (Fock matrix) for        *
*     frozen and inactive orbitals).                                   *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : array of real*8                                        *
*               MO-coefficients                                        *
*     PUVX    : array of real*8                                        *
*               two-electron integrals (pu!vx)                         *
*     D       : array of real*8                                        *
*               averaged one-body density matrix                       *
*     D1A     : array of real*8                                        *
*               active one body density matrix in AO-basis             *
*     FI      : array of real*8                                        *
*               inactive Fock matrix                                   *
*     FA      : array of real*8                                        *
*               active Fock matrix                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated for MOLCAS version 2                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1991               *
*     - updated for MOLCAS version 3                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1993               *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='FMAT    ')
#include "WrkSpc.fh"

      Dimension CMO(*) , PUVX(*) , D(*) , D1A(*) , FI(*) , FA(*)

      Call qEnter ('Fmat')
C Local print level (if any)
      IPRLEV=IPRLOC(4)
      !iPrLev=DEBUG-1
      If ( iPrLev.ge.DEBUG ) then
        write(6,*) ('*',i=1,65)
        write(6,*) 'Entering FMAT routine called by MSCTL!'
        write(6,*) ('*',i=1,65)
        write(6,*) 'printing input matrices :'
        write(6,*) ('*',i=1,65)
        Write(LF,*)
        Write(LF,*) ' CMOs in FMAT'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
           iBas = nBas(iSym)
           call wrtmat(CMO(ioff),iBas,iBas, iBas, iBas)
           iOff = iOff + iBas*iBas
         End Do

         Write(LF,*)
         Write(LF,*) ' PUVX in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         call wrtmat(PUVX,1,nFint, 1, nFint)

         Write(LF,*)
         Write(LF,*) ' ---------------------'
        CALL TRIPRT('Averaged one-body density matrix D, in MO in FMAT',
     &              ' ',D,NAC)

         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(D1A(iOff),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End DO

         Write(LF,*)
         Write(LF,*) ' FI in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FI(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do

         Write(LF,*)
         Write(LF,*) ' FA in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FA(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do
       End If

*************************************************************
* Here we should start the real work!
*************************************************************
C Local print level (if any)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*     create FA in AO basis
      Call GetMem('Scr1','Allo','Real',iTmp1,nTot1)
      Call Fold(nSym,nBas,D1A,Work(iTmp1))
      If(KSDFT.ne.'SCF'
     & .and. KSDFT.ne.'TLSDA'
     & .and. KSDFT.ne.'FTLSDA'
     & .and. KSDFT.ne.'TPBE'
     & .and. KSDFT.ne.'FTPBE'
     & .and. KSDFT.ne.'TBLYP'
     & .and. KSDFT.ne.'FTBLYP'
     & .and. KSDFT.ne.'TREVPBE'
     & .and. KSDFT.ne.'FTREVPBE') NewFock=0

      If (NewFock.eq.0) Then
         nBMX=0
         Do iSym=1,nSym
            nBMX=Max(nBMX,nBas(iSym))
         End Do
         Call FZero(FA,nTot1)
!         write(6,*) 'ExFac : ', ExFac
         Call FTwo_Drv(nSym,nBas,nAsh,nSkipX,
     &                    Work(iTmp1),D1A,FA,nTot1,
     &                    ExFac,nTot2,nBMX,CMO)
      End If

*     Inactive-active contribution to ECAS
      VIA=dDot_(nTot1,FI,1,Work(iTmp1),1)
      ECAS=EMY+VIA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A,E20.10)') ' Total core energy:            ',EMY
        Write(LF,'(A,E20.10)') ' inactive-active interaction:  ',VIA
        Write(LF,'(A,E20.10)') ' CAS energy (core+interaction):',ECAS
      End If
      Call GetMem('Scr1','Free','Real',iTmp1,nTot1)

*     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

*     transform FI from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FI(iOff1),Work(iTmp1),1,iBas,iBas)
C        Call MXMA(Work(iTmp1),1,iBas,
C     &            CMO(iOff2+(iFro*iBas)),1,iBas,
C     &            Work(iTmp2),1,iBas,
C     &            iBas,iBas,iOrb)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Work(iTmp1),iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Work(iTmp2),iBas)
        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FI(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

*     transform FA from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FA(iOff1),Work(iTmp1),1,iBas,iBas)
C        Call MXMA(Work(iTmp1),1,iBas,
C     &            CMO(iOff2+(iFro*iBas)),1,iBas,
C     &            Work(iTmp2),1,iBas,
C     &            iBas,iBas,iOrb)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Work(iTmp1),iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Work(iTmp2),iBas)
        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FA(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

c**************************************************************************
c              Add DFT part to Fock matrix:                               *
c**************************************************************************
      If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM'.and.
     & (KSDFT(1:5).ne.'TLSDA'.and. !GLM
     &  KSDFT(1:5).ne.'TBLYP'.and.
     &  KSDFT(1:6).ne.'FTLSDA'.and.
     &  KSDFT(1:5).ne.'FTPBE'.and.
     &  KSDFT(1:7).ne.'TREVPBE'.and.
     &  KSDFT(1:8).ne.'FTREVPBE'.and.
     &  KSDFT(1:6).ne.'FTBLYP'.and.
     &  KSDFT(1:4).ne.'TPBE')) Then
        ipTmpFckI=-99999
        ipTmpFckA=-99999
        Call Get_dExcdRa(ipTmpFck,nTmpFck)
        If(nTmpFck.eq.NTOT1) Then
           ipTmpFckI=ipTmpFck
        Else If(nTmpFck.eq.2*NTOT1) Then
           ipTmpFckI=ipTmpFck
           ipTmpFckA=ipTmpFck+nTot1
        Else
           Write(LF,*) ' Somethings wrong in dim. DFT',nTmpFck
           Call Abend()
        End If
        Call GetMem('ScrD1a','Allo','Real',iTmpD1A,nTot1)
        Call Fold(nSym,nBas,D1A,Work(iTmpD1A))
        VIA_DFT=dDot_(nTot1,Work(ipTmpFckI),1,Work(iTmpD1A),1)
        Call GetMem('ScrD1a','Free','Real',iTmpD1A,nTot1)
*
*          Transform alpha density from AO to MO
*
        iOff1 = 1
        iOff2 = 1
        iOff3 = 1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          iOrb = nOrb(iSym)
          iFro = nFro(iSym)
          Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
          Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
          Call Square(Work(ipTmpFckI+iOff1-1),
     &                Work(iTmp1),1,iBas,iBas)
C          Call MXMA(Work(iTmp1),1,iBas,
C     &              CMO(iOff2+(iFro*iBas)),1,iBas,
C     &              Work(iTmp2),1,iBas,
C     &              iBas,iBas,iOrb)
          Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Work(iTmp1),iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Work(iTmp2),iBas)
          Call MXMT(Work(iTmp2),iBas,1,
     &              CMO(iOff2+(iFro*iBas)),1,iBas,
     &              Work(ipTmpFckI+iOff3-1),
     &              iOrb,iBas)
          Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
          Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
          iOff1 = iOff1 + (iBas*iBas+iBas)/2
          iOff2 = iOff2 + iBas*iBas
          iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
        End Do
*
*          Transform Active DFT Fock from AO to MO
*
        If(ipTmpFckA.ne.-99999) Then
        iOff1 = 1
        iOff2 = 1
        iOff3 = 1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          iOrb = nOrb(iSym)
          iFro = nFro(iSym)
          Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
          Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
          Call Square(Work(ipTmpFckA+iOff1-1),
     &                Work(iTmp1),1,iBas,iBas)
C          Call MXMA(Work(iTmp1),1,iBas,
C     &              CMO(iOff2+(iFro*iBas)),1,iBas,
C     &              Work(iTmp2),1,iBas,
C     &              iBas,iBas,iOrb)
          Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Work(iTmp1),iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Work(iTmp2),iBas)
          Call MXMT(Work(iTmp2),iBas,1,
     &              CMO(iOff2+(iFro*iBas)),1,iBas,
     &              Work(ipTmpFckA+iOff3-1),
     &              iOrb,iBas)
          Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
          Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
          iOff1 = iOff1 + (iBas*iBas+iBas)/2
          iOff2 = iOff2 + iBas*iBas
          iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
        End Do
        End If
*
c        If(DFTFOCK(1:4).ne.'ROKS') Then
c          Write(LF,*) ' Just add a,b to FA,FI',DFTFOCK(1:4)
c        Else
c          Write(LF,*) ' ROKS formula',DFTFOCK(1:4)
c        End If
*
        If(DFTFOCK(1:4).ne.'ROKS') Then
          call daxpy_(NTOT1,1.0D0,Work(ipTmpFckI),1,FI,1)
          If(ipTmpFckA.ne.-99999)
     &    call daxpy_(NTOT1,1.0D0,Work(ipTmpFckA),1,FA,1)
        Else If (DFTFOCK(1:4).eq.'ROKS') Then
           iOff1 = 0
           Do iSym = 1,nSym
              Do iOrb=1,nOrb(iSym)
                 Do jOrb=1,iOrb
                    ij=iOff1+iOrb*(iOrb-1)/2+jOrb
                    If(iOrb.le.nIsh(iSym)) Then
                      FI(ij)=FI(ij)+0.5d0*
     &                  (Work(ipTmpFckI+ij-1)+Work(ipTmpFckA+ij-1))
                    End If
                    If (iOrb.gt.nIsh(iSym).and.
     &                  iOrb.le.nIsh(iSym)+nAsh(iSym)) Then
                       If (jOrb.le.nIsh(iSym)) Then
                          FI(ij)=FI(ij)+Work(ipTmpFckA+ij-1)
                       Else
                          FI(ij)=FI(ij)+0.5d0*(Work(ipTmpFckI+ij-1)+
     &                                         Work(ipTmpFckA+ij-1))
                       End If
                    End If
                    If (iOrb.gt.nIsh(iSym)+nAsh(iSym)) Then
                       If(jOrb.gt.nIsh(iSym).and.
     &                    jOrb.le.nIsh(iSym)+nAsh(iSym)) Then
                          FI(ij)=FI(ij)+Work(ipTmpFckI+ij-1)
                       Else
                          FI(ij)=FI(ij)+0.5d0*(Work(ipTmpFckI+ij-1)+
     &                                         Work(ipTmpFckA+ij-1))
                       End If

                    End If
                 End Do
              End Do
              iOff1 = iOff1 + (nOrb(iSym)*nOrb(iSym)+nOrb(iSym))/2
           End Do
        Else
           Write(LF,*) " Not implemented yet"
        End If
        Call Free_Work(ipTmpFck)
      End If
***************************************************************************
*     update Fock matrix
      If (NewFock.eq.1) Call Upd_FA_m(PUVX,FA,D,ExFac)

*     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

      Call qExit('Fmat')

      Return
      End
