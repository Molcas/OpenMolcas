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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CheckAllIntegrals(QuitOnError,Silent,Mode,tau_,
     &                               MaxErr,MaxAbsErr,AverageErr,RMSErr,
     &                               MaxErr_OffDiagonal,
     &                               MaxAbsErr_OffDiagonal,
     &                               DiffNorm,DiffSum,ConvNorm,LDFNorm,
     &                               ConvSum,LDFSum)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compare two-electron integrals computed from LDF
C              coefficients (robust or nonrobust representation) to
C              conventional integrals. "All" means that all integrals
C              involving pairs of atom pairs, (AB|CD), are compared.
C
C     QuitOnError=T: quit if accuracy or CS violations are detected.
C     Silent=T: no printing (except for Cauchy-Schwarz and accuracy
C               violations).
C     Silent=F: print info for each pair of APs (AB|CD) as well as
C               total statistics (which is also returned to caller).
C     Mode=1: robust LDF integrals (2- and 3-index integrals)
C     Mode=2: nonrobust LDF integrals (2-index integrals)
C     Mode=3: half-and-half (3-index integrals)
C
      Implicit None
      Logical QuitOnError
      Logical Silent
      Integer Mode
      Real*8  tau_
      Real*8  MaxErr, MaxAbsErr, AverageErr, RMSErr
      Real*8  MaxErr_OffDiagonal, MaxAbsErr_OffDiagonal
      Real*8  DiffNorm, DiffSum
      Real*8  ConvNorm, LDFNorm
      Real*8  ConvSum, LDFSum
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"

      Character*21 SecNam
      Parameter (SecNam='LDF_CheckAllIntegrals')

      Real*8 Tol, Tol2
      Parameter (Tol=1.0d-12)

      Logical  LDF_IntegralPrescreeningInfoIsSet, isSymmetric
      External LDF_IntegralPrescreeningInfoIsSet, isSymmetric

      Integer  LDF_nBas_Atom, LDF_AtomPair_DiagDim
      External LDF_nBas_Atom, LDF_AtomPair_DiagDim

      Logical Verbose
      Logical IPI_set_here

      Integer AB, CD
      Integer ABmax, CDmax
      Integer nAB, nCD, nABCD
      Integer ip_Diff, l_Diff
      Integer ip_TOC, l_TOC
      Integer ip_TOCS, l_TOCS
      Integer ip_TOCSI, l_TOCSI
      Integer ip_DMax, l_DMax
      Integer ip_Map, l_Map
      Integer ip, l
      Integer k, ip0
      Integer nDiagErr
      Integer i2CF, nD2C, iRow, iCol
      Integer nRow, nCol

      Real*8  tau
      Real*8  TotDim, xABCD
      Real*8  ME, MAE, RMSE, DS
      Real*8  CNrm, LNrm
      Real*8  CSm, LSm
      Real*8  Tst, DABCD
      Real*8  nCS, nCSI
      Real*8  PairType(4)
      Real*8  MxD2C

      Integer i, j
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      Real*8  DMax
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      DMax(i)=Work(ip_DMax-1+i)

      Verbose=.not.Silent
      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
         Call LDF_SetIntegralPrescreeningInfo()
         IPI_set_here=.True.
      Else
         IPI_set_here=.False.
      End If

      l_TOC=NumberOfAtomPairs
      Call GetMem('TOC','Allo','Inte',ip_TOC,l_TOC)
      Call iZero(iWork(ip_TOC),l_TOC)
      l_TOCS=l_TOC
      Call GetMem('TOCS','Allo','Inte',ip_TOCS,l_TOCS)
      Call iZero(iWork(ip_TOCS),l_TOCS)
      l_TOCSI=l_TOC
      Call GetMem('TOCSI','Allo','Inte',ip_TOCSI,l_TOCSI)
      Call iZero(iWork(ip_TOCSI),l_TOCSI)

      l_DMax=NumberOfAtomPairs
      Call GetMem('Dmax','Allo','Real',ip_DMax,l_DMax)
      ip0=ip_DMax-1
      Do AB=1,NumberOfAtomPairs
         l=LDF_AtomPair_DiagDim(AB)
         ip=iWork(ip_AP_Diag-1+AB)
         Work(ip0+AB)=Work(ip)
         Do k=1,l-1
            Work(ip0+AB)=max(Work(ip0+AB),Work(ip+k))
         End Do
      End Do

      tau=max(tau_,0.0d0)

      MaxErr=0.0d0
      MaxAbsErr=0.0d0
      AverageErr=0.0d0
      RMSErr=0.0d0
      MaxErr_OffDiagonal=0.0d0
      MaxAbsErr_OffDiagonal=0.0d0
      DiffNorm=0.0d0
      DiffSum=0.0d0
      ConvNorm=0.0d0
      LDFNorm=0.0d0
      ConvSum=0.0d0
      LDFSum=0.0d0
      nCS=0.0d0
      nCSI=0.0d0
      PairType(1)=0.0d0
      PairType(2)=0.0d0
      PairType(3)=0.0d0
      PairType(4)=0.0d0
      ABmax=0
      CDmax=0
      nDiagErr=0

      TotDim=0.0d0
      Do AB=1,NumberOfAtomPairs
         nAB=LDF_nBas_Atom(AP_Atoms(1,AB))
     &      *LDF_nBas_Atom(AP_Atoms(2,AB))
         Do CD=1,AB
            nCD=LDF_nBas_Atom(AP_Atoms(1,CD))
     &         *LDF_nBas_Atom(AP_Atoms(2,CD))
            nABCD=nAB*nCD
            If (nABCD.gt.0) Then
               xABCD=dble(nABCD)
               TotDim=TotDim+xABCD
               l_Diff=nABCD
               Call GetMem('Diff','Allo','Real',ip_Diff,l_Diff)
               Call LDF_DiffIntegrals(Mode,tau,AB,CD,
     &                                l_Diff,Work(ip_Diff),
     &                                CNrm,LNrm,CSm,LSm)
               If (AB.eq.CD) Then
                  If (.not.isSymmetric(Work(ip_Diff),nAB,Tol)) Then
                     Call WarningMessage(2,
     &                               SecNam//': DiffInt not symmetric!')
                     Write(6,'(A,2I10,A)')
     &               'AB,CD=',AB,CD,'  DiffInt:'
                     Call Cho_Output(Work(ip_Diff),1,nAB,1,nAB,nAB,nAB,
     &                               1,6)
                     Call LDF_Quit(1)
                  End If
               End If
               If (AP_2CFunctions(1,AB).gt.0) Then
                  nRow=AP_2CFunctions(1,AB)
                  If (AP_Atoms(1,AB).eq.AP_Atoms(2,AB)) Then
                     nCol=2
                  Else
                     nCol=1
                  End If
                  l_Map=nRow*nCol
                  Call GetMem('ChkIMap','Allo','Inte',ip_Map,l_Map)
                  Call LDF_Map2CF(AB,nRow,nCol,iWork(ip_Map))
                  Tol2=min(Thr_Accuracy*1.0d-2,Tol)
                  nD2C=0
                  Do i2CF=1,AP_2CFunctions(1,AB)
                     iRow=iWork(ip_Map-1+i2CF)-1
                     MxD2C=0.0d0
                     Do k=0,nCD-1
                        MxD2C=max(MxD2C,abs(Work(ip_Diff+nAB*k+iRow)))
                     End Do
                     If (MxD2C.gt.Tol2) Then
                        Write(6,'(A,2I10,A,I10,A,1P,D20.10)')
     &                  'AB,CD=',AB,CD,
     &                  'Integrals corresponding to AB 2C function',
     &                  i2CF,' are not exact. Max abs error=',
     &                  MxD2C
                        nD2C=nD2C+1
                     End If
                  End Do
                  If (AP_Atoms(1,AB).eq.AP_Atoms(2,AB)) Then
                     Do i2CF=1,AP_2CFunctions(1,AB)
                        iRow=iWork(ip_Map-1+nRow+i2CF)-1
                        MxD2C=0.0d0
                        Do k=0,nCD-1
                           MxD2C=max(MxD2C,
     &                               abs(Work(ip_Diff+nAB*k+iRow)))
                        End Do
                        If (MxD2C.gt.Tol2) Then
                           Write(6,'(A,2I10,A,I10,A,1P,D20.10)')
     &                     'AB,CD=',AB,CD,
     &                     'Integrals corresponding to AB 2C function',
     &                     i2CF,' are not exact. Max abs error=',
     &                     MxD2C
                           nD2C=nD2C+1
                        End If
                     End Do
                  End If
                  Call GetMem('ChkIMap','Free','Inte',ip_Map,l_Map)
                  If (nD2C.ne.0) Then
                     Call WarningMessage(2,SecNam//
     &                    ': Integrals corresponding to 2C functions'//
     &                    ' are not exact!')
                     Call LDF_Quit(1)
                  End If
               End If
               If (AB.ne.CD .and. AP_2CFunctions(1,CD).gt.0) Then
                  nRow=AP_2CFunctions(1,CD)
                  If (AP_Atoms(1,CD).eq.AP_Atoms(2,CD)) Then
                     nCol=2
                  Else
                     nCol=1
                  End If
                  l_Map=nRow*nCol
                  Call GetMem('ChkIMap','Allo','Inte',ip_Map,l_Map)
                  Call LDF_Map2CF(CD,nRow,nCol,iWork(ip_Map))
                  Tol2=min(Thr_Accuracy*1.0d-2,Tol)
                  nD2C=0
                  Do i2CF=1,AP_2CFunctions(1,CD)
                     iCol=iWork(ip_Map-1+i2CF)-1
                     MxD2C=0.0d0
                     Do k=0,nAB-1
                        MxD2C=max(MxD2C,abs(Work(ip_Diff+nAB*iCol+k)))
                     End Do
                     If (MxD2C.gt.Tol2) Then
                        Write(6,'(A,2I10,A,I10,A,1P,D20.10)')
     &                  'AB,CD=',AB,CD,
     &                  'Integrals corresponding to AB 2C function',
     &                  i2CF,' are not exact. Max abs error=',
     &                  MxD2C
                        nD2C=nD2C+1
                     End If
                  End Do
                  If (AP_Atoms(1,CD).eq.AP_Atoms(2,CD)) Then
                     Do i2CF=1,AP_2CFunctions(1,CD)
                        iCol=iWork(ip_Map-1+nRow+i2CF)-1
                        MxD2C=0.0d0
                        Do k=0,nAB-1
                           MxD2C=max(MxD2C,
     &                               abs(Work(ip_Diff+nAB*iCol+k)))
                        End Do
                        If (MxD2C.gt.Tol2) Then
                           Write(6,'(A,2I10,A,I10,A,1P,D20.10)')
     &                     'AB,CD=',AB,CD,
     &                     'Integrals corresponding to AB 2C function',
     &                     i2CF,' are not exact. Max abs error=',
     &                     MxD2C
                           nD2C=nD2C+1
                        End If
                     End Do
                  End If
                  Call GetMem('ChkIMap','Free','Inte',ip_Map,l_Map)
                  If (nD2C.ne.0) Then
                     Call WarningMessage(2,SecNam//
     &                    ': Integrals corresponding to 2C functions'//
     &                    ' are not exact!')
                     Call LDF_Quit(1)
                  End If
               End If
               ME=0.0d0
               MAE=0.0d0
               RMSE=0.0d0
               DS=0.0d0
               ip0=ip_Diff-1
               Do k=1,nABCD
                  ME=max(ME,Work(ip0+k))
                  MAE=max(MAE,abs(Work(ip0+k)))
                  RMSE=RMSE+Work(ip0+k)**2
                  DS=DS+Work(ip0+k)
               End Do
               If (AB.ne.CD) Then
                  MaxErr_OffDiagonal=max(MaxErr_OffDiagonal,ME)
                  MaxAbsErr_OffDiagonal=max(MaxAbsErr_OffDiagonal,ME)
               End If
               MaxErr=max(MaxErr,ME)
               If (MAE.gt.MaxAbsErr) Then
                  MaxAbsErr=MAE
                  ABmax=AB
                  CDmax=CD
               End If
               AverageErr=AverageErr+DS
               RMSErr=RMSErr+RMSE
               DiffNorm=DiffNorm+RMSE
               DiffSum=DiffSum+DS
               ConvNorm=ConvNorm+CNrm**2
               LDFNorm=LDFNorm+LNrm**2
               ConvSum=ConvSum+CSm
               LDFSum=LDFSum+LSm
               If (Verbose) Then
                  Write(6,'(A,2I10,1X,A,4I10,A)')
     &            'LDF integral error analysis for AB,CD=',AB,CD,
     &            '(A,B,C,D=',AP_Atoms(1,AB),AP_Atoms(2,AB),
     &             AP_Atoms(1,CD),AP_Atoms(2,CD),'):'
                  Write(6,'(3X,A,10X,I10)')
     &            'Dimension...........',nABCD
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Conventional norm...',CNrm,
     &            'Conventional sum....',CSm
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'LDF norm............',LNrm,
     &            'LDF sum.............',LSm
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Difference norm.....',sqrt(RMSE),
     &            'Difference sum......',DS
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Max Error...........',ME,
     &            'Max Abs Error.......',MAE
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Average Error.......',DS/xABCD,
     &            'RMS Error...........',sqrt(RMSE/xABCD)
                  Call xFlush(6)
               End If
               If (AB.eq.CD) Then
                  Tst=MAE-DMax(AB)
                  If (abs(Tst).gt.1.0d-8) Then
                     Write(6,'(A,2I10,1P,3(1X,A,D20.10))')
     &               'AB,CD=',AB,CD,'MAE=',MAE,'DMax=',DMax(AB),
     &               'Difference=',Tst
                     Call xFlush(6)
                     nDiagErr=nDiagErr+1
                  End If
                  DABCD=DMax(AB)
               Else
                  DABCD=sqrt(DMax(AB)*DMax(CD))
               End If
               Tst=MAE-DABCD
               If (Tst.gt.0.0d0) Then
                  If (DABCD.gt.1.0d-8) Then ! DABCD>1d-8
                     If ((Tst/DABCD).gt.1.0d-2) Then ! 1 pct allowed
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [SIGNIFICANT]'
                        nCS=nCS+1.0d0
                        iWork(ip_TOCS-1+AB)=iWork(ip_TOCS-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCS-1+CD)=iWork(ip_TOCS-1+CD)+1
                        End If
                     Else
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [INSIGNIFICANT]'
                        nCSI=nCSI+1.0d0
                        iWork(ip_TOCSI-1+AB)=iWork(ip_TOCSI-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCSI-1+CD)=iWork(ip_TOCSI-1+CD)+1
                        End If
                     End If
                  Else If (DABCD.gt.1.0d-12) Then ! 1d-8>=DABCD>1d-12
                     If ((Tst/DABCD).gt.5.0d-1) Then ! 50 pct allowed
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [SIGNIFICANT]'
                        nCS=nCS+1.0d0
                        iWork(ip_TOCS-1+AB)=iWork(ip_TOCS-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCS-1+CD)=iWork(ip_TOCS-1+CD)+1
                        End If
                     Else
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [INSIGNIFICANT]'
                        nCSI=nCSI+1.0d0
                        iWork(ip_TOCSI-1+AB)=iWork(ip_TOCSI-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCSI-1+CD)=iWork(ip_TOCSI-1+CD)+1
                        End If
                     End If
                  Else If (DABCD.gt.1.0d-15) Then ! 1d-12>=DABCD>1d-15
                     If (Tst.gt.1.0d-10) Then
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [SIGNIFICANT]'
                        nCS=nCS+1.0d0
                        iWork(ip_TOCS-1+AB)=iWork(ip_TOCS-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCS-1+CD)=iWork(ip_TOCS-1+CD)+1
                        End If
                     Else
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [INSIGNIFICANT]'
                        nCSI=nCSI+1.0d0
                        iWork(ip_TOCSI-1+AB)=iWork(ip_TOCSI-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCSI-1+CD)=iWork(ip_TOCSI-1+CD)+1
                        End If
                     End If
                  Else ! 1d-15>=DABCD
                     If (Tst.gt.1.0d-11) Then
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [SIGNIFICANT]'
                        nCS=nCS+1.0d0
                        iWork(ip_TOCS-1+AB)=iWork(ip_TOCS-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCS-1+CD)=iWork(ip_TOCS-1+CD)+1
                        End If
                     Else
                        Write(6,
     &                  '(3X,A,2I10,1P,3(1X,A,D20.10),A)')
     &                  'CS violation: AB,CD=',AB,CD,
     &                  'MAE=',MAE,'CS upper bound=',DABCD,
     &                  'Difference=',Tst,' [INSIGNIFICANT]'
                        nCSI=nCSI+1.0d0
                        iWork(ip_TOCSI-1+AB)=iWork(ip_TOCSI-1+AB)+1
                        If (CD.ne.AB) Then
                           iWork(ip_TOCSI-1+CD)=iWork(ip_TOCSI-1+CD)+1
                        End If
                     End If
                  End If
                  Call xFlush(6)
               End If
               If (LDF2) Then
                  If (MAE.gt.Thr_Accuracy) Then
                     Write(6,'(3X,A,2I10,A,1P,D20.10,A,A,D20.10)')
     &               'AB,CD=',AB,CD,' Max Abs Error=',MAE,' > ',
     &               'target accuracy=',Thr_Accuracy
                     Call xFlush(6)
                     iWork(ip_TOC-1+AB)=iWork(ip_TOC-1+AB)+1
                     If (CD.ne.AB) Then
                        iWork(ip_TOC-1+CD)=iWork(ip_TOC-1+CD)+1
                     End If
                     If (AP_2CFunctions(1,AB).gt.0 .or.
     &                   AP_2CFunctions(1,CD).gt.0) Then
                        If (AP_1CLinDep(1,AB).gt.0 .or.
     &                      AP_1CLinDep(1,CD).gt.0) Then
                           PairType(4)=PairType(4)+1.0d0
                        Else
                           PairType(3)=PairType(3)+1.0d0
                        End If
                     Else
                        If (AP_1CLinDep(1,AB).gt.0 .or.
     &                      AP_1CLinDep(1,CD).gt.0) Then
                           PairType(2)=PairType(2)+1.0d0
                        Else
                           PairType(1)=PairType(1)+1.0d0
                        End If
                     End If
                  End If
               End If
               Call GetMem('Diff','Free','Real',ip_Diff,l_Diff)
            Else
               If (Verbose) Then
                  Write(6,'(A,2I10,1X,A,4I10,A)')
     &            'LDF integral error analysis for AB,CD=',AB,CD,
     &            '(A,B,C,D=',AP_Atoms(1,AB),AP_Atoms(2,AB),
     &             AP_Atoms(1,CD),AP_Atoms(2,CD),'):'
                  Write(6,'(3X,A,10X,I10)')
     &            'Dimension...........',nABCD
                  xABCD=0.0d0
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Conventional norm...',xABCD,
     &            'Conventional sum....',xABCD
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'LDF norm............',xABCD,
     &            'LDF sum.............',xABCD
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Difference norm.....',xABCD,
     &            'Difference sum......',xABCD
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Max Error...........',xABCD,
     &            'Max Abs Error.......',xABCD
                  Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &            'Average Error.......',xABCD,
     &            'RMS Error...........',xABCD
                  Call xFlush(6)
               End If
            End If
         End Do
      End Do

      If (TotDim.gt.0.0d0) Then
         AverageErr=AverageErr/TotDim
         RMSErr=sqrt(RMSErr/TotDim)
      End If
      DiffNorm=sqrt(DiffNorm)
      ConvNorm=sqrt(ConvNorm)
      LDFNorm=sqrt(LDFNorm)

      If (Verbose) Then
         Call Cho_Head('LDF Integral Error Analysis','-',80,6)
         If (Mode.eq.1) Then
            Write(6,'(3X,A)') 'LDF integral representation: ROBUST'
         Else If (Mode.eq.2) Then
            Write(6,'(3X,A)') 'LDF integral representation: NONROBUST'
         Else If (Mode.eq.3) Then
            Write(6,'(3X,A)')
     &      'LDF integral representation: HALF-AND-HALF'
         Else
            Write(6,'(3X,A)') 'LDF integral representation: UNKNOWN'
         End If
         Write(6,'(3X,A,1P,D20.10)')
     &   'LDF integral prescreening threshold:',tau
         Write(6,'(3X,A,1P,D20.10)')
     &   'Dimension...........',TotDim
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &   'Conventional norm...',ConvNorm,
     &   'Conventional sum....',ConvSum
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &   'LDF norm............',LDFNorm,
     &   'LDF sum.............',LDFSum
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &   'Difference norm.....',DiffNorm,
     &   'Difference sum......',DiffSum
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10,A,2I10,A)')
     &   'Max Error...........',MaxErr,
     &   'Max Abs Error.......',MaxAbsErr,'  (AB,CD=',ABmax,CDmax,')'
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10,A,2I10,A)')
     &   'Max OffD Error......',MaxErr_OffDiagonal,
     &   'Max Abs OffD Error..',MaxAbsErr_OffDiagonal
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &   'Average Error.......',AverageErr,
     &   'RMS Error...........',RMSErr
         Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &   'CS Violations.......',nCS,
     &   'CS Near-Violations..',nCSI
         If (LDF2) Then
            Write(6,'(3X,A,1P,D20.10)')
     &      'Accuracy violations.',PairType(1)+PairType(2)
     &                            +PairType(3)+PairType(4)
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &      '1C only, no LinDep..',PairType(1),
     &      '1C only, LinDep.....',PairType(2)
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &      '1C+2C, no LinDep....',PairType(3),
     &      '1C+2C, LinDep.......',PairType(4)
         End If
         Call xFlush(6)
      Else
         If (nCS.gt.0.0d0 .or.
     &       (PairType(1)+PairType(2)+PairType(3)+PairType(4)).gt.0.0d0)
     &   Then
            Call Cho_Head('LDF Integral Error Analysis','-',80,6)
            If (Mode.eq.1) Then
               Write(6,'(3X,A)') 'LDF integral representation: ROBUST'
            Else If (Mode.eq.2) Then
               Write(6,'(3X,A)')
     &         'LDF integral representation: NONROBUST'
            Else If (Mode.eq.3) Then
               Write(6,'(3X,A)')
     &         'LDF integral representation: HALF-AND-HALF'
            Else
               Write(6,'(3X,A)') 'LDF integral representation: UNKNOWN'
            End If
            Write(6,'(3X,A,1P,D20.10)')
     &      'LDF integral prescreening threshold:',tau
            Write(6,'(3X,A,1P,D20.10)')
     &      'Dimension...........',TotDim
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &      'CS Violations.......',nCS,
     &      'CS Near-Violations..',nCSI
            If (LDF2) Then
               Write(6,'(3X,A,1P,D20.10)')
     &         'Accuracy violations.',PairType(1)+PairType(2)
     &                               +PairType(3)+PairType(4)
               Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &         '1C only, no LinDep..',PairType(1),
     &         '1C only, LinDep.....',PairType(2)
               Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)')
     &         '1C+2C, no LinDep....',PairType(3),
     &         '1C+2C, LinDep.......',PairType(4)
            End If
         End If
      End If

      If (nCSI.gt.0.0d0) Then
         Do AB=1,NumberOfAtomPairs
            If (iWork(ip_TOCSI-1+AB).gt.0) Then
               Write(6,'(3X,A,I10,A,I10,A)')
     &         'Atom pair',AB,' is involved in',iWork(ip_TOCSI-1+AB),
     &         ' insignificant Cauchy-Schwarz violations'
            End If
         End Do
         Write(6,*)
         Call xFlush(6)
      End If

      If (nCS.gt.0.0d0) Then
         Do AB=1,NumberOfAtomPairs
            If (iWork(ip_TOCS-1+AB).gt.0) Then
               Write(6,'(3X,A,I10,A,I10,A)')
     &         'Atom pair',AB,' is involved in',iWork(ip_TOCS-1+AB),
     &         ' Cauchy-Schwarz violations'
            End If
         End Do
         Write(6,*)
         Call xFlush(6)
      End If

      If ((PairType(1)+PairType(2)+PairType(3)+PairType(4)).gt.0.0d0)
     & Then
         Do AB=1,NumberOfAtomPairs
            If (iWork(ip_TOC-1+AB).gt.0) Then
               Write(6,'(3X,A,I10,A,I10,A)')
     &         'Atom pair',AB,' is involved in',iWork(ip_TOC-1+AB),
     &         ' accuracy violations'
            End If
         End Do
         Write(6,*)
         Call xFlush(6)
      End If

      Call GetMem('Dmax','Free','Real',ip_DMax,l_DMax)
      Call GetMem('TOCSI','Free','Inte',ip_TOCSI,l_TOCSI)
      Call GetMem('TOCS','Free','Inte',ip_TOCS,l_TOCS)
      Call GetMem('TOC','Free','Inte',ip_TOC,l_TOC)
      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      If (QuitOnError) Then
         If ((PairType(1)+PairType(2)+PairType(3)+PairType(4)).gt.0.0d0)
     &   Then
            Call WarningMessage(2,
     &                  SecNam//': Accuracy violations detected!')
            Call LDF_Quit(1)
         End If
         If (nCS.gt.0.0d0) Then
            Call WarningMessage(2,
     &                  SecNam//': Cauchy-Schwarz violations detected!')
            Call LDF_Quit(1)
         End If
         If (ABmax.ne.CDmax) Then
            Call WarningMessage(2,SecNam//
     &         ': Max Abs Error should be on the diagonal - it is not!')
            Call LDF_Quit(1)
         End If
         If (nDiagErr.ne.0) Then
            Call WarningMessage(2,SecNam//': off-diagonal MAEs found!')
            Call LDF_Quit(1)
         End If
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_DiffIntegrals(Mode,tau,AB,CD,l_DiffInt,DiffInt,
     &                             ConvNorm,LDFNorm,ConvSum,LDFSum)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute difference between conventional integrals and
C              two-electron integrals computed from LDF coefficients
C              (robust or nonrobust representation).
C
      Implicit None
      Integer Mode
      Real*8  tau
      Integer AB, CD
      Integer l_DiffInt
      Real*8  DiffInt(l_DiffInt)
      Real*8  ConvNorm, LDFNorm
      Real*8  ConvSum, LDFSum
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*17 SecNam
      Parameter (SecNam='LDF_DiffIntegrals')

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      real*8 ddot_
      external ddot_

      Logical IPI_set_here

      Integer nABCD
      Integer ip_LDFInt, l_LDFInt
      Integer k, ip0

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      nABCD=LDF_nBas_Atom(AP_Atoms(1,AB))
     &     *LDF_nBas_Atom(AP_Atoms(2,AB))
     &     *LDF_nBas_Atom(AP_Atoms(1,CD))
     &     *LDF_nBas_Atom(AP_Atoms(2,CD))
      If (nABCD.gt.0) Then
         If (nABCD.gt.l_DiffInt) Then
            Call WarningMessage(2,
     &                         SecNam//': insufficient array dimension')
            Call LDF_Quit(1)
         End If
         Call LDF_ComputeValenceIntegrals(AB,CD,l_DiffInt,DiffInt)
         ConvNorm=sqrt(dDot_(nABCD,DiffInt,1,DiffInt,1))
         ConvSum=DiffInt(1)
         Do k=2,nABCD
            ConvSum=ConvSum+DiffInt(k)
         End Do
         If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
            Call LDF_SetIntegralPrescreeningInfo()
            IPI_set_here=.True.
         Else
            IPI_set_here=.False.
         End If
         l_LDFInt=nABCD
         Call GetMem('LDFInt','Allo','Real',ip_LDFInt,l_LDFInt)
         Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,CD,
     &                                         l_LDFInt,Work(ip_LDFInt))
         LDFNorm=sqrt(dDot_(nABCD,Work(ip_LDFInt),1,Work(ip_LDFInt),1))
         ip0=ip_LDFInt-1
         LDFSum=Work(ip0+1)
         Do k=2,nABCD
            LDFSum=LDFSum+Work(ip0+k)
         End Do
         Call dAXPY_(nABCD,-1.0d0,Work(ip_LDFInt),1,DiffInt,1)
         Call GetMem('LDFInt','Free','Real',ip_LDFInt,l_LDFInt)
         If (IPI_set_here) Then
            Call LDF_UnsetIntegralPrescreeningInfo()
         End If
      Else
         ConvNorm=0.0d0
         LDFNorm=0.0d0
         ConvSum=0.0d0
         LDFSum=0.0d0
      End If

      End
