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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CheckPSD_Full(Diagonalization,Mode,
     &                             UseExactIntegrals,tau,irc)
C
C     Thomas Bondo Pedersen, January 2012.
C
C     Purpose: Check if full integral matrix is positive semidefinite.
C
C              Diagonalization: .True.: eigenvalues, .False.: CD
C              Mode=0: conventional integrals
C                  =1: robust LDF integrals
C                  =2: nonrobust LDF integrals
C                  =3: half-and-half LDF integrals
C              UseExactIntegrals=0: do not use exact integrals (Mode>0)
C                  =1: use exact integrals in diagonal blocks
C                  =2: use exact integrals in off-diagonal blocks
C              tau: LDF prescreening threshold (not used when Mode=0)
C              irc: return code (0 if PSD)
C
C     WARNING: requires a LOT of CPU and memory....
C
      Implicit None
      Logical Diagonalization
      Integer Mode
      Integer UseExactIntegrals
      real*8  tau
      Integer irc
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*17 SecNam
      Parameter (SecNam='LDF_CheckPSD_Full')

      Real*8 Thr, Tol
      Parameter (Thr=1.0d-10, Tol=1.0d-12)

      Logical  LDF_IntegralPrescreeningInfoIsSet, isSymmetric
      External LDF_IntegralPrescreeningInfoIsSet, isSymmetric
      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical IPI_set_here, doDiagonalization

      Integer AB, CD
      Integer A, B
      Integer nAB, nCD
      Integer ip_Int, l_Int
      Integer ip_I, l_I
      Integer ip_Indx, l_Indx
      Integer ip_Stat, l_Stat
      Integer iShell
      Integer m, N

      Real*8 x

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      doDiagonalization=Diagonalization
      Call Cho_Head('Enter '//SecNam,'=',80,6)
      Write(6,'(A,I3)') 'Integral mode:',Mode
      If (doDiagonalization) Then
          Write(6,'(A)') 'Using diagonalization'
      Else
          Write(6,'(A)') 'Using Cholesky decomposition'
      End If
      Call xFlush(6)

      ! Prescreening info
      If (Mode.gt.0 .and. Mode.le.3) Then
         If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
            Call LDF_SetIntegralPrescreeningInfo()
            IPI_set_here=.True.
         Else
            IPI_set_here=.False.
         End If
      Else If (Mode.eq.0) Then
         IPI_set_here=.False.
      Else
         Call WarningMessage(2,SecNam//': illegal Mode')
         Call LDF_Quit(1)
         IPI_set_here=.False.
      End If

      ! Allocations
      N=nBas_Valence*(nBas_Valence+1)/2
      x=dble(N)**2
      l_Int=int(x)
      If (l_Int.lt.0) Then
         Call WarningMessage(2,SecNam//': integer overflow (?)')
         Call LDF_Quit(1)
      End If
      Call GetMem('PSDInt','Allo','Real',ip_Int,l_Int)

      l_I=0
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         l_I=max(l_I,LDF_nBas_Atom(A)*LDF_nBas_Atom(B))
      End Do
      l_I=2*l_I**2
      Call GetMem('PSD_I','Allo','Real',ip_I,l_I)

      l_Indx=nShell_Valence
      Call GetMem('PSDIndx','Allo','Inte',ip_Indx,l_Indx)
      m=0
      Do iShell=1,nShell_Valence
         iWork(ip_Indx-1+iShell)=m
         m=m+nBasSh(iShell)
      End Do

      ! Compute integrals
      Call Cho_dZero(Work(ip_Int),l_Int)
      If (Mode.eq.0) Then ! exact integrals
         Do AB=1,NumberOfAtomPairs
            nAB=LDF_nBas_Atom(AP_Atoms(1,AB))
     &         *LDF_nBas_Atom(AP_Atoms(2,AB))
            Do CD=1,AB-1
               nCD=LDF_nBas_Atom(AP_Atoms(1,CD))
     &            *LDF_nBas_Atom(AP_Atoms(2,CD))
               Call LDF_ComputeValenceIntegrals(AB,CD,l_I,Work(ip_I))
               Call LDF_CheckPSD_Full_D(AB,CD,iWork(ip_Indx),
     &                                  nAB,nCD,Work(ip_I),
     &                                  N,Work(ip_Int))
               Call Trnsps(nAB,nCD,Work(ip_I),Work(ip_I+nAB*nCD))
               Call LDF_CheckPSD_Full_D(CD,AB,iWork(ip_Indx),
     &                                  nCD,nAB,Work(ip_I+nAB*nCD),
     &                                  N,Work(ip_Int))
            End Do
            Call LDF_ComputeValenceIntegrals(AB,AB,l_I,Work(ip_I))
            Call LDF_CheckPSD_Full_D(AB,AB,iWork(ip_Indx),
     &                               nAB,nAB,Work(ip_I),
     &                               N,Work(ip_Int))
         End Do
      Else
         Do AB=1,NumberOfAtomPairs
            nAB=LDF_nBas_Atom(AP_Atoms(1,AB))
     &         *LDF_nBas_Atom(AP_Atoms(2,AB))
            Do CD=1,AB-1
               nCD=LDF_nBas_Atom(AP_Atoms(1,CD))
     &            *LDF_nBas_Atom(AP_Atoms(2,CD))
               If (UseExactIntegrals.eq.2) Then
                  Call LDF_ComputeValenceIntegrals(AB,CD,l_I,Work(ip_I))
               Else
                  Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,CD,
     &                                                  l_I,Work(ip_I))
               End If
               Call LDF_CheckPSD_Full_D(AB,CD,iWork(ip_Indx),
     &                                  nAB,nCD,Work(ip_I),
     &                                  N,Work(ip_Int))
               Call Trnsps(nAB,nCD,Work(ip_I),Work(ip_I+nAB*nCD))
               Call LDF_CheckPSD_Full_D(CD,AB,iWork(ip_Indx),
     &                                  nCD,nAB,Work(ip_I+nAB*nCD),
     &                                  N,Work(ip_Int))
            End Do
            If (UseExactIntegrals.eq.1) Then
               Call LDF_ComputeValenceIntegrals(AB,AB,l_I,Work(ip_I))
            Else
               Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,AB,
     &                                               l_I,Work(ip_I))
            End If
            Call LDF_CheckPSD_Full_D(AB,AB,iWork(ip_Indx),
     &                               nAB,nAB,Work(ip_I),
     &                               N,Work(ip_Int))
         End Do
      End If

      ! Deallocation
      Call GetMem('PSDIndx','Free','Inte',ip_Indx,l_Indx)
      Call GetMem('PSD_I','Free','Real',ip_I,l_I)

      ! Check symmetry
      If (.not.isSymmetric(Work(ip_Int),N,Tol)) Then
         Call WarningMessage(2,
     &                        SecNam//': integral matrix not symmetric')
         Write(6,'(A,1P,D20.10)') 'Tolerance for symmetry=',Tol
         Call LDF_Quit(1)
      End If

      ! Check PSD
      irc=-1
      If (doDiagonalization) Then ! diagonalization
         l_Stat=11
         Call GetMem('PSDStat','Allo','Real',ip_Stat,l_Stat)
         Call LDF_CheckPSD_Full_Diag(N,Work(ip_Int),Work(ip_Stat),irc)
         If (irc.lt.0) Then
            Call WarningMessage(2,
     &                           SecNam//': irc<0 from diagonalization')
            Call LDF_Quit(1)
         Else If (irc.gt.0) Then
            Call WarningMessage(0,
     &                         SecNam//': full integral matrix not PSD')
            Write(6,'(I10,A)') irc,' negative eigenvalues:'
            Do i=1,irc
               Write(6,'(A,1X,I10,1X,1P,D20.10)')
     &         'Negative eigenvalue no.',i,Work(ip_Int-1+i)
            End Do
            Call xFlush(6)
         Else
            Call WarningMessage(0,SecNam//': full integral matrix PSD')
         End If
         Write(6,'(A,10X,I10)')
     &   'Number of eigenvalues..................',N
         Write(6,'(A,10X,I10)')
     &   'Number of large negative eigenvalues...',int(Work(ip_Stat+8))
         Write(6,'(A,10X,I10)')
     &   'Number of slightly neg. eigenvalues....',int(Work(ip_Stat+9))
         Write(6,'(A,10X,I10)')
     &   'Number of zero or pos. eigenvalues.....',int(Work(ip_Stat+10))
         Write(6,'(A,1P,D20.10)')
     &   'Minimum eigenvalue.....................',Work(ip_Stat)
         Write(6,'(A,1P,D20.10)')
     &   'Maximum................................',Work(ip_Stat+1)
         Write(6,'(A,1P,D20.10)')
     &   'Sum....................................',Work(ip_Stat+2)
         Write(6,'(A,1P,D20.10)')
     &   'Norm...................................',Work(ip_Stat+3)
         Write(6,'(A,1P,D20.10)')
     &   'Average................................',Work(ip_Stat+4)
         Write(6,'(A,1P,D20.10)')
     &   'Standard deviation.....................',Work(ip_Stat+5)
         Write(6,'(A,1P,D20.10)')
     &   'Skewness...............................',Work(ip_Stat+6)
         Write(6,'(A,1P,D20.10)')
     &   'Kurtosis...............................',Work(ip_Stat+7)
         Call xFlush(6)
         Call GetMem('PSDStat','Free','Real',ip_Stat,l_Stat)
      Else ! Cholesky decomposition
         Call LDF_CheckPSD_Full_CD(N,Work(ip_Int),Thr,irc)
         If (irc.lt.0) Then
            Call WarningMessage(2,
     &                    SecNam//': irc<0 from Cholesky decomposition')
            Call LDF_Quit(1)
         Else If (irc.eq.1) Then
            Call WarningMessage(0,
     &                         SecNam//': full integral matrix not PSD')
            x=9.9d9
            Do i=1,N
               x=min(x,Work(ip_Int-1+N*(i-1)+i))
            End Do
            Write(6,'(A,1X,1P,D20.10)')
     &      'Smallest diagonal element in CD:',x
            Call xFlush(6)
         Else If (irc.gt.1) Then
            Call WarningMessage(2,
     &                    SecNam//': irc>1 from Cholesky decomposition')
            Call LDF_Quit(1)
         Else
            Call WarningMessage(0,SecNam//': full integral matrix PSD')
         End If
      End If
      If (irc.ne.0) irc=1

      ! Deallocation
      Call GetMem('PSDInt','Free','Real',ip_Int,l_Int)

      ! Unset prescreening info (if set here)
      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPSD_Full_D(AB,CD,iSB,
     &                               nAB,nCD,Block,
     &                               N,Full)
      Implicit None
      Integer AB, CD
      Integer iSB(*)
      Integer nAB, nCD
      Real*8  Block(nAB,nCD)
      Integer N
      Real*8  Full(N,N)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer C, D
      Integer nShellC, nShellD
      Integer iS, jS
      Integer ipi, ipj
      Integer iShell, jShell
      Integer iOffi, iOffj
      Integer j_
      Integer ij, ij0
      Integer iOffB

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms, iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)
      nShellC=LDF_nShell_Atom(C)
      nShellD=LDF_nShell_Atom(D)
      ipi=LDF_lShell_Atom(C)-1
      ipj=LDF_lShell_Atom(D)-1

      iOffB=0
      Do jS=1,nShellD
         jShell=iWork(ipj+jS)
         iOffj=iSB(jShell)
         Do iS=1,nShellC
            iShell=iWork(ipi+iS)
            iOffi=iSB(iShell)
            Do j=1,nBasSh(jShell)
               j_=iOffj+j
               ij0=iOffB+nBasSh(iShell)*(j-1)
               Do i=1,nBasSh(iShell)
                  ij=iTri(iOffi+i,j_)
                  Call LDF_Block2Full_Packed(AB,Full(1,ij),iSB,
     &                                       Block(1,ij0+i))
               End Do
            End Do
            iOffB=iOffB+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPSD_Full_Diag(N,X,Stat,irc)
      Implicit None
      Integer N
      Real*8  X(N,N)
      Real*8  Stat(11)
      Integer irc
#include "WrkSpc.fh"

      Character*22 SecNam
      Parameter (SecNam='LDF_CheckPSD_Full_Diag')

      Integer iDummy
      Parameter (iDummy=-987654)

      Real*8 Dummy
      Real*8 ThrNeg
      Parameter (Dummy=-9.87654d0)
      Parameter (ThrNeg=-1.0d-10)

      Real*8   dLAMCH_
      External dLAMCH_

      Integer nFound
      Integer ip_EigVal, l_EigVal
      Integer ip_EigVec, l_EigVec
      Integer ip, l
      Integer ip_W, l_W, l_WT
      Integer ip_iW, l_iW, l_iWT

      Real*8 AbsTol

      If (N.lt.1) Then
         irc=0
      Else
         AbsTol=dLAMCH_('Safe minimum')
         l_EigVal=N
         Call GetMem('EigVal','Allo','Real',ip_EigVal,l_EigVal)
         l_EigVec=N
         Call GetMem('EigVec','Allo','Real',ip_EigVec,l_EigVec)
         l=2*N
         Call GetMem('iSuppZ','Allo','Inte',ip,l)
         l_WT=1
         Call GetMem('WA','Allo','Real',ip_W,l_WT)
         l_iWT=1
         Call GetMem('iWA','Allo','Inte',ip_iW,l_iWT)
         call dsyevr_('N','A','L',N,X,N,Dummy,Dummy,iDummy,iDummy,
     &                AbsTol,nFound,Work(ip_EigVal),Work(ip_EigVec),1,
     &                iWork(ip),Work(ip_W),-1,iWork(ip_iW),-1,irc)
         If (irc.ne.0) Then
            Call WarningMessage(2,SecNam//': nonzero rc from dSYEVR[0]')
            Call LDF_Quit(1)
         End If
         l_W=int(Work(ip_W))
         l_iW=iWork(ip_iW)
         Call GetMem('iWA','Free','Inte',ip_iW,l_iWT)
         Call GetMem('WA','Free','Real',ip_W,l_WT)
         Call GetMem('WA','Allo','Real',ip_W,l_W)
         Call GetMem('iWA','Allo','Inte',ip_iW,l_iW)
         call dsyevr_('N','A','L',N,X,N,Dummy,Dummy,iDummy,iDummy,
     &                AbsTol,nFound,Work(ip_EigVal),Work(ip_EigVec),1,
     &                iWork(ip),Work(ip_W),l_W,iWork(ip_iW),l_iW,irc)
         If (irc.ne.0) Then
            Call WarningMessage(2,SecNam//': nonzero rc from dSYEVR')
            Call LDF_Quit(1)
         Else
            If (nFound.ne.N) Then
               Call WarningMessage(2,SecNam//': nFound != N')
               Call LDF_Quit(1)
            End If
C-tbp:
      write(6,*) 'ERI eigenvalues:'
      write(6,'(1P,6D20.12)') (work(ip_eigval+irc),irc=0,N-1)
            Call dCopy_(N,Work(ip_EigVal),1,X(1,1),1) ! save eigenvalues
            Call LDF_CheckPSD_Full_Stat(N,X(1,1),ThrNeg,Stat)
            irc=int(Stat(9))
         End If
         Call GetMem('iWA','Free','Inte',ip_iW,l_iW)
         Call GetMem('WA','Free','Real',ip_W,l_W)
         Call GetMem('iSuppZ','Free','Inte',ip,l)
         Call GetMem('EigVec','Free','Real',ip_EigVec,l_EigVec)
         Call GetMem('EigVal','Free','Real',ip_EigVal,l_EigVal)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPSD_Full_CD(N,X,Thr,irc)
      Implicit None
      Integer N
      Real*8  X(N,N)
      Real*8  Thr
      Integer irc
#include "WrkSpc.fh"

      Integer ip, l, nVec

      If (N.lt.1) Then
         irc=0
      Else
         l=N**2
         Call GetMem('PSDCDV','Allo','Real',ip,l)
         Call CD_InCore(X,N,Work(ip),N,nVec,Thr,irc)
         If (irc.eq.101) irc=1
         Call GetMem('PSDCDV','Free','Real',ip,l)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPSD_Full_Stat(N,X,ThrNeg,Stat)
C1 Stat=(min,
C2       max,
C3       sum,
C4       norm,
C5       average,
C6       std dev wrt average,
C7       skewness,
C8       kurtosis
C9       #x for which x < ThrNeg
C10      #x for which ThrNeg <= x < 0
C11      #x for which x >= 0
      Implicit None
      Integer N
      Real*8  X(N)
      Real*8  Stat(*)
      Real*8  ThrNeg

      Integer nStat
      Parameter (nStat=11)

      real*8 ddot_
      external ddot_

      Real*8 Nx, mu, m1, m2, m3, m4, y
      Real*8 Tol

      Integer i, n1, n2

      If (ThrNeg.ge.0.0d0) Then
         Tol=-ThrNeg
      Else
         Tol=ThrNeg
      End If

      If (N.gt.0) Then
         Nx=dble(N)
         ! (1)=min value
         Stat(1)=X(1)
         Do i=2,N
            Stat(1)=min(Stat(1),X(i))
         End Do
         ! (2)=max value
         Stat(2)=X(1)
         Do i=2,N
            Stat(2)=max(Stat(2),X(i))
         End Do
         ! (3)=sum
         Stat(3)=X(1)
         Do i=2,N
            Stat(3)=Stat(3)+X(i)
         End Do
         ! compute average and moments
         mu=Stat(3)/Nx
         m1=0.0d0
         m2=0.0d0
         m3=0.0d0
         m4=0.0d0
         Do i=1,N
            y=X(i)-mu
            m1=m1+y
            m2=m2+y**2
            m3=m3+y**3
            m4=m4+y**4
         End Do
         m1=m1/Nx
         m2=m2/Nx
         m3=m3/Nx
         m4=m4/Nx
         ! (4)=norm
         Stat(4)=sqrt(dDot_(N,X,1,X,1))
         ! (5)=average
         Stat(5)=mu
         ! (6)=std dev (wrt average)
         Stat(6)=sqrt(m2)
         ! (7)=skewness
         If (m2.ne.0.0d0) Then
            Stat(7)=m3/sqrt(m2**3)
         Else
            Stat(7)=0.0d0
         End If
         ! (8)=kurtosis
         If (m2.ne.0.0d0) Then
            Stat(8)=m4/(m2**2)-3.0d0
         Else
            Stat(8)=0.0d0
         End If
         ! (9)=#x for which x < ThrNeg
         ! (10)=#x for which ThrNeg <= x < 0
         ! (11)=#x for which x >= 0
         n1=0
         n2=0
         Do i=1,N
            If (sign(1.0d0,X(i)).lt.0.0d0) Then
               If (X(i).lt.Tol) Then
                  n1=n1+1
               Else
                  n2=n2+1
               End If
            End If
         End Do
         Stat(9)=dble(n1)
         Stat(10)=dble(n2)
         Stat(11)=dble(N-n1-n2)
      Else
         Call dZero(Stat,nStat)
      End If

      End
