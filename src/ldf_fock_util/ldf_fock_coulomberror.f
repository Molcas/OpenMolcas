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
      Subroutine LDF_Fock_CoulombErrorAnalysis(ComputeF,
     &                                         Mode,PackedD,PackedF,
     &                                         nD,FactC,ip_D,ip_F)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: analyze Coulomb error
C
C              F(uv)-Ftilde(uv)= FactC * sum_kl { (uv|kl)
C                                                -[uv|kl] } * D(kl)
C
C              where [uv|kl] are the LDF integrals, and compare to
C              the upper bound.
C
C     If ComputF: the LDF Fock matrix is computed here
C     Else: on input, ip_F should point to the Fock matrix computed from
C     LDF integrals (replaced with the error on exit)!
C
      Implicit None
      Logical ComputeF
      Integer Mode
      Logical PackedD
      Logical PackedF
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_Fock_CoulombErrorAnalysis')

      real*8 ddot_
      external ddot_

      Logical PrintNorm
      Parameter (PrintNorm=.False.)

      Real*8 Stat(7,3)
      Real*8 RMS1, RMS2, RMS3

      Logical Add
      Logical Packed_myF

      Integer ip_myF, l_myF
      Integer ipF, lF
      Integer iD
      Integer i

      If (nD.lt.1) Return

      ! Compute error
      Call LDF_Fock_CoulombError(PrintNorm,ComputeF,
     &                           Mode,PackedD,PackedF,
     &                           nD,FactC,ip_D,ip_F)

      ! Compute upper bound
      Add=.False.
      Packed_myF=PackedF
      l_myF=nD
      Call GetMem('CEAmyFP','Allo','Inte',ip_myF,l_myF)
      If (Packed_myF) Then
         lF=nBas_Valence*(nBas_Valence+1)/2
      Else
         lF=nBas_Valence**2
      End If
      Do iD=1,nD
         Call GetMem('CEAmyF','Allo','Real',ipF,lF)
         iWork(ip_myF-1+iD)=ipF
      End Do
      Call LDF_Fock_CoulombUpperBound_Full(PrintNorm,
     &                                     Add,PackedD,Packed_myF,
     &                                     nD,FactC,ip_D,iWork(ip_myF))

      ! Analysis
      Call Cho_Head('Coulomb Error','-',80,6)
      Do iD=1,nD
         Call Statistics(Work(iWork(ip_myF-1+iD)),lF,Stat(1,1),
     &                   1,2,3,4,5,6,7)
         RMS1=dDot_(lF,Work(iWork(ip_myF-1+iD)),1,
     &                Work(iWork(ip_myF-1+iD)),1)
         Call Statistics(Work(ip_F(iD)),lF,Stat(1,2),1,2,3,4,5,6,7)
         RMS2=dDot_(lF,Work(ip_F(iD)),1,
     &                Work(ip_F(iD)),1)
         Do i=1,lF
            Work(iWork(ip_myF-1+iD)-1+i)=Work(iWork(ip_myF-1+iD)-1+i)
     &                                  -abs(Work(ip_F(iD)-1+i))
         End Do
         Call Statistics(Work(iWork(ip_myF-1+iD)),lF,Stat(1,3),
     &                   1,2,3,4,5,6,7)
         RMS3=dDot_(lF,Work(iWork(ip_myF-1+iD)),1,
     &                Work(iWork(ip_myF-1+iD)),1)
         Write(6,'(/,2X,A,I10,A)')
     &   'Coulomb error for density',iD,' (Upper bound,Actual,Diff):'
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Average error......',Stat(1,1),Stat(1,2),Stat(1,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Abs average error..',Stat(2,1),Stat(2,2),Stat(2,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Min error..........',Stat(3,1),Stat(3,2),Stat(3,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Max error..........',Stat(4,1),Stat(4,2),Stat(4,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Max abs error......',Stat(5,1),Stat(5,2),Stat(5,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Variance...........',Stat(6,1),Stat(6,2),Stat(6,3)
         Write(6,'(2X,A,1P,3D20.10)')
     &   'Norm...............',sqrt(RMS1),sqrt(RMS2),sqrt(RMS3)
         If (lF.gt.0) Then
            RMS1=sqrt(RMS1/dble(lF))
            RMS2=sqrt(RMS2/dble(lF))
            RMS3=sqrt(RMS3/dble(lF))
         Else
            RMS1=0.0d0
            RMS2=0.0d0
            RMS3=0.0d0
         End If
         Write(6,'(2X,A,1P,3D20.10)')
     &   'RMS error..........',RMS1,RMS2,RMS3
         Call xFlush(6)
         If ((Stat(5,1)-Stat(5,2)).lt.0.0d0) Then
            If (abs(Stat(5,1)-Stat(5,2)).gt.1.0d-6) Then
               Call WarningMessage(2,
     &           SecNam//': max abs error is greater than upper bound!')
               Call LDF_Quit(1)
            End If
         End If
      End Do

      ! Deallocations
      Do iD=1,nD
         ipF=iWork(ip_myF-1+iD)
         Call GetMem('CEAmyF','Free','Real',ipF,lF)
      End Do
      Call GetMem('CEAmyFP','Free','Inte',ip_myF,l_myF)

      End
      Subroutine LDF_Fock_CoulombError(PrintNorm,ComputeF,
     &                                 Mode,PackedD,PackedF,
     &                                 nD,FactC,ip_D,ip_F)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute Coulomb error
C
C              F(uv)-Ftilde(uv)= FactC * sum_kl { (uv|kl)
C                                                -[uv|kl] } * D(kl)
C
C              where [uv|kl] are the LDF integrals.
C
C     If ComputF: the LDF Fock matrix is computed here
C     Else: on input, ip_F should point to the Fock matrix computed from
C     LDF integrals (replaced with the error on exit)!
C
      Implicit None
      Logical PrintNorm
      Logical ComputeF
      Integer Mode
      Logical PackedD
      Logical PackedF
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      real*8 ddot_
      external ddot_

      Logical Timing
      Logical Add

      Integer IntegralOption
      Integer ip_myF, l_myF
      Integer ipF, lF
      Integer iD

      Real*8  ThrPS(2)

      If (ComputeF) Then
         IntegralOption=0
         Timing=.False.
         ThrPS(1)=0.0d0
         ThrPS(2)=0.0d0
         Add=.False.
         Call LDF_Fock_CoulombOnly(IntegralOption,
     &                             Timing,Mode,ThrPS,
     &                             Add,PackedD,PackedF,
     &                             nD,FactC,ip_D,ip_F)
      End If

      If (PackedF) Then
         lF=nBas_Valence*(nBas_Valence+1)/2
      Else
         lF=nBas_Valence**2
      End If
      l_myF=nD
      Call GetMem('myFPtr','Allo','Inte',ip_myF,l_myF)
      Do iD=1,nD
         Call GetMem('myF','Allo','Real',ipF,lF)
         iWork(ip_myF-1+iD)=ipF
      End Do
      IntegralOption=222 ! use conventional integrals
      Timing=.False.
      ThrPS(1)=0.0d0
      ThrPS(2)=0.0d0
      Add=.False.
      Call LDF_Fock_CoulombOnly(IntegralOption,
     &                          Timing,Mode,ThrPS,
     &                          Add,PackedD,PackedF,
     &                          nD,FactC,ip_D,iWork(ip_myF))
      Do iD=1,nD
         ipF=iWork(ip_myF-1+iD)
         Call dAXPY_(lF,-1.0d0,Work(ipF),1,Work(ip_F(iD)),1)
         Call dScal_(lF,-1.0d0,Work(ip_F(iD)),1)
      End Do
      Do iD=1,nD
         ipF=iWork(ip_myF-1+iD)
         Call GetMem('myF','Free','Real',ipF,lF)
      End Do
      Call GetMem('myFPtr','Free','Inte',ip_myF,l_myF)

      If (PrintNorm) Then
         Do iD=1,nD
            Write(6,'(A,I10,A,1P,D20.10)')
     &      'Norm of Coulomb error for density',iD,':',
     &      sqrt(dDot_(lF,Work(ip_F(iD)),1,Work(ip_F(iD)),1))
         End Do
         Call xFlush(6)
      End If

      End
