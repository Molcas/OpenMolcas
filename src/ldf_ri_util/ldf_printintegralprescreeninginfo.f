************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine LDF_PrintIntegralPrescreeningInfo(tau)
      Implicit None
      Real*8  tau
#include "WrkSpc.fh"
#include "localdf.fh"
#include "localdf_int.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_atom_pair_info.fh"

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet

      Integer  LDF_nAtom, LDF_nShell_Atom, LDF_nAuxShell_Atom
      Integer  LDF_GlobalToAtomicShell
      External LDF_nAtom, LDF_nShell_Atom, LDF_nAuxShell_Atom
      External LDF_GlobalToAtomicShell

      Real*8 tau2
      Real*8 Total, Calc

      Integer A, B, AB, CD
      Integer nSA, nSB, iSA, iSB, ip
      Integer nnShl, ijS, iShell, jShell

      Integer i, j
      Real*8  Imax
      Real*8  Gmax_1C
      Real*8  Gmax_2C
      Integer ip_Imax
      Integer ip_Gmax_1C
      Integer ip_Gmax_2C
      Integer AP_Atoms
      Integer AP_2CFunctions
      Imax(i)=Work(ip_IDiag_Mx-1+i)
      Gmax_1C(i)=Work(ip_GDiag_1C_Mx-1+i)
      Gmax_2C(i)=Work(ip_GDiag_2C_Mx-1+i)
      ip_Imax(i)=iWork(ip_IDiag+2*(i-1)+1)
      ip_Gmax_1C(i)=iWork(ip_GDiag_1C+2*(i-1)+1)
      ip_Gmax_2C(i)=iWork(ip_GDiag_2C+2*(i-1)+1)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Return

      Call Cho_Head('Max integral diagonals, atom pair based','=',80,6)
      Write(6,'(/,A)')
     & 'Atom Pair      max(uv|uv)'
      Write(6,'(A)')
     & '-------------------------'
      Do AB=1,NumberOfAtomPairs
         Write(6,'(I9,1X,1P,D15.6)') AB,Imax(AB)
      End Do
      Write(6,'(A)')
     & '-------------------------'

      Call Cho_Head('Max integral diagonals, shell based','=',80,6)
      Write(6,'(/,A)')
     & '   Atom A    Atom B   Shell A   Shell B      max(uv|uv)'
      Write(6,'(A)')
     & '-------------------------------------------------------'
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nSA=LDF_nShell_Atom(A)
         nSB=LDF_nShell_Atom(B)
         ip=ip_Imax(AB)-1
         Do iSB=1,nSB
            Do iSA=1,nSA
               Write(6,'(4(I9,1X),1P,D15.6)')
     &         A,B,iSA,iSB,Work(ip+nSA*(iSB-1)+iSA)
            End Do
         End Do
      End Do
      Write(6,'(A)')
     & '-------------------------------------------------------'

      Call Cho_Head('Max 1C G diagonals, atom based','=',80,6)
      Write(6,'(/,A)')
     & '     Atom    max(J_A|J_A)'
      Write(6,'(A)')
     & '-------------------------'
      Do A=1,LDF_nAtom()
         Write(6,'(I9,1X,1P,D15.6)') A,Gmax_1C(A)
      End Do
      Write(6,'(A)')
     & '-------------------------'

      Call Cho_Head('Max 1C G diagonals, shell based','=',80,6)
      Write(6,'(/,A)')
     & '     Atom  AuxShell    max(J_A|J_A)'
      Write(6,'(A)')
     & '-----------------------------------'
      Do A=1,LDF_nAtom()
         nSA=LDF_nAuxShell_Atom(A)
         ip=ip_Gmax_1C(A)-1
         Do iSA=1,nSA
            Write(6,'(2(I9,1X),1P,D15.6)')
     &      A,iSA,Work(ip+iSA)
         End Do
      End Do
      Write(6,'(A)')
     & '-----------------------------------'

      If (LDF2) Then
         Call Cho_Head('Max 2C G diagonals, atom pair based','=',80,6)
         Write(6,'(/,A)')
     &   'Atom Pair  max(J_AB|J_AB)'
         Write(6,'(A)')
     &   '-------------------------'
         Do AB=1,NumberOfAtomPairs
            If (AP_2CFunctions(1,AB).gt.0) Then
               Write(6,'(I9,1X,1P,D15.6)') AB,Gmax_2C(AB)
            End If
         End Do
         Write(6,'(A)')
     &   '-------------------------'
         Call Cho_Head('Max 2C G diagonals, shell based','=',80,6)
         Write(6,'(/,A)')
     &   '   Atom A    Atom B   Shell A   Shell B  max(J_AB|J_AB)'
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Do AB=1,NumberOfAtomPairs
            If (AP_2CFunctions(1,AB).gt.0) Then
               Call LDF_SetIndxG(AB)
               nnShl=iWork(ip_GDiag_2C+2*(AB-1))
               If (nnShl.ne.l_2CList_2) Then
                  Call WarningMessage(2,
     &          'LDF_PrintIntegraPrescreeningInfo: nnShl != l_2CList_2')
                  Call LDF_Quit(1)
               End If
               ip=ip_Gmax_2C(AB)
               A=AP_Atoms(1,AB)
               B=AP_Atoms(2,AB)
               Do ijS=0,nnShl-1
                  iShell=iWork(ip_2CList+3*ijS)
                  iSA=LDF_GlobalToAtomicShell(A,iShell)
                  jShell=iWork(ip_2CList+3*ijS+1)
                  iSB=LDF_GlobalToAtomicShell(B,jShell)
                  Write(6,'(4(I9,1X),1P,D15.6)')
     &            A,B,iSA,iSB,Work(ip+ijS)
               End Do
               Call LDF_UnsetIndxG()
            End If
         End Do
         Write(6,'(A)')
     &   '-------------------------------------------------------'
      End If

      tau2=tau**2
      Call Cho_Head('Atom/Pair block screening statistics','=',80,6)
      Write(6,'(/,A,1P,D15.6)') 'Screening threshold: ',tau
      Calc=0.0d0
      Do AB=1,NumberOfAtomPairs
         Do A=1,LDF_nAtom()
            If (Imax(AB)*Gmax_1C(A).gt.tau2) Then
               Calc=Calc+1.0d0
            End If
         End Do
      End Do
      Total=dble(NumberOfAtomPairs)*dble(LDF_nAtom())
      Write(6,'(A,F7.2,A)')
     & 'Screening percentage for (u_A v_B|J_C): ',
     & 1.0d2*(Total-Calc)/Total,' %'
      If (LDF2) Then
         Total=0.0d0
         Calc=0.0d0
         Do AB=1,NumberOfAtomPairs
            Do CD=1,NumberOfAtomPairs
               If (AP_2CFunctions(1,CD).gt.0) Then
                  Total=Total+1.0d0
                  If (Imax(AB)*Gmax_2C(CD).gt.tau2) Then
                     Calc=Calc+1.0d0
                  End If
               End If
            End Do
         End Do
         If (Total.gt.0.0d0) Then
            Write(6,'(A,F7.2,A)')
     &      'Screening percentage for (u_A v_B|J_CD):',
     &      1.0d2*(Total-Calc)/Total,' %'
         End If
      End If
      Calc=0.0d0
      Do A=1,LDF_nAtom()
         Do B=1,A
            If (Gmax_1C(A)*Gmax_1C(B).gt.tau2) Then
               Calc=Calc+1.0d0
            End If
         End Do
      End Do
      Total=dble(LDF_nAtom())*(dble(LDF_nAtom())+1.0d0)/2.0d0
      Write(6,'(A,F7.2,A)')
     & 'Screening percentage for (J_A|J_B):     ',
     & 1.0d2*(Total-Calc)/Total,' %'
      If (LDF2) Then
         Total=0.0d0
         Calc=0.0d0
         Do A=1,LDF_nAtom()
            Do CD=1,NumberOfAtomPairs
               If (AP_2CFunctions(1,CD).gt.0) Then
                  Total=Total+1.0d0
                  If (Gmax_1C(A)*Gmax_2C(CD).gt.tau2) Then
                     Calc=Calc+1.0d0
                  End If
               End If
            End Do
         End Do
         If (Total.gt.0.0d0) Then
            Write(6,'(A,F7.2,A)')
     &      'Screening percentage for (J_A|J_CD):    ',
     &      1.0d2*(Total-Calc)/Total,' %'
         End If
         Total=0.0d0
         Calc=0.0d0
         Do AB=1,NumberOfAtomPairs
            If (AP_2CFunctions(1,AB).gt.0) Then
               Do CD=1,AB
                  If (AP_2CFunctions(1,CD).gt.0) Then
                     Total=Total+1.0d0
                     If (Gmax_2C(AB)*Gmax_2C(CD).gt.tau2) Then
                        Calc=Calc+1.0d0
                     End If
                  End If
               End Do
            End If
         End Do
         If (Total.gt.0.0d0) Then
            Write(6,'(A,F7.2,A)')
     &      'Screening percentage for (J_AB|J_CD):   ',
     &      1.0d2*(Total-Calc)/Total,' %'
         End If
      End If

      Call xFlush(6)

      End
