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
      Subroutine LDF_VerifyFit_Drv(irc)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: check fit by verifying that fitting equations are
C              are fulfilled for each atom pair.
C
      Implicit None
      Integer irc
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"

      Character*17 SecNam
      Parameter (SecNam='LDF_VerifyFit_Drv')

      Real*8 RMSTol
      Parameter (RMSTol=1.0d-10)

      Logical  LDF_ConstraintInfoIsSet
      External LDF_ConstraintInfoIsSet

      Integer  iPrintLevel, LDF_nBas_Atom, LDF_nBasAux_Pair_wLD
      Integer  LDF_nBasAux_Pair
      External iPrintLevel, LDF_nBas_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBasAux_Pair

      Logical LinDepRemoved
      Logical Silent
      Logical ConstraintInfoSetHere

      Integer AB
      Integer ip_C, l_C, l

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (NumberOfAtomPairs.lt.1) Then
         irc=0
         Return
      End If

      If (LDF_Constraint.lt.-1 .or. LDF_Constraint.gt.0) Then
         Call WarningMessage(2,SecNam//': unknown constraint')
         Write(6,'(A,I10)') 'Constraint=',LDF_Constraint
         Call LDF_Quit(1)
      End If

      If (LDF_ConstraintInfoIsSet(LDF_Constraint)) Then
         ConstraintInfoSetHere=.False.
      Else
         Call LDF_SetConstraint(LDF_Constraint)
         ConstraintInfoSetHere=.True.
      End If

      l_C=LDF_nBas_Atom(AP_Atoms(1,1))
     &   *LDF_nBas_Atom(AP_Atoms(2,1))
     &   *LDF_nBasAux_Pair_wLD(1)
      Do AB=2,NumberOfAtomPairs
         l_C=max(l_C,LDF_nBas_Atom(AP_Atoms(1,AB))
     &              *LDF_nBas_Atom(AP_Atoms(2,AB))
     &              *LDF_nBasAux_Pair_wLD(AB))
      End Do
      Call GetMem('VFC','Allo','Real',ip_C,l_C)

      LinDepRemoved=.False.
      Silent=iPrintLevel(-1).lt.3
      irc=0
      AB=0
      Do While (AB.lt.NumberOfAtomPairs .and. irc.eq.0)
         AB=AB+1
         If (LDF_Constraint.eq.0) Then ! get lambda and overlap
            l=LDF_nBas_Atom(AP_Atoms(1,AB))
     &       *LDF_nBas_Atom(AP_Atoms(2,AB))
     &       *LDF_nBasAux_Pair(AB)
            Call LDF_ReadUnconstrainedCoefficients(AB,l,Work(ip_C),irc)
            If (irc.eq.-1) Then
               Call WarningMessage(2,
     &         SecNam//': unconstrained coefficients not found on disk')
               Call LDF_Quit(1)
            Else If (irc.ne.0) Then
               Call WarningMessage(2,SecNam//
     &  ': non-zero return code from LDF_ReadUnconstrainedCoefficients')
               Write(6,'(A,I10)') 'irc=',irc
               Call LDF_Quit(1)
            End If
            Call LDF_AddChargeConstraintCorrection(AB,l,Work(ip_C))
         End If
         l=LDF_nBas_Atom(AP_Atoms(1,AB))
     &    *LDF_nBas_Atom(AP_Atoms(2,AB))
     &    *LDF_nBasAux_Pair_wLD(AB)
         Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l)
         Call LDF_VerifyFit(LinDepRemoved,Silent,LDF_Constraint,
     &                      RMSTol,AB,l,Work(ip_C),irc)
         If (irc.ne.0) Then
            Write(6,'(A,A,I10)')
     &      SecNam,': LDF_VerifyFit returned code',irc
            Write(6,'(A)')
     &      'Parameters passed to LDF_VerifyFit:'
            Write(6,'(3X,A,L1)')
     &      'LinDepRemoved=',LinDepRemoved
            Write(6,'(3X,A,L1)')
     &      'Silent=',Silent
            Write(6,'(3X,A,1P,D20.10)')
     &      'RMSTol=',RMSTol
            Write(6,'(3X,A,I10)')
     &      'AB=',AB
            Write(6,'(3X,A,I10)')
     &      'l=',l
         End If
      End Do

      Call GetMem('VFC','Free','Real',ip_C,l_C)

      If (ConstraintInfoSetHere) Then
         Call LDF_UnsetConstraint(LDF_Constraint)
      End If

      End
      Subroutine LDF_VerifyFit(LinDepRemoved,Silent,Constraint,
     &                         RMSTol,AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: Check fit by verifying that fitting equations are
C              fulfilled for atom pair AB.
C
      Implicit None
      Logical LinDepRemoved
      Logical Silent
      Integer Constraint
      Real*8  RMSTol
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc

      If (LinDepRemoved) Then
         Call LDF_VerifyFit_1(Silent,Constraint,RMSTol,AB,l_C,C,irc)
      Else
         Call LDF_VerifyFit_2(Silent,Constraint,RMSTol,AB,l_C,C,irc)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_VerifyFit_1(Silent,Constraint,RMSTol,AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: Check fit by verifying that fitting equations are
C              fulfilled for atom pair AB.
C              Linearly dependent functions removed from C.
C
      Implicit None
      Logical Silent
      Integer Constraint
      Real*8  RMSTol
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_charge_constraint_info.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_VerifyFit_1')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Real*8   dDot_, LDF_AtomicDistance, Cho_dSumElm
      external ddot_, LDF_AtomicDistance, Cho_dSumElm

      Character*4 Label1, Label2

      Integer n, M
      Integer ip_Int, l_Int
      Integer ip_G, l_G
      Integer ip_Stat, l_Stat
      Integer ip_JInt, l_JInt

      Real*8  NrmI, Nrm, RMS
      Real*8  SmI, Sm

      Integer i, j
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Set dimensions
      n=LDF_nBas_Atom(AP_Atoms(1,AB))*LDF_nBas_Atom(AP_Atoms(2,AB))
      M=LDF_nBasAux_Pair(AB)
      If (n.lt.1 .or. M.lt.1) Then
         irc=0
         Return
      End If
      If (l_C.lt.n*M) Then
         irc=-1
         Return
      End If

      ! Set index arrays for G matrix
      Call LDF_SetIndxG(AB)

      ! Allocate and compute 3-index integrals
      ! Add constraint term if needed
      l_Int=n*M
      Call GetMem('VFInt','Allo','Real',ip_Int,l_Int)
      Call LDF_ComputeIntegrals_uvJ(AB,l_Int,Work(ip_Int))
      If (Constraint.eq.0) Then ! charge constraint
         l_JInt=M
         Call GetMem('VFJInt','Allo','Real',ip_JInt,l_JInt)
         ! Get n(J)=\int J(r)dr for this pair with 1C LinDep removed
         ! and including 2C functions (copied from overlap)
         Call LDF_CC_GetRHS(AB,n,Work(ip_CC_Overlap),M,Work(ip_JInt))
         ! Add correction to integrals: (uv|J) + lambda(uv)*n(J)
         Call dGeR(n,M,1.0d0,Work(ip_CC_lambda),1,Work(ip_JInt),1,
     &             Work(ip_Int),n)
         Call GetMem('VFJInt','Free','Real',ip_JInt,l_JInt)
      Else If (Constraint.ne.-1) Then ! unknown constraint
         Call WarningMessage(2,SecNam//': unknown constraint')
         Write(6,'(A,I10)') 'Constraint=',Constraint
         Call LDF_Quit(1)
      End If
      NrmI=sqrt(dDot_(l_Int,Work(ip_Int),1,Work(ip_Int),1))
      SmI=Cho_dSumElm(Work(ip_Int),l_Int)

      ! Allocate and compute G matrix
      l_G=M**2
      Call GetMem('VFG','Allo','Real',ip_G,l_G)
      Call LDF_ComputeGMat(AB,M,Work(ip_G))

      ! Compute difference
      Call dGeMM_('N','N',n,M,M,
     &            -1.0d0,C,n,Work(ip_G),M,
     &            1.0d0,Work(ip_Int),n)

      ! Compute difference norm and evaluate irc
      Nrm=dDot_(l_Int,Work(ip_Int),1,Work(ip_Int),1)
      RMS=sqrt(Nrm/dble(l_Int))
      Nrm=sqrt(Nrm)
      If (RMS.gt.RMSTol) Then
         irc=1
      Else
         irc=0
      End If
      Sm=Cho_dSumElm(Work(ip_Int),l_Int)

      ! Compute and print statistics
      If (.not.Silent) Then
         Call LDF_SetAtomicLabels()
         Call LDF_GetAtomicLabel(AP_Atoms(1,AB),Label1)
         Call LDF_GetAtomicLabel(AP_Atoms(2,AB),Label2)
         l_Stat=7
         Call GetMem('VFStat','Allo','Real',ip_Stat,l_Stat)
         Call Statistics(Work(ip_Int),l_Int,Work(ip_Stat),1,2,3,4,5,6,7)
         Call Cho_Head(SecNam//': fit verification info','-',80,6)
         Write(6,'(2X,A,10X,I10,2X,A,2I10,2X,A,1X,A)')
     &   'Atom pair...........',AB,
     &   'Atoms...............',AP_Atoms(1,AB),AP_Atoms(2,AB),
     &   Label1,Label2
         Write(6,'(2X,A,10X,I10,2X,A,1P,D20.10)')
     &   'Auxiliary basis dim.',LDF_nBasAux_Pair(AB),
     &   'Atomic distance.....',LDF_AtomicDistance(AP_Atoms(1,AB),
     &                                             AP_Atoms(2,AB))
         Write(6,'(2X,A,I10,A,2X,A,10X,I10)')
     &   '1C LinDep...........',AP_1CLinDep(1,AB),'   (Excl.)',
     &   '2C Functions........',AP_2CFunctions(1,AB)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Integral norm.......',NrmI,
     &   'Difference norm.....',Nrm
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Integral sum........',SmI,
     &   'Difference sum......',Sm
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Average.............',Work(ip_Stat),
     &   'Abs Average.........',Work(ip_Stat+1)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Max Difference......',Work(ip_Stat+3),
     &   'Max Abs Difference..',Work(ip_Stat+4)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Variance............',Work(ip_Stat+5),
     &   'Unbiased Variance...',Work(ip_Stat+6)
         Write(6,'(2X,A,1P,D20.10)')
     &   'RMS.................',RMS
         Call xFlush(6)
         Call GetMem('VFStat','Free','Real',ip_Stat,l_Stat)
         Call LDF_UnsetAtomicLabels()
      End If

      ! Deallocation
      Call GetMem('VFG','Free','Real',ip_G,l_G)
      Call GetMem('VFInt','Free','Real',ip_Int,l_Int)

      ! Unset index arrays for G matrix
      Call LDF_UnsetIndxG()

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_VerifyFit_2(Silent,Constraint,RMSTol,AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: Check fit by verifying that fitting equations are
C              fulfilled for atom pair AB.
C              Linearly dependent functions included in C.
C
      Implicit None
      Logical Silent
      Integer Constraint
      Real*8  RMSTol
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_charge_constraint_info.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_VerifyFit_2')

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet

      Integer  LDF_nBasAux_Pair
      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Atom
      External LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Atom

      Real*8   dDot_, LDF_AtomicDistance, Cho_dSumElm
      external ddot_, LDF_AtomicDistance, Cho_dSumElm

      Character*4 Label1, Label2

      Logical IPI_set_here

      Integer A, B
      Integer n, M, MA, MB, MAB
      Integer ip_Int, l_Int
      Integer ip_G, l_G
      Integer ip_Stat, l_Stat
      Integer ipI, ipJ, l
      Integer ip_2CInt, l_2CInt
      Integer nRow_Map, nCol_Map
      Integer ip_Map, l_Map
      Integer ip2, ipS
      Integer i2C

      Real*8  NrmI, Nrm, RMS
      Real*8  SmI, Sm

      Real*8 tau
      Parameter (tau=0.0d0)

      Integer i, j
      Integer Map
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      Map(i)=iWork(ip_Map-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)

      ! Set prescreening info (if not already done)
      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
         Call LDF_SetIntegralPrescreeningInfo()
         IPI_set_here=.True.
      Else
         IPI_set_here=.False.
      End If

      ! Get atoms
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)

      ! Set and check dimensions
      n=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      M=LDF_nBasAux_Pair_wLD(AB)
      MA=LDF_nBasAux_Atom(A)
      MB=LDF_nBasAux_Atom(B)
      MAB=AP_2CFunctions(1,AB)
      If (n.lt.1 .or. M.lt.1) Then
         irc=0
         Return
      End If
      If (l_C.lt.n*M) Then
         irc=-1
         Return
      End If

      ! Allocate G matrix blocks
      l_G=max(MA,MB,MAB)*max(MA,MB,MAB)
      Call GetMem('VFG','Allo','Real',ip_G,l_G)

      ! Allocate 3-index integrals
      l_Int=n*M
      Call GetMem('VFInt','Allo','Real',ip_Int,l_Int)

      ! Compute 3-index integrals (u_A v_B | J_A;J_B;J_AB)
      l=n*MA
      ipI=ip_Int
      Call LDF_Compute3IndexIntegrals_1(AB,A,tau,l,Work(ipI))
      ipI=ipI+l
      If (B.ne.A) Then
         l=n*MB
         Call LDF_Compute3IndexIntegrals_1(AB,B,tau,l,Work(ipI))
         ipI=ipI+l
      End If
      If (MAB.gt.0) Then
         l=n*MAB
         Call LDF_Compute3IndexIntegrals_2(AB,AB,tau,l,Work(ipI))
      End If
      ! Add constraint term if needed
      If (Constraint.eq.0) Then ! charge constraint
         ! (uv|J) + lambda(uv)*n(J), n(J)=\int J(r)dr
         ipI=ip_Int
         ipJ=iWork(ip_CC_AuxIntVec-1+A)
         Call dGeR(n,MA,1.0d0,Work(ip_CC_lambda),1,Work(ipJ),1,
     &             Work(ipI),n)
         ipI=ipI+n*MA
         If (B.ne.A) Then
            ipJ=iWork(ip_CC_AuxIntVec-1+B)
            Call dGeR(n,MB,1.0d0,Work(ip_CC_lambda),1,Work(ipJ),1,
     &                Work(ipI),n)
            ipI=ipI+n*MB
         End If
         If (MAB.gt.0) Then
            l_2CInt=MAB
            Call GetMem('VF2CInt','Allo','Real',ip_2CInt,l_2CInt)
            nRow_Map=MAB
            nCol_Map=1
            l_Map=nRow_Map*nCol_Map
            Call GetMem('VFMap','Allo','Inte',ip_Map,l_Map)
            Call LDF_Map2CF(AB,nRow_Map,nCol_Map,iWork(ip_Map))
            ip2=ip_2CInt-1
            ipS=ip_CC_Overlap-1
            Do i2C=1,MAB
               Work(ip2+i2C)=Work(ipS+Map(i2C))
            End Do
            Call dGeR(n,MAB,1.0d0,Work(ip_CC_lambda),1,Work(ip_2CInt),1,
     &                Work(ipI),n)
            Call GetMem('VFMap','Free','Inte',ip_Map,l_Map)
            Call GetMem('VF2CInt','Free','Real',ip_2CInt,l_2CInt)
         End If
      Else If (Constraint.ne.-1) Then ! unknown constraint
         Call WarningMessage(2,SecNam//': unknown constraint')
         Write(6,'(A,I10)') 'Constraint=',Constraint
         Call LDF_Quit(1)
      End If
      NrmI=sqrt(dDot_(l_Int,Work(ip_Int),1,Work(ip_Int),1))
      SmI=Cho_dSumElm(Work(ip_Int),l_Int)

      ! Compute contributions from 2-index integrals (J_A|K_A)
      l=MA**2
      Call LDF_Compute2IndexIntegrals_11(A,A,tau,l,Work(ip_G))
      Call dGeMM_('N','N',n,MA,MA,
     &            -1.0d0,C,n,Work(ip_G),max(MA,1),
     &            1.0d0,Work(ip_Int),n)

      ! Compute contributions from 2-index integrals (J_A|K_B)
      ! and (J_B|K_B)
      If (B.ne.A) Then
         l=MA*MB
         Call LDF_Compute2IndexIntegrals_11(A,B,tau,l,Work(ip_G))
         Call dGeMM_('N','N',n,MB,MA,
     &               -1.0d0,C,n,Work(ip_G),max(MA,1),
     &               1.0d0,Work(ip_Int+n*MA),n)
         Call dGeMM_('N','T',n,MA,MB,
     &               -1.0d0,C(n*MA+1),n,Work(ip_G),max(MA,1),
     &               1.0d0,Work(ip_Int),n)
         l=MB**2
         Call LDF_Compute2IndexIntegrals_11(B,B,tau,l,Work(ip_G))
         Call dGeMM_('N','N',n,MB,MB,
     &               -1.0d0,C(n*MA+1),n,Work(ip_G),max(MB,1),
     &               1.0d0,Work(ip_Int+n*MA),n)
      End If

      ! Compute contributions from 2-index integrals (J_A|K_AB)
      ! and (J_B|K_AB) and (J_AB|K_AB)
      If (MAB.gt.0) Then
         l=MA*MAB
         Call LDF_Compute2IndexIntegrals_12(A,AB,tau,l,Work(ip_G))
         If (B.eq.A) Then
            Call dGeMM_('N','N',n,MAB,MA,
     &                  -1.0d0,C,n,Work(ip_G),max(MA,1),
     &                  1.0d0,Work(ip_Int+n*MA),n)
            Call dGeMM_('N','T',n,MA,MAB,
     &                  -1.0d0,C(n*MA+1),n,Work(ip_G),max(MA,1),
     &                  1.0d0,Work(ip_Int),n)
         Else
            Call dGeMM_('N','N',n,MAB,MA,
     &                  -1.0d0,C,n,Work(ip_G),max(MA,1),
     &                  1.0d0,Work(ip_Int+n*(MA+MB)),n)
            Call dGeMM_('N','T',n,MA,MAB,
     &                  -1.0d0,C(n*(MA+MB)+1),n,Work(ip_G),max(MA,1),
     &                  1.0d0,Work(ip_Int),n)
         End If
         If (B.ne.A) Then
            l=MB*MAB
            Call LDF_Compute2IndexIntegrals_12(B,AB,tau,l,Work(ip_G))
            Call dGeMM_('N','N',n,MAB,MB,
     &                  -1.0d0,C(n*MA+1),n,Work(ip_G),max(MB,1),
     &                  1.0d0,Work(ip_Int+n*(MA+MB)),n)
            Call dGeMM_('N','T',n,MB,MAB,
     &                  -1.0d0,C(n*(MA+MB)+1),n,Work(ip_G),max(MB,1),
     &                  1.0d0,Work(ip_Int+n*MA),n)
         End If
         l=MAB**2
         Call LDF_Compute2IndexIntegrals_22(AB,AB,tau,l,Work(ip_G))
         If (B.eq.A) Then
            Call dGeMM_('N','N',n,MAB,MAB,
     &                  -1.0d0,C(n*MA+1),n,Work(ip_G),MAB,
     &                  1.0d0,Work(ip_Int+n*MA),n)
         Else
            Call dGeMM_('N','N',n,MAB,MAB,
     &                  -1.0d0,C(n*(MA+MB)+1),n,Work(ip_G),MAB,
     &                  1.0d0,Work(ip_Int+n*(MA+MB)),n)
         End If
      End If

      ! Compute difference norm and evaluate irc
      Nrm=dDot_(l_Int,Work(ip_Int),1,Work(ip_Int),1)
      RMS=sqrt(Nrm/dble(l_Int))
      Nrm=sqrt(Nrm)
      If (RMS.gt.RMSTol) Then
         irc=1
      Else
         irc=0
      End If
      Sm=Cho_dSumElm(Work(ip_Int),l_Int)

      ! Compute and print statistics
      If (.not.Silent) Then
         Call LDF_SetAtomicLabels()
         Call LDF_GetAtomicLabel(AP_Atoms(1,AB),Label1)
         Call LDF_GetAtomicLabel(AP_Atoms(2,AB),Label2)
         l_Stat=7
         Call GetMem('VFStat','Allo','Real',ip_Stat,l_Stat)
         Call Statistics(Work(ip_Int),l_Int,Work(ip_Stat),1,2,3,4,5,6,7)
         Call Cho_Head(SecNam//': fit verification info','-',80,6)
         Write(6,'(2X,A,10X,I10,2X,A,2I10,2X,A,1X,A)')
     &   'Atom pair...........',AB,
     &   'Atoms...............',AP_Atoms(1,AB),AP_Atoms(2,AB),
     &   Label1,Label2
         Write(6,'(2X,A,10X,I10,2X,A,1P,D20.10)')
     &   'Auxiliary basis dim.',LDF_nBasAux_Pair(AB),
     &   'Atomic distance.....',LDF_AtomicDistance(AP_Atoms(1,AB),
     &                                             AP_Atoms(2,AB))
         Write(6,'(2X,A,I10,A,2X,A,10X,I10)')
     &   '1C LinDep...........',AP_1CLinDep(1,AB),'   (Incl.)',
     &   '2C Functions........',AP_2CFunctions(1,AB)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Integral norm.......',NrmI,
     &   'Difference norm.....',Nrm
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Integral sum........',SmI,
     &   'Difference sum......',Sm
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Average.............',Work(ip_Stat),
     &   'Abs Average.........',Work(ip_Stat+1)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Max Difference......',Work(ip_Stat+3),
     &   'Max Abs Difference..',Work(ip_Stat+4)
         Write(6,'(2X,A,1P,D20.10,2X,A,D20.10)')
     &   'Variance............',Work(ip_Stat+5),
     &   'Unbiased Variance...',Work(ip_Stat+6)
         Write(6,'(2X,A,1P,D20.10)')
     &   'RMS.................',RMS
         Call xFlush(6)
         Call GetMem('VFStat','Free','Real',ip_Stat,l_Stat)
         Call LDF_UnsetAtomicLabels()
      End If

      ! Deallocation
      Call GetMem('VFInt','Free','Real',ip_Int,l_Int)
      Call GetMem('VFG','Free','Real',ip_G,l_G)

      ! Unset integral prescreening info (if done here)
      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      End
