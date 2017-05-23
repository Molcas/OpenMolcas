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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_ComputeCoulombIntermediates(Timing,nD,ip_DBlocks,
     &                                                    ip_V,ip_CNorm)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Compute Coulomb intermediates
C
C              V(J) = sum_uv C(uv,J)*D(uv)
C
C              using LDF fitting coefficients.
C              If ip_CNorm>0 on entry, Frobenius norms of the fitting
C              coefficients are computed as a byproduct, stored as:
C              Work(ip_CNorm-1+4*(AB-1)+1)
C                  =sum_uAvBJ sqrt[C(uAvB,J)**2]       {all J}
C              Work(ip_CNorm-1+4*(AB-1)+2)
C                  =sum_uAvBJA sqrt[C(uAvB,JA)**2]     {J on A}
C              Work(ip_CNorm-1+4*(AB-1)+3)
C                  =sum_uAvBJB sqrt[C(uAvB,JB)**2]     {J on B}
C              Work(ip_CNorm-1+4*(AB-1)+4)
C                  =sum_uAvBJAB sqrt[C(uAvB,JAB)**2]   {J on AB (2CF)}
C
C     FOR BLOCKED VERSION: Call LDF_ComputeCoulombIntermediates0
C
      Implicit None
      Logical Timing
      Integer nD
      Integer ip_DBlocks(nD)
      Integer ip_V(nD)
      Integer ip_CNorm
#if defined (_MOLCAS_MPP_)
#include "para_info.fh"
#endif
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*31 SecNam
      Parameter (SecNam='LDF_ComputeCoulombIntermediates')

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer  LDF_nBas_Atom, LDF_nBasAux_Atom
      Integer  LDF_nBasAux_Pair_wLD, LDF_nAtom
      External LDF_nBas_Atom, LDF_nBasAux_Atom
      External LDF_nBasAux_Pair_wLD, LDF_nAtom

      real*8 ddot_
      external ddot_

      Logical doNorm

      Real*8 tC1, tC2
      Real*8 tW1, tW2

      Integer TaskListID
      Integer iD
      Integer ip_C, l_C
      Integer iAtomPair
      Integer iAtom, jAtom, nAtom
      Integer nuv, M
      Integer ipD, ipV, ipC

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      If (Timing) Call CWTIme(tC1,tW1)

      ! Initialize V arrays
      Do iD=1,nD
         Call LDF_ZeroAuxBasVector(ip_V(iD))
      End Do

      ! Allocate array for storing coefficients
      l_C=0
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nuv=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
         M=LDF_nBasAux_Pair_wLD(iAtomPair)
         l_C=max(l_C,nuv*M)
      End Do
      Call GetMem('LDFCBlk','Allo','Real',ip_C,l_C)

      ! compute norm?
      doNorm=ip_CNorm.gt.0
#if defined (_MOLCAS_MPP_)
      ! Init norm array
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         If (doNorm) Then
            Call Cho_dZero(Work(ip_CNorm),4*NumberOfAtomPairs)
         End If
      End If
#endif

      ! Compute V
      nAtom=LDF_nAtom()
      Call Init_Tsk(TaskListID,NumberOfAtomPairs)
      Do While (Rsv_Tsk(TaskListID,iAtomPair))
         Call LDF_CIO_ReadC_wLD(iAtomPair,Work(ip_C),l_C)
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nuv=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
         ipC=ip_C
         M=LDF_nBasAux_Atom(iAtom)
         If (doNorm) Then
            Work(ip_CNorm+4*(iAtomPair-1))
     &                  =sqrt(dDot_(nuv*LDF_nBasAux_Pair_wLD(iAtomPair),
     &                                      Work(ip_C),1,Work(ip_C),1))
            Work(ip_CNorm+4*(iAtomPair-1)+1)
     &                       =sqrt(dDot_(nuv*M,Work(ipC),1,Work(ipC),1))
         End If
         Do iD=1,nD
            ipD=iWork(ip_DBlocks(iD)-1+iAtomPair)
            ipV=iWork(ip_V(iD)-1+iAtom)
            Call dGeMV_('T',nuv,M,
     &                 1.0d0,Work(ipC),nuv,
     &                 Work(ipD),1,1.0d0,Work(ipV),1)
         End Do
         If (jAtom.ne.iAtom) Then
            ipC=ipC+nuv*M
            M=LDF_nBasAux_Atom(jAtom)
            If (doNorm) Then
               Work(ip_CNorm+4*(iAtomPair-1)+2)
     &                       =sqrt(dDot_(nuv*M,Work(ipC),1,Work(ipC),1))
            End If
            Do iD=1,nD
               ipD=iWork(ip_DBlocks(iD)-1+iAtomPair)
               ipV=iWork(ip_V(iD)-1+jAtom)
               Call dGeMV_('T',nuv,M,
     &                    1.0d0,Work(ipC),nuv,
     &                    Work(ipD),1,1.0d0,Work(ipV),1)
            End Do
         Else
            If (doNorm) Then
               Work(ip_CNorm+4*(iAtomPair-1)+2)
     &                                 =Work(ip_CNorm+4*(iAtomPair-1)+1)
            End If
         End If
         If (AP_2CFunctions(1,iAtomPair).gt.0) Then
            ipC=ipC+nuv*M
            M=AP_2CFunctions(1,iAtomPair)
            If (doNorm) Then
               Work(ip_CNorm+4*(iAtomPair-1)+3)
     &                       =sqrt(dDot_(nuv*M,Work(ipC),1,Work(ipC),1))
            End If
            Do iD=1,nD
               ipD=iWork(ip_DBlocks(iD)-1+iAtomPair)
               ipV=iWork(ip_V(iD)-1+nAtom+iAtomPair)
               Call dGeMV_('T',nuv,M,
     &                    1.0d0,Work(ipC),nuv,
     &                    Work(ipD),1,1.0d0,Work(ipV),1)
            End Do
         Else
            If (doNorm) Then
               Work(ip_CNorm+4*(iAtomPair-1)+3)=0.0d0
            End If
         End If
      End Do
      Call Free_Tsk(TaskListID)
      If (Timing) Then
         Call CWTIme(tC2,tW2)
         Write(6,'(A,2(1X,F12.2),A)')
     &   'Time spent computing Coulomb (V) intermediates:   ',
     &   tC2-tC1,tW2-tW1,' seconds'
      End If

      ! Deallocation
      Call GetMem('LDFCBlk','Free','Real',ip_C,l_C)

#if defined (_MOLCAS_MPP_)
      ! Get complete V and norm on all nodes
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         If (Timing) Call CWTime(tC1,tW1)
         Do iD=1,nD
            Call LDF_P_AddAuxBasVector(ip_V(iD))
         End Do
         If (doNorm) Then
            Call GAdGOp(Work(ip_CNorm),4*NumberOfAtomPairs,'+')
         End If
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Parallel overhead for Coulomb (V) intermediates:  ',
     &      tC2-tC1,tW2-tW1,' seconds'
         End If
      End If
#endif

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ComputeCoulombIntermediates0(nD,
     &                                            ip_DBlocks,ip_VBlocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute Coulomb intermediates
C
C              V(J) = sum_uv C(uv,J)*D(uv)
C
C              using LDF fitting coefficients.
C
C     BLOCKED VERSION
C
      Implicit None
      Integer nD
      Integer ip_DBlocks(nD)
      Integer ip_VBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer TaskListID
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer ip_C, l_C
      Integer iD
      Integer ipV, ipD
      Integer nuv, M

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

#if defined (_MOLCAS_MPP_)
      ! Initialize V arrays
      Do iD=1,nD
         Call LDF_ZeroBlockVector(ip_VBlocks(iD))
      End Do
#endif

      ! Allocate array for storing coefficients
      l_C=0
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nuv=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
         M=LDF_nBasAux_Pair(iAtomPair)
         l_C=max(l_C,nuv*M)
      End Do
      Call GetMem('LDFCBlk','Allo','Real',ip_C,l_C)

      ! Compute atom pair blocks of V
      Call Init_Tsk(TaskListID,NumberOfAtomPairs)
      Do While (Rsv_Tsk(TaskListID,iAtomPair))
         ! Get dimensions
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nuv=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
         M=LDF_nBasAux_Pair(iAtomPair)
         ! Read coefficients for this atom pair
         Call LDF_CIO_ReadC(iAtomPair,Work(ip_C),l_C)
         ! Compute V for this atom pair and for each density
         Do iD=1,nD
            ipD=iWork(ip_DBlocks(iD)-1+iAtomPair)
            ipV=iWork(ip_VBlocks(iD)-1+iAtomPair)
            Call dGeMV_('T',nuv,M,
     &                 1.0d0,Work(ip_C),nuv,
     &                 Work(ipD),1,0.0d0,Work(ipV),1)
         End Do
      End Do
      Call Free_Tsk(TaskListID)

      ! Deallocation
      Call GetMem('LDFCBlk','Free','Real',ip_C,l_C)

#if defined (_MOLCAS_MPP_)
      ! Get complete V on all nodes
      Do iD=1,nD
         Call LDF_P_AddBlockVector(ip_VBlocks(iD))
      End Do
#endif

      End
