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
      Subroutine LDF_SetChargeConstraint()
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: set info for charge constraint.
C              To unset, call LDF_UnsetChargeConstraint()
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_charge_constraint_info.fh"
#include "ldf_atom_pair_info.fh"

      Character*23 SecNam
      Parameter (SecNam='LDF_SetChargeConstraint')

      Integer  LDF_nAtom, LDF_nBasAux_Atom, LDF_nBas_Atom
      External LDF_nAtom, LDF_nBasAux_Atom, LDF_nBas_Atom

      Character*8 Label

      Integer nAtom
      Integer l
      Integer ip0
      Integer A, AB
      Integer ip

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (ChargeConstraintSet) Return

      ! Allocate memory for storing integrals over 1C aux functions
      nAtom=LDF_nAtom()
      l=nAtom
      Call GetMem('CCAIV_BP','Allo','Inte',ip_CC_AuxIntVec,l)
      l=0
      ip0=ip_CC_AuxIntVec-1
      Do A=1,nAtom
         iWork(ip0+A)=l
         l=l+LDF_nBasAux_Atom(A)
      End Do
      Call GetMem('CCAuxInt','Allo','Real',ip,l)
      Do A=1,nAtom
         iWork(ip0+A)=iWork(ip0+A)+ip
      End Do

      ! Compute integrals over 1C aux functions
      Label='Mltpl  0'
      Call LDF_SetOneEl(Label)
      Do A=1,nAtom
         l=LDF_nBasAux_Atom(A)
         ip=iWork(ip0+A)
         Call LDF_ComputeAuxInt_1(A,l,Work(ip))
      End Do
      Call LDF_UnsetOneEl(Label)

      ! Allocate memory for largest overlap and multiplier blocks
      l=0
      Do AB=1,NumberOfAtomPairs
         l=max(l,LDF_nBas_Atom(AP_Atoms(1,AB))
     &          *LDF_nBas_Atom(AP_Atoms(2,AB)))
      End Do
      l_CC_Overlap=l
      Call GetMem('CLDFOv','Allo','Real',ip_CC_Overlap,l_CC_Overlap)
      l_CC_lambda=l
      Call GetMem('CLDFla','Allo','Real',ip_CC_lambda,l_CC_lambda)

      ChargeConstraintSet=.True.

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetChargeConstraint()
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: unset info for charge constraint.
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_charge_constraint_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Integer nAtom
      Integer l
      Integer A
      Integer ip

      If (.not.ChargeConstraintSet) Return

      ! Deallocate memory for storing integrals over 1C aux functions
      nAtom=LDF_nAtom()
      l=0
      Do A=1,nAtom
         l=l+LDF_nBasAux_Atom(A)
      End Do
      ip=iWork(ip_CC_AuxIntVec)
      Call GetMem('CCAuxInt','Free','Real',ip,l)
      l=nAtom
      Call GetMem('CCAIV_BP','Free','Inte',ip_CC_AuxIntVec,l)
      ip_CC_AuxIntVec=0

      ! Deallocate largest overlap and multiplier blocks
      Call GetMem('CLDFOv','Free','Real',ip_CC_Overlap,l_CC_Overlap)
      l_CC_Overlap=0
      ip_CC_Overlap=0
      Call GetMem('CLDFla','Free','Real',ip_CC_lambda,l_CC_lambda)
      l_CC_lambda=0
      ip_CC_lambda=0

      ChargeConstraintSet=.False.

      End
