!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
      Subroutine Mk_EOrb()
      Use SCF_Arrays, only: FockAO, EOrb, CMO
      use InfSCF, only: nSym, nBas, nOrb, nD
      Implicit None
!
      Integer nFck, nEOrb, nCMO, iD
!
      nFck =SIZE(FockAO,1)
      nEOrb=SIZE(EOrb,1)
      nCMO =SIZE(CMO,1)

      Do iD = 1, nD
         Call MkEorb_(FockAO(:,iD),nFck,CMO(:,iD),nCMO,EOrb(:,iD),nEorb,nSym,nBas,nOrb)
         If (iD==1) Then
            Call Put_darray('OrbE',   Eorb(:,iD),SIZE(EOrb,1))
         Else
            Call Put_darray('OrbE_ab',Eorb(:,iD),SIZE(EOrb,1))
         End If
      End Do
!
      Return
      End Subroutine Mk_EOrb
      Subroutine MkEorb_(FockAO,nFck,CMO,nCMO,Eorb,nEorb,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
!  This routine calculates the diagonal elements of the MO Fock matrix *
!  (orbital energies).                                                 *
!                                                                      *
!  Input:                                                              *
!    FockAO  Fock matrix in AO basis                                   *
!    CMO     Orbitals                                                  *
!                                                                      *
!  Output:                                                             *
!    Eorb    Orbital energies.                                         *
!                                                                      *
!***********************************************************************
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
      Integer nFck, nCMO, nEOrb
      Real*8 FockAO(nFck)
      Real*8 CMO(nCMO)
      Real*8 EOrb(nEOrb)
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Real*8  t
      Real*8, Dimension(:), Allocatable:: FckSqr
      Integer iSym
      Integer iBas
      Integer jBas
      Integer iOrb
      Integer MaxTri
      Integer MaxSqr
      Integer iOffTri
      Integer iOffCMO
      Integer npFckSqr
      Integer indE
      Integer indF
      Integer indx
      Integer jndx
!----------------------------------------------------------------------*
! Some preliminary setup.                                              *
!----------------------------------------------------------------------*
      MaxTri = 0
      MaxSqr = 0
      Do iSym=1,nSym
         MaxTri=Max(MaxTri,nBas(iSym)*(nBas(iSym)+1)/2)
         MaxSqr=Max(MaxSqr,nBas(iSym)*nBas(iSym))
      End Do
!----------------------------------------------------------------------*
! Allocate matrices.                                                   *
!----------------------------------------------------------------------*
      npFckSqr=MaxSqr
      Call mma_allocate(FckSqr,npFckSqr,Label='FckSqr')
!----------------------------------------------------------------------*
! Do compute orbital energies                                          *
!----------------------------------------------------------------------*
      iOffTri=0
      iOffCMO=0
      indE=1
      Do iSym=1,nSym
         If(nOrb(iSym)>0) Then
            Call Square(FockAO(1+iOffTri),FckSqr,1,nBas(iSym),nBas(iSym))
            Do iOrb=1,nOrb(iSym)
               t=Zero
               indF=1
               Do iBas=1,nBas(iSym)
                  Do jBas=1,nBas(iSym)
                     indx=iBas+(iOrb-1)*nBas(iSym)+iOffCMO
                     jndx=jBas+(iOrb-1)*nBas(iSym)+iOffCMO
                     t=t+CMO(indx)*CMO(jndx)*FckSqr(indF)
                     indF=indF+1
                  End Do
               End Do
               EOrb(indE)=t
               indE=indE+1
            End Do
         End If
         iOffTri=iOffTri+nBas(iSym)*(nBas(iSym)+1)/2
         iOffCMO=iOffCMO+nBas(iSym)*nOrb(iSym)
      End Do
!----------------------------------------------------------------------*
! Debug print                                                          *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      indE=1
      Do iSym=1,nSym
         If(nOrb(iSym).gt.0) Then
            Call RecPrt('mkeor: Orbital energies','(20F10.4)',Eorb(indE),1,nOrb(iSym))
         End If
         indE=indE+nOrb(iSym)
      End Do
#endif
!----------------------------------------------------------------------*
! Deallocate matrices.                                                 *
!----------------------------------------------------------------------*
      Call mma_deallocate(FckSqr)
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
      Return
      End Subroutine MkEorb_
