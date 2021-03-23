!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Annihil_rho(Dmat,nBas)

      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8 Dmat(*)
      Integer nBas
      Character*(LENIN4) Name(mxBas)
      Integer, Allocatable:: nBas_per_Atom(:), nBas_Start(:)
      Real*8, Allocatable:: Charge_B(:)
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Get_iScalar('Unique atoms',nAtoms)

      If (nAtoms.lt.1) Then
         Write(6,'(A,I9)') 'nUniqAt =',nAtoms
         Call Abend()
      End If
!
      Call mma_allocate(nBas_per_Atom,nAtoms,Label='nBpA')
      Call mma_allocate(nBas_Start,nAtoms,Label='nB_Start')

      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas)

      Call BasFun_Atom(nBas_per_Atom,nBas_Start,                        &
     &                 Name,nBas,nAtoms,.false.)
!
      Call mma_allocate(Charge_B,nAtoms,Label='Charge_B')
      Call Get_dArray('Nuclear charge',Charge_B,nAtoms)
!
      ZA=0.0d0
      iAt=1
      Do while (iAt.le.nAtoms .and. ZA.eq.0.0d0)
         ZA = Charge_B(iAt)
         iAt = iAt + 1
      End Do
      iAt_B=iAt-1  ! start of atoms of subsystem B
      Call mma_deallocate(Charge_B)
!
      If (iAt_B.eq.1) Then ! subsystem B comes first
         ZB=1.0d0
         nAt_B=1
         Do while (nAt_B.le.nAtoms .and. ZB.gt.0.0d0)
            ZB = Charge_B(nAt_B)
            nAt_B = nAt_B + 1
         End Do
         nAt_B=nAt_B-1  ! end of atoms of subsystem B
         nBas_B = nBas_Start(nAt_B) - 1
         Do j=nBas_B,nBas-1
            jj=j*(j+1)/2
            Do i=1,j
               ijj=i+jj
               Dmat(ijj)=0.0d0
            End Do
         End Do
      Else
         nBas_A = nBas_Start(iAt_B) - 1
         nAA=nBas_A*(nBas_A+1)/2
         Call FZero(Dmat,nAA)
         Do j=nBas_A,nBas-1
            jj=j*(j+1)/2
            Do i=1,nBas_A
               ijj=i+jj
               Dmat(ijj)=0.0d0
            End Do
         End Do
      EndIf
!
      Call mma_deallocate(nBas_Start)
      Call mma_deallocate(nBas_per_Atom)
!
!  Annihilated density written to runfile for use in Coulomb gradients
!
      Length=nBas*(nBas+1)/2
      Call Put_D1ao_Var(Dmat,Length)

      Return
      End
