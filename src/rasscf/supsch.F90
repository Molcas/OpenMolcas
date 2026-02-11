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
! Copyright (C) 1997, Luis Serrano-Andres                              *
!***********************************************************************
      SUBROUTINE SUPSCH(SMAT,CMOO,CMON)
!
!     Program RASSCF
!
!     Objective: To check the order of the input of natural orbitals
!                to obtain the right labels for the Supersymmetry matrix.
!
!     Called from ortho, neworb, fckpt2, and natorb.
!
!     Luis Serrano-Andres
!     University of Lund, Sweden, 1997
!     **** Molcas-4 *** Release 97 04 01 **********
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use general_data, only: NSYM,NBAS

      IMPLICIT None
      Real*8 CMOO(*),CMON(*),SMAT(*)

      Real*8, Allocatable:: Temp1(:), Temp2(:)
      Integer, Allocatable:: IxSym2(:)
      Integer :: iSym, nOrb_Tot, nOrbMx
!
!
      nOrbMX=0
      nOrb_tot=0
      Do iSym=1,nSym
         nOrbMX=Max(nOrbMX,nBas(iSym))
         nOrb_tot=nOrb_tot+nBas(iSym)
      End Do
!
      Call mma_allocate(Temp1,nOrbMX*nOrbMX,Label='Temp1')
      Call mma_allocate(Temp2,nOrbMX*nOrbMX,Label='Temp2')
      Call mma_allocate(IxSym2,nOrb_tot,Label='IxSym2')
!
      Call SUPSCH_(SMAT,CMOO,CMON,Temp1,Temp2,nOrbMX,IxSym2,nOrb_tot)
!
      Call mma_deallocate(IxSym2)
      Call mma_deallocate(Temp2)
      Call mma_deallocate(Temp1)
!
!
      End SUBROUTINE SUPSCH
