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
      Subroutine NatOrb_MCLR(Dens,CMOO,CMON,OCCN)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One
      use MCLR_Data, only: ipCM, ipMat, nDens2
      use input_mclr, only: nSym,nBas,kPrint
      Implicit None
      Real*8 Dens(*),CMOO(*),CMON(*),OCCN(*)

      Real*8, Allocatable:: EVal(:), EVec(:)
      Integer iO, iS, ij, i, j, ii, iSt, iEnd

      Call mma_allocate(EVec,ndens2,Label='EVec')
      Call mma_allocate(EVal,ndens2,Label='EVal')
!
!         Diagonalize the density matrix and transform orbitals
!
      If (iAnd(kprint,8).eq.8) Then
         Write(6,*)
         Write(6,*) '           Effective natural population '
         Write(6,*) '           ============================ '
         Write(6,*)
      End If
      io=0
      Do is=1,nsym
         ij=0
         Do i=0,nbas(is)-1
            Do j=0,i
               ij=ij+1
               Eval(ij)=Dens(ipMat(is,is)+i+j*nbas(is))
            End DO
         End DO
         EVec(:)=Zero
         Call dCopy_(nBas(is),[One],0,EVec,nbas(is)+1)
         CALL JACOB(EVal,EVec,nbas(is),nbas(is))
         ii=0
         DO i=1,nbas(is)
            ii=ii+i
            OCCN(io+i)=Eval(ii)
         END DO
         IST=IO+1
         IEND=IO+NBAS(is)
         If (iAnd(kprint,2).eq.2)                                       &
     &      Write (6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')           &
     &             'sym',iS,':',(OCCN(I),I=IST,IEND)
         If (nBas(is).ge.1)                                             &
     &      CALL DGEMM_('N','N',                                        &
     &                  NBAS(is),NBAS(is),NBAS(is),                     &
     &                  One,CMOO(ipCM(is)),NBAS(is),                    &
     &                  EVec,NBAS(is),                                  &
     &                  Zero,CMON(ipCM(is)),NBAS(is))
         io=io+nbas(is)
      End DO

      Call mma_deallocate(EVec)
      Call mma_deallocate(Eval)

      End Subroutine NatOrb_MCLR
