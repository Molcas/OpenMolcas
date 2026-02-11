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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
      Subroutine CalcFckO(CMO,FI,FA,FckO)
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: NAC
      use general_data, only: NTOT1,NTOT2,NSYM,NASH,NBAS,NFRO,NISH
      Implicit None

!*****Input
      Real*8,DIMENSION(NTOT1)::FI,FA
      Real*8,Dimension(NTOT2)::CMO
!*****Output
      Real*8 ::FckO(NAC,NAC)
!*****Auxiliary quantities
      INTEGER NB,NA,NI,IOff1,IOff2,IOff3
      INTEGER IBas,JBas,ISym,IOrb,JOrb
      Real*8, Allocatable:: FIAAO(:,:), Scr(:,:), FckOt(:,:)

      FckO(:,:)=0.0d0

      IOff1=0
      IOff2=1
      IOff3=0
      DO ISym=1,NSym
        NB=NBas(ISym)
        NA=NASH(ISym)
        NI=NISH(ISym)+NFro(ISym)
       IF(NASH(ISym).gt.0) THEN
        Call mma_allocate(FIAAO,nB,nB,Label='FIAAO')
        Call mma_allocate(Scr,nB,nA,Label='Scr')
        Call mma_allocate(FckOt,NA,NA,Label='FckOt')
        FckOt(:,:)=0.0D0
!        write(6,*)'Print FI Matrix'
!        Do IBas=1,NB
!         write(6,*)(FI(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
!        End Do
!        write(6,*)'Print FA Matrix'
!        Do IBas=1,NB
!         write(6,*)(FA(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
!        End Do
!        write(6,*)'Active CMO mat for sym',ISym
!        CALL RecPrt(' ',' ',CMO(IOff2+NI*NB),NB,NA)
        Do IBas=1,NB
         do JBas=1,IBas
          FIAAO(jBas,iBas)= FI(IOff1+(IBas-1)*IBas/2+JBas)+             &
     &                      FA(IOff1+(IBas-1)*IBas/2+JBas)
          FIAAO(iBas,jBas)=FIAAO(jBas,iBas)
         end do
        End Do
!        write(6,*)'FIA mat for sym',ISym
!        CALL RecPrt(' ',' ',FIAAO,NB,NB)
        CALL DGEMM_('n','n',NB,NA,NB,1.0D0,FIAAO,                       &
     &       NB,CMO(IOff2+NI*NB),NB,0.0D0,Scr,NB)
        CALL DGEMM_('t','n',NA,NA,NB,1.0D0,CMO(IOff2+NI*NB),            &
     &       NB,Scr,NB,0.0D0,FckOt,NA)
!        write(6,*)'FckO mat for sym',ISym
!        CALL RecPrt(' ',' ',FckOt,NA,NA)
        Do IOrb=1,NA
         do JOrb=1,NA
          FckO(IOrb+IOff3,JOrb+IOff3)= FckOt(jOrb,iOrb)
         end do
        End Do
        Call mma_deallocate(FIAAO)
        Call mma_deallocate(Scr)
        Call mma_deallocate(FckOt)
       END IF
       IOff1=IOff1+NB*(NB+1)/2
       IOff2=IOff2+NB**2
       IOff3=IOff3+NA
      END DO

      END Subroutine CalcFckO
