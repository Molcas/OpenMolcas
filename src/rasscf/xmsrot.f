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
      Subroutine XMSRot(CMO,FI,FA)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on May. 21, 2020, created this file.               *
! ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: LROOTS, NAC
      use general_data, only: NTOT1,NTOT2

      Implicit None

!*****Input
      Real*8,DIMENSION(NTOT1):: FI,FA
      Real*8,Dimension(NTOT2)::CMO
!*****Auxiliary quantities
      Real*8,DIMENSION(:,:),Allocatable::FckO
!*****FckO:  Fock matrix for MO
      Real*8,DIMENSION(:,:),Allocatable::FckS,EigVec
!*****FckS:  Fock matrix for states
      Real*8,DIMENSION(:,:,:),Allocatable::GDMat
!*****GDMat: density matrix or transition density matrix

!     Allocating Memory
      CALL mma_allocate(GDMat,lRoots*(lRoots+1)/2,NAC,NAC)
      CALL mma_allocate(FckO,NAC,NAC)
      CALL mma_allocate(FckS,lRoots,lRoots)
      CALL mma_allocate(EigVec,lRoots,lRoots)
!
      CALL CalcFckO(CMO,FI,FA,FckO)

      CALL GetGDMat(GDMAt)
!
      CALL CalcFckS(FckO,GDMat,FckS)
!
      CALL CalcEigVec(FckS,lRoots,EigVec)
!
      call printmat('ROT_VEC','XMS-PDFT',eigvec,lroots,lroots,7,8,'N')
!     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(FckO)
      CALL mma_deallocate(FckS)
      CALL mma_deallocate(EigVec)

      End Subroutine XMSRot

!*****************************************************

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

!*****************************************************

      Subroutine GetGDMat(GDMat)
      use lucia_data, only: DStmp, Dtmp
      use stdalloc, only: mma_allocate, mma_deallocate
      use Lucia_Interface, only: Lucia_Util
      use rasscf_global, only: lRoots, nAc, iAdr15
      use general_data, only: JOBIPH,NCONF
      Implicit None

!     Output
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat

!     Auxiliary qunatities
      INTEGER CIDisk1,CIDisk2
      INTEGER NIJ2, jRoot, kRoot, iOrb, jOrb
      Real*8, Allocatable:: SDtmp(:), TmpD(:)
      Real*8, Allocatable:: VecL(:), VecR(:)


      Call mma_allocate(VecL,NConf,Label='VecL')
      Call mma_allocate(VecR,NConf,Label='VecR')
      Call mma_allocate(TmpD,NAC**2,Label='TmpD')
      Call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
      SDtmp(:)=DStmp(:)
      DStmp(:)=0.0D0
      TmpD(:)=Dtmp(:)
      Dtmp(:)=0.0D0
      CIDisk1=IADR15(4)
      Do jRoot=1,lRoots
       Call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
       CIDisk2=IADR15(4)
       Do kRoot=1,jRoot
        Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
!        write(6,*) 'VecL and VecR for states',jRoot,kRoot
!        write(6,*)(VecL(I),I=0,NConf-1)
!        write(6,*)(VecR(I),I=0,NConf-1)
        Call Lucia_Util('Densi',                                        &
     &                  CI_Vector=VecL(:),                              &
     &                  RVEC=VECR(:))
!        write(6,*)'GDMat for states',jRoot,kRoot
         dO IOrb=1,NAC
          do JOrb=1,NAC
          NIJ2=jRoot*(jRoot-1)/2+kRoot
          GDMat(NIJ2,JOrb,IOrb)=Dtmp(JOrb+(IOrb-1)*NAC)
          end do
!          write(6,'(10(F8.4,2X))')(GDMat(NIJ2,IOrb,JOrb),JOrb=1,NAC)
         eND dO
       End Do
      End DO
      DStmp(:)=SDtmp(:)
      DTmp(:)=TmpD(:)
      Call mma_deallocate(SDtmp)
      Call mma_deallocate(TmpD)
      Call mma_deallocate(VecL)
      Call mma_deallocate(VecR)

      END Subroutine GetGDMat

!*****************************************************

      Subroutine CalcFckS(FckO,GDMat,FckS)
      use rasscf_global, only: lRoots, nAc
      Implicit None


!*****Input
      Real*8,DIMENSION(NAC,NAC)::FckO
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat
!*****Output
      Real*8,DIMENSION(lRoots,lRoots)::FckS
!*****Auxiliary variables
      INTEGER IState,JState, iOrb, jOrb

      FckS(:,:)=0.0d0

      DO IState=1,lRoots
       Do JState=1,IState
        dO IOrb=1,NAC
         do JOrb=1,NAC
          FckS(IState,JState)=FckS(IState,JState)+FckO(IOrb,JOrb)*      &
     &GDMat(IState*(IState-1)/2+JState,IOrb,JOrb)
         end do
        eND DO
        FckS(JState,IState)=FckS(IState,JState)
        End Do
      END DO

!      CALL PrintMat('XMS_Mat','test',FckS,LRoots,LRoots,0,4,'N')

      END Subroutine CalcFckS

!*****************************************************

      Subroutine CalcEigVec(Matrix,NDIM,EigVec)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None

!*****Input
      INTEGER NDim
!*****Input & Output
      Real*8,DIMENSION(NDIM,NDIM)::Matrix,EigVec
!*****Calculating rotation matrix
      Real*8, Allocatable:: Mat(:),Val(:,:), Scr(:)
      INTEGER NScr,INFO
      Real*8,DIMENSION(2)::WGRONK
!*****Auxiliary quantities
      INTEGER NElem ! NElem=NDim**2
      INTEGER IRow,ICol,IRIC,NI
      LOGICAL UseJacob

      UseJacob=.true.
      EigVec(:,:)=0.0d0

      IF(UseJacob) THEN
       NElem=NDim*(NDim+1)/2
       Call mma_allocate(Mat,nElem,Label='Mat')
       Call mma_allocate(Val,nDim,nDim,Label='Val')
       IRIC=0
       DO IRow=1,NDIM
        Do ICol=1,IRow
         IRIC=IRIC+1
         Mat(IRIC)=Matrix(IRow,ICol)
        End Do
       END DO
       Val(:,:)=0.0D0
       DO NI=1,NDim
        Val(NI,NI)=1.0D0
       END DO
!       write(6,*)'eigenvector matrix before diag'
!       CALL RECPRT(' ',' ',Val,NDIM,NDIM)
!       write(6,*)'matrix to be diagonalized'
!       CALL TriPrt(' ',' ',Mat,NDIM)
       CALL JACOB(Mat,Val,NDim,NDim)
!       write(6,*)'eigenvector matrix'
!       CALL RECPRT(' ',' ',Val,NDIM,NDIM)
!       DO IRow=1,NDIM
!         write(6,*) (EigVec(IRow,ICol), ICol=1,NDim)
!       END DO
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=Val(iCol,iRow)
        End Do
       END DO
       Call mma_deallocate(Val)
       Call mma_deallocate(Mat)

      ELSE
       NElem=NDim**2
       Call mma_allocate(Mat,nElem,Label='Mat')
       Call mma_allocate(Val,nDim,nDim,Label='Val')
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         Mat((ICol-1)*NDIM+IRow)=Matrix(IRow,ICol)
        End Do
       END DO
       Val(:,:)=0.0D0
       Call Dsyev_('V','U',NDIM,Mat,NDIM,Val,WGRONK,-1,INFO)
       NScr=Int(WGRONK(1))
       Call mma_allocate(Scr,nScr,Label='Scr')
       Call Dsyev_('V','U',NDIM,Mat,NDIM,Val,Scr,NScr,INFO)
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=Mat((IRow-1)*NDIM+ICol)
        End Do
       END DO
       Call mma_deallocate(Scr)
       Call mma_deallocate(Val)
       Call mma_deallocate(Mat)
       END IF

      End Subroutine CalcEigVec

!*****************************************************

      Subroutine PrintMat(FileName,MatInfo,Matrix,NRow,NCol,            &
     &                    LenName,LenInfo,Trans)
      Implicit None

      INTEGER NRow,NCol,LenName,LenInfo
      CHARACTER(Len=LenName)::FileName
      CHARACTER(Len=LenInfo)::MatInfo
      CHARACTER(Len=1)::Trans
      CHARACTER(Len=80)::PrtFmt
      Real*8,DIMENSION(NRow,NCol)::Matrix

      INTEGER LU,IsFreeUnit,IRow,ICol
      External IsFreeUnit


      IF(LenName.gt.0) THEN
      LU=100
      LU=IsFreeUnit(LU)
      CALL Molcas_Open(LU,FileName)
      ELSE
      LU=6
      END IF

      IF(Trans.eq.'N') THEN
       WRITE(PrtFmt,'(A,I5,A)')                                         &
     & '(',NCol,'(ES24.14E4,1X))'
       DO IRow=1,NRow
        write(LU,PrtFmt)                                                &
     &  (Matrix(IRow,ICol),ICol=1,NCol)
       END DO
      ELSE
       WRITE(PrtFmt,'(A,I5,A)')                                         &
     & '(',NRow,'(ES24.14E4,1X))'
       DO ICol=1,NCol
        write(LU,PrtFmt)                                                &
     &  (Matrix(IRow,ICol),IRow=1,NRow)
       END DO
      END IF
      WRITE(LU,*)MatInfo
      IF(LenName.gt.0) THEN
       Close(LU)
      END IF
      End Subroutine PrintMat

!*****************************************************

      Subroutine ReadMat(FileName,MatInfo,Matrix,NRow,NCol,             &
     &                   LenName,LenInfo,Trans)
      Implicit None

      INTEGER NRow,NCol,LenName,LenInfo
      CHARACTER(Len=LenName)::FileName
      CHARACTER(Len=LenInfo)::MatInfo
      CHARACTER(Len=1)::Trans
      Real*8,DIMENSION(NRow,NCol)::Matrix

      INTEGER LU,IsFreeUnit,IRow,ICol
      External IsFreeUnit

      IF(LenName.gt.0) THEN
      LU=100
      LU=IsFreeUnit(LU)
      CALL Molcas_Open(LU,FileName)
      ELSE
      LU=6
      END IF
      IF(Trans.eq.'N') THEN
       DO IRow=1,NRow
        Read(LU,*) (Matrix(IRow,ICol),ICol=1,NCol)
       END DO
      ELSE
       DO ICol=1,NCol
        Read(LU,*) (Matrix(IRow,ICol),IRow=1,NRow)
       END DO
      END IF
      Read(LU,*)MatInfo
      IF(LenName.gt.0) THEN
       Close(LU)
      END IF
      End Subroutine ReadMat
