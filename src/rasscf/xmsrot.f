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
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
      Subroutine XMSRot(CMO,FI,FA)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on May. 21, 2020, created this file.               *
* ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
#include "rasscf_lucia.fh"

******Input
      Real*8,DIMENSION(NTOT1):: FI,FA
      Real*8,Dimension(NTOT2)::CMO
******Auxillary quantities
      Real*8,DIMENSION(:,:),Allocatable::FckO
******FckO:  Fock matrix for MO
      Real*8,DIMENSION(:,:),Allocatable::FckS,EigVec
******FckS:  Fock matrix for states
      Real*8,DIMENSION(:,:,:),Allocatable::GDMat
******GDMat: density matrix or transition density matrix

C     Allocating Memory
      CALL mma_allocate(GDMat,lRoots*(lRoots+1)/2,NAC,NAC)
      CALL mma_allocate(FckO,NAC,NAC)
      CALL mma_allocate(FckS,lRoots,lRoots)
      CALL mma_allocate(EigVec,lRoots,lRoots)
C
      CALL CalcFckO(CMO,FI,FA,FckO)

      CALL GetGDMat(GDMAt)
C
      CALL CalcFckS(FckO,GDMat,FckS)
C
      CALL CalcEigVec(FckS,lRoots,EigVec)
C
      call printmat('ROT_VEC','XMS-PDFT',eigvec,lroots,lroots,7,8,'N')
C     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(FckO)
      CALL mma_deallocate(FckS)
      CALL mma_deallocate(EigVec)

      RETURN
      End Subroutine

***********************************************************************

******************************************************
      Subroutine CalcFckO(CMO,FI,FA,FckO)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
#include "rasscf_lucia.fh"
******Input
      Real*8,DIMENSION(NTOT1)::FI,FA
      Real*8,Dimension(NTOT2)::CMO
******Output
      Real*8,DIMENSION(NAC,NAC)::FckO
******Auxillary quantities
      INTEGER LFIAAO,LScr,NB,NA,NI,IOff1,IOff2,IOff3,LFckOt
      INTEGER IBas,JBas,ISym,IOrb,JOrb

      DO IOrb=1,NAC
       Do JOrb=1,NAC
        FckO(JOrb,IOrb)=0.0d0
       End Do
      END DO

      IOff1=0
      IOff2=1
      IOff3=0
      DO ISym=1,NSym
        NB=NBas(ISym)
        NA=NASH(ISym)
        NI=NISH(ISym)+NFro(ISym)
       IF(NASH(ISym).gt.0) THEN
        Call GetMem('FIAAO','ALLO','REAL',LFIAAO,NB**2)
        Call GetMem('Scra' ,'ALLO','REAL',LScr,  NB*NA)
        Call GetMem('FckOt','ALLO','REAL',LFckOt,NA*NA)
        CALL FZero(WORK(LFckOt),NA**2)
C        write(6,*)'Print FI Matrix'
C        Do IBas=1,NB
C         write(6,*)(FI(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
C        End Do
C        write(6,*)'Print FA Matrix'
C        Do IBas=1,NB
C         write(6,*)(FA(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
C        End Do
C        write(6,*)'Active CMO mat for sym',ISym
C        CALL RecPrt(' ',' ',CMO(IOff2+NI*NB),NB,NA)
        Do IBas=1,NB
         do JBas=1,IBas
          WORK(LFIAAO+(IBas-1)*NB+JBas-1)=
     & FI(IOff1+(IBas-1)*IBas/2+JBas)+
     & FA(IOff1+(IBas-1)*IBas/2+JBas)
          WORK(LFIAAO+(JBas-1)*NB+IBas-1)=
     & WORK(LFIAAO+(IBas-1)*NB+JBas-1)
         end do
        End Do
C        write(6,*)'FIA mat for sym',ISym
C        CALL RecPrt(' ',' ',Work(LFIAAO),NB,NB)
        CALL DGEMM_('n','n',NB,NA,NB,1.0D0,Work(LFIAAO),
     &       NB,CMO(IOff2+NI*NB),NB,0.0D0,Work(LScr),NB)
        CALL DGEMM_('t','n',NA,NA,NB,1.0D0,CMO(IOff2+NI*NB),
     &       NB,Work(LScr),NB,0.0D0,Work(LFckOt),NA)
C        write(6,*)'FckO mat for sym',ISym
C        CALL RecPrt(' ',' ',Work(LFckOt),NA,NA)
        Do IOrb=1,NA
         do JOrb=1,NA
          FckO(IOrb+IOff3,JOrb+IOff3)=
     &Work(LFckOt+(IOrb-1)*NA+JOrb-1)
         end do
        End Do
        Call GetMem('FIAAO','FREE','REAL',LFIAAO,NB**2)
        Call GetMem('Scra' ,'FREE','REAL',LScr,  NB*NA)
        Call GetMem('FckOt','FREE','REAL',LFckOt,NA*NA)
       END IF
       IOff1=IOff1+NB*(NB+1)/2
       IOff2=IOff2+NB**2
       IOff3=IOff3+NA
      END DO

      RETURN
      END Subroutine

******************************************************

******************************************************
      Subroutine GetGDMat(GDMat)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
#include "rasscf_lucia.fh"
*     Output
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat
*     Auxillary qunatities
      INTEGER CIDisk1,CIDisk2,iVecL,iVecR,iDummy
      INTEGER tlw6,tlw7,ldtmp,lsdtmp,NIJ2
      Dimension Dummy(1)
      tlw6=lw6
      tlw7=lw7
      Call GetMem('LVEC','ALLO','REAL',iVecL,NConf)
      Call GetMem('RVEC','ALLO','REAL',iVecR,NConf)
      Call GetMem('Dtmp','ALLO','REAL',ldtmp,NAC**2)
      Call GetMem('SDtmp','ALLO','REAL',lsdtmp,NAC**2)
      lw6=ldtmp
      lw7=lsdtmp
      CIDisk1=IADR15(4)
      Do jRoot=1,lRoots
       Call DDafile(JOBIPH,2,Work(iVecL),nConf,CIDisk1)
       C_Pointer=iVecL
       CIDisk2=IADR15(4)
       Do kRoot=1,jRoot
        Call DDafile(JOBIPH,2,Work(iVecR),nConf,CIDisk2)
C        write(6,*) 'VecL and VecR for states',jRoot,kRoot
C        write(6,*)(WORK(iVecL+I),I=0,NConf-1)
C        write(6,*)(WORK(iVecR+I),I=0,NConf-1)
        Call Lucia_Util('Densi',iVecR,iDummy,Dummy)
C        Call Lucia_Util('Densi',0,iDummy,rdum)
C        write(6,*)'GDMat for states',jRoot,kRoot
         dO IOrb=1,NAC
          do JOrb=1,NAC
          NIJ2=jRoot*(jRoot-1)/2+kRoot
          GDMat(NIJ2,IOrb,JOrb)=WORK(LW6+JOrb-1+(IOrb-1)*NAC)
          end do
C          write(6,'(10(F8.4,2X))')(GDMat(NIJ2,IOrb,JOrb),JOrb=1,NAC)
         eND dO
       End Do
      End DO
      lw6=tlw6
      lw7=tlw7
      Call GetMem('LVEC','FREE','REAL',iVecL,NConf)
      Call GetMem('RVEC','FREE','REAL',iVecR,NConf)
      Call GetMem('Dtmp','FREE','REAL',ldtmp,NAC**2)
      Call GetMem('SDtmp','Free','REAL',lsdtmp,NAC**2)
      RETURN
      END Subroutine

******************************************************
      Subroutine CalcFckS(FckO,GDMat,FckS)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
#include "rasscf_lucia.fh"

******Input
      Real*8,DIMENSION(NAC,NAC)::FckO
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat
******Output
      Real*8,DIMENSION(lRoots,lRoots)::FckS
******Auxillary variables
      INTEGER IState,JState

      DO IState=1,lRoots
       Do JState=1,IState
        FckS(IState,JState)=0.0d0
       End Do
      END DO

      DO IState=1,lRoots
       Do JState=1,IState
        dO IOrb=1,NAC
         do JOrb=1,NAC
          FckS(IState,JState)=FckS(IState,JState)+FckO(IOrb,JOrb)*
     &GDMat(IState*(IState-1)/2+JState,IOrb,JOrb)
         end do
        eND DO
        FckS(JState,IState)=FckS(IState,JState)
        End Do
      END DO

C      CALL PrintMat('XMS_Mat','test',FckS,LRoots,LRoots,0,4,'N')


      RETURN
      END SUBROUTINE
******************************************************


******************************************************
      Subroutine CalcEigVec(Matrix,NDIM,EigVec)

#include "WrkSpc.fh"
******Input
      INTEGER NDim
******Input & Output
      Real*8,DIMENSION(NDIM,NDIM)::Matrix,EigVec
******Calculating rotation matrix
      INTEGER LMat,LVal,NScr,INFO
      Real*8,DIMENSION(2)::WGRONK
******Auxillary quantities
      INTEGER NElem ! NElem=NDim**2
      INTEGER IRow,ICol,IRIC
      LOGICAL UseJacob

      UseJacob=.true.

      DO ICol=1,NDIM
       Do IRow=1,NDIM
        EigVec(IRow,ICol)=0.0d0
       End Do
      END DO
      IF(UseJacob) THEN
       NElem=NDim*(NDim+1)/2
       Call GetMem('Mat','ALLO','REAL',LMat,NElem)
       Call GetMem('EVa','ALLO','REAL',LVal,NDIM**2)
       IRIC=0
       DO IRow=1,NDIM
        Do ICol=1,IRow
         WORK(LMat+IRIC)=Matrix(IRow,ICol)
         IRIC=IRIC+1
        End Do
       END DO
       CALL FZero(WORK(LVal),NDIM**2)
       IRIC=0
       DO NI=1,NDim
        WORK(LVal+IRIC)=1.0D0
        IRIC=IRIC+NDIM+1
       END DO
C       write(6,*)'eigenvector matrix before diag'
C       CALL RECPRT(' ',' ',WORK(LVal),NDIM,NDIM)
C       write(6,*)'matrix to be diagonalized'
C       CALL TriPrt(' ',' ',WORK(LMat),NDIM)
       CALL JACOB(WORK(LMat),WORK(LVal),NDim,NDim)
C       write(6,*)'eigenvector matrix'
C       CALL RECPRT(' ',' ',WORK(LVal),NDIM,NDIM)
C       DO IRow=1,NDIM
C         write(6,*) (EigVec(IRow,ICol), ICol=1,NDim)
C       END DO
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=WORK(LVal+(ICol-1)*NDIM+IRow-1)
        End Do
       END DO
       Call GetMem('EVa','FREE','REAL',LVal,NElem)
       Call GetMem('Mat','FREE','REAL',LMat,NElem)

      ELSE
       NElem=NDim**2
       Call GetMem('Mat','ALLO','REAL',LMat,NElem)
       Call GetMem('EVa','ALLO','REAL',LVal,NDIM**2)
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         WORK(LMat+(ICol-1)*NDIM+IRow-1)=Matrix(IRow,ICol)
        End Do
       END DO
       CALL FZERO(WORK(LVal),NElem)
       Call Dsyev_('V','U',NDIM,Work(LMat),NDIM,Work(LVal),
     &             WGRONK,-1,INFO)
       NScr=Int(WGRONK(1))
       Call GetMem('Scr','Allo','Real',LScr,NScr)
       Call Dsyev_('V','U',NDIM,Work(LMat),NDIM,Work(LVal),
     &                Work(LScr),NScr,INFO)
       DO ICol=1,NDIM
        Do IRow=1,NDIM
         EigVec(IRow,ICol)=WORK(LMat+(ICol-1)*NDIM+IRow-1)
        End Do
       END DO
       Call GetMem('Scr','Free','Real',LScr,NScr)
       Call GetMem('EVa','FREE','REAL',LVal,NElem)
       Call GetMem('Mat','FREE','REAL',LMat,NElem)
       END IF

      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine PrintMat(FileName,MatInfo,Matrix,NRow,NCol,
     &LenName,LenInfo,Trans)

      INTEGER NRow,NCol,LenName
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
        write(LU,*) (Matrix(IRow,ICol),ICol=1,NCol)
       END DO
      ELSE
       DO ICol=1,NCol
        write(LU,*) (Matrix(IRow,ICol),IRow=1,NRow)
       END DO
      END IF
      WRITE(LU,*)MatInfo
      IF(LenName.gt.0) THEN
       Close(LU)
      END IF
      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine ReadMat(FileName,MatInfo,Matrix,NRow,NCol,
     &LenName,LenInfo,Trans)

      INTEGER NRow,NCol,LenName
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
      RETURN
      End Subroutine
******************************************************


******************************************************
       Subroutine MatToWork2DRR(Mat,L,M,LWork,Trans)
#include "WrkSpc.fh"
       INTEGER L, M, LWork
       Real*8,DIMENSION(L,M)::Mat
       INTEGER I,J
       CHARACTER(Len=1)::Trans
       IF(Trans.ne.'T') THEN
        DO I=1,L
         Do J=1,M
          WORK(LWork+J-1+(I-1)*M)=Mat(J,I)
         End Do
        END DO
       ELSE
        DO I=1,L
         Do J=1,M
          WORK(LWork+J-1+(I-1)*M)=Mat(I,J)
         End Do
        END DO
       END IF
       RETURN
       End Subroutine
******************************************************

******************************************************
       Subroutine WorkToMat2DRR(Mat,L,M,LWork,Trans)
#include "WrkSpc.fh"
       INTEGER L, M, LWork
       Real*8,DIMENSION(L,M)::Mat
       CHARACTER(Len=1)::Trans
       INTEGER I,J
       IF(Trans.ne.'T') THEN
        DO I=1,L
         Do J=1,M
          Mat(J,I)=WORK(LWork+J-1+(I-1)*M)
         End Do
        END DO
       ELSE
        DO I=1,L
         Do J=1,M
          Mat(I,J)=WORK(LWork+J-1+(I-1)*M)
         End Do
        END DO
       END IF
       RETURN
       End Subroutine
******************************************************
