************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      PROGRAM ASC2JOB
      use rasscf_global, only: BName, Header, IADR15, iPT2, iRoot,
     &                         IROOT, lRoots, NACPAR, NACPR2, NORBT,
     &                         nRoots, NTOT3, PotNuc, Title, Weight
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      INTEGER FMTIPH, I, IAD15, ISYM, J, lSym, NASHT, NFOCK, nHeader,
     &        nName, nTitle
      INTEGER, EXTERNAL :: isFreeUnit

#include "rasdim.fh"
#include "general.fh"
      REAL*8, ALLOCATABLE:: ADR1(:), ADR2(:), ADR(:)

      CALL INIT_LINALG
      CALL PrgmInit('Asc2Job')

      JOBIPH=15
      JOBIPH=isFreeUnit(JOBIPH)
      CALL DANAME(JOBIPH,'JOBIPH')

      FMTIPH=JOBIPH+1
      FMTIPH=isFreeUnit(FMTIPH)
      CALL MOLCAS_OPEN(FMTIPH,'FMTIPH')

      READ(FMTIPH,*)
      READ(FMTIPH,'(15I10)') (IADR15(I), I=1,30)
      IAD15=0
      CALL IDAFILE(JOBIPH,1,IADR15,30,IAD15)

      nName=LENIN8*mxOrb
      nHeader=144
      nTitle=4*18*mxTit

 100  FORMAT(I10)
 200  FORMAT(2X,10I10)
 300  FORMAT(A)
 400  FORMAT(E20.12)
 500  FORMAT(2X,5E20.12)
 510  FORMAT(2X,6E20.12)
      READ(FMTIPH,*)
      READ(FMTIPH,100) nActEl
      READ(FMTIPH,*)
      READ(FMTIPH,100) iSpin
      READ(FMTIPH,*)
      READ(FMTIPH,100) nSym
      READ(FMTIPH,*)
      READ(FMTIPH,100) lSym
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nFro(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nISh(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nASh(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nDel(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nBas(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,100) nConf
      READ(FMTIPH,*)
      READ(FMTIPH,'(36A)') (Header(I), I=1,72)
      READ(FMTIPH,*)
      READ(FMTIPH,'(80A)') (Title(I), I=1,9)
      READ(FMTIPH,*)
      READ(FMTIPH,400) PotNuc
      READ(FMTIPH,*)
      READ(FMTIPH,100) lRoots
      READ(FMTIPH,*)
      READ(FMTIPH,100) nRoots
      READ(FMTIPH,*)
      READ(FMTIPH,200) (iRoot(I), I=1,MxRoot)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nRS1(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nRS2(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,200) (nRS3(I), I=1,MxSym)
      READ(FMTIPH,*)
      READ(FMTIPH,100) nHole1
      READ(FMTIPH,*)
      READ(FMTIPH,100) nElec3
      READ(FMTIPH,*)
      READ(FMTIPH,100) iPT2
      READ(FMTIPH,*)
      READ(FMTIPH,500) (Weight(I), I=1,MxRoot)

      NTOT=0
      NTOT2=0
      NASHT=0
      NTOT3=0
      NFOCK=0
      DO ISYM=1,NSYM
         NTOT=NTOT+nBas(ISYM)
         NTOT2=NTOT2+nBas(ISYM)**2
         NASHT=NASHT+nASh(ISYM)
         nOrb(iSym)=nBas(iSym)-nFro(iSym)-nDel(iSym)
         NTOT3=NTOT3+nOrb(iSym)*(nOrb(iSym)+1)/2
         NFOCK=NFOCK+(nISh(ISYM)+nASh(ISYM))**2
      END DO
      NACPAR=NASHT*(NASHT+1)/2
      NACPR2=NACPAR*(NACPAR+1)/2

      READ(FMTIPH,*)
      READ(FMTIPH,300) (BName(I), I=1,NTOT)

      IAD15=IADR15(1)
      CALL WR_RASSCF_Info(JOBIPH,1,IAD15,NACTEL,ISPIN,NSYM,LSYM,
     &            NFRO,NISH,NASH,NDEL,NBAS,MxSym,
     &            BName,nName,NCONF,HEADER,nHeader,
     &            TITLE,nTitle,POTNUC,LROOTS,NROOTS,
     &            IROOT,MxRoot,NRS1,NRS2,NRS3,
     &            NHOLE1,NELEC3,IPT2,WEIGHT)

      IAD15=IADR15(2)
      Call mma_allocate(ADR1,NTOT2,Label='ADR1')
      Call mma_allocate(ADR2,NTOT ,Label='ADR2')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR1(I), I=1,NTOT2)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR2(I), I=1,NTOT)
      CALL DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
      CALL DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
      Call mma_deallocate(ADR1)
      Call mma_deallocate(ADR2)

      IAD15=IADR15(3)
      Call mma_allocate(ADR1,NACPAR,Label='ADR1')
      Call mma_allocate(ADR2,NACPR2,Label='ADR2')
      READ(FMTIPH,*)
      Do I = 1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR1(J), J=1,NACPAR)
        Call DDafIle(JOBIPH,1,ADR1,NACPAR,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR1(J), J=1,NACPAR)
        Call DDafIle(JOBIPH,1,ADR1,NACPAR,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR2(J), J=1,NACPR2)
        Call DDafIle(JOBIPH,1,ADR2,NACPR2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR2(J), J=1,NACPR2)
        Call DDafIle(JOBIPH,1,ADR2,NACPR2,IAD15)
      ENDDO
      Call mma_deallocate(ADR1)
      Call mma_deallocate(ADR2)

      IAD15=IADR15(4)
      Call mma_allocate(ADR,NCONF,LABEL='ADR')
      READ(FMTIPH,*)
      Do I = 1,LROOTS
      READ(FMTIPH,*)
      READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR(J), J=1,NCONF)
        Call DDafIle(JOBIPH,1,ADR,NCONF,IAD15)
      ENDDO
      Call mma_deallocate(ADR)

      IAD15=IADR15(5)
      Call mma_allocate(ADR,NFOCK,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,NFOCK)
      CALL DDAFILE(JOBIPH,1,ADR,NFOCK,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(6)
      Call mma_allocate(ADR,mxRoot*mxIter,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,mxRoot*mxIter)
      CALL DDAFILE(JOBIPH,1,ADR,mxRoot*mxIter,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(7)
      Call mma_allocate(ADR,6*mxIter,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,510) (ADR(I), I=1,6*mxIter)
      CALL DDAFILE(JOBIPH,1,ADR,6*mxIter,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(9)
      Call mma_allocate(ADR,NTOT2,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,NTOT2)
      CALL DDAFILE(JOBIPH,1,ADR,NTOT2,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(10)
      Call mma_allocate(ADR,NTOT3,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,NTOT3)
      CALL DDAFILE(JOBIPH,1,ADR,NTOT3,IAD15)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,NTOT3)
      CALL DDAFILE(JOBIPH,1,ADR,NTOT3,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(11)
      Call mma_allocate(ADR,NORBT,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(I), I=1,NORBT)
      CALL DDAFILE(JOBIPH,1,ADR,NORBT,IAD15)
      Call mma_deallocate(ADR)

      IAD15=IADR15(12)
      Call mma_allocate(ADR1,NTOT2,Label='ADR1')
      Call mma_allocate(ADR2,NTOT ,Label='ADR2')
      READ(FMTIPH,*)
      DO I = 1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR1(J), J=1,NTOT2)
        CALL DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR2(J), J=1,NTOT)
        CALL DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
      ENDDO
      Call mma_deallocate(ADR1)
      Call mma_deallocate(ADR2)

      IAD15=IADR15(14)
      Call mma_allocate(ADR1,NTOT2,Label='ADR1')
      Call mma_allocate(ADR2,NTOT ,Label='ADR2')
      READ(FMTIPH,*)
      Do I=1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR1(J), J=1,NTOT2)
        CALL DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (ADR2(J), J=1,NTOT)
        CALL DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
      End Do
      Call mma_deallocate(ADR1)
      Call mma_deallocate(ADR2)

      IAD15=IADR15(17)
      Call mma_allocate(ADR,LROOTS**2,LABEL='ADR')
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (ADR(J), J=1,LROOTS**2)
      CALL DDAFILE(JOBIPH,1,ADR,LROOTS**2,IAD15)
      Call mma_deallocate(ADR)

      CALL DACLOS(JOBIPH)
      CLOSE(FMTIPH)

      END
