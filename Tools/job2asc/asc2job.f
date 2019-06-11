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
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER FMTIPH

#include "rasdim.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "general.fh"

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
      READ(FMTIPH,300) (Name(I), I=1,NTOT)

      IAD15=IADR15(1)
      CALL WR_RASSCF_Info(JOBIPH,1,IAD15,NACTEL,ISPIN,NSYM,LSYM,
     &            NFRO,NISH,NASH,NDEL,NBAS,MxSym,
     &            NAME,nName,NCONF,HEADER,nHeader,
     &            TITLE,nTitle,POTNUC,LROOTS,NROOTS,
     &            IROOT,MxRoot,NRS1,NRS2,NRS3,
     &            NHOLE1,NELEC3,IPT2,WEIGHT)

      IAD15=IADR15(2)
      CALL GETMEM('ADR1','ALLO','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','ALLO','REAL',LADR2,NTOT)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR1+I-1), I=1,NTOT2)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR2+I-1), I=1,NTOT)
      CALL DDAFILE(JOBIPH,1,WORK(LADR1),NTOT2,IAD15)
      CALL DDAFILE(JOBIPH,1,WORK(LADR2),NTOT,IAD15)
      CALL GETMEM('ADR1','FREE','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','FREE','REAL',LADR2,NTOT)

      IAD15=IADR15(3)
      CALL GETMEM('ADR1','ALLO','REAL',LADR1,NACPAR)
      CALL GETMEM('ADR2','ALLO','REAL',LADR2,NACPR2)
      READ(FMTIPH,*)
      Do I = 1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR1+J-1), J=1,NACPAR)
        Call DDafIle(JOBIPH,1,WORK(LADR1),NACPAR,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR1+J-1), J=1,NACPAR)
        Call DDafIle(JOBIPH,1,WORK(LADR1),NACPAR,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR2+J-1), J=1,NACPR2)
        Call DDafIle(JOBIPH,1,WORK(LADR2),NACPR2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR2+J-1), J=1,NACPR2)
        Call DDafIle(JOBIPH,1,WORK(LADR2),NACPR2,IAD15)
      ENDDO
      CALL GETMEM('ADR1','FREE','REAL',LADR1,NACPAR)
      CALL GETMEM('ADR2','FREE','REAL',LADR2,NACPR2)

      IAD15=IADR15(4)
      CALL GETMEM('ADR','ALLO','REAL',LADR,NCONF)
      READ(FMTIPH,*)
      Do I = 1,LROOTS
      READ(FMTIPH,*)
      READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR+J-1), J=1,NCONF)
        Call DDafIle(JOBIPH,1,WORK(LADR),NCONF,IAD15)
      ENDDO
      CALL GETMEM('ADR','FREE','REAL',LADR,NCONF)

      IAD15=IADR15(5)
      CALL GETMEM('ADR','ALLO','REAL',LADR,NFOCK)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,NFOCK)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),NFOCK,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,NFOCK)

      IAD15=IADR15(6)
      CALL GETMEM('ADR','ALLO','REAL',LADR,mxRoot*mxIter)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,mxRoot*mxIter)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),mxRoot*mxIter,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,mxRoot*mxIter)

      IAD15=IADR15(7)
      CALL GETMEM('ADR','ALLO','REAL',LADR,6*mxIter)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,510) (WORK(LADR+I-1), I=1,6*mxIter)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),6*mxIter,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,6*mxIter)

      IAD15=IADR15(9)
      CALL GETMEM('ADR','ALLO','REAL',LADR,NTOT2)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,NTOT2)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),NTOT2,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,NTOT2)

      IAD15=IADR15(10)
      CALL GETMEM('ADR','ALLO','REAL',LADR,NTOT3)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,NTOT3)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),NTOT3,IAD15)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,NTOT3)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),NTOT3,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,NTOT3)

      IAD15=IADR15(11)
      CALL GETMEM('ADR','ALLO','REAL',LADR,NORBT)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+I-1), I=1,NORBT)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),NORBT,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,NORBT)

      IAD15=IADR15(12)
      CALL GETMEM('ADR1','ALLO','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','ALLO','REAL',LADR2,NTOT)
      READ(FMTIPH,*)
      DO I = 1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR1+J-1), J=1,NTOT2)
        CALL DDAFILE(JOBIPH,1,WORK(LADR1),NTOT2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR2+J-1), J=1,NTOT)
        CALL DDAFILE(JOBIPH,1,WORK(LADR2),NTOT,IAD15)
      ENDDO
      CALL GETMEM('ADR1','FREE','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','FREE','REAL',LADR2,NTOT)

      IAD15=IADR15(14)
      CALL GETMEM('ADR1','ALLO','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','ALLO','REAL',LADR2,NTOT)
      READ(FMTIPH,*)
      Do I=1,LROOTS
      READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR1+J-1), J=1,NTOT2)
        CALL DDAFILE(JOBIPH,1,WORK(LADR1),NTOT2,IAD15)
        READ(FMTIPH,*)
        READ(FMTIPH,*)
        READ(FMTIPH,500) (WORK(LADR2+J-1), J=1,NTOT)
        CALL DDAFILE(JOBIPH,1,WORK(LADR2),NTOT,IAD15)
      End Do
      CALL GETMEM('ADR1','FREE','REAL',LADR1,NTOT2)
      CALL GETMEM('ADR2','FREE','REAL',LADR2,NTOT)

      IAD15=IADR15(17)
      CALL GETMEM('ADR','ALLO','REAL',LADR,LROOTS**2)
      READ(FMTIPH,*)
      READ(FMTIPH,*)
      READ(FMTIPH,500) (WORK(LADR+J-1), J=1,LROOTS**2)
      CALL DDAFILE(JOBIPH,1,WORK(LADR),LROOTS**2,IAD15)
      CALL GETMEM('ADR','FREE','REAL',LADR,LROOTS**2)

      CALL DACLOS(JOBIPH)
      CLOSE(FMTIPH)

      END
