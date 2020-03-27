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
      SUBROUTINE PROPCT()
      IMPLICIT REAL*8 (A-H,O-Z)
*PAM04      DIMENSION HWork(*)
      CHARACTER*8 FNAME,LABEL
c      CHARACTER*8 REMARK
      CHARACTER*30 REMARK
      CHARACTER*100 REALNAME
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      DIMENSION DUMMY(1),IDUMMY(7,8)
      REAL*8, ALLOCATABLE :: PROP(:,:,:)
*
      CALL QENTER('PROPCT')
*      LCMO=LPRP
*      LCNO=LCMO+NCMO
*      LOCC=LCNO+NCMO
*      LTDAO=LOCC+NBAST
*      LDAO=LTDAO
*      LAFOLD=LDAO+NBAST**2
*      LSFOLD=LAFOLD+NBTRI
*      LPINT =LSFOLD+NBTRI
*      LSCR  =LPINT+NBTRI+4
*      NSCR=MAX(NBTRI,NBMAX**2)
*      LTOP  =LSCR+NSCR-1
      CALL GETMEM('CMO','ALLO','REAL',LCMO,NCMO)
      CALL GETMEM('CNO','ALLO','REAL',LCNO,NCMO)
      CALL GETMEM('OCC','ALLO','REAL',LOCC,NBAST)
      CALL GETMEM('DAO','ALLO','REAL',LDAO,NBAST**2)
      LTDAO=LDAO
      CALL GETMEM('AFOLD','ALLO','REAL',LAFOLD,NBTRI)
      CALL GETMEM('SFOLD','ALLO','REAL',LSFOLD,NBTRI)
      CALL GETMEM('PINT','ALLO','REAL',LPINT,NBTRI+4)
      NSCR=MAX(NBTRI,NBMAX**2)
      CALL GETMEM('SCR','ALLO','REAL',LSCR,NSCR)
C LOOP OVER OPERATORS:
      IOPT=8
      NPROP=0
      DO 98 I=1,100
C PICK UP OPERATOR LABELS FROM ONE-ELECTRON FILE:
        LABEL='UNDEF'
        CALL iRDONE(IRTC,1+IOPT,LABEL,IPC,IDUMMY,ISYMLB)
        IF(IRTC.NE.0) GOTO 99
        IOPT=16
CPAM96        IF(IAND(1,ISYMLB).EQ.0) GOTO 98
        IF(MOD(ISYMLB,2).EQ.0) GOTO 98
        NPROP=NPROP+1
        PNAME(NPROP)=LABEL
        IPCOMP(NPROP)=IPC
        PTYPE(NPROP)='HERM'
        IF(LABEL.EQ.'VELOCITY') PTYPE(NPROP)='ANTI'
        IF(LABEL.EQ.'ANGMOM  ') PTYPE(NPROP)='ANTI'
98    CONTINUE
99    CONTINUE
      CALL MMA_ALLOCATE(PROP,NRROOT,NRROOT,NPROP,'PROP')
      CALL DCOPY_(NRROOT*NRROOT*NPROP,[0.0D00],0,PROP,1)
      IDDMO=0
      DO 100 ISTATE=1,NRROOT
C PICK UP DMO
*PAM04        CALL dDAFILE(LUEIG,2,HWork(LDMO),NBTRI,IDDMO)
        CALL dDAFILE(LUEIG,2,Work(LDMO),NBTRI,IDDMO)
C PICK UP CMO
        IDISK=ITOC17(1)
        CALL dDAFILE(LUONE,2,Work(LCMO),NCMO,IDISK)
C COMPUTE & WRITE NATURAL ORBITALS
*PAM04        CALL NATORB_MRCI(HWork(LCMO),HWork(LDMO),HWork(LCNO),
        CALL NATORB_MRCI(Work(LCMO),Work(LDMO),Work(LCNO),
     *                   Work(LOCC),Work(LSCR))
        WRITE(FNAME,'(A5,I2.2)')'CIORB',ISTATE
        REALNAME=FNAME
*PAM04        IF(ISTATE.EQ.1) Call Add_Info('CI_DENS1',HWork(LDMO),1,5)
        IF(ISTATE.EQ.1) Call Add_Info('CI_DENS1',Work(LDMO),1,5)
        REMARK='* MRCI  '
c** Gusarov , include 1st root acpf energy to CiOrb file:
c         IF(ICPF.EQ.1) REMARK='* ACPF  '
        IF(ICPF.EQ.1)
     &  write(REMARK,'(8H* ACPF  ,f22.16)') ESMALL(1)+ESHIFT
c
        CALL WRVEC(REALNAME,LUVEC,'CO',NSYM,NBAS,NBAS,
     &  Work(LCNO), Work(LOCC), Dummy, iDummy,REMARK)
        WRITE(6,*)
        WRITE(6,'(A,I2)')' NATURAL ORBITALS OF STATE NR. ',ISTATE
        WRITE(6,*)' FULL SET OF ORBITALS ARE SAVED ON FILE ',REALNAME
        CALL PRORB(Work(LCNO),Work(LOCC))
        WRITE(6,*)' ',('*',I=1,70)
C CREATE DAO
        CALL MKDAO(Work(LCNO),Work(LOCC),Work(LDAO))
C CALL PMATEL TO CALCULATE CHARGES AND PROPERTIES.
C PUT PROPERTIES INTO APPROPRIATE MATRICES.
        CALL PMATEL (ISTATE,ISTATE,PROP,Work(LPINT),Work(LSCR),
     *               Work(LCNO),Work(LOCC),
     *               Work(LSFOLD),Work(LAFOLD),Work(LDAO))
100   CONTINUE
C ENERGIES SAVED FROM PREVIOUS OUTPUT REPEATED HERE FOR CONVENIENCE:
      WRITE(6,*)
      WRITE(6,*)' SUMMARY OF ENERGIES:'
      DO 987 ISTA=1,NRROOT,4
        IEND=MIN(ISTA+3,NRROOT)
        WRITE(6,'(1X,A,I8,3I16)')
     *  '               ROOT:',(I,I=ISTA,IEND)
        WRITE(6,'(1X,A,4F16.8)')
     *  '       TOTAL ENERGY:',(ENGY(I,1),I=ISTA,IEND)
        IF(ICPF.EQ.0) THEN
          WRITE(6,'(1X,A,4F16.8)')
     *  'DAVIDSON CORRECTION:',(ENGY(I,2),I=ISTA,IEND)
          WRITE(6,'(1X,A,4F16.8)')
     *  '    ACPF CORRECTION:',(ENGY(I,3),I=ISTA,IEND)
      END IF
      WRITE(6,*)
987   CONTINUE
* ---------------------------------------------------
*PAM Grep-able energy output for convenience:
      WRITE(6,*)
      WRITE(6,*)' Energies, machine-readable format:'
      IF (ICPF.EQ.0) THEN
        DO I=1,NRROOT
         WRITE(6,'(1X,A,I3,3(5X,A,F16.8))')
     &     ' CI State ',I,'Total energy:',ENGY(I,1),
     &     'QDav:',ENGY(I,2),'QACPF:',ENGY(I,3)
        END DO
      ELSE
        DO I=1,NRROOT
         WRITE(6,'(1X,A,I3,5X,A,F16.8)')
     &     ' ACPF State ',I,'Total energy:',ENGY(I,1)
        END DO
      END IF
      WRITE(6,*)
* ---------------------------------------------------
      IF(NPROP.GT.0) THEN
C WRITE EXPECTATION VALUES:
        WRITE(6,*)
        WRITE(6,*)' EXPECTATION VALUES OF VARIOUS OPERATORS:'
        WRITE(6,*)
     *  '(Note: Electronic multipoles include a negative sign.)'
        DO 110 IPROP=1,NPROP
          IF(PTYPE(IPROP).EQ.'ANTI') GOTO 110
          DO 105 ISTA=1,NRROOT,4
            IEND=MIN(ISTA+3,NRROOT)
            WRITE(6,*)
            WRITE(6,'(1X,A,A8,A,I4)')
     *    '   PROPERTY :',PNAME(IPROP),
     *    '   COMPONENT:',IPCOMP(IPROP)
            WRITE(6,'(1X,A,3F16.8)')
     *    '    GAUGE ORIGIN:',(PORIG(I,IPROP),I=1,3)
            WRITE(6,'(1X,A,I8,3I16)')
     *    '            ROOT:',(I,I=ISTA,IEND)
            WRITE(6,'(1X,A,4F16.8)')
     *    '      ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
            WRITE(6,'(1X,A,4F16.8)')
     *    '         NUCLEAR:',(PNUC(IPROP),I=ISTA,IEND)
            WRITE(6,'(1X,A,4F16.8)')
     *    '           TOTAL:',(PNUC(IPROP)+PROP(I,I,IPROP),I=ISTA,IEND)
105       CONTINUE
110     CONTINUE
        WRITE(6,*)
      END IF
      IF(ITRANS.EQ.0) GOTO 1000
      DO 200 ISTATE=2,NRROOT
        DO 200 JSTATE=1,ISTATE-1
C PICK UP TDMA
*PAM04        CALL dDAFILE(LUEIG,2,HWork(LTDMO),NBAST**2,IDDMO)
        CALL dDAFILE(LUEIG,2,Work(LTDMO),NBAST**2,IDDMO)
C CREATE TDAO
*PAM04        CALL MKTDAO(HWork(LCMO),HWork(LTDMO),HWork(LTDAO),HWork(LSCR))
        CALL MKTDAO(Work(LCMO),Work(LTDMO),Work(LTDAO),Work(LSCR))
C CALL PMATEL TO CALCULATE TRANSITION PROPERTIES
C PUT PROPERTIES INTO APPROPRIATE MATRICES.
        IF(NPROP.EQ.0) GOTO 200
        CALL PMATEL (ISTATE,JSTATE,PROP,Work(LPINT),Work(LSCR),
     *               Work(LCNO),Work(LOCC),
     *               Work(LSFOLD),Work(LAFOLD),Work(LTDAO))
200   CONTINUE
      IF(NPROP.EQ.0) GOTO 1000
C WRITE PROPERTY MATRICES.
      WRITE(6,*)
      WRITE(6,*)' MATRIX ELEMENTS OF VARIOUS OPERATORS:'
      WRITE(6,*)' (INCLUDING ANY NUCLEAR CONTRIBUTIONS)'
      DO 299 IPROP=1,NPROP
        DO 299 I=1,NRROOT
          PROP(I,I,IPROP)=PROP(I,I,IPROP)+PNUC(IPROP)
299   CONTINUE
      DO 300 IPROP=1,NPROP
        DO 300 ISTA=1,NRROOT,4
          IEND=MIN(ISTA+3,NRROOT)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     *    '   PROPERTY :',PNAME(IPROP),
     *    '   COMPONENT:',IPCOMP(IPROP)
          WRITE(6,'(1X,A,3F16.8)')
     *    '    GAUGE ORIGIN:',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,3I16)')
     *    '            ROOT:',(I,I=ISTA,IEND)
          DO 310 J=1,NRROOT
            WRITE(6,'(15X,I2,4F16.8)')
     *        J,(PROP(J,I,IPROP),I=ISTA,IEND)
310       CONTINUE
300   CONTINUE
      WRITE(6,*)
1000  CONTINUE
      CALL GETMEM('CMO','FREE','REAL',LCMO,NCMO)
      CALL GETMEM('CNO','FREE','REAL',LCNO,NCMO)
      CALL GETMEM('OCC','FREE','REAL',LOCC,NBAST)
      CALL GETMEM('DAO','FREE','REAL',LDAO,NBAST**2)
      CALL GETMEM('AFOLD','FREE','REAL',LAFOLD,NBTRI)
      CALL GETMEM('SFOLD','FREE','REAL',LSFOLD,NBTRI)
      CALL GETMEM('PINT','FREE','REAL',LPINT,NBTRI+4)
      CALL GETMEM('SCR','FREE','REAL',LSCR,NSCR)
      CALL MMA_DEALLOCATE(PROP)
      CALL QEXIT('PROPCT')
      RETURN
      END
