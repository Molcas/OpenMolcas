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
      SUBROUTINE TRAONE(CMO)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      DIMENSION CMO(NCMO)
      DIMENSION nBasXX(8),Keep(8)
      Logical iSquar, Found

c Objective: Transformation of one-electron integrals
c (effective one electron Hamiltonian) for CASPT2.

      CALL QENTER('TRAONE')
#ifdef _DEBUG_
      IFTEST=1
#else
      IFTEST=0
#endif

      Call GetOrd(IRC,iSquar,nSymXX,nBasXX,Keep)
      IF ( IPRGLB.GE.VERBOSE ) THEN
        If(iSquar)      WRITE(6,*) 'TRAONE OrdInt status: squared'
        If(.not.iSquar) WRITE(6,*) 'TRAONE OrdInt status: non-squared'
      ENDIF
      IERR=0
      DO ISYM=1,NSYM
        IF (NBAS(ISYM).NE.NBASXX(ISYM)) IERR=1
      END DO
      IF(IERR.NE.0) THEN
        WRITE(6,*)'     *** ERROR IN SUBROUTINE TRAONE ***'
        WRITE(6,*)'          INCOMPATIBLE BASIS DATA'
        WRITE(6,*)
        WRITE(6,*)' JOBIPH NR OF SYMM:', NSYM
        WRITE(6,*)' JOBIPH NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
        WRITE(6,*)
        WRITE(6,*)' ORDINT NR OF SYMM:', NSYMXX
        WRITE(6,*)' ORDINT NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBASXX(I),I=1,NSYMXX)
        CALL ERRTRA
        CALL ABEND()
      END IF
c Allocate FLT,DLT, and DSQ.
      CALL GETMEM('WFLT','ALLO','REAL',LWFLT,NBTRI)
c Read nuclear repulsion energy:
      IRC=-1
      IOPT=0
      ICOMP=0
      ISYLBL=1
      IF ( IFTEST.NE.0 ) WRITE(6,*)' GET POTNUC FROM RUNFILE'
      Call Get_dScalar('PotNuc',PotNuc)
      IF ( IFTEST.NE.0 ) WRITE(6,*)' POTNUC:',POTNUC
c Read one-electron hamiltonian matrix into FLT.
      IRC=-1
      IOPT=6
      ICOMP=1
      ISYLBL=1
      IF ( IFTEST.NE.0 ) WRITE(6,*)' CALLING RDONE (ONEHAM)'
      CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,WORK(LWFLT),ISYLBL)
      IF ( IFTEST.NE.0 ) WRITE(6,*)' BACK FROM RDONE'
      IF(IRC.NE.0) THEN
        WRITE(6,*)'TRAONE Error: RDONE failed reading OneHam.'
        Call Quit(_RC_IO_ERROR_READ_)
      END IF

      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)'     TEST PRINTS FROM TRAONE.'
        WRITE(6,*)'     NAKED 1-EL HAMILTONIAN IN AO BASIS'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NBAS(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',WORK(LWFLT+ISTLT),NBAS(ISYM))
            ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
          END IF
        END DO
      END IF

c If this is a perturbative reaction field calculation then
c modifiy the one-electron Hamiltonian by the reaction field and
c the nuclear attraction by the cavity self-energy

      If ( RFpert ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
*
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call Get_dScalar('RF Self Energy',ERFSelf)
         Call Get_dArray('Reaction field',Work(lTemp),nTemp)
         If (Found) Call NameRun('RUNFILE')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nTemp,1.0D0,Work(lTemp),1,WORK(LWFLT),1)
*
         Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
         IF ( IFTEST.NE.0 ) THEN
           WRITE(6,*)' 1-EL HAMILTONIAN INCLUDING REACTION FIELD'
           ISTLT=0
           DO ISYM=1,NSYM
             IF ( NBAS(ISYM).GT.0 ) THEN
               WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
               CALL TRIPRT(' ',' ',WORK(LWFLT+ISTLT),NBAS(ISYM))
               ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
             END IF
           END DO
         END IF
      End If

      EONE=0.0d0
      ETWO=0.0d0
c The following section is needed for frozen orbitals:
      IF(NFROT.EQ.0) GOTO 300
      CALL GETMEM('WDLT','ALLO','REAL',LWDLT,NBTRI)
      CALL GETMEM('WDSQ','ALLO','REAL',LWDSQ,NBSQT)
c Compute the density matrix of the frozen orbitals
c The DLT matrix contains the same data as DSQ, but
c with symmetry blocks in lower triangular format, and
c with non-diagonal elements doubled.
      CALL DCOPY_(NBTRI,[0.0D0],0,WORK(LWDLT),1)
      CALL DCOPY_(NBSQT,[0.0D0],0,WORK(LWDSQ),1)
      ISTMO=1
      ISTSQ=LWDSQ
      ISTLT=LWDLT
      DO 100 ISYM=1,NSYM
        NF=NFRO(ISYM)
        NB=NBAS(ISYM)
        IF(NB.EQ.0) GOTO 100
        IF(NF.EQ.0) GOTO 110
        CALL DGEMM_('N','T',NB,NB,NF,2.0D0,CMO(ISTMO),NB,
     &             CMO(ISTMO),NB,0.0D0,WORK(ISTSQ),NB)
        IJ=ISTLT-1
        DO 130 IB=1,NB
          DO 140 JB=1,IB
            IJ=IJ+1
            WORK(IJ)=2.0D0*WORK(ISTSQ+JB-1+(IB-1)*NB)
140       CONTINUE
          WORK(IJ)=0.5D0*WORK(IJ)
130     CONTINUE
110     CONTINUE
        ISTMO=ISTMO+NB*NB
        ISTSQ=ISTSQ+NB*NB
        ISTLT=ISTLT+NB*(NB+1)/2
100   CONTINUE

c One-electron contribution to the core energy.
c Note that FLT still contains only the naked
c  one-electron hamiltonian.
      EONE=DDOT_(NBTRI,WORK(LWDLT),1,WORK(LWFLT),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate Fock-matrix for frozen orbitals
*     and compute the total core energy
*     Look out-- we temporarily allocate all available memory.
*
      ExFac=1.0D0

         Call FTwo_Drv(nSym,nBas,nFro,KEEP,
     &                 WORK(LWDLT),WORK(LWDSQ),WORK(LWFLT),NBTRI,
     &                 ExFac,nBSQT,nBMX,CMO)

*                                                                      *
************************************************************************
*                                                                      *
c Compute the two-electron contribution to the core energy
      ETWO=0.5D0*(DDOT_(NBTRI,WORK(LWDLT),1,WORK(LWFLT),1)-EONE)
      CALL GETMEM('WDSQ','FREE','REAL',LWDSQ,NBSQT)
      CALL GETMEM('WDLT','FREE','REAL',LWDLT,NBTRI)
c Previous section was bypassed if NFROT.EQ.0.
 300  CONTINUE
      ECORE=POTNUC+EONE+ETWO
      IF ( IFTEST.NE.0 ) THEN
         WRITE(6,'(6X,A,E20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
         WRITE(6,'(6X,A,E20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
         WRITE(6,'(6X,A,E20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
         WRITE(6,'(6X,A,E20.10)') '       TOTAL CORE ENERGY:',ECORE
      ENDIF

c Allocate FMO, TMP:
      NWTMP=2*NBMX**2
      CALL GETMEM('WFMO','ALLO','REAL',LWFMO,notri)
      CALL GETMEM('WTMP','ALLO','REAL',LWTMP,NWTMP)

c Transform one-electron effective Hamiltonian:
      CALL DCOPY_(notri,[0.0D0],0,WORK(LWFMO),1)
      CALL DCOPY_(NWTMP,[0.0D0],0,WORK(LWTMP),1)
      ICMO=1
      IAO =LWFLT
      IMO =LWFMO
      DO 200 ISYM=1,NSYM
         ICMO=ICMO+NBAS(ISYM)*NFRO(ISYM)
         IOFF=LWTMP+NBAS(ISYM)*NBAS(ISYM)
         IF(NORB(ISYM).GT.0) THEN
           CALL SQUARE(WORK(IAO),WORK(LWTMP),1,NBAS(ISYM),NBAS(ISYM))

           CALL DGEMM_('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),
     &                  1.0d0,CMO(ICMO),NBAS(ISYM),WORK(LWTMP),
     &                  NBAS(ISYM),0.0d0,WORK(IOFF),NORB(ISYM))

           CALL MXMT(WORK(IOFF),    1,NORB(ISYM),
     &             CMO(ICMO),     1,NBAS(ISYM),
     &             WORK(IMO),
     &             NORB(ISYM),NBAS(ISYM))
         END IF
         ICMO=ICMO+NBAS(ISYM)*(NORB(ISYM)+NDEL(ISYM))
         IAO =IAO +NBAS(ISYM)*(NBAS(ISYM)+1)/2
         IMO =IMO +NORB(ISYM)*(NORB(ISYM)+1)/2
200   CONTINUE

      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)'      EFFECTIVE 1-EL HAMILTONIAN IN MO BASIS'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',WORK(LWFMO+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
      END IF
      IDISK=IEOF1M
      IAD1M(3)=IDISK
      CALL DDAFILE(LUONEM,1,WORK(LWFMO),notri,IDISK)
      IEOF1M=IDISK
      CALL DCOPY_(NOTRI,WORK(LWFMO),1,WORK(LHONE),1)
      CALL GETMEM('WTMP','FREE','REAL',LWTMP,NWTMP)
      CALL GETMEM('WFMO','FREE','REAL',LWFMO,notri)
      CALL GETMEM('WFLT','FREE','REAL',LWFLT,NBTRI)

      CALL QEXIT('TRAONE')

      RETURN
      End
