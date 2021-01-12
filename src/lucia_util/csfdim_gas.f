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
      SUBROUTINE CSFDIM_GAS(IOCCLS,NOCCLS,ISYM,IPRCSF)
*
* Initializing routine for CSF-DET expansions
*
* information about the number of dets,csf's for
* each symmetry. CI space is defined by the NOCCLS
* occupation classes IOCCLS
*
*
* DETERMINE BASE ADRESSES
*             KDFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
*             KCFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
*             KDTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
*             KICONF_OCC(I) : SPACE FOR STORING  NCNSM
*                        CONFIGURATION EXPANSIONS
* Local memory requirements : IWORK(NACTOB)
*
* ( Spin signaled by PSSIGN in CIINFO)
*
* Adapted for GAS calculations and LUCIA, Dec. 2001

#include "implicit.fh"
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "glbbas.fh"
#include "cgas.fh"
#include "spinfo_lucia.fh"
#include "warnings.fh"
#include "gasstr.fh"
* Input type of occupation classes
      INTEGER IOCCLS(NGAS,*)
      INTEGER IDUM_ARR(1)
      INTEGER TMP_CNF(MXPORB+1),HEXS_CNF(MXPORB+1),
     &        maxingas(N_ELIMINATED_GAS)
*
      IDUM = 0
      IDUM_ARR=0
*
      NTEST = 00
      NTEST = MAX(IPRCSF,NTEST)
      IF(NTEST.GE.10) WRITE(6,*) '  PSSIGN : ', PSSIGN
      IF(NTEST.GE.10) WRITE(6,*) ' MULTS, MS2 = ', MULTS,MS2
      NELEC = IELSUM(IOCCLS(1,1),NGAS)
*
*.. Define parameters in SPINFO
*
*. Allowed number of open orbitals
      MINOP = ABS(MS2)
      CALL MAX_OPEN_ORB(MAXOP,IOCCLS,NGAS,NOCCLS,NOBPT)
      IF( NTEST .GE. 6 )
     &WRITE(6,*) ' MINOP MAXOP ',MINOP,MAXOP
C
C.. Number of prototype sd's and csf's per configuration prototype
C
           ITP = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
*. Unpaired electrons :
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN .AND. IAEL-IBEL .EQ. MS2 .AND.
     &            IAEL .GE. 0 .AND. IBEL .GE. 0) THEN
          IF(PSSIGN.EQ. 0.0D0 .OR. IOPEN .EQ. 0 ) THEN
*. Number of determinants is in general set to number of combinations
            NPDTCNF(ITP) = IBION_LUCIA(IOPEN,IAEL)
            NPCMCNF(ITP) = NPDTCNF(ITP)
          ELSE
            NPDTCNF(ITP) = IBION_LUCIA(IOPEN,IAEL)/2
            NPCMCNF(ITP) = NPDTCNF(ITP)
          END IF
          IF(IOPEN .GE. MULTS-1) THEN
            NPCSCNF(ITP) = IWEYLF(IOPEN,MULTS)
          ELSE
            NPCSCNF(ITP) = 0
          END IF
        ELSE
          NPDTCNF(ITP) = 0
          NPCMCNF(ITP) = 0
          NPCSCNF(ITP) = 0
        END IF
      END DO
*
      IF(NTEST.GE.5) THEN
      IF(PSSIGN .EQ. 0.0D0 ) THEN
        WRITE(6,*) '  (Combinations = Determinants ) '
      ELSE
        WRITE(6,*) '  (Spin combinations in use ) '
      END IF
      WRITE(6,'(/A)') ' Information about prototype configurations '
      WRITE(6,'( A)') ' ========================================== '
      WRITE(6,'(/A)')
     &'  Open orbitals   Combinations    CSFs '
      DO IOPEN = MINOP,MAXOP,2
        WRITE(6,'(5X,I3,10X,I6,7X,I6)')
     &  IOPEN,NPCMCNF(IOPEN+1),NPCSCNF(IOPEN+1)
      END DO
*
      END IF
C
C.. Number of Configurations per occupation type
C
      If (NOCCLS.gt.MXPCSM) Then
         WRITE(6,*)' A known bug has reoccurred -- It seems that'
         WRITE(6,*)' the named constant MXPCSM must be increased'
         WRITE(6,*)' from its current value MXPCSM=',MXPCSM
         WRITE(6,*)' to AT LEAST NOCCLS=',NOCCLS
         WRITE(6,*)' This parameter is found in the file'
         WRITE(6,*)'  <molcas>/src/lucia_util/mxpdim.fh'
         WRITE(6,*)' Change it. Then ''cd'' to molcas root'
         WRITE(6,*)' directory and give command ''make''.'
         WRITE(6,*)' But this may also be a bug. Please tell the'
         WRITE(6,*)' molcas developers!'
         Call Quit(_RC_INTERNAL_ERROR_)
      End If
*MGD : max occupation in removed GAS spaces
      Do i=1,N_ELIMINATED_GAS
        iGAS=IELIMINATED_IN_GAS(i)
        maxingas(i)=0
        DO JOCCLS = 1, NOCCLS
          maxingas(i)=max(IOCCLS(iGAS,JOCCLS),maxingas(i))
        End Do
      End Do
      Do i=1,maxop+1
        HEXS_CNF(i)=0
        NCONF_PER_OPEN(i,ISYM)=0
      End Do

*
      DO JOCCLS = 1, NOCCLS
        IF(JOCCLS.EQ.1) THEN
          INITIALIZE_CONF_COUNTERS = 1
        ELSE
          INITIALIZE_CONF_COUNTERS = 0
        END IF
*
        IDOREO = 0
        IF(JOCCLS.EQ.1) THEN
          NCONF_ALL_SYM_PREV = 0
        ELSE
          NCONF_ALL_SYM_PREV = NCONF_ALL_SYM
        END IF
        Do i=1,maxop+1
          TMP_CNF(i)=0
        End Do
        CALL GEN_CONF_FOR_OCCLS(IOCCLS(1,JOCCLS),
     &                             IDUM,
     &                          INITIALIZE_CONF_COUNTERS,
     &                             NGAS,   ISYM,  MINOP,  MAXOP,  NSMST,
     &                                1,
*
     &                            NOCOB,
     &                            NOBPT,
     &                          TMP_CNF,
     &                          NCONF_OCCLS,
     &                          IB_CONF_REO,
*
     &                          IB_CONF_OCC,
     &                          IDUM_ARR,IDOREO,IDUM_ARR,NCONF_ALL_SYM,
     &                          idum_arr,
     &                          nconf_tot)
         Do i=1,maxop+1
           NCONF_PER_OPEN(i,ISYM)=NCONF_PER_OPEN(i,ISYM)+TMP_CNF(i)
         End Do
*MGD add to hexs_cnf only if the configuration does not have max occupation
* in the selected GAS space
         If (N_ELIMINATED_GAS.gt.0) Then
           ielim=0
           Do j=1,N_ELIMINATED_GAS
             jGAS=IELIMINATED_IN_GAS(j)
             If (IOCCLS(jGAS,JOCCLS).eq.maxingas(j)) ielim=1
             IF(I_ELIMINATE_GAS.GT.1) THEN
               If (IOCCLS(jGAS,JOCCLS).eq.maxingas(j)-1) ielim=1
             END IF
           End Do
           If (ielim.eq.0) Then
             Do i=1,maxop+1
               HEXS_CNF(i)=HEXS_CNF(i)+TMP_CNF(i)
             End Do
           EndIf
         EndIf
c.. testing
c      write(6,*)'nconf_per_open after first call of gen_conf_for_occls'
c      call iwrtma(nconf_per_open,1,4,1,4)
c
*. NCONF_ALL_SYM is accumulated, so
           IF(JOCCLS.EQ.1) THEN
             NCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = NCONF_ALL_SYM
             IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = 1
           ELSE
* PAM2009: It was discovered that these two arrays could be overrun.
* The arrays are declared in spinfo_lucia.fh, and their dimension
* is MXPCSM, which is set in mxpdim.fh -- both included above.
* So MXPCSM is now increased from 20 to 40 -- if this is not a final
* solution remains to be discovered:
             NCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = NCONF_ALL_SYM
     &   -   NCONF_ALL_SYM_PREV
*
             IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS) =
     &       IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS-1)
     &   +    NCONF_ALL_SYM_FOR_OCCLS(JOCCLS-1)

          END IF
      END DO
*. Number of CSF's in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPCSCNF,
     &                   NCSF)
*. Number of SD's in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPDTCNF,
     &                    NSD)
*. Number of combinations in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPCMCNF,
     &                    NCMB)
*MGD
      nCSF_HEXS=0
      If (N_ELIMINATED_GAS.gt.0) Then
        CALL NCNF_TO_NCOMP(MAXOP,HEXS_CNF,NPCSCNF,
     &                   NCSF_HEXS)
      EndIf
*
      NCSF_PER_SYM(ISYM) = NCSF
      NSD_PER_SYM(ISYM) = NSD
      NCONF_PER_SYM(ISYM) = IELSUM(NCONF_PER_OPEN(1,ISYM),MAXOP+1)
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Number of CSFs  ', NCSF
        WRITE(6,*) ' Number of SDs   ', NSD
        WRITE(6,*) ' Number of Confs ', NCONF_PER_SYM(ISYM)
        WRITE(6,*) ' Number of CMBs  ', NCMB
      END IF
C
*. Total number of configurations and length of configuration list
C     INFO_CONF_LIST(NCONF_PER_OPEN,MAXOP,NEL,
C    &                          LENGTH_LIST,NCONF_TOT,IB_REO,IB_OCC)
      CALL INFO_CONF_LIST(NCONF_PER_OPEN(1,ISYM),
     &                    MAXOP,NELEC,LENGTH_LIST,NCONF_TOT,IB_CONF_REO,
     &                    IB_CONF_OCC)
C.. Permanent and local memory for csf routines
C
C    memory for CSDTMT arrays.
C    Largest block of proto type combinations .
C    Largest number of prototype csf's
C
      LIDT = 0
      LICS = 0
      LDTOC = 0
      MXPTBL = 0
      MXDT = 0
      LCONF = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
        LIDT = LIDT + NPCMCNF(ITP) * IOPEN
        LICS = LICS + NPCSCNF(ITP) * IOPEN
        LDTOC= LDTOC + NPCSCNF(ITP)*NPCMCNF(ITP)
        MXDT =   MAX(MXDT,NPCMCNF(ITP) )
        MXPTBL = MAX(NPCMCNF(ITP)*IOPEN,MXPTBL)
      END DO
C. Memory needed to store ICONF array
      LCONF = 0
      ILCNF = 0
*
C     LDET = NSD
      LLCONF = 0
      ILLCNF = 0
      DO IOPEN = 0, MAXOP
        ITYP = IOPEN + 1
        ICL = ( NELEC-IOPEN)/2
        LLCONF = LLCONF + NCONF_PER_OPEN(ITYP,ISYM)*(IOPEN+ICL)
        ILLCNF = ILLCNF + NCONF_PER_OPEN(ITYP,ISYM)
      END DO
C?    WRITE(6,*) ' MEMORY FOR HOLDING CONFS OF SYM... ',ISYM,LLCONF
      LCONF = MAX(LCONF,LLCONF)
      ILCNF = MAX(ILCNF,ILLCNF)
*
       IF(NTEST.GE.5) THEN
       WRITE(6,'(/A,I8)')
     & '  Memory for holding list of configurations ',LCONF
       WRITE(6,'(/A,I8)')
     & '  Size of CI expansion (combinations)',NSD
       WRITE(6,'(/A,I8)')
     & '  Size of CI expansion (confs)',ILCNF
       END IF
C
C. permanent memory for csf proto type arrays
C
      CALL GETMEM('DFTP  ','ALLO','INTE',KDFTP,LIDT)
      CALL GETMEM('CFTP  ','ALLO','INTE',KCFTP,LICS)
      CALL GETMEM('D_TO_C','ALLO','REAL',KDTOC,LDTOC)
C     CALL MEMADD(KDFTP,LIDT,KFREE,1)
C     CALL MEMADD(KCFTP,LICS,KFREE,1)
C     CALL MEMADD(KDTOC,LDTOC,KFREE,2)
      IF( NTEST .GE.6 ) THEN
        WRITE(6,*) ' MEMORY ALLOCATION IN CSFDIM : '
        WRITE(6,*) ' POINTERS FOR PROTOTYPE INFORMATION '
        WRITE(6,'(A,3I7)') ' KDFTP KCFTP KDTOC ',KDFTP,KCFTP,KDTOC
      END IF
C
C. PERMANENT ARRAYS FOR
C. HOLDING CONFIGURATION EXPANSIONS AND REORDER ARRAYS
C
*. Occupation of configurations
      CALL GETMEM('CNFOCC','ALLO','INTE',KICONF_OCC(ISYM),LCONF)
*. Reorder array for configurations
c     CALL MEMMAN(KICONF_REO(ISYM),NCONF_ALL_SYM,'ADDL  ',1,'CNFREO')
      CALL GETMEM('CNFREO','ALLO','INTE',KICONF_REO(ISYM),NCONF_TOT)
*. Reorder array for determinants, index and sign
      CALL GETMEM('SDREOI','ALLO','INTE',KSDREO_I(ISYM),NSD)
*. Sign array is kept real for compatibility
c.. Giovanni and Dongxia decided to delete 'SDREOS' array since it is not used
c     CALL MEMMAN(KSDREO_S(ISYM),NSD,'ADDL  ',1,'SDREOS')
c..
COLD  CALL MEMADD(KICONF(ICNSM),LCONF,KFREE,1)
COLD  CALL MEMADD(KICTS(ICNSM),LDET ,KFREE,1)
COLD  CALL MEMADD(KSCTS(ICNSM),LDET ,KFREE,2)
c     write(6,*)'check me!! dimensions:'
c     write(6,*)'nconf_all_sym =',nconf_all_sym
c     write(6,*)'nconf_tot =', nconf_tot
C
c     ntest=10
      IF( NTEST .GE. 6 ) THEN
        WRITE(6,*) ' KICONF_OCC,  kiconf_reo = ',
     &  KICONF_OCC(ISYM), kiconf_reo(isym)
        WRITE(6,*) ' KSDREO_I',KSDREO_I(ISYM)
      END IF
*
* Arrays for addressing prototype determinants for each number of orbitals
*
      DO IOPEN = MINOP, MAXOP
        ITYP = IOPEN + 1
*
        IALPHA = (IOPEN+MS2)/2
        LZ = IOPEN*IALPHA
        LPTDT = IBION_LUCIA(IOPEN,IALPHA)
        CALL GETMEM('Z_PTDT','ALLO','INTE',KZ_PTDT(ITYP),LZ)
        CALL GETMEM('RE_PTD','ALLO','INTE',KREO_PTDT(ITYP),LPTDT)
      END DO
*
* Array giving first determinant with given number of electrons
* in list of determinants ordered  according to the number of open orbitals
*
      IZERO = 0
      CALL ISETVC(IB_SD_FOR_OPEN,IZERO,MAXOP+1)
      IB = 1
      DO IOPEN = MINOP, MAXOP
        IB_SD_FOR_OPEN(IOPEN+1) = IB
        IF(MOD(IOPEN-MS2,2).EQ.0) THEN
          IB = IB + NCONF_PER_OPEN(IOPEN+1,ISYM)*NPCMCNF(IOPEN+1)
        END IF
      END DO
*
      RETURN
      END

      SUBROUTINE CSFDIM_FREE(ISYM)
* Free resources allocated by CSFDIM_GAS

#include "implicit.fh"
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "glbbas.fh"
#include "cgas.fh"
#include "spinfo_lucia.fh"
#include "warnings.fh"

      DO IOPEN = MINOP, MAXOP
        ITYP = IOPEN + 1
*
        IALPHA = (IOPEN+MS2)/2
        LZ = IOPEN*IALPHA
        LPTDT = IBION_LUCIA(IOPEN,IALPHA)
        CALL GETMEM('Z_PTDT','FREE','INTE',KZ_PTDT(ITYP),LZ)
        CALL GETMEM('RE_PTD','FREE','INTE',KREO_PTDT(ITYP),LPTDT)
      END DO

c     LDET = NSD_PER_SYM(ISYM)
c     LCONF = 0
c     LLCONF = 0
c     DO IOPEN = 0, MAXOP
c       ITYP = IOPEN + 1
* FIXME: NELEC is undefined
c       ICL = ( NELEC-IOPEN)/2
c       LLCONF = LLCONF + NCONF_PER_OPEN(ITYP,ISYM)*(IOPEN+ICL)
c     END DO
c     LCONF = MAX(LCONF,LLCONF)

      CALL GETMEM('DFTP  ','FREE','INTE',KDFTP,LIDT)
      CALL GETMEM('CFTP  ','FREE','INTE',KCFTP,LICS)
      CALL GETMEM('D_TO_C','FREE','REAL',KDTOC,LDTOC)

      CALL GETMEM('CNFOCC','FREE','INTE',KICONF_OCC(ISYM),LCONF)
      CALL GETMEM('CNFREO','FREE','INTE',KICONF_REO(ISYM),NCONF_TOT)
      CALL GETMEM('SDREOI','FREE','INTE',KSDREO_I(ISYM),LDET)
      END
