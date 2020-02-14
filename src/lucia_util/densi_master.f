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
      SUBROUTINE densi_master(rvec)
*
* Controls the calculation of the densities, when Lucia is called
* from Molcas Rasscf.
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "clunit.fh"
#include "glbbas.fh"
#include "orbinp.fh"
#include "lucinp.fh"
#include "spinfo_lucia.fh"
#include "cstate.fh"
#include "rasscf_lucia.fh"
#include "io_util.fh"
      integer rvec
      logical iPack,tdm
      dimension dummy(1)
*
* Put CI-vector from RASSCF on luc
*
*     if rvec>0, it should be a pointer to a second CI vector
*     and a one-particle transition density matrix will be computed
      tdm = rvec.gt.0
      CALL GETMEM('LSCR1 ','ALLO','REAL',LSCR1,NSD_PER_SYM(IREFSM))
      CALL GETMEM('LSCR2 ','ALLO','REAL',LSCR2,NSD_PER_SYM(IREFSM))
      CALL COPVEC(WORK(C_POINTER),WORK(LSCR1),NCSF_PER_SYM(IREFSM))
      ITMP_POINTER = C_POINTER
      Call GetMem('lrec','allo','inte',ivlrec,MXNTTS)
      IF (tdm) THEN
         CALL GETMEM('LSCR3 ','ALLO','REAL',LSCR3,NSD_PER_SYM(IREFSM))
         CALL GETMEM('LSCR4 ','ALLO','REAL',LSCR4,NSD_PER_SYM(IREFSM))
         CALL COPVEC(WORK(rvec),WORK(LSCR3),NCSF_PER_SYM(IREFSM))
         CALL CSDTVC(WORK(LSCR3),WORK(LSCR4),1,WORK(KDTOC_POINTER),
     &               iWORK(KSDREO_POINTER), IREFSM, 1)
         C_POINTER = LSCR3
         CALL CPCIVC(LUHC, MXNTTS, IREFSM, 1,iwork(ivlrec))
      END IF
      C_POINTER = LSCR1
      CALL CSDTVC(WORK(LSCR1),WORK(LSCR2),1,WORK(KDTOC_POINTER),
     &     iWORK(KSDREO_POINTER), IREFSM, 1)
      CALL CPCIVC(LUC, MXNTTS, IREFSM, 1,iwork(ivlrec))
      Call GetMem('lrec','free','inte',ivlrec,MXNTTS)

*
* Determine length of arrays VEC1 and VEC2
*
c      IF(ISIMSYM.EQ.0) THEN
         LBLOCK = MXSOOB
c      ELSE
c         LBLOCK = MXSOOB_AS
c      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
*PAM06      LBLOCK = MAX(XISPSM(IREFSM,1),DBLE(MXSOOB))
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = 2*INT(XISPSM(IREFSM,1))
*
* Allocate arrays
*
      IDUM=0
*     CALL MEMMAN(IDUM, IDUM, 'MARK', IDUM, 'DENS_M')
      CALL GETMEM('VEC1  ','ALLO','REAL',KVEC1,LBLOCK)
      CALL GETMEM('KC2   ','ALLO','REAL',KVEC3,kvec3_length)
*
* Copy Sigma-vector from disc to core
*
      CALL GETMEM('VEC2  ','ALLO','REAL',KVEC2,LBLOCK)
       IF (iSigma_on_disk .ne. 0) THEN
          Call GetMem('lvec','Allo','inte',ivlrec,MXNTTS)
          CALL cpsivc(lusc34, mxntts, work(kvec2),iWork(ivlrec))
          Call GetMem('lvec','Free','inte',ivlrec,MXNTTS)
       ELSE
         Do i = 1, Lblock
            work(kvec2+i-1) = 0.0d0
         Enddo
       ENDIF
*
* Information needed on file handling
*
      LBLK = - 1
*
* Copy vector on file LUC to LUSC1 and LUHC
*
      IDISK(LUC)=0
      IDISK(LUSC1)=0
      CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,LBLK)
      IF (.not.tdm) CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,LBLK)
*
* Calculate one- and two-body densities
*
      IPACK = .TRUE.
      DUMMY = 0.0D0
      IF (tdm) THEN
         CALL densi2_lucia(1,work(lw6),dummy,dummy,dummy,
     &   work(kvec1),work(kvec2),lusc1,luhc,exps2,1,work(lw7),IPACK)
      ELSE
         CALL densi2_lucia(2,work(krho1),dummy,Work(lw8),Work(lw9),
     &   work(kvec1),work(kvec2),lusc1,luhc,exps2,1,work(ksrho1),IPACK)
      END IF

*
* Explanation of calling parameters
*
C      2      : DONE!!! - Calculate both one and two body densities.
C      krho1  : DONE!!! - Output - include in glbbas.fh.
C      krho2  : DONE!!! - Output - include in glbbas.fh.
C      kvec1  : DONE!!! - CI-vector
C      kvec2  : DONE!!! - Sigma-vector
C      lusc1  : DONE!!! - file pointer
C      luhc   : DONE!!! - file pointer
C      exps2  : DONE!!! - Output - expectation value of S**2.
C      1      : DONE!!! - Calculate spin density
C      ksrho1 : DONE!!! - Comming with glbbas.fh.
*
      IF (.not.tdm) THEN
*        Save densities in trigonal format for use in Molcas
*
         CALL TriPak(work(krho1), work(lw6), 1, ntoob, ntoob)
         CALL TriPak(work(ksrho1), work(lw7), 1, ntoob, ntoob)
      END IF
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
*
      CALL CSDTVC(work(lscr1),work(lscr2),2,work(kdtoc_pointer),
     &     iwork(KSDREO_POINTER), iRefSm, 1)
      C_POINTER = iTmp_pointer
*
*     CALL MEMMAN(IDUM, IDUM, 'FLUSM', IDUM, 'DENS_M')
      CALL GETMEM('LSCR1 ','FREE','REAL',LSCR1,NSD_PER_SYM(IREFSM))
      CALL GETMEM('LSCR2 ','FREE','REAL',LSCR2,NSD_PER_SYM(IREFSM))
      CALL GETMEM('VEC1  ','FREE','REAL',KVEC1,LBLOCK)
      CALL GETMEM('KC2   ','FREE','REAL',KVEC3,kvec3_length)
      CALL GETMEM('VEC2  ','FREE','REAL',KVEC2,LBLOCK)
      IF (tdm) THEN
         CALL GETMEM('LSCR3 ','FREE','REAL',LSCR3,NSD_PER_SYM(IREFSM))
         CALL GETMEM('LSCR4 ','FREE','REAL',LSCR4,NSD_PER_SYM(IREFSM))
      END IF
*
      RETURN
      END
*
