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
      SUBROUTINE DENSCT(AREF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
      DIMENSION AREF(*)
      DIMENSION IDC(MXROOT)
*
      IDREST=0
      IDDMO=0
      CALL GETMEM('CI','ALLO','REAL',LCI,NCONF)
      CALL GETMEM('SGM','ALLO','REAL',LSGM,NCONF)
      CALL GETMEM('ASCR2','ALLO','REAL',LASCR2,NVMAX**2)
      CALL GETMEM('BSCR2','ALLO','REAL',LBSCR2,NVMAX**2)
      CALL GETMEM('FSCR2','ALLO','REAL',LFSCR2,NVSQ)
      DO 10 I=1,NRROOT
        IDC(I)=IDREST
        CALL dDAFILE(LUREST,2,Work(LCI),NCONF,IDREST)
        CALL FZERO(Work(LDMO),NBTRI)
*PAM04        IF(ICPF.NE.0) CALL DCORR(HWork(LJREFX),HWork(LAREF),
*PAM04     *                                 HWork(LCSPCK),
*PAM04     *                           HWork(LINTSY),HWork(LINDX),
*PAM04     *                           HWork(LDMO))
        IF(ICPF.NE.0) CALL DCORR(IWork(LJREFX),AREF,
     *                                 IWork(LCSPCK),
     *                           IWork(LINTSY),IWork(LINDX),
     *                           Work(LDMO))
*PAM04        CALL FIJD (HWork(LINTSY),HWork(LINDX),HWork(LCI),
*PAM04     *             HWork(LDMO),HWork(LJREFX),HWork(LAREF))
        CALL FIJD (IWork(LINTSY),IWork(LINDX),Work(LCI),
     *             Work(LDMO),IWork(LJREFX),AREF)
*PAM04        CALL AID *(HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LDMO),
        CALL AID (IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LDMO),
     *            Work(LASCR2),Work(LBSCR2),Work(LFSCR2))
*PAM04        CALL ABD (HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),
        CALL ABD (IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),
*PAM04     *                  HWork(LCI),HWork(LDMO),
     *                  Work(LCI),Work(LDMO),
     *            Work(LASCR2),Work(LBSCR2),Work(LFSCR2),
*PAM04     *            HWork(LJREFX))
     *            IWork(LJREFX))
*PAM04        CALL dDAFILE(LUEIG,1,HWork(LDMO),NBTRI,IDDMO)
        CALL dDAFILE(LUEIG,1,Work(LDMO),NBTRI,IDDMO)
10    CONTINUE
      IF(ITRANS.EQ.0) GOTO 100
      DO 20 I=2,NRROOT
        IDREST=IDC(I)
        CALL dDAFILE(LUREST,2,Work(LCI),NCONF,IDREST)
        DO 21 J=1,I-1
          IDREST=IDC(J)
          CALL dDAFILE(LUREST,2,Work(LSGM),NCONF,IDREST)
*PAM04          CALL FZERO(HWork(LTDMO),NBAST**2)
          CALL FZERO(Work(LTDMO),NBAST**2)
*PAM04          CALL FIJTD (HWork(LINTSY),HWork(LINDX),HWork(LCI),
*PAM04     *                HWork(LSGM),HWork(LTDMO))
          CALL FIJTD (IWork(LINTSY),IWork(LINDX),Work(LCI),
     *                Work(LSGM),Work(LTDMO))
*PAM04          CALL AITD (HWork(LINTSY),HWork(LINDX),HWork(LCI),
          CALL AITD (IWork(LINTSY),IWork(LINDX),Work(LCI),
     *                Work(LSGM),
*PAM04     *               HWork(LTDMO),HWork(LASCR2),HWork(LBSCR2),
     *               Work(LTDMO),Work(LASCR2),Work(LBSCR2),
     *               Work(LFSCR2))
*PAM04          CALL ABTD (HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),
          CALL ABTD (IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),
     *                     Work(LCI),Work(LSGM),
*PAM04     *               HWork(LTDMO),HWork(LASCR2),HWork(LBSCR2),
     *               Work(LTDMO),Work(LASCR2),Work(LBSCR2),
     *               Work(LFSCR2))
*PAM04        CALL dDAFILE(LUEIG,1,HWork(LTDMO),NBAST**2,IDDMO)
          CALL dDAFILE(LUEIG,1,Work(LTDMO),NBAST**2,IDDMO)
21      CONTINUE
20    CONTINUE
100   CONTINUE
      CALL GETMEM('CI','FREE','REAL',LCI,NCONF)
      CALL GETMEM('SGM','FREE','REAL',LSGM,NCONF)
      CALL GETMEM('ASCR2','FREE','REAL',LASCR2,NVMAX**2)
      CALL GETMEM('BSCR2','FREE','REAL',LBSCR2,NVMAX**2)
      CALL GETMEM('FSCR2','FREE','REAL',LFSCR2,NVSQ)
      RETURN
      END
