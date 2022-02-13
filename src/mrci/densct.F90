!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine DENSCT(AREF)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: AREF(*)
#include "mrci.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: I, IDC(MXROOT), IDDMO, IDREST, J !IFG

IDREST = 0
IDDMO = 0
call GETMEM('CI','ALLO','REAL',LCI,NCONF)
call GETMEM('SGM','ALLO','REAL',LSGM,NCONF)
call GETMEM('ASCR2','ALLO','REAL',LASCR2,NVMAX**2)
call GETMEM('BSCR2','ALLO','REAL',LBSCR2,NVMAX**2)
call GETMEM('FSCR2','ALLO','REAL',LFSCR2,NVSQ)
do I=1,NRROOT
  IDC(I) = IDREST
  call dDAFILE(LUREST,2,Work(LCI),NCONF,IDREST)
  call FZERO(Work(LDMO),NBTRI)
  !PAM04 if (ICPF /= 0) call DCORR(HWork(LJREFX),HWork(LAREF),HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),HWork(LDMO))
  if (ICPF /= 0) call DCORR(IWork(LJREFX),AREF,IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),Work(LDMO))
  !PAM04 call FIJD (HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LDMO),HWork(LJREFX),HWork(LAREF))
  call FIJD(IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LDMO),IWork(LJREFX),AREF)
  !PAM04 call AID (HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LDMO),
  call AID(IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LDMO),Work(LASCR2),Work(LBSCR2),Work(LFSCR2))
  !PAM04 call ABD (HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LDMO),HWork(LJREFX))
  call ABD(IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LDMO),Work(LASCR2),Work(LBSCR2),Work(LFSCR2),IWork(LJREFX))
  !PAM04 call dDAFILE(LUEIG,1,HWork(LDMO),NBTRI,IDDMO)
  call dDAFILE(LUEIG,1,Work(LDMO),NBTRI,IDDMO)
end do
if (ITRANS /= 0) then
  do I=2,NRROOT
    IDREST = IDC(I)
    call dDAFILE(LUREST,2,Work(LCI),NCONF,IDREST)
    do J=1,I-1
      IDREST = IDC(J)
      call dDAFILE(LUREST,2,Work(LSGM),NCONF,IDREST)
      !PAM04 call FZERO(HWork(LTDMO),NBAST**2)
      call FZERO(Work(LTDMO),NBAST**2)
      !PAM04 call FIJTD (HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LSGM),HWork(LTDMO))
      call FIJTD(IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LSGM),Work(LTDMO))
      !PAM04 call AITD (HWork(LINTSY),HWork(LINDX),HWork(LCI),HWork(LTDMO),HWork(LASCR2),HWork(LBSCR2),
      call AITD(IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LSGM),Work(LTDMO),Work(LASCR2),Work(LBSCR2),Work(LFSCR2))
      !PAM04 call ABTD (HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),HWork(LTDMO),HWork(LASCR2),HWork(LBSCR2),
      call ABTD(IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),Work(LCI),Work(LSGM),Work(LTDMO),Work(LASCR2),Work(LBSCR2),Work(LFSCR2))
      !PAM04 call dDAFILE(LUEIG,1,HWork(LTDMO),NBAST**2,IDDMO)
      call dDAFILE(LUEIG,1,Work(LTDMO),NBAST**2,IDDMO)
    end do
  end do
end if
call GETMEM('CI','FREE','REAL',LCI,NCONF)
call GETMEM('SGM','FREE','REAL',LSGM,NCONF)
call GETMEM('ASCR2','FREE','REAL',LASCR2,NVMAX**2)
call GETMEM('BSCR2','FREE','REAL',LBSCR2,NVMAX**2)
call GETMEM('FSCR2','FREE','REAL',LFSCR2,NVSQ)

return

end subroutine DENSCT
