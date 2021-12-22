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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine CI_SELECT(INDOUT,ICAD,IBUFL,L0,L1,L2,L3,KBUF,NTPB,NBINS)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: INDOUT(*), ICAD(*), IBUFL(*)
integer(kind=iwp), intent(in) :: L0(*), L1(*), L2(*), L3(*), KBUF, NTPB, NBINS
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
integer(kind=iwp) :: ID, IDIAG, IFIN, II, IID, IND1, IND2, IST, IT1, IT2, ITIM, JJ, JJD, JTYP, KBUF2, M2MIN, NA, NB, NC, ND, NSA, &
                     NSABC, NSAV1, NSAV2, NSAVE, NSB, NSC, NSCD, NSD

call JTIME(IST)
call CI_SELECT_INTERNAL(INDOUT,IBUFL)

return

! This is to allow type punning without an explicit interface
contains

subroutine CI_SELECT_INTERNAL(INDOUT,IBUFL)

  integer(kind=iwp), target, intent(_OUT_) :: INDOUT(*), IBUFL(*)
  integer(kind=iwp) :: IREC, ITT, M1, M2, M3, M4
  real(kind=wp), pointer :: dINDOUT(:), dIBUFL(:)

  call c_f_pointer(c_loc(INDOUT),dINDOUT,[1])
  call c_f_pointer(c_loc(IBUFL),dIBUFL,[1])

  KBUF2 = (RTOI+1)*KBUF+2
  ID = 0
  do IREC=1,NBINS
    IBUFL(IREC) = 0
    ICAD(IREC) = ID
    INDOUT(ID+KBUF2) = -1
    ID = ID+KBUF2
  end do
  II = 0
  IID = 0
  JJ = 0
  JJD = 0
  JTYP = 0
  ! DIAGONAL ELEMENTS
  IADD11 = 0
  IAD10(3) = IADD10
  IOUT = 0
  NMAT = 0
  IDIAG = 1
  do M3=1,LN
    do M4=1,M3
      do M1=M3,LN
        M2MIN = 1
        if (M1 == M3) M2MIN = M4
        do M2=M2MIN,M1
          if ((M1 /= M3) .or. (M2 /= M4)) GO TO 310
          if (M1 == M2) GO TO 310
          NA = ICH(M1)
          NB = ICH(M2)
          if (NA >= NB) GO TO 131
          NSAVE = NA
          NA = NB
          NB = NSAVE
131       NC = ICH(M3)
          ND = ICH(M4)
          if (NC >= ND) GO TO 130
          NSAVE = NC
          NC = ND
          ND = NSAVE
130       if (NA > NC) GO TO 132
          if (NA == NC) GO TO 133
          NSAV1 = NA
          NSAV2 = NB
          NA = NC
          NB = ND
          NC = NSAV1
          ND = NSAV2
          GO TO 132
133       if (NB > ND) GO TO 132
          NSAVE = NB
          NB = ND
          ND = NSAVE
132       call INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
310       continue
        end do
      end do
    end do
  end do
  call AIAI(dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB,NBINS)
  call EMPTY(dIBUFL,IBUFL(RtoI*kBuf+1),ICAD,dINDOUT,KBUF,NTPB)
  ICOP1(nCOP+1) = IOUT
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+IOUT
  ICOP1(nCOP+1) = -1
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  write(IW,601) NMAT
601 format(/6X,'COEFFICIENTS FOR DIAG',I9)
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,701) ITIM
701 format(6X,'TIME FOR DIAG',I17)
  IAD10(4) = IADD10
  IST = IFIN
  JTYP = 1
  !FUE start modification
  !if (IFIRST == 0) call AI(JTYP,INDOUT,L0,L1,L2,L3)
  if (IFIRST == 0) then
    call AI(JTYP,INDOUT,L0,L1,L2,L3)
  else
    ICOP1(nCOP+1) = 0
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    ICOP1(nCOP+1) = -1
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  end if
  !FUE end modification
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,702) ITIM
702 format(6X,'TIME FOR ABCI',I17)
  IAD10(5) = IADD10
  IST = IFIN
  IOUT = 0
  NMAT = 0
  IDIAG = 0
  do M3=1,LN
    NSC = NSM(ICH(M3))
    do M4=1,M3
      NSD = NSM(ICH(M4))
      NSCD = MUL(NSC,NSD)
      do M1=M3,LN
        NSA = NSM(ICH(M1))
        NSABC = MUL(NSA,NSCD)
        M2MIN = 1
        if (M1 == M3) M2MIN = M4
        do M2=M2MIN,M1
          NA = ICH(M1)
          NB = ICH(M2)
          NSB = NSM(NB)
          if (NSB /= NSABC) GO TO 510
          if (NA >= NB) GO TO 231
          NSAVE = NA
          NA = NB
          NB = NSAVE
231       NC = ICH(M3)
          ND = ICH(M4)
          if (NC >= ND) GO TO 230
          NSAVE = NC
          NC = ND
          ND = NSAVE
230       if (NA > NC) GO TO 232
          if (NA == NC) GO TO 233
          NSAV1 = NA
          NSAV2 = NB
          NA = NC
          NB = ND
          NC = NSAV1
          ND = NSAV2
          GO TO 232
233       if (NB > ND) GO TO 232
          NSAVE = NB
          NB = ND
          ND = NSAVE
232       IOUT = IOUT+1
          ICOP1(IOUT) = 0
          if (IOUT < NBUF) GO TO 460
          ICOP1(nCOP+1) = NBUF
          call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
          call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
          NMAT = NMAT+NBUF
          IOUT = 0
460       IOUT = IOUT+1
          !IND1 = NA+2**8*NB
          !IND2 = IND1+2**16*NC
          !ICOP1(IOUT) = IND2+2**24*ND
          IND1 = ior(NA,ishft(NB,8))
          IND2 = ior(IND1,ishft(NC,16))
          ICOP1(IOUT) = ior(IND2,ishft(ND,24))
          if (IOUT < NBUF) GO TO 511
          ICOP1(nCOP+1) = NBUF
          call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
          call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
          NMAT = NMAT+NBUF
          IOUT = 0
511       if (NA /= NC) GO TO 520
          if (NB == ND) GO TO 525
          if (NB == NC) GO TO 524
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT61(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
            call INT62(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
524       do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
525       if (NA == NB) GO TO 510
          call INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
          GO TO 510
520       if (NB /= ND) GO TO 530
          if (NC == ND) GO TO 529
          call INT5(ND,NC,NA)
          GO TO 510
529       do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
530       if (NB /= NC) GO TO 535
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT4(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
535       if (NA == NB) GO TO 546
          if (NC /= ND) GO TO 547
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT8(NB,NA,NC,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
546       if (NC == ND) GO TO 510
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT8(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
547       if (NB > ND) GO TO 540
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT3(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
540       if (NB > NC) GO TO 545
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT2(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
          GO TO 510
545       do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT1(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
510       continue
        end do
      end do
    end do
  end do
  ICOP1(nCOP+1) = IOUT
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+IOUT
  ICOP1(nCOP+1) = -1
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  write(IW,600) NMAT
600 format(/6X,'COEFFICIENTS FOR IJKL',I9)
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,703) ITIM
703 format(6X,'TIME FOR IJKL',I17)
  IAD10(6) = IADD10
  IST = IFIN
  NMAT = 0
  call AIBJ(L0,L1,L2,L3,INDOUT)
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,704) ITIM
704 format(6X,'TIME FOR AIBJ',I17)
  IAD10(7) = IADD10
  IST = IFIN
  call AIJK(INDOUT,L0,L1,L2,L3)
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,705) ITIM
705 format(6X,'TIME FOR AIJK',I17)
  IAD10(8) = IADD10
  IST = IFIN
  call ONEEL_GUGA()
  IAD10(9) = IADD10
  JTYP = 0
  call AI(JTYP,INDOUT,L0,L1,L2,L3)
  call JTIME(IFIN)
  ITIM = IFIN-IST
  write(IW,706) ITIM
706 format(6X,'TIME FOR ONEEL',I16)
  nullify(dINDOUT,dIBUFL)

  return

end subroutine CI_SELECT_INTERNAL

end subroutine CI_SELECT
