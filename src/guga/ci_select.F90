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

subroutine CI_SELECT(L0,L1,L2,L3,KBUF,NTPB,NBINS,LW1)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use guga_global, only: IADD10, IADD11, ICH, IFIRST, ILIM, IOUT, LN, Lu_10, MXVERT, NBUF, NMAT, NSM
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6, RtoI

implicit none
integer(kind=iwp), intent(in) :: L0(4*MXVERT), L1(4*MXVERT), L2(4*MXVERT), L3(4*MXVERT), KBUF, NTPB, NBINS, LW1
integer(kind=iwp) :: ID, IDIAG, IFIN, II, IID, IND1, IND2, IREC, IST, IT1, IT2, ITIM, ITT, JJ, JJD, JTYP, KBUF2, M1, M2, M2MIN, &
                     M3, M4, NA, NB, NC, ND, NSA, NSABC, NSAV1, NSAV2, NSAVE, NSB, NSC, NSCD, NSD
integer(kind=iwp), allocatable :: ICAD(:), IBUFL(:)
integer(kind=iwp), allocatable, target :: INDOUT(:)
real(kind=wp), allocatable :: BUFL(:)
real(kind=wp), pointer :: dINDOUT(:)

call JTIME(IST)
KBUF2 = (RtoI+1)*KBUF+2
call mma_allocate(ICAD,NBINS,label='ICAD')
call mma_allocate(IBUFL,max(NBINS,KBUF+2),label='IBUFL')
call mma_allocate(BUFL,KBUF,label='BUFL')
call mma_allocate(INDOUT,LW1,label='INDOUT')

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
call c_f_pointer(c_loc(INDOUT),dINDOUT,[1])
do M3=1,LN
  do M4=1,M3
    do M1=M3,LN
      M2MIN = 1
      if (M1 == M3) M2MIN = M4
      do M2=M2MIN,M1
        if ((M1 /= M3) .or. (M2 /= M4)) cycle
        if (M1 == M2) cycle
        NA = ICH(M1)
        NB = ICH(M2)
        if (NA < NB) then
          NSAVE = NA
          NA = NB
          NB = NSAVE
        end if
        NC = ICH(M3)
        ND = ICH(M4)
        if (NC < ND) then
          NSAVE = NC
          NC = ND
          ND = NSAVE
        end if
        if (NA <= NC) then
          if (NA /= NC) then
            NSAV1 = NA
            NSAV2 = NB
            NA = NC
            NB = ND
            NC = NSAV1
            ND = NSAV2
          else if (NB <= ND) then
            NSAVE = NB
            NB = ND
            ND = NSAVE
          end if
        end if
        call INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
      end do
    end do
  end do
end do
call AIAI(dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB,NBINS)
call EMPTY(BUFL,IBUFL,ICAD,dINDOUT,KBUF,NTPB)
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(u6,601) NMAT
call JTIME(IFIN)
ITIM = IFIN-IST
write(u6,701) ITIM
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
write(u6,702) ITIM
IAD10(5) = IADD10
IST = IFIN
IOUT = 0
NMAT = 0
IDIAG = 0
do M3=1,LN
  NSC = NSM(ICH(M3))
  do M4=1,M3
    NSD = NSM(ICH(M4))
    NSCD = Mul(NSC,NSD)
    do M1=M3,LN
      NSA = NSM(ICH(M1))
      NSABC = Mul(NSA,NSCD)
      M2MIN = 1
      if (M1 == M3) M2MIN = M4
      do M2=M2MIN,M1
        NA = ICH(M1)
        NB = ICH(M2)
        NSB = NSM(NB)
        if (NSB /= NSABC) cycle
        if (NA < NB) then
          NSAVE = NA
          NA = NB
          NB = NSAVE
        end if
        NC = ICH(M3)
        ND = ICH(M4)
        if (NC < ND) then
          NSAVE = NC
          NC = ND
          ND = NSAVE
        end if
        if (NA <= NC) then
          if (NA /= NC) then
            NSAV1 = NA
            NSAV2 = NB
            NA = NC
            NB = ND
            NC = NSAV1
            ND = NSAV2
          else if (NB <= ND) then
            NSAVE = NB
            NB = ND
            ND = NSAVE
          end if
        end if
        IOUT = IOUT+1
        ICOP1(IOUT) = 0
        if (IOUT >= NBUF) then
          ICOP1(nCOP+1) = NBUF
          call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
          call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
          NMAT = NMAT+NBUF
          IOUT = 0
        end if
        IOUT = IOUT+1
        !IND1 = NA+2**8*NB
        !IND2 = IND1+2**16*NC
        !ICOP1(IOUT) = IND2+2**24*ND
        IND1 = ior(NA,ishft(NB,8))
        IND2 = ior(IND1,ishft(NC,16))
        ICOP1(IOUT) = ior(IND2,ishft(ND,24))
        if (IOUT >= NBUF) then
          ICOP1(nCOP+1) = NBUF
          call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
          call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
          NMAT = NMAT+NBUF
          IOUT = 0
        end if
        if (NA /= NC) then
          if (NB /= ND) then
            if (NB /= NC) then
              if (NA == NB) then
                if (NC /= ND) then
                  do ITT=1,ILIM
                    IT1 = (ITT-1)*MXVERT
                    IT2 = IT1
                    call INT8(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
                  end do
                end if
              else if (NC /= ND) then
                if (NB > ND) then
                  if (NB > NC) then
                    do ITT=1,ILIM
                      IT1 = (ITT-1)*MXVERT
                      IT2 = IT1
                      call INT1(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
                    end do
                  else
                    do ITT=1,ILIM
                      IT1 = (ITT-1)*MXVERT
                      IT2 = IT1
                      call INT2(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
                    end do
                  end if
                else
                  do ITT=1,ILIM
                    IT1 = (ITT-1)*MXVERT
                    IT2 = IT1
                    call INT3(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
                  end do
                end if
              else
                do ITT=1,ILIM
                  IT1 = (ITT-1)*MXVERT
                  IT2 = IT1
                  call INT8(NB,NA,NC,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
                end do
              end if
            else
              do ITT=1,ILIM
                IT1 = (ITT-1)*MXVERT
                IT2 = IT1
                call INT4(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
              end do
            end if
          else if (NC == ND) then
            do ITT=1,ILIM
              IT1 = (ITT-1)*MXVERT
              IT2 = IT1
              call INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
            end do
          else
            call INT5(ND,NC,NA)
          end if
        else if (NB == ND) then
          if (NA /= NB) call INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
        else if (NB == NC) then
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
        else
          do ITT=1,ILIM
            IT1 = (ITT-1)*MXVERT
            IT2 = IT1
            call INT61(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
            call INT62(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,L0,L1,L2,L3)
          end do
        end if
      end do
    end do
  end do
end do
nullify(dINDOUT)
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(u6,600) NMAT
call JTIME(IFIN)
ITIM = IFIN-IST
write(u6,703) ITIM
IAD10(6) = IADD10
IST = IFIN
NMAT = 0
call AIBJ(L0,L1,L2,L3,INDOUT)
call JTIME(IFIN)
ITIM = IFIN-IST
write(u6,704) ITIM
IAD10(7) = IADD10
IST = IFIN
call AIJK(INDOUT,L0,L1,L2,L3)
call JTIME(IFIN)
ITIM = IFIN-IST
write(u6,705) ITIM
IAD10(8) = IADD10
IST = IFIN
call ONEEL_GUGA()
IAD10(9) = IADD10
JTYP = 0
call AI(JTYP,INDOUT,L0,L1,L2,L3)
call JTIME(IFIN)
ITIM = IFIN-IST
write(u6,706) ITIM

call mma_deallocate(ICAD)
call mma_deallocate(IBUFL)
call mma_deallocate(BUFL)
call mma_deallocate(INDOUT)

return

600 format(/6X,'COEFFICIENTS FOR IJKL',I9)
601 format(/6X,'COEFFICIENTS FOR DIAG',I9)
701 format(6X,'TIME FOR DIAG',I17)
702 format(6X,'TIME FOR ABCI',I17)
703 format(6X,'TIME FOR IJKL',I17)
704 format(6X,'TIME FOR AIBJ',I17)
705 format(6X,'TIME FOR AIJK',I17)
706 format(6X,'TIME FOR ONEEL',I16)

end subroutine CI_SELECT
