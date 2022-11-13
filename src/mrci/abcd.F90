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

subroutine ABCD(INTSYM,indx,ISAB,C,S,ACBDS,ACBDT,BUFIN)

use mrci_global, only: IPASS, IRC, IROW, JJS, KBUFF1, LN, LSYM, Lu_80, NSM, NSYM, NVIR, NVIRP, NVIRT, SQ2, SQ2INV
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), indx(*), ISAB(NVIRT,NVIRT)
real(kind=wp), intent(inout) :: C(*), S(*)
real(kind=wp), intent(_OUT_) :: ACBDS(*), ACBDT(*), BUFIN(*)
integer(kind=iwp) :: IAC, IACMAX, IACMIN, IAD16, IFIN1, IFIN2, ILOOP, IN1, INB, INDA, INPS, INPT, INS, INSB, INSIN, INUMB, ISAC, &
                     IST, IST1, IST2, ISTEP, ISYM, ITAIL, NA, NC, NDMAX, NOV, NSAC, NSACL, NSC, NVT
real(kind=wp) :: TERM
real(kind=wp), external :: DDOT_

!vv hand-made loop unrolling to fix a bug in GCC 3.x
IAD16 = 0
INSIN = KBUFF1
call CSCALE(indx,INTSYM,C,SQ2)
call CSCALE(indx,INTSYM,S,SQ2INV)
NVT = IROW(NVIRT+1)
NOV = (NVT-1)/IPASS+1
IACMAX = 0
!do ISTEP=1,IPASS
ISTEP = 1
if (IPASS >= 1) then

  do
    IACMIN = IACMAX+1
    IACMAX = IACMAX+NOV
    if (IACMAX > NVT) IACMAX = NVT
    if (IACMIN <= IACMAX) then
      !do ISYM=1,NSYM
      ISYM = 1
      if (NSYM >= 1) then
        do
          IST1 = IRC(3)+JJS(ISYM+9)+1
          IFIN1 = IRC(3)+JJS(ISYM+10)
          INPS = IFIN1-IST1+1
          IST2 = IRC(2)+JJS(ISYM)+1
          IFIN2 = IRC(2)+JJS(ISYM+1)
          INPT = IFIN2-IST2+1
          ITAIL = INPS+INPT
          if (ITAIL /= 0) then
            IN1 = -NVIRT
            !do NA=1,NVIRT
            NA = 1
            if (NVIRT >= 1) then
              do
                IN1 = IN1+NVIRT
                do NC=1,NA
                  IAC = IROW(NA)+NC
                  if (IAC < IACMIN) cycle
                  if (IAC > IACMAX) cycle
                  if (NA == 1) cycle
                  NSAC = MUL(NSM(LN+NA),NSM(LN+NC))
                  NSACL = MUL(NSAC,LSYM)
                  if (NSACL /= ISYM) cycle
                  ISAC = ISAB(NA,NC)
                  NSC = NSM(LN+NC)
                  NDMAX = NVIRP(NSC)+NVIR(NSC)
                  if (NDMAX > NA) NDMAX = NA
                  INS = ISAB(NA,NDMAX)
                  ! MOVE INS ITEMS FROM FILE, UNIT 16, VIA BUFFER, INTO ACBDS,
                  ! AND THEN INTO ACBDT:
                  ILOOP = 0
                  do
                    INSB = INS
                    do
                      if (INSIN >= KBUFF1) then
                        ! INSB ITEMS REMAIN TO MOVE.
                        ! INSIN ITEMS HAVE ALREADY BEEN MOVED FROM THE BUFFER.
                        call dDAFILE(Lu_80,2,BUFIN,KBUFF1,IAD16)
                        INSIN = 0
                      end if
                      INB = KBUFF1-INSIN
                      ! INB FRESH ITEMS ARE STILL REMAINING IN BUFFER.
                      INUMB = min(INSB,INB)
                      ! MOVE INUMB ITEMS.
                      IST = INS-INSB+1
                      if (ILOOP == 0) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDS(IST),1)
                      if (ILOOP == 1) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDT(IST),1)
                      INSIN = INSIN+INUMB
                      INSB = INSB-INUMB
                      if (INSB <= 0) exit
                    end do
                    ILOOP = ILOOP+1
                    if (ILOOP /= 1) exit
                  end do
                  ! INS ITEMS HAVE BEEN TRANSFERRED TO ACBDS AND TO ACBDT.
                  if (INPS /= 0) then
                    INDA = IST1
                    if (IFIN1 >= IST1) then
                      !do INDA=IST1,IFIN1
                      do
                        TERM = DDOT_(INS,C(indx(INDA)+1),1,ACBDS,1)
                        S(indx(INDA)+ISAC) = S(indx(INDA)+ISAC)+TERM
                        call DAXPY_(INS,C(indx(INDA)+ISAC),ACBDS,1,S(indx(INDA)+1),1)
                        !end do
                        INDA = INDA+1
                        if (INDA > IFIN1) exit
                      end do
                    end if
                  end if
                  if ((INPT == 0) .or. (NA == NC)) cycle
                  INDA = IST2
                  if (IFIN2 >= IST2) then
                    !do INDA=IST2,IFIN2
                    do
                      TERM = DDOT_(INS,C(indx(INDA)+1),1,ACBDT,1)
                      S(indx(INDA)+ISAC) = S(indx(INDA)+ISAC)+TERM
                      call DAXPY_(INS,C(indx(INDA)+ISAC),ACBDT,1,S(indx(INDA)+1),1)
                      !end do
                      INDA = INDA+1
                      if (INDA > IFIN2) exit
                    end do
                  end if
                end do
                !vv end of unrolling loop
                !NC = NC+1
                !if (NC == NA) exit
                !end do
                NA = NA+1
                if (NA > NVIRT) exit
              end do
            end if
          end if
          !end do
          ISYM = ISYM+1
          if (ISYM > NSYM) exit
        end do
      end if
    end if

    !end do
    ISTEP = ISTEP+1
    if (ISTEP > IPASS) exit
  end do
end if
call CSCALE(indx,INTSYM,C,SQ2INV)
call CSCALE(indx,INTSYM,S,SQ2)

return

end subroutine ABCD
