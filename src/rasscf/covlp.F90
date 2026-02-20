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

subroutine COVLP(C1IN,C2IN,DIA,PA,SXN,C1,C2,X,OVL)
! RASSCF program version IBM-3090: SX section
!
! Purpose:Calculation of the overlap between two super-CI
! vectors C1IN and C2IN. The result is given in OVL.
! C1,C2, and X are scratch areas.
!
! ********** IBM-3090 Release 89 01 25 **********
!PAM01 Added: replace correct overlap by adding a diagonal
!PAM01 quantity to the overlap of brillouin states.

use Index_Functions, only: iTri
use rasscf_global, only: NROOT, NSXS
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NASH, NISH, NSSH, NSYM
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: C1IN(*), C2IN(*), DIA(*), PA(*), SXN(NSXS)
real(kind=wp), intent(out) :: C1(NSXS), C2(NSXS), X(NSXS), OVL
integer(kind=iwp) :: iAshI, iAshJ, iC1, iC2, iPrLev, ISTBM, ISTC2, iSTIA, ISYM, JSYM, NAE, NAEJ, NAO, NAOJ, NEO, NIA, NIAJ, NIO, &
                     NIOJ, NP, NQ, NT, NTT, NTUT, NTUVX, NU, NUT, NV, NVT, NVXT, NX, NXT
real(kind=wp) :: C1C2, FAC, OVLADD, PRQS, TERM
real(kind=wp), external :: DDot_

IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering COVLP'
!PAM02 Note structure of SX-vectors: First NROOT elements are special.
!PAM02 Elements NROOT+1,..,NROOT+NSXS contain the usual SX elements.
!PAM02 NROOT=1 always right now. Part of the code is prepared for using
!PAM02 several roots, so most of the code must use the general case.
OVL = sum(C1IN(1:NROOT)*C2IN(1:NROOT))

!PAM01 Adding overlap from small shift of SX overlaps:
OVL = OVL+1.0e-6_wp*DDOT_(NSXS,C1IN(NROOT+1),1,C2IN(NROOT+1),1)

! renormalize the C vector (simple element-by-element scaling).

C1(:) = SXN(:)*C1IN(NROOT+1:NROOT+NSXS)
C2(:) = SXN(:)*C2IN(NROOT+1:NROOT+NSXS)

ISTIA = 0
ISTBM = 0
IASHI = 0
do ISYM=1,NSYM
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NIA = NIO+NAO
  NEO = NSSH(ISYM)
  NAE = NAO+NEO
  if ((NIA /= 0) .and. (NAE /= 0)) then

    ! p is secondary (p = q = a)

    if (NEO /= 0) then
      call DGEMM_('N','N',NIA,NEO,NIA,One,DIA(ISTIA+1),NIA,C1(ISTBM+1+NIA*NAO),NIA,Zero,X,NIA)
      OVLADD = DDOT_(NIA*NEO,X,1,C2(ISTBM+1+NIA*NAO),1)
      OVL = OVL+OVLADD
    end if
  end if
  ISTIA = ISTIA+NIA**2
  ISTBM = ISTBM+NIA*NAE
  IASHI = IASHI+NAO
end do

! A very long loop over  symmetry
ISTIA = 0
ISTBM = 0
IASHI = 0
do ISYM=1,NSYM
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NIA = NIO+NAO
  NEO = NSSH(ISYM)
  NAE = NAO+NEO

  ! r is inactive (r = s = i); p and q are active

  if ((NIO /= 0) .and. (NAO /= 0)) then
    IC1 = ISTBM
    do NP=NIO+1,NIA
      IC2 = ISTBM
      do NQ=NIO+1,NIA
        C1C2 = sum(C1(IC1+1:IC1+NIO)*C2(IC2+1:IC2+NIO))
        FAC = -DIA(ISTIA+NIA*(NP-1)+NQ)
        if (NP == NQ) FAC = FAC+Two
        OVLADD = C1C2*FAC
        OVL = OVL+OVLADD
        IC2 = IC2+NIA
      end do
      IC1 = IC1+NIA
    end do
  end if

  ! r,s active and p,q active  (p,r=t,u; q,s=v,x)

  do NT=2,NAO
    NTT = NT+IASHI
    do NU=1,NT-1
      NUT = NU+IASHI
      NTUT = iTri(NTT,NUT)

      IASHJ = 0
      TERM = Zero
      ISTC2 = 0
      do JSYM=1,NSYM
        NAOJ = NASH(JSYM)
        NIOJ = NISH(JSYM)
        NIAJ = NIOJ+NAOJ
        NAEJ = NAOJ+NSSH(JSYM)
        if (NAOJ > 1) then
          if (JSYM == ISYM) then
            !--------
            do NV=2,NAOJ
              NVT = NV+IASHJ
              do NX=1,NV-1
                NXT = NX+IASHJ
                NVXT = iTri(NVT,NXT)
                NTUVX = iTri(NTUT,NVXT)
                PRQS = -Four*PA(NTUVX)
                if (NU == NX) PRQS = PRQS+DIA(ISTIA+NIA*(NT+NIO-1)+NV+NIO)
                if (NT == NV) PRQS = PRQS+DIA(ISTIA+NIA*(NU+NIO-1)+NX+NIO)
                if (NT == NX) PRQS = PRQS-DIA(ISTIA+NIA*(NU+NIO-1)+NV+NIO)
                if (NU == NV) PRQS = PRQS-DIA(ISTIA+NIA*(NT+NIO-1)+NX+NIO)
                TERM = TERM+PRQS*C2(ISTC2+NIAJ*(NV-1)+NIOJ+NX)
              end do
            end do
            !--------
          else
            do NV=2,NAOJ
              NVT = NV+IASHJ
              do NX=1,NV-1
                NXT = NX+IASHJ
                NVXT = iTri(NVT,NXT)
                NTUVX = iTri(NTUT,NVXT)
                PRQS = -Four*PA(NTUVX)
                TERM = TERM+PRQS*C2(ISTC2+NIAJ*(NV-1)+NIOJ+NX)
              end do
            end do
          end if
          !--------
        end if
        ISTC2 = ISTC2+NIAJ*NAEJ
        IASHJ = IASHJ+NAOJ
      end do
      OVLADD = C1(ISTBM+NIA*(NT-1)+NIO+NU)*TERM
      OVL = OVL+OVLADD
    end do
  end do

  ISTIA = ISTIA+NIA**2
  ISTBM = ISTBM+NIA*NAE
  IASHI = IASHI+NAO

  ! End of very long loop over  symmetry
end do

if (IPRLEV >= DEBUG) write(u6,'(1X,A,F15.9)') ' OVERLAP IN COVLP:',OVL

end subroutine COVLP
