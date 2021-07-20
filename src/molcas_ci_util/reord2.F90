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
! Copyright (C) 1990,1996, Markus P. Fuelscher                         *
!               1990, Jeppe Olsen                                      *
!***********************************************************************

subroutine Reord2(NORB,NEL,IREFSM,IMODE,ICONF,ISPIN,CIOLD,CINEW,KCNF)
!***********************************************************************
!                                                                      *
!     Rearrange CI-vectors                                             *
!     iMode=0 --> from SGA to split graph GUGA order                   *
!     iMode=1 --> from split graph GUGA to SGA order                   *
!                                                                      *
!     calling arguments:                                               *
!     nOrb    : integer                                                *
!               total number of active orbitals                        *
!     nEl     : integer                                                *
!               total number of active electrons                       *
!     iRefSm  : integer                                                *
!               state symmetry                                         *
!     iMode   : integer                                                *
!               switch selecting reordering mode (see above)           *
!     iConf   : array of integer                                       *
!               string information                                     *
!     iSpin   : array of integer                                       *
!               spin coupling information                              *
!     nSm     : array of integer                                       *
!               symmetry per active orbital                            *
!     CIold   : array of real*8                                        *
!               incoming CI vector                                     *
!     CInew   : array of real*8                                        *
!               outgoing CI vector                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and J. Olsen                                      *
!     University of Lund, Sweden, 1990                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NORB, NEL, IREFSM, IMODE, ICONF(*), ISPIN(*)
real(kind=wp), intent(in) :: CIOLD(*)
real(kind=wp), intent(_OUT_) :: CINEW(*)
integer(kind=iwp), intent(out) :: KCNF(NEL)
#include "rasdim.fh"
#include "spinfo.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
integer(kind=iwp) :: i, IC, ICL, ICNBS, ICNBS0, ICSBAS, ICSFJP, IIBCL, IIBOP, IICSF, IOPEN, IP, IPBAS, IPRLEV, ISG, ITYP, &
                     IWALK(mxAct), JOCC, KOCC, KORB, LPRINT
integer(kind=iwp), external :: IPHASE, ISGNUM

IPRLEV = IPRLOC(3)
! LOOP OVER CONFIGURATIONS TYPES

ICSFJP = 0
ICNBS0 = 0
IPBAS = 0
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1
  ICL = (NEL-IOPEN)/2
  ! BASE ADRESS FOR CONFIGURATION OF THIS TYPE
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
  end if
  ! BASE ADRESS FOR PROTOTYPE SPIN COUPLINGS
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NCSFTP(ITYP-1)*(IOPEN-1)
  end if

  !LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
  !SPIN COUPLINGS

  do IC=1,NCNFTP(ITYP,IREFSM)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    do IICSF=1,NCSFTP(ITYP)
      ICSFJP = ICSFJP+1
      ICSBAS = IPBAS+(IICSF-1)*IOPEN
      ! Obtain configuration in standard RASSCF form
      IIBOP = 1
      IIBCL = 1
      JOCC = ICL+IOPEN
      do KOCC=0,JOCC-1
        KORB = ICONF(ICNBS+KOCC)
        if (KORB < 0) then
          ! Doubly occupied orbital
          KCNF(IIBCL) = abs(KORB)
          IIBCL = IIBCL+1
        else
          ! Singly occupied orbital
          KCNF(ICL+IIBOP) = KORB
          IIBOP = IIBOP+1
        end if
      end do

      ! COMPUTE STEP VECTOR
      call STEPVEC(KCNF(1),KCNF(ICL+1),ICL,IOPEN,ISPIN(ICSBAS),NORB,IWALK)
      ! GET SPLIT GRAPH ORDERING NUMBER
      ISG = ISGNUM(IWORK(LDOWN),IWORK(LUP),IWORK(LDAW),IWORK(LRAW),IWORK(LUSGN),IWORK(LLSGN),IWALK)
      ! GET PHASE PHASE FACTOR
      IP = IPHASE(IWORK(LDRT),IWORK(LUP),IWALK)
      if (IMODE == 0) then
        CINEW(ISG) = CIOLD(ICSFJP)
        if (IP < 0) CINEW(ISG) = -CIOLD(ICSFJP)
      else
        CINEW(ICSFJP) = CIOLD(ISG)
        if (IP < 0) CINEW(ICSFJP) = -CIOLD(ISG)
      end if
    end do
  end do
end do

if (IPRLEV >= DEBUG) then
  LPRINT = min(200,ICSFJP)
  write(u6,*)
  write(u6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
  write(u6,'(10F12.8)') (CIOLD(I),I=1,LPRINT)
  write(u6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
  write(u6,'(10F12.8)') (CINEW(I),I=1,LPRINT)
  write(u6,*)
end if

return

end subroutine Reord2
