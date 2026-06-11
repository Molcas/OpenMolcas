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

subroutine SG_ReOrd(SGS,EXS,IREFSM,IMODE,ICONF,ISPIN,nConf,CIOLD,CINEW)
!***********************************************************************
!                                                                      *
!     Rearrange CI-vectors                                             *
!     iMode=0 --> from SGA to split graph GUGA order                   *
!     iMode=1 --> from split graph GUGA to SGA order                   *
!                                                                      *
!     calling arguments:                                               *
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
!     CIold   : array of real                                          *
!               incoming CI vector                                     *
!     CInew   : array of real                                          *
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

use sguga, only: EXStruct, SGStruct
use output_ras, only: IPRLOC
use spinfo, only: MINOP, NCNFTP, NCSFTP, NTYP
use PrintLevel, only: DEBUG
use Molcas, only: MxAct
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
type(SGStruct), intent(in) :: SGS
type(EXStruct), intent(in) :: EXS
integer(kind=iwp), intent(in) :: IREFSM, IMODE, ICONF(*), ISPIN(*), nConf
real(kind=wp), intent(in) :: CIOLD(nConf)
real(kind=wp), intent(_OUT_) :: CINEW(nConf)

integer(kind=iwp) :: i, IC, ICL, ICNBS, ICNBS0, ICSBAS, ICSFJP, IIBCL, IIBOP, IICSF, IOPEN, IP, IPBAS, IPRLEV, ISG, ITYP, &
                     IWALK(mxAct), JOCC, KOCC, KORB, LPRINT, nOrb, nEl
integer(kind=iwp), external :: SG_PHASE, SG_NUM
integer(kind=iwp) :: KCNF(MxAct)

IPRLEV = IPRLOC(3)

nOrb=SGS%nLev
nEl=SGS%nActEl

ICSFJP = 0
ICNBS0 = 0 ! dummy initialize
IPBAS = 0 ! dummy initialize
! LOOP OVER CONFIGURATIONS TYPES
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1
  ICL = (NEL-IOPEN)/2
  ! BASE ADDRESS FOR CONFIGURATION OF THIS TYPE
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
  end if
  ! BASE ADDRESS FOR PROTOTYPE SPIN COUPLINGS
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NCSFTP(ITYP-1)*(IOPEN-1)
  end if

  ! LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
  ! SPIN COUPLINGS

  do IC=1,NCNFTP(ITYP,IREFSM)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    do IICSF=1,NCSFTP(ITYP)
      ICSFJP = ICSFJP+1
      ICSBAS = IPBAS+(IICSF-1)*IOPEN
      KCNF(:)=0
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
      call STEPVEC(KCNF(1:ICL),KCNF(ICL+1),ICL,IOPEN,ISPIN(ICSBAS),NORB,IWALK)

      ! GET SPLIT GRAPH ORDERING NUMBER
      ISG = SG_NUM(SGS,EXS,IWALK)
      ! GET PHASE PHASE FACTOR
      IP = SG_PHASE(SGS,IWALK)
      if (IMODE == 0) then
        if (IP < 0) then
          CINEW(ISG) = -CIOLD(ICSFJP)
        else
          CINEW(ISG) = CIOLD(ICSFJP)
        end if
      else
        if (IP < 0) then
          CINEW(ICSFJP) = -CIOLD(ISG)
        else
          CINEW(ICSFJP) = CIOLD(ISG)
        end if
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

end subroutine SG_Reord
