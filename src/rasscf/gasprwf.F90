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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine gasprwf(NORB,NEL,IREFSM,ICONF,ISPIN,CICOEF,kcnf)
!***********************************************************************
!                                                                      *
!     PURPOSE: PRINT THE WAVEFUNCTION FOR GAS                          *
!                                                                      *
!***********************************************************************
!                                                                      *
!     calling arguments:                                               *
!     nOrb    : integer                                                *
!               total number of active orbitals                        *
!     nEl     : integer                                                *
!               total number of active electrons                       *
!     iRefSm  : integer                                                *
!               state symmetry                                         *
!     iConf   : array of integer                                       *
!               string information                                     *
!     iSpin   : array of integer                                       *
!               spin coupling information                              *
!     nSm     : array of integer                                       *
!               symmetry per active orbital                            *
!     CiCoef   : array of real*8                                       *
!               incoming CI vector                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                      *
!***********************************************************************

use rasscf_global, only: PrwThr, nSm
use output_ras, only: LF
use spinfo, only: NTYP, MINOP, NCNFTP, NCSFTP
use Molcas, only: MxAct

implicit none
integer nOrb, nEl
integer ICONF(*), ISPIN(*)
real*8 CICOEF(*)
integer KCNF(NEL)
integer IWALK(mxAct)
character(len=120) Line
integer iRefSM, IC, ICL, ICNBS, ICNBS0, iCSBAS, ICSFJP, IIBCL, IIBOP, IICSF, iOff, iOpen, iOrb, ipBas, iSym, iTyp, jOCC, kOCC, kOrb
real*8 COEF

! print headline

Line(1:16) = '      Conf/sym  '
iOff = 16
iSym = nSm(1)
do iorb=1,norb
  if (nsm(iorb) /= isym) ioff = ioff+1
  write(line(ioff+iorb:),'(I1)') nsm(iorb)
  if (nsm(iorb) /= isym) isym = nsm(iorb)
end do
iOff = iOff+norb+3
Line(iOff:iOff+15) = '   Coeff Weight'
write(LF,'(A)') Line(1:iOff+15)
Line = ' '

! Loop over configuration types

ICSFJP = 0
ICNBS0 = 0
IPBAS = 0
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
      ! Obtain configuration in standard RASSCF form
      IIBOP = 1
      IIBCL = 1
      JOCC = ICL+IOPEN
      do KOCC=0,JOCC-1
        KORB = ICONF(ICNBS+KOCC)
        if (KORB < 0) then
          ! Doubly occupied orbitals
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
      ! SKIP IT OR PRINT IT?
      COEF = CICOEF(ICSFJP)
      if (abs(COEF) < PRWTHR) cycle
      ! PRINT IT
      write(Line(1:),'(I8)') icsfjp
      iOff = 10
      iSym = nSm(1)
      do iorb=1,norb
        if (nSm(iorb) /= iSym) iOff = iOff+1
        if (iwalk(iorb) == 3) then
          write(Line(iOff+iorb:),'(A1)') '2'
        else if (iwalk(iorb) == 2) then
          write(Line(iOff+iorb:),'(A1)') 'd'
        else if (iwalk(iorb) == 1) then
          write(Line(iOff+iorb:),'(A1)') 'u'
        else if (iwalk(iorb) == 0) then
          write(Line(iOff+iorb:),'(A1)') '0'
        end if
        if (nSm(iorb) /= iSym) iSym = nSm(iorb)
      end do
      iOff = iOff+norb+3
      write(Line(iOff:),'(2F8.5)') COEF,COEF**2
      write(LF,'(6X,A)') Line(1:iOff+15)
      Line = ' '

    end do
  end do
end do

end subroutine gasprwf
