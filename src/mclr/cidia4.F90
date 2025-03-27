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
! Copyright (C) 1991,1994, Jeppe Olsen                                 *
!***********************************************************************

subroutine CIDIA4(NAEL,IASTR,NBEL,IBSTR,NORB,DIAG,NSMST,H,ISMOST,IBLTP,XA,XB,SCR,RJ,RK,NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,ISSOA, &
                  ISSOB,LUDIA,ECORE,PLSIGN,PSSIGN,NTOOB,ICISTR)
! Calculate determinant diagonal
! Turbo-ras version
!
! ========================
! General symmetry version
! ========================
!
! Jeppe Olsen, Winter of 1991
! K => J - K moved outside, April 1994

use Constants, only: Zero, One, Half
use Definitions, only: wp

implicit real*8(A-H,O-Z)
! General input
dimension NSSOA(NOCTPA,*), NSSOB(NOCTPB,*)
dimension ISSOA(NOCTPA,*), ISSOB(NOCTPB,*)
dimension IASTR(NAEL,*), IBSTR(NBEL,*)
dimension H(NORB)
! Specific input
dimension IOCOC(NOCTPA,NOCTPB)
dimension ISMOST(*), IBLTP(*)
! Scratch
dimension RJ(NTOOB,NTOOB), RK(NTOOB,NTOOB)
dimension XA(NORB), XB(NORB), SCR(2*NORB)
! Output
dimension DIAG(*)
dimension IDUM(1)

if (PSSIGN == -One) then
  XADD = 1.0e6_wp
else
  XADD = Zero
end if

!*3 Diagonal elements according to Handys formulae
!   (corrected for error)
!
!   DIAG(IDET) = HII*(NIA+NIB)
!              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
!              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
!              +         J(I,J) * NIA*NJB

! K goes to J - K

IDET = 0
ITDET = 0
if (LUDIA /= 0) rewind LUDIA
do IASM=1,NSMST
  IBSM = ISMOST(IASM)
  if ((IBSM == 0) .or. (IBLTP(IASM) == 0)) cycle
  if (IBLTP(IASM) == 2) then
    IREST1 = 1
  else
    IREST1 = 0
  end if

  do IATP=1,NOCTPA
    if (IREST1 == 1) then
      MXBTP = IATP
    else
      MXBTP = NOCTPB
    end if
    do IBTP=1,MXBTP
      if (IOCOC(IATP,IBTP) == 0) cycle
      IBSTRT = ISSOB(IBTP,IBSM)
      IBSTOP = IBSTRT+NSSOB(IBTP,IBSM)-1
      do IB=IBSTRT,IBSTOP
        IBREL = IB-IBSTRT+1

        ! Terms depending only on IB

        XB(:) = Zero
        HB = Zero
        RJBB = Zero

        do IEL=1,NBEL
          IBEL = IBSTR(IEL,IB)
          HB = HB+H(IBEL)

          do JEL=1,NBEL
            RJBB = RJBB+RK(IBSTR(JEL,IB),IBEL)
          end do

          do IORB=1,NORB
            XB(IORB) = XB(IORB)+RJ(IORB,IBEL)
          end do
        end do
        EB = HB+Half*RJBB+ECORE

        if ((IREST1 == 1) .and. (IATP == IBTP)) then
          IASTRT = ISSOA(IATP,IASM)-1+IBREL
        else
          IASTRT = ISSOA(IATP,IASM)
        end if
        IASTOP = ISSOA(IATP,IASM)+NSSOA(IATP,IASM)-1
        do IA=IASTRT,IASTOP
          IDET = IDET+1
          ITDET = ITDET+1
          X1 = EB
          X2 = Zero
          do IEL=1,NAEL
            IAEL = IASTR(IEL,IA)
            X1 = X1+(H(IAEL)+XB(IAEL))
            do JEL=1,NAEL
              X2 = X2+RK(IASTR(JEL,IA),IAEL)
            end do
          end do
          DIAG(IDET) = X1+Half*X2
          if (IB == IA) DIAG(IDET) = DIAG(IDET)+XADD
        end do
      end do
      ! Yet a RAS block of the diagonal has been constructed
      if (ICISTR >= 2) then
        IDUM(1) = IDET
        call ITODS(IDUM,1,-1,LUDIA)
        call TODSC_MCLR(DIAG,IDET,-1,LUDIA)
        IDET = 0
      end if
    end do
  end do

end do

if (ICISTR >= 2) then
  IDUM(1) = -1
  call ITODS(IDUM,1,-1,LUDIA)
end if

return
! Avoid unused argument lines
if (.false.) then
  call Unused_real_array(XA)
  call Unused_real_array(SCR)
  call Unused_real(PLSIGN)
end if

end subroutine CIDIA4
