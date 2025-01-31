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
! Copyright (C) 2011, Jeppe Olsen                                      *
!               2011, Giovanni Li Manni                                *
!***********************************************************************

subroutine ADDDIA_TERMS(NAEL,IASTR,NBEL,IBSTR,NORB,CVEC,SVEC,NSMST,H,XA,XB,SCR,RJ,RK,NSSOA,NSSOB,ECORE,IPRNT,NTOOB,RJKAA,IASPGP, &
                        IASM,IBSPGP,IBSM,FACTOR)
!. Update Sigma vector with diagonal terms for a given block
!     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
!                         + FACTOR*DIAG(IASPGP,IBSPGP)CVEC(IASPGP,IBSPGP)
! ========================
! General symmetry version
! ========================
!
! Jeppe Olsen and Giovanni Li Manni, September 2011
!
! I12 = 1 => only one-body part
!     = 2 =>      one+two-body part

implicit real*8(A-H,O-Z)
! Input
dimension NSSOA(NSMST,*), NSSOB(NSMST,*)
dimension H(NORB)
dimension CVEC(*)
! Scratch
dimension RJ(NTOOB,NTOOB), RK(NTOOB,NTOOB)
dimension XA(NORB), XB(NORB), SCR(2*NORB)
dimension IASTR(NAEL,*), IBSTR(NBEL,*)
dimension RJKAA(*)
dimension IDUM(1)
! Output
dimension SVEC(*)

NTEST = 0
NTEST = max(NTEST,IPRNT)
IDUM(1) = 0
I12 = 2
!if (LUIN > 0) rewind LUIN
!if (LUOUT > 0) rewind LUOUT

if (NTEST >= 20) then
  write(6,*) ' ======================'
  write(6,*) ' ADDDIA_TERMS in action'
  write(6,*) ' ======================'
  write(6,*)
  write(6,*) ' IASM, IASPGP, IBSM, IBSPGP = ',IASM,IASPGP,IBSM,IBSPGP
end if

if (NTEST >= 1000) then
  write(6,*) ' Diagonal one electron integrals'
  call WRTMAT(H,1,NORB,1,NORB)
  if (I12 == 2) then
    write(6,*) ' Coulomb and exchange integrals'
    call WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
    write(6,*)
    call WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
  end if
  write(6,*) ' FACTOR = ',FACTOR
end if

!*3 Diagonal elements according to Handys formulae
!   (corrected for error)
!
!   DIAG(IDET) = HII*(NIA+NIB)
!              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,A)*N(J,A)
!              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,B)*N(J,B)
!              +         J(I,J) * N(I,A)*N(J,B)
! N(X) are occupation numbers

! K goes to J - K
if (I12 == 2) call VECSUM(RK,RK,RJ,-1.0d0,+1.0d0,NTOOB**2)

! Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
!                          0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
!
! Obtain alpha strings of sym IASM and type IASPGP
call GETSTR_TOTSM_SPGP(1,IASPGP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM,IDUM)

NIA = NSSOA(IASM,IASPGP)

if (NTEST >= 1000) then
  write(6,*) ' After GETSTR for A strings'
  write(6,*) ' alpha strings obtained'
  call IWRTMA(IASTR,NAEL,NIA,NAEL,NIA)
end if

do IA=1,NIA
  EAA = 0.0d0
  do IEL=1,NAEL
    IAEL = IASTR(IEL,IA)
    EAA = EAA+H(IAEL)
    if (I12 == 2) then
      do JEL=1,NAEL
        EAA = EAA+0.5d0*RK(IASTR(JEL,IA),IAEL)
      end do
    end if
  end do
  RJKAA(IA) = EAA
end do
! Obtain alpha strings of sym IBSM and type IBTP
call GETSTR_TOTSM_SPGP(2,IBSPGP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM,IDUM)
NIB = NSSOB(IBSM,IBSPGP)
IDET = 0
do IB=1,NIB
  ! Terms depending only on IB
  HB = 0.0d0
  RJBB = 0.0d0
  call SETVEC(XB,0.0d0,NORB)

  do IEL=1,NBEL
    IBEL = IBSTR(IEL,IB)
    HB = HB+H(IBEL)

    if (I12 == 2) then
      do JEL=1,NBEL
        RJBB = RJBB+RK(IBSTR(JEL,IB),IBEL)
      end do

      do IORB=1,NORB
        XB(IORB) = XB(IORB)+RJ(IORB,IBEL)
      end do
    end if
  end do

  EB = HB+0.5d0*RJBB+ECORE
  do IA=1,NSSOA(IASM,IASPGP)
    IDET = IDET+1
    X = EB+RJKAA(IA)
    do IEL=1,NAEL
      X = X+XB(IASTR(IEL,IA))
    end do
    SVEC(IDET) = SVEC(IDET)+CVEC(IDET)*(X+FACTOR)
  end do ! IA
end do ! IB

if (NTEST >= 1000) then
  write(6,*) ' Input and output vectord, ADDDIA_TERMS'
  call WRTMAT(CVEC,1,IDET,1,IDET)
  call WRTMAT(SVEC,1,IDET,1,IDET)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(XA)
  call Unused_real_array(SCR)
end if

end subroutine ADDDIA_TERMS
