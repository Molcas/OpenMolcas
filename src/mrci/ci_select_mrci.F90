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

subroutine CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)

implicit real*8(A-H,O-Z)
dimension AREF(NREF,NREF), CISEL(NREF,*), IROOT(NRROOT)
dimension PLEN(NREF)

if (NSEL == 0) return
! SELECTION BY PROJECTION ONTO SPACE SPANNED BY CISEL VECTORS. IROOT()
! IS SET TO SELECT THE NRROOT VECTORS WITH MAX PROJECTED LENGTH.
do J=1,NREF
  SUM = 0.0d00
  do ISEL=1,NSEL
    SUM1 = 0.0d00
    do I=1,NREF
      SUM1 = SUM1+AREF(I,J)*CISEL(I,ISEL)
    end do
    SUM = SUM+SUM1**2
  end do
  PLEN(J) = SUM+J*1.0D-12
end do
! SELECT BY MAGNITUDE OF PLEN:
do J=1,NRROOT
  PMAX = PLEN(1)
  JMAX = 1
  do JJ=2,NREF
    if (PMAX >= PLEN(JJ)) goto 150
    PMAX = PLEN(JJ)
    JMAX = JJ
150 continue
  end do
  PLEN(JMAX) = -PMAX
end do
I = 0
do IR=1,NREF
  PL = PLEN(IR)
  if (PL < 0.0d00) then
    I = I+1
    IROOT(I) = IR
    PL = -PL
  end if
  PLEN(IR) = PL-IR*1.0D-12
end do

return

end subroutine CI_SELECT_MRCI
