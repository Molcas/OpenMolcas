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

subroutine GT1DIS(H1DIA,IREOTS,IPNT,H,ISMFTO,IBSO,NACOB)
! diagonal of one electron integrals over active orbitals

implicit real*8(A-H,O-Z)
! Input
integer IREOTS(*), IPNT(*), ISMFTO(*), IBSO(*)
dimension H(*)
! Output
dimension H1DIA(*)

do IIOB=1,NACOB
  IOB = IREOTS(IIOB)
  ISM = ISMFTO(IIOB)
  IOBREL = IOB-IBSO(ISM)+1
  !write(6,*) ' IIOB IOB ISM IOBREL'
  !write(6,*) IIOB,IOB,ISM,IOBREL
  H1DIA(IIOB) = H(IPNT(ISM)-1+IOBREL*(IOBREL+1)/2)
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' Diagonal one electron integrals'
  call WRTMAT(H1DIA,1,NACOB,1,NACOB)
end if

end subroutine GT1DIS
