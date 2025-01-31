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

subroutine ICOPMT(MATI,NRI,NCI,MATO,NRO,NCO)
! Copy integer matrix MATI to MATO

! Input
integer MATI(NRI,NCI)
! Output
integer MATO(NRO,NCO)

NREFF = min(NRI,NRO)
NCEFF = min(NCI,NCO)

do IC=1,NCEFF
  do IR=1,NREFF
    MATO(IR,IC) = MATI(IR,IC)
  end do
end do

end subroutine ICOPMT
