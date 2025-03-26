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

subroutine CMP2CN(ICNL,NCLL,NOPL,ICNR,NCLR,NOPR,ISCR,NORB,NDIFF)
! Number of differences in occupation of two configurations

dimension ICNL(*), ICNR(*)
dimension ISCR(*)

! Length of Scratch: Number of orbitals

ISCR(1:NORB) = 0
do ICL=1,NCLL
  ISCR(ICNL(ICL)) = 2
end do
do IOP=1,NOPL
  ISCR(ICNL(NCLL+IOP)) = 1
end do

NDIFF = 0
do ICL=1,NCLR
  NDIFF = NDIFF+2-ISCR(ICNR(ICL))
end do
do IOP=1,NOPR
  if (ISCR(ICNR(NCLR+IOP)) == 0) NDIFF = NDIFF+1
end do

return

end subroutine CMP2CN
