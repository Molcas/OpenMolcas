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
! Copyright (C) 1997,1998, Jeppe Olsen                                 *
!***********************************************************************

subroutine GETH1(H,ISM,ITP,JSM,JTP)
! One-electron integrals over orbitals belonging to
! given OS class
!
! The orbital symmetries  are used to obtain the total
! symmetry of the one-electron integrals.
! It is therefore assumed that ISM, JSM represents a correct symmetry block
! of the integrals
!
! Jeppe Olsen, Version of fall 97
!              Summer of 98 : CC options added

use lucia_data, only: NOBPTS

implicit none
integer ISM, ITP, JSM, JTP
! Output
real*8 H(*)
integer NI, NJ, IJ, J, I, NTEST
real*8, external :: GETH1E

NI = NOBPTS(ITP,ISM)
NJ = NOBPTS(JTP,JSM)

! Normal one-electron integrals

IJ = 0
do J=1,NJ
  do I=1,NI
    IJ = IJ+1
    H(IJ) = GETH1E(I,ITP,ISM,J,JTP,JSM)
  end do
end do

NTEST = 0
if (NTEST /= 0) then
  write(6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
  call WRTMAT(H,NI,NJ,NI,NJ)
end if

end subroutine GETH1
