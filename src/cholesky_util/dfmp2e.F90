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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine DfMP2E(NVar,NOcc,NCore,E,EMin,EMax)
!-----------------------------------------------------------------------
! Function : Define MP2 energy
!-----------------------------------------------------------------------

implicit real*8(A-H,O-Z)
parameter(TWO=2.0D+00)
real*8 E(NVar)
integer IDX(NVar)

do I=1,NVar
  IDX(I) = I
end do
!call OrdV(E,IDX,NVar)

EMin = TWO*(E(IDX(NCore+NOcc+1))-E(IDX(NCore+NOcc)))
EMax = TWO*(E(IDX(NVar))-E(IDX(NCore+1)))

return

end subroutine DfMP2E
