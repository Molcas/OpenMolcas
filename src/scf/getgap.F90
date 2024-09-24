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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

subroutine GetGap(Eorb,nData,nAufb,Gap,Efermi)
!***********************************************************************
!                                                                      *
! This routine figure out the homo lumo gap.                           *
!                                                                      *
!***********************************************************************

use Constants, only: Half

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
integer nData
integer nAufb
real*8 Eorb(nData)
real*8 Gap
real*8 Efermi
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer i, j, k
real*8 tmp

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Sort array                                                           *
!----------------------------------------------------------------------*
!write(6,*) 'Unsorted orbitals energies'
!write(6,'(10F12.6)') Eorb
do i=1,nData-1
  k = i
  do j=i+1,nData
    if (Eorb(j) < Eorb(k)) k = j
  end do
  tmp = Eorb(k)
  Eorb(k) = Eorb(i)
  Eorb(i) = tmp
end do
!Write(6,*) 'Sorted orbitals energies'
!Write(6,'(10F12.6)') Eorb
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (nAufb <= 0) then
  Gap = 1.0d3
  Efermi = Eorb(1)
else if (nAufb >= nData) then
  Gap = 1.0d3
  Efermi = Eorb(nData)+1.0d-3
else
  Gap = Eorb(nAufb+1)-Eorb(nAufb)
  Efermi = half*(Eorb(nAufb+1)+Eorb(nAufb))
end if
!write(6,*) 'Gap:',Gap
!write(6,*) 'Efermi:',Efermi
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine GetGap
