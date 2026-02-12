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
! Copyright (C) 2006, Per-Olof Widmark                                 *
!***********************************************************************

subroutine Scram(CMO,nSym,nBas,nOrb,ScrFac)
!***********************************************************************
!                                                                      *
! This routine scrambles start orbitals in order to introduce symmetry *
! breaking where it is desirable.                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
! Written: September 2006                                              *
!                                                                      *
!***********************************************************************

use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
real(kind=wp), intent(in) :: ScrFac
integer(kind=iwp) :: iBas, indx, iOff, iOrb, iSeed = 13, iSym, jndx, jOrb
real(kind=wp) :: p, q, u, v
real(kind=wp), external :: Random_Molcas

!----------------------------------------------------------------------*
! Do small rotations                                                   *
!----------------------------------------------------------------------*
iOff = 0
do iSym=1,nSym
  !write(u6,*) 'Scrambling irrep',iSym
  do iOrb=1,nOrb(iSym)-1
    jOrb = iOrb+1
    q = ScrFac*(Two*Random_Molcas(iSeed)-One)
    p = sqrt(One-q*q)
    !write(u6,*) 'q=',q
    do iBas=1,nBas(iSym)
      indx = iOff+(iOrb-1)*nBas(iSym)+iBas
      jndx = iOff+(jOrb-1)*nBas(iSym)+iBas
      u = p*CMO(indx)-q*CMO(jndx)
      v = q*CMO(indx)+p*CMO(jndx)
      CMO(indx) = u
      CMO(jndx) = v
    end do
  end do
  iOff = iOff+nBas(iSym)*nOrb(iSym)
end do
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Scram
