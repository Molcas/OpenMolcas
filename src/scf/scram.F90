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

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 CMO(*)
integer nSym
integer nBas(nSym)
integer nOrb(nSym)
real*8 ScrFac
!----------------------------------------------------------------------*
! External references                                                  *
!----------------------------------------------------------------------*
real*8 Random_Molcas
external Random_Molcas
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer iSeed
save iSeed
integer iSym
integer iOrb
integer jOrb
integer iBas
integer iOff
integer indx
integer jndx
real*8 p
real*8 q
real*8 u
real*8 v
data iSeed/13/

!----------------------------------------------------------------------*
! Do small rotations                                                   *
!----------------------------------------------------------------------*
iOff = 0
do iSym=1,nSym
  !write(6,*) 'Scrambling irrep',iSym
  do iOrb=1,nOrb(iSym)-1
    jOrb = iOrb+1
    q = ScrFac*(Two*Random_Molcas(iSeed)-One)
    p = sqrt(One-q*q)
    !write(6,*) 'q=',q
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
