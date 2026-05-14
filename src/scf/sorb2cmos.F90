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

!#define _DEBUGPRINT_
subroutine Sorb2CMOs(CMO,nCMO,nD,Occ,nnB,nBas,nOrb,nSym,OrbType)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMO, nD, nnB, nSym, nBas(nSym), nOrb(nSym)
real(kind=wp), intent(inout) :: CMO(nCMO,nD), Occ(nnB,nD)
integer(kind=iwp), intent(inout) :: OrbType(nnB,nD)
integer(kind=iwp) :: iD, iOff1, iOff2, iOrb, iSym, iTmp, jOrb, kOrb, nOcc
real(kind=wp) :: Occ_i, Occ_j
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iOff, jOff

do iD=1,nD
  iOff = 1
  jOff = 1
  do iSym=1,nSym
    call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
    call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),nBas(iSym),nOrb(iSym))
    iOff = iOff+nOrb(iSym)
    jOff = jOff+nBas(iSym)*nOrb(iSym)
  end do
end do
#endif

! Sort orbitals into
! 1) occupied and virtual orbitals
! 2) within each block sort according to the orbital energy

do iD=1,nD
  !write(u6,*)
  !write(u6,*) 'iD=',iD
  !write(u6,*)
  iOff1 = 0
  iOff2 = 1
  do iSym=1,nSym

    nOcc = 0

    ! Sort first the orbitals according to the occupation numbers.

    do iOrb=1,nOrb(iSym)-1

      Occ_i = Occ(iOff1+iOrb,iD)
      !write(u6,*) 'Occ_i,iOrb=',Occ_i,iOrb
      kOrb = 0
      do jOrb=iOrb+1,nOrb(iSym)
        Occ_j = Occ(iOff1+jOrb,iD)
        !write(u6,*) 'Occ_j,jOrb=',Occ_j,jOrb
        if ((Occ_i == Zero) .and. (Occ_j > Occ_i)) then
          Occ_i = Occ_j
          kOrb = jOrb
        end if
      end do
      !write(u6,*) 'kOrb=',kOrb
      if (kOrb /= 0) then
        iTmp = OrbType(iOff1+iOrb,iD)
        OrbType(iOff1+iOrb,iD) = OrbType(iOff1+kOrb,iD)
        OrbType(iOff1+kOrb,iD) = iTmp

        Occ_i = Occ(iOff1+iOrb,iD)
        Occ(iOff1+iOrb,iD) = Occ(iOff1+kOrb,iD)
        Occ(iOff1+kOrb,iD) = Occ_i
        call DSwap_(nBas(iSym),CMO(iOff2+(iOrb-1)*nBas(iSym),iD),1,CMO(iOff2+(kOrb-1)*nBas(iSym),iD),1)
      end if

      if (Occ(iOff1+iOrb,iD) /= Zero) nOcc = nOcc+1

    end do

    iOff1 = iOff1+nOrb(iSym)
    iOff2 = iOff2+nBas(iSym)*nOrb(iSym)
  end do    ! iSym
end do      ! iD
#ifdef _DEBUGPRINT_
do iD=1,nD
  iOff = 1
  jOff = 1
  do iSym=1,nSym
    call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
    call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),nBas(iSym),nOrb(iSym))
    iOff = iOff+nOrb(iSym)
    jOff = jOff+nBas(iSym)*nOrb(iSym)
  end do
end do
#endif

return

end subroutine Sorb2CMOs
