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
! Copyright (C) 1992, Martin Schuetz                                   *
!               1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2017,2022, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SOiniH()
!***********************************************************************
!                                                                      *
!     purpose: generate initial Hessian (diagonal) from                *
!              orbital energies (for second order update)              *
!                                                                      *
!***********************************************************************

use InfSCF, only: FockMO, HDiag, nFro, nOcc, nOrb, nSym, OrbType
use Constants, only: Zero, One, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iD, iOcc, iOff_F, iOff_H, ioffs, iSym, iVir, jOcc, jVir, nD, nOccmF, nOrbmF
real(kind=wp) :: Hii
real(kind=wp), pointer :: Fock(:,:)
real(kind=wp), parameter :: Hii_Max = One, Hii_Min = 0.05_wp

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Set the array to silly large values. In the case of UHF these
! will remain but should not make any difference. They are actully
! needed to make the rs-rfo code work.

! Compute the diagonal values of the Fock matrix, stored in EOrb.
! Call Mk_EOrb(CMO_Ref,Size(CMO_Ref,1),Size(CMO_Ref,2))

nD = size(FockMO,2)
HDiag(:) = Zero

#ifdef _DEBUGPRINT_
write(u6,*) 'nD=',nD
do iD=1,nD
  write(u6,*) 'iD=',iD
  write(u6,'(A,8I3)') 'nOcc',(nOcc(iSym,iD),iSym=1,nSym)
end do
write(u6,'(A,8I3)') 'nFro',(nFro(iSym),iSym=1,nSym)
write(u6,'(A,8I3)') 'nOrb',(nOrb(iSym),iSym=1,nSym)
call RecPrt('SOIniH: FockMO',' ',FockMO,1,size(FockMO))
call NrmClc(FockMO,size(FockMO),'SOIniH','FockMO')
#endif
iOff_H = 1
do iD=1,nD

  iOffs = 0
  iOff_F = 0
  do iSym=1,nSym

    ! loop over all occ orbitals in sym block

    iOffs = iOffs+nFro(iSym)
    ! number of Occupied, excluding frozen
    nOccmF = nOcc(iSym,iD)-nFro(iSym)
    ! number of Orbitals, excluding frozen
    nOrbmF = nOrb(iSym)-nFro(iSym)

    if (nOrb(iSym) < 1) cycle
    Fock(1:nOrb(iSym),1:nOrb(iSym)) => FockMO(iOff_F+1:iOff_F+nOrb(iSym)**2,iD)

    do iOcc=ioffs+1,ioffs+nOccmF
      jOcc = iOcc-iOffs

      ! loop over all virt orbitals in sym block

      do iVir=ioffs+nOccmF+1,ioffs+nOrbmF
        jVir = iVir-iOffs

        if (OrbType(iVir,iD) == OrbType(iOcc,iD)) then

          Hii = Four*(Fock(jVir,jVir)-Fock(jOcc,jOcc))/real(nD,kind=wp)

          !write(u6,*) 'Hii, iOff_H=', Hii, iOff_H
          !write(u6,*) 'Fock(jVir,jVir), jVir=',Fock(jVir,jVir),jVir
          !write(u6,*) 'Fock(jOcc,jOcc), jOcc=',Fock(jOcc,jOcc),jOcc
          if (Hii < Zero) then
            !write(u6,*) 'SOIniH: Hii<0.0, Hii=',Hii
            Hii = max(Hii_Max,abs(Hii))
            !write(u6,*) '        Reset to Hii=',Hii
          else if (abs(Hii) < Hii_Min) then
            !write(u6,*) 'SOIniH: Abs(Hii)<Hii_Min, Hii=',Hii
            !write(u6,*) 'jVir,jOcc=',jVir,jOcc
            !write(u6,*) 'Fock(jOcc,jOcc)=',Fock(jOcc,jOcc)
            !write(u6,*) 'Fock(jVir,jVir)=',Fock(jVir,jVir)
            Hii = Hii_Min
            !write(u6,*) 'SOIniH: Reset to          Hii=',Hii
          end if
          HDiag(iOff_H) = Hii
        end if

        iOff_H = iOff_H+1

      end do  ! iVir

    end do    ! iOcc

    nullify(Fock)
    iOff_F = iOff_F+nOrb(iSym)**2
    iOffs = iOffs+nOrbmF

  end do ! iSym
end do ! iD

#ifdef _DEBUGPRINT_
call RecPrt('HDiag',' ',HDiag(:),1,size(HDiag))
call NrmClc(HDiag,size(HDiag),'SOIniH','HDiag')
#endif

return

end subroutine SOiniH
