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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine yHx(X,Y,nXY)
!***********************************************************************
!                                                                      *
!     purpose: multiply an approximation of the orbital Hessian times  *
!              a trial vector, X.                                      *
!                                                                      *
!     Y = H(approximate) x  X                                          *
!                                                                      *
!***********************************************************************

use InfSCF, only: FockMO, nFro, nOcc, nOrb, nSym, OrbType
use Constants, only: Zero, One, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nXY
real(kind=wp), target, intent(in) :: X(nXY)
real(kind=wp), target, intent(out) :: Y(nXY)
integer(kind=iwp) :: iD, iOcc, iOff_F, iOff_XY, iSym, iVir, jOcc, jVir, nD, nOccmF, nOrbmF
real(kind=wp) :: Hij, Tmp
real(kind=wp), parameter :: Hii_Max = One, Hii_Min = 0.05_wp
real(kind=wp), pointer :: Fock(:,:), XP(:,:), YP(:,:)

!----------------------------------------------------------------------*

!write(u6,*)
!call NrmClc(FockMO,SIZE(FockMO),'yHx','FockMO(:,:)')
!call NrmClc(X,SIZE(X),'yHx','X(:)')
!call RecPrt('yHx: FockMO',' ',FockMO,Size(FockMO,1),Size(FockMO,2))
!call RecPrt('yHx: X',' ',X,1,Size(X))

nD = size(FockMO,2)
iOff_XY = 0
do iD=1,nD

  iOff_F = 0
  do iSym=1,nSym

    ! loop over all occ orbitals in sym block

    ! number of Occupied, excluding frozen
    nOccmF = nOcc(iSym,iD)-nFro(iSym)
    ! number of Orbitals, excluding frozen
    nOrbmF = nOrb(iSym)-nFro(iSym)

    Fock(1:nOrb(iSym),1:nOrb(iSym)) => FockMO(iOff_F+1:iOff_F+nOrb(iSym)**2,iD)
    XP(nOccmF+1:nOrbmF,1:nOccmF) => X(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))
    YP(nOccmF+1:nOrbmF,1:nOccmF) => Y(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))

    do iOcc=1,nOccmF
      do iVir=nOccmF+1,nOrbmF

        Tmp = Zero
        do jOcc=1,nOccmF
          do jVir=nOccmF+1,nOrbmF

            Hij = Zero
            if ((OrbType(iVir,iD) == OrbType(iOcc,iD)) .and. (OrbType(jVir,iD) == OrbType(jOcc,iD)) .and. &
                (OrbType(iVir,iD) == OrbType(jOcc,iD))) then

              if ((iVir == jVir) .and. (iOcc == jOcc)) then
                Hij = (Four*(Fock(iVir,jVir)-Fock(iOcc,jOcc))/real(nD,kind=wp))

                if (Hij < Zero) then
                  !write (u6,*) 'Hii<0.0, Hii=',Hij
                  Hij = max(Hii_Max,abs(Hij))
                else if (abs(Hij) < Hii_Min) then
                  !write(u6,*) 'Abs(Hii)<0.05, Hii=',Hij
                  !write(u6,*) 'jVir,jOcc=',jVir,jOcc
                  !write(u6,*) 'Fock(jOcc,jOcc)=',Fock(jOcc,jOcc)
                  !write(u6,*) 'Fock(jVir,jVir)=',Fock(jVir,jVir)
                  Hij = Hii_Min
                end if

              else if ((iVir == jVir) .and. (iOcc /= jOcc)) then
                Hij = (Four*(-Fock(iOcc,jOcc))/real(nD,kind=wp))
              else if ((iOcc == jOcc) .and. (iVir /= jVir)) then
                Hij = (Four*(Fock(iVir,jVir))/real(nD,kind=wp))
              end if
              Tmp = Tmp+Hij*XP(jVir,jOcc)

            end if

          end do  ! jVir
        end do    ! jOcc
        YP(iVir,iOcc) = Tmp

      end do  ! iVir
    end do    ! iOcc

    nullify(Fock,XP,YP)
    iOff_XY = iOff_XY+nOccmF*(nOrbmF-nOccmF)
    iOff_F = iOff_F+nOrb(iSym)**2

  end do ! iSym
end do ! iD
#ifdef _DEBUGPRINT_
call NrmClc(Y,size(Y),'yHx','Y(:)')
call RecPrt('yHx: Y',' ',Y,1,size(Y))
#endif

#undef _DEBUGPRINT_
end subroutine yHx
