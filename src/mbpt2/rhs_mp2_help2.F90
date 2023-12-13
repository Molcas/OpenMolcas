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

subroutine rhs_mp2_help2(iSymA,iSymB,iSymI,iSymJ,Int1,Int2,Scr1)

#include "intent.fh"

use MBPT2_Global, only: Density, Mp2Lagr
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymA, iSymB, iSymI, iSymJ
real(kind=wp), intent(_OUT_) :: Int1(*), Int2(*), Scr1(*)
integer(kind=iwp) :: iA, iAA, iAK, iB, iBB, iC, iCC, iCI, iI, iIC, iJ, iK, iKA, nB, nJ
real(kind=wp) :: Fac, xacbi, xaibc, xiajk, xikja
#include "corbinf.fh"

!-----------------------------------------------------------------------
!write(u6,*) 'Starting subroutine with symmetries: '
!write(u6,*) 'I,J,A,B',I,J,iA,iB
!write(u6,*) 'NewSym',NewSym
!write(u6,*) 'IpoDens',IpoDens
!-----------------------------------------------------------------------

! Reads the <..|AB> - integrals for the L3-term

nB = nExt(iSymB)+nDel(iSymB)
do iA=1,nExt(iSymA)+nDel(iSymA)
  iAA = iA+nFro(iSymA)+nOcc(iSymA)
  ! Decide i the AB-matrix is symmetric, if so
  ! only loop over lower triangle
  if (iSymA == iSymB) nB = iA
  do iB=1,nB
    iBB = iB+nFro(iSymB)+nOcc(iSymB)

    call Exch(iSymI,iSymA,iSymJ,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Int1,Scr1)
    if (iSymA /= iSymB) then
      call Exch(iSymJ,iSymA,iSymI,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Int2,Scr1)
    end if
    !write(u6,*) 'iOrbA',iA,'iOrbB',iB
    !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymA)+nDel(iSymA),nOrb(iSymB)+nDel(iSymA))
    !if (iSymA /= iSymB) then
    !  call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSymB)+nDel(iSymB),nOrb(iSymA)+nDel(iSymA))
    !end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Lagrangian term 3
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if ((iSymA == iSymI) .and. (iSymB == iSymJ)) then
      if ((iA == iB) .and. (iSymA == iSymB)) then
        Fac = Half
      else
        Fac = One
      end if
      iIC = (nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
      iCI = nFro(iSymB)+nOcc(iSymB)
      do iC=1,nExt(iSymB)+nDel(iSymB)
        iCC = iC+nFro(iSymB)+nOcc(iSymB)
        do iI=1,nOcc(iSymA)+nFro(iSymA)
          xaibc = Int1(iIC+iI+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (iSymA == iSymB) then
            xacbi = Int1(iCI+iC+(iI-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xacbi = Int2(iCI+iC+(iI-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          Mp2Lagr%SB(iSymA)%A2(iI,iA) = Mp2Lagr%SB(iSymA)%A2(iI,iA)-Fac*Density%SB(iSymB)%A2(iCC,iBB)*(Two*xaibc-xacbi)
          if (iSymA == iSymB) then
            Mp2Lagr%SB(iSymA)%A2(iI,iB) = Mp2Lagr%SB(iSymA)%A2(iI,iB)-Fac*Density%SB(iSymB)%A2(iCC,iAA)*(Two*xacbi-xaibc)
          end if
          !write(u6,*) 'AIBC',iA,iI,iB,iC
          !write(u6,*) 'AIBC2',iB,iI,iA,iC
          !write(u6,*) 'icab',xaibc
          !write(u6,*) 'ibac',xacbi
          !write(u6,*) 'Dens',Density%SB(iSymB)%A2(iCC,iBB)
          !write(u6,*) 'Dens2',Density%SB(iSymB)%A2(iCC,iAA)
          !write(u6,*) 'A',(Two*xaibc-xacbi)
          !write(u6,*) 'A2',(Two*xacbi-xaibc)
          !write(u6,*) 'Bidrag',Density%SB(iSymB)%A2(iCC,iBB)*(Two*xaibc-xacbi)
          !write(u6,*) 'Bidrag2',Density%SB(iSymB)%A2(iCC,iAA)*(Two*xacbi-xaibc)

        end do
      end do

      ! And again we need another loop for the symmetries that isnt written
      ! explicitly. all references to A and B need to change too.

      if (iSymA /= iSymB) then
        iIC = (nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
        iCI = nFro(iSymA)+nOcc(iSymA)
        do iC=1,nExt(iSymA)+nDel(iSymA)
          iCC = iC+nFro(iSymA)+nOcc(iSymA)
          do iI=1,nFro(iSymB)+nOcc(iSymB)
            xaibc = Int1(iCI+iC+(iI-1)*(nOrb(iSymA)+nDel(iSymA)))
            xacbi = Int2(iIC+iI+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            Mp2Lagr%SB(iSymB)%A2(iI,iB) = Mp2Lagr%SB(iSymB)%A2(iI,iB)-Fac*Density%SB(iSymA)%A2(iCC,iAA)*(Two*xaibc-xacbi)
          end do
        end do
      end if
    end if

  end do
end do
! Reads one block of <ij|..> at the time.

nJ = nFro(iSymJ)+nOcc(iSymJ)
do iI=1,nFro(iSymI)+nOcc(iSymI)
  ! Decide if the IJ-block i symmetric and
  ! if so only read a triangle
  if (iSymI == iSymJ) nJ = iI
  !***** Check if the IJ-matrix is symmetric, if so
  !      only loop over lower triangle.
  do iJ=1,nJ
    ! If we are only looping over a triangular matrix
    ! we need to count each offdiagonal element twice.
    ! Copy one AB-block of integrals to the workspace
    call Exch(iSymA,iSymI,iSymB,iSymJ,iI,iJ,Int1,Scr1)
    if (iSymI /= iSymJ) then
      call Exch(iSymB,iSymI,iSymA,iSymJ,iI,iJ,Int2,Scr1)
    end if
    !write(u6,*)
    !write(u6,*) ' *  i,j = ',iI,iJ
    !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymA),nOrb(iSymB))
    !if (iSymI /= iSymJ) then
    !  call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSymB),nOrb(iSymA))
    !end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Lagrangian term 4
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if ((iSymI == iSymA) .and. (iSymJ == iSymB)) then
      if ((iI == iJ) .and. (iSymI == iSymJ)) then
        Fac = Half
      else
        Fac = One
      end if

      iKA = (nFro(iSymJ)+nOcc(iSymJ))*(nOrb(iSymI)+nDel(iSymI))
      iAK = nFro(iSymJ)+nOcc(iSymJ)
      do iK=1,nFro(iSymI)+nOcc(iSymI)
        do iA=1,nExt(iSymJ)+nDel(iSymJ)
          xikja = Int1(iKA+iK+(iA-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymI == iSymJ) then
            xiajk = Int1(iAK+iA+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xiajk = Int2(iAK+iA+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          Mp2Lagr%SB(iSymJ)%A2(iJ,iA) = Mp2Lagr%SB(iSymJ)%A2(iJ,iA)-Fac*Density%SB(iSymI)%A2(iK,iI)*(Two*xikja-xiajk)
          if (iSymI == iSymJ) then
            Mp2Lagr%SB(iSymJ)%A2(iI,iA) = Mp2Lagr%SB(iSymJ)%A2(iI,iA)-Fac*Density%SB(iSymI)%A2(iK,iJ)*(Two*xiajk-xikja)
          end if
          !write(u6,*) 'xikja',xikja
          !write(u6,*) 'xiajk',xiajk
        end do
      end do

      if (iSymI /= iSymJ) then
        iKA = (nFro(iSymI)+nOcc(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))
        iAK = nFro(iSymI)+nOcc(iSymI)
        do iK=1,nFro(iSymJ)+nOcc(iSymJ)
          do iA=1,nExt(iSymI)+nDel(iSymI)
            xikja = Int1(iAK+iA+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
            xiajk = Int2(iKA+iK+(iA-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            Mp2Lagr%SB(iSymI)%A2(iI,iA) = Mp2Lagr%SB(iSymI)%A2(iI,iA)-Fac*Density%SB(iSymJ)%A2(iK,iJ)*(Two*xikja-xiajk)
          end do
        end do
      end if
    end if
  end do
end do

return

end subroutine rhs_mp2_help2
