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

subroutine rhs_mp2_help2(iSymA,iSymB,iSymI,iSymJ)

use MBPT2_Global, only: ip_Density, ip_Mp2Lagr, ipInt1, ipInt2, ipScr1
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymA, iSymB, iSymI, iSymJ
integer(kind=iwp) :: iA, iAK, iB, iC, iCI, iI, iIC, iJ, iK, iKA, nB, nJ
real(kind=wp) :: Fac, xacbi, xaibc, xiajk, xikja
#include "WrkSpc.fh"
#include "corbinf.fh"
! statement functions
integer(kind=iwp) :: i, j, k, iVirVir, iOccOcc, iMp2Lagr
iVirVir(i,j,k) = nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iOccOcc(i,j,k) = j-1+(nOrb(k)+nDel(k))*(i-1)
iMp2Lagr(i,j,k) = ip_Mp2Lagr(k)+j-1+(nOcc(k)+nFro(k))*(i-1)

!-----------------------------------------------------------------------
!write(u6,*) 'Starting subroutine with symmetries: '
!write(u6,*) 'I,J,A,B',I,J,iA,iB
!write(u6,*) 'NewSym',NewSym
!write(u6,*) 'IpoDens',IpoDens
!-----------------------------------------------------------------------

! Reads the <..|AB> - integrals for the L3-term

nB = nExt(iSymB)+nDel(iSymB)
do iA=1,nExt(iSymA)+nDel(iSymA)
  ! Decide i the AB-matrix is symmetric, if so
  ! only loop over lower triangle
  if (iSymA == iSymB) nB = iA
  do iB=1,nB

    call Exch(iSymI,iSymA,iSymJ,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Work(ipInt1),Work(ipScr1))
    if (iSymA /= iSymB) then
      call Exch(iSymJ,iSymA,iSymI,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Work(ipInt2),Work(ipScr1))
    end if
    !write(u6,*) 'iOrbA',iA,'iOrbB',iB
    !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymA)+nDel(iSymA),nOrb(iSymB)+nDel(iSymA))
    !if (iSymA /= iSymB) then
    !  call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymB)+nDel(iSymB),nOrb(iSymA)+nDel(iSymA))
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
        do iI=1,nOcc(iSymA)+nFro(iSymA)
          xaibc = Work(ipInt1+iIC+iI-1+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (iSymA == iSymB) then
            xacbi = Work(ipInt1+iCI+iC-1+(iI-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xacbi = Work(ipInt2+iCI+iC-1+(iI-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          Work(iMp2Lagr(iA,iI,iSymA)) = Work(iMp2Lagr(iA,iI,iSymA))- &
                                        Fac*Work(ip_Density(iSymB)+iVirVir(iB,iC,iSymB))*(Two*xaibc-xacbi)
          if (iSymA == iSymB) then
            Work(iMp2Lagr(iB,iI,iSymA)) = Work(iMp2Lagr(iB,iI,iSymA))- &
                                          Fac*Work(ip_Density(iSymB)+iVirVir(iA,iC,iSymB))*(Two*xacbi-xaibc)
          end if
          !write(u6,*) 'AIBC',iA,iI,iB,iC
          !write(u6,*) 'AIBC2',iB,iI,iA,iC
          !write(u6,*) 'icab',xaibc
          !write(u6,*) 'ibac',xacbi
          !write(u6,*) 'Dens',Work(ip_Density(iSymB)+iVirVir(iB,iC,iSymB))
          !write(u6,*) 'Dens2',Work(ip_Density(iSymB)+iVirVir(iA,iC,iSymB))
          !write(u6,*) 'A',(Two*xaibc-xacbi)
          !write(u6,*) 'A2',(Two*xacbi-xaibc)
          !write(u6,*) 'Bidrag',Work(ip_Density(iSymB)+iVirVir(iB,iC,iSymB))*(Two*xaibc-xacbi)
          !write(u6,*) 'Bidrag2',Work(ip_Density(iSymB)+iVirVir(iA,iC,iSymB))*(Two*xacbi-xaibc)

        end do
      end do

      ! And again we need another loop for the symmetries that isnt written
      ! explicitly. all references to A and B need to change too.

      if (iSymA /= iSymB) then
        iIC = (nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
        iCI = nFro(iSymA)+nOcc(iSymA)
        do iC=1,nExt(iSymA)+nDel(iSymA)
          do iI=1,nFro(iSymB)+nOcc(iSymB)
            xaibc = Work(ipInt1+iCI+iC-1+(iI-1)*(nOrb(iSymA)+nDel(iSymA)))
            xacbi = Work(ipInt2+iIC+iI-1+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            Work(iMp2Lagr(iB,iI,iSymB)) = Work(iMp2Lagr(iB,iI,iSymB))- &
                                          Fac*Work(ip_Density(iSymA)+iVirVir(iA,iC,iSymA))*(Two*xaibc-xacbi)
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
    call Exch(iSymA,iSymI,iSymB,iSymJ,iI,iJ,Work(ipInt1),Work(ipScr1))
    if (iSymI /= iSymJ) then
      call Exch(iSymB,iSymI,iSymA,iSymJ,iI,iJ,Work(ipInt2),Work(ipScr1))
    end if
    !write(u6,*)
    !write(u6,*) ' *  i,j = ',iI,iJ
    !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymA),nOrb(iSymB))
    !if (iSymI /= iSymJ) then
    !  call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymB),nOrb(iSymA))
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
          xikja = Work(ipInt1+iKA+iK-1+(iA-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymI == iSymJ) then
            xiajk = Work(ipInt1+iAK+iA-1+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xiajk = Work(ipInt2+iAK+iA-1+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          Work(iMp2Lagr(iA,iJ,iSymJ)) = Work(iMp2Lagr(iA,iJ,iSymJ))- &
                                        Fac*Work(ip_Density(iSymI)+iOccOcc(iI,iK,iSymI))*(Two*xikja-xiajk)
          if (iSymI == iSymJ) then
            Work(iMp2Lagr(iA,iI,iSymJ)) = Work(iMp2Lagr(iA,iI,iSymJ))- &
                                          Fac*Work(ip_Density(iSymI)+iOccOcc(iJ,iK,iSymI))*(Two*xiajk-xikja)
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
            xikja = Work(ipInt1+iAK+iA-1+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
            xiajk = Work(ipInt2+iKA+iK-1+(iA-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            Work(iMp2Lagr(iA,iI,iSymI)) = Work(iMp2Lagr(iA,iI,iSymI))- &
                                          Fac*Work(ip_Density(iSymJ)+iOccOcc(iJ,iK,iSymJ))*(Two*xikja-xiajk)
          end do
        end do
      end if
    end if
  end do
end do

return

end subroutine rhs_mp2_help2
