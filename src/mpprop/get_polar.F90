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

subroutine get_polar(nPrim,nBas,nAtoms,nCenters,NOCOB,OENE,nOrb,OCOF,RCHC,LNearestAtom,LFirstRun)
!EB  OENE,ONUM,nOrb,OCOF,RCHC,LNearestAtom)

use Constants, only: Zero, One, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPrim, nBas, nAtoms, nCenters, NOCOB, nOrb, OCOF(nBas,nPrim), RCHC(3,nBas)
real(kind=wp), intent(in) :: OENE(nOrb)
logical(kind=iwp), intent(in) :: LNearestAtom, LFirstRun
#include "WrkSpc.fh"
#include "MpData.fh"
#include "MolProp.fh"
integer(kind=iwp) :: i, iA, iPBas, iStdOut, j, K, KK, L, LL, nA, nB
real(kind=wp) :: FOE, FracA, FracB, PAX, PAY, PAZ, Pol(6,nAtoms,nAtoms), Pd(3,nAtoms), R, RA, RB, RIJX, RIJY, RIJZ, Smallest

iStdOut = u6
iA = 0 ! Added by EB

write(iStdOut,*) '  '
write(iStdOut,*) ' CALCULATE THE POLARIZATION TENSOR '
write(iStdOut,*) '  '

! CONSTRUCT THE POLARIZATION CONTRIBUTION FROM EACH PAIR OF ATOMS

! NBOND assigns the array value where bondvaluearray starts
do nA=1,nAtoms
  do nB=1,nAtoms
    do i=1,6
      POL(i,nA,nB) = Zero
    end do
  end do
end do
write(iStdOut,*)
write(iStdOut,*) 'No occupied orbitals',NOCOB
write(iStdOut,*) 'No orbitals',nBas
write(iStdOut,*)

do i=1,NOCOB
  do j=NOCOB+1,nBas

    ! ORBITAL ENERGIES

    FOE = Four/(OENE(J)-OENE(I))

    ! CALCULATE EXPANSION CENTER FOR THE TWO ORBITALS

    RIJX = Half*(RCHC(1,I)+RCHC(1,J))
    RIJY = Half*(RCHC(2,I)+RCHC(2,J))
    RIJZ = Half*(RCHC(3,I)+RCHC(3,J))
    do nA=1,nAtoms
      PAX = Zero
      PAY = Zero
      PAZ = Zero
      do iPBas=1,nAtomPBas(NA)
        K = iAtPrTab(iPBas,NA)
        do L=1,nPrim
          if (L > K) then
            KK = L
            LL = K
          else
            KK = K
            LL = L
          end if
          PAX = PAX+OCOF(I,K)*OCOF(J,L)*(work(iwork(iMltPlAd(1))+kk*(kk-1)/2+ll-1)+work(iwork(iMltPlAd(0))+ &
                kk*(kk-1)/2+ll-1)*(CordMltPl(1,1)-RIJX))

          PAY = PAY+OCOF(I,K)*OCOF(J,L)*(work(iwork(iMltPlAd(1)+1)+kk*(kk-1)/2+ll-1)+work(iwork(iMltPlAd(0))+ &
                kk*(kk-1)/2+ll-1)*(CordMltPl(2,1)-RIJY))

          PAZ = PAZ+OCOF(I,K)*OCOF(J,L)*(work(iwork(iMltPlAd(1)+2)+kk*(kk-1)/2+ll-1)+work(iwork(iMltPlAd(0))+ &
                kk*(kk-1)/2+ll-1)*(CordMltPl(3,1)-RIJZ))
        end do
      end do
      PD(1,nA) = PAX
      PD(2,nA) = PAY
      PD(3,nA) = PAZ
    end do
    do nA=1,nAtoms
      do nB=1,nAtoms
        POL(1,nA,nB) = POL(1,nA,nB)+PD(1,nA)*PD(1,nB)*FOE
        POL(2,nA,nB) = POL(2,nA,nB)+PD(1,nA)*PD(2,nB)*FOE
        POL(3,nA,nB) = POL(3,nA,nB)+PD(1,nA)*PD(3,nB)*FOE
        POL(4,nA,nB) = POL(4,nA,nB)+PD(2,nA)*PD(2,nB)*FOE
        POL(5,nA,nB) = POL(5,nA,nB)+PD(2,nA)*PD(3,nB)*FOE
        POL(6,nA,nB) = POL(6,nA,nB)+PD(3,nA)*PD(3,nB)*FOE
      end do
    end do
  end do
end do
do nA=1,nAtoms
  do i=1,6
    Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1) = POL(i,nA,nA)
    Work(iAtPolAd+nAtoms*(i-1)+nA-1) = POL(i,nA,nA)
  end do
  do nB=1,nA-1
    do i=1,6
      Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1) = POL(i,nA,nB)+POL(i,nB,nA)
    end do
  end do
end do
do nA=1,nAtoms
  do nB=1,nA-1
    FracA = Frac(nA,nB)
    FracB = One-FracA
    ! Find the closest atom
    if (LNearestAtom .and. (.not. BondMat(nA,nB))) then
      Smallest = huge(Smallest)
      do i=1,nAtoms
        R = sqrt((Cor(1,nA,nB)-Cor(1,i,i))**2+(Cor(2,nA,nB)-Cor(2,i,i))**2+(Cor(3,nA,nB)-Cor(3,i,i))**2)
        if (R < Smallest) then
          Smallest = R
          iA = i
        end if
      end do
      RA = sqrt((Cor(1,nA,nB)-Cor(1,nA,nA))**2+(Cor(2,nA,nB)-Cor(2,nA,nA))**2+(Cor(3,nA,nB)-Cor(3,nA,nA))**2)
      RB = sqrt((Cor(1,nA,nB)-Cor(1,nB,nB))**2+(Cor(2,nA,nB)-Cor(2,nB,nB))**2+(Cor(3,nA,nB)-Cor(3,nB,nB))**2)
      if (((iA /= nA) .and. (iA /= nB)) .and. ((Smallest < RA) .and. (Smallest < RB))) then
        iA = iA
        FracA = One
        FracB = Zero
        !write(iStdOut,*)
        !write(iStdOut,*) ' Moving polaizability between the atoms', nA, nB
        !write(iStdOut,*) ' to the atom                  ', iA
        !write(iStdOut,*)
      else
        iA = nA
      end if
    else
      iA = nA
    end if
    if (BondMat(nA,nB)) then
      do i=1,6
        Work(iAtPolAd+nAtoms*(i-1)+iA-1) = Work(iAtPolAd+nAtoms*(i-1)+iA-1)+Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
        Work(iAtPolAd+nAtoms*(i-1)+nB-1) = Work(iAtPolAd+nAtoms*(i-1)+nB-1)+Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
      end do
    else
      do i=1,6
        Work(iAtPolAd+nAtoms*(i-1)+iA-1) = Work(iAtPolAd+nAtoms*(i-1)+iA-1)+Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
        Work(iAtPolAd+nAtoms*(i-1)+nB-1) = Work(iAtPolAd+nAtoms*(i-1)+nB-1)+Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
        Work(iAtBoPolAd+nCenters*(i-1)+iA*(iA+1)/2-1) = Work(iAtBoPolAd+nCenters*(i-1)+iA*(iA+1)/2-1)+ &
                                                        Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
        Work(iAtBoPolAd+nCenters*(i-1)+nB*(nB+1)/2-1) = Work(iAtBoPolAd+nCenters*(i-1)+nB*(nB+1)/2-1)+ &
                                                        Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
      end do
    end if
  end do
end do
! If a UHF system is used, correct for the double counting of alpha and beta electrons
if (.not. LFirstRun) then
  FracA = Half
  do nA=1,nAtoms
    do i=1,6
      Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1) = Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1)*FracA
      Work(iAtPolAd+nAtoms*(i-1)+nA-1) = Work(iAtPolAd+nAtoms*(i-1)+nA-1)*FracA
    end do
    do nB=1,nA-1
      if (BondMat(nA,nB)) then
        do i=1,6
          Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1) = Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
        end do
      end if
    end do
  end do
end if

return

end subroutine get_polar
