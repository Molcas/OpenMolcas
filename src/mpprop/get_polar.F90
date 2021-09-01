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

subroutine get_polar(nPrim,nBas,nAtoms,NOCOB,OENE,nOrb,OCOF,RCHC,LNearestAtom,LFirstRun)
!EB  OENE,ONUM,nOrb,OCOF,RCHC,LNearestAtom)

use MPProp_globals, only: AtPol, AtBoPol, BondMat, Cor, CordMltPl, Frac, iAtPrTab, nAtomPBas, MltPl
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPrim, nBas, nAtoms, NOCOB, nOrb
real(kind=wp), intent(in) :: OENE(nOrb), OCOF(nBas,nPrim), RCHC(3,nBas)
logical(kind=iwp), intent(in) :: LNearestAtom, LFirstRun
integer(kind=iwp) :: i, iA, iPBas, iStdOut, j, K, KK, kl, L, LL, nA, nB
real(kind=wp) :: FOE, FracA, FracB, PAX, PAY, PAZ, R, RA, RB, RIJX, RIJY, RIJZ, Smallest
real(kind=wp), allocatable :: Pol(:,:,:), Pd(:,:)

iStdOut = u6
iA = 0 ! Added by EB

write(iStdOut,*) '  '
write(iStdOut,*) ' CALCULATE THE POLARIZATION TENSOR '
write(iStdOut,*) '  '

! CONSTRUCT THE POLARIZATION CONTRIBUTION FROM EACH PAIR OF ATOMS

call mma_allocate(Pd,3,nAtoms,label='Pd')
call mma_allocate(Pol,6,nAtoms,nAtoms,label='Pol')
Pol(:,:,:) = Zero

! NBOND assigns the array value where bondvaluearray starts
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
          kl = kk*(kk-1)/2+ll
          PAX = PAX+OCOF(I,K)*OCOF(J,L)*(MltPl(1)%A(kl,1)+MltPl(0)%A(kl,1)*(CordMltPl(1,1)-RIJX))
          PAY = PAY+OCOF(I,K)*OCOF(J,L)*(MltPl(1)%A(kl,2)+MltPl(0)%A(kl,1)*(CordMltPl(2,1)-RIJY))
          PAZ = PAZ+OCOF(I,K)*OCOF(J,L)*(MltPl(1)%A(kl,3)+MltPl(0)%A(kl,1)*(CordMltPl(3,1)-RIJZ))
        end do
      end do
      PD(1,nA) = PAX
      PD(2,nA) = PAY
      PD(3,nA) = PAZ
    end do
    do nA=1,nAtoms
      do nB=1,nAtoms
        Pol(1,nA,nB) = Pol(1,nA,nB)+Pd(1,nA)*Pd(1,nB)*FOE
        Pol(2,nA,nB) = Pol(2,nA,nB)+Pd(1,nA)*Pd(2,nB)*FOE
        Pol(3,nA,nB) = Pol(3,nA,nB)+Pd(1,nA)*Pd(3,nB)*FOE
        Pol(4,nA,nB) = Pol(4,nA,nB)+Pd(2,nA)*Pd(2,nB)*FOE
        Pol(5,nA,nB) = Pol(5,nA,nB)+Pd(2,nA)*Pd(3,nB)*FOE
        Pol(6,nA,nB) = Pol(6,nA,nB)+Pd(3,nA)*Pd(3,nB)*FOE
      end do
    end do
  end do
end do
do nA=1,nAtoms
  AtBoPol(:,nA*(nA+1)/2) = Pol(:,nA,nA)
  AtPol(:,nA) = Pol(:,nA,nA)
  do nB=1,nA-1
    AtBoPol(:,nA*(nA-1)/2+nB) = Pol(:,nA,nB)+Pol(:,nB,nA)
  end do
end do
call mma_deallocate(Pd)
call mma_deallocate(Pol)
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
      AtPol(:,iA) = AtPol(:,iA)+AtBoPol(:,nA*(nA-1)/2+nB)*FracA
      AtPol(:,nB) = AtPol(:,nB)+AtBoPol(:,nA*(nA-1)/2+nB)*FracB
    else
      AtPol(:,iA) = AtPol(:,iA)+AtBoPol(:,nA*(nA-1)/2+nB)*FracA
      AtPol(:,nB) = AtPol(:,nB)+AtBoPol(:,nA*(nA-1)/2+nB)*FracB
      AtBoPol(:,iA*(iA+1)/2) = AtBoPol(:,iA*(iA+1)/2)+AtBoPol(:,nA*(nA-1)/2+nB)*FracA
      AtBoPol(:,nB*(nB+1)/2) = AtBoPol(:,nB*(nB+1)/2)+AtBoPol(:,nA*(nA-1)/2+nB)*FracB
    end if
  end do
end do
! If a UHF system is used, correct for the double counting of alpha and beta electrons
if (.not. LFirstRun) then
  FracA = Half
  do nA=1,nAtoms
    AtBoPol(:,nA*(nA+1)/2) = AtBoPol(:,nA*(nA+1)/2)*FracA
    AtPol(:,nA) = AtPol(:,nA)*FracA
    do nB=1,nA-1
      if (BondMat(nA,nB)) then
        AtBoPol(:,nA*(nA-1)/2+nB) = AtBoPol(:,nA*(nA-1)/2+nB)*FracA
      end if
    end do
  end do
end if

return

end subroutine get_polar
