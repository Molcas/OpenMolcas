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

subroutine rhs_mp2_help1(iSymA,iSymB,iSymI,iSymJ,Int1,Int2,Scr1)

#include "intent.fh"

use MBPT2_Global, only: Density, EMP2, EOcc, EVir, mAdDel, mAdFro, mAdOcc, mAdVir, Mp2Lagr, VECL2, WDensity
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymA, iSymB, iSymI, iSymJ
real(kind=wp), intent(_OUT_) :: Int1(*), Int2(*), Scr1(*)
integer(kind=iwp) :: iA, iAA, iAB, iAC, iB, iBB, iBA, iBC, iBDel, iBK, iC, iCA, iCB, iCJ, iI, iII, iIJ, iIK, iJ, iJC, iJFroz, iJI, &
                     iJJ, iJK, iK, iKB, iKI, iKJ, nB, nJ
real(kind=wp) :: EDenom, EDiff, EDiffac, EDiffbc, EDiffik, EDiffjk, fac_ab, fac_ij, T_ab, T_ba, T_ij, T_ji, Tij, Tji, xacbj, &
                 xaibj, xajbc, xajbi, xbcaj, xbiaj, xiajb, xiajc, xiakb, xibja, xibjc, xibjk, xicja, xicjb, xikjb, xjakb, xkaib, &
                 xkajb
#include "corbinf.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! First we read in blocks of all integrals of the form <I.|J.>
! If the I and J is of the same symmetry we only need to
! consider blocks with OrbJ < OrbI
nJ = nOcc(iSymJ)
do iI=1,nOcc(iSymI)
  ! Check if we need to consider all combinations of orbitals
  if (iSymI == iSymJ) nJ = iI
  do iJ=1,nJ
    !write(u6,*) 'iSymI, iSymJ',iSymI,iSymJ
    !write(u6,*) 'iI, iJ',iI,iJ
    ! If only blocks with OrbJ < OrbI is used we need
    ! to multiply the result with two for offdiag terms.
    fac_ij = One
    if ((iI /= iJ) .and. (iSymI == iSymJ)) fac_ij = Two
    ! If I is not equal J there is a symmetry J I B A
    ! that is identical and not stored.
    if (iSymI /= iSymJ) fac_ij = Two
    ! Needed?

    ! Copy one AB-block of integrals to the workspace
    call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nFro(iSymI),iJ+nFro(iSymJ),Int1,Scr1)
    if (iSymI /= iSymJ) then
      call Exch(iSymB,iSymI,iSymA,iSymJ,iI+nFro(iSymI),iJ+nFro(iSymJ),Int2,Scr1)
    end if
    !write(u6,*) 'Rhs_mp2_help1'
    !write(u6,*) ' *  i,j = ',iI,iJ
    !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymA)+nDel(iSymA),nOrb(iSymB)+nDel(iSymB))
    !if (iSymI /= iSymJ) then
    !  call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSymB)+nDel(iSymB),nOrb(iSymA)+nDel(iSymA))
    !end if
    !write(u6,*) ''

    ! Calculates the EMP2-energy as a test
    iAB = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
    iBA = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))

    do iA=1,nExt(iSymA)
      nB = nExt(iSymB)
      if (iSymA == iSymB) nB = iA
      do iB=1,nB
        EDiff = One/(EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1))

        xiajb = Int1(iAB+iA+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (iSymI == iSymJ) then
          xibja = Int1(iAB+iB+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xibja = Int2(iBA+iB+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        end if
        EMP2 = EMP2+Ediff*Fac_ij*(xiajb*(Two*xiajb-xibja))
        VECL2 = VECL2+Ediff**2*Fac_ij*(xiajb*(Two*xiajb-xibja))
        if (((iA /= iB) .and. (iSymA == iSymB)) .or. (iSymI /= iSymJ)) then
          EMP2 = EMP2+Ediff*Fac_ij*(xibja*(Two*xibja-xiajb))
          VECL2 = VECL2+Ediff**2*Fac_ij*(xibja*(Two*xibja-xiajb))
        end if
      end do !iOrbB
    end do   !iOrbA
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! P Virtual-Virtual
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Calculate the virtual-virtual perturbation to the
    ! P-matrix
    ! In this and all other things we calculate from now on
    ! we need to do an extra loop over symmetries where J is
    ! not equal to I since only symmetries where I > J is saved.
    ! In the extra loops we do some special handling to get the
    ! results from the symmetries with I < J.

    iAC = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
    iCA = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
    iBC = iAC
    iCB = iCA
    do iC=1,nExt(iSymB)
      do iA=1,nExt(iSymA)
        iAA = iA+nFro(iSymA)+nOcc(iSymA)
        xiajc = Int1(iAC+iA+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (ISymI == iSymJ) then
          xicja = Int1(iCA+iC+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xicja = Int2(iCA+iC+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        end if
        EDiffac = EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iC-1)
        T_ij = (Two*xiajc-xicja)/EDiffac
        T_ji = (Two*xicja-xiajc)/EDiffac
        do iB=1,nExt(iSymA)
          iBB = iB+nFro(iSymA)+nOcc(iSymA)
          xibjc = Int1(iBC+iB+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (ISymI == iSymJ) then
            xicjb = Int1(iCB+iC+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xicjb = Int2(iCB+iC+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          EDiffbc = EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymB)+iC-1)-EVir(mAdVir(iSymA)+iB-1)
          Density%SB(iSymA)%A2(iBB,iAA) = Density%SB(iSymA)%A2(iBB,iAA)+Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
          WDensity%SB(iSymA)%A2(iAA,iBB) = WDensity%SB(iSymA)%A2(iAA,iBB)-Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
        end do
        !write(u6,*) 'IJABC',iI,iJ,iA,iB,iC
        !write(u6,*) 'Symm',iSymI,iSymJ,iSymA,iSymB,iSymC
        !write(u6,*) 'EDiffBC',EDiffbc
        !write(u6,*) 'EDiffAC',EDiffac
        !write(u6,*) 'xicja',xicja
        !write(u6,*) 'xicjb',xicjb
        !write(u6,*) 'xiajc',xiajc
        !write(u6,*) 'xibjc',xibjc
        do iBDel=1,nDel(iSymA)
          iBB = iBDel+nFro(iSymA)+nOcc(iSymA)+nExt(iSymA)
          xibjc = Int1(iBC+iBDel+nExt(iSymA)+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (ISymI == iSymJ) then
            xicjb = Int1(iBC+iC+(iBDel+nExt(iSymA)-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xicjb = Int2(iCB+iC+(iBDel+nExt(iSymA)-1)*(nOrb(iSymB)+nDel(iSymB)))

          end if
          EDiffbc = EVir(mAdVir(iSymA)+iA-1)-EVir(mAdDel(iSymA)+iBDel-1)
          Density%SB(iSymA)%A2(iAA,iBB) = Density%SB(iSymA)%A2(iAA,iBB)+Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
          WDensity%SB(iSymA)%A2(iAA,iBB) = WDensity%SB(iSymA)%A2(iAA,iBB)-Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
          !-------------------------------------------------------------
          !write(u6,*) 'xicjb',xicjb
          !write(u6,*) 'xibjc',xibjc
          !write(u6,*) 'EDiffbc',EDiffbc
          !write(u6,*) 'E_a',EVir(mAdVir(iSymA)+iA-1)
          !write(u6,*) 'E_B',EVir(mAdDel(iSymA)+iBDel-1)
          !write(u6,*) 'Tij',T_ij
          !write(u6,*) 'Tji',T_ji
          !write(u6,*) 'iB, iC',iBDel+nExt(iSymA),iC
          !write(u6,*) 'iI, iJ',iI,iJ
          !-------------------------------------------------------------
        end do
        !---------------------------------------------------------------
        !write(u6,*) 'The Denom is',EDenom
        !write(u6,*) 'xicja',xicja
        !write(u6,*) 'xicjb',xicjb
        !write(u6,*) 'xiajc',xiajc
        !write(u6,*) 'xibjc',xibjc
        !write(u6,*) 'this is with the unmirrored'
        !---------------------------------------------------------------

      end do !iOrbB
    end do   !iOrbA

    ! Here is the extra loop.
    if (iSymA /= iSymB) then
      iAC = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
      iCA = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
      iBC = iAC
      iCB = iCA
      do iC=1,nExt(iSymA)
        do iA=1,nExt(iSymB)
          iAA = iA+nFro(iSymB)+nOcc(iSymB)
          EDiffac = EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymB)+iA-1)-EVir(mAdVir(iSymA)+iC-1)
          xiajc = Int1(iCA+iC+(iA-1)*(nOrb(iSymA)+nDel(iSymA)))
          xicja = Int2(iAC+iA+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
          T_ij = (Two*xiajc-xicja)/EDiffac
          T_ji = (Two*xicja-xiajc)/EDiffac
          do iB=1,nExt(iSymB)
            iBB = iB+nFro(iSymB)+nOcc(iSymB)
            EDiffbc = EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymA)+iC-1)-EVir(mAdVir(iSymB)+iB-1)
            xibjc = Int1(iCB+iC+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
            xicjb = Int2(iBC+iB+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            Density%SB(iSymB)%A2(iBB,iAA) = Density%SB(iSymB)%A2(iBB,iAA)+Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
            WDensity%SB(iSymB)%A2(iAA,iBB) = WDensity%SB(iSymB)%A2(iAA,iBB)-Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))

          end do
          do iBDel=1,nDel(iSymB)
            iBB = iBDel+nFro(iSymB)+nOcc(iSymB)+nExt(iSymB)
            xibjc = Int1(iBC+iBDel+nExt(iSymB)+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            xicjb = Int2(iCB+iC+(iBDel+nExt(iSymB)-1)*(nOrb(iSymA)+nDel(iSymA)))
            EDiffbc = EVir(mAdVir(iSymB)+iA-1)-EVir(mAdDel(iSymB)+iBDel-1)
            Density%SB(iSymB)%A2(iAA,iBB) = Density%SB(iSymB)%A2(iAA,iBB)+Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
            WDensity%SB(iSymB)%A2(iAA,iBB) = WDensity%SB(iSymB)%A2(iAA,iBB)-Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
          end do
          !-------------------------------------------------------------
          !write(u6,*) 'The Denom is',EDenom
          !write(u6,*) 'xicja',xicja
          !write(u6,*) 'xicjb',xicjb
          !write(u6,*) 'xiajc',xiajc
          !write(u6,*) 'xibjc',xibjc
          !write(u6,*) 'this is with the mirrored symmetry'
          !-------------------------------------------------------------

        end do !iOrbB
      end do   !iOrbA
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Lagrangian term 2
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    iAB = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
    iBA = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
    iKB = (nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
    iBK = nFro(iSymB)+nOcc(iSymB)
    do iA=1,nExt(iSymA)
      iAA = iA+nFro(iSymA)+nOcc(iSymA)
      do iB=1,nExt(iSymB)
        EDenom = (EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1))
        xiajb = Int1(iAB+iA+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (isymI == iSymJ) then
          xibja = Int1(iBA+iB+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xibja = Int2(iBA+iB+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        end if
        Tij = (Two*xiajb-xibja)/EDenom
        Tji = (Two*xibja-xiajb)/EDenom
        do iK=1,nFro(iSymA)+nOcc(iSymA)
          xikjb = Int1(iKB+iK+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (isymI == iSymJ) then
            xibjk = Int1(iBK+iB+(iK-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xibjk = Int2(iBK+iB+(iK-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          Mp2Lagr%SB(iSymA)%A2(iK,iA) = Mp2Lagr%SB(iSymA)%A2(iK,iA)+Fac_ij*(Tij*Xikjb+Tji*Xibjk)
          WDensity%SB(iSymA)%A2(iAA,iK) = WDensity%SB(iSymA)%A2(iAA,iK)-Two*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
          !-------------------------------------------------------------
          !write(u6,*) 'The Denom is',EDenom
          !write(u6,*) 'Tij  ',Tij
          !write(u6,*) 'Tji  ',Tji
          !write(u6,*) 'xiajb',xiajb
          !write(u6,*) 'xiajb',xiajb
          !write(u6,*) 'xibja',xibja
          !write(u6,*) 'xikjb',xikjb
          !write(u6,*) 'xibjk',xibjk
          !-------------------------------------------------------------
        end do !OrbK
      end do   !OrbA
    end do     !OrbB

    ! In the case that the symmetries A and B are not the same
    ! only the case where A is greater than B is written on disk
    ! so we must do an extra loop for where B is greater than A

    if (iSymA /= iSymB) then
      iAB = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
      iBA = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymB))
      iKB = (nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))
      iBK = nFro(iSymA)+nOcc(iSymA)
      do iA=1,nExt(iSymB)
        iAA = iA+nFro(iSymB)+nOcc(iSymB)
        do iB=1,nExt(iSymA)
          EDenom = (EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymB)+iA-1)-EVir(mAdVir(iSymA)+iB-1))
          xiajb = Int1(iBA+iB+(iA-1)*(nOrb(iSymA)+nDel(iSymA)))
          xibja = Int2(iAB+iA+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          Tij = (Two*xiajb-xibja)/EDenom
          Tji = (Two*xibja-xiajb)/EDenom
          do iK=1,nFro(iSymB)+nOcc(iSymB)
            xikjb = Int1(iBK+iB+(iK-1)*(nOrb(iSymA)+nDel(iSymA)))

            xibjk = Int2(iKB+iK+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))

            Mp2Lagr%SB(iSymB)%A2(iK,iA) = Mp2Lagr%SB(iSymB)%A2(iK,iA)+Fac_ij*(Tij*Xikjb+Tji*Xibjk)
            WDensity%SB(iSymB)%A2(iAA,iK) = WDensity%SB(iSymB)%A2(iAA,iK)-Two*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
            !-----------------------------------------------------------
            !write(u6,*) 'The Denom is',EDenom
            !write(u6,*) 'Tij  ',Tij
            !write(u6,*) 'Tji  ',Tji
            !write(u6,*) 'xijab',xijab
            !write(u6,*) 'xijba',xijba
            !write(u6,*) 'xijkb',xijkb
            !write(u6,*) 'xijbk',xijbk
            !write(u6,*) 'this is with the mirrored symmetry'
            !------------------------------------------------------------
          end do !OrbK
        end do   !OrbA
      end do     !OrbB
    end if

    ! When we are done with a block we free the memory again.

  end do         !OrbI
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Terms that need <A.|B.>-integrals
!                                                                      *
!***********************************************************************
!                                                                      *
nB = nExt(iSymB)
do iA=1,nExt(iSymA)
  ! Decide i the AB-matrix is symmetric, if so
  ! only loop over lower triangle
  if (iSymA == iSymB) nB = iA

  do iB=1,nB
    !write(u6,*) 'iSymA, iSymB',iSymA,iSymB
    !write(u6,*) 'iA, iB',iA,iB
    fac_ab = One
    if ((iA /= iB) .and. (iSymA == iSymB)) fac_ab = Two
    ! If I is not equal J there is a symmetry J I B A
    ! that is identical and not stored.
    if (iSymA /= iSymB) fac_ab = Two

    call Exch(iSymI,iSymA,iSymJ,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Int1,Scr1)
    if (iSymA /= iSymB) then
      call Exch(iSymJ,iSymA,iSymI,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Int2,Scr1)
    end if
    !write(u6,*)
    !write(u6,*) ' *  A,B = ',iA,iB
    !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymI)+nDel(iSymI),nOrb(iSymJ)+nDel(iSymJ))
    !if (iSymA /= iSymB) then
    !   call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSymJ) + nDel(iSymJ),nOrb(iSymI) + nDel(iSymI))
    !end if
    !write(u6,*)

    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Calculate Pij
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! We start by taking i and j in symI and k in symJ
    iIK = nFro(iSymI)+nFro(iSymJ)*(nOrb(iSymI)+nDel(iSymI))
    iKI = nFro(iSymJ)+nFro(iSymI)*(nOrb(iSymJ)+nDel(iSymJ))
    iJK = iIK
    iKJ = iKI
    do iK=1,nOcc(iSymJ)
      do iI=1,nOcc(iSymI)
        iII = iI+nFro(iSymI)
        EDiffik = EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iK-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1)
        xiakb = Int1(iIK+iI+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
        if (iSymA == iSymB) then
          xkaib = Int1(iKI+iK+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
        else
          xkaib = Int2(iKI+iK+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
        end if
        T_ab = (Two*xiakb-xkaib)/EDiffik
        T_ba = (Two*xkaib-xiakb)/EDiffik
        do iJ=1,nOcc(iSymI)
          iJJ = iJ+nFro(iSymI)
          ! Calculating the denominator

          EDiffjk = EOcc(mAdOcc(iSymI)+iJ-1)+EOcc(mAdOcc(iSymJ)+iK-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1)

          xjakb = Int1(iJK+iJ+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymA == iSymB) then
            xkajb = Int1(iKJ+iK+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xkajb = Int2(iKJ+iK+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          !-------------------------------------------------------------
          !write(u6,*) 'The Denom is',EDenom
          !write(u6,*) 'xikab',xikab
          !write(u6,*) 'xjkab',xjkab
          !write(u6,*) 'xkiab',xkiab
          !write(u6,*) 'xkjab',xkjab
          !-------------------------------------------------------------
          Density%SB(iSymI)%A2(iJJ,iII) = Density%SB(iSymI)%A2(iJJ,iII)-Fac_ab*(xjakb*T_ab+xkajb*T_ba)/EDiffjk
          WDensity%SB(iSymI)%A2(iJJ,iII) = WDensity%SB(iSymI)%A2(iJJ,iII)-Fac_ab*(xjakb*T_ab+xkajb*T_ba)

        end do
        do iJFroz=1,nFro(iSymI)
          xjakb = Int1(iJFroz+(iK+nFro(iSymJ)-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymI == iSymJ) then
            xkajb = Int1(iK+nFro(iSymJ)+(iJFroz-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xkajb = Int2(iK+nFro(iSymJ)+(iJFroz-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          EDiffjk = EOcc(mAdOcc(iSymI)+iI-1)-EOcc(mAdFro(iSymI)+iJFroz-1)
          Density%SB(iSymI)%A2(iJFroz,iII) = Density%SB(iSymI)%A2(iJFroz,iII)+Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))/EDiffjk
          WDensity%SB(iSymI)%A2(iJFroz,iII) = WDensity%SB(iSymI)%A2(iJFroz,iII)-Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))
        end do
      end do
    end do
    ! And then we do the same thing with i and j in symJ and k in symI
    if (iSymI /= iSymJ) then
      iIK = nFro(iSymJ)+nFro(iSymI)*(nOrb(iSymJ)+nDel(iSymJ))
      iKI = nFro(iSymI)+nFro(iSymJ)*(nOrb(iSymI)+nDel(iSymI))
      iJK = iIK
      iKJ = iKI
      do iK=1,nOcc(iSymI)
        do iI=1,nOcc(iSymJ)
          iII = iI+nFro(iSymJ)
          EDiffik = EOcc(mAdOcc(iSymJ)+iI-1)+EOcc(mAdOcc(iSymI)+iK-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1)
          xiakb = Int1(iKI+iK+(iI-1)*(nOrb(iSymI)+nDel(iSymI)))
          xkaib = Int2(iIK+iI+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          T_ab = (Two*xiakb-xkaib)/EDiffik
          T_ba = (Two*xkaib-xiakb)/EDiffik
          do iJ=1,nOcc(iSymJ)
            iJJ = iJ+nFro(iSymJ)
            ! Calculating the denominator

            EDiffjk = EOcc(mAdOcc(iSymJ)+iJ-1)+EOcc(mAdOcc(iSymI)+iK-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1)

            xjakb = Int1(iKJ+iK+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
            xkajb = Int2(iJK+iJ+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))

            Density%SB(iSymJ)%A2(iJJ,iII) = Density%SB(iSymJ)%A2(iJJ,iII)-Fac_ab*(xjakb*T_ab+xkajb*T_ba)/EDiffjk
            WDensity%SB(iSymJ)%A2(iJJ,iII) = WDensity%SB(iSymJ)%A2(iJJ,iII)-Fac_ab*(xjakb*T_ab+xkajb*T_ba)

            !-----------------------------------------------------------
            !write(u6,*) 'The Denom is',EDenom
            !write(u6,*) 'xikab',xikab
            !write(u6,*) 'xjkab',xjkab
            !write(u6,*) 'xkiab',xkiab
            !write(u6,*) 'xkjab',xkjab
            !-----------------------------------------------------------
          end do
          do iJFroz=1,nFro(iSymJ)
            xjakb = Int1(iK+nFro(iSymI)+(iJFroz-1)*(nOrb(iSymI)+nDel(iSymI)))
            xkajb = Int2(iJFroz+(iK+nFro(iSymI)-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            EDiffjk = EOcc(mAdOcc(iSymJ)+iI-1)-EOcc(mAdFro(iSymJ)+iJFroz-1)
            Density%SB(iSymJ)%A2(iJFroz,iII) = Density%SB(iSymJ)%A2(iJFroz,iII)+Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))/EDiffjk
            WDensity%SB(iSymJ)%A2(iJFroz,iII) = WDensity%SB(iSymJ)%A2(iJFroz,iII)-Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))
          end do
        end do
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Lagrangian term 1
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    iCJ = nFro(iSymI)+nOcc(iSymI)+nFro(iSymJ)*(nOrb(iSymI)+nDel(iSymI))
    iJC = nFro(iSymJ)+(nFro(iSymI)+nOcc(iSymI))*(nOrb(iSymJ)+nDel(iSymJ))
    iIJ = nFro(iSymI)+nFro(iSymJ)*(nOrb(iSymI)+nDel(iSymI))
    IJI = nFro(iSymJ)+nFro(iSymI)*(nOrb(iSymJ)+nDel(iSymJ))
    do iC=1,nExt(iSymI)+nDel(iSymI)
      do iI=1,nOcc(iSymI)
        iII = iI+nFro(iSymI)
        do iJ=1,nOcc(iSymJ)
          EDenom = One/(EOcc(mAdOcc(iSymI)+iI-1)+EOcc(mAdOcc(iSymJ)+iJ-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1))
          xaibj = Int1(iIJ+iI+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
          xacbj = Int1(iCJ+iC+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymA == iSymB) then
            xbiaj = Int1(iJI+iJ+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xbcaj = Int1(iJC+iJ+(iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xbiaj = Int2(iJI+iJ+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xbcaj = Int2(iJC+iJ+(iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          Mp2Lagr%SB(iSymI)%A2(iII,iC) = Mp2Lagr%SB(iSymI)%A2(iII,iC)+ &
                                         EDenom*Fac_ab*(-Two*(xaibj*xacbj+xbiaj*xbcaj)+xaibj*xbcaj+xbiaj*xacbj)

          !-------------------------------------------------------------
          !write(u6,*) 'EDenom',EDenom
          !write(u6,*) 'Fac_ab',Fac_ab
          !write(u6,*) 'xaibj',xaibj
          !write(u6,*) 'xacbj',xacbj
          !write(u6,*) 'xbiaj',xbiaj
          !write(u6,*) 'xbcaj',xbcaj
          !-------------------------------------------------------------
        end do
      end do
    end do
    ! And now we switch symmetries and go again
    if (iSymA /= iSymB) then
      iCJ = nFro(iSymJ)+nOcc(iSymJ)+nFro(iSymI)*(nOrb(iSymJ)+nDel(iSymJ))
      iJC = nFro(iSymI)+(nFro(iSymJ)+nOcc(iSymJ))*(nOrb(iSymI)+nDel(iSymI))
      iIJ = nFro(iSymJ)+nFro(iSymI)*(nOrb(iSymJ)+nDel(iSymJ))
      IJI = nFro(iSymI)+nFro(iSymJ)*(nOrb(iSymI)+nDel(iSymI))
      do iC=1,nExt(iSymJ)+nDel(iSymJ)
        do iI=1,nOcc(iSymJ)
          iII = iI+nFro(iSymJ)
          do iJ=1,nOcc(iSymI)
            EDenom = One/(EOcc(mAdOcc(iSymJ)+iI-1)+EOcc(mAdOcc(iSymI)+iJ-1)-EVir(mAdVir(iSymA)+iA-1)-EVir(mAdVir(iSymB)+iB-1))
            xaibj = Int1(iJI+iJ+(iI-1)*(nOrb(iSymI)+nDel(iSymI)))
            xacbj = Int1(iJC+iJ+(iC-1)*(nOrb(iSymI)+nDel(iSymI)))
            xajbi = Int2(iIJ+iI+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xajbc = Int2(iCJ+iC+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            Mp2Lagr%SB(iSymJ)%A2(iII,iC) = Mp2Lagr%SB(iSymJ)%A2(iII,iC)+ &
                                           EDenom*Fac_ab*(-Two*(xaibj*xacbj+xajbi*xajbc)+xaibj*xajbc+xajbi*xacbj)

            !-----------------------------------------------------------
            !write(u6,*) 'EDenom',EDenom
            !write(u6,*) 'Fac_ab',Fac_ab
            !write(u6,*) 'xijab',xijab
            !write(u6,*) 'xcjab',xcjab
            !write(u6,*) 'xjiab',xjiab
            !write(u6,*) 'xjcab',xjcab
            !-----------------------------------------------------------
          end do
        end do
      end do
    end if

  end do
end do

return

end subroutine rhs_mp2_help1
