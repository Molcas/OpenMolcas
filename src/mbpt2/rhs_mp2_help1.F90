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

subroutine rhs_mp2_help1(iSymA,iSymB,iSymI,iSymJ)

implicit real*8(a-h,o-z)
!defining One etc.
#include "real.fh"
#include "trafo.fh"
#include "WrkSpc.fh"
#include "files_mbpt2.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
!                                                                      *
!***********************************************************************
!                                                                      *
! Statement function
iVirVir(i,j,k) = nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iOccOcc(i,j,k) = nFro(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iVirDel(i,j,k) = nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)+nExt(k)-1)
iOccFro(i,j,k) = j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iVirOcc(i,j,k) = j-1+nFro(k)+nOcc(k)+(nOrb(k)+nDel(k))*(i-1)
iMp2Lagr(i,j,k) = ip_Mp2Lagr(k)+j-1+(nOcc(k)+nFro(k))*(i-1)

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
    !write(6,*) 'iSymI, iSymJ',iSymI,iSymJ
    !write(6,*) 'iI, iJ',iI,iJ
    ! If only blocks with OrbJ < OrbI is used we need
    ! to multiply the result with two for offdiag terms.
    fac_ij = One
    if ((iI /= iJ) .and. (iSymI == iSymJ)) fac_ij = Two
    ! If I is not equal J there is a symmetry J I B A
    ! that is identical and not stored.
    if (iSymI /= iSymJ) fac_ij = Two
    ! Needed?

    ! Copy one AB-block of integrals to the workspace
    call Exch(iSymA,iSymI,iSymB,iSymJ,iI+nFro(iSymI),iJ+nFro(iSymJ),Work(ipInt1),Work(ipScr1))
    if (iSymI /= iSymJ) then
      call Exch(iSymB,iSymI,iSymA,iSymJ,iI+nFro(iSymI),iJ+nFro(iSymJ),Work(ipInt2),Work(ipScr1))
    end if
    !write(6,*) 'Rhs_mp2_help1'
    !write(6,*) ' *  i,j = ',iI,iJ
    !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymA)+nDel(iSymA),nOrb(iSymB)+nDel(iSymB))
    !if (iSymI /= iSymJ) then
    !  call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymB)+nDel(iSymB),nOrb(iSymA)+nDel(iSymA))
    !end if
    !write(6,*) ''

    ! Calculates the EMP2-energy as a test
    iAB = nFro(iSymA)+nOcc(iSymA)+(nFro(iSymB)+nOcc(iSymB))*(nOrb(iSymA)+nDel(iSymA))
    iBA = nFro(iSymB)+nOcc(iSymB)+(nFro(iSymA)+nOcc(iSymA))*(nOrb(iSymB)+nDel(iSymB))

    do iA=1,nExt(iSymA)
      nB = nExt(iSymB)
      if (iSymA == iSymB) nB = iA
      do iB=1,nB
        EDiff = One/(work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1))

        xiajb = work(ipInt1+iAB+iA-1+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (iSymI == iSymJ) then
          xibja = work(ipInt1+iAB+iB-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xibja = work(ipInt2+iBA+iB-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
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
        xiajc = work(ipInt1+iAC+iA-1+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (ISymI == iSymJ) then
          xicja = work(ipInt1+iCA+iC-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xicja = work(ipInt2+iCA+iC-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        end if
        EDiffac = work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iC-1)
        T_ij = (Two*xiajc-xicja)/EDiffac
        T_ji = (Two*xicja-xiajc)/EDiffac
        do iB=1,nExt(iSymA)
          xibjc = work(ipInt1+iBC+iB-1+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (ISymI == iSymJ) then
            xicjb = work(ipInt1+iCB+iC-1+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xicjb = work(ipInt2+iCB+iC-1+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          EDiffbc = work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymB)+iC-1)-work(mAdVir(iSymA)+iB-1)
          Work(ip_Density(iSymA)+iVirVir(iA,iB,iSymA)) = Work(ip_Density(iSymA)+iVirVir(iA,iB,iSymA))+ &
                                                         Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
          Work(ip_WDensity(iSymA)+iVirVir(iB,iA,iSymA)) = Work(ip_WDensity(iSymA)+iVirVir(iB,iA,iSymA))- &
                                                          Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
        end do
        !write(6,*) 'IJABC',iI,iJ,iA,iB,iC
        !write(6,*) 'Symm',iSymI,iSymJ,iSymA,iSymB,iSymC
        !write(6,*) 'EDiffBC',EDiffbc
        !write(6,*) 'EDiffAC',EDiffac
        !write(6,*) 'xicja',xicja
        !write(6,*) 'xicjb',xicjb
        !write(6,*) 'xiajc',xiajc
        !write(6,*) 'xibjc',xibjc
        do iBDel=1,nDel(iSymA)
          xibjc = work(ipInt1+iBC+iBDel-1+nExt(iSymA)+(iC-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (ISymI == iSymJ) then
            xicjb = work(ipInt1+iBC+iC-1+(iBDel+nExt(iSymA)-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xicjb = work(ipInt2+iCB+iC-1+(iBDel+nExt(iSymA)-1)*(nOrb(iSymB)+nDel(iSymB)))

          end if
          EDiffbc = Work(madVir(iSymA)+iA-1)-Work(mAdDel(iSymA)+iBDel-1)
          Work(ip_Density(iSymA)+iVirDel(iBDel,iA,iSymA)) = Work(ip_Density(iSymA)+iVirDel(iBDel,iA,iSymA))+ &
                                                            Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
          Work(ip_WDensity(iSymA)+iVirDel(iBDel,iA,iSymA)) = Work(ip_WDensity(iSymA)+iVirDel(iBDel,iA,iSymA))- &
                                                             Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
          !-------------------------------------------------------------
          !write(6,*) 'index Dens',iVirDel(iBDel,iA,iSymA)
          !write(6,*) 'xicjb',xicjb
          !write(6,*) 'xibjc',xibjc
          !write(6,*) 'EDiffbc',EDiffbc
          !write(6,*) 'E_a',Work(madVir(iSymA)+iA-1)
          !write(6,*) 'E_B',Work(mAdDel(iSymA)+iBDel-1)
          !write(6,*) 'Tij',T_ij
          !write(6,*) 'Tji',T_ji
          !write(6,*) 'iB, iC',iBDel+nExt(iSymA),iC
          !write(6,*) 'iI, iJ',iI,iJ
          !-------------------------------------------------------------
        end do
        !Work(indexW) = Work(indexW)-Fac_ij*EDiffbc*(Two*(xijac*xijbc+xijca*xijcb)-xijac*xijcb-xijca*xijbc)
        !---------------------------------------------------------------
        !write(6,*) 'The Denom is',EDenom
        !write(6,*) 'xicja',xicja
        !write(6,*) 'xicjb',xicjb
        !write(6,*) 'xiajc',xiajc
        !write(6,*) 'xibjc',xibjc
        !write(6,*) 'this is with the unmirrored'
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
          EDiffac = work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymB)+iA-1)-work(mAdVir(iSymA)+iC-1)
          xiajc = work(ipInt1+iCA+iC-1+(iA-1)*(nOrb(iSymA)+nDel(iSymA)))
          xicja = work(ipInt2+iAC+iA-1+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
          T_ij = (Two*xiajc-xicja)/EDiffac
          T_ji = (Two*xicja-xiajc)/EDiffac
          do iB=1,nExt(iSymB)
            EDiffbc = work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymA)+iC-1)-work(mAdVir(iSymB)+iB-1)
            xibjc = work(ipInt1+iCB+iC-1+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
            xicjb = work(ipInt2+iBC+iB-1+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            Work(ip_Density(iSymB)+iVirVir(iA,iB,iSymB)) = Work(ip_Density(iSymB)+iVirVir(iA,iB,iSymB))+ &
                                                           Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
            Work(ip_WDensity(iSymB)+iVirVir(iB,iA,iSymB)) = Work(ip_WDensity(iSymB)+iVirVir(iB,iA,iSymB))- &
                                                            Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))

          end do
          do iBDel=1,nDel(iSymB)
            xibjc = work(ipInt1+iBC+iBDel-1+nExt(iSymB)+(iC-1)*(nOrb(iSymB)+nDel(iSymB)))
            xicjb = work(ipInt2+iCB+iC-1+(iBdel+nExt(iSymB)-1)*(nOrb(iSymA)+nDel(iSymA)))
            EDiffbc = Work(madVir(iSymB)+iA-1)-Work(mAdDel(iSymB)+iBDel-1)
            Work(ip_Density(iSymA)+iVirDel(iBDel,iA,iSymB)) = Work(ip_Density(iSymB)+iVirDel(iBDel,iA,iSymB))+ &
                                                              Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))/EDiffbc
            Work(ip_WDensity(iSymA)+iVirDel(iBDel,iA,iSymB)) = Work(ip_WDensity(iSymB)+iVirDel(iBDel,iA,iSymB))- &
                                                               Fac_ij*((xibjc*T_ij)+(xicjb*T_ji))
          end do
          !-------------------------------------------------------------
          !write(6,*) 'The Denom is',EDenom
          !write(6,*) 'xicja',xicja
          !write(6,*) 'xicjb',xicjb
          !write(6,*) 'xiajc',xiajc
          !write(6,*) 'xibjc',xibjc
          !write(6,*) 'this is with the mirrored symmetry'
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
      do iB=1,nExt(iSymB)
        EDenom = (work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1))
        xiajb = work(ipInt1+iAB+iA-1+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
        if (isymI == iSymJ) then
          xibja = work(ipInt1+iBA+iB-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        else
          xibja = work(ipInt2+iBA+iB-1+(iA-1)*(nOrb(iSymB)+nDel(iSymB)))
        end if
        Tij = (Two*xiajb-xibja)/EDenom
        Tji = (Two*xibja-xiajb)/EDenom
        do iK=1,nFro(iSymA)+nOcc(iSymA)
          xikjb = work(ipInt1+iKB+iK-1+(iB-1)*(nOrb(iSymA)+nDel(iSymA)))
          if (isymI == iSymJ) then
            xibjk = work(ipInt1+iBK+iB-1+(iK-1)*(nOrb(iSymB)+nDel(iSymB)))
          else
            xibjk = work(ipInt2+iBK+iB-1+(iK-1)*(nOrb(iSymB)+nDel(iSymB)))
          end if
          Work(iMp2Lagr(iA,iK,iSymA)) = Work(iMp2Lagr(iA,iK,iSymA))+Fac_ij*(Tij*Xikjb+Tji*Xibjk)
          Work(ip_WDensity(iSymA)+iVirOcc(iK,iA,iSymA)) = Work(ip_WDensity(iSymA)+iVirOcc(iK,iA,iSymA))- &
                                                          2.0d0*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
          !-------------------------------------------------------------
          !write(6,*) 'The Denom is',EDenom
          !write(6,*) 'Tij  ',Tij
          !write(6,*) 'Tji  ',Tji
          !write(6,*) 'xiajb',xiajb
          !write(6,*) 'xiajb',xiajb
          !write(6,*) 'xibja',xibja
          !write(6,*) 'xikjb',xikjb
          !write(6,*) 'xibjk',xibjk
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
        do iB=1,nExt(iSymA)
          EDenom = (work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymB)+iA-1)-work(mAdVir(iSymA)+iB-1))
          xiajb = work(ipInt1+iBA+iB-1+(iA-1)*(nOrb(iSymA)+nDel(iSymA)))
          xibja = work(ipInt2+iAB+iA-1+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))
          Tij = (Two*xiajb-xibja)/EDenom
          Tji = (Two*xibja-xiajb)/EDenom
          do iK=1,nFro(iSymB)+nOcc(iSymB)
            xikjb = work(ipInt1+iBK+iB-1+(iK-1)*(nOrb(iSymA)+nDel(iSymA)))

            xibjk = work(ipInt2+iKB+iK-1+(iB-1)*(nOrb(iSymB)+nDel(iSymB)))

            Work(iMp2Lagr(iA,iK,iSymB)) = Work(iMp2Lagr(iA,iK,iSymB))+Fac_ij*(Tij*Xikjb+Tji*Xibjk)
            Work(ip_WDensity(iSymB)+iVirOcc(iK,iA,iSymB)) = Work(ip_WDensity(iSymB)+iVirOcc(iK,iA,iSymB))- &
                                                            2.0d0*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
            !-----------------------------------------------------------
            !write(6,*) 'The Denom is',EDenom
            !write(6,*) 'Tij  ',Tij
            !write(6,*) 'Tji  ',Tji
            !write(6,*) 'xijab',xijab
            !write(6,*) 'xijba',xijba
            !write(6,*) 'xijkb',xijkb
            !write(6,*) 'xijbk',xijbk
            !write(6,*) 'this is with the mirrored symmetry'
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
    !write(6,*) 'iSymA, iSymB',iSymA,iSymB
    !write(6,*) 'iA, iB',iA,iB
    fac_ab = One
    if ((iA /= iB) .and. (iSymA == iSymB)) fac_ab = Two
    ! If I is not equal J there is a symmetry J I B A
    ! that is identical and not stored.
    if (iSymA /= iSymB) fac_ab = Two

    call Exch(iSymI,iSymA,iSymJ,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Work(ipInt1),Work(ipScr1))
    if (iSymA /= iSymB) then
      call Exch(iSymJ,iSymA,iSymI,iSymB,iA+nOcc(iSymA)+nFro(iSymA),iB+nOcc(iSymB)+nFro(iSymB),Work(ipInt2),Work(ipScr1))
    end if
    !write(6,*)
    !write(6,*) ' *  A,B = ',iA,iB
    !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymI)+nDel(iSymI),nOrb(iSymJ)+nDel(iSymJ))
    !if (iSymA /= iSymB) then
    !   call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymJ) + nDel(iSymJ),nOrb(iSymI) + nDel(iSymI))
    !end if
    !write(6,*)

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
        EDiffik = work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iK-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1)
        xiakb = Work(ipInt1+iIK+iI-1+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
        if (iSymA == iSymB) then
          xkaib = Work(ipInt1+iKI+iK-1+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
        else
          xkaib = Work(ipInt2+iKI+iK-1+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
        end if
        T_ab = (Two*xiakb-xkaib)/EDiffik
        T_ba = (Two*xkaib-xiakb)/EDiffik
        do iJ=1,nOcc(iSymI)
          ! Calculating the denominator

          EDiffjk = work(mAdOcc(iSymI)+iJ-1)+work(mAdOcc(iSymJ)+iK-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1)

          xjakb = Work(ipInt1+iJK+iJ-1+(iK-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymA == iSymB) then
            xkajb = Work(ipInt1+iKJ+iK-1+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xkajb = Work(ipInt2+iKJ+iK-1+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          !-------------------------------------------------------------
          !write(6,*) 'The Denom is',EDenom
          !write(6,*) 'xikab',xikab
          !write(6,*) 'xjkab',xjkab
          !write(6,*) 'xkiab',xkiab
          !write(6,*) 'xkjab',xkjab
          !-------------------------------------------------------------
          Work(ip_Density(iSymI)+iOccOcc(iI,iJ,iSymI)) = Work(ip_Density(iSymI)+iOccOcc(iI,iJ,iSymI))- &
                                                         Fac_ab*(xjakb*T_ab+xkajb*T_ba)/EDiffjk
          Work(ip_WDensity(iSymI)+iOccOcc(iI,iJ,iSymI)) = Work(ip_WDensity(iSymI)+iOccOcc(iI,iJ,iSymI))- &
                                                          Fac_ab*(xjakb*T_ab+xkajb*T_ba)

        end do
        do iJFroz=1,nFro(iSymI)
          xjakb = Work(ipInt1+iJFroz-1+(iK+nFro(iSymJ)-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymI == iSymJ) then
            xkajb = Work(ipInt1+iK-1+nFro(iSymJ)+(iJFroz-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xkajb = Work(ipInt2+iK-1+nFro(iSymJ)+(iJFroz-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          EDiffjk = Work(madOcc(iSymI)+iI-1)-Work(mAdFro(iSymI)+iJFroz-1)
          Work(ip_Density(iSymI)+iOccFro(iI,iJFroz,iSymI)) = Work(ip_Density(iSymI)+iOccFro(iI,iJFroz,iSymI))+ &
                                                             Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))/EDiffjk
          Work(ip_WDensity(iSymI)+iOccFro(iI,iJFroz,iSymI)) = Work(ip_WDensity(iSymI)+iOccFro(iI,iJFroz,iSymI))- &
                                                              Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))
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
          EDiffik = work(mAdOcc(iSymJ)+iI-1)+work(mAdOcc(iSymI)+iK-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1)
          xiakb = Work(ipInt1+iKI+iK-1+(iI-1)*(nOrb(iSymI)+nDel(iSymI)))
          xkaib = Work(ipInt2+iIK+iI-1+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          T_ab = (Two*xiakb-xkaib)/EDiffik
          T_ba = (Two*xkaib-xiakb)/EDiffik
          do iJ=1,nOcc(iSymJ)
            ! Calculating the denominator

            EDiffjk = work(mAdOcc(iSymJ)+iJ-1)+work(mAdOcc(iSymI)+iK-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1)

            xjakb = Work(ipInt1+iKJ+iK-1+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
            xkajb = Work(ipInt2+iJK+iJ-1+(iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))

            Work(ip_Density(iSymJ)+iOccOcc(iI,iJ,iSymJ)) = Work(ip_Density(iSymJ)+iOccOcc(iI,iJ,iSymJ))- &
                                                           Fac_ab*(xjakb*T_ab+xkajb*T_ba)/EDiffjk
            Work(ip_WDensity(iSymJ)+iOccOcc(iI,iJ,iSymJ)) = Work(ip_WDensity(iSymJ)+iOccOcc(iI,iJ,iSymJ))- &
                                                            Fac_ab*(xjakb*T_ab+xkajb*T_ba)

            !-----------------------------------------------------------
            !write(6,*) 'The Denom is',EDenom
            !write(6,*) 'xikab',xikab
            !write(6,*) 'xjkab',xjkab
            !write(6,*) 'xkiab',xkiab
            !write(6,*) 'xkjab',xkjab
            !-----------------------------------------------------------
          end do
          do iJFroz=1,nFro(iSymJ)
            xjakb = Work(ipInt1+iK+nFro(iSymI)-1+(iJFroz-1)*(nOrb(iSymI)+nDel(iSymI)))
            xkajb = Work(ipInt2+iJFroz-1+(iK+nFro(iSymI)-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            EDiffjk = Work(madOcc(iSymJ)+iI-1)-Work(mAdFro(iSymJ)+iJFroz-1)
            Work(ip_Density(iSymJ)+iOccFro(iI,iJFroz,iSymJ)) = Work(ip_Density(iSymJ)+iOccFro(iI,iJFroz,iSymJ))+ &
                                                               Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))/EDiffjk
            Work(ip_WDensity(iSymJ)+iOccFro(iI,iJFroz,iSymJ)) = Work(ip_WDensity(iSymJ)+iOccFro(iI,iJFroz,iSymJ))- &
                                                                Fac_ab*((xjakb*T_ab)+(xkajb*T_ba))
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
        do iJ=1,nOcc(iSymJ)
          EDenom = One/(work(mAdOcc(iSymI)+iI-1)+work(mAdOcc(iSymJ)+iJ-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1))
          xaibj = Work(ipInt1+iIJ+iI-1+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
          xacbj = Work(ipInt1+iCJ+iC-1+(iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
          if (iSymA == iSymB) then
            xbiaj = Work(ipInt1+iJI+iJ-1+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xbcaj = Work(ipInt1+iJC+iJ-1+(iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          else
            xbiaj = Work(ipInt2+iJI+iJ-1+(iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xbcaj = Work(ipInt2+iJC+iJ-1+(iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
          end if
          Work(iMp2Lagr(iC,iI+nFro(iSymI),iSymI)) = Work(iMp2Lagr(iC,iI+nFro(iSymI),iSymI))+ &
                                                    EDenom*Fac_ab*(-Two*(xaibj*xacbj+xbiaj*xbcaj)+xaibj*xbcaj+xbiaj*xacbj)

          !-------------------------------------------------------------
          !write(6,*) 'EDenom',EDenom
          !write(6,*) 'Fac_ab',Fac_ab
          !write(6,*) 'xaibj',xaibj
          !write(6,*) 'xacbj',xacbj
          !write(6,*) 'xbiaj',xbiaj
          !write(6,*) 'xbcaj',xbcaj
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
          do iJ=1,nOcc(iSymI)
            EDenom = One/(work(mAdOcc(iSymJ)+iI-1)+work(mAdOcc(iSymI)+iJ-1)-work(mAdVir(iSymA)+iA-1)-work(mAdVir(iSymB)+iB-1))
            xaibj = Work(ipInt1+iJI+iJ-1+(iI-1)*(nOrb(iSymI)+nDel(iSymI)))
            xacbj = Work(ipInt1+iJC+iJ-1+(iC-1)*(nOrb(iSymI)+nDel(iSymI)))
            xajbi = Work(ipInt2+iIJ+iI-1+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            xajbc = Work(ipInt2+iCJ+iC-1+(iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
            Work(iMp2Lagr(iC,iI+nFro(iSymJ),iSymJ)) = Work(iMp2Lagr(iC,iI+nFro(iSymJ),iSymJ))+ &
                                                      EDenom*Fac_ab*(-Two*(xaibj*xacbj+xajbi*xajbc)+xaibj*xajbc+xajbi*xacbj)

            !-----------------------------------------------------------
            !write(6,*) 'EDenom',EDenom
            !write(6,*) 'Fac_ab',Fac_ab
            !write(6,*) 'xijab',xijab
            !write(6,*) 'xcjab',xcjab
            !write(6,*) 'xjiab',xjiab
            !write(6,*) 'xjcab',xjcab
            !-----------------------------------------------------------
          end do
        end do
      end do
    end if

  end do
end do

return

end subroutine rhs_mp2_help1
