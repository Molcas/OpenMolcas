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

subroutine Diff_Numerical(nAt,nB,MP,nij,EC,iANr,Ttot,Ttot_Inv,lMax,TP,dLimmo,Thrs1,Thrs2,nThrs,iPrint,ThrsMul,Pot_Expo,Pot_Point, &
                          Pot_Fac,Diffed)

use Index_Functions, only: nTri3_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt, nB, nij, iANr(nAt), lMax, nThrs, iPrint
real(kind=wp), intent(in) :: MP(nij,*), EC(3,nij), Ttot(nB,nB), Ttot_Inv(nB,nB), TP(nAt), dLimmo(2), Thrs1, Thrs2, ThrsMul
real(kind=wp), intent(out) :: Pot_Expo(nij*2), Pot_Point(nij), Pot_Fac(nij*4)
logical(kind=iwp), intent(out) :: Diffed(nij*2)
integer(kind=iwp) :: iAtom, iDC, ij, iK, irc, jAtom, k, kaunt, kauntA, kComp, l, nAbove, nEPP, nK, nPick
real(kind=wp) :: A(2), BS, Chi2B, chPoint, dM, dMag, ThrsMul_Clever
logical(kind=iwp) :: AboveMul(2), AboveThr
character(len=50) :: UtChar
character(len=10) :: OneFile
integer(kind=iwp), allocatable :: Center(:), Pick(:)
real(kind=wp), allocatable :: dMullig(:), DPick(:), EPCo(:,:), Potte(:)
real(kind=wp), external :: vdwRad
interface
  subroutine Diff_Aux1(nEPotPoints,EPCo,nB,OneFile)
    import :: wp, iwp
    integer(kind=iwp), intent(out) :: nEPotPoints
    real(kind=wp), allocatable, intent(out) :: EPCo(:,:)
    integer(kind=iwp), intent(in) :: nB
    character(len=10), intent(in) :: OneFile
  end subroutine Diff_Aux1
end interface

! Pick up some auxiliary stuff.

write(OneFile,'(A)') 'ONEINT'
call Diff_Aux1(nEPP,EPCo,nB,OneFile)
call mma_allocate(Center,nB,label='BasIndCent')
call Get_iArray('Center Index',Center,nB)
call mma_allocate(Pick,nEPP,label='PickPoints')
call mma_allocate(DPick,nEPP,label='DistPick')
call mma_allocate(dMullig,nTri3_Elem1(lMax),label='dMullig')

!-- Do a 'clever' determination of the threshold for the multipole magnitude.

!call Diff_ThrsMul(MP,ThrsMul,ThrsMul_Clever,nAt,nij)
ThrsMul_Clever = ThrsMul

! Run a numerical fit for each active centre.

kauntA = 1
nAbove = 0
do iAtom=1,nAt
  do jAtom=1,iAtom
    ij = iAtom*(iAtom-1)/2+jAtom

    ! Pick up the nuclei+core charge.

    if (iAtom == jAtom) then
      chPoint = TP(iAtom)
    else
      chPoint = Zero
    end if

    ! Pick out the multipole, the prefactors. If none is above a
    ! certain threshold, then it is not meaningful to make them
    ! diffuse. Also check which individual multipoles that should
    ! be made diffuse.

    kaunt = 1
    AboveThr = .false.
    do l=0,lMax
      kComp = (l+1)*(l+2)/2
      dMag = Zero
      do k=1,kComp
        dM = MP(kauntA,kaunt)
        dMullig(kaunt) = dM
        dMag = dMag+dM**2
        kaunt = kaunt+1
      end do
      dMag = sqrt(dMag)
      if ((dMag > ThrsMul_Clever) .and. (l <= 1)) then
        AboveThr = .true.
        AboveMul(l+1) = .true.
      else if ((dMag <= ThrsMul_Clever) .and. (l <= 1)) then
        AboveMul(l+1) = .false.
      end if
    end do

    if (AboveThr) then
      nAbove = nAbove+1

      ! Select the potential points which should be used for this centre.

      !BS = Half*(Bragg_Slater(iANr(iAtom))+Bragg_Slater(iANr(jAtom)))
      BS = Half*(vdWRad(iANr(iAtom))+vdWRad(iANr(jAtom)))
      call PickPoints(nPick,Pick,DPick,nEPP,EPCo,EC(:,ij),dLimmo,BS)

      ! Compute the true potential from the density assigned to this centre.

      call mma_allocate(Potte,nPick,label='Potential')
      call EPotPoint(Potte,nPick,Pick,DPick,Ttot,Ttot_Inv,iANr(iAtom),nB,iAtom,jAtom,Center)

      ! Print the true potential for given centre if requested.

      if (iPrint >= 5) then
        write(UtChar,'(A,2I3)') 'Partial density potential, centre',iAtom,jAtom
        call RecPrt(UtChar,' ',Potte,nPick,1)
      end if

      ! All hail to the chiefs, i.e. Levenberg and Marquardt.

      call LevMarquart(Potte,nPick,Pick,EPCo,EC(:,ij),dMullig,lMax,A,iAtom,jAtom,chPoint,Thrs1,Thrs2,nThrs,Chi2B,iPrint,AboveMul)
      call mma_deallocate(Potte)
    end if

    ! Store things for later.

    kaunt = 1
    Pot_Point(kauntA) = chPoint
    do iDC=1,2
      nK = iDC*(iDC+1)/2
      do iK=1,nK
        Pot_Fac(4*(kauntA-1)+kaunt) = dMullig(kaunt)
        kaunt = kaunt+1
      end do
      if (.not. AboveThr) then
        Diffed(2*(kauntA-1)+iDC) = .false.
      else
        if ((A(iDC) < Three) .and. AboveMul(iDC)) then
          Diffed(2*(kauntA-1)+iDC) = .true.
          Pot_Expo(2*(kauntA-1)+iDC) = A(iDC)
        else
          Diffed(2*(kauntA-1)+iDC) = .false.
          Pot_Expo(2*(kauntA-1)+iDC) = Ten !Dummy
        end if
      end if
    end do

    ! Step the atom-pair counter.

    kauntA = kauntA+1
  end do
end do

! Deallocations.

call mma_deallocate(Center)
call mma_deallocate(Pick)
call mma_deallocate(DPick)
call mma_deallocate(EPCo)
call mma_deallocate(dMullig)
irc = -1
call ClsOne(irc,0)

return

end subroutine Diff_Numerical
