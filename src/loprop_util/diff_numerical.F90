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

subroutine Diff_Numerical(nAt,nB,ipMP,ipC,nij,EC,iANr,ip_Ttot,ip_Ttot_Inv,lMax,iTP,dLimmo,Thrs1,Thrs2,nThrs,iPrint,ThrsMul, &
                          Pot_Expo,Pot_Point,Pot_Fac,Diffed)

use Constants, only: Zero, Three, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAt, nB, ipMP, ipC, nij, iANr(nAt), ip_Ttot, ip_Ttot_Inv, lMax, iTP, nThrs, iPrint
real(kind=wp) :: EC(3,nij), dLimmo(2), Thrs1, Thrs2, ThrsMul, Pot_Expo(nij*2), Pot_Point(nij), Pot_Fac(nij*4)
logical(kind=iwp) :: Diffed(nij*2)
integer(kind=iwp) :: iAtom, iDC, ij, iK, ip_Center, ipDPick, ipEPCo, iPotte, ipPick, irc, jAtom, k, kaunt, kauntA, kComp, l, &
                     nAbove, nEPP, nK, nPick
real(kind=wp) :: A(2), BS, Chi2B, chPoint, dM, dMag, dMullig((lMax*(lMax**2+6*lMax+11)+6)/6), ThrsMul_Clever
logical(kind=iwp) :: AboveMul(2), AboveThr
character(len=50) :: UtChar
character(len=10) :: OneFile
real(kind=wp), external :: vdwRad
#include "WrkSpc.fh"

! Pick up some auxiliary stuff.

write(OneFile,'(A)') 'ONEINT'
call Diff_Aux1(nEPP,ipEPCo,nB,OneFile)
call GetMem('BasIndCent','Allo','Inte',ip_Center,nB)
call Get_iArray('Center Index',iWork(ip_Center),nB)
call GetMem('PickPoints','Allo','Inte',ipPick,nEPP)
call GetMem('DistPick','Allo','Real',ipDPick,nEPP)

!-- Do a 'clever' determination of the threshold for the multipole magnitude.

!call Diff_ThrsMul(ipMP,ThrsMul,ThrsMul_Clever,nAt,nij,lMax)
ThrsMul_Clever = ThrsMul

! Run a numerical fit for each active centre.

kauntA = 0
nAbove = 0
do iAtom=1,nAt
  do jAtom=1,iAtom
    ij = iAtom*(iAtom-1)/2+jAtom

    ! Pick up the nuclei+core charge.

    if (iAtom == jAtom) then
      chPoint = Work(iTP+iAtom-1)
    else
      chPoint = Zero
    end if

    ! Pick out the multipole, the prefactors. If none is above a
    ! certain threshold, then it is not meaningful to make them
    ! diffuse. Also check which individual multipoles that should
    ! be made diffuse.

    kaunt = 0
    AboveThr = .false.
    do l=0,lMax
      kComp = (l+1)*(l+2)/2
      dMag = Zero
      do k=1,kComp
        dM = Work(ipMP+nij*kaunt+kauntA)
        dMullig(kaunt+1) = dM
        dMag = dMag+dM**2
        kaunt = kaunt+1
      end do
      dMag = sqrt(dMag)
      if ((dMag > ThrsMul_Clever) .and. (l <= 1)) then
        AboveThr = .true.
        AboveMul(l+1) = .true.
      elseif ((dMag <= ThrsMul_Clever) .and. (l <= 1)) then
        AboveMul(l+1) = .false.
      end if
    end do

    if (AboveThr) then
      nAbove = nAbove+1

      ! Select the potential points which should be used for this centre.

      !BS = Half*(Bragg_Slater(iANr(iAtom))+Bragg_Slater(iANr(jAtom)))
      BS = Half*(vdWRad(iANr(iAtom))+vdWRad(iANr(jAtom)))
      call PickPoints(nPick,ipPick,ipDPick,nEPP,ipEPCo,EC(1,ij),dLimmo,BS)

      ! Compute the true potential from the density assigned to this centre.

      call GetMem('Potential','Allo','Real',iPotte,nPick)
      call EPotPoint(iPotte,nPick,ipPick,ipDPick,nEPP,ip_Ttot,ip_Ttot_Inv,iANr(iAtom),nB,iAtom,jAtom,ip_Center)

      ! Print the true potential for given centre if requested.

      if (iPrint >= 5) then
        write(UtChar,'(A,2I3)') 'Partial density potential, centre',iAtom,jAtom
        call RecPrt(UtChar,' ',Work(iPotte),nPick,1)
      end if

      ! All hail to the chiefs, i.e. Levenberg and Marquardt.

      call LevMarquart(iPotte,nPick,ipPick,nEPP,ipEPCo,EC(1,ij),dMullig,lMax,A,iAtom,jAtom,chPoint,Thrs1,Thrs2,nThrs,Chi2B,iPrint, &
                       AboveMul)
      call GetMem('Potential','Free','Real',iPotte,nPick)
    end if

    ! Store things for later.

    kaunt = 0
    Pot_Point(kauntA+1) = chPoint
    do iDC=1,2
      nK = iDC*(iDC+1)/2
      do iK=1,nK
        kaunt = kaunt+1
        Pot_Fac(4*kauntA+kaunt) = dMullig(kaunt)
      end do
      if (.not. AboveThr) then
        Diffed(2*kauntA+iDC) = .false.
      else
        if ((A(iDC) < Three) .and. AboveMul(iDC)) then
          Diffed(2*kauntA+iDC) = .true.
          Pot_Expo(2*kauntA+iDC) = A(iDC)
        else
          Diffed(2*kauntA+iDC) = .false.
          Pot_Expo(2*kauntA+iDC) = Ten !Dummy
        end if
      end if
    end do

    ! Step the atom-pair counter.

    kauntA = kauntA+1
  end do
end do

! Deallocations.

call GetMem('BasIndCent','Free','Inte',ip_Center,nB)
call GetMem('PickPoints','Free','Inte',ipPick,nEPP)
call GetMem('DistPick','Free','Real',ipDPick,nEPP)
call GetMem('PotPointCoord','Free','Real',ipEPCo,3*nEPP)
irc = -1
call ClsOne(irc,0)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(ipC)

end subroutine Diff_Numerical
