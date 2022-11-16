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

subroutine No_ESPF(natom,Forces,DoTinker)

use Basis_Info
use Center_Info
use external_centers
use Gateway_global, only: Primitive_pass
use Gateway_Info, only: PotNuc
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
#include "espf.fh"
real*8 A(3), B(3), RB(3)
integer iDCRR(0:7), jCoSet(8,8), iStb(0:7)
logical EQ, NoLoop, Forces, DoTinker, DoRys
character*180 Line, Get_Ln
external Get_Ln

iPL = iPL_espf()
NoLoop = .true.

! Nothing to do in a gradient calc

if (Forces) return

! Update the nuclear repulsion energy

call Get_dScalar('PotNuc',RepNuc)
RepNuc_old = RepNuc

! Read the MM contribution to the total energy and add it
! to the Nuclear Repulsion term

if (DoTinker) then
  Tke = Zero
  ITkQMMM = IsFreeUnit(30)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'MMEnergy ') /= 0) call Get_F1(2,TkE)
  end do
  close(ITkQMMM)
  TkE = TkE*ToHartree
  RepNuc = RepNuc+TkE
  if (iPL >= 3) write(6,1000) RepNuc_old,TkE,RepNuc
end if

if (allocated(XF) .and. (nOrd_XF >= 0)) then
  write(6,*) 'Here we are!!'

  DoRys = .true.
  nDiff = 0
  call GetInf(DoRys,nDiff)
  Primitive_Pass = .true.

  ! Add contribution for interaction external field and nuclear
  ! charges. Here we will have charge-charge, and charge-dipole
  ! interaction.

  ZA = Zero
  DAx = Zero
  DAy = Zero
  DAz = Zero
  Qxx = Zero
  Qxy = Zero
  Qxz = Zero
  Qyy = Zero
  Qyz = Zero
  Qzz = Zero

  PNX = Zero
  iDum = 0
  do iFd=1,nXF
    if (nOrd_XF == 0) then
      ZA = XF(4,iFd)
      NoLoop = ZA == Zero
    else if (nOrd_XF == 1) then
      ZA = XF(4,iFd)
      DAx = XF(5,iFd)
      DAy = XF(6,iFd)
      DAz = XF(7,iFd)
      NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero)
    else if (nOrd_XF == 2) then
      ZA = XF(4,iFd)
      DAx = XF(5,iFd)
      DAy = XF(6,iFd)
      DAz = XF(7,iFd)
      Qxx = XF(8,iFd)
      Qxy = XF(9,iFd)
      Qxz = XF(10,iFd)
      Qyy = XF(11,iFd)
      Qyz = XF(12,iFd)
      Qzz = XF(13,iFd)
      NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero) .and. (Qxx == Zero) .and. (Qxy == Zero) &
               .and. (Qxz == Zero) .and. (Qyy == Zero) .and. (Qyz == Zero) .and. (Qzz == Zero)
    else
      call WarningMessage(2,'Option not implemented yet!')
      call Quit_OnUserError()
    end if
    if (NoLoop) Go To 102
    A(1:3) = XF(1:3,iFd)
    iChxyz = iChAtm(A)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    ndc = 0
    do jCnttp=1,nCnttp
      ZB = dbsc(jCnttp)%Charge
      if (dbsc(jCnttp)%pChrg) Go To 202
      if (ZB == Zero) Go To 202
      if (dbsc(jCnttp)%Frag) Go To 202
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        temp0 = Zero
        temp1 = Zero
        temp2 = Zero
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          ! The index A=RB is illegal.
          if (.not. EQ(A,RB)) then
            ABx = A(1)-RB(1)
            ABy = A(2)-RB(2)
            ABz = A(3)-RB(3)
            r12 = sqrt(ABx**2+ABy**2+ABz**2)

            fab = One
            if (dbsc(jCnttp)%ECP) then
              ! Add contribution from M1 operator
              do iM1xp=1,dbsc(jCnttp)%nM1
                Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gamma*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(jCnttp)%nM2
                Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gamma*r12**2)
              end do
            end if
            temp0 = temp0+fab/r12
            if (nOrd_XF >= 1) temp1 = temp1-fab*(DAx*ABx+DAy*ABy+DAz*ABz)/r12**3
            if (nOrd_XF >= 2) then

              temp2 = temp2+fab*0.5d0*(3.0d0*(Qxx*ABx**2+2.0d0*Qxy*ABx*ABy+2.0d0*Qxz*ABx*ABz+Qyy*ABy**2+ &
                      2.0d0*Qyz*ABy*ABz+Qzz*ABz**2)/r12**5-One/r12**3*(Qxx+Qyy+Qzz))
            end if

          end if
        end do
        PNX = PNX+((ZAZB*temp0+ZB*(temp1+temp2))*dble(nIrrep))/dble(LmbdR)

      end do
202   continue
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
102 continue
  end do

  if (iPL >= 3) write(6,1100) RepNuc,PNX,RepNuc+PNX

  PotNuc = PotNuc+PNX
  call Put_dScalar('PotNuc',RepNuc)
end if

! Update the 1-e integrals

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(natom)

1000 format(/,' RepNuc + MM = ',F13.8,' + ',F13.8,' = ',F13.8)
1100 format(/,' RepNuc + Point charges = ',F13.8,' + ',F13.8,' = ',F13.8)

end subroutine No_ESPF
