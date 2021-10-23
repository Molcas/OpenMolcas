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

subroutine CHO_eval_waxy(irc,Scr,ChoV1,ChoV2,W_PWXY,nAorb,JSYM,NUMV,DoTraInt,CMO)

use Data_structures, only: SBA_Type, twxy_Type, DSBA_Type

implicit real*8(a-h,o-z)
real*8 W_PWXY(*)
integer nAorb(*)
type(DSBA_Type) CMO
type(SBA_Type) ChoV1, ChoV2
type(twxy_type) Scr
integer off_PWXY(8,8,8)
logical DoTraInt
#include "real.fh"
#include "Molcas.fh"
#include "general.fh"
#include "wadr.fh"
!************************************************
!Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1
!************************************************

if (NumV < 1) return

! Computing the integrals
! -----------------------------------------------------
! (wa|xy)  <-  (wa|xy)  +  sum_J  L(wa,#J) * L(xy,#J)
!======================================================

do iSymy=1,nSym

  iSymx = MulD2h(iSymy,JSYM)

  Nxy = size(ChoV2%SB(iSymx)%A2,1)

  if ((iSymx <= iSymy) .and. (Nxy > 0)) then

    do iSyma=1,nSym

      iSymw = MulD2h(iSyma,JSYM)

      Nwa = size(ChoV1%SB(iSymw)%A3,1)*size(ChoV1%SB(iSymw)%A3,2)

      if (Nwa <= 0) cycle

      call DGEMM_('N','T',Nwa,Nxy,NumV,ONE,ChoV1%SB(iSymw)%A3,Nwa,ChoV2%SB(iSymx)%A2,Nxy,ONE,Scr%SB(iSymw,iSymx)%A,Nwa)

    end do

  end if

end do

! MO-transformation to build the (pw|xy) integrals
! is performed (p is a general index).
! The storage of the result is the one required by
! the RASSCF program : LT-storage in xy, with sym(x)<=sym(y)
!
!       (pw|xy)  =  sum_a  C(a,p) * (wa|xy)
!
!---------------------------------------------------------
if (DoTraInt) then

  ! generate offsets to (pw|xy)
  !
  ! Important: the way GET_TUVX is written requires
  !            a strict loop structure for the definition
  !            of the offsets

  iStack = 0
  do iSymp=1,nSym
    iOrb = nOrb(iSymp)
    do iSymw=1,nSym
      jAsh = nAorb(iSymw)
      ijSym = MulD2h(iSymp,iSymw)
      do iSymy=1,nSym
        kAsh = nAorb(iSymy)
        iSymx = MulD2h(ijSym,iSymy)
        if (iSymx <= iSymy) then
          lAsh = nAorb(iSymx)
          kl_Orb_pairs = kAsh*lAsh+min(0,ijSym-2)*kAsh*(lAsh-1)/2
          off_PWXY(iSymp,iSymw,iSymx) = iStack
          iStack = iStack+iOrb*jAsh*kl_Orb_pairs
        end if
      end do
    end do
  end do
  nPWXY = iStack
  !
  ! Reordering and MO-transformation
  !
  do iSymy=1,nSym

    iSymx = MulD2h(iSymy,JSYM)

    if (iSymx <= iSymy) then

      Nxy = nAorb(iSymx)*nAorb(iSymy)+min(0,JSYM-2)*nAorb(iSymx)*(nAorb(iSymy)-1)/2

      do iSymw=1,nSym

        iSyma = MulD2h(iSymw,JSYM) ! =iSymp

        Nwa = nAorb(iSymw)*nBas(iSyma)
        Npw = nOrb(iSyma)*nAorb(iSymw)

        iS = 1+nBas(iSyma)*nFro(iSyma)

        do ixy=1,Nxy

          ipMpw = off_PWXY(iSyma,iSymw,iSymx)+1+Npw*(ixy-1)

          ! --------------------------------------------------------
          !       M(p,w)[xy]  =  sum_a  C(a,p) * M(w,a)[xy]
          ! --------------------------------------------------------
          nBas_a = max(nBas(iSyma),1)
          nAob_w = max(nAorb(iSymw),1)
          nOrb_a = max(nOrb(iSyma),1)

          call DGEMM_('T','T',nOrb(iSyma),nAorb(iSymw),nBas(iSyma),ONE,CMO%SB(iSyma)%A1(iS:),nBas_a,Scr%SB(iSymw,iSymx)%A(:,ixy), &
                      nAob_w,ZERO,W_PWXY(ipMpw),nOrb_a)

        end do

      end do

    end if

  end do

end if

irc = 0

return

end subroutine CHO_eval_waxy
