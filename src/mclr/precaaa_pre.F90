!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine Precaaa_Pre(ActInt,A_J,Scr)

use MCLR_Data, only: nA
use input_mclr, only: nAsh, nBas, nIsh, nSym, ntAsh
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: ActInt(ntAsh,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(_OUT_) :: A_J(*), Scr(*)
integer(kind=iwp) :: iA, iAabs, iAtot, iB, iBabs, iBtot, iC, iCabs, iCtot, iD, iDabs, iDtot, iSym, jSym
real(kind=wp) :: Val

do iSym=1,nSym
  do iA=1,nAsh(iSym)
    iAabs = iA+nA(iSym)
    iAtot = iA+nIsh(iSym)
    do iB=1,nAsh(iSym) ! iA
      iBabs = iB+nA(iSym)
      iBtot = iB+nIsh(iSym)
      do jSym=1,iSym
        call Coul(iSym,iSym,jSym,jSym,iAtot,iBtot,A_J,Scr)
        do iC=1,nAsh(jSym)
          iCabs = iC+nA(jSym)
          iCtot = iC+nIsh(jSym)
          do iD=1,nAsh(jSym) ! iC
            iDabs = iD+nA(jSym)
            iDtot = iD+nIsh(jSym)
            Val = A_J(iCtot+nBas(jSym)*(iDtot-1))
            ActInt(iAabs,iBabs,iCabs,iDabs) = Val
            !ActInt(iBabs,iAabs,iCabs,iDabs) = Val
            !ActInt(iAabs,iBabs,iDabs,iCabs) = Val
            !ActInt(iBabs,iAabs,iDabs,iCabs) = Val
            !ActInt(iCabs,iDabs,iAabs,iBabs) = Val
            !ActInt(iCabs,iDabs,iBabs,iAabs) = Val
            !ActInt(iDabs,iCabs,iAabs,iBabs) = Val
            !ActInt(iDabs,iCabs,iBabs,iAabs) = Val
          end do
        end do
      end do
    end do
  end do
end do

end subroutine Precaaa_Pre
