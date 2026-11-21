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
! ActInt(a,b,c,d) = (ab|cd) for all active orbitals a, b, c, d, with the
! active orbitals in symmetry order (index nA(iSym)+iAsh). Used by Precaaa.
!
! All symmetry combinations with Mul(Mul(sym(a),sym(b)),Mul(sym(c),sym(d))) = 1
! are needed (Precaaa uses Coulomb- as well as exchange-type elements);
! the remaining elements vanish by symmetry and are set to zero.

use Symmetry_Info, only: Mul
use MCLR_Data, only: nA
use Orbital_Spaces, only: nAsh, nIsh, nSym
use input_mclr, only: nOrb, ntAsh
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: ActInt(ntAsh,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(_OUT_) :: A_J(*), Scr(*)
integer(kind=iwp) :: iA, iAabs, iB, iBabs, iC, iCabs, iCtot, iD, iDabs, iDtot, iSa, iSb, iSc, iSd

ActInt(:,:,:,:) = Zero
do iSc=1,nSym
  do iSd=1,nSym
    do iC=1,nAsh(iSc)
      iCabs = iC+nA(iSc)
      iCtot = iC+nIsh(iSc)
      do iD=1,nAsh(iSd)
        iDabs = iD+nA(iSd)
        iDtot = iD+nIsh(iSd)
        do iSa=1,nSym
          iSb = Mul(Mul(iSa,iSc),iSd)
          if ((nAsh(iSa) == 0) .or. (nAsh(iSb) == 0)) cycle
          ! A_J(p,q) = (pq|cd), p in iSa, q in iSb, leading dimension nOrb(iSa)
          call Coul(iSa,iSb,iSc,iSd,iCtot,iDtot,A_J,Scr)
          do iB=1,nAsh(iSb)
            iBabs = iB+nA(iSb)
            do iA=1,nAsh(iSa)
              iAabs = iA+nA(iSa)
              ActInt(iAabs,iBabs,iCabs,iDabs) = A_J(iA+nIsh(iSa)+nOrb(iSa)*(iB+nIsh(iSb)-1))
            end do
          end do
        end do
      end do
    end do
  end do
end do

end subroutine Precaaa_Pre
