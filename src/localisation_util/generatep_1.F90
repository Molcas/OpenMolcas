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
! Copyright (C) Yannick Carissan                                       *
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GenerateP_1(Ovlp,cMO,Sbar,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Reduce operation count and use BLAS.

use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: nBasis, nOrb2Loc, nAtoms, nBas_per_Atom(*), nBas_Start(*)
real(kind=wp) :: Ovlp(nBasis,nBasis), cMO(nBasis,*), Sbar(nBasis,nOrb2Loc), PA(nOrb2Loc,nOrb2Loc,nAtoms)
character(len=LenIn8) :: BName(*)
logical(kind=iwp) :: Debug
#include "WrkSpc.fh"
#include "real.fh"
integer(kind=iwp) :: iAt, iMO_s, iMO_t
real(kind=wp) :: PAst, PAts
character(len=LenIn8) :: PALbl

! Compute Sbar(mu,s) = sum_{nu} Ovlp(mu,nu) * cMO(nu,s)

call DGEMM_('N','N',nBasis,nOrb2Loc,nBasis,One,Ovlp,nBasis,cMO,nBasis,Zero,Sbar,nBasis)

do iAt=1,nAtoms

  ! Compute MA(s,t) = sum_{mu_in_A} cMO(mu,s) * Sbar(mu,t)

  call DGEMM_('T','N',nOrb2Loc,nOrb2Loc,nBas_per_Atom(iAt),One,cMO(nBas_Start(iAt),1),nBasis,Sbar(nBas_Start(iAt),1),nBasis,Zero, &
              PA(1,1,iAt),nOrb2Loc)

  ! Compute <s|PA|t> by symmetrization of MA.

  do iMO_s=1,nOrb2Loc
    do iMO_t=iMO_s+1,nOrb2Loc
      PAst = PA(iMO_s,iMO_t,iAt)
      PAts = PA(iMO_t,iMO_s,iAt)
      PA(iMO_s,iMO_t,iAt) = Half*(PAst+PAts)
      PA(iMO_t,iMO_s,iAt) = PA(iMO_s,iMO_t,iAt)
    end do !iMO_t
  end do !iMO_s

end do !iAt

if (Debug) then
  write(u6,*) 'In GenerateP'
  write(u6,*) '------------'
  do iAt=1,nAtoms
    PALbl = 'PA__'//BName(nBas_Start(iAt))(1:LenIn)
    call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
  end do
end if

return

end subroutine GenerateP_1
