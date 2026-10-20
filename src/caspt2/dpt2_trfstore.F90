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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DPT2_TrfStore(Scal,NBSQT,DPT2q,DPT2n,Trf,WRK)

use caspt2_module, only: NBAS, NDEL, NORB, NSYM
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT
real(kind=wp), intent(in) :: Scal, DPT2q(NBSQT), Trf(NBSQT)
real(kind=wp), intent(inout) :: DPT2n(NBSQT)
real(kind=wp), intent(out) :: WRK(NBSQT)
integer(kind=iwp) :: iMO, iSym, nOrbI

iMO = 1
do iSym=1,nSym
  if (nOrb(iSym) > 0) then
    nOrbI = nBas(iSym)-nDel(iSym)
    !! Quasi-canonical -> natural transformation of DPT2
    call DGemm_('N','N',nOrbI,nOrbI,nOrbI,One,Trf(iMO),nOrbI,DPT2q(iMO),nOrbI,Zero,WRK,nOrbI)
    call DGemm_('N','T',nOrbI,nOrbI,nOrbI,Scal,WRK,nOrbI,Trf(iMO),nOrbI,One,DPT2n(iMO),nOrbI)
  end if
  iMO = iMO+nOrbI*nOrbI
end do

return

end subroutine DPT2_TrfStore
