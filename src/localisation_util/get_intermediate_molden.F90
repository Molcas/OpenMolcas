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
!***********************************************************************

subroutine get_intermediate_molden(CMO,nBasis,nOrb2Loc,nIter)
! call this subroutine to view intermediately localized orbitals at one (!) iteration of interest

! creates an additional MD_LOC (MD_LOIM) file, with extension .imlocal.molden

! the code is mostly copied from localisation.F90

use Localisation_globals, only: EOrb, Ind, nBas, nOrb, nSym, Occ, Silent
use Definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc, nIter
real(kind=wp), intent(in) :: CMO(nBasis,norb2Loc)
integer(kind=iwp) :: IndT(7,8), iSym, iUHF, j, k, kIndT, LU_
character(len=80) :: Title
character(len=9) :: Filename
integer(kind=iwp), external :: isFreeUnit

! Write TMPORB file.
! ------------------

Title = ''
write(Title,'(A)') 'Intermediate Orbitals'
LU_ = isFreeUnit(12)
j = 0
IndT(:,:) = 0
do iSym=1,nSym
  do k=1,nOrb(iSym)
    kIndT = Ind(j+k)
    if ((kIndT > 0) .and. (kIndT <= 7)) then
      IndT(kIndT,iSym) = IndT(kIndT,iSym)+1
    else
      call WarningMessage(2,'Localisation: Illegal orbital type')
      write(u6,'(A,I6,A,I2,A,I9)') 'Orbital',k,' of sym.',iSym,' has illegal type:',kIndT
      call Abend()
    end if
  end do
  j = j+nBas(iSym)
end do
call WrVec_Localisation('TMPORB',LU_,'COEI',nSym,nBas,nBas,CMO,Occ,EOrb,IndT,Title)
if (.not. Silent) then
  write(u6,'(1X,A)') 'The TMPORB file has been written.'
  write(u6,'(3(A),I5)') ' Title=',trim(Title),' LU_=',LU_
end if

! Write MOLDEN file.
! ------------------

iUHF = 0
write(Filename,'(A,I0)') 'MD_LOC.',nIter
call Molden_Interface(iUHF,'TMPORB',Filename)
if (.not. Silent) write(u6,'(1X,A)') 'The MOLDEN file has been written.'

end subroutine get_intermediate_molden
