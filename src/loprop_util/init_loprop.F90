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

subroutine Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,LP_context,nSize,nBas1,nBas2,nBasMax)

use loprop_arrays, only: LP_context_type
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nSym, nBas(8), nOrb(8), nAtoms, nSize, nBas1, nBas2, nBasMax
real(kind=wp), intent(out) :: CoC(3)
type(LP_context_type), intent(out) :: LP_context
integer(kind=iwp) :: i, iDum, iSym
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iBas
#endif
real(kind=wp) :: DET
logical(kind=iwp) :: lOrb
integer(kind=iwp), parameter :: Occ = 1, Vir = 0

!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call qpg_iarray('nOrb',lOrb,iDum)
if (lOrb) then
  call Get_iArray('nOrb',nOrb,nSym)
else
  nOrb(1:nSym) = nBas(1:nSym)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nSize = 0
nBas1 = 0
nBas2 = 0
nBasMax = 0
do iSym=1,nSym
  nSize = nSize+nBas(iSym)*(nBas(iSym)+1)/2
  nBas1 = nBas1+nBas(iSym)
  nBas2 = nBas2+nBas(iSym)**2
  nBasMax = max(nBasMax,nBas(iSym))
end do
nSize = nSize+4
!                                                                      *
!***********************************************************************
!                                                                      *
! Center of Charge
call Get_dArray('Center of Charge',CoC,3)

! List coordinates of Coordinates
call Get_iScalar('LP_nCenter',nAtoms)
call mma_allocate(LP_context%C,3,nAtoms,label='C')
call Get_dArray('LP_Coor',LP_context%C,3*nAtoms)

! Effective charge at each center
call mma_allocate(LP_context%Q_Nuc,nAtoms,label='nAtoms')
call Get_dArray('LP_Q',LP_context%Q_Nuc,nAtoms)

! Atom number of each center
call mma_allocate(LP_context%ANr,nAtoms,label='ANr')
call Get_iArray('LP_A',LP_context%ANr,nAtoms)

! Pick up information of orbital type. Occ/Vir
call mma_allocate(LP_context%otype,nbas1,label='otype')
call Get_iArray('Orbital Type',LP_context%otype,nBas1)
do i=1,nBas1
  if ((LP_context%otype(i) /= Occ) .and. (LP_context%otype(i) /= Vir)) then
    write(u6,*) 'Orbital type vector is corrupted!'
    call Abend()
  end if
end do

! Pick up index array of which center a basis function belong.
call mma_allocate(LP_context%center,nbas1,label='center')
call Get_iArray('Center Index',LP_context%center,nBas1)

#ifdef _DEBUGPRINT_
write(u6,*) '******* LoProp Debug Info *******'
call RecPrt('Coordinates',' ',LP_context%C,3,nAtoms)
call RecPrt('Charges',' ',LP_context%Q_Nuc,1,nAtoms)
write(u6,*) 'Atom Nr:',LP_context%ANr(:)
do iBas=1,nBas1
  if (LP_context%otype(iBas) == Occ) then
    write(u6,'(A,I3,A)') 'Basis function ',iBas,' is occupied'
  else
    write(u6,'(A,I3,A)') 'Basis function ',iBas,' is virtual'
  end if
end do
write(u6,*)
do iBas=1,nBas1
  write(u6,'(A,I3,A,I3)') 'Basis function ',iBas,' belongs to center ',LP_context%center(iBas)
end do
write(u6,*) '*********************************'
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of symmetry we need the desymmetrization matrix

if (nSym /= 1) then
  call mma_allocate(LP_context%P,nbas1,nbas1,label='P')
  call mma_allocate(LP_context%PInv,nbas1,nbas1,label='PInv')
  call Get_dArray('SM',LP_context%P,nbas1**2)
# ifdef _DEBUGPRINT_
  call RecPrt('SM',' ',LP_context%P,nbas1,nbas1)
# endif
  call MINV(LP_context%P,LP_context%PInv,DET,nBas1)
# ifdef _DEBUGPRINT_
  call RecPrt('SMInv',' ',LP_context%PInv,nbas1,nbas1)
# endif
  call DGeTMi(LP_context%PInv,nbas1,nbas1)
else
  call mma_allocate(LP_context%P,0,0,label='P')
  call mma_allocate(LP_context%PInv,0,0,label='PInv')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Init_LoProp
