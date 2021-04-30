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

subroutine FCIN(FLT,nFLT,DLT,FSQ,DSQ,EMY,CMO)

#include "intent.fh"

use motra_global, only: Debug, iPrint, nBas, nFro, nSym, nTot1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: nFLT
real(kind=wp), intent(inout) :: FLT(nFLT)
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: DLT(*), FSQ(*), DSQ(*)
real(kind=wp), intent(out) :: EMY
integer(kind=iwp) :: ISTLTT, ISYM, LBUF, n_Bas, NB, NPQ, NTFRO
real(kind=wp) :: EONE, ETWO
real(kind=wp), allocatable :: Temp(:), W1(:), W2(:)
logical(kind=iwp) :: DoCholesky
integer(kind=iwp), external :: mma_avmem

! Construct the one-electron density matrix for frozen space

call DONEI(DLT,DSQ,CMO)

! Compute the one-electron energy contribution to EMY

EONE = Zero
do NPQ=1,NTOT1
  EONE = EONE+DLT(NPQ)*FLT(NPQ)
end do
EMY = EONE
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A,E20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
end if

! Quit here if there are no frozen orbitals

NTFRO = 0
n_Bas = 0
do ISYM=1,NSYM
  NTFRO = NTFRO+NFRO(ISYM)
  n_Bas = max(n_Bas,nBas(iSym))
end do
if (NTFRO == 0) then
  return
end if

! Compute the two-electron contribution to the Fock matrix

call mma_allocate(Temp,nFlt,label='Temp')
Temp(:) = Zero

call DecideOnCholesky(DoCholesky)

if (DoCholesky) then
  call Cho_Fock_MoTra(nSym,nBas,nFro,DLT,DSQ,FLT,nFLT,FSQ,One)

  ! Print the Fock-matrix

  if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
    write(u6,'(6X,A)') 'Fock matrix in AO basis'
    ISTLTT = 1
    do ISYM=1,NSYM
      NB = NBAS(ISYM)
      if (NB > 0) then
        write(u6,'(6X,A,I2)') 'symmetry species:',ISYM
        call TRIPRT(' ',' ',Temp(ISTLTT),NB) ! ??? Temp is zero
        ISTLTT = ISTLTT+NB*(NB+1)/2
      end if
    end do
  end if

else

  call mma_allocate(W2,n_Bas**2,label='FCIN2')
  LBUF = int(mma_avmem()*0.9_wp,kind=iwp)/RtoB
  call mma_allocate(W1,LBUF,label='FCIN1')

  call FTWOI(DLT,DSQ,Temp,nFlt,FSQ,LBUF,W1,W2)

  call mma_deallocate(W2)
  call mma_deallocate(W1)
end if

call DaXpY_(nFlt,One,Temp,1,Flt,1)
call mma_deallocate(Temp)

! Add the two-electron contribution to EMY
ETWO = -EONE
do NPQ=1,NTOT1
  ETWO = ETWO+DLT(NPQ)*FLT(NPQ)
end do
EMY = EONE+Half*ETWO
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A,E20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
end if

! Exit

return

end subroutine FCIN
