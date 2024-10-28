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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine MODens()
!***********************************************************************
!                                                                      *
!     purpose: Compute density matrix in molecular orbital basis       *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use InfSCF, only: CMO, Dens, DMOMax, MaxBas, MaxBXO, MaxOrb, nBas, nD, nDens, nOcc, nOrb, nSym, Ovrlp, TEEE, TimFld
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use InfSCF, only: nBO, nBT
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: i, iD, iiBO, iiBT, iOvl, iSym, iT, j, jD
real(kind=wp) :: CPU1, CPU2, Tim1, Tim2, Tim3
real(kind=wp), allocatable :: Aux1(:), Aux2(:), DMoO(:), DnsS(:), OvlS(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Timing(Cpu1,Tim1,Tim2,Tim3)

! Allocate memory for squared density matrix
call mma_allocate(DnsS,MaxBas**2,Label='DnsS')

! Allocate memory for squared overlap matrix
call mma_allocate(OvlS,MaxBas**2,Label='OvlS')

! Allocate memory for density matrix in MO basis
call mma_allocate(DMoO,nTri_Elem(MaxOrb),Label='DMoO')

! Allocate memory for auxiliary matrix
call mma_allocate(Aux1,MaxBxO,Label='Aux1')

! Allocate memory for auxiliary matrix
call mma_allocate(Aux2,MaxBxO,Label='Aux2')

DMOMax = Zero
do jD=1,nD

# ifdef _DEBUGPRINT_
  call NrmClc(Dens(1,jD,nDens),nBT,'MoDens','D in AO   ')
  call NrmClc(Ovrlp,nBT,'MoDens','Overlap   ')
  call NrmClc(CMO(1,jD),nBO,'MoDens','CMOs      ')
  write(u6,*) 'nOcc=',(nOcc(i,jD),i=1,nSym)
  !write(u6,'(F16.8)') DXot(MaxBxO,CMO(1,jD),1,CMO(1,jD),1)
# endif
  it = 1
  id = 1
  iOvl = 1
  do iSym=1,nSym

    iiBO = nBas(iSym)*nOrb(iSym)
    iiBT = nTri_Elem(nBas(iSym))

    if ((nOcc(iSym,jD) > 0) .or. (Teee .and. (nBas(iSym) > 0))) then
      call DSq(Dens(id,jD,nDens),DnsS,1,nBas(iSym),nBas(iSym))
      call Square(Ovrlp(iOvl),OvlS,1,nBas(iSym),nBas(iSym))

      call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym), &
                  One,OvlS,nBas(iSym), &
                  CMO(it,jD),nBas(iSym), &
                  Zero,Aux1,nBas(iSym))
      call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym), &
                  One,DnsS,nBas(iSym), &
                  Aux1,nBas(iSym), &
                  Zero,Aux2,nBas(iSym))
      call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym), &
                  One,OvlS,nBas(iSym), &
                  Aux2,nBas(iSym), &
                  Zero,Aux1,nBas(iSym))
      call DGEMM_Tri('T','N',nOrb(iSym),nOrb(iSym),nBas(iSym), &
                     One,CMO(it,jD),nBas(iSym), &
                     Aux1,nBas(iSym), &
                     Zero,DMoO,nOrb(iSym))

      !call TriPrt('D(mo)','(8F12.6)',DMoO,nOrb(iSym))
      do i=nOcc(iSym,jD)+1,nOrb(iSym)
        do j=1,nOcc(iSym,jD)
          DMOMax = max(DMOMax,real(nD,kind=wp)*abs(DMoO(iTri(i,j))))
        end do
      end do
    end if

    it = it+iiBO
    id = id+iiBT
    iOvl = iOvl+iiBT
  end do

end do

! Deallocate memory
call mma_deallocate(Aux2)
call mma_deallocate(Aux1)
call mma_deallocate(DMoO)
call mma_deallocate(OvlS)
call mma_deallocate(DnsS)

#ifdef _DEBUGPRINT_
write(u6,*) ' DMOMax in MODens',DMOMax
#endif
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(14) = TimFld(14)+(Cpu2-Cpu1)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine MODens
