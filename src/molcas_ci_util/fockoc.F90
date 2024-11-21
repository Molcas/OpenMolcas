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

subroutine FOCKOC(FOCC,F,CMO)
! RASSCF program version IBM-3090: SX section
!
! PURPOSE: Construct a fock matrix for the occupied orbitals and
!          write it on the job interphase for later use in
!          the RASSCF gradient programs.
!          The fock matrix is in MO basis with the frozen orbitals
!          excluded. It is symmetry blocked in contrast to earlier
!          versions.
!
! called from FOCK if IFINAL=1.
!
! ********** IBM-3090 MOLCAS Release: 90 02 22 **********

use Index_Functions, only: nTri_Elem
use wadr, only: FockOcc
use rasscf_global, only: IADR15, NO2M
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: FOCC(*)
real(kind=wp), intent(in) :: F(*), CMO(*)
#include "rasdim.fh"
#include "general.fh"
real(kind=wp), allocatable :: SCR1(:), SCR2(:)
integer(kind=iwp) :: IAD15, iBas, iCMO, ij, IPQ, ISTFCK, ISYM, jBas, jFock, kl, lk, NAO, NIO, NO, NOO, NP, NQ

ISTFCK = 0
IPQ = 0
! Loop over symmetry:
do ISYM=1,NSYM
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NOO = NIO+NAO
  NO = NORB(ISYM)
  if (NOO /= 0) then
    do NP=1,NOO
      do NQ=1,NOO
        IPQ = IPQ+1
        FOCC(IPQ) = F(ISTFCK+NO*(NQ-1)+NP)
      end do
    end do
  end if
  ISTFCK = ISTFCK+NO**2
end do

IAD15 = IADR15(5)
call DDAFILE(JOBIPH,1,FOCC,IPQ,IAD15)

! Fock matrices added -- R L 921008.

! Start processing the all occupied indices Fock matrix

call mma_allocate(Scr1,no2m,Label='Scr1')
call mma_allocate(Scr2,no2m,Label='Scr2')
FockOcc(:) = Zero
#ifdef _DEBUGPRINT_
write(u6,*) 'nTot1=',nTot1
#endif

! Construct the occupied part of the Fock matrix in SO/AO basis.

ISTFCK = 1
jFock = 1
iCMo = 1
! A long loop over symmetry:
do ISYM=1,NSYM
  ! Hmm maybe you should check so I use correct nbas/norb
  NO = nOrb(isym)
  if (NO /= 0) then

    ! Transform to SO/AO basis.

#   ifdef _DEBUGPRINT_
    write(u6,*) 'iSym=',iSym
    call RecPrt('F(iStFck)',' ',F(iStFck),nOrb(iSym),nOrb(iSym))
    call RecPrt('CMO(iCMO)',' ',CMO(iCMO),nBas(iSym),nBas(iSym))
#   endif
    call DGEMM_('N','N', &
                nBas(iSym),nOrb(isym),nOrb(isym), &
                One,CMO(iCMo),nBas(iSym), &
                F(iStFck),nOrb(isym), &
                Zero,Scr1,nBas(iSym))
    call DGEMM_('N','T', &
                nBas(iSym),nBas(iSym),nOrb(isym), &
                One,Scr1,nBas(iSym), &
                CMO(iCMo),nBas(iSym), &
                Zero,Scr2,nBas(iSym))
    ij = jFock
    do iBas=1,nBas(iSym)
      do jBas=1,iBas-1
        kl = nBas(iSym)*(jBas-1)+iBas
        lk = nBas(iSym)*(iBas-1)+jBas
        FockOcc(ij) = Scr2(kl)+Scr2(lk)
        ij = ij+1
      end do
      kl = nBas(iSym)*(iBas-1)+iBas
      if (ij-jFock+1 > nTot1) then
        write(u6,*) ij,jFock,nTot1
        call Abend()
      end if
      FockOcc(ij) = Scr2(kl)
      ij = ij+1
    end do
  end if
# ifdef _DEBUGPRINT_
  call TriPrt('FAO',' ',FockOcc(jFock),nBas(iSym))
# endif
  jFock = jFock+nTri_Elem(nBas(iSym))
  iCMo = iCMo+nBas(iSym)**2
  ISTFCK = ISTFCK+NO**2
  ! End of long loop over symmetry
end do
call mma_deallocate(Scr2)
call mma_deallocate(Scr1)

! End of addition, R L 921008.

end subroutine FOCKOC
