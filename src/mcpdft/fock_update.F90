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
!  fock_update
!
!> @brief construct PDFT generalized fock matrix
!***********************************************************************

subroutine FOCK_update(F,FI,FP,D,P,Q,NQ,FINT,CMO)
! This subroutine is supposed to add the dft portions of the mcpdft fock
! matrix to the Fock matrix pieces that have already been built for the
! CASSCF portion.

! Calculation of the MCSCF fock matrix F(eq.(7) in I.J.Q.C.S14,175)
! FP is the matrix FI+FA (FP is FA at entrance)
! F is stored as a symmetry blocked square matrix, by columns.
! Note that F contains all elements, also the zero elements
! occurring when the first index is secondary.
!
! FockOcc is further constructed as well

use Index_Functions, only: iTri, nTri_Elem
use PrintLevel, only: DEBUG
use mcpdft_output, only: iPrLoc
use rasscf_global, only: ISTORD, ISTORP, nacpar, nFint, nTot3, nTot4
use general_data, only: nash, nbas, nish, norb, nsym, nTot1, nTot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: iwp, wp, u6

implicit none
real(kind=wp), intent(inout) :: F(nTot4), FP(nTot1)
real(kind=wp), intent(in) :: FI(nTot1), D(nacpar), P(*), FINT(nFint), CMO(nTot2)
integer(kind=iwp), intent(in) :: NQ
real(kind=wp), intent(out) :: Q(NQ)
integer(kind=iwp) :: ipFMCSCF, iprlev, ISTAV(8), ISTD, ISTFCK, ISTFP, ISTP, ISTSQ(8), iSym, JSTF, NAO, NI, NIO, NM, NO, NO2, NOR, &
                     NP, NT, NTM, NTV, NUVX, NV, NVI, NVM
real(kind=wp) :: QNTM
real(kind=wp), allocatable :: TF(:)

IPRLEV = IPRLOC(4)

call mma_allocate(TF,NTOT4,Label='TF')
TF(:) = Zero

ISTSQ(1) = 0
ISTAV(1) = 0
do iSym=2,nSym
  ISTSQ(iSym) = ISTSQ(iSym-1)+nBas(iSym-1)**2
  ISTAV(iSym) = ISTAV(iSym-1)+nBas(iSym-1)*nAsh(iSym-1)
end do

! add FI to FA to obtain FP
FP(:ntot3) = FP(:ntot3)+FI(:ntot3)

ISTFCK = 0
ISTFP = 0
ISTD = 0
! A long loop over symmetry
do ISYM=1,NSYM
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NO = NORB(ISYM)
  NO2 = nTri_Elem(NO)
  if (NO == 0) cycle

  ! First index in F is inactive
  if (NIO /= 0) then
    do NP=1,NO
      do NI=1,NIO
        TF(ISTFCK+NO*(NP-1)+NI) = Two*FP(ISTFP+iTri(NP,NI))
      end do
    end do
  end if

  ! first index in F active
  if (NAO /= 0) then

    ISTP = ISTORP(ISYM)+1
    JSTF = ISTORD(ISYM)+1
    NUVX = (ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

    ! first compute the Q-matrix (equation (19))
    !
    ! Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
    !
    ! P is packed in xy and pre-multiplied by 2
    !                   and reordered
    call DGEMM_('N','N', &
                NO,NAO,NUVX, &
                One,FINT(JSTF),NO, &
                P(ISTP),NUVX, &
                Zero,Q,NO)

    ! Now Q should contain the additional 2-electron part of the fock matrix
    ! for mcpdft, for the active region, at least.
    !
    ! We should also have contributions from terms like FI and FA, too.
    ! FA takes care of the 1-RDM/2e- integral terms?
    ! FI takes care of the one-body hamiltonian and the occ/occ and occ/act
    ! contributions.

    ! Fock matrix
    NTM = 0
    do NT=1,NAO
      do NM=1,NO
        NTM = NTM+1
        QNTM = Q(NTM)
        do NV=1,NAO
          NVI = NV+NIO
          NTV = iTri(NT,NV)+ISTD
          NVM = iTri(NVI,NM)+ISTFP
          QNTM = QNTM+D(NTV)*FI(NVM)
        end do
        TF(ISTFCK+NO*(NM-1)+NT+NIO) = QNTM
      end do
    end do
  end if

  ! End of long loop over symmetry
  ISTFCK = ISTFCK+NO**2
  ISTFP = ISTFP+NO2
  ISTD = ISTD+nTri_Elem(NAO)
end do

F(:ntot4) = F(:ntot4)+TF(:ntot4)

! Calculate Fock matrix for occupied orbitals.

! I am going to add the Fock matrix temporarily to the Runfile.
! I don't want to construct it again in MCLR in the case of gradients.
if (iPrLev >= DEBUG) then
  write(u6,'(A)') ' MCSCF Fock-matrix in MO-basis'
  ipFMCSCF = 1
  do iSym=1,nSym
    nOr = nOrb(iSym)
    call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
    ipFMCSCF = ipFMCSCF+nOr*nOr
  end do
end if

! This will set FockOcc in wadr
call FOCKOC(Q,F,CMO)

call mma_deallocate(TF)

end subroutine FOCK_update
