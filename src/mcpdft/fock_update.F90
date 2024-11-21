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

!> @brief construct PDFT generalized fock matrix
SUBROUTINE FOCK_update(F,FI,FP,D,P,Q,FINT,CMO)
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
  use definitions,only:iwp,wp,u6
  use constants,only:zero,one,two
  use stdalloc,only:mma_allocate,mma_deallocate
  use printlevel,only:debug
  use mcpdft_output,only:iPrLoc
  use rasscf_global,only:nTot3,nTot4,ISTORP,iTri,ISTORD
  use general_data,only:nsym,nbas,nash,nish,norb

  implicit none
  real(kind=wp) :: FI(*),FP(*),D(*),P(*),Q(*),FINT(*),F(*),CMO(*)
  integer(kind=iwp) :: ISTSQ(8),ISTAV(8)

  real(kind=wp),allocatable :: TF(:)

  real(kind=wp) :: QNTM
  integer(kind=iwp) :: ipFMCSCF,ISTD,ISTFCK,ISTFP,ISTP,iprlev
  integer(kind=iwp) :: iSym,JSTF,N1,N2,NAO,NI,NIO,NM,NO,NO2
  integer(kind=iwp) :: NOR,NP,NT,NTM,NTV,NUVX,NV,NVI,NVM

  IPRLEV = IPRLOC(4)

  Call mma_allocate(TF,NTOT4,Label='TF')
  TF(:) = zero

  ISTSQ(1) = 0
  ISTAV(1) = 0
  DO iSym = 2,nSym
    ISTSQ(iSym) = ISTSQ(iSym-1)+nBas(iSym-1)**2
    ISTAV(iSym) = ISTAV(iSym-1)+nBas(iSym-1)*nAsh(iSym-1)
  EndDo

  ! add FI to FA to obtain FP
  FP(:ntot3) = FP(:ntot3)+FI(:ntot3)

  ISTFCK = 0
  ISTFP = 0
  ISTD = 0
  ! A long loop over symmetry
  DO ISYM = 1,NSYM
    NIO = NISH(ISYM)
    NAO = NASH(ISYM)
    NO = NORB(ISYM)
    NO2 = (NO**2+NO)/2
    IF(NO == 0) then
      cycle
    endif

    ! First index in F is inactive
    IF(NIO /= 0) THEN
      DO NP = 1,NO
        DO NI = 1,NIO
          N1 = MAX(NP,NI)
          N2 = MIN(NP,NI)
          TF(ISTFCK+NO*(NP-1)+NI) = two*FP(ISTFP+(N1**2-N1)/2+N2)
        ENDDO
      ENDDO
    ENDIF

    ! first index in F active
    IF(NAO /= 0) THEN

      ISTP = ISTORP(ISYM)+1
      JSTF = ISTORD(ISYM)+1
      NUVX = (ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

      ! first compute the Q-matrix (equation (19))
      !
      ! Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
      !
      ! P is packed in xy and pre-multiplied by 2
      !                   and reordered
      CALL DGEMM_('N','N', &
                  NO,NAO,NUVX, &
                  one,FINT(JSTF),NO, &
                  P(ISTP),NUVX, &
                  zero,Q,NO)

      ! Now Q should contain the additional 2-electron part of the fock matrix
      ! for mcpdft, for the active region, at least.
      !
      ! We should also have contributions from terms like FI and FA, too.
      ! FA takes care of the 1-RDM/2e- integral terms?
      ! FI takes care of the one-body hamiltonian and the occ/occ and occ/act
      ! contributions.

      ! Fock matrix
      NTM = 0
      DO NT = 1,NAO
        DO NM = 1,NO
          NTM = NTM+1
          QNTM = Q(NTM)
          DO NV = 1,NAO
            NVI = NV+NIO
            NTV = ITRI(MAX(NT,NV))+MIN(NT,NV)+ISTD
            NVM = ITRI(MAX(NVI,NM))+MIN(NVI,NM)+ISTFP
            QNTM = QNTM+D(NTV)*FI(NVM)
          ENDDO
          TF(ISTFCK+NO*(NM-1)+NT+NIO) = QNTM
        ENDDO
      ENDDO
    ENDIF

    ! End of long loop over symmetry
    ISTFCK = ISTFCK+NO**2
    ISTFP = ISTFP+NO2
    ISTD = ISTD+(NAO**2+NAO)/2
  ENDDO

  F(:ntot4) = F(:ntot4)+TF(:ntot4)

  ! Calculate Fock matrix for occupied orbitals.

  ! I am going to add the Fock matrix temporarily to the Runfile.  I don't
  ! want to construct it again in MCLR in the case of gradients.
  If(iPrLev >= DEBUG) then
    Write(u6,'(A)') ' MCSCF Fock-matrix in MO-basis'
    ipFMCSCF = 1
    Do iSym = 1,nSym
      nOr = nOrb(iSym)
      Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
      ipFMCSCF = ipFMCSCF+nOr*nOr
    EndDo
  EndIf

  ! This will set FockOcc in wadr
  CALL FOCKOC(Q,F,CMO)

  Call mma_deallocate(TF)

ENDSUBROUTINE FOCK_update
