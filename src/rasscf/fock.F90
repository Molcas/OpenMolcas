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

subroutine FOCK(F,BM,FI,FP,D,P,Q,FINT,IFINAL,CMO)
! RASSCF program version IBM-3090: SX section
!
!*****************************************************************************
!   BM    : output
!   F     : used only in this routine
!   Q     : used only in this routine
!   FP    : in Input it is the FA matrix in MO. in Output it is FI + FA in MO.
!   FI    : input. FI in MO basis
!   D     : input. Active 1RDM in MO basis
!   P     : input. Active 2RDM in MO basis
!   FINT  : input. Active two-electron integrals in MO basis.
!   IFINAL: input
!   CMO   : input
!*****************************************************************************
! Calculation of the MCSCF fock matrix F(eq.(7) in I.J.Q.C.S14,175)
! FP is the matrix FI+FA (FP is FA at entrance)
! F is stored as a symmetry blocked square matrix, by columns.
! Note that F contains all elements, also the zero elements
! occurring when the first index is secondary.
! F is used to construct the Brillouin elements and the first row
! of the super-CI Hamiltonian, while FP is used as the effective
! one-electron operator in the construction of the super-CI
! interaction matrix.
!
!      ********** IBM-3090 MOLCAS Release: 90 02 22 **********

use Index_Functions, only: iTri, nTri_Elem
use Fock_util_global, only: ALGO, DoCholesky
use rasscf_global, only: CBLB, CBLBM, E2act, ECAS, HalfQ1, IBLB, IBLBM, ISTORD, ISTORP, iSymBB, ixSym, iZROT, JBLB, JBLBM, KSDFT, &
                         NTOT3, via_DFT
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NASH, NBAS, NFRO, NISH, NORB, NSSH, NSYM
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: F(*), BM(*), FI(*), FP(*), D(*), P(*), Q(*), FINT(*), CMO(*)
integer(kind=iwp) :: iFinal
integer(kind=iwp) :: ipBM, ipFMCSCF, ipMOs, ipQs, iPrLev, ISTAV(8), ISTBM, ISTD, ISTFCK, ISTFP, ISTP, ISTSQ(8), ISTZ, iSym, IX, &
                     IX1, JSTF, N1, N2, NAO, NAS, NEO, NI, NIA, NIO, NIS, NM, NO, NO2, nOr, NP, NPQ, NQ, nSs, NT, NTM, NTT, NTU, &
                     NTV, NU, NUVX, NV, NVI, NVM
real(kind=wp) :: CASDFT_En, CSX, ECAS0, QNTM

IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering FOCK'

ISTSQ(1) = 0
ISTAV(1) = 0
do iSym=2,nSym
  ISTSQ(iSym) = ISTSQ(iSym-1)+nBas(iSym-1)**2
  ISTAV(iSym) = ISTAV(iSym-1)+nBas(iSym-1)*nAsh(iSym-1)
end do

!***********************************************************************
! add FI to FA to obtain FP
!***********************************************************************
FP(1:NTOT3) = FP(1:NTOT3)+FI(1:NTOT3)

!***********************************************************************
! Loop over symmetry blocks
!***********************************************************************
ISTFCK = 0
ISTFP = 0
ISTD = 0
ISTBM = 0
IX1 = 0
ISTZ = 0
E2act = Zero
HalfQ1 = Zero

do ISYM=1,NSYM
  IX = IX1+NFRO(ISYM)
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NEO = NSSH(ISYM)
  NIA = NIO+NAO
  NO = NORB(ISYM)
  NO2 = nTri_Elem(NO)
  CSX = Zero
  N1 = 0
  N2 = 0

  ! If no orbitals in this symmetry block go to next symm
  if (NO == 0) GO TO 90

  !************************
  ! clear F_gen matrix
  !************************
  F(ISTFCK+1:ISTFCK+NO**2) = Zero

  !*********************************************************************
  ! first index in F is inactive: F_gen is twice FP=FI+FA (Eq.10.8.27 MEST)
  !*********************************************************************
  if (NIO /= 0) then
    do NP=1,NO
      do NI=1,NIO
        F(ISTFCK+NO*(NP-1)+NI) = 2*FP(ISTFP+iTri(NP,NI))
      end do
    end do
  end if

  !*********************************************************************
  ! first index in F active
  !*********************************************************************
  if (NAO /= 0) then
    ISTP = ISTORP(ISYM)+1
    JSTF = ISTORD(ISYM)+1
    NUVX = (ISTORP(ISYM+1)-ISTORP(ISYM))/NAO

    if ((.not. DoCholesky) .or. (ALGO == 1)) then
      !*****************************************************************
      ! Compute the Q-matrix (Eq. (19) in IJQC S14 175 1980 or Eq. 10.8.31 MEST)
      ! Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
      ! P is packed in xy and pre-multiplied by 2 and reordered
      !*****************************************************************
      if (IPRLEV >= DEBUG) then
        write(u6,*) 'PUVX integrals in FOCK'
        call wrtmat(FINT(JSTF),NO,NUVX,NO,NUVX)
        write(u6,*) 'two-elec density mat OR DMAT*DMAT in FOCK'
        call wrtmat(P(ISTP),NAO,NUVX,NAO,NUVX)
      end if

      call DGEMM_('N','N',NO,NAO,NUVX,One,FINT(JSTF),NO,P(ISTP),NUVX,Zero,Q,NO)

    else if (ALGO == 2) then
      !*****************************************************************
      ! --- the Q-matrix has been already computed as Q(a,v)
      ! --- where a is an AO index and v is an active index
      ! --- Transform the 1st index to MOs (one symmetry at the time)
      ! --- Q(m,v) = C(a,m) * Q(a,v)
      !*****************************************************************
      ipQS = 1+ISTAV(iSym)
      ipMOs = 1+ISTSQ(iSym)+nBas(iSym)*nFro(iSym)

      call DGEMM_('T','N',nOrb(iSym),nAsh(iSym),nBas(iSym),One,CMO(ipMOs),nBas(iSym),Q(ipQS),nBas(iSym),Zero,Q(1),nOrb(iSym))
    else
      write(u6,*) 'FOCK: illegal Cholesky parameter ALGO= ',ALGO
      call abend()
    end if

    if (IPRLEV >= DEBUG) call recprt('Q-mat in fock.f',' ',Q(1),NO,NAO)

    !*******************************************************************
    ! active-active interaction energy term in the RASSCF energy: trace of Q
    !*******************************************************************
    ECAS0 = ECAS
    do NT=1,NAO
      NTT = (NT-1)*NO+NIO+NT
      ECAS = ECAS+Half*Q(NTT)
      HALFQ1 = HALFQ1+Half*Q(NTT)
      E2act = E2act+Half*Q(NTT)
    end do
    if (IPRLEV >= DEBUG) then
      write(u6,*) 'Two-electron contribution (Q term):',ECAS-ECAS0
      write(u6,*) 'ECAS aft adding Q in fock.f :',ECAS
    end if

    !*******************************************************************
    ! Continue... Fock matrix for first index active
    !*******************************************************************
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
        F(ISTFCK+NO*(NM-1)+NT+NIO) = QNTM
      end do
    end do

  end if
  !*********************************************************************
  !^ End if related to first index in F active
  !*********************************************************************

  !*********************************************************************
  !       The Brillouin matrix BM(pq)=F(qp)-F(pq)
  !*********************************************************************
  NPQ = ISTBM
  do NP=NIO+1,NO
    do NQ=1,NIA
      NPQ = NPQ+1
      BM(NPQ) = F(ISTFCK+NO*(NP-1)+NQ)-F(ISTFCK+NO*(NQ-1)+NP)

      ! Set zeroes to BM elements corresponding to rotations not allowed
      ! as controlled by the array IZROT.

      if ((NP <= NIA) .and. (NQ > NIO)) then
        NT = NP-NIO
        NU = NQ-NIO
        if (NT <= NU) then
          BM(NPQ) = Zero
        else
          NTU = ISTZ+nTri_Elem(NT-2)+NU
          if (IZROT(NTU) /= 0) BM(NPQ) = Zero
          if (IXSYM(IX+NP) /= IXSYM(IX+NQ)) BM(NPQ) = Zero
        end if
      end if

      ! check for largest Brillouin matrix element

      if (abs(BM(NPQ)) < abs(CSX)) GO TO 20
      if (IXSYM(IX+NP) /= IXSYM(IX+NQ)) GO TO 20
      CSX = BM(NPQ)
      N1 = NQ+NFRO(ISYM)
      N2 = NP+NFRO(ISYM)
20    continue
    end do
  end do

90 continue
  ISTFCK = ISTFCK+NO**2
  ISTFP = ISTFP+NO2
  ISTD = ISTD+nTri_Elem(NAO)
  ISTBM = ISTBM+(NIO+NAO)*(NAO+NEO)
  IX1 = IX1+NBAS(ISYM)
  ISTZ = ISTZ+nTri_Elem(NAO-1)
  CBLB(ISYM) = CSX
  IBLB(ISYM) = N1
  JBLB(ISYM) = N2
end do
!***********************************************************************
!^ End of long loop over symmetry
!***********************************************************************

if (iPrLev >= DEBUG) then
  CASDFT_En = Zero
  if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_En)
  write(u6,'(A,2F22.16)') ' RASSCF energy: ',ECAS+CASDFT_En,VIA_DFT
end if
if (iPrLev >= DEBUG) then
  write(u6,'(A)') ' MCSCF Fock-matrix in MO-basis'
  ipFMCSCF = 1
  do iSym=1,nSym
    nOr = nOrb(iSym)
    call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
    ipFMCSCF = ipFMCSCF+nOr*nOr
  end do
end if
if (iPrLev >= DEBUG) then
  write(u6,'(A)') ' Brillouin matrix'
  ipBM = 1
  do iSym=1,nSym
    nIs = nIsh(iSym)
    nAs = nAsh(iSym)
    nSs = nSsh(iSym)
    call RecPrt(' ',' ',BM(ipBM),(nAs+nIs),(nAs+nSs))
    ipBM = ipBM+(nAs+nIs)*(nAs+nSs)
  end do
end if

! Maximum BLB matrix element all symmetries

CBLBM = Zero
ISYMBB = 0
do ISYM=1,NSYM
  if (abs(CBLB(ISYM)) < abs(CBLBM)) GO TO 150
  CBLBM = CBLB(ISYM)
  IBLBM = IBLB(ISYM)
  JBLBM = JBLB(ISYM)
  ISYMBB = ISYM
150 continue
end do

!***********************************************************************
! Calculate Fock matrix for occupied orbitals.
!***********************************************************************
if (iFinal == 1) call FOCKOC(Q,F,CMO)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' >>> Exit Fock <<< '
  write(u6,*)
end if

end subroutine FOCK
