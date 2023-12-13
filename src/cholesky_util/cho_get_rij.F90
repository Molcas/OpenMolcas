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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!  CHO_get_Rij
!
!> @author F. Aquilante
!>
!> @details
!> Computes the \f$ R \f$ matrix used for the
!> maximization of Edmiston--Ruedenberg functional
!>
!> \f[ \mathit{ER} = \sum_i (ii|ii) = \mathrm{Tr}(R) \f]
!>
!> for a given set of MOs.
!>
!> \f$ R \f$ is defined from the two-electron integrals
!> computed from the MO-transformed Cholesky vectors
!>
!> \f[ R_{ij} = (ij|jj) = \sum_K L_{ij,K} L_{jj,K} \f]
!>
!> and the condition for the maximization of the ER-functional
!> is given by
!>
!> \f[ \mathrm{grad}(\mathit{ER})_{ij} = 4(R_{ij} - R_{ji}) = 0 \quad (\forall i,j) \f]
!>
!> @note
!> Requires initialization of the Cholesky information.
!>
!> @param[out] irc     Return code
!> @param[in]  MO      type DSBA_Type of block of the MO matrix, stored as \p C(k,a)
!> @param[in]  nOcc    Number of orbitals to be localized in each symmetry
!> @param[out] Rij     \p nOcc &times; \p nOcc symmetry blocked matrix \f$  R_{ij} = (ij|jj) \f$
!> @param[in]  timings Switch on/off timings printout
!***********************************************************************

subroutine CHO_get_Rij(irc,MO,nOcc,Rij,timings)

use Cholesky, only: InfVec, nBas, nDimRS, nSym, NumCho
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
type(DSBA_Type), intent(in) :: MO(1)
integer(kind=iwp), intent(in) :: nOcc(*)
real(kind=wp), intent(_OUT_) :: Rij(*)
logical(kind=iwp), intent(in) :: timings
integer(kind=iwp) :: iBatch, iE, iLoc, iOcc(8), iOcs(8), IREDC, iS, iSkip(8), iSwap, IVEC2, iVrs, JNUM, jpR, JRED, JRED_, JRED1, &
                     JRED2, JSYM, jVEC, kMOs, kS, kSym, lj, LREAD, LWORK, Maj, Mneed, Mocc, MUSED, n1, nBatch, nMOs, nOcs, nRS, &
                     NUMV, nVec, nVrs
real(kind=wp) :: TCI1, TCI2, TCR1, TCR2, TCT1, TCT2, tintg(2), tmotr(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, &
                 tread(2), TWI1, TWI2, TWR1, TWR2, TWT1, TWT2
type(SBA_Type) :: Laq(1)
real(kind=wp), allocatable, target :: Lab(:)
real(kind=wp), pointer :: pLab(:,:,:), pLjj(:)
logical(kind=iwp), parameter :: DoRead = .false.
character(len=*), parameter :: SECNAM = 'CHO_get_Rij'

IREDC = -1

JSYM = 1
if (NumCho(JSYM) < 1) then
  write(u6,*) SECNAM//': No total symmetric vectors present'
  irc = 77
  return
end if

do kS=1,nSym
  if ((nOcc(kS) > nBas(kS)) .or. (nOcc(kS) < 0)) then
    write(u6,*) SECNAM//': Wrong nOcc in symmetry ',kS
    write(u6,*) 'nOcc(',kS,')= ',nOcc(kS)
    write(u6,*) 'nBas(',kS,')= ',nBas(kS)
    irc = 79
    return
  end if
end do

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = Zero  !time for reading the vectors
tmotr(:) = Zero  !time for MO transformation of the vectors
tintg(:) = Zero  !time for computing the functional

! compute some offsets and other quantities
iOcc(1) = 0
iOcs(1) = 0
Mocc = nOcc(1)
do kSym=2,nSym
  iOcc(kSym) = iOcc(kSym-1)+nOcc(kSym-1)
  iOcs(kSym) = iOcs(kSym-1)+nOcc(kSym-1)**2
  Mocc = max(Mocc,nOcc(kSym))
end do

nOcs = iOcs(nSym)+nOcc(nSym)**2

Rij(1:nOcs) = Zero

do kS=1,nSym
  iSkip(kS) = min(nOcc(kS),1) ! initialize skipping flags
end do

! Memory need for 1 of the half-transformed vectors : L(aj)
Maj = sum(nBas(1:nSym)*nOcc(1:nSym))

iLoc = 3 ! use scratch location in reduced index arrays

JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

do JRED=JRED1,JRED2

  ! Memory management section -----------------------------

  call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

  if (nVrs == 0) cycle

  if (nVrs < 0) then
    write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs < 0. STOP!!'
    call abend()
  end if

  call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
  if (irc /= 0) then
    write(u6,*) SECNAM//': cho_X_setred non-zero return code. rc= ',irc
    call abend()
  end if

  IREDC = JRED

  nRS = nDimRS(JSYM,JRED)

  Mneed = max(nRS,Mocc**2+1) ! mem. for Lab and (Lij + Ljj)

  call mma_maxDBLE(LWORK)

  nVec = min(LWORK/(Maj+Mneed),nVrs)

  if (nVec < 1) then
    write(u6,*) SECNAM//': Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',Maj+Mneed
    irc = 33
    call Abend()
    nBatch = -9999  ! dummy assignment
  end if

  LREAD = nRS*nVec

  call mma_allocate(Lab,Mneed*nVec,Label='Lab')

  ! BATCH over the vectors in JSYM=1 ----------------------------

  nBatch = (nVrs-1)/nVec+1

  do iBatch=1,nBatch

    if (iBatch == nBatch) then
      JNUM = nVrs-nVec*(nBatch-1)
    else
      JNUM = nVec
    end if
    iSwap = 2  ! LiK,b are returned
    call Allocate_DT(Laq(1),nOcc,nBas,JNUM,JSYM,nSym,iSwap)

    JVEC = nVec*(iBatch-1)+iVrs
    IVEC2 = JVEC-1+JNUM

    call CWTIME(TCR1,TWR1)

    JRED_ = JRED
    call CHO_VECRD(Lab,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED_,MUSED)

    if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
      irc = 77
      return
    end if

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

    ! First half-transformation of the vectors : Lab,J --> LiJ,b
    !-----------------------------------------------------------

    kMOs = 1
    nMOs = 1

    call CWTIME(TCT1,TWT1)

    call CHO_X_getVtra(irc,Lab,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,MO(1),Laq(1),DoRead)

    if (irc /= 0) return

    call CWTIME(TCT2,TWT2)
    tmotr(1) = tmotr(1)+(TCT2-TCT1)
    tmotr(2) = tmotr(2)+(TWT2-TWT1)

    do kSym=1,nSym

      n1 = nOcc(kSym)
      iS = 1
      iE = n1*JNUM*n1
      pLab(1:n1,1:JNUM,1:n1) => Lab(iS:iE)
      iS = iE+1
      iE = iE+JNUM
      pLjj(1:JNUM) => Lab(iS:iE)

      if (iSkip(kSym) /= 0) then

        call CWTIME(TCT1,TWT1)
        !---------------------------------------------------------------
        ! Second half-transformation  L(iK,j) = sum_b  L(iK,b) * C(j,b)
        !---------------------------------------------------------------

        call DGEMM_('N','T',nOcc(kSym)*JNUM,nOcc(kSym),nBas(kSym),One,Laq(1)%SB(kSym)%A3,nOcc(kSym)*JNUM,MO(1)%SB(kSym)%A2, &
                    nOcc(kSym),Zero,pLab,nOcc(kSym)*JNUM)

        call CWTIME(TCT2,TWT2)
        tmotr(1) = tmotr(1)+(TCT2-TCT1)
        tmotr(2) = tmotr(2)+(TWT2-TWT1)

        call CWTIME(TCI1,TWI1)

        do lj=1,nOcc(kSym)

          pLjj(1:JNUM) = pLab(lj,1:JNUM,lj)

          !-------------------------------------------------------------
          ! Compute   R(i,j) = sum_K  L(i,K)[j] * L(K)[j]
          !-------------------------------------------------------------

          jpR = iOcs(kSym)+nOcc(kSym)*(lj-1)+1

          call DGEMV_('N',nOcc(kSym),JNUM,One,pLab(:,:,lj),nOcc(kSym),pLjj,1,One,Rij(jpR),1)

        end do

        call CWTIME(TCI2,TWI2)
        tintg(1) = tintg(1)+(TCI2-TCI1)
        tintg(2) = tintg(2)+(TWI2-TWI1)

      end if
      nullify(pLjj)
      nullify(pLab)

    end do

    call Deallocate_DT(Laq(1))
  end do  !end batch loop

  ! free memory
  call mma_deallocate(Lab)

end do   ! loop over red sets

! sync Rij

call GAdGOp(Rij,nOcs,'+')

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  write(u6,*)
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*) 'Timing from ',SECNAM,'           CPU      WALL  '
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'TRANSFORM VECTORS                         ',tmotr(1),tmotr(2)
  write(u6,'(2x,A26,2f10.2)') 'COMPUTE Rij = (ij|jj)                     ',tintg(1),tintg(2)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

irc = 0

return

end subroutine CHO_get_Rij
