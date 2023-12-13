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
! Copyright (C) Thomas Bondo Pedersen                                  *
!               2020,2021, Roland Lindh                                *
!***********************************************************************
!  Cho_X_CalculateGMat
!
!> @brief
!>   Calculate Cholesky \f$ G \f$ matrix (metric matrix)
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Calculate the metric matrix from Cholesky vectors (i.e. exact).
!>
!> \f[ G_{IJ} = \sum_{K=1}^{\min(I,J)} L_{IK} L_{JK} \f]
!>
!> The matrix is symmetric and stored on disk (file ``AVECXX``) in
!> triangular storage. The file is opened and closed in this
!> routine using routines ::DAName_MF_WA and ::DAClos, respectively
!> (i.e. \c FileName is opened as a word addressable multifile).
!> The calculation failed if \p irc is different from zero on exit.
!>
!> @note
!> This routine should *not* be used with DF.
!>
!> @param[out] irc Return code
!***********************************************************************

subroutine Cho_X_CalculateGMat(irc)

use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: iiBstR, InfVec, nnBstR, nSym, NumCho
use Cholesky_procedures, only: Cho_CGM_InfVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: i, iDisk, idRS2RS, iI, iJ, iLoc, iOpt, iRed, iRedC, iSym, j, K, kG_IJ, KK, KK1, KKK, kOffV, l_G, l_iRS2RS, &
                     l_Wrk, lUnit, mUSed, nVRead
real(kind=wp) :: V_J
logical(kind=iwp) :: isDF
character(len=6) :: FileName
integer(kind=iwp), pointer :: InfVcT(:,:,:)
integer(kind=iwp), allocatable :: iRS2RS(:), NVT(:)
real(kind=wp), allocatable :: G(:), Wrk(:)
#ifdef _DEBUGPRINT_
real(kind=wp), external :: ddot_
#endif

! Set return code.
! ----------------

irc = 0

! Refuse to perform calculation for DF.
! -------------------------------------

call DecideOnDF(isDF)
if (isDF) then
  irc = -1
  return
end if

! Scratch location in index arrays.
! ---------------------------------

iLoc = 3 ! do NOT change (used implicitly by reading routine)

! Get pointer to InfVec array for all vectors (needed for parallel
! runs) and the total number of vectors.
! ----------------------------------------------------------------

call mma_allocate(NVT,nSym,Label='NVT')
call Cho_CGM_InfVec(InfVcT,NVT,size(NVT))

! Copy rs1 to location 2.
! -----------------------

call Cho_X_RSCopy(irc,1,2)
if (irc /= 0) then
  irc = 1
  call Finish_this()
  return
end if

! Calculate triangular G matrix.
! G(IJ)=sum_K L(I,K)*L(J,K).
! ------------------------------

iRedC = -1
do iSym=1,nSym

  ! Open file.
  ! ----------

  write(FileName,'(A4,I2.2)') 'AVEC',iSym-1
  lUnit = 7
  call DAName_MF_WA(lUnit,FileName)
  iDisk = 0

  l_iRS2RS = nnBstR(iSym,1)
  call mma_allocate(iRS2RS,l_iRS2RS,Label='iRS2RS')
  iRS2RS(:) = 0
  l_G = nTri_Elem(NVT(iSym))
  call mma_allocate(G,l_G,Label='G')
  call mma_MaxDBLE(l_Wrk)
  call mma_allocate(Wrk,l_Wrk,Label='Wrk')
  Wrk(:) = Zero
  G(:) = Zero

  idRS2RS = -2
  KK1 = 1
  do while (KK1 <= NumCho(iSym))
    nVRead = 0
    mUsed = 0
    call Cho_X_VecRd(Wrk,size(Wrk),KK1,NumCho(iSym),iSym,nVRead,iRedC,mUsed)
    if (nVRead < 1) then
      irc = 2
      ! exit after deallocation
      call Finish_this()
      return
    end if
    kOffV = 0
    do KKK=0,nVRead-1
      KK = KK1+KKK
      iRed = InfVec(KK,2,iSym)
      if (iRedC /= iRed) then
        call Cho_X_SetRed(irc,iLoc,iRed)
        if (irc /= 0) then
          irc = 3
          ! exit after deallocation
          call Finish_this()
          return
        end if
        iRedC = iRed
      end if
      if (idRS2RS /= iRedC) then
        call Cho_RS2RS(iRS2RS,size(iRS2RS),2,iLoc,iRedC,iSym)
        idRS2RS = iRedC
      end if
      K = InfVec(KK,5,iSym)
      do J=K,NVT(iSym)
        iJ = iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
        V_J = Wrk(kOffV+iJ)
        do I=K,J
          iI = iRS2RS(InfVcT(I,1,iSym)-iiBstR(iSym,1))
          kG_IJ = iTri(I,J)
          G(kG_IJ) = G(kG_IJ)+Wrk(kOffV+iI)*V_J
        end do
      end do
      kOffV = kOffV+nnBstR(iSym,iLoc)
    end do
    KK1 = KK1+nVRead
  end do
  call Cho_GADGOP(G,size(G),'+')
  iOpt = 1
  call DDAFile(lUnit,iOpt,G,size(G),iDisk)
# ifdef _DEBUGPRINT_
  call TriPrt('G-matix',' ',G,NVT(iSym))
  write(u6,'(A,I2,A,ES16.7)') 'G matrix, sym.',iSym,': Norm = ',sqrt(dDot_(size(G),G,1,G,1))
# endif
  call mma_deallocate(Wrk)
  call mma_deallocate(G)
  call mma_deallocate(iRS2RS)

  ! Close file
  ! ----------
  call DAClos(lUnit)
end do

call Finish_this()

contains

subroutine Finish_this()

  call mma_deallocate(NVT)

end subroutine Finish_this

end subroutine Cho_X_CalculateGMat
