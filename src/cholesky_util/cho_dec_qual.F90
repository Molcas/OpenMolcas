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
!  Cho_Dec_Qual
!
!> @author F. Aquilante
!>
!> @details
!> The symmetric positive (semi-)definite matrix
!>
!> \f[ Q_{ab,cd} = Q_{ab,cd} - \sum_J L_{ab,J} L_{cd,J} \f]
!>
!> is evaluated and the diagonal elements are returned in \p QDiag.
!> Subsequently, \p Q is Cholesky decomposed to yield \p NumV vectors
!> \f$ K_{ab,I} \f$ as
!>
!> \f[ Q_{ab,cd} = \sum_I K_{ab,I} K_{cd,I} \f]
!>
!> @note
!> \p Q is destroyed and cannot be used after this routine.
!>
!> @param[in,out] Diag  Array of diagonal integrals stored as in the first reduced set
!> @param[in]     Lab   Array of Cholesky vectors (symmetry blocked) for the qualified shell pairs
!> @param[in,out] Q     Array (symmetry blocked) containing the matrix of the qualified integrals
!> @param[out]    Kab   Array symmetry blocked of Cholesky vectors resulting from the decomposition
!>                      of the qualified integral matrix
!> @param[out]    iD    Index array. On exit \p iD(k,jSym) contains the index of the qualified diagonals
!>                      from which the \p k -th vector of \p jSym was generated
!> @param[out]    NumV  Array of the number of vectors in each symmetry from the decomposition of the
!>                      qualified integral matrix
!> @param[out]    QDiag Array (symmetry blocked) containing the updated qualified diagonals
!***********************************************************************

subroutine Cho_Dec_Qual(Diag,Lab,Q,Kab,iD,NumV,QDiag)

use Cholesky, only: Cho_1Center, nQual, nSym, Span, ThrCom
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: Diag(*), Q(*)
real(kind=wp), intent(in) :: Lab(*)
real(kind=wp), intent(_OUT_) :: Kab(*), QDiag(*)
integer(kind=iwp), intent(_OUT_) :: iD(*), NumV(*)
integer(kind=iwp) :: ipD, ipQ, ipQD, ipSQK, ipSQL, iQ, irc, jSym, mQ, nVecG(8)
real(kind=wp) :: Dmax(8), Thr
logical(kind=iwp) :: Sync
character(len=*), parameter :: SecNam = 'Cho_Dec_Qual'

irc = 0

! Determine the max diagonal (qualified excluded)
! For 1-center decomposition, this is done later.
if (Cho_1Center) then
  Dmax(1:nSym) = Zero
else
  Sync = .false.
  call Cho_P_MaxDX(Diag,Sync,Dmax)
end if

! Set the number of previous Cholesky vectors (global)
call Cho_P_GetGV(nVecG,nSym)

ipD = 1
ipQ = 1
ipSQL = 1
ipSQK = 1
ipQD = 0
do jSym=1,nSym

  mQ = max(nQual(jSym),1)

  ! Do the subtraction  Q({ab}|{cd}) -= sum_J  L({ab},J) * L({cd},J)
  ! -----------------------------------------------------------------
  call DGEMM_('N','T',nQual(jSym),nQual(jSym),nVecG(jSym),-One,Lab(ipSQL),mQ,Lab(ipSQL),mQ,One,Q(ipQ),mQ)

  ! Extract diagonal of updated Q({ab}|{cd})
  ! ----------------------------------------
  do iQ=1,nQual(jSym)
    QDiag(ipQD+iQ) = Q(ipQ-1+nQual(jSym)*(iQ-1)+iQ)
  end do

  ! For 1-Center decomposition, find max. diagonal among qualified
  ! --------------------------------------------------------------
  if (Cho_1Center) then
    do iQ=1,nQual(jSym)
      Dmax(jSym) = max(Dmax(jSym),QDiag(ipQD+iQ))
    end do
  end if

  ! Cholesky decompose Q({ab}|{cd})
  ! Note: Q is destroyed in CD_InCore_p
  ! -----------------------------------
  Thr = max(ThrCom,Dmax(jSym)*Span)
  call CD_InCore_p(Q(ipQ),nQual(jSym),Kab(ipSQK),nQual(jSym),iD(ipD),NumV(jSym),Thr,irc)

  if (irc /= 0) then
    write(u6,*) SecNam,' non-zero rc on exit from CD_InCore_p: ',irc
    call Cho_Quit('Decomposition error in '//SecNam,104)
  end if

  ipSQL = ipSQL+nQual(jSym)*nVecG(jSym)
  ipSQK = ipSQK+nQual(jSym)**2
  ipQ = ipQ+nQual(jSym)**2
  ipD = ipD+nQual(jSym)
  ipQD = ipQD+nQual(jSym)

end do

return

end subroutine Cho_Dec_Qual
