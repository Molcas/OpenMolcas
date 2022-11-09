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
!*****************************************************************
!
! This subroutine substitutes the previous version of the RdOrd
! to allow the generation of two-electron integrals from a set
! of Cholesky vectors
!
!     Data declarations:
!
!        First  : logical variable for identifying the 1st call
!                 to this subroutine
!
!   F. Aquilante
!
!*****************************************************************
!***********************************************************************
!  RdOrd
!
!> @brief
!>   Driver for the actual ::RdOrd_ routines and the Cholesky vector integral generator
!> @author F. Aquilante
!>
!> @note
!> The integrals are returned in the order \f$ (lk|ji) \f$.
!>
!> @param[out] rc   return code
!> @param[in]  iOpt option code (= ``1``: start reading at first shell distribution
!>                  in the given product symmetry. = ``2``: continue reading)
!> @param[in]  iSym irred. representation of "first" symmetry label
!> @param[in]  jSym irred. representation of "second" symmetry label
!> @param[in]  kSym irred. representation of "third" symmetry label
!> @param[in]  lSym irred. representation of "fourth" symmetry label
!> @param[out] Buf  contains on output the integrals
!> @param[in]  lBuf length of the integral buffer
!> @param[out] nMat number of submatrices read in
!***********************************************************************

subroutine RdOrd(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc, nMat
integer(kind=iwp), intent(in) :: iOpt, iSym, jSym, kSym, lSym, lBuf
real(kind=wp), intent(_OUT_) :: Buf(*)
logical(kind=iwp) :: First = .true., DoCholesky = .false.

if (First) then
  call DecideOnCholesky(DoCholesky)
  if (DoCholesky) call INIT_GETINT(rc) ! initialize information
  First = .false.
end if

if (DoCholesky) then
  call Get_Int(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
else
  call RdOrd_(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
end if

! Debug printing
!
!call get_iscalar('nSym',nsym)
!call get_iarray('nbas',nbas,nsym)
!if (kSym == lSym) then
!  Nkl = nTri_Elem(nBas(kSym))
!else
!  Nkl = nBas(kSym)*nBas(lSym)
!end if
!write(u6,*)
!write(u6,*)'Integrals: symm(ij|kl), nMat ',iSym,jSym,kSym,lSym,nMat
!write(u6,*)'Nbask,Nbasl,Nkl= ',nBas(kSym),nBas(lSym),Nkl
!xnrm = sqrt(ddot_(Nkl*nMat,Buf,1,Buf,1))
!write(u6,*)'norm= ',xnrm
!call cho_output(Buf,1,nMat,1,Nkl,nMat,Nkl,1,u6)
!call xflush(u6)

return

end subroutine RdOrd
