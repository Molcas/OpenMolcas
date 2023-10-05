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

subroutine CHO_X_GET_PARDIAG(jSym,iSO_ab)
!***********************************************************************
! Returns an array of the "a" and "b" indices that give rise to the
! parent diagonal from which a given J-index has been originated
! by the (molecular) Cholesky decomposition procedure
!
!   ia=iSO_ab(1,numcho(jSym))  contains the index of the basis "a"
!                              within its symm. block.
!                              (Note: there is an integer function
!                              "cho_isao(ia)" that returns the
!                              irrep to which "a" belongs)
!
!   "Location 2" of iSO_ab returns the same kind of info for "b"
!
!
! Author: F. Aquilante
!
!***********************************************************************

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Cholesky, only: iiBstR, InfVec, iRS2F, NumCho
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: jSym
integer(kind=iwp), intent(out) :: iSO_ab(2,*)
integer(kind=iwp) :: iOff, jRab, jv, kRab
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iRank, kv, nTot
integer(kind=iwp), allocatable :: ip_List(:), ip_Aux(:,:)
#endif

iOff = 0

#ifdef _MOLCAS_MPP_
call mma_allocate(ip_List,[0,nProcs-1],Label='ip_List')
ip_List(:) = 0
ip_List(MyRank) = NumCho(jSym)
call Cho_GAIGOP(ip_List,nProcs,'+')
nTot = 0
do iRank=1,nProcs
  nTot = nTot+ip_List(iRank-1)
  if (iRank == MyRank) iOff = nTot
end do
call mma_deallocate(ip_List)
#endif

do jv=1,NumCho(jSym)

  jRab = InfVec(jv,1,jSym) ! addr. 1st red set within jSym

  kRab = iiBstr(jSym,1)+jRab ! global addr. 1st red set

  iSO_ab(1,jv+iOff) = iRS2F(1,kRab) !global addr. of a in SO list
  iSO_ab(2,jv+iOff) = iRS2F(2,kRab) !same for b, always (a >= b)
end do

#ifdef _MOLCAS_MPP_
call mma_allocate(ip_Aux,2,nTot,Label='ip_Aux')
ip_Aux(:,:) = 0
do jv=1,NumCho(jSym)
  kv = InfVec(jv,5,jSym)
  ip_Aux(1,kv) = iSO_ab(1,jv+iOff)
  ip_Aux(2,kv) = iSO_ab(2,jv+iOff)
end do
iSO_ab(:,1:nTot) = ip_Aux(:,1:nTot)
call mma_deallocate(ip_Aux)
call Cho_GAIGOP(iSO_ab,2*nTot,'+')
#endif

return

end subroutine CHO_X_GET_PARDIAG
