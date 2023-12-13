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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_CheckDiagFromZ(irc,NVT,l_NVT,nBlock,l_nBlock,nV,l_nV1,l_nV2,iV1,l_iV11,l_iV12,ip_Z,l_Z1,l_Z2,Z,l_Z,Report)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Check integral diagonal from Z vectors:
!
!          (J|J) = sum[K=1,J] Z(J,K)*Z(J,K)
!
!          Return codes
!          irc=0: all is fine
!          irc<0: too negative diagonals encountered, but otherwise
!                 calculation seems converged
!          irc>0: calculation not converged

use Index_Functions, only: iTri
use Cholesky, only: LuPri, nnBstRT, nSym, ThrCom, ThrNeg, TOONEG, WARNEG
use Cholesky_procedures, only: Cho_X_GetIP_InfVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_NVT, NVT(l_NVT), l_nBlock, nBlock(l_nBlock), l_nV1, l_nV2, nV(l_NV1,l_NV2), l_iV11, l_iV12, &
                                 iV1(l_iV11,l_iV12), l_Z1, l_Z2, ip_Z(l_Z1,l_Z2), l_Z
real(kind=wp), intent(in) :: Z(l_Z)
logical(kind=iwp), intent(in) :: Report
integer(kind=iwp) :: iD, iSym, J_inBlock, jBlock, J, K_inBlock, kblock, kOffZ, n1, n2, n3, n4, n5, nTot
real(kind=wp) :: Damax, Damin, Dmax, Dmin
integer(kind=iwp), pointer :: InfVct(:,:,:)
real(kind=wp), allocatable :: IntDia(:)
character(len=*), parameter :: SecNam = 'Cho_CheckDiagFromZ'

! Get pointer to global InfVec array
call Cho_X_getIP_InfVec(InfVcT)

! Allocate memory for exact integral diagonal
call mma_allocate(IntDia,nnBstRT(1),Label='IntDia')

! Read diagonal
call Cho_IODiag(IntDia,2)

! Subtract Z vector contributions
do iSym=1,nSym
  do kBlock=1,nBlock(iSym)
    do K_inBlock=1,nV(kBlock,iSym)
      kOffZ = ip_Z(iTri(kBlock,kBlock),iSym)-1
      do J_inBlock=K_inBlock,nV(kBlock,iSym)
        J = iV1(kBlock,iSym)+J_inBlock-1
        iD = InfVcT(J,1,iSym)
        IntDia(iD) = IntDia(iD)-Z(kOffZ+iTri(J_inBlock,K_inBlock))**2
      end do
    end do
    do jBlock=kBlock+1,nBlock(iSym)
      do K_inBlock=1,nV(kBlock,iSym)
        kOffZ = ip_Z(iTri(jBlock,kBlock),iSym)-1+nV(jBlock,iSym)*(K_inBlock-1)
        do J_inBlock=1,nV(jBlock,iSym)
          J = iV1(jBlock,iSym)+J_inBlock-1
          iD = InfVcT(J,1,iSym)
          IntDia(iD) = IntDia(iD)-Z(kOffZ+J_inBlock)**2
        end do
      end do
    end do
  end do
end do

! Total number of vectors
! ...should equal number of converged diagonals
nTot = sum(NVT(1:nSym))

! Count diagonals smaller than threshold
n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0
Dmax = -9.0e9_wp
Damax = Zero
Dmin = 9.0e9_wp
Damin = 9.0e9_wp
do iSym=1,nSym
  do J=1,NVT(iSym)
    iD = InfVcT(J,1,iSym)
    Dmax = max(Dmax,IntDia(iD))
    Damax = max(Damax,abs(IntDia(iD)))
    Dmin = min(Dmin,IntDia(iD))
    Damin = min(Damin,abs(IntDia(iD)))
    if (IntDia(iD) <= ThrCom) n1 = n1+1
    if (IntDia(iD) < Zero) n2 = n2+1
    if (IntDia(iD) < ThrNeg) n3 = n3+1
    if (IntDia(iD) < WarNeg) n4 = n4+1
    if (IntDia(iD) < TooNeg) n5 = n5+1
  end do
end do

! Write a report if requested
if (Report) then
  call Cho_Head(SecNam//': Report on (J|J) Diagonal from Z','=',80,LuPri)
  write(LuPri,'(/,A,I8)') 'Total dimension of diagonal............',nnBstRT(1)
  write(LuPri,'(A,I8)') 'Number of Cholesky vectors.............',nTot
  write(LuPri,'(A,I8)') 'Converged diagonals....................',n1
  write(LuPri,'(A,I8)') 'Unconverged diagonals..................',nTot-n1
  write(LuPri,'(A,I8)') 'Negative diagonals.....................',n2
  write(LuPri,'(A,I8)') 'Neg. diag. that would be zeroed........',n3
  write(LuPri,'(A,I8)') 'Neg. diag. that would cause warning....',n4
  write(LuPri,'(A,I8)') 'Neg. diag. that would cause crash......',n5
  write(LuPri,'(A,ES15.6)') 'Max diagonal...........................',Dmax
  write(LuPri,'(A,ES15.6)') 'Min diagonal...........................',Dmin
  write(LuPri,'(A,ES15.6)') 'Max abs diagonal.......................',Damax
  write(LuPri,'(A,ES15.6)') 'Min abs diagonal.......................',Damin
end if

! Deallocation
call mma_deallocate(IntDia)

! Set return code and return
if (n1 == nTot) then
  if (n5 /= 0) then
    irc = -10
  else
    irc = 0
  end if
else
  irc = 10
end if

end subroutine Cho_CheckDiagFromZ
