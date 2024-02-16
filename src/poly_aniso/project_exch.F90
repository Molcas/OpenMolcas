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

subroutine project_exch(N1,N2,S1,S2,M1,M2,E1,E2,HEXCH,Jpar,Jc)
! this function determines the local pseudospins and rotates the hamiltonian
! to the local pseudospin basis

use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, u6

implicit none
integer N1, N2
real(kind=8) :: E1(N1), E2(N2) ! spin-orbit energies on each site
complex(kind=8) :: S1(3,N1,N1), S2(3,N2,N2) ! spin matrices on each site
complex(kind=8) :: M1(3,N1,N1), M2(3,N2,N2) ! magnetic moment matrices on each site
complex(kind=8) :: HEXCH(N1,N1,N2,N2) ! exchange hamiltonian
complex(kind=8) :: HEXCH2(N1,N1,N2,N2) ! exchange hamiltonian
complex(kind=8) :: HEXCH3(N1,N1,N2,N2) ! exchange hamiltonian
!--  local variables --
integer :: ns1, ns2, is1, is2, iprint !, i1, j1, i2, j2
!integer :: k1, k2, q1, q2, js1, js2, ms1, ms2
real(kind=8) :: Ethr, gtens(2,3), maxes(2,3,3), Jc(3,3)
complex(kind=8) :: Z1(N1,N1), Z2(N2,N2)
complex(kind=8) :: SR1(3,N1,N1), MR1(3,N1,N1)
complex(kind=8) :: SR2(3,N2,N2), MR2(3,N2,N2)
complex(kind=8) :: TMP(N1,N1)
!complex(kind=8) :: DIP_O1(N1,N1)
!complex(kind=8) :: DIP_W1(N1,N1)
!complex(kind=8) :: DIP_O2(N2,N2)
!complex(kind=8) :: DIP_W2(N2,N2)
!complex(kind=8) :: SP_MOW1,SP_MOW2
!complex(kind=8) :: QMAT(N1,N1,N2,N2) !, trace
!real(kind=8) :: WCG ! Clebsch-Gordan Coefficients
!logical DBG
!external WCG
complex(kind=8) :: Jpar(N1-1,-N1+1:N1-1,N2-1,-N2+1:N2-1)

!DBG = .false.

!write(u6,'(A)') 'J parameters in the initial ab intio basis:'
!call JKQPar(N1,N2,HEXCH,Jpar)
!call tensor2cart(1,1,Jpar(1,-1:1,1,-1:1),Jc)

! determine the pseudospin on each site (Z1 and Z2):
! threshold for determination of the local pseudospin main anisotropy axis
Ethr = 0.2_wp
ns1 = 0
ns2 = 0
do is1=1,N1
  if (E1(is1) < Ethr) then
    ns1 = ns1+1
  end if
end do
do is1=1,N2
  if (E2(is1) < Ethr) then
    ns2 = ns2+1
  end if
end do
write(u6,'(A,i3)') 'size of local pseudospin, site 1  =',ns1
write(u6,'(A,i3)') 'size of local pseudospin, site 2  =',ns2
gtens(:,:) = Zero
maxes(:,:,:) = Zero
call atens(M1(1:3,1:ns1,1:ns1),ns1,gtens(1,:),maxes(1,:,:),2)
call atens(M2(1:3,1:ns2,1:ns2),ns2,gtens(2,:),maxes(2,:,:),2)
! rotate the magnetic moment to the coordinate system of main magnetic axes on Ln
call rotmom2(S1,N1,maxes(1,:,:),SR1)
call rotmom2(M1,N1,maxes(1,:,:),MR1)
call rotmom2(S2,N2,maxes(2,:,:),SR2)
call rotmom2(M2,N2,maxes(2,:,:),MR2)
iprint = 1
call pseudospin(MR1,N1,Z1,3,1,iprint)
call pseudospin(MR2,N2,Z2,3,1,iprint)
! rewrite the exchange matrix in the basis of local pseudospins:
do is1=1,N2
  do is2=1,N2
    call ZGEMM_('C','N',N1,N1,N1,cOne,Z1(1:N1,1:N1),N1,HEXCH(1:N1,1:N1,is1,is2),N1,cZero,TMP(1:N1,1:N1),N1)
    call ZGEMM_('N','N',N1,N1,N1,cOne,TMP(1:N1,1:N1),N1,Z1(1:N1,1:N1),N1,cZero,HEXCH2(1:N1,1:N1,is1,is2),N1)
  end do
end do
do is1=1,N1
  do is2=1,N1
    call ZGEMM_('C','N',N2,N2,N2,cOne,Z2(1:N2,1:N2),N2,HEXCH2(is1,is2,1:N2,1:N2),N2,cZero,TMP(1:N2,1:N2),N2)
    call ZGEMM_('N','N',N2,N2,N2,cOne,TMP(1:N2,1:N2),N2,Z2(1:N2,1:N2),N2,cZero,HEXCH3(is1,is2,1:N2,1:N2),N2)
  end do
end do
call JKQPar(N1,N2,HEXCH3,Jpar)
call tensor2cart(Jpar(1,-1:1,1,-1:1),Jc)

return

end subroutine project_exch
