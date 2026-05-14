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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine FckAcc_NoSym(iCmp,jCmp,kCmp,lCmp,Shijij,iShell,nijkl,AOInt,FMat,DMat,nDens,iAO,iAOst,iBas,jBas,kBas,lBas,ExFac)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmetry adapted Fock matrix.                        *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          In addition to this complication we have that the order of  *
!          indices in the integrals are not ordered canonically but    *
!          rather in an order such that the contraction step will be   *
!          optimal. Hence, special care has to be taken when tracing   *
!          the density with the integrals so that both entities have   *
!          the same order.                                             *
!                                                                      *
!          The Fock matrix is computed in lower triangular form.       *
!                                                                      *
!          The density matrix is not folded if the shell indices and   *
!          the angular indices are identical.                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. February '93                            *
!***********************************************************************

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use Gateway_Info, only: CutInt
use Constants, only: Zero, One, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCmp, jCmp, kCmp, lCmp, iShell(4), nijkl, nDens, iAO(4), iAOst(4), iBas, jBas, kBas, lBas
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(in) :: AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), DMat(nDens), ExFac
real(kind=wp), intent(inout) :: FMat(nDens)
integer(kind=iwp) :: i, i1, i2, i3, i4, ij, ik, il, iSO, j, jk, jl, jSO, k, kl, kSO, l, lSO, nijkl_
real(kind=wp) :: AOijkl, D_jk, D_jl, D_kl, F_jk, F_jl, F_kl, Fac, Fac_C, Fac_E
logical(kind=iwp) :: Shij, Shkl

!write(u6,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)
!if (iPrint >= 49) &
!  write(u6,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One,0)

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

Shij = (iShell(1) == iShell(2))
Shkl = (iShell(3) == iShell(4))

Fac = One
if (Shij) Fac = Fac*Half
if (Shkl) Fac = Fac*Half
if (Shijij) Fac = Fac*Half
Fac_C = Four*Fac
Fac_E = -Fac*ExFac

do i1=1,iCmp
  iSO = iAOtSO(iAO(1)+i1,0)+iAOst(1)
  do i2=1,jCmp
    jSO = iAOtSO(iAO(2)+i2,0)+iAOst(2)
    do i3=1,kCmp
      kSO = iAOtSO(iAO(3)+i3,0)+iAOst(3)
      do i4=1,lCmp
        lSO = iAOtSO(iAO(4)+i4,0)+iAOst(4)

        nijkl_ = 0
        do l=lSO,lSO+lBas-1
          do k=kSO,kSO+kBas-1
            kl = iTri(k,l)
            D_kl = DMat(kl)*Fac_C
            F_kl = Zero

            do j=jSO,jSO+jBas-1
              jk = iTri(j,k)
              jl = iTri(j,l)
              D_jl = DMat(jl)*Fac_E
              D_jk = DMat(jk)*Fac_E
              F_jl = Zero
              F_jk = Zero
              do i=iSO,iSO+iBas-1
                nijkl_ = nijkl_+1
                AOijkl = AOInt(nijkl_,i1,i2,i3,i4)
                if (abs(AOijkl) < CutInt) cycle
                ij = iTri(i,j)
                ik = iTri(i,k)
                il = iTri(i,l)

                FMat(ij) = FMat(ij)+AOijkl*D_kl
                F_kl = F_kl+AOijkl*DMat(ij)

                FMat(ik) = FMat(ik)+AOijkl*D_jl
                F_jl = F_jl+AOijkl*DMat(ik)
                FMat(il) = FMat(il)+AOijkl*D_jk
                F_jk = F_jk+AOijkl*DMat(il)

              end do
              FMat(jl) = FMat(jl)+F_jl*Fac_E
              FMat(jk) = FMat(jk)+F_jk*Fac_E
            end do
            FMat(kl) = FMat(kl)+F_kl*Fac_C
          end do
        end do

      end do
    end do
  end do
end do

end subroutine FckAcc_NoSym
