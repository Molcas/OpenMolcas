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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DesymP(iAng,iCmp,jCmp,kCmp,lCmp,Shijij,iShll,iShell,iAO,kOp,ijkl,Aux,nAux,PAO,PSO,nPSO)
!***********************************************************************
!                                                                      *
!  Object: to transform the integrals in AO basis to symmetry adapted  *
!          orbitals , SO. This is done by accumulating the AO inte-    *
!          grals onto the SO integrals.                                *
!          Observe that one of the operator is the Unit operation      *
!          performed on center A. However, since we scramble the order *
!          we do not really know which center is which. However, the   *
!          Unit operator will always give a factor of one. Hence, this *
!          is a convenient way to do the symmetry transformation with- *
!          out having to know the order.                               *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to desymmetrization of the second order density *
!             matrix, January '92.                                     *
!***********************************************************************

use Index_Functions, only: nTri3_Elem
use Basis_Info, only: Shells
use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep, Prmt
use SOAO_Info, only: iAOtSO
use Real_Spherical, only: iSphCr
use Constants, only: Zero, Eight, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAng(4), iCmp, jCmp, kCmp, lCmp, iShll(4), iShell(4), iAO(4), kOp(4), ijkl, nAux, nPSO
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(out) :: Aux(nAux), PAO(ijkl,iCmp,jCmp,kCmp,lCmp)
real(kind=wp), intent(in) :: PSO(ijkl,nPSO)
integer(kind=iwp) :: i1, i2, i3, i4, iAux, iChBs, ii, is, iSym(0:7), j, j1, j12, j123, j2, j3, j4, jChBs, jj, js, jSym(0:7), &
                     kChBs, kk, ks, kSym(0:7), lChBs, ll, ls, lSym(0:7), MemSO2, niSym, njSym, nkSym, nlSym
real(kind=wp) :: Fact, FactNs, pa, pb, pc, Xa, Xb, Xg
logical(kind=iwp) :: Shij, Shkl
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
real(kind=wp), external :: DDot_
#endif

Shij = (iShell(1) == iShell(2))
Shkl = (iShell(3) == iShell(4))
MemSO2 = 1

#ifdef _DEBUGPRINT_
call RecPrt(' In DesymP: PSO ',' ',PSO,ijkl,nPSO)
call WrCheck(' In DesymP: PSO ',PSO,ijkl*nPSO)
write(u6,*) 'iCmp,jCmp,kCmp,lCmp,nPSO=',iCmp,jCmp,kCmp,lCmp,nPSO
write(u6,*) Shij,Shkl,Shijij
write(u6,*) 'kOp=',kOp
#endif
Fact = Eight
if (Shij) Fact = Fact*Half
if (Shkl) Fact = Fact*Half
if (Shijij) Fact = Fact*Half

! Initialize second order density matrix in AO basis

PAO(:,:,:,:,:) = Zero

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in PAO in a
! sequential way.

ii = nTri3_Elem(iAng(1))
jj = nTri3_Elem(iAng(2))
kk = nTri3_Elem(iAng(3))
ll = nTri3_Elem(iAng(4))
do i1=1,iCmp
  iChBs = iChBas(ii+i1)
  if (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
  pa = Prmt(iOper(kOp(1)),iChBs)
  niSym = 0
  do j=0,nIrrep-1
    if (iAOtSO(iAO(1)+i1,j) > 0) then
      iSym(niSym) = j
      niSym = niSym+1
    end if
  end do
  do i2=1,jCmp
    jChBs = iChBas(jj+i2)
    if (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
    pb = Prmt(iOper(kOp(2)),jChBs)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do
    do i3=1,kCmp
      kChBs = iChBas(kk+i3)
      if (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
      pc = Prmt(iOper(kOp(3)),kChBs)
      nkSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(3)+i3,j) > 0) then
          kSym(nkSym) = j
          nkSym = nkSym+1
        end if
      end do
      do i4=1,lCmp
        lChBs = iChBas(ll+i4)
        if (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
        ! Parity factor due to symmetry operations applied to the
        ! angular part of the basis functions.
        FactNs = pa*pb*pc*Prmt(iOper(kOp(4)),lChBs)
        nlSym = 0
        do j=0,nIrrep-1
          if (iAOtSO(iAO(4)+i4,j) > 0) then
            lSym(nlSym) = j
            nlSym = nlSym+1
          end if
        end do

        ! Loop over irreps which are spanned by the basis functions.
        ! Again, the loop structure is restricted to ensure unique
        ! integrals.

        if (nIrrep == 1) then
          !FactNs = 1
          PAO(:,i1,i2,i3,i4) = PAO(:,i1,i2,i3,i4)+Fact*PSO(:,MemSO2)
          MemSO2 = MemSO2+1
          cycle
        end if

        iAux = 0
        do is=0,niSym-1
          j1 = iSym(is)
          Xa = real(iChTbl(j1,kOp(1)),kind=wp)*FactNs
          do js=0,njSym-1
            j2 = jSym(js)
            Xb = real(iChTbl(j2,kOp(2)),kind=wp)*Xa
            j12 = ieor(j1,j2)
            do ks=0,nkSym-1
              j3 = kSym(ks)
              Xg = real(iChTbl(j3,kOp(3)),kind=wp)*Xb
              j123 = ieor(j12,j3)
              do ls=0,nlSym-1
                j4 = lSym(ls)
                if (j123 == j4) then
                  iAux = iAux+1
                  Aux(iAux) = real(iChTbl(j4,kOp(4)),kind=wp)*Xg*Fact
                  exit
                end if
              end do
            end do
          end do
        end do

        if (iAux /= 0) then
#         ifdef _DEBUGPRINT_
          call RecPrt(' Aux',' ',Aux,iAux,1)
#         endif
          if (iAux /= 1) then
            call DNaXpY(iAux,ijkl,Aux,1,PSO(1,MemSO2),1,ijkl,PAO(1,i1,i2,i3,i4),1,0)
          else
            PAO(:,i1,i2,i3,i4) = PAO(:,i1,i2,i3,i4)+Aux(1)*PSO(:,MemSO2)
          end if
          MemSO2 = MemSO2+iAux
        end if

      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' On exit from DesymP: PAO ',' ',PAO,ijkl,iCmp*jCmp*kCmp*lCmp)
do i=1,ijkl
  write(u6,*) DDot_(iCmp*jCmp*kCmp*lCmp,PAO(i,1,1,1,1),ijkl,PAO(i,1,1,1,1),ijkl)
end do
call WrCheck('DesymP: PAO ',PAO,ijkl*iCmp*jCmp*kCmp*lCmp)
#endif

end subroutine DesymP
