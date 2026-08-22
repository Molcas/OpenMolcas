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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine DWder_MCLR(mode,idsym,der1,nder1,der2,nder2,DWOut)

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use input_mclr, only: ERASSCF, iToc, nAsh, nBas, nIsh, nRoots, nSym, ntAsh, weight
use DWSol, only: DWSol_der
use MCLR_Data, only: ipCM, LUJOB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: mode, idSym, nder1, nder2
real(kind=wp), intent(in) :: der1(nder1), der2(nder2)
real(kind=wp), intent(inout) :: DWOut(nRoots,nRoots)
integer(kind=iwp) :: i, iA, ij1, ip1, iS, j, jA, jS, jDisk, kA, kl1, lA, ng1, ng2
real(kind=wp) :: OMGDER, rdum(1), Scal
real(kind=wp), allocatable :: G1q(:), G2q(:), DERHII(:), DEROMG(:)

if (nRoots == 1) return

! derivative of the dynamical weight
! I think this implementation uses the following constraint condition
! w_I^Z * (<I|H|I> - E_I)
! where E is used as a parameter, and w_I^Z is the undetermined multiplier (not weight itself)
!
! For DW-MCSCF, partial derivative of the Lagrangian wrt H_{IJ} can be evaluated with either orbital or CI rotations.
! At present, it is evaluated in ISR_TimesE2.
! For DW solvation, it is only with orbital rotations (here)?

ng1 = nTri_Elem(ntash)
ng2 = nTri_Elem(ng1)
call mma_allocate(G1q,ng1,Label='G1q')
call mma_allocate(G2q,ng2,Label='G2q')
call mma_allocate(DEROMG,nRoots,Label='DEROMG')
DEROMG = Zero

jdisk = itoc(3)
do j=1,nRoots
  call dDaFile(LUJOB,2,G1q,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,2,G2q,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)

  OMGDER = Zero
  do iS=1,nSym
    if (nAsh(iS) > 0) then
      jS = Mul(iS,idSym)
      do iA=1,nAsh(is)
        do jA=1,nAsh(js)
          ij1 = iTri(ia,ja)
          ip1 = nIsh(iS)+iA-1+nBas(jS)*(nIsh(js)+jA-1)+ipCM(is)
          OMGDER = OMGDER+der1(ip1)*G1q(ij1)
          if (mode == 1) then
            !! DW-MCSCF
            do kA=1,nAsh(is)
              do lA=1,nAsh(js)
                kl1 = iTri(ka,la)
                !! dimension of rmoaa is itri(nnA**2,nnA**2), but the actual array is similar to G2q
                !! dimension of G2q  : itri(itri(nnA,nnA),itri(nnA,nnA))
                scal = One
                if ((ij1 >= kl1) .and. (kA == lA)) then
                  scal = Two
                else if ((ij1 < kl1) .and. (iA == jA)) then
                  scal = Two
                end if
                OMGDER = OMGDER+scal*der2(itri(ij1,kl1))*G2q(itri(ij1,kl1))*Half
              end do
            end do
          end if
        end do
      end do
    end if
  end do
  DEROMG(j) = OMGDER
end do

call mma_allocate(DERHII,nRoots,Label='DERHII')
DERHII(:) = Zero

!! weight is optional
call DWSol_Der(mode,DEROMG,DERHII,ERASSCF,weight)

!! this contributions should not be treated as multipliers, I guess
do i=1,nRoots
  DWOut(i,i) = DWOut(i,i)+DERHII(i)*Half
end do
!write(u6,*) 'dwout'
!call sqprt(dwout,nroots)

call mma_deallocate(G1q)
call mma_deallocate(G2q)
call mma_deallocate(DERHII)
call mma_deallocate(DEROMG)

return

end subroutine DWder_MCLR
