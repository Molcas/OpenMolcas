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

subroutine fetch_data_RunFile_all(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,angmom,eDmom,amfi,HSO,eso_au,esfs_au)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, Onei, auTocm, gElectron
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nss, nstate
integer(kind=iwp), intent(out) :: multiplicity(nstate)
real(kind=wp), intent(out) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), eDmom(3,nstate,nstate), amfi(3,nstate,nstate), &
                              eso_au(nss), esfs_au(nstate)
complex(kind=wp), intent(out) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss), DM(3,nss,nss)
complex(kind=wp), intent(inout) :: HSO(nss,nss)
integer(kind=iwp) :: i, i1, info, ipar, iss, ist, j, j1, jst, l, mult, multI, multJ, mxjob, ndata, njob
real(kind=wp) :: diff
logical(kind=iwp) :: found_amfi, found_edmom, found_hsoi, found_hsor
integer(kind=iwp), allocatable :: ibas(:,:), jbnum(:), mltplt(:), nstat(:)
real(kind=wp), allocatable :: tmpI(:,:), tmpR(:,:), W(:)
complex(kind=wp), allocatable :: M_tmp(:,:), tmp(:,:), u1(:,:)
complex(kind=wp), external :: Spin
real(kind=wp), parameter :: g_e = -gElectron, thr_deg = 0.2e-13_wp ! a.u. = 0.2e-13*auTocm = 4.38949263e-9 cm-1

! get basic sizes:
njob = 0
mxjob = 0
call get_iScalar('NJOB_SINGLE',njob)
call get_iScalar('MXJOB_SINGLE',mxjob)
! allocate temporary memory:
call mma_allocate(jbnum,nstate,'jbnum')
call mma_allocate(mltplt,mxjob,'mltplt')
call mma_allocate(nstat,mxjob,'nstat')
! get the information from RUNFILE:
mltplt = 0
jbnum = 0
nstat = 0
call get_iArray('MLTP_SINGLE',MLTPLT,MXJOB)
call get_iArray('JBNUM_SINGLE',JBNUM,NSTATE)
call get_iArray('NSTAT_SINGLE',NSTAT,MXJOB)
! computing the multiplicity of each state:
do i=1,nstate
  multiplicity(i) = mltplt(jbnum(i))
end do
call mma_deallocate(jbnum)
call mma_deallocate(mltplt)
call mma_deallocate(nstat)

! fetch the spin-orbit energies:
call get_dArray('ESO_SINGLE',eso,nss)
call get_dArray('ESO_LOW',eso_au,nss)

! fetch the spin-free energies:
call get_dArray('ESFS_SINGLE',esfs,nstate)
call get_dArray('ESFS_SINGLEAU',esfs_au,nstate)

! fetch the U matrix:
call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')

call get_dArray('UMATR_SINGLE',tmpR,nss*nss)
call get_dArray('UMATI_SINGLE',tmpI,nss*nss)
U(:,:) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)

! fetch the angular momentum integrals:
call get_dArray('ANGM_SINGLE',angmom,3*nstate*nstate)

! fetch the electric dipole moment integrals:
call qpg_dArray('DIP1_SINGLE',FOUND_EDMOM,NDATA)
if (found_edmom) then
  call get_dArray('DIP1_SINGLE',edmom,3*nstate*nstate)
else
  edmom(:,:,:) = Zero
end if

! fetch the amfi integrals:
call qpg_dArray('AMFI_SINGLE',FOUND_AMFI,NDATA)
if (found_amfi) then
  call get_dArray('AMFI_SINGLE',amfi,3*nstate*nstate)
else
  amfi(:,:,:) = Zero
end if

! fetch the spin-orbit hamiltonian
call qpg_dArray('HAMSOR_SINGLE',FOUND_HSOR,NDATA)
call qpg_dArray('HAMSOI_SINGLE',FOUND_HSOI,NDATA)
if (FOUND_HSOR .and. FOUND_HSOI) then
  call get_dArray('HAMSOR_SINGLE',tmpR,nss*nss)
  call get_dArray('HAMSOI_SINGLE',tmpI,nss*nss)
  HSO(:,:) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
  !---------------------------------------------------------------------
  ! if HSO is found, proceed to diagonalize it
  call mma_allocate(W,nss,'W')
  call mma_allocate(U1,nss,nss,'U1')
  info = 0
  call diag_c2(hso,nss,info,W,U1)
  ! correct for numerical degeneracies:
  do i=1,nss-1
    !wtmp = W(i)
    do j=i+1,nss
      diff = abs(w(i)-w(j))
      if (diff < thr_deg) w(j) = w(i)
    end do
  end do

  ESO(:) = (W(:)-W(1))*auTocm

  call mma_deallocate(W)
  call mma_deallocate(U1)
  !---------------------------------------------------------------------
end if
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

!-----
! generate a local indexing table:
call mma_allocate(ibas,[1,nstate],[-50,50],label='ibas')
iss = 0
ibas(:,:) = 0
ipar = mod(multiplicity(1),2)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) cycle
    Iss = Iss+1
    Ibas(Ist,I) = Iss
  end do ! i
end do ! Ist

!-----
! expand the spin free basis to the spin-orbit basis:
MM(:,:,:) = cZero
ML(:,:,:) = cZero
MS(:,:,:) = cZero
DM(:,:,:) = cZero
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) cycle
    do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
      if ((Ipar == 0) .and. (J == 0)) cycle
      i1 = Ibas(Ist,I)
      j1 = Ibas(Ist,J)
      do l=1,3
        MM(l,i1,j1) = -Spin(l,Mult,I,J)*g_e
        MS(l,i1,j1) = Spin(l,Mult,I,J)
      end do ! l
    end do ! J
  end do ! I
end do ! Ist

do Ist=1,nstate
  MultI = Multiplicity(Ist)
  do Jst=1,nstate
    MultJ = Multiplicity(Jst)
    if (MultI == MultJ) then
      do I=-(MultI-Ipar)/2,(MultI-Ipar)/2
        if ((Ipar == 0) .and. (I == 0)) cycle
        i1 = Ibas(Ist,I)
        j1 = Ibas(Jst,I)
        MM(:,i1,j1) = MM(:,i1,j1)-Angmom(:,Ist,Jst)*Onei
        ML(:,i1,j1) = ML(:,i1,j1)+Angmom(:,Ist,Jst)*Onei
        DM(:,i1,j1) = DM(:,i1,j1)+eDmom(:,Ist,Jst)*cOne
      end do ! I
    end if
  end do ! Jst
end do ! Ist

call mma_deallocate(ibas)

! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(M_tmp,nss,nss,'M_tmp')
call mma_allocate(tmp,nss,nss,'tmp')
do L=1,3
  ! spin moment
  M_tmp(:,:) = MS(L,:,:)
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,M_tmp,nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,M_tmp,nss)
  MS(L,:,:) = M_tmp(:,:)
  ! orbital moment
  M_tmp(:,:) = ML(L,:,:)
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,M_tmp,nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,M_tmp,nss)
  ML(L,:,:) = M_tmp(:,:)
  ! magnetic moment
  M_tmp(:,:) = MM(L,:,:)
  call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,M_tmp,nss,cZero,TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,M_tmp,nss)
  MM(L,:,:) = M_tmp(:,:)

  if (found_EDMOM) then
    ! electric dipole moment
    M_tmp(:,:) = DM(L,:,:)
    call ZGEMM_('C','N',nss,nss,nss,cOne,U,nss,M_tmp,nss,cZero,TMP,nss)
    call ZGEMM_('N','N',nss,nss,nss,cOne,TMP,nss,U,nss,cZero,M_tmp,nss)
    DM(L,:,:) = M_tmp(:,:)
  end if
end do ! L
call mma_deallocate(M_tmp)
call mma_deallocate(tmp)

! check the commutation rules of spin:
call check_commutation(nss,MS,.false.)

return

end subroutine fetch_data_RunFile_all
