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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer :: nss, nstate
integer :: multiplicity(nstate)
real(kind=8) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), eDmom(3,nstate,nstate), amfi(3,nstate,nstate), eso_au(nss), &
                esfs_au(nstate)
complex(kind=8) :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
complex(kind=8) :: DM(3,nss,nss)
complex(kind=8) :: U(nss,nss), HSO(nss,nss)
! local variables:
integer :: njob, mxjob, ndata, iss, ibas(nstate,-50:50)
integer :: i, j, i1, j1, ist, jst, mult, multI, multJ
integer :: l, ipar, info
real(kind=8) :: g_e, au2cm, thr_deg, diff
! allocatable local arrays:
integer, allocatable :: mltplt(:), jbnum(:), nstat(:) !,lroot(:)
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:), W(:)
complex(kind=8), allocatable :: tmp(:,:), u1(:,:)
complex(kind=8) :: Spin
external :: Spin
logical :: found_edmom, found_amfi, found_hsor, found_hsoi

g_e = 2.00231930437180_wp
au2cm = 219474.6313702_wp
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
multiplicity = 0
do i=1,nstate
  multiplicity(i) = mltplt(jbnum(i))
end do
call mma_deallocate(jbnum)
call mma_deallocate(mltplt)
call mma_deallocate(nstat)

! fetch the spin-orbit energies:
eso = 0.0_wp
eso_au = 0.0_wp
call get_dArray('ESO_SINGLE',eso,nss)
call get_dArray('ESO_LOW',eso_au,nss)

! fetch the spin-free energies:
esfs = 0.0_wp
esfs_au = 0.0_wp
call get_dArray('ESFS_SINGLE',esfs,nstate)
call get_dArray('ESFS_SINGLEAU',esfs_au,nstate)

! fetch the U matrix:
call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')

tmpR = 0.0_wp
tmpI = 0.0_wp
call get_dArray('UMATR_SINGLE',tmpR,nss*nss)
call get_dArray('UMATI_SINGLE',tmpI,nss*nss)
U = (0.0_wp,0.0_wp)
do i=1,nss
  do j=1,nss
    U(i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
  end do
end do

! fetch the angular momentum integrals:
angmom = 0.0_wp
call get_dArray('ANGM_SINGLE',angmom,3*nstate*nstate)

! fetch the electric dipole moment integrals:
edmom = 0.0_wp
found_EDMOM = .false.
call qpg_dArray('DIP1_SINGLE',FOUND_EDMOM,NDATA)
if (found_edmom) call get_dArray('DIP1_SINGLE',edmom,3*nstate*nstate)

! fetch the amfi integrals:
amfi = 0.0_wp
found_AMFI = .false.
call qpg_dArray('AMFI_SINGLE',FOUND_AMFI,NDATA)
if (found_amfi) call get_dArray('AMFI_SINGLE',amfi,3*nstate*nstate)

! fetch the spin-orbit hamiltonian
FOUND_HSOR = .false.
FOUND_HSOI = .false.
call qpg_dArray('HAMSOR_SINGLE',FOUND_HSOR,NDATA)
call qpg_dArray('HAMSOI_SINGLE',FOUND_HSOI,NDATA)
if (FOUND_HSOR .and. FOUND_HSOI) then
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  call get_dArray('HAMSOR_SINGLE',tmpR,nss*nss)
  call get_dArray('HAMSOI_SINGLE',tmpI,nss*nss)
  call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,HSO,1)
  do i=1,nss
    do j=1,nss
      HSO(i,j) = cmplx(tmpR(i,j),tmpI(i,j),wp)
    end do
  end do
  !---------------------------------------------------------------------
  ! if HSO is found, proceed to diagonalize it
  call mma_allocate(W,nss,'W')
  call mma_allocate(U1,nss,nss,'U1')
  call dcopy_(nss,[0.0_wp],0,W,1)
  call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,U1,1)
  info = 0
  call diag_c2(hso,nss,info,W,U1)
  ! correct for numerical degeneracies:
  thr_deg = 0.2D-13 ! a.u. = 0.2D-13*au2cm = 4.38949263E-09 cm-1
  do i=1,nss-1
    !wtmp = W(i)
    do j=i+1,nss
      diff = abs(w(i)-w(j))
      if (diff < thr_deg) w(j) = w(i)
    end do
  end do

  do i=1,nss
    ESO(i) = (W(i)-W(1))*au2cm
  end do

  call mma_deallocate(W)
  call mma_deallocate(U1)
  !---------------------------------------------------------------------
end if
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

!-----
! generate a local indexing table:
iss = 0
ibas = 0
ipar = mod(multiplicity(1),2)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) go to 310
    Iss = Iss+1
    Ibas(Ist,I) = Iss
310 continue
  end do ! i
end do ! ist

!-----
! expand the spin free basis to the spin-orbit basis:
call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,MM,1)
call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,ML,1)
call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,MS,1)
do Ist=1,nstate
  Mult = Multiplicity(Ist)
  do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
    if ((Ipar == 0) .and. (I == 0)) go to 301
    do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
      if ((Ipar == 0) .and. (J == 0)) go to 302
      do l=1,3
        i1 = Ibas(Ist,I)
        j1 = Ibas(Ist,J)
        MM(l,i1,j1) = -Spin(l,Mult,I,J)*g_e
        MS(l,i1,j1) = Spin(l,Mult,I,J)
      end do ! l
302   continue
    end do ! J
301 continue
  end do ! I
end do ! Ist

do Ist=1,nstate
  MultI = Multiplicity(Ist)
  do Jst=1,nstate
    MultJ = Multiplicity(Jst)
    if (MultI == MultJ) then
      do I=-(MultI-Ipar)/2,(MultI-Ipar)/2
        if ((Ipar == 0) .and. (I == 0)) go to 303
        do l=1,3
          i1 = Ibas(Ist,I)
          j1 = Ibas(Jst,I)
          MM(l,i1,j1) = MM(l,i1,j1)-cmplx(0.0_wp,Angmom(l,Ist,Jst),wp)
          ML(l,i1,j1) = ML(l,i1,j1)+cmplx(0.0_wp,Angmom(l,Ist,Jst),wp)
          DM(l,i1,j1) = DM(l,i1,j1)+cmplx(eDmom(l,Ist,Jst),0.0_wp,wp)
        end do   ! l
303     continue
      end do   ! I
    end if
  end do   ! Jst
end do   ! Ist

! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(tmp,nss,nss,'tmp')
do L=1,3
  TMP = (0.0_wp,0.0_wp)
  ! spin moment
  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,MS(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),MS(L,:,:),nss)
  ! orbital moment
  TMP = (0.0_wp,0.0_wp)
  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,ML(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),ML(L,:,:),nss)
  ! magnetic moment
  TMP = (0.0_wp,0.0_wp)
  call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,MM(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
  call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),MM(L,:,:),nss)

  if (found_EDMOM) then
    ! electric dipole moment
    TMP = (0.0_wp,0.0_wp)
    call ZGEMM_('C','N',nss,nss,nss,(1.0_wp,0.0_wp),U,nss,DM(L,:,:),nss,(0.0_wp,0.0_wp),TMP,nss)
    call ZGEMM_('N','N',nss,nss,nss,(1.0_wp,0.0_wp),TMP,nss,U,nss,(0.0_wp,0.0_wp),DM(L,:,:),nss)
  end if
end do ! L
call mma_deallocate(tmp)

! check the commutation rules of spin:
call check_commutation(nss,MS(1:3,1:nss,1:nss),.false.)

return

end subroutine fetch_data_RunFile_all
