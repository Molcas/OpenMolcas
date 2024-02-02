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

subroutine generate_isotrop_site(nss,nsfs,nexch,nLoc,gtens_input,riso,D,EoverD,E,M,S)

implicit none
#include "warnings.h"
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: nLoc
integer, intent(inout) :: nexch
integer, intent(inout) :: nss, nsfs
real(kind=8), intent(in) :: gtens_input(3)
real(kind=8), intent(in) :: riso(3,3)
real(kind=8), intent(in) :: D, EoverD ! ZFS factors
! spin-orbit energy states, starting from 0
real(kind=8), intent(out) :: E(nExch)
complex(kind=8), intent(out) :: M(3,nExch,nExch)
complex(kind=8), intent(out) :: S(3,nExch,nExch)
! local variables:
integer :: i, j, l
integer :: info
!complex(kind=8) :: spin
complex(kind=8) :: redme
complex(kind=8), allocatable :: HZFS(:,:), S2(:,:), Wc(:,:)
complex(kind=8), allocatable :: SX2(:,:), SY2(:,:), SZ2(:,:)
complex(kind=8), allocatable :: Z(:,:), tmp(:,:)
complex(kind=8), allocatable :: MTMP(:,:,:), STMP(:,:,:)
real(kind=8), allocatable :: W(:), gtens(:), maxes(:,:)
real(kind=8) :: dznrm2_, RM, RS, dnrm2_
real(kind=8) :: g(3), ma(3,3)
external :: spin, dznrm2_, dnrm2_
logical :: dbg
dbg = .false.
!----------------------------------------------------------------------|

nsfs = 1
nss = nexch

if (dbg) then
  write(6,*) 'GENERATE_SITE:  nss   = ',nss
  write(6,*) 'GENERATE_SITE:  nsfs  = ',nsfs
  write(6,*) 'GENERATE_SITE:  nLoc  = ',nLoc
  write(6,*) 'GENERATE_SITE:  gfact = ',(gtens_input(l),l=1,3)
  write(6,*) 'GENERATE_SITE:  EoverD= ',EoverD
  write(6,*) 'GENERATE_SITE:  D     = ',D
  write(6,*) 'GENERATE_SITE:  riso  = ',((riso(i,j),i=1,3),j=1,3)
end if

do i=1,nExch
  E(i) = 0.0_wp
end do

call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,S,1)
call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,M,1)

if (nss >= 2) then
  call mma_allocate(Wc,nss,nss,'Wc')
  call ESO(nss,1,1,S(1,1:nss,1:nss),S(2,1:nss,1:nss),redME)
  call ESO(nss,1,0,S(3,1:nss,1:nss),Wc(1:nss,1:nss),redME)
  call mma_deallocate(Wc)
  call zcopy_(3*nss*nss,S,1,M,1)
  call zdscal_(nss*nss,-gtens_input(1),M(1,1:nss,1:nss),1)
  call zdscal_(nss*nss,-gtens_input(2),M(2,1:nss,1:nss),1)
  call zdscal_(nss*nss,-gtens_input(3),M(3,1:nss,1:nss),1)

  if (dbg) then
    call prmom('GENERATE_SITE:     SPIN MOMENT:',S,nss)
    call prmom('GENERATE_SITE: MAGNETIC MOMENT:',M,nss)
  end if

  RM = dznrm2_(3*nss*nss,M(1:3,1:nss,1:nss),1)
  RS = dznrm2_(3*nss*nss,S(1:3,1:nss,1:nss),1)
  if (dbg) write(6,'(A,2ES22.14)') 'Norms of M and S:',RM,RS

  if ((RM > 0.0_wp) .and. (RS > 0.0_wp)) then
    ! rotate the spin and magnetic moment to the general coordinate system:
    call mma_allocate(MTMP,3,nExch,nExch,'MTMP')
    call mma_allocate(STMP,3,nExch,nExch,'STMP')
    call mma_allocate(Z,nExch,nExch,'Z')
    call mma_allocate(tmp,nExch,nExch,'tmp')

    call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,MTMP,1)
    call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,STMP,1)

    call zcopy_(3*nExch*nExch,M,1,MTMP,1)
    call zcopy_(3*nExch*nExch,S,1,STMP,1)
    call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,M,1)
    call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,S,1)
    call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,Z,1)

    if (dbg) write(6,'(A,ES20.10)') 'GENERATE_SITE: Norm of riso:',dnrm2_(9,riso,1)

    call zcopy_(3*nExch*nExch,MTMP,1,M,1)
    call zcopy_(3*nExch*nExch,STMP,1,S,1)
    call rotmom(STMP,nExch,riso,S)
    call rotmom(MTMP,nExch,riso,M)

    if (dbg) then
      call mma_allocate(gtens,3,'gtens')
      call mma_allocate(maxes,3,3,'maxes')
      call dcopy_(3,[0.0_wp],0,gtens,1)
      call dcopy_(3*3,[0.0_wp],0,maxes,1)
      call atens(M,nExch,gtens,maxes,2)
      call mma_deallocate(gtens)
      call mma_deallocate(maxes)
    end if

    if (abs(D) > 0.0_wp) then
      ! compute the ZFS
      call mma_allocate(HZFS,nExch,nExch,'HZFS')
      call mma_allocate(SX2,nExch,nExch,'SX2')
      call mma_allocate(SY2,nExch,nExch,'SY2')
      call mma_allocate(SZ2,nExch,nExch,'SZ2')
      call mma_allocate(S2,nExch,nExch,'S2')
      call mma_allocate(W,nExch,'W')

      call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,HZFS,1)
      call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SX2,1)
      call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SY2,1)
      call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SZ2,1)
      call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,S2,1)
      call dcopy_(nExch,[0.0_wp],0,W,1)

      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),S(1,:,:),nEXCH,S(1,:,:),nEXCH,(0.0_wp,0.0_wp),SX2(:,:),nEXCH)
      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),S(2,:,:),nEXCH,S(2,:,:),nEXCH,(0.0_wp,0.0_wp),SY2(:,:),nEXCH)
      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),S(3,:,:),nEXCH,S(3,:,:),nEXCH,(0.0_wp,0.0_wp),SZ2(:,:),nEXCH)

      S2(:,:) = SX2(:,:)+SY2(:,:)+SZ2(:,:)

      if (dbg) then
        call pa_prMat('GENERATE_SITE: SX2',SX2,nexch)
        call pa_prMat('GENERATE_SITE: SY2',SY2,nexch)
        call pa_prMat('GENERATE_SITE: SZ2',SZ2,nexch)
        call pa_prMat('GENERATE_SITE: S2 ',S2,nexch)
      end if

      HZFS(:,:) = D*(SZ2(:,:)-S2(:,:)/3.0_wp)+D*EoverD*(SX2(:,:)-SY2(:,:))

      if (dbg) call print_ZFS('GENERATE_SITE: ZFS matrix:',HZFS,nExch)

      info = 0
      call diag_c2(HZFS,nExch,info,W,Z)

      do i=1,nExch
        E(i) = W(i)-W(1)
      end do
      if (dbg) then
        do i=1,nExch
          write(6,'(A,i2,A,F14.8)') 'ZFS  E(',i,') = ',E(i)
        end do
      end if
      ! rotate the spin and magnetic moment:

      do L=1,3
        call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,TMP,1)
        call zgemm_('C','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),Z,nEXCH,M(L,:,:),nEXCH,(0.0_wp,0.0_wp),TMP,nEXCH)
        call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,M(L,:,:),1)
        call zgemm_('N','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),TMP,nEXCH,Z,nEXCH,(0.0_wp,0.0_wp),M(L,:,:),nEXCH)

        call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,TMP,1)
        call zgemm_('C','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),Z,nEXCH,S(L,:,:),nEXCH,(0.0_wp,0.0_wp),TMP,nEXCH)
        call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,S(L,:,:),1)
        call zgemm_('N','N',nEXCH,nEXCH,nEXCH,(1.0_wp,0.0_wp),TMP,nEXCH,Z,nEXCH,(0.0_wp,0.0_wp),S(L,:,:),nEXCH)
      end do  ! L

      call mma_deallocate(SX2)
      call mma_deallocate(SY2)
      call mma_deallocate(SZ2)
      call mma_deallocate(S2)
      call mma_deallocate(W)
      call mma_deallocate(HZFS)
    end if ! ZFS is defined

    call mma_deallocate(MTMP)
    call mma_deallocate(STMP)
    call mma_deallocate(Z)
    call mma_deallocate(tmp)
  end if ! dznrm2_ M and S

  if (dbg) then
    write(6,'(A)') 'g tensor at the end of GENERATE_SPIN'
    g = 0.0_wp
    ma = 0.0_wp
    call atens(M,nExch,g,ma,2)
  end if
end if ! nss>=2

return

end subroutine generate_isotrop_site
