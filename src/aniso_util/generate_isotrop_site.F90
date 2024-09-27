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

subroutine generate_isotrop_site(nss,nsfs,nexch,gtens_input,riso,D,EoverD,E,M,S)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, cZero, cOne
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: nss, nsfs, nexch
real(kind=wp), intent(in) :: gtens_input(3), riso(3,3), D, EoverD ! ZFS factors
real(kind=wp), intent(out) :: E(nExch) ! spin-orbit energy states, starting from 0
complex(kind=wp), intent(out) :: M(3,nExch,nExch), S(3,nExch,nExch)
integer(kind=iwp) :: info, l
real(kind=wp) :: RM, RS
complex(kind=wp) :: redme
real(kind=wp), allocatable :: W(:)
complex(kind=wp), allocatable :: HZFS(:,:), MTMP(:,:,:), S2(:,:), STMP(:,:,:), SX2(:,:), SY2(:,:), SZ2(:,:), tmp(:,:), tmp2(:,:), &
                                 Z(:,:)
real(kind=wp), external :: dznrm2_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, j
real(kind=wp) :: g(3), gtens(3), ma(3,3), maxes(3,3)
real(kind=wp), external :: dnrm2_
#endif
!----------------------------------------------------------------------|

nsfs = 1
nss = nexch

#ifdef _DEBUGPRINT_
write(u6,*) 'GENERATE_SITE:  nss   = ',nss
write(u6,*) 'GENERATE_SITE:  nsfs  = ',nsfs
write(u6,*) 'GENERATE_SITE:  gfact = ',(gtens_input(l),l=1,3)
write(u6,*) 'GENERATE_SITE:  EoverD= ',EoverD
write(u6,*) 'GENERATE_SITE:  D     = ',D
write(u6,*) 'GENERATE_SITE:  riso  = ',((riso(i,j),i=1,3),j=1,3)
#endif

E(:) = Zero
S(:,:,:) = cZero
M(:,:,:) = cZero

if (nss >= 2) then
  call mma_allocate(tmp,nss,nss,'tmp')
  call mma_allocate(tmp2,nss,nss,'tmp2')
  call ESO(nss,1,1,tmp,tmp2,redME)
  S(1,:,:) = tmp(:,:)
  S(2,:,:) = tmp2(:,:)
  call ESO(nss,1,0,tmp,tmp2,redME)
  S(3,:,:) = tmp(:,:)
  call mma_deallocate(tmp2)
  M(1,:,:) = -gtens_input(1)*S(1,:,:)
  M(2,:,:) = -gtens_input(2)*S(2,:,:)
  M(3,:,:) = -gtens_input(3)*S(3,:,:)

# ifdef _DEBUGPRINT_
  call prmom('GENERATE_SITE:     SPIN MOMENT:',S,nss)
  call prmom('GENERATE_SITE: MAGNETIC MOMENT:',M,nss)
# endif

  RM = dznrm2_(3*nss*nss,M,1)
  RS = dznrm2_(3*nss*nss,S,1)
# ifdef _DEBUGPRINT_
  write(u6,'(A,2ES22.14)') 'Norms of M and S:',RM,RS
# endif

  if ((RM > Zero) .and. (RS > Zero)) then
    ! rotate the spin and magnetic moment to the general coordinate system:
    call mma_allocate(MTMP,3,nExch,nExch,'MTMP')
    call mma_allocate(STMP,3,nExch,nExch,'STMP')
    call mma_allocate(Z,nExch,nExch,'Z')

    MTMP(:,:,:) = M(:,:,:)
    STMP(:,:,:) = S(:,:,:)
    Z(:,:) = Zero

#   ifdef _DEBUGPRINT_
    write(u6,'(A,ES20.10)') 'GENERATE_SITE: Norm of riso:',dnrm2_(9,riso,1)
#   endif

    call rotmom(STMP,nExch,riso,S)
    call rotmom(MTMP,nExch,riso,M)

#   ifdef _DEBUGPRINT_
    gtens(:) = Zero
    maxes(:,:) = Zero
    call atens(M,nExch,gtens,maxes,2)
#   endif

    if (abs(D) > Zero) then
      ! compute the ZFS
      call mma_allocate(HZFS,nExch,nExch,'HZFS')
      call mma_allocate(SX2,nExch,nExch,'SX2')
      call mma_allocate(SY2,nExch,nExch,'SY2')
      call mma_allocate(SZ2,nExch,nExch,'SZ2')
      call mma_allocate(S2,nExch,nExch,'S2')
      call mma_allocate(W,nExch,'W')

      W(:) = Zero

      call mma_allocate(tmp2,nEXCH,nEXCH,label='tmp2')
      tmp2(:,:) = S(1,:,:)
      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,cOne,tmp2,nEXCH,tmp2,nEXCH,cZero,SX2,nEXCH)
      tmp2(:,:) = S(2,:,:)
      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,cOne,tmp2,nEXCH,tmp2,nEXCH,cZero,SY2,nEXCH)
      tmp2(:,:) = S(3,:,:)
      call zgemm_('C','N',nEXCH,nEXCH,nEXCH,cOne,tmp2,nEXCH,tmp2,nEXCH,cZero,SZ2,nEXCH)

      S2(:,:) = SX2(:,:)+SY2(:,:)+SZ2(:,:)

#     ifdef _DEBUGPRINT_
      call pa_prMat('GENERATE_SITE: SX2',SX2,nexch)
      call pa_prMat('GENERATE_SITE: SY2',SY2,nexch)
      call pa_prMat('GENERATE_SITE: SZ2',SZ2,nexch)
      call pa_prMat('GENERATE_SITE: S2 ',S2,nexch)
#     endif

      HZFS(:,:) = D*(SZ2(:,:)-S2(:,:)/Three)+D*EoverD*(SX2(:,:)-SY2(:,:))

#     ifdef _DEBUGPRINT_
      call print_ZFS('GENERATE_SITE: ZFS matrix:',HZFS,nExch)
#     endif

      info = 0
      call diag_c2(HZFS,nExch,info,W,Z)

      E(:) = W(:)-W(1)
#     ifdef _DEBUGPRINT_
      do i=1,nExch
        write(u6,'(A,i2,A,F14.8)') 'ZFS  E(',i,') = ',E(i)
      end do
#     endif
      ! rotate the spin and magnetic moment:

      do L=1,3
        tmp2(:,:) = M(L,:,:)
        call zgemm_('C','N',nEXCH,nEXCH,nEXCH,cOne,Z,nEXCH,tmp2,nEXCH,cZero,TMP,nEXCH)
        call zgemm_('N','N',nEXCH,nEXCH,nEXCH,cOne,TMP,nEXCH,Z,nEXCH,cZero,tmp2,nEXCH)
        M(L,:,:) = tmp2(:,:)

        tmp2(:,:) = S(L,:,:)
        call zgemm_('C','N',nEXCH,nEXCH,nEXCH,cOne,Z,nEXCH,tmp2,nEXCH,cZero,TMP,nEXCH)
        call zgemm_('N','N',nEXCH,nEXCH,nEXCH,cOne,TMP,nEXCH,Z,nEXCH,cZero,tmp2,nEXCH)
        S(L,:,:) = tmp2(:,:)
      end do  ! L

      call mma_deallocate(SX2)
      call mma_deallocate(SY2)
      call mma_deallocate(SZ2)
      call mma_deallocate(S2)
      call mma_deallocate(W)
      call mma_deallocate(HZFS)
      call mma_deallocate(tmp2)
    end if ! ZFS is defined

    call mma_deallocate(MTMP)
    call mma_deallocate(STMP)
    call mma_deallocate(Z)
  end if ! dznrm2_ M and S
  call mma_deallocate(tmp)

# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'g tensor at the end of GENERATE_SPIN'
  g(:) = Zero
  ma(:,:) = Zero
  call atens(M,nExch,g,ma,2)
# endif
end if ! nss>=2

return

end subroutine generate_isotrop_site
