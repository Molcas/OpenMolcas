************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine generate_isotrop_site( nss, nsfs, nexch, nLoc,
     &                                  gtens_input,riso,D,EoverD,
     &                                  E, M, S )

      Implicit none
#include "warnings.fh"
#include "stdalloc.fh"
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in)           :: nLoc
      Integer, intent(inout)        :: nexch
      Integer, intent(inout)        :: nss, nsfs
      Real(kind=wp), intent(in)     :: gtens_input(3)
      Real(kind=wp), intent(in)     :: riso(3,3)
      Real(kind=wp), intent(in)     :: D, EoverD ! ZFS factors
      ! spin-orbit energy states, starting from 0
      Real(kind=wp), intent(out)    :: E(nExch)
      Complex(kind=wp), intent(out) :: M(3,nExch,nExch)
      Complex(kind=wp), intent(out) :: S(3,nExch,nExch)

      ! local variables:
      Integer         :: ib(-30:30), i, j, l
      Integer         :: info
!      Complex(kind=wp):: spin
      Complex(kind=wp) :: redme
      Complex(kind=wp), allocatable :: HZFS(:,:), S2(:,:), Wc(:,:)
      Complex(kind=wp), allocatable :: SX2(:,:), SY2(:,:), SZ2(:,:)
      Complex(kind=wp), allocatable :: Z(:,:), tmp(:,:)
       Complex(kind=wp), allocatable :: MTMP(:,:,:), STMP(:,:,:)
      Real(kind=wp), allocatable    :: W(:),gtens(:),maxes(:,:)
      Real(kind=wp)                 :: dznrm2_, RM, RS, dnrm2_
      Real(kind=wp)                 :: g(3), ma(3,3)
      External                      :: spin, dznrm2_, dnrm2_
      Logical                       :: dbg
      dbg=.false.
!----------------------------------------------------------------------|
      Call qEnter('generate_isotrop_site')

      Do j=-30,30
         ib(j)=0
      End Do

      nsfs=1
      nss  = nexch

      If(dbg) Write(6,*) 'GENERATE_SITE:  nss   = ', nss
      If(dbg) Write(6,*) 'GENERATE_SITE:  nsfs  = ', nsfs
      If(dbg) Write(6,*) 'GENERATE_SITE:  nLoc  = ', nLoc
      If(dbg) Write(6,*) 'GENERATE_SITE:  gfact = ',
     &                                           (gtens_input(l),l=1,3)
      If(dbg) Write(6,*) 'GENERATE_SITE:  EoverD= ', EoverD
      If(dbg) Write(6,*) 'GENERATE_SITE:  D     = ', D
      If(dbg) Write(6,*) 'GENERATE_SITE:  riso  = ',
     &                                        ((riso(i,j),i=1,3),j=1,3)

      Do i=1,nExch
         E(i)=0.0_wp
      End Do

      Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,S,1)
      Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,M,1)

      Call mma_allocate(Wc,nss,nss,'Wc')
      Call ESO(nss,1,1,S(1,1:nss,1:nss),S(2,1:nss,1:nss),redME)
      Call ESO(nss,1,0,S(3,1:nss,1:nss),Wc( 1:nss,1:nss),redME)
      Call mma_deallocate(Wc)
      Call zcopy_(3*nss*nss,S,1,M,1)
      Call zdscal_(nss*nss,-gtens_input(1),M(1,1:nss,1:nss),1)
      Call zdscal_(nss*nss,-gtens_input(2),M(2,1:nss,1:nss),1)
      Call zdscal_(nss*nss,-gtens_input(3),M(3,1:nss,1:nss),1)

      If(dbg) Call prmom('GENERATE_SITE:     SPIN MOMENT:',S,nss)
      If(dbg) Call prmom('GENERATE_SITE: MAGNETIC MOMENT:',M,nss)

      RM=dznrm2_(3*nss*nss,M(1:3,1:nss,1:nss),1)
      RS=dznrm2_(3*nss*nss,S(1:3,1:nss,1:nss),1)
      If(dbg) Write(6,'(A,2ES22.14)') 'Norms of M and S:', RM, RS

      If( (RM>0.0_wp).AND.(RS>0.0_wp) ) Then
!        rotate the spin and magnetic moment to the general coordinate system:
         Call mma_allocate(MTMP,3,nExch,nExch,'MTMP')
         Call mma_allocate(STMP,3,nExch,nExch,'STMP')
         Call mma_allocate(Z,nExch,nExch,'Z')
         Call mma_allocate(tmp,nExch,nExch,'tmp')

         Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,MTMP,1)
         Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,STMP,1)

         Call zcopy_(3*nExch*nExch,M,1,MTMP,1)
         Call zcopy_(3*nExch*nExch,S,1,STMP,1)
         Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,M,1)
         Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,S,1)
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,Z,1)


         If(dbg) Write(6,'(A,ES20.10)') 'GENERATE_SITE: Norm of riso:',
     &                                   dnrm2_(9,riso,1)
         If(dbg) Call xFlush(6)

         Call zcopy_(3*nExch*nExch,MTMP,1,M,1)
         Call zcopy_(3*nExch*nExch,STMP,1,S,1)
         Call rotmom( STMP, nExch, riso, S )
         Call rotmom( MTMP, nExch, riso, M )

         If(dbg) Then
            Call mma_allocate(gtens,3,'gtens')
            Call mma_allocate(maxes,3,3,'maxes')
            Call dcopy_(3,[0.0_wp],0,gtens,1)
            Call dcopy_(3*3,[0.0_wp],0,maxes,1)
            Call atens(M, nExch, gtens, maxes, 2)
            Call mma_deallocate(gtens)
            Call mma_deallocate(maxes)
         End If

         If( abs(D)>0.0_wp ) Then
         ! compute the ZFS
            Call mma_allocate(HZFS,nExch,nExch,'HZFS')
            Call mma_allocate(SX2,nExch,nExch,'SX2')
            Call mma_allocate(SY2,nExch,nExch,'SY2')
            Call mma_allocate(SZ2,nExch,nExch,'SZ2')
            Call mma_allocate(S2,nExch,nExch,'S2')
            Call mma_allocate(W,nExch,'W')

            Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,HZFS,1)
            Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SX2,1)
            Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SY2,1)
            Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,SZ2,1)
            Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,S2,1)
            Call dcopy_(nExch     , [0.0_wp],        0,W,1)

            SX2(:,:)=MATMUL(S(1,:,:),S(1,:,:))
            SY2(:,:)=MATMUL(S(2,:,:),S(2,:,:))
            SZ2(:,:)=MATMUL(S(3,:,:),S(3,:,:))
            S2(:,:)=SX2(:,:)+SY2(:,:)+SZ2(:,:)

            If(dbg) Call pa_prMat('GENERATE_SITE: SX2',SX2,nexch)
            If(dbg) Call pa_prMat('GENERATE_SITE: SY2',SY2,nexch)
            If(dbg) Call pa_prMat('GENERATE_SITE: SZ2',SZ2,nexch)
            If(dbg) Call pa_prMat('GENERATE_SITE: S2 ',S2, nexch)

            HZFS(:,:)=D          * ( SZ2(:,:)-S2(:,:)/3.0_wp) +
     &                D * EoverD * ( SX2(:,:)-SY2(:,:)      )

            If(dbg) Call print_ZFS('GENERATE_SITE: ZFS matrix:',HZFS,
     &                              nExch)

            info=0
            Call diag_c2(HZFS,nExch,info,W,Z)

            Do i=1,nExch
               E(i)=W(i)-W(1)
            End Do
            If(dbg) Then
              Do i=1,nExch
                Write(6,'(A,i2,A,F14.8)') 'ZFS  E(',i,') = ',E(i)
              End Do
            End If
            ! rotate the spin and magnetic moment:

            Do L=1,3
              Call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,TMP,1)
              Call zgemm_('C','N',nEXCH,nEXCH,nEXCH,
     &                   (1.0_wp,0.0_wp),Z, nEXCH,
     &                                   M(L,:,:), nEXCH,
     &                   (0.0_wp,0.0_wp),TMP, nEXCH )
              Call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,M(L,:,:),1)
              Call zgemm_('N','N',nEXCH,nEXCH,nEXCH,
     &                   (1.0_wp,0.0_wp),TMP,nEXCH,
     &                                     Z,nEXCH,
     &                   (0.0_wp,0.0_wp), M(L,:,:), nEXCH )

              Call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,TMP,1)
              Call zgemm_('C','N',nEXCH,nEXCH,nEXCH,
     &                   (1.0_wp,0.0_wp),Z,nEXCH,
     &                                   S(L,:,:), nEXCH,
     &                   (0.0_wp,0.0_wp),TMP,nEXCH )
              Call zcopy_(nexch*nexch,[(0.0_wp,0.0_wp)],0,S(L,:,:),1)
              Call zgemm_('N','N',nEXCH,nEXCH,nEXCH,
     &                   (1.0_wp,0.0_wp),TMP,nEXCH,
     &                                     Z,nEXCH,
     &                   (0.0_wp,0.0_wp),S(L,:,:),nEXCH )
            End Do  ! L

            Call mma_deallocate(SX2)
            Call mma_deallocate(SY2)
            Call mma_deallocate(SZ2)
            Call mma_deallocate(S2)
            Call mma_deallocate(W)
            Call mma_deallocate(HZFS)
         End If ! ZFS is defined

         Call mma_deallocate(MTMP)
         Call mma_deallocate(STMP)
         Call mma_deallocate(Z)
         Call mma_deallocate(tmp)
      End If ! dznrm2_ M and S

      If(dbg) Then
         Write(6,'(A)') 'g tensor at the end of GENERATE_SPIN'
         g=0.0_wp
         ma=0.0_wp
         Call atens(M, nExch, g, ma, 2)
      End If


      Call qExit('generate_isotrop_site')
      Return
      End subroutine generate_isotrop_site
