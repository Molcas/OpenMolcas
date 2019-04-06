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
      Subroutine exchctl( exch, nneq, neqv, neq, nexch, nmax,
     &                    lmax, npair, i_pair, MxRank1, MxRank2,
     &                    imaxrank,
     &                    Jex, JAex, JAex9,
     &                    JDMex, JITOexR, JITOexI, eso, SM,  MM,
     &                    coord, rot, rlg, riso,
     &                    tpar, upar, lant, itype,
     &                    Dipol, AnisoLines1, AnisoLines3, AnisoLines9,
     &                    KE, KEOPT, DM_exchange, JITO_exchange,
     &                    W,  Z,  S,   M, iPrint, mem  )


!     this Subroutine is a control Subroutine for the exchange interaction,
!     diagonalization of total hamiltonian and computation of matrix elements
!     of magnetic and spin moment
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
      ! global variables:
      Integer, intent(in)           :: nneq
      Integer, intent(in)           :: neqv ! max of neq(nneq)
      Integer, intent(in)           :: neq(nneq)
      Integer, intent(in)           :: nexch(nneq)
      Integer, intent(in)           :: nmax
      Integer, intent(in)           :: lmax
      Integer, intent(in)           :: npair
      Integer, intent(in)           :: i_pair(npair,2)
      Integer, intent(in)           :: exch
!     ( takes values from 1-7 for Gd-Yb respectively
      Integer, intent(in)           :: lant
      Integer, intent(in)           :: iPrint
      Integer, intent(in)           :: mem ! memory allocated so far
      Integer, intent(in)           :: MxRank1, MxRank2
      Integer, intent(in)           :: imaxrank(npair,2)
      Character(1), intent(in)      :: itype(nneq)

      Real(kind=wp), intent(in)     :: eso(nneq,nmax)
      Real(kind=wp), intent(in)     :: Jex(npair)
      Real(kind=wp), intent(in)     :: JAex(npair,3)
      Real(kind=wp), intent(in)     :: JDMex(npair,3)
      Real(kind=wp), intent(in)     :: JAex9(npair,3,3)
      Real(kind=wp), intent(in)     ::
     &                          JITOexR(nPair,MxRank1,-MxRank1:MxRank1,
     &                                        MxRank2,-MxRank2:MxRank2)
      Real(kind=wp), intent(in)     ::
     &                          JITOexI(nPair,MxRank1,-MxRank1:MxRank1,
     &                                        MxRank2,-MxRank2:MxRank2)
      Real(kind=wp), intent(in)     :: coord(nneq,3)
      Real(kind=wp), intent(in)     :: rot(nneq,neqv,3,3)
      Real(kind=wp), intent(in)     :: rlg(nneq,neqv,3,3)
      Real(kind=wp), intent(in)     :: riso(nneq,3,3)
      Real(kind=wp), intent(in)     :: tpar
      Real(kind=wp), intent(in)     :: upar

      Complex(kind=wp), intent(inout)  :: SM(nneq,3,nmax,nmax)
      Complex(kind=wp), intent(inout)  :: MM(nneq,3,nmax,nmax)

      Logical, intent(in)           :: AnisoLines1
      Logical, intent(in)           :: AnisoLines3
      Logical, intent(in)           :: AnisoLines9
      Logical, intent(in)           :: Dipol
      Logical, intent(in)           :: KE
      Logical, intent(in)           :: DM_exchange
      Logical, intent(in)           :: JITO_exchange




      Real(kind=wp), intent(out)    :: W(exch)

      Complex(kind=wp), intent(out) :: Z(exch,exch)
      Complex(kind=wp), intent(out) :: S(3,exch,exch)
      Complex(kind=wp), intent(out) :: M(3,exch,exch)
!------------------------------------------------------------------
      ! local variables
      Integer                       :: i,j,l,lp,lb1,lb2,nb,nb1,nb2,
     &                                 isite,is1,js1,i1,i2,j1,j2,k,
     &                                 ibuf,q1,q2,k1,k2,n1,n2
      Integer                       :: norder
      Integer                       :: NmaxPop
      Integer, allocatable          :: intc(:)   !  intc(lmax)
      Integer, allocatable          :: ibas(:,:) !  ibas(exch,lmax)
      Integer, allocatable          :: icoord(:) !  icoord(lmax)
      Integer, allocatable          :: nind(:,:) !  nind(lmax,2)

      Real(kind=wp)                 :: vect(3)
      Real(kind=wp)                 :: dist
      Real(kind=wp), allocatable    :: wlin(:) ! wlin(exch)
      Real(kind=wp), allocatable    :: wlin1(:)! wlin1(exch)
      Real(kind=wp), allocatable    :: wlin3(:)! wlin3(exch)
      Real(kind=wp), allocatable    :: wlin9(:)! wlin9(exch)
      Real(kind=wp), allocatable    :: wdip(:) ! wdip(exch)
      Real(kind=wp), allocatable    :: wkex(:) ! wkex(exch)
      Real(kind=wp), allocatable    :: wdmo(:) ! wdmo(exch)
      Real(kind=wp), allocatable    :: wito(:) ! wito(exch)

      Complex(kind=wp), allocatable :: S1(:,:,:) ! S1(3,nmax,nmax)
      Complex(kind=wp), allocatable :: M1(:,:,:) ! M1(3,nmax,nmax)
      Complex(kind=wp), allocatable :: S2(:,:,:) ! S2(3,nmax,nmax)
      Complex(kind=wp), allocatable :: M2(:,:,:) ! M2(3,nmax,nmax)
      Complex(kind=wp), allocatable :: ZA1(:,:), ZA2(:,:)
      Complex(kind=wp), allocatable :: SM1(:,:,:) ! SM1(3,nmax,nmax)
      Complex(kind=wp), allocatable :: MM1(:,:,:) ! MM1(3,nmax,nmax)
      Complex(kind=wp), allocatable :: SM2(:,:,:) ! SM2(3,nmax,nmax)
      Complex(kind=wp), allocatable :: MM2(:,:,:) ! MM2(3,nmax,nmax)
      Complex(kind=wp), allocatable :: HLIN1(:,:,:,:,:)
!                                      HLIN1(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HLIN3(:,:,:,:,:)
!                                      HLIN3(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HLIN9(:,:,:,:,:)
!                                      HLIN9(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HDIP(:,:,:,:,:)
!                                      HDIP(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HKEX(:,:,:,:,:)
!                                      HKEX(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HDMO(:,:,:,:,:)
!                                      HDMO(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: HITO(:,:,:,:,:)
!                                      HITO(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), allocatable :: tmp(:,:) ! tmp(exch,exch)
c two options for KE:
      Integer                       :: KEOPT
      Integer, parameter            :: exchR=8
      Integer                       :: nmaxR
      Integer                       :: nsta
      Integer, allocatable          :: nexchR(:)  ! nexchR(nneq)
      Integer, allocatable          :: ibasR(:,:) ! ibasR(exchR,lmax)
      Integer, allocatable          :: intcR(:)   ! intcR(lmax)
      Real(kind=wp), allocatable    :: WR(:)  ! WR(exchR)
      Real(kind=wp), allocatable    :: rotR(:,:,:,:)
!                                      rotR(nneq,neqv,3,3)
      Complex(kind=wp), allocatable :: ZR(:,:) ! ZR(exchR,exchR)
      Complex(kind=wp), allocatable :: HKEXR(:,:,:,:,:)
!                                      HKEXR(npair,2,2,2,2)
      Complex(kind=wp), allocatable :: MR(:,:,:) ! MR(3,exchR,exchR)
      Complex(kind=wp), allocatable :: SR(:,:,:) ! SR(3,exchR,exchR)
      Complex(kind=wp), allocatable :: SMR(:,:,:,:) ! SMR(nneq,3,2,2)
      Complex(kind=wp), allocatable :: MMR(:,:,:,:) ! MMR(nneq,3,2,2)
c      Complex(kind=wp) ::  JAllDip(npair,nmax,-nmax:nmax,nmax,-nmax:nmax)
c      Complex(kind=wp) ::  JAllEx( npair,nmax,-nmax:nmax,nmax,-nmax:nmax)
c      Real(kind=wp) ::  J1Dip(npair,3,3)
c      Real(kind=wp) ::  J1Ex(npair,3,3)
      Real(kind=wp)    :: mg1(3,3), mg2(3,3)
      Integer          :: CtoB, RtoB, ItoB, mem_local
      Logical          :: DBG !, testlines
      Real(kind=wp)    :: dnrm2_
      External         :: norder, dnrm2_  !,ilaenv

      Call qEnter('PA_exchctl')
      DBG=.false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      If(DBG) Then
         Write(6,'(A)') 'Enter EXCHCTL'
         Write(6,'(A,  i8)') 'exch   = ',exch
         Write(6,'(A,  i8)') 'nneq   = ',nneq
         Write(6,'(A,  i8)') 'neqv   = ',neqv
         Write(6,'(A,  i8)') 'nmax   = ',nmax
         Write(6,'(A,  i8)') 'lmax   = ',lmax
         Write(6,'(A,  i8)') 'nPair  = ',nPair
         Write(6,'(A,  i8)') 'nPair  = ',nPair
         Write(6,'(A,  i8)') 'MxRank1= ',MxRank1
         Write(6,'(A,  i8)') 'MxRank2= ',MxRank2
         Write(6,'(A,10i4)') 'neq()  = ',(neq(i),i=1,nneq)
         Write(6,'(A,10i4)') 'nexch()= ',(nexch(i),i=1,nneq)
         Write(6,'(A,  L2)') 'AnisoLines1  = ',AnisoLines1
         Write(6,'(A,  L2)') 'AnisoLines3  = ',AnisoLines3
         Write(6,'(A,  L2)') 'AnisoLines9  = ',AnisoLines9
         Write(6,'(A,  L2)') 'DM_exchange  = ',DM_exchange
         Write(6,'(A,  L2)') 'Dipol        = ',Dipol
         Write(6,'(A,  L2)') 'JITO_exchange= ',JITO_exchange
         If ( AnisoLines1 ) Then
            Do i=1, nPair
              Write(6,'(A,2I4,F10.5)') 'LIN1', i_pair(i,1), i_pair(i,2),
     &                                Jex(i)
            End Do
         End If

         If ( AnisoLines3 ) Then
           Do i=1, nPair
             Write(6,'(A,2I4,3F10.5)') 'LIN3', i_pair(i,1), i_pair(i,2),
     &                              (JAex(i,j),j=1,3)
           End Do
         End If

         If( AnisoLines9 ) Then
           Do i=1, nPair
             Write(6,'(A,2I4,9F10.5)') 'LIN9', i_pair(i,1), i_pair(i,2),
     &                              ((JAex9(i,j,k),j=1,3),k=1,3)
           End Do
         End If

         If(DM_exchange) Then
            Do i=1,nPair
           Write(6,'(A,2I4,3F10.5)') 'DMEX', i_pair(i,1), i_pair(i,2),
     &                              (JDMex(i,j),j=1,3)
            End Do
         End If

         If(Dipol) Then
         Write(6,'(A)') 'COORD(i):'
           Do i=1,nneq
              Write(6,'(i3,3F14.8)') i, (coord(i,j),j=1,3)
           End Do
         End If

         If(JITO_exchange) Then
            Do i=1,nPair
           Write(6,'(A,4I4,3F10.5)') 'ITOJ', i_pair(i,1), i_pair(i,2),
     &                          imaxrank(i,1),imaxrank(i,2)
              Do k1=1,imaxrank(i,1),2
               Do q1=-k1,k1
                Do k2=1,imaxrank(i,2),2
                 Do q2=-k2,k2
                  Write(6,'(4I3,2x,2ES21.14)') k1,q1,k2,q2,
     &                JITOexR(i,k1,q1,k2,q2), JITOexI(i,k1,q1,k2,q2)
                 End Do
                End Do
               End Do
              End Do
            End Do ! ipair
         End If

         Write(6,'(A)') 'ESO(i):'
         Do i=1,nneq
            Write(6,'(i3,90F14.8)') i, (eso(i,j),j=1,nmax)
         End Do
         Write(6,'(90A)') 'itype()=', (itype(i),' ',i=1,nneq)

         Write(6,'(A,  i8)') 'neqv   = ',neqv
         Do i=1, nneq
         Write(6,'(A,  i8)') 'riso( site=',i,'):'
            Do j=1,3
              Write(6,'(3ES22.14)') (riso(i,j,k),k=1,3)
            End Do
         End Do
      End If ! dbg
!-----------------------------------------------------------------------
! allocate memory for this function:
      ItoB=8
      RtoB=8
      CtoB=16
      mem_local=0
      If(lmax>0) Then
        ! exchange energy spectrum
        Call mma_allocate(intc,lmax,'intc')
        Call mma_allocate(icoord,lmax,'icoord')
        Call mma_allocate(nind,lmax,2,'nind')
        Call icopy(lmax,[0],0,intc,1)
        Call icopy(lmax,[0],0,icoord,1)
        Call icopy(2*lmax,[0],0,nind,1)
        mem_local=mem_local+4*lmax*ItoB
        If(exch>0) Then
          Call mma_allocate(ibas,exch,lmax,'ibas')
          Call icopy(exch*lmax,[0],0,ibas,1)
          mem_local=mem_local+exch*lmax*ItoB
        End If
      End If
      If(exch>0) Then
        Call mma_allocate(wlin ,exch,'wlin ')
        Call mma_allocate(wlin1,exch,'wlin1')
        Call mma_allocate(wlin3,exch,'wlin3')
        Call mma_allocate(wlin9,exch,'wlin9')
        Call mma_allocate(wdip ,exch,'wdip ')
        Call mma_allocate(wkex ,exch,'wkex ')
        Call mma_allocate(wdmo ,exch,'wdmo ')
        Call mma_allocate(wito ,exch,'wito ')
        Call dcopy_(exch,[0.0_wp],0,wlin ,1)
        Call dcopy_(exch,[0.0_wp],0,wlin1,1)
        Call dcopy_(exch,[0.0_wp],0,wlin3,1)
        Call dcopy_(exch,[0.0_wp],0,wlin9,1)
        Call dcopy_(exch,[0.0_wp],0,wdip ,1)
        Call dcopy_(exch,[0.0_wp],0,wkex ,1)
        Call dcopy_(exch,[0.0_wp],0,wdmo ,1)
        Call dcopy_(exch,[0.0_wp],0,wito ,1)
        mem_local=mem_local+8*exch*RtoB
      End If
      If(nmax>0) Then
        Call mma_allocate( S1,3,nmax,nmax,' S1')
        Call mma_allocate( M1,3,nmax,nmax,' M1')
        Call mma_allocate( S2,3,nmax,nmax,' S2')
        Call mma_allocate( M2,3,nmax,nmax,' M2')
        Call mma_allocate( ZA1, nmax,nmax,' Z1')
        Call mma_allocate( ZA2, nmax,nmax,' Z2')
        Call mma_allocate(SM1,3,nmax,nmax,'SM1')
        Call mma_allocate(SM2,3,nmax,nmax,'SM2')
        Call mma_allocate(MM1,3,nmax,nmax,'MM1')
        Call mma_allocate(MM2,3,nmax,nmax,'MM2')
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S1,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, M1,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S2,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, M2,1)
        Call zcopy_(  nmax*nmax,[(0.0_wp,0.0_wp)],0, ZA1,1)
        Call zcopy_(  nmax*nmax,[(0.0_wp,0.0_wp)],0, ZA2,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,SM1,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,SM2,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,MM1,1)
        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,MM2,1)
        mem_local=mem_local+8*3*nmax*nmax*CtoB

        If(npair>0) Then
          ibuf=npair*nmax*nmax*nmax*nmax
          Call mma_allocate(HLIN1,npair,nmax,nmax,nmax,nmax,'HLIN1')
          Call mma_allocate(HLIN3,npair,nmax,nmax,nmax,nmax,'HLIN3')
          Call mma_allocate(HLIN9,npair,nmax,nmax,nmax,nmax,'HLIN9')
          Call mma_allocate(HDIP,npair,nmax,nmax,nmax,nmax,'HDIP')
          Call mma_allocate(HKEX,npair,nmax,nmax,nmax,nmax,'HKEX')
          Call mma_allocate(HDMO,npair,nmax,nmax,nmax,nmax,'HDMO')
          Call mma_allocate(HITO,npair,nmax,nmax,nmax,nmax,'HITO')
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN1,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN3,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN9,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDIP,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HKEX,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDMO,1)
          Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HITO,1)
          mem_local=mem_local+7*ibuf*CtoB
        End If
      End If

      If(exch>0) Then
        Call mma_allocate(tmp,exch,exch,'tmp')
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,tmp,1)
        mem_local=mem_local+exch*exch*CtoB
      End If


      If(nneq>0) Then
        Call mma_allocate(nexchR,nneq,'nexchR')
        Call icopy( nneq,[0],0,nexchR,1)
        mem_local=mem_local+nneq*ItoB

        Call mma_allocate(SMR,nneq,3,2,2,'SMR')
        Call mma_allocate(MMR,nneq,3,2,2,'MMR')
        Call zcopy_(nneq*3*2*2,[(0.0_wp,0.0_wp)],0,SMR,1)
        Call zcopy_(nneq*3*2*2,[(0.0_wp,0.0_wp)],0,MMR,1)
        mem_local=mem_local+2*nneq*3*2*2*CtoB

        If(neqv>0) Then
          Call mma_allocate(rotR,nneq,neqv,3,3,'rotR')
          Call dcopy_(nneq*neqv*3*3,[0.0_wp],0,rotR,1)
          mem_local=mem_local+nneq*neqv*3*3*RtoB
        End If
      End If

      If(exchR>0) Then
        If(lmax>0) Then
          Call mma_allocate(ibasR,nneq,lmax,'ibasR')
          Call icopy( nneq*lmax,[0],0,ibasR,1)
          mem_local=mem_local+nneq*lmax*ItoB
        End If
        Call mma_allocate(WR,exchR,'WR')
        Call dcopy_(exchR,[0.0_wp],0,WR,1)
        mem_local=mem_local+exchR*RtoB

        Call mma_allocate(ZR,exchR,exchR,'ZR')
        Call mma_allocate(MR,3,exchR,exchR,'MR')
        Call mma_allocate(SR,3,exchR,exchR,'SR')
        Call zcopy_(  exchR*exchR,[(0.0_wp,0.0_wp)],0,ZR,1)
        Call zcopy_(3*exchR*exchR,[(0.0_wp,0.0_wp)],0,MR,1)
        Call zcopy_(3*exchR*exchR,[(0.0_wp,0.0_wp)],0,SR,1)
        mem_local=mem_local+7*exchR*exchR*CtoB
      End If

      If(npair>0) Then
        Call mma_allocate(HKEXR,npair,2,2,2,2,'HKEXR')
        Call zcopy_(npair*2*2*2*2,[(0.0_wp,0.0_wp)],0,HKEXR,1)
        mem_local=mem_local+npair*2*2*2*2*CtoB
      End If

      If(lmax>0) Then
        Call mma_allocate(intcR,lmax,'intcR')
        Call icopy(lmax,[0],0,intcR,1)
        mem_local=mem_local+lmax*ItoB
      End If
      If(dbg) Write(6,*) 'EXCHCTL:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'EXCHCTL:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local
!-----------------------------------------------------------------------
      l=0
      Do i=1,nneq
        Do j=1,neq(i)
          l=l+1
          nind(l,1)=i
          nind(l,2)=j
        End Do
      End Do
      intc(1)=1
      If (lmax.gt.1) Then
        Do i=2,lmax
          isite=nind(i-1, 1)
          intc(i)=intc(i-1)*nexch(isite)
        End Do
      End If
      Do nb=1,exch
        nb1=nb-1
        Do i=1,lmax
          ibas(nb, lmax-i+1)= nb1 / intc(lmax-i+1)
          nb1=nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
        End Do
      End Do
      If(dbg) Then
        Write(6,'(34x,A,1x,20i3)')  'site Nr.', (i,i=1,lmax)
        Do nb=1,exch
          Write(6,'(A,i5,A,20i3)') 'COUPLING: basis set:  ibas(',nb,
     &                             ' ,isite) = ',(ibas(nb,i)+1,i=1,lmax)
        End Do
      End If ! dbg

!-----------------------------------------------------------------------
! Lines model of magnetic couping  -- 1 parameter
      If ( AnisoLines1 ) Then
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN1,1)
       Do lp=1, npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)

        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, .true. )
        End If
        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, .true. )
        End If

        If(dbg) Call prMom('SM(i1) bf Lines1',SM(i1,1:3,1:n1,1:n1),n1)
        If(dbg) Call prMom('SM(i2) bf Lines1',SM(i2,1:3,1:n2,1:n2),n2)

        ! build the Lines exchange matrix:
        Call Lines_Exchange( Jex(lp), n1, n2,
     &                       SM(i1,1:3,1:n1,1:n1),
     &                       SM(i2,1:3,1:n2,1:n2),
     &              HLIN1(lp,1:n1,1:n1,1:n2,1:n2) )

        If(dbg) Then
         Call prMom('SM(i1) af Lines1',SM(i1,1:3,1:n1,1:n1),n1)
         Call prMom('SM(i2) af Lines1',SM(i2,1:3,1:n2,1:n2),n2)
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Write(6,'(A,i2)') 'in local pseudospin basis'
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HLIN1(',i,',',j,',',k,',',l,') = ',HLIN1(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If! dbg

       End Do
       End If ! nPair
      End If
!-----------------------------------------------------------------------





      ! Anisotropic Lines model of magnetic couping -- 3 parameters)
      ! Jxx, Jyy, Jzz
      If (AnisoLines3 ) Then
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN3,1)
       Do lp=1, npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)

        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, dbg )
        End If

        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, dbg )
        End If

        Call Aniso_Lines_Exchange3( JAex(lp,1:3), n1, n2,
     &                              SM(i1,1:3,1:n1,1:n1),
     &                              SM(i2,1:3,1:n2,1:n2),
     &                     HLIN3(lp,1:n1,1:n1,1:n2,1:n2) )

        If(dbg) Then
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HLIN3(',i,',',j,',',k,',',l,') = ',HLIN3(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If ! dbg

       End Do ! lp
       End If ! nPair
      End If
!-----------------------------------------------------------------------




      ! Anisotropic Lines model of magnetic couping -- 9 parameters)
      ! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
      If ( AnisoLines9 ) Then
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN9,1)
       Do lp=1, npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)

        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, dbg )
        End If
        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, dbg )
        End If
        ! rotate the input J matrix by the mg1 and mg2
        If(dbg) Call prMom('SM(i1) bf Lines9',SM(i1,1:3,1:n1,1:n1),n1)
        If(dbg) Call prMom('SM(i2) bf Lines9',SM(i2,1:3,1:n2,1:n2),n2)


        Call Aniso_Lines_Exchange9( JAex9(lp,1:3,1:3), n1, n2,
     &                              SM(i1,1:3,1:n1,1:n1),
     &                              SM(i2,1:3,1:n2,1:n2),
     &                     HLIN9(lp,1:n1,1:n1,1:n2,1:n2) )


        If(dbg) Then
         Call prMom('SM(i1) af Lines9',SM(i1,1:3,1:n1,1:n1),n1)
         Call prMom('SM(i2) af Lines9',SM(i2,1:3,1:n2,1:n2),n2)
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HLIN9(',i,',',j,',',k,',',l,') = ',HLIN9(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If! dbg

       End Do
       End If ! nPair
      End If



!-----------------------------------------------------------------------
!     dipolar couping
      If(Dipol) Then
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDIP,1)
       Do lp=1,npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)

        vect(1:3)=0.0_wp
        dist=0.0_wp
        Call dirvect( coord(i1,1:3), rlg(i1,j1,1:3,1:3),
     &                coord(i2,1:3), rlg(i2,j2,1:3,1:3),
     &                vect(1:3), dist )
        If(DBG) Write(6,'(A,i3,3ES20.12,2x,ES20.12)')
     &                  'EXCHCTL: DIPOL: lp, vect, R:',
     &                                   lp,(vect(i),i=1,3), dist
        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, dbg )
        End If
        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, dbg )
        End If

        Call Dipol_Exchange( n1, n2, vect(1:3), dist,
     &                       MM(i1,1:3,1:n1,1:n1),
     &                       MM(i2,1:3,1:n2,1:n2),
     &               HDIP(lp,1:n1,1:n1,1:n2,1:n2) )

        If(dbg) Then
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HDIP (',i,',',j,',',k,',',l,') = ',HDIP(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If ! dbg

       End Do ! lp
       End If
      End If
!-----------------------------------------------------------------------





!     Dzyaloshinsky-Morya antisymmetric couping
      If(DM_exchange) Then
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDMO,1)
       Do lp=1,npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)
        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, dbg )
        End If
        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, dbg )
        End If

        Call Dzyaloshinsky_Morya_Exchange( JDMex(lp,1:3), n1, n2,
     &                                     SM(i1,1:3,1:n1,1:n1),
     &                                     SM(i2,1:3,1:n2,1:n2),
     &                             HDMO(lp,1:n1,1:n1,1:n2,1:n2) )
        If(dbg) Then
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HDMO (',i,',',j,',',k,',',l,') = ',HDMO(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If! dbg

       End Do
       End If ! nPair
      End If
!-----------------------------------------------------------------------




!     JITO exchange interaction
      If(JITO_exchange) Then
       If(dbg) Write(6,'(A)') 'EXCHCTL:  Enterring  JITO_exchange'
       If(nPair>0) Then
       Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HITO,1)
       Do lp=1,npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)
        ! find local pseudospin and rotate the spin and magnetic moment
        ! to the local pseudospin basis
        If(itype(i1)=='A') Then
         Call prep_mom_exchange( n1, rot(i1,j1,1:3,1:3),
     &                            SM(i1,1:3,1:n1,1:n1),
     &                            MM(i1,1:3,1:n1,1:n1), mg1, dbg )
        End If
        If(itype(i2)=='A') Then
         Call prep_mom_exchange( n2, rot(i2,j2,1:3,1:3),
     &                            SM(i2,1:3,1:n2,1:n2),
     &                            MM(i2,1:3,1:n2,1:n2), mg2, dbg )
        End If

        ! using Naoya's ITO:, in general coordinate system
        Call JITO_Exchange_Int( MxRank1, MxRank2, imaxrank(lp,1:2),
     &                          n1, n2,
     &                    JITOexR(lp,1:MxRank1,-MxRank1:MxRank1,
     &                               1:MxRank2,-MxRank2:MxRank2 ),
     &                    JITOexI(lp,1:MxRank1,-MxRank1:MxRank1,
     &                               1:MxRank2,-MxRank2:MxRank2 ),
     &                    HITO(lp,1:n1,1:n1,1:n2,1:n2) )
        If(dbg) Then
         Write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
         Do i=1,n1
          Do j=1,n1
           Do k=1,n2
            Do l=1,n2
              Write(6,'(4(a,i2),a,2ES24.14)')
     &         'HITO (',i,',',j,',',k,',',l,') = ',HITO(lp,i,j,k,l)
            End Do
           End Do
          End Do
         End Do
        End If! dbg

       End Do
       End If ! nPair
       If(dbg) Write(6,'(A)') 'EXCHCTL:  Exiting JITO_exchange'
      End If
!-----------------------------------------------------------------------











      If (KE) Then
        If (nPair>0) Then
        HKEX=(0.0_wp,0.0_wp)
        HKEX=(0.0_wp,0.0_wp)
        HKEXR=(0.0_wp,0.0_wp)
        MMR=(0.0_wp,0.0_wp)
        SMR=(0.0_wp,0.0_wp)
        Do lp=1,npair
          lb1=i_pair(lp,1)
          lb2=i_pair(lp,2)
          i1=nind(lb1,1) ! indices of non-equivalent sites
          i2=nind(lb2,1) ! indices of non-equivalent sites
          j1=nind(lb1,2) ! indices of equivalent sites
          j2=nind(lb2,2) ! indices of equivalent sites

          n1=nexch(i1)
          n2=nexch(i2)
          Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0, S1,1)
          Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0, S2,1)
          Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0, M1,1)
          Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0, M2,1)
          Call rotmom2( MM(i1,1:3,1:n1,1:n1), n1, rot(i1,j1,1:3,1:3),
     &                  M1(1:3,1:n1,1:n1) )
          Call rotmom2( SM(i1,1:3,1:n1,1:n1), n1, rot(i1,j1,1:3,1:3),
     &                  S1(1:3,1:n1,1:n1) )
          Call rotmom2( MM(i2,1:3,1:n2,1:n2), n2, rot(i2,j2,1:3,1:3),
     &                  M2(1:3,1:n2,1:n2) )
          Call rotmom2( SM(i2,1:3,1:n2,1:n2), n2, rot(i2,j2,1:3,1:3),
     &                  S2(1:3,1:n2,1:n2) )
          ! KEOPT=1 ! FULL
          ! KEOPT=2 ! Full, 1/U
          ! KEOPT=3 ! FULL  + reduced form
          ! KEOPT=4 ! Full, 1/U + reduced form
          MM1=(0.0_wp,0.0_wp)
          SM1=(0.0_wp,0.0_wp)
          MM2=(0.0_wp,0.0_wp)
          SM2=(0.0_wp,0.0_wp)
!IFG: this call to KE_Exchange does not match at all its definition, please fix
          Call WarningMessage(2,'There is surely a bug here')
          If (.False.) Call Unused_real(tpar)
          If (.False.) Call Unused_real(upar)
          If (.False.) Call Unused_integer(lant)
!         Call KE_Exchange(n1,n2,
!    &                     M1( 1:3, 1:n1, 1:n1 ),
!    &                     S1( 1:3, 1:n1, 1:n1 ),
!    &                     M2( 1:3, 1:n2, 1:n2 ),
!    &                     S2( 1:3, 1:n2, 1:n2 ),
!    &                     eso(i1,1:n1),eso(i2,1:n2),
!    &                     tpar, upar, lant, KEOPT,
!    &       HKEX(lp,1:n1,1:n1,1:n2,1:n2),
!    &                     MM1(1:3, 1:n1, 1:n1 ),
!    &                     SM1(1:3, 1:n1, 1:n1 ),
!    &                     MM2(1:3, 1:n2, 1:n2 ),
!    &                     SM2(1:3, 1:n2, 1:n2 ) )

          If((KEOPT.eq.1).OR.(KEOPT.eq.2)) Then
            Do is1=1,n2
              Do js1=1,n2
                HKEXR(lp,1,1,is1,js1)=HKEX(lp, 1, 1,is1,js1)
                HKEXR(lp,1,2,is1,js1)=HKEX(lp, 1, 2,is1,js1)
                HKEXR(lp,2,1,is1,js1)=HKEX(lp, 2, 1,is1,js1)
                HKEXR(lp,2,2,is1,js1)=HKEX(lp, 2, 2,is1,js1)
              End Do
            End Do
            Do l=1,3
              Do is1=1,2
                Do js1=1,2
                  MMR(i1,l,is1,js1)=MM1(l,is1,js1)
                  MMR(i2,l,is1,js1)=MM2(l,is1,js1)
                  SMR(i1,l,is1,js1)=SM1(l,is1,js1)
                  SMR(i2,l,is1,js1)=SM2(l,is1,js1)
                End Do
              End Do
            End Do

          Else If((KEOPT.eq.3).OR.(KEOPT.eq.4)) Then
            Do is1=1,n2
              Do js1=1,n2
                HKEXR(lp,1,1,is1,js1)=HKEX(lp, 1, 1,is1,js1)
                HKEXR(lp,1,2,is1,js1)=HKEX(lp, 1,n1,is1,js1)
                HKEXR(lp,2,1,is1,js1)=HKEX(lp,n1, 1,is1,js1)
                HKEXR(lp,2,2,is1,js1)=HKEX(lp,n1,n1,is1,js1)
              End Do
            End Do

            Do l=1,3
              MMR(i1,l,1,1)=MM1(l, 1, 1)
              MMR(i1,l,1,2)=MM1(l, 1,n1)
              MMR(i1,l,2,1)=MM1(l,n1, 1)
              MMR(i1,l,2,2)=MM1(l,n1,n1)

              MMR(i2,l,1,1)=MM2(l, 1, 1)
              MMR(i2,l,1,2)=MM2(l, 1,n2)
              MMR(i2,l,2,1)=MM2(l,n2, 1)
              MMR(i2,l,2,2)=MM2(l,n2,n2)

              SMR(i1,l,1,1)=SM1(l, 1, 1)
              SMR(i1,l,1,2)=SM1(l, 1,n1)
              SMR(i1,l,2,1)=SM1(l,n1, 1)
              SMR(i1,l,2,2)=SM1(l,n1,n1)

              SMR(i2,l,1,1)=SM2(l, 1, 1)
              SMR(i2,l,1,2)=SM2(l, 1,n2)
              SMR(i2,l,2,1)=SM2(l,n2, 1)
              SMR(i2,l,2,2)=SM2(l,n2,n2)
            End Do

            If(DBG) Then
              Write(6,'(A)')'site Ln'
              Do i=1,2
                Do j=1,2
                  Write(6,'(3(A,i1),A,3(2F20.14,2x))')
     &              'MMR(',i1,',L,',i,',',j,')= ',(MMR(i1,l,i,j),l=1,3)
                End Do
              End Do
              Write(6,'(/)')
              Do i=1,2
                Do j=1,2
                  Write(6,'(3(A,i1),A,3(2F20.14,2x))')
     &              'SMR(',i1,',L,',i,',',j,')= ',(SMR(i1,l,i,j),l=1,3)
                End Do
              End Do
              Write(6,'(/)')
              Write(6,'(A)')'site Radical'
              Do i=1,2
                Do j=1,2
                  Write(6,'(3(A,i1),A,3(2F20.14,2x))')
     &              'MMR(',i2,',L,',i,',',j,')= ',(MMR(i2,l,i,j),l=1,3)
                End Do
              End Do
              Write(6,'(/)')
              Do i=1,2
                Do j=1,2
                  Write(6,'(3(A,i1),A,3(2F20.14,2x))')
     &              'SMR(',i2,',L,',i,',',j,')= ',(SMR(i2,l,i,j),l=1,3)
                End Do
              End Do
            End If ! DBG

          End If ! KEOPT
        End Do ! lp

!       in case of KEOPT=3 or KEOPT=4 Then we need to compute the spectrum and the properties
!       in the reduced form, where nexch(i1)=2 ( ground Doublet only).
!       exchnew=8:
        If(KEOPT.le.4) Then
          nmaxR=2
          nexchR(1)=2
          nexchR(2)=2
          Do i=1,nneq
            nexchR(i)=2
          End Do
          intcR(:)=0
          ibasR(:,:)=0
          intcR(1)=1
          If (lmax.gt.1) Then
            Do i=2,lmax
              isite=nind(i-1, 1)
              intcR(i)=intcR(i-1)*nexchR(isite)
            End Do
          End If
          Do nb=1,exchR
            nb1=nb-1
            Do i=1,lmax
              ibasR(nb, lmax-i+1)= nb1 / intcR(lmax-i+1)
              nb1=nb1-ibasR(nb,lmax-i+1)*intcR(lmax-i+1)
            End Do
          End Do
cccc----------------------------------------------------------------------------------
          HLIN1=(0.0_wp,0.0_wp)
          HLIN3=(0.0_wp,0.0_wp)
          HLIN9=(0.0_wp,0.0_wp)
          HDIP=(0.0_wp,0.0_wp)
          WLIN=0.0_wp
          WLIN1=0.0_wp
          WLIN3=0.0_wp
          WLIN9=0.0_wp
          WDIP=0.0_wp
          WKEX=0.0_wp
          WR=0.0_wp
          ZR=(0.0_wp,0.0_wp)
          ! print the Exchange Hamiltonian:
          Call pa_prham( exchR, npair, i_pair,  nneq,   neq, nexchR,
     &                   nmaxR,
     &                lmax, eso(1:nneq,1:nmaxR),
     &                HLIN1(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HLIN3(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HLIN9(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HDIP(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &               HKEXR(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HDMO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HITO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                Dipol, AnisoLines1, AnisoLines3, AnisoLines9, KE,
     &                DM_exchange, .False. )
          ! diagonalize the Hamiltonian:
          Call pa_diagham( exchR, npair, i_pair, nneq, neq, nexchR,
     &                     nmaxR,
     &                lmax, eso(1:nneq,1:nmaxR),
     &                HLIN1(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HLIN3(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HLIN9(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HDIP(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &               HKEXR(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HDMO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                HITO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),
     &                Dipol, .False.,
     &                AnisoLines1, AnisoLines3, AnisoLines9, KE,.False.,
     &                WLIN1(1:exchR), WLIN3(1:exchR), WLIN9(1:exchR),
     &                WLIN(1:exchR), WDIP(1:exchR), WKEX(1:exchR),
     &                WDMO(1:exchR), WITO(1:exchR),
     &                  WR(1:exchR), ZR(1:exchR,1:exchR) )
          ! print the resulting eigenstates:
          Call pa_preigen( exchR,  lmax,   ibasR,  Dipol,
     &                     AnisoLines1, AnisoLines3, AnisoLines9, KE,
     &                     .False., .False.,
     &                     WLIN(1:exchR), WDIP(1:exchR),
     &                     WKEX(1:exchR), WDMO(1:exchR), WITO(1:exchR),
     &                     WR(1:exchR), ZR(1:exchR,1:exchR), 0 )
          ! compute the moments:
          rotR=0.0_wp
          rotR(1,1,1,1)=1.0_wp
          rotR(1,1,2,2)=1.0_wp
          rotR(1,1,3,3)=1.0_wp
          rotR(1,2,1,1)=1.0_wp
          rotR(1,2,2,2)=1.0_wp
          rotR(1,2,3,3)=1.0_wp
          rotR(2,1,1,1)=1.0_wp
          rotR(2,1,2,2)=1.0_wp
          rotR(2,1,3,3)=1.0_wp
          MR=(0.0_wp,0.0_wp)
          SR=(0.0_wp,0.0_wp)
          Do L=1,3
            Do isite=1,lmax
              Do nb1=1,exchR
                Do lp=1,lmax
                  icoord(lp)=ibasR(nb1,lp)
                End Do
                i1=nind(isite,1)
                j1=nind(isite,2)
               is1=ibasR(nb1,isite)+1

                Do js1=1,nexchR(i1)
                  icoord(isite)=js1-1
                  nb2=norder(icoord,intcR,lmax)
                  MR( l, nb1, nb2 ) = MR( l, nb1, nb2 ) +
     &                  rotR( i1, j1, l, 1 ) * MMR( i1, 1, is1, js1 )
     &                 +rotR( i1, j1, l, 2 ) * MMR( i1, 2, is1, js1 )
     &                 +rotR( i1, j1, l, 3 ) * MMR( i1, 3, is1, js1 )
                  SR( l, nb1, nb2 ) = SR( l, nb1, nb2 ) +
     &                  rotR( i1, j1, l, 1 ) * SMR( i1, 1, is1, js1 )
     &                 +rotR( i1, j1, l, 2 ) * SMR( i1, 2, is1, js1 )
     &                 +rotR( i1, j1, l, 3 ) * SMR( i1, 3, is1, js1 )

                End Do  ! js1
              End Do  ! nb1
            End Do  ! isite
            TMP(:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,
     &                  (1.0_wp,0.0_wp), ZR,EXCHR,
     &                                   MR(L,:,:), EXCHR,
     &                  (0.0_wp,0.0_wp), TMP, EXCHR )
            MR(L,:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,
     &                  (1.0_wp,0.0_wp),TMP,EXCHR,
     &                                   ZR, EXCHR,
     &                  (0.0_wp,0.0_wp), MR(L,:,:), EXCHR )
            TMP(:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,
     &                  (1.0_wp,0.0_wp), ZR, EXCHR,
     &                                   SR(L,:,:), EXCHR,
     &                  (0.0_wp,0.0_wp), TMP, EXCHR )
            SR(L,:,:)=(0.0_wp,0.0_wp)
            Call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,
     &                  (1.0_wp,0.0_wp),TMP, EXCHR,
     &                                   ZR, EXCHR,
     &                  (0.0_wp,0.0_wp), SR(L,:,:), EXCHR )
            Do is1=1,8
              Do js1=1,8
                Write(6,'(3(A,i1),A,2F20.14)') 'MR(',l,',',
     &                     is1,',',js1,') = ',MR(l,is1,js1)
              End Do
            End Do
          End Do  ! L
          ! print the localized moments on sites:
          nsta=exchR
          ! assuming max 10 equivalent magnetic sites
          Call momloc2(nsta, nmaxR, nneq, neq, neqv,
     &                rotR(1:nneq,1:10,:,:), lmax, nexchR,
     &                wR(1:nsta),zR(1:nsta,1:nsta),
     &                MR(1:3,1:nsta,1:nsta),
     &                SR(1:3,1:nsta,1:nsta),
     &                MMR(1:nneq,1:3,1:nmaxR,1:nmaxR),
     &                SMR(1:nneq,1:3,1:nmaxR,1:nmaxR) )

          Call barrier(exchR,MR(1:3,1:exchR,1:exchR),WR(1:exchR),1,2)
        End If !KEOPT

        End If ! npair>0, index lp
      End If !KE


c--------------------------------------------------------------
c ALL exchange couplings for all pairs are now known.
c printout the Hamiltonians:
      If(DBG) Then
        Write(6,'(A,i5)') 'exch  = ',exch
        Write(6,'(A,i5)') 'npair = ',npair
        Write(6,'(A,i5)') 'nneq  = ',nneq
        Write(6,'(A,i5)') 'nmax  = ',nmax
        Write(6,'(A,i5)') 'lmax  = ',lmax
        Write(6,'(A,i5)') 'iPrint= ',iPrint
        Write(6,*) 'AnisoLines1 = ',AnisoLines1
        Write(6,*) 'AnisoLines3 = ',AnisoLines3
        Write(6,*) 'AnisoLines9 = ',AnisoLines9
        Write(6,*) 'Dipol = ',Dipol
        Write(6,*) 'JITO  = ',JITO_exchange
        Write(6,*) 'KE    = ',KE
        Do i=1,npair
          Write(6,'(A,i2,A,i3,3x,A,i2,A,i3)')
     &                       'i_pair(',i,',1) = ',i_pair(i,1),
     &                       'i_pair(',i,',2) = ',i_pair(i,2)
        End Do
        Do i=1,nneq
          Write(6,'(A,i2,A,i3)')    '   neq(',i,') = ',neq(i)
        End Do
        Do i=1,nneq
          Write(6,'(A,i2,A,i3)')    ' nexch(',i,') = ',neq(i)
        End Do
        Do i=1,nneq
          Write(6,'(A,i2,A,100F10.3)') '   eso(',i,') = ',
     &                                (eso(i,j),j=1,nexch(i))
        End Do
        Call xFlush(6)
      End If


      If((iPrint > 2).OR.(DBG)) Then
      ! print the Exchange Hamiltonian:
      Call pa_prham( exch, npair, i_pair, nneq,
     &               neq, nexch, nmax, lmax, eso,
     &               HLIN1, HLIN3, HLIN9, HDIP, HKEX, HDMO, HITO,
     &               Dipol, AnisoLines1, AnisoLines3, AnisoLines9,
     &               KE, DM_exchange, JITO_exchange )
      End If

      ! diagonalize the Hamiltonian:
      Call pa_diagham( exch, npair, i_pair,  nneq,   neq, nexch, nmax,
     &              lmax,   eso,   HLIN1, HLIN3, HLIN9,  HDIP,
     &              HKEX, HDMO, HITO, Dipol, DM_exchange, AnisoLines1,
     &              AnisoLines3, AnisoLines9, KE, JITO_exchange,
     &              WLIN1, WLIN3, WLIN9, WLIN, WDIP, WKEX, WDMO, WITO,
     &              W, Z )

      ! printout the resulting eigenstates:
      Call pa_preigen( exch,  lmax,  ibas, Dipol, AnisoLines1,
     &                 AnisoLines3, AnisoLines9, KE, DM_exchange,
     &                 JITO_exchange,
     &                 WLIN,  WDIP, WKEX, WDMO, WITO, W, Z, iPrint )
      !Z =  exchange eigenstates:
      NmaxPop=500
      If(NmaxPop.gt.exch) Then
        NmaxPop=exch
      End If
      Call PopAnalysis( nneq, neq, exch, nexch, nmax, lmax, NmaxPop,Z)
      Do i=exch,1,-1
        w(i)=w(i)-w(1)
      End Do

      !some verification
      If(dnrm2_(exch,WLIN,1).gt.1.0d-13)
     &   Call Add_Info('EXCHCTL::  WLIN',WLIN(1:NmaxPop),NmaxPop,8)
      If(dnrm2_(exch,WDIP,1).gt.1.0d-13)
     &   Call Add_Info('EXCHCTL::  WDIP',WDIP(1:NmaxPop),NmaxPop,8)
      If(dnrm2_(exch,WKEX,1).gt.1.0d-13)
     &   Call Add_Info('EXCHCTL::  WKEX',WKEX(1:NmaxPop),NmaxPop,8)
      If(dnrm2_(exch,W,1).gt.1.0d-13)
     &   Call Add_Info('EXCHCTL::     W',W(1:exch),exch,8)
c compute the moments:
      Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,M,1)
      Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,S,1)
      If(dbg) Then
       Write(6,'(A)') 'Magnetic moments before the build of '//
     &                'coupled M and S matrices'
       If(nPair>0) Then
        Do lp=1,npair
          lb1=i_pair(lp,1)
          lb2=i_pair(lp,2)
          i1=nind(lb1,1) ! indices of non-equivalent sites
          i2=nind(lb2,1) ! indices of non-equivalent sites
          n1=nexch(i1)
          n2=nexch(i2)
          Call prMom('EXCHCTL,Before M ans S, SM(i1):',
     &                SM(i1,1:3,1:n1,1:n1), n1 )
          Call prMom('EXCHCTL,Before M ans S, SM(i2):',
     &                SM(i2,1:3,1:n2,1:n2), n2 )
          Call prMom('EXCHCTL,Before M ans S, MM(i1):',
     &                MM(i1,1:3,1:n1,1:n1), n1 )
          Call prMom('EXCHCTL,Before M ans S, MM(i2):',
     &                MM(i2,1:3,1:n2,1:n2), n2 )
        End Do
       End If
      End If

      Do L=1,3
        Do isite=1,lmax
          Do nb1=1,exch
            Do lp=1,lmax
              icoord(lp)=ibas(nb1,lp)
            End Do
              i1=nind(isite,1)
              j1=nind(isite,2)
            is1=ibas(nb1,isite)+1

            Do js1=1,nexch(i1)
              icoord(isite)=js1-1
              nb2=norder(icoord,intc,lmax)
              M( l, nb1, nb2 ) = M( l, nb1, nb2 ) +
     &              rot( i1, j1, l, 1 ) * MM( i1, 1, is1, js1 )
     &             +rot( i1, j1, l, 2 ) * MM( i1, 2, is1, js1 )
     &             +rot( i1, j1, l, 3 ) * MM( i1, 3, is1, js1 )
              S( l, nb1, nb2 ) = S( l, nb1, nb2 ) +
     &              rot( i1, j1, l, 1 ) * SM( i1, 1, is1, js1 )
     &             +rot( i1, j1, l, 2 ) * SM( i1, 2, is1, js1 )
     &             +rot( i1, j1, l, 3 ) * SM( i1, 3, is1, js1 )

            End Do  ! js1
          End Do  ! nb1
        End Do  ! isite

        ! magnetic moment
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,TMP,1)
        Call zgemm_('C','N',EXCH,EXCH,EXCH,
     &             (1.0_wp,0.0_wp),Z, EXCH,
     &                             M(L,:,:), EXCH,
     &             (0.0_wp,0.0_wp),TMP, EXCH )
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,M(L,:,:),1)
        Call zgemm_('N','N',EXCH,EXCH,EXCH,
     &             (1.0_wp,0.0_wp),TMP,EXCH,
     &                               Z,EXCH,
     &             (0.0_wp,0.0_wp), M(L,:,:), EXCH )
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,TMP,1)
        ! spin moment
        Call zgemm_('C','N',EXCH,EXCH,EXCH,
     &             (1.0_wp,0.0_wp),Z,EXCH,
     &                             S(L,:,:), EXCH,
     &             (0.0_wp,0.0_wp),TMP,EXCH )
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,S(L,:,:),1)
        Call zgemm_('N','N',EXCH,EXCH,EXCH,
     &             (1.0_wp,0.0_wp),TMP,EXCH,
     &                               Z,EXCH,
     &             (0.0_wp,0.0_wp),S(L,:,:),EXCH )
      End Do  ! L


      If(npair >0 ) Then
      ! ITO decomposition of the exchange and dipolar interactions:
      Call pr_ito_int( npair, i_pair, lmax, nexch, nneq, neqv, itype,
     &                 neq, nmax, eso(1:nneq,1:nmax),
     &                  MM(1:nneq,1:3,1:nmax,1:nmax),
     &                  SM(1:nneq,1:3,1:nmax,1:nmax), rot,
     &                 Dipol, AnisoLines1, AnisoLines3, AnisoLines9,
     &                 DM_exchange, JITO_exchange,
     &                 HLIN1, HLIN3, HLIN9, HDIP, HDMO, HITO )
      End If

      ! projection on the Ising Hamiltonian
      ! accessible format
c      Write(6,'(/)')
c      Write(6,'(A)')'Complete decomposition of the exchange interaction'
c      Write(6,'(A)')




!-----------------------------------------------------------------------
! deallocate memory for this function:
      If(lmax>0) Then
        ! exchange energy spectrum
        Call mma_deallocate(intc)
        Call mma_deallocate(icoord)
        Call mma_deallocate(nind)
        If(exch>0) Then
          Call mma_deallocate(ibas)
        End If
      End If
      If(exch>0) Then
        Call mma_deallocate(wlin)
        Call mma_deallocate(wlin1)
        Call mma_deallocate(wlin3)
        Call mma_deallocate(wlin9)
        Call mma_deallocate(wdip)
        Call mma_deallocate(wkex)
        Call mma_deallocate(wdmo)
        Call mma_deallocate(wito)
      End If
      If(nmax>0) Then
        Call mma_deallocate(S1)
        Call mma_deallocate(M1)
        Call mma_deallocate(S2)
        Call mma_deallocate(M2)
        Call mma_deallocate(ZA1)
        Call mma_deallocate(ZA2)
        Call mma_deallocate(SM1)
        Call mma_deallocate(SM2)
        Call mma_deallocate(MM1)
        Call mma_deallocate(MM2)
        If(npair>0) Then
          Call mma_deallocate(HLIN1)
          Call mma_deallocate(HLIN3)
          Call mma_deallocate(HLIN9)
          Call mma_deallocate(HDIP)
          Call mma_deallocate(HKEX)
          Call mma_deallocate(HDMO)
          Call mma_deallocate(HITO)
        End If
      End If

      If(exch>0) Then
        Call mma_deallocate(tmp)
      End If

      If(nneq>0) Then
        Call mma_deallocate(nexchR)
        Call mma_deallocate(SMR)
        Call mma_deallocate(MMR)
        If(neqv>0) Then
          Call mma_deallocate(rotR)
        End If
      End If

      If(exchR>0) Then
        If(lmax>0) Then
          Call mma_deallocate(ibasR)
        End If
        Call mma_deallocate(WR)
        Call mma_deallocate(ZR)
        Call mma_deallocate(MR)
        Call mma_deallocate(SR)
      End If

      If(npair>0) Then
        Call mma_deallocate(HKEXR)
      End If

      If(lmax>0) Then
        Call mma_deallocate(intcR)
      End If




c  results of projection of the exchange interaction on the Ising Hamiltonian:




c      open (76,file='eigenstates.txt')
c      Do i=1,EXCH
c        Do j=1,EXCH
c      Write(76,'(2E24.14)') Z(i,j)
c        End Do
c      End Do
c      close(76)
c
c      open (77,file='eigenstates_full.txt')
c      ZZR=0.0_wp
c      ZZI=0.0_wp
c      Do i=1,EXCH
c        Do j=1,EXCH
c      read(77,'(2E24.14)') ZZR(i,j),ZZI(i,j)
c        End Do
c      End Do
c      close(77)
c      ZF=(0.0_wp,0.0_wp)
c      Do i=1,EXCH
c        Do j=1,EXCH
c        ZF(i,j)=cmplx( ZZR(i,j), ZZI(i,j), 8)
c        End Do
c      End Do
c      OVLP=(0.0_wp,0.0_wp)
c      Call ZGEMM_('C','N',EXCH,  EXCH,  EXCH, (1.0_wp,0.0_wp),
c     &             Z(  1:EXCH,1:EXCH), EXCH,
c     &            ZF(  1:EXCH,1:EXCH), EXCH, (0.0_wp,0.0_wp),
c     &          OVLP(  1:EXCH,1:EXCH), EXCH )
c      Do i=1,1
c        Do j=1,EXCH
c        If(ABS(OVLP(j,i)).gt.0.0000001_wp ) Then
c      Write(6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)')
c     & 'OVLP-1(',j,',',i,') = ',OVLP(j,i),
c     & 'ABS = ', ABS(OVLP(j,i)), 'ABS^2 = ',
c     & ABS(OVLP(j,i))*ABS(OVLP(j,i))
c       End If
c        End Do
c      End Do
c      OVLP=(0.0_wp,0.0_wp)
c      Call ZGEMM_('C','N',EXCH,  EXCH,  EXCH, (1.0_wp,0.0_wp),
c     &            ZF(  1:EXCH,1:EXCH), EXCH,
c     &             Z(  1:EXCH,1:EXCH), EXCH, (0.0_wp,0.0_wp),
c     &          OVLP(  1:EXCH,1:EXCH), EXCH )
c      Do i=1,1
c        Do j=1,EXCH
c        If(ABS(OVLP(j,i)).gt.0.0000001_wp ) Then
c      Write(6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)')
c     & 'OVLP-2(',i,',',j,') = ',OVLP(i,j),
c     & 'ABS = ', ABS(OVLP(i,j)), 'ABS^2 = ',
c     & ABS(OVLP(i,j))*ABS(OVLP(i,j))
c        End If
c        End Do
c      End Do
c 199  Continue
      Call qExit('PA_exchctl')
      Return
      End



      Subroutine prep_mom_exchange( n, R, S, M, mg, dbg)

      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)             :: n
      Real(kind=wp), intent(in)       :: R(3,3)
      Real(kind=wp), intent(out)      :: mg(3,3)
      Complex(kind=wp), intent(inout) :: S(3,n,n), M(3,n,n)
      Logical                         :: dbg
      ! local data:
      Real(kind=wp)                   :: g(3)
      Complex(kind=wp), allocatable   :: Mt(:,:,:), St(:,:,:), Z(:,:)

      Call qEnter('PA_prep_mom_exch')
!-----------------------------------------------------------------------
      Call mma_allocate(Mt,3,n,n,'Mt')
      Call mma_allocate(St,3,n,n,'St')
      Call mma_allocate(Z,n,n,'Z')

      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,Mt,1)
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,St,1)
      Call dcopy_(3  ,[0.0_wp],0,  g,1)
      Call dcopy_(3*3,[0.0_wp],0, mg,1)
      ! make a local backup of the data:
      Call zcopy_(3*n*n,M,1,Mt,1)
      Call zcopy_(3*n*n,S,1,St,1)

      If(dbg) Call prMom('PA_prep_mom_exch, input S',St,n)
      If(dbg) Call prMom('PA_prep_mom_exch, input M',Mt,n)


      ! rotate the momentum using the R rotation matrix --
      ! to the local axes for a symmetric compound:
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
      Call rotmom2( St, n, R, S)
      Call rotmom2( Mt, n, R, M)
      ! back-up again:
      Call zcopy_(3*n*n,M,1,Mt,1)
      Call zcopy_(3*n*n,S,1,St,1)



      ! find local magnetic axes:
      Call atens( M, n, g, mg, 2)
      ! rotate the momentum using the  mg  rotation matrix --
      ! to the local magnetic axes:
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
      Call rotmom2( St, n, mg, S)
      Call rotmom2( Mt, n, mg, M)



      ! find local pseudospin:
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Z,1)
      Call pseudospin( M, n, Z, 3,1, 1)
      If(dbg) Call pa_prmat('PA_prep_mom_exch, Z:',Z,n)

      ! Transform the moment into their local pseudospins
      Call UTMU2( n, n, Z, S )
      Call UTMU2( n, n, Z, M )
      If(dbg) Call prMom('PA_prep_mom_exch, S:', S, n)
      If(dbg) Call prMom('PA_prep_mom_exch, M:', M, n)
      ! back-up again:
      Call zcopy_(3*n*n,M,1,Mt,1)
      Call zcopy_(3*n*n,S,1,St,1)



      ! rotate back the moment, so that we preserve the
      ! original coordinate system of the computed molecule
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,M,1)
      Call zcopy_(3*n*n,[(0.0_wp,0.0_wp)],0,S,1)
      Call rotmom( St, n, mg, S)
      Call rotmom( Mt, n, mg, M)

!-----------------------------------------------------------------------

      Call mma_deallocate(Mt)
      Call mma_deallocate(St)
      Call mma_deallocate(Z)

      Call qExit('PA_prep_mom_exch')

!-----------------------------------------------------------------------
! old preparation of the data for Lines exchange
!        ! rotate the moments to the general coordinate system
!        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S1,1)
!        Call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0, S2,1)
!
!        Call rotmom2( SM(i1,1:3,1:n1,1:n1), n1, rot(i1,j1,1:3,1:3),
!     &                S1(1:3,1:n1,1:n1) )
!        Call rotmom2( SM(i2,1:3,1:n2,1:n2), n2, rot(i2,j2,1:3,1:3),
!     &                S2(1:3,1:n2,1:n2) )
!        Call Lines_Exchange( Jex(lp), n1, n2,
!     &                       S1(1:3,1:n1,1:n1),
!     &                       S2(1:3,1:n2,1:n2),
!     &              HLIN1(lp,1:n1,1:n1,1:n2,1:n2) )
      Return
      End
