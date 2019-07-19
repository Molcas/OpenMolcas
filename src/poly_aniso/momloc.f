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
      Subroutine momloc2( N, NL, nneq,  neq, neqv, r_rot, nsites,
     &                    nexch, W, Z, dipexch, s_exch, dipso, s_so )
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer          :: N
      Integer          :: NL
      Integer          :: nneq
      Integer          :: neq(nneq)
      Integer          :: neqv  ! neqv = MAXVAL(neq(:))
      Integer          :: nsites
      Integer          :: nexch(nneq)
      Real(kind=wp)    :: W(N)
!     assuming 10 equivalent magnetic sites, which is too much for many cases
      Real(kind=wp)    :: R_rot(NNEQ,neqv,3,3)
      Complex(kind=wp) :: dipexch(3,N,N)
      Complex(kind=wp) ::  s_exch(3,N,N)
      Complex(kind=wp) :: dipso(nneq,3,NL,NL)
      Complex(kind=wp) ::  s_so(nneq,3,NL,NL)
      Complex(kind=wp) :: Z(N,N)
#include "stdalloc.fh"
c local variables:
      Integer          :: L, i, j, m, k
      Integer          :: nmult, isite
      Integer          :: icod(nsites)
      Integer          :: ib(N,nsites)
      Integer          :: nind(nsites,2)
      Integer          :: l_exch, jEnd
      Integer          :: i1,j1,iss1,jss1,nb1,nb2,iss
      Integer          :: icoord(nsites)
      Integer          :: norder
      Real(kind=wp)    :: gtens(3)
      Real(kind=wp)    :: maxes(3,3)
      Real(kind=wp)    :: st(3)
      Real(kind=wp)    :: H
      Real(kind=wp)    :: E_thres
      Real(kind=wp)    :: zJ
      Real(kind=wp)    :: g_e
      Character(60)    :: fmtline
      logical          :: DBG
      Real(kind=wp),allocatable     :: WM(:)     ! WM(N)
      Real(kind=wp), allocatable    :: MM(:,:,:) ! MM(nsites,3,N)
      Real(kind=wp), allocatable    :: LM(:,:,:) ! LM(nsites,3,N)
      Real(kind=wp), allocatable    :: SM(:,:,:) ! SM(nsites,3,N)
      Real(kind=wp), allocatable    :: JM(:,:,:) ! JM(nsites,3,N)
      Complex(kind=wp), allocatable :: ZM(:,:)  ! ZM(N,N)
      Complex(kind=wp), allocatable :: VL(:,:)  ! VL(N,N)
      Complex(kind=wp), allocatable :: TMP(:,:) ! TMP(N,N)
      ! temporary data for ZEEM:
      Real(kind=wp), allocatable :: RWORK(:)
      Complex(kind=wp), allocatable :: HZEE(:), WORK(:), W_c(:)
c
      Call qEnter('PA_momloc2')
      DBG=.false.
      g_e=2.0023193043718_wp

      Call mma_allocate(WM,N,'WM')
      Call mma_allocate(MM,nsites,3,N,'MM')
      Call mma_allocate(LM,nsites,3,N,'LM')
      Call mma_allocate(SM,nsites,3,N,'SM')
      Call mma_allocate(JM,nsites,3,N,'JM')
      Call mma_allocate(ZM,N,N,'ZM')
      Call mma_allocate(VL,N,N,'VL')
      Call mma_allocate(TMP,N,N,'TMP')
      ! temporary arrays used in ZEEM_SA:
      Call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
      Call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
      Call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
      Call mma_allocate(W_c,N,'ZEEM_W_c')

      ! zero everything:
      Call dcopy_(N,[0.0_wp],0,WM,1)
      Call dcopy_(nsites*3*N,[0.0_wp],0,MM,1)
      Call dcopy_(nsites*3*N,[0.0_wp],0,LM,1)
      Call dcopy_(nsites*3*N,[0.0_wp],0,SM,1)
      Call dcopy_(nsites*3*N,[0.0_wp],0,JM,1)
      Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,ZM,1)
      Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
      Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,TMP,1)

      Call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
      Call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
      Call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
      Call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)

c  initialisations:
      isite=0
      nind=0
      Do i=1,nneq
         Do j=1,neq(i)
            isite=isite+1
            nind(isite,1)=i
            nind(isite,2)=j
         End Do
      End Do
      isite=0
      icod=0
      icod(1)=1
      Do i=2,nsites
         isite=nind(i-1,1)
         icod(i) = icod(i-1)*nexch(isite)
      End Do
      ib=0
      Do i=1,N
         j=i-1
         Do isite=1,nsites
            k=nsites-isite+1
            ib(i,k) = j/icod(k)
            j = j - ib(i,k) * icod(k)
         End Do
      End Do
c  find the multiplicity of the low-lying group of states,
c  energy threshold E_thres = 1.d-2 cm-1;
      nmult=0
      E_thres=1.d-3
      Do i=1,N
        If( ABS(W(i)-W(1)).lt.E_thres) Then
           nmult=nmult+1
        End If
      End Do
      If(nmult<2) nmult=2 !minimum value needed to compute g-tensor
c  find the main magnetic axes of this manIfold:
      gtens=0.0_wp
      maxes=0.0_wp
      Call atens( dipexch(1:3,1:nmult,1:nmult), nmult, gtens, maxes, 1)
      If(DBG) Then
      Write(6,'(A)') 'MOMLOC2:  g tensor of the ground manIfold:'
      Write(6,*)
      Write(6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(1),
     & ' axis X: ',(maxes(j,1),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(2),
     & ' axis Y: ',(maxes(j,2),j=1,3)
      Write(6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(3),
     & ' axis Z: ',(maxes(j,3),j=1,3)
      End If
c  construct the Zeeman matrix in the lowest N exchange states  and diagonalize it
      st=0.0_wp
      H=1.d-4 ! Tesla
      zJ=0.0_wp !absence of intermolecular interaction
      Call zeem_sa( N, H, maxes(1,3),maxes(2,3),maxes(3,3),
     &           W(1:N), dipexch(1:3,1:N,1:N),
     &                    s_exch(1:3,1:N,1:N),
     &           ST, zJ,  WM(1:N), ZM(1:N,1:N),
     &           DBG, RWORK, HZEE, WORK, W_c )
ccc  eigenvectors
      If (DBG) Then
         Write(6,*)
         Write(6,'(100a)') (('%'),i=1,96)
         Write(6,'(10x,a)') 'MOMLOC2:  '//
     &       'EigenVectors of the Total Magnetic Interaction'
         Write(6,'(100a)') (('%'),i=1,96)

         Write(6,'(A,16A)') '-',('---',m=1,nsites),'|',
     &       ('---------------------------|',i=1,4)
         Write(fmtline,'(A,i2,A)') '(1x,',nsites,'A,39x,A,38x,A)'
         Write(6,fmtline) ('   ',m=1,nsites),
     &        'eigenvectors of the exchange matrix','|'
         Do J=1,N,4
            jEnd=MIN(N,J+3)
            Write(6,'(A,16A)') '-',('---',m=1,nsites),'|',
     &       ('---------------------------|',i=j,jEnd)

            Write(6,'(A,5A)') 'Exch.  |',
     &       ('      exchange state       |',i=j,jEnd)

            Write(6,'(A,6A)') 'basis  |',
     &       ('                           |',i=j,jEnd)
            Write(6,'(A,6(a,i5,a))') 'on site|',
     &       ('         ',i,'             |',i=j,jEnd)
              Write(fmtline,'(A,i2,A)')
     &              '(A,',nsites,'i3,a,6(a,f19.12,2x,a))'
            Write(6,fmtline) ' ',(i,i=1,nsites),'|',
     &       ('  E =',wm(i)-wm(1),' |',i=j,jEnd)

            Write(6,'(A,16A)') '-',('---',m=1,nsites),'|',
     &       ('----- Real ------- Imag ---|',i=j,jEnd)

               Write(fmtline,'(A,i2,A)')
     &              '(A,',nsites,'i3,A,5(2F13.9,1x,a))'
            Do iss=1,N
               Write(6,fmtline) '<',
     &              ( ib(iss,m)+1,m=1,nsites),'|',
     &              ( ZM(iss,i),'|',i=j,jEnd)
            End Do  !iss
            Write(6,'(A,16A)') '-',('---',m=1,nsites),'|',
     &          ('---------------------------|',i=j,jEnd)
            Write(6,*)
         End Do ! j
      End If !DBG


c  generate the exchange basis:
      Do isite=1,nsites
         Do L=1,3
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Do nb1=1,N
                Do l_exch=1,nsites
                icoord(l_exch)=ib(nb1,l_exch)
                End Do
            i1=nind(isite,1)
            j1=nind(isite,2)
          iss1=ib(nb1,isite)+1

               Do jss1=1,nexch(i1)
               icoord(isite)=jss1-1
               nb2=norder(icoord,icod,nsites)
               VL(nb1, nb2) = VL(nb1, nb2) +
     &               r_rot( i1, j1, l, 1 ) * dipso( i1, 1, iss1, jss1 )
     &              +r_rot( i1, j1, l, 2 ) * dipso( i1, 2, iss1, jss1 )
     &              +r_rot( i1, j1, l, 3 ) * dipso( i1, 3, iss1, jss1 )
               End Do  ! jss1
            End Do  ! nb1
c rotate this matrix to exchange basis:
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,TMP,1)
            Call ZGEMM_('C','N',N,N,N,
     &                   (1.0_wp,0.0_wp),   Z, N,
     &                                     VL, N,
     &                   (0.0_wp,0.0_wp), TMP, N )
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Call ZGEMM_('N','N',N,N,N,
     &                   (1.0_wp,0.0_wp),TMP, N,
     &                                     Z, N,
     &                   (0.0_wp,0.0_wp), VL, N )
c rotate this matrix to Zeeman basis:
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,TMP,1)
            Call ZGEMM_('C','N',N,N,N,(1.0_wp,0.0_wp),
     &                 ZM(1:N,1:N), N,
     &                 VL(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N )
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Call ZGEMM_('N','N',N,N,N, (1.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N,
     &                 ZM(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                 VL(1:N,1:N), N )
            Do i=1,N
               MM(isite,L,i)=dble(VL(i,i))
            End Do
c spin moment
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Do nb1=1,N
                Do l_exch=1,nsites
                icoord(l_exch)=ib(nb1,l_exch)
                End Do
            i1=nind(isite,1)
            j1=nind(isite,2)
          iss1=ib(nb1,isite)+1

               Do jss1=1,nexch(i1)
               icoord(isite)=jss1-1
               nb2=norder(icoord,icod,nsites)
               VL(nb1, nb2) = VL(nb1, nb2) +
     &               r_rot( i1, j1, l, 1 ) *  s_so( i1, 1, iss1, jss1 )
     &              +r_rot( i1, j1, l, 2 ) *  s_so( i1, 2, iss1, jss1 )
     &              +r_rot( i1, j1, l, 3 ) *  s_so( i1, 3, iss1, jss1 )
               End Do  ! jss1
            End Do  ! nb1
c rotate this matrix to exchange basis and to Zeeman basis
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,TMP,1)
            Call ZGEMM_('C','N',N,N,N,  (1.0_wp,0.0_wp),
     &                  Z(1:N,1:N), N,
     &                 VL(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N )
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Call ZGEMM_('N','N',N,N,N,  (1.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N,
     &                  Z(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                 VL(1:N,1:N), N )
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,TMP,1)
            Call ZGEMM_('C','N',N,N,N,  (1.0_wp,0.0_wp),
     &                 ZM(1:N,1:N), N,
     &                 VL(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N )
            Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,VL,1)
            Call ZGEMM_('N','N',N,N,N,  (1.0_wp,0.0_wp),
     &                TMP(1:N,1:N), N,
     &                 ZM(1:N,1:N), N, (0.0_wp,0.0_wp),
     &                 VL(1:N,1:N), N )
            Do i=1,N
               SM(isite,L,i)=dble(VL(i,i))
            End Do
         End Do  ! L
      End Do  ! isite
c compute and print the calculated expectation values:
      Write(6,*)
      Write(6,'(A)') 'EXPECTATION VALUES'
      Write(6,'(5A)') '--------|----|',
     & ('------------------------------|',i=1,4)
      Write(6,'(5A)') 'Exchange|Site|',
     & '   MAGNETIC MOMENT (M=-L-2S)  |',
     & '        SPIN MOMENT (S)       |',
     & '      ORBITAL MOMENT (L)      |',
     & '     TOTAL MOMENT (J=L+S)     |'
      Write(6,'(5A)') ' state  | Nr.|',
     & ('     X         Y         Z    |',i=1,4)
      Write(6,'(5A)') '--------|----|',
     & ('------------------------------|',i=1,4)
      Do i=1,N
c  we proceed to compute expectation values for this nb1 exchange state
         Do isite=1,nsites
            Do L=1,3
               LM(isite,L,i)=-MM(isite,L,i)-g_e*SM(isite,L,i)
               JM(isite,L,i)= LM(isite,L,i)+    SM(isite,L,i)
            End Do
         End Do

         Do isite=1,nsites
         If(isite.eq.int((nsites+1)/2)) Then
      Write(6,'(i5,3x,A,1x,i2,1x,A,4(3(F9.5,1x),A))')
     & i,'|',isite,'|',
     & MM(isite,1,i),MM(isite,2,i),MM(isite,3,i),'|',
     & SM(isite,1,i),SM(isite,2,i),SM(isite,3,i),'|',
     & LM(isite,1,i),LM(isite,2,i),LM(isite,3,i),'|',
     & JM(isite,1,i),JM(isite,2,i),JM(isite,3,i),'|'
         Else
      Write(6,'(8x,A,1x,i2,1x,A,4(3(F9.5,1x),A))')
     & '|',isite,'|',
     & MM(isite,1,i),MM(isite,2,i),MM(isite,3,i),'|',
     & SM(isite,1,i),SM(isite,2,i),SM(isite,3,i),'|',
     & LM(isite,1,i),LM(isite,2,i),LM(isite,3,i),'|',
     & JM(isite,1,i),JM(isite,2,i),JM(isite,3,i),'|'
         End If
         End Do
      Write(6,'(5A)') '--------|----|',
     & ('------------------------------|',i1=1,4)
      End Do

      ! deallocate temporary arrays:
      Call mma_deallocate(WM)
      Call mma_deallocate(MM)
      Call mma_deallocate(LM)
      Call mma_deallocate(SM)
      Call mma_deallocate(JM)
      Call mma_deallocate(ZM)
      Call mma_deallocate(VL)
      Call mma_deallocate(TMP)
      Call mma_deallocate(RWORK)
      Call mma_deallocate(HZEE)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W_c)


      Call qExit('PA_momloc2')
      Return
      End Subroutine momloc2
