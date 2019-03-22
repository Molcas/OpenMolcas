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
      Subroutine pr_ito_int( npair, i_pair, lmax, nexch, nneq, neqv,
     &                       itype, neq, nmax, soe, MM, SM, rot, Dipol,
     &                       AnisoLines1, AnisoLines3, AnisoLines9,
     &                       DM_exchange, JITO_exchange,
     &                       HLIN1, HLIN3, HLIN9, HDIP, HDMO, HITO )
!     this function prints the parameters of the exchange interaction in an
!     accessible format
      Implicit None
#include "stdalloc.fh"
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: npair
      Integer, intent(in)          :: nneq
      Integer, intent(in)          :: neqv
      Integer, intent(in)          :: lmax
      Integer, intent(in)          :: nmax
      Integer, intent(in)          :: nexch(nneq)
      Integer, intent(in)          :: neq(nneq)
      Integer, intent(in)          :: i_pair(npair,2)
      Real(kind=wp), intent(in)    :: rot(nneq,neqv,3,3)
      Real(kind=wp), intent(in)    :: soe(nneq,nmax)
      Complex(kind=wp), intent(in) :: MM(nneq,3,nmax,nmax)
      Complex(kind=wp), intent(in) :: SM(nneq,3,nmax,nmax)
      Complex(kind=wp), intent(in) :: HLIN1(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), intent(in) :: HLIN3(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), intent(in) :: HLIN9(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), intent(in) :: HDIP(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), intent(in) :: HDMO(npair,nmax,nmax,nmax,nmax)
      Complex(kind=wp), intent(in) :: HITO(npair,nmax,nmax,nmax,nmax)
      Logical, intent(in)          :: Dipol
      Logical, intent(in)          :: AnisoLines1
      Logical, intent(in)          :: AnisoLines3
      Logical, intent(in)          :: AnisoLines9
      Logical, intent(in)          :: DM_exchange
      Logical, intent(in)          :: JITO_exchange
      Character(1), intent(in)     :: itype(nneq)
      ! local variables
      Integer       ::   i,j,l,k,lp,i1,i2,j1,j2,lb1,lb2,iopt,ibuf,
     &                   is1,is2,js1,js2,k1,k2,q1,q2,n1,n2,nsize
      Integer       ::   nind(lmax,2),l1(2),l2(2),l3(2),l4(2)
      Real(kind=wp) ::   J1C(3,3), J1Cr(3,3) !, J1C_trans(3,3)
      Complex(kind=wp), allocatable :: JN(:,:,:,:)
      Complex(kind=wp), allocatable :: JB(:,:,:,:)
      Complex(kind=wp), allocatable :: JS(:,:,:,:)
      Real(kind=wp)    :: dznrm2_,RL1,RL3,RL9,RDI,RDM,RIT
      Real(kind=wp)    :: g1(3),g2(3),mg1(3,3),mg2(3,3)
      External         :: dznrm2_
      Real(kind=wp)    :: cm_to_MHz
      logical DBG

      cm_to_MHz=29979.2458_wp
      DBG=.false.
c some initializations:
      nind(:,:)=0
      l=0
      Do i=1,nneq
        Do j=1,neq(i)
          l=l+1
          nind(l,1)=i
          nind(l,2)=j
        End Do
      End Do

      ibuf=npair*nmax*nmax*nmax*nmax
      If(ibuf==0) Then
         Write(6,'(A)') 'in UTMU:   ibuf=0 !!!'
         Write(6,*) 'npair= ', npair
         Write(6,*) 'nmax = ', nmax
         Call xFlush(6)
         Call xquit(128)
      End If
      RL1=dznrm2_(ibuf,HLIN1,1)
      RL3=dznrm2_(ibuf,HLIN3,1)
      RL9=dznrm2_(ibuf,HLIN9,1)
      RDI=dznrm2_(ibuf,HDIP,1)
      RDM=dznrm2_(ibuf,HDMO,1)
      RIT=dznrm2_(ibuf,HITO,1)
cccccccccccccccccccccccccccccccc
      If(DBG) Then
        Write(6,'(A,i6)') ' nmax =', nmax
        Write(6,'(A,i6)') ' lmax =', lmax
        Write(6,'(A,i6)') 'npair =', nmax
        Write(6,'(A,i6)') ' nneq =', nneq
        Write(6,'(20A)') ' itype =', (itype(i),i=1,nneq)
        Write(6,'(A)') 'i_pair'
        Do i=1,npair
          Write(6,'(2I5)') i_pair(i,1),i_pair(i,2)
        End Do
        Write(6,'(A)') 'equivalent sites'
        Do i=1, nneq
          Write(6,'(A,i2,A,i4)') 'neq(',i,') =',neq(i)
        End Do
        Write(6,'(A)') 'exchange basis'
        Do i=1, nneq
          Write(6,'(A,i2,A,i4)') 'nexch(',i,') =',nexch(i)
        End Do
        Write(6,'(A)') 'rotation matrix'
        Do i=1,nneq
          Write(6,'(A,i3)') 'site', i
          Do j=1,neq(i)
            Do l=1,3
              Write(6,'(3F12.6)') (rot(i,j,l,i1),i1=1,3)
            End Do
          End Do
        End Do
        Write(6,'(A,i3)') 'site', i
        Write(6,'(A)' ) 'magnetic moment, initial'
        Do i=1,nneq
          Write(6,'(A ,i3)') 'site', i
          Call prMom('pr_ito_int:: magnetic moment, initial',
     &                mm(i,:,:,:),nmax)
        End Do
        Write(6,'(A)' ) 'spin-orbit energies'
        Do i=1,nmax
          Write(6,'(1 5F14.5)') (soe(j,i),j=1,nneq)
        End Do

        Do lp=1,npair
         lb1=i_pair(lp,1)
         lb2=i_pair(lp,2)
          i1=nind(lb1,1) ! indices of non-equivalent sites
          i2=nind(lb2,1) ! indices of non-equivalent sites
          j1=nind(lb1,2) ! indices of equivalent sites
          j2=nind(lb2,2) ! indices of equivalent sites
          If (AnisoLines1.AND.(RL1>0.0_wp)) Then
            Write(6,'(A,i5)') 'HLIN1,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HLIN1(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

          If (AnisoLines3.AND.(RL3>0.0_wp)) Then
            Write(6,'(A,i5)') 'HLIN3,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HLIN3(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

          If (AnisoLines9.AND.(RL9>0.0_wp)) Then
            Write(6,'(A,i5)') 'HLIN9,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HLIN9(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

          If (Dipol.AND.(RDI>0.0_wp)) Then
            Write(6,'(A,i5)') 'HDIP,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HDIP(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

          If (DM_exchange.AND.(RDM>0.0_wp)) Then
            Write(6,'(A,i5)') 'HDMO,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HDMO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

          If (JITO_exchange.AND.(RIT>0.0_wp)) Then
            Write(6,'(A,i5)') 'HITO,  interacting pair ',lp
            Do is1=1,nexch(i1)
              Do is2=1,nexch(i1)
                Do js1=1,nexch(i2)
                  Write(6,'(10(2F10.6,2x))')
     &                     (HITO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
                End Do
              End Do
            End Do
          End If

        End Do !lp
        Call prMom('SM(i1) bf Lines1',SM(i1,1:3,1:n1,1:n1),n1)
        Call prMom('SM(i2) bf Lines1',SM(i2,1:3,1:n2,1:n2),n2)
      End If !DBG


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Write(6,*)
      Write(6,'(100a)') (('%'),j=1,100)
      If ( (.not.AnisoLines1).and.(.not.AnisoLines3).and.
     &     (.not.AnisoLines9).and.(.not.Dipol).and.
     &     (.not.DM_exchange).and.(.not.JITO_exchange) ) Then
        Write(6,'(20x,A)') 'ITO decomposition of exchange and/or '//
     &                     'dipolar couplings.'
        Write(6,'(A)') 'AnisoLines1 =  FALSE.'
        Write(6,'(A)') 'AnisoLines3 =  FALSE.'
        Write(6,'(A)') 'AnisoLines9 =  FALSE.'
        Write(6,'(A)') 'Dipol       =  FALSE.'
        Write(6,'(A)') 'JITO_exch   =  FALSE.'
        Write(6,'(A)') 'DM_exchange =  FALSE.'
        Write(6,'(A)') 'Nothing to Do.'
        Return

      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &         (.not.Dipol).and.(.not.DM_exchange).and.
     &         (.not.JITO_exchange) ) Then
        Write(6,'(20x,A)') 'ITO decomposition of the Lines exchange '//
     &                     'interactions.'

      Else If((.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9))
     &        .and.Dipol.and.(.not.DM_exchange).and.
     &        (.not.JITO_exchange)) Then
        Write(6,'(20x,A)') 'ITO decomposition of the dipolar '//
     &                     'interactions.'

      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and.(.not.DM_exchange).and.
     &          (.not.JITO_exchange) ) Then
        Write(6,'(20x,A)') 'ITO decomposition of exchange and/or '//
     &                     'dipolar couplings.'
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.(.not.Dipol).and.(.not.DM_exchange).and.
     &          JITO_exchange ) Then
        Write(6,'(20x,A)') 'ITO decomposition of anisotropic '//
     &                     'exchange interaction. '
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &         .and.Dipol.and.(.not.DM_exchange).and.
     &          JITO_exchange ) Then
        Write(6,'(20x,A)') 'ITO decomposition of anisotropic '//
     &                     'exchange interaction and dipolar '//
     &                     'couplings.'

      End If
      Write(6,'(100a)') (('%'),j=1,100)
      Write(6,*)

c decompose the exchange interaction in products of ITO
c first rotate the magnetic moments to the general coordinate system:
      Do lp=1,npair
        lb1=i_pair(lp,1)
        lb2=i_pair(lp,2)
        i1=nind(lb1,1) ! indices of non-equivalent sites
        i2=nind(lb2,1) ! indices of non-equivalent sites
        j1=nind(lb1,2) ! indices of equivalent sites
        j2=nind(lb2,2) ! indices of equivalent sites

        n1=nexch(i1)
        n2=nexch(i2)
        Write(6,'(A)') 'PART 1: Magnetic exchange is written in the '//
     &                 'coordinate systems of the LOCAL main '//
     &                 'magnetic axes of the interacting sites.'
        iopt=1
        Write(6,'(A)')
        Write(6,'(100A)') ('-',i=1,100)
        Write(6,'(A,i2)') 'Interacting pair',lp

        ! JN= exch. parameters in Naoya's ITO operators
        ! JL= exch. parameters in Liviu's ITO operators
        ! JS= exch. parameters in Stevens ESO operators
        l1(1)=  1
        l1(2)=  n1-1
        l2(1)=-(n1-1)
        l2(2)= (n1-1)
        l3(1)=  1
        l3(2)=  n2-1
        l4(1)=-(n2-1)
        l4(2)= (n2-1)
        nsize=(n1-1)*(n2-1)*(2*(n1-1)+1)*(2*(n2-1)+1)
        Call mma_allocate(JN,l1,l2,l3,l4,'JN')
        Call mma_allocate(JB,l1,l2,l3,l4,'JB')
        Call mma_allocate(JS,l1,l2,l3,l4,'JS')

!====================================================================
        If( AnisoLines1.AND.(RL1>0.0_wp) ) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HLIN1(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )


          ! re-write the first rank tensor in cartesian representation:
          ! using Naoya's ITO parameters
          J1C=0.0_wp
          J1Cr=0.0_wp
          Call tensor2cart( JB(1,-1:1,1,-1:1), J1C )

!         rotate cartesian J1C matrix by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          mg1=0.0_wp; mg2=0.0_wp; g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )
          Do i=1,3
           Do j=1,3
             Do l=1,3
               Do k=1,3
                 J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
               End Do
             End Do
           End Do
          End Do

        Write(6,'(A)')
        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: LINES-1'
        Write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') 'To be used directly in exchange Hamiltonian'//
     &                 ' of the kind:'
        Write(6,'(A)') ' H = -J * S1 * S2'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES22.14)') ( -J1Cr(i,j),j=1,3)
        End Do

        ! print out the data:
        Write(6,'(A)')
        Write(6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
        Write(6,'( 5x,A)') 'with absolute values larger than:  0.5d-14 '
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|           Line'//
     &                 's  Exchange  Interaction        |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'

        Do k1=1,n1-1,2
         Do q1=-k1,k1
          Do k2=1,n2-1,2
           Do q2=-k2,k2
            Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &            k1,'|',q1,'|',k2,'|',q2,'|', JN(k1,q1,k2,q2),'|'
           End Do
          End Do
         End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'


!      !------------------------------------------------------------------------
!      ! verify the back transform for HAM!:
!        Call mma_allocate(HAM,n1,n1,n2,n2,'HAM')
!        Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,HAM,1)
!
!        S1a=( 0.0_wp, 0.0_wp)
!        S2a=( 0.0_wp, 0.0_wp)
!        S1b=( 0.0_wp, 0.0_wp)
!        S2b=( 0.0_wp, 0.0_wp)
!
!        Call ESO(n1,1,1,S1b(1,1:n1,1:n1),S1b(2,1:n1,1:n1),redME)
!        Call ESO(n1,1,0,S1b(3,1:n1,1:n1),W1b(  1:n1,1:n1),redME)
!        Call zcopy_(3*n1*n1,S1b,1,S2b,1)
!
!        Call prmom('SM(i1) bf recover HAM',SM(i1,1:3,1:n1,1:n1),n1)
!        Call prmom('SM(i2) bf recover HAM',SM(i2,1:3,1:n2,1:n2),n2)
!
!
!        J1C_trans=0.0_wp
!        Do i=1,3
!         Do j=1,3
!          J1C_trans(i,j)=-J1Cr(i,j)
!         End Do
!        End Do
!        Write(6,'(/)')
!        Call Aniso_Lines_Exchange9( J1C_trans, n1, n2,
!     &          S1b(1:3,1:n1,1:n1),
!     &          S2b(1:3,1:n2,1:n2), HAM )
!        Call Aniso_Lines_Exchange9( J1C_trans, n1, n2,
!     &          SM(i1,1:3,1:n1,1:n1),
!     &          SM(i2,1:3,1:n2,1:n2), HAM )
!
!        Write(6,'(A,i5)') 'HLIN1: ORIG, REGEN, DIFF:'
!        Do is1=1,n1
!         Do is2=1,n1
!          Do js1=1,n2
!           Do js2=1,n2
!           Write(6,'(4I3,3(2ES22.13,3x))') is1,is2,js1,js2,
!     &               HLIN1(lp,is1,is2,js1,js2),
!     &                    HAM(is1,is2,js1,js2),
!     &      HLIN1(lp,is1,is2,js1,js2)-HAM(is1,is2,js1,js2)
!           End Do
!          End Do
!         End Do
!        End Do
!
!        Call mma_deallocate(HAM)
        End If



!====================================================================
        If( AnisoLines3.AND.(RL3>0.0_wp) ) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HLIN3(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )

          J1C=0.0_wp
          J1Cr=0.0_wp
          ! using Liviu's ITO parameters
          Call tensor2cart( JB( 1,-1:1, 1,-1:1), J1C )

!         rotate cartesian JLinC1 by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          mg1=0.0_wp;mg2=0.0_wp;g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )
          Do i=1,3
           Do j=1,3
             Do l=1,3
               Do k=1,3
                 J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
               End Do
             End Do
           End Do
          End Do

          Write(6,'(A)')
          Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                   '(rank-1) exchange interaction: LINES-3'
          Write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
          Write(6,'(A)') 'LOCAL AXES:::'
          Write(6,'(A)') 'To be used directly in exchange Hamiltonian'//
     &                   'of the kind:'
          Write(6,'(A)') ' H = -J * S1 * S2'
          Write(6,'(A)') '     (  xx   xy  xz  )  '
          Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
          Write(6,'(A)') '     (  zx   zy  zz  )  '
          Do i=1,3
            Write(6,'(3ES22.14)') ( -J1Cr(i,j),j=1,3)
          End Do

        ! print out the data:
        Write(6,'(A)')
        Write(6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
        Write(6,'( 5x,A)') 'with absolute values larger than:  0.1d-20 '
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|        Lines-3'//
     &                 '  Exchange  Interaction         |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'
        Do k1=1,N1-1,2
          Do q1=-k1,k1
            Do k2=1,N2-1,2
              Do q2=-k2,k2
!                If(ABS(JLin3(lp,k1,q1,k2,q2)) .gt. 0.1e-20_wp ) Then
                   Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &                   k1,'|',q1,'|',k2,'|',q2,'|',
     &                   JN(k1,q1,k2,q2),'|'
!                End If
              End Do
            End Do
          End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'
        End If

!====================================================================
        If( AnisoLines9.AND.(RL9>0.0_wp) ) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HLIN9(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )


          J1C=0.0_wp
          J1Cr=0.0_wp
          Call tensor2cart( JB(1,-1:1, 1,-1:1), J1C)

!         rotate cartesian JLinC1 by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          mg1=0.0_wp; mg2=0.0_wp; g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )
          Do i=1,3
           Do j=1,3
             Do l=1,3
               Do k=1,3
                 J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
               End Do
             End Do
           End Do
          End Do

        Write(6,'(A)')
        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: LINES-9'
        Write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') 'To be used directly in exchange Hamiltonian'//
     &                 'of the kind:'
        Write(6,'(A)') ' H = -J * S1 * S2'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
        End Do

        ! print out the data:
        Write(6,'(A)')
        Write(6,'(10x,A)') 'Parameters of the ITOs:'
        Write(6,'( 5x,A)') 'with absolute values larger than:  0.1e-20 '
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|        Lines-9'//
     &                 '  Exchange  Interaction         |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'
        Do k1=1,N1-1,2
         Do q1=-k1,k1
          Do k2=1,N2-1,2
           Do q2=-k2,k2
            Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &            k1,'|',q1,'|',k2,'|',q2,'|',
     &            JN(k1,q1,k2,q2),'|'
           End Do
          End Do
         End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'
        End If

!====================================================================


        If (Dipol.AND.(RDI>0.0_wp)) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HDIP(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )

          J1C=0.0_wp
          J1Cr=0.0_wp
          Call tensor2cart( JB(1,-1:1,1,-1:1), J1C(1:3,1:3) )

!         rotate cartesian JLinC1 by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          mg1=0.0_wp; mg2=0.0_wp; g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )
          Do i=1,3
           Do j=1,3
             Do l=1,3
               Do k=1,3
                 J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
               End Do
             End Do
           End Do
          End Do

        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: DIPOL'
        Write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') 'To be used directly in exchange Hamiltonian'//
     &                 'of the kind:'
        Write(6,'(A)') ' H = -J * S1 * S2'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
        End Do
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|      Dipolar  '//
     &                 'Exchange  Interaction           |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'
        Do k1=1,N1-1,2
          Do q1=-k1,k1
            Do k2=1,N2-1,2
              Do q2=-k2,k2
                   Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &                   k1,'|',q1,'|',k2,'|',q2,'|',
     &                   JN(k1,q1,k2,q2),'|'
              End Do
            End Do
          End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'
        End If

!====================================================================


        If (DM_exchange.AND.(RDM>0.0_wp)) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HDMO(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )

          J1C=0.0_wp
          J1Cr=0.0_wp
          Call tensor2cart( JB(1,-1:1,1,-1:1), J1C )

          mg1=0.0_wp; mg2=0.0_wp; g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )

!         rotate cartesian JLinC1 by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          Do i=1,3
           Do j=1,3
             Do l=1,3
               Do k=1,3
                 J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
               End Do
             End Do
           End Do
          End Do

        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: '//
     &                 'Dzyaloshinski-Morya'
        Write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') 'To be used directly in exchange Hamiltonian'//
     &                 'of the kind:'
        Write(6,'(A)') ' H = -J * S1 * S2'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
        End Do
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|     Dzyaloshin'//
     &                 'sky - Morya Interaction         |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'
        Do k1=1,N1-1,2
         Do q1=-k1,k1
          Do k2=1,N2-1,2
           Do q2=-k2,k2
            Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &            k1,'|',q1,'|',k2,'|',q2,'|',
     &            JN(k1,q1,k2,q2),'|'
           End Do
          End Do
         End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'
        End If

!====================================================================


        If (JITO_exchange.AND.(RIT>0.0_wp)) Then

          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)
          Call zcopy_(nsize,(0.0_wp,0.0_wp),0,
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),1)

          Call newjkqpar(n1,n2,HITO(lp,1:n1,1:n1,1:n2,1:n2),
     &            JN(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JB(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)),
     &            JS(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) )

          J1C(1:3,1:3)=0.0_wp
          Call tensor2cart( JB(1,-1:1,1,-1:1), J1C)

          mg1=0.0_wp; mg2=0.0_wp; g1=0.0_wp; g2=0.0_wp;
          Call atens( MM(i1,1:3,1:n1,1:n1), n1, g1, mg1, 2 )
          Call atens( MM(i2,1:3,1:n2,1:n2), n2, g2, mg2, 2 )

!         rotate cartesian JLinC1 by mg1 and mg2, in order to represent
!         the interaction matrix in the original coordinate system:
          J1Cr(1:3,1:3)=0.0_wp
          Do i=1,3
           Do j=1,3
            Do l=1,3
             Do k=1,3
              J1Cr(i,j)= J1Cr(i,j) + mg1(i,l)*mg2(j,k)*J1C(l,k)
             End Do
            End Do
           End Do
          End Do

        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: '//
     &                 'Anisotropic ITO exchange'
        Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
        End Do
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '--------------------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|     Anisotropi'//
     &                 'c ITO Exchange Interaction      |'
        Write(6,'(A)') '------|------|------|------|---------- Real'//
     &                 ' ----------------- Imag --------|'
        Do k1=1,N1-1,2
          Do q1=-k1,k1
            Do k2=1,N2-1,2
              Do q2=-k2,k2
                 Write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)')
     &                 k1,'|',q1,'|',k2,'|',q2,'|',
     &                 JN(k1,q1,k2,q2),'|'
              End Do
            End Do
          End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 ' -------------------------------|'
        End If

!====================================================================




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        Write(6,'(100A)') ('-',i=1,100)
c        Write(6,*)
c        Write(6,'(A)') 'PART 2: Magnetic exchange is written in the '//
c     &                 'GENERAL coordinate systems of the computed '//
c     &                 'system.'
c        Write(6,'(A)') 'The cartesian Z axis of the system is chosen'//
c     &                 'the quantisation axis for the local '//
c     &                 'pseudospins of the two interacting sites.'
c        Write(6,'(A)')
c
c       iopt=2
c       unity(1:3,1:3)=0.0_wp
c       unity(1,1)=1.0_wp
c       unity(2,2)=1.0_wp
c       unity(3,3)=1.0_wp
c
c       If(Lines.or.AnisoLines) Then
c
c         Call transHam( n1, n2, unity(1:3,1:3), unity(1:3,1:3),
c     &                  MM(i1,1:3,1:n1,1:n1), MM(i2,1:3,1:n2,1:n2),
c     &                  itype(i1), itype(i2),
c     &                  HLIN(lp, 1:n1,1:n1, 1:n2,1:n2),
c     &                 HLIN3(lp, 1:n1,1:n1, 1:n2,1:n2), iopt )
c         JLinG(lp, 1:(n1-1), -(n1-1):(n1-1),
c     &             1:(n2-1), -(n2-1):(n2-1) )=(0.0_wp,0.0_wp)
c         Call JKQPar( n1, n2, HLIN3(lp,1:n1,1:n1,1:n2,1:n2),
c     &                JLinG( lp,1:(n1-1), -(n1-1):(n1-1),
c     &                          1:(n2-1), -(n2-1):(n2-1)) )
c         JLinCG(lp,1:3,1:3)=0.0_wp
c         Call tensor2cart(1,1,JLinG(lp,1,-1:1,1,-1:1),
c     &                       JLinCG(lp,1:3,1:3) )
c        Write(6,'(A)')
c        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
c     &                 '(rank-1) exchange interaction: LINES'
c        Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
c        Write(6,'(A)') 'GENERAL COORD:::'
c        Write(6,'(A)') '     (  xx   xy  xz  )  '
c        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
c        Write(6,'(A)') '     (  zx   zy  zz  )  '
c        Do i=1,3
c          Write(6,'(3ES22.14)') (JLinCG(lp,i,j),j=1,3)
c        End Do
c       End If

c       If (Dipol) Then
c         Call transHam( n1, n2, unity(1:3,1:3), unity(1:3,1:3),
c     &                  MM(i1,1:3,1:n1,1:n1), MM(i2,1:3,1:n2,1:n2),
c     &                  itype(i1), itype(i2),
c     &                  HDIP(lp, 1:n1,1:n1, 1:n2,1:n2),
c     &                 HDIP3(lp, 1:n1,1:n1, 1:n2,1:n2), iopt )
c          JDipG(lp, 1:(n1-1), -(n1-1):(n1-1),
c     &              1:(n2-1), -(n2-1):(n2-1) )=(0.0_wp,0.0_wp)
c         Call JKQPar( n1, n2, HDIP3(lp,1:n1,1:n1,1:n2,1:n2),
c     &                JDipG( lp,1:(n1-1), -(n1-1):(n1-1),
c     &                          1:(n2-1), -(n2-1):(n2-1) ) )
c         JDipCG(lp,1:3,1:3)=0.0_wp
c         Call tensor2cart(1,1,JDipG(lp,1,-1:1,1,-1:1),
c     &                       JDipCG(lp,1:3,1:3) )
c        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
c     &                 '(rank-1) exchange interaction: DIPOL'
c        Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
c        Write(6,'(A)') 'GENERAL COORD:::'
c        Write(6,'(A)') '     (  xx   xy  xz  )  '
c        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
c        Write(6,'(A)') '     (  zx   zy  zz  )  '
c        Do i=1,3
c          Write(6,'(3ES22.14)') (JDipCG(lp,i,j),j=1,3)
c        End Do
c       End If




c        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
c     &                 '(rank-1) exchange interaction: LINES+DIPOL'
c        Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
c        Write(6,'(A)') 'GENERAL COORD:::'
c        Write(6,'(A)') '     (  xx   xy  xz  )  '
c        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
c        Write(6,'(A)') '     (  zx   zy  zz  )  '
c        Do i=1,3
c          Write(6,'(3F24.14)')
c     &            ((JLinCG(lp,i,j)+JDipCG(lp,i,j))*cm_to_MHz,j=1,3)
c        End Do


c print out the data:
c      Write(6,'(A)')
c      Write(6,'(10x,A)') 'Parameters of the ITOs:'
c      Write(6,'( 5x,A)') 'with absolute values larger than:  0.5d-14 '
c      If((Lines.or.AnisoLines).and.(.not.Dipol)) Then
c      Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
c     & '----------------------|'
c      Write(6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exch'//
c     & 'ange  Interaction     |'
c      Write(6,'(A)') '------|------|------|------|------ Real ---'//
c     & '-------- Imag --------|'
c      Do k1=1,N1-1,2
c        Do q1=-k1,k1
c          Do k2=1,N2-1,2
c            Do q2=-k2,k2
c                  If(ABS(JLinG(lp,k1,q1,k2,q2)) .gt. 0.5d-14 ) Then
c      Write(6,'(4(i4,2x,A),2(1x,E17.10),1x,A)')
c     & k1,'|',q1,'|',k2,'|',q2,'|',JLinG(lp,k1,q1,k2,q2),'|'
c                  End If
c            End Do
c          End Do
c        End Do
c      End Do
c      Write(6,'(A)') '------|------|------|------|---------------'//
c     & '----------------------|'
c      Else If((.not.(Lines.or.AnisoLines)).and.Dipol) Then
c      Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
c     & '----------------------|'
c      Write(6,'(A)') ' rank | proj.| rank | proj.|  Dipolar  Exch'//
c     & 'ange  Interaction     |'
c      Write(6,'(A)') '------|------|------|------|------ Real ---'//
c     & '-------- Imag --------|'
c      Do k1=1,N1-1,2
c        Do q1=-k1,k1
c          Do k2=1,N2-1,2
c            Do q2=-k2,k2
c                  If(ABS(JDipG(lp,k1,q1,k2,q2)) .gt. 0.5d-14 ) Then
c      Write(6,'(4(i4,2x,A),2(1x,E17.10),1x,A)')
c     & k1,'|',q1,'|',k2,'|',q2,'|',JDipG(lp,k1,q1,k2,q2),'|'
c                  End If
c            End Do
c          End Do
c        End Do
c      End Do
c      Write(6,'(A)') '------|------|------|------|---------------'//
c     & '----------------------|'
c      Else If((Lines.or.AnisoLines).and.Dipol) Then
c      Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
c     & '----------------------|-------------------------------------|'
c      Write(6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exch'//
c     & 'ange  Interaction     |  Dipolar  Exchange  Interaction     |'
c      Write(6,'(A)') '------|------|------|------|------ Real ---'//
c     & '-------- Imag --------|------ Real ----------- Imag --------|'
c      Do k1=1,N1-1,2
c        Do q1=-k1,k1
c          Do k2=1,N2-1,2
c            Do q2=-k2,k2
c                  If ( (ABS(JLinG(lp,k1,q1,k2,q2)) .gt. 0.5d-14) .OR.
c     &                 (ABS(JDipG(lp,k1,q1,k2,q2)) .gt. 0.5d-14) )  Then
c      Write(6,'(4(i4,2x,A),2(1x,E17.10),1x,A,2(1x,E17.10),1x,A)')
c     & k1,'|',q1,'|',k2,'|',q2,'|',JLinG(lp,k1,q1,k2,q2),'|',
c     & JDipG(lp,k1,q1,k2,q2),'|'
c                  End If
c            End Do
c          End Do
c        End Do
c      End Do
c      Write(6,'(A)') '------|------|------|------|---------------'//
c     & '----------------------|-------------------------------------|'
c      End If
c
c
c
c
c
c
!      Write(6,'(A)')
!      Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
!     & '(rank-1) exchange interaction'
!      Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
!      Write(6,'(A)') '     (  xx   xy  xz  )  '
!      Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
!      Write(6,'(A)') '     (  zx   zy  zz  )  '
!
!      express the rank-1 tensors in the sum of
!      Isotrop part    C * unit matrix
!      Symmetric part
!      Antisimmetric part
!      print out the data:
!
!      JDipCG(lp,:,:)=0.0_wp
!      If(Lines.or.AnisoLines) Then
!      Write(6,'(A)') 'Lines matrix:'
!      Do i=1,3
!      Write(6,'(3F12.6)') (JLinCG(lp,i,j),j=1,3)
!      End Do
!
!
!      ELin(lp,:,:)=0.0_wp
!      ALin(lp,:,:)=0.0_wp
!      SLin(lp,:,:)=0.0_wp
!      ELin(lp,1,1)=(JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
!      ELin(lp,2,2)=(JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
!      ELin(lp,3,3)=(JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
!      Do is1=1,3
!       Do is2=1,3
!       ALin(lp,is1,is2)=(JLinCG(lp,is1,is2)-JLinCG(lp,is2,is1))/2.0_wp
!       End Do
!      End Do
!      tsum=0.0_wp
!      Do l=1,3
!      tsum=tsum+JLinCG(lp,l,l)
!      End Do
!      test=0.0_wp
!      Do is1=1,3
!      test(is1,is1)=2.0_wp*tsum/3.0_wp
!      End Do
!
!      Do is1=1,3
!       Do is2=1,3
!       SLin(lp,is1,is2)=( JLinCG(lp,is1,is2)+JLinCG(lp,is2,is1)
!     &                   -test(is1,is2) )/2.0_wp
!       End Do
!      End Do
!      Write(6,'(A)') 'Elin * unit'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (ELin(lp,is1,is2),is2=1,3)
!      End Do
!      Write(6,'(A)') '-2/3 * test'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (-2.0_wp*test(is1,is2)/3.0_wp,is2=1,3)
!      End Do
!
!      Write(6,'(A)') 'ALin:'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (ALin(lp,is1,is2),is2=1,3)
!      End Do
!      Write(6,'(A)') 'SLin:'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (SLin(lp,is1,is2),is2=1,3)
!      End Do
!      End If
!
!      If(Dipol) Then
!      Write(6,'(A)') 'Dipolar exchange matrix:'
!      Do i=1,3
!      Write(6,'(3F12.6)') (JDipCG(lp,i,j),j=1,3)
!      End Do
!
!      EDip(lp,:,:)=0.0_wp
!      ADip(lp,:,:)=0.0_wp
!      SDip(lp,:,:)=0.0_wp
!      EDip(lp,1,1)=(JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
!      EDip(lp,2,2)=(JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
!      EDip(lp,3,3)=(JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
!      Do is1=1,3
!       Do is2=1,3
!       ADip(lp,is1,is2)=(JDipCG(lp,is1,is2)-JDipCG(lp,is2,is1))/2.0_wp
!       End Do
!      End Do
!      tsum=0.0_wp
!      Do l=1,3
!      tsum=tsum+JDipCG(lp,l,l)
!      End Do
!      test=0.0_wp
!      Do is1=1,3
!      test(is1,is1)=2.0_wp*tsum/3.0_wp
!      End Do
!      Do is1=1,3
!       Do is2=1,3
!       SDip(lp,is1,is2)=( JDipCG(lp,is1,is2)+JDipCG(lp,is2,is1)
!     &                   -test(is1,is2) )/2.0_wp
!       End Do
!      End Do
!      Write(6,'(A)') 'EDip * unit'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (EDip(lp,is1,is2),is2=1,3)
!      End Do
!      Write(6,'(A)') '-2/3 * test'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (-2.0_wp*test(is1,is2)/3.0_wp,is2=1,3)
!      End Do
!      Write(6,'(A)') 'ADip:'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (ADip(lp,is1,is2),is2=1,3)
!      End Do
!      Write(6,'(A)') 'SDip:'
!      Do is1=1,3
!      Write(6,'(3F12.6)') (SDip(lp,is1,is2),is2=1,3)
!      End Do
!      End If
!
        Call mma_deallocate(JN)
        Call mma_deallocate(JB)
        Call mma_deallocate(JS)
      End Do ! lp, interacting pairs


      Return
      End
