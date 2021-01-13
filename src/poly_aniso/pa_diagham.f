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
       Subroutine pa_diagham( exch, npair, i_pair,  nneq,   neq, nexch,
     &                        nmax,  lmax,    eso,  HLIN1, HLIN3, HLIN9,
     &                        HDIP,  HKEX, HDMO, HITO, Dipol,
     &                        DM_exchange, AnisoLines1, AnisoLines3,
     &                        AnisoLines9, KE, JITO_exchange,
     &                        WLIN1, WLIN3, WLIN9, WLIN, WDIP,
     &                        WKEX, WDMO, WITO,  W,     Z  )
c this function builds and diagonalizes the interaction Hamiltonians
      Implicit None
#include "stdalloc.fh"
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: exch
      Integer, intent(in)           :: npair
      Integer, intent(in)           :: i_pair(npair,2)
      Integer, intent(in)           :: nneq
      Integer, intent(in)           :: neq(nneq)
      Integer, intent(in)           :: nexch(nneq)
      Integer, intent(in)           :: nmax
      Integer, intent(in)           :: lmax
      Logical, intent(in)           :: Dipol
      Logical, intent(in)           :: KE
      Logical, intent(in)           :: DM_exchange
      Logical, intent(in)           :: AnisoLines1
      Logical, intent(in)           :: AnisoLines3
      Logical, intent(in)           :: AnisoLines9
      Logical, intent(in)           :: JITO_exchange
      Real(kind=8), intent(in)     :: eso(nneq,nmax)
      Complex(kind=8), intent(in)  :: HLIN1(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HLIN3(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HLIN9(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HDIP(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HKEX(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HDMO(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in)  :: HITO(npair,nmax,nmax,nmax,nmax)
      ! output data:
      Real(kind=8), intent(out)    :: wlin(exch) ! total 1+3+9
      Real(kind=8), intent(out)    :: wlin1(exch)
      Real(kind=8), intent(out)    :: wlin3(exch)
      Real(kind=8), intent(out)    :: wlin9(exch)
      Real(kind=8), intent(out)    :: wdip(exch)
      Real(kind=8), intent(out)    :: wkex(exch)
      Real(kind=8), intent(out)    :: wdmo(exch)
      Real(kind=8), intent(out)    :: wito(exch)

      Real(kind=8), intent(out)    :: w(exch)
      Complex(kind=8), intent(out) :: Z(exch,exch)
c local variables
      Complex(kind=8), allocatable :: HTOT(:,:)
      Integer, allocatable :: nind(:,:), intc(:), ibas(:,:), icoord(:)
      Integer  :: nb1, nb2, lb1, lb2, i1, i2, is1, is2,
     &            js1, js2, nb, i, j, l, lp, lb
      Integer  :: norder
      external :: norder
c diag:
      Integer          :: info, lwork
      Real(kind=8), allocatable :: rwork(:) !rwork(3*exch-2)
      Complex(kind=8), allocatable :: work(:) !work(2*exch-1)
c allocate memory and initialize variables:
      If( exch >= 0 ) Then
         Call mma_allocate(HTOT,exch,exch,'HTOT')
         Call mma_allocate(WORK,(2*exch-1),'WORK')
         If( lmax >= 0 ) Then
            Call mma_allocate(ibas,exch,lmax,'ibas')
         End If
         Call mma_allocate(rwork,(3*exch-2),'rwork')
      End If

      If( lmax >= 0 ) Then
         Call mma_allocate(nind,lmax,2,'nind')
         Call mma_allocate(intc,lmax,'intc')
         Call mma_allocate(icoord,lmax,'icoord')
      End If
      Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,Z,1)
      Call dcopy_(exch,[0.0_wp],0,w,1)
      Call dcopy_(exch,[0.0_wp],0,wlin,1)
      Call dcopy_(exch,[0.0_wp],0,wlin1,1)
      Call dcopy_(exch,[0.0_wp],0,wlin3,1)
      Call dcopy_(exch,[0.0_wp],0,wlin9,1)
      Call dcopy_(exch,[0.0_wp],0,wdip,1)
      Call dcopy_(exch,[0.0_wp],0,wkex,1)
      Call dcopy_(exch,[0.0_wp],0,wdmo,1)
      Call dcopy_(exch,[0.0_wp],0,wito,1)

c generate the tables:
      l=0
      Call icopy(2*lmax,[0],0,nind,1)
      Do i=1,nneq
        Do j=1,neq(i)
          l=l+1
          nind(l,1)=i
          nind(l,2)=j
        End Do
      End Do

      Call icopy(lmax,[0],0,intc,1)
      intc(1)=1
      If (lmax.gt.1) Then
        Do i=2,lmax
          i1=nind(i-1, 1)
          intc(i)=intc(i-1)*nexch(i1)
        End Do
      End If

      Call icopy(exch*lmax,[0],0,ibas,1)
      Do nb=1,exch
        nb1=nb-1
        Do i=1,lmax
          ibas(nb, lmax-i+1)= nb1 / intc(lmax-i+1)
          nb1=nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
        End Do
      End Do
c build the interaction Hamiltonians
!----------------------------------------------------------------------!
      If( AnisoLines1 ) Then
        Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HLIN1(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp

          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                HTOT(nb1,nb1)=HTOT(nb1,nb1) +
     &                        cmplx( eso(i,ibas(nb1,lb)+1), 0.0_wp, wp )
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wlin1,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If( AnisoLines3 ) Then
        Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HLIN3(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp

          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                HTOT(nb1,nb1)=HTOT(nb1,nb1) +
     &                        cmplx( eso(i,ibas(nb1,lb)+1), 0.0_wp, wp )
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wlin3,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If( AnisoLines9 ) Then
        Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HLIN9(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp

          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                HTOT(nb1,nb1)=HTOT(nb1,nb1) +
     &                        cmplx( eso(i,ibas(nb1,lb)+1), 0.0_wp, wp )
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wlin9,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If( AnisoLines1.OR.AnisoLines3.OR.AnisoLines9 ) Then
        Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HLIN1(lp,is1,js1,is2,js2)
     &                                      + HLIN3(lp,is1,js1,is2,js2)
     &                                      + HLIN9(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp

          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                HTOT(nb1,nb1)=HTOT(nb1,nb1) +
     &                        cmplx( eso(i,ibas(nb1,lb)+1), 0.0_wp, wp )
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wlin,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If(Dipol) Then
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HDIP(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp
          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                ! kind=8, complex double precision
                HTOT(nb1,nb1)=HTOT(nb1,nb1)+
     &                        cmplx(eso(i,ibas(nb1,lb)+1),0.0_wp,wp)
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wdip,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If( DM_exchange ) Then
        Call zcopy_(exch*exch, [(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HDMO(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp

          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                HTOT(nb1,nb1)=HTOT(nb1,nb1) +
     &                        cmplx( eso(i,ibas(nb1,lb)+1), 0.0_wp, wp )
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wdmo,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If(KE) Then
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HKEX(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wkex,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      If(JITO_exchange) Then
        Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,HTOT,1)
        Do nb1 = 1,exch
          Do lp = 1,npair
            Call icopy(lmax,[0],0,icoord,1)
            Do i = 1,lmax
              icoord(i) = ibas(nb1,i)
            End Do

            lb1 = i_pair(lp ,1)
            lb2 = i_pair(lp ,2)
            i1  =   nind(lb1,1)
            i2  =   nind(lb2,1)
            is1 = icoord(lb1) + 1
            is2 = icoord(lb2) + 1
            Do js1 = 1, nexch(i1)
              icoord(lb1) = js1 - 1
              Do js2 = 1, nexch(i2)
                icoord(lb2) = js2 - 1
                nb2 = norder(icoord,intc,lmax)
                HTOT(nb1,nb2)=HTOT(nb1,nb2) + HITO(lp,is1,js1,is2,js2)
              End Do ! js2
            End Do ! js1
          End Do ! lp
        End Do !nb1
        ! diagonalize
        info=0
        lwork=0
        lwork=2*exch-1
        Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
        Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
        Call zheev('n','u',exch,htot,exch,wito,work,lwork,rwork,info)
      End If



!----------------------------------------------------------------------!
      !cccccccc  total Hamiltonian cccccccc
      Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,HTOT,1)
      Do nb1 = 1,exch
        Do lp = 1,npair
          Call icopy(lmax,[0],0,icoord,1)
          Do i = 1,lmax
            icoord(i) = ibas(nb1,i)
          End Do

          lb1 = i_pair(lp ,1)
          lb2 = i_pair(lp ,2)
          i1  =   nind(lb1,1)
          i2  =   nind(lb2,1)
          is1 = icoord(lb1) + 1
          is2 = icoord(lb2) + 1
          Do js1 = 1, nexch(i1)
            icoord(lb1) = js1 - 1
            Do js2 = 1, nexch(i2)
              icoord(lb2) = js2 - 1
              nb2 = norder(icoord,intc,lmax)

              HTOT(nb1,nb2)= HTOT(nb1,nb2)
     &                     + HLIN1(lp,is1,js1,is2,js2)
     &                     + HLIN3(lp,is1,js1,is2,js2)
     &                     + HLIN9(lp,is1,js1,is2,js2)
     &                     +  HDIP(lp,is1,js1,is2,js2)
     &                     +  HKEX(lp,is1,js1,is2,js2)
     &                     +  HITO(lp,is1,js1,is2,js2)
            End Do ! js2
          End Do ! js1
        End Do ! lp

        If(.not.KE) Then
          l=0
          lb=1
          Do i = 1, nneq
            Do j = 1, neq(i)
              l = l + 1
              If( l.eq.lb ) Then
                !kind=8, complex double precision
                HTOT(nb1,nb1)=HTOT(nb1,nb1)+
     &                        cmplx(eso(i,ibas(nb1,lb)+1),0.0_wp,wp)
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End If ! .not. KE
      End Do !nb1
c diagonalize
      info =0
      lwork=0
      lwork=2*exch-1
      Call dcopy_((3*exch-2),[0.0_wp],0,rwork,1)
      Call zcopy_((2*exch-1),[(0.0_wp,0.0_wp)],0,WORK,1)
      Call zheev('v','u',exch,htot,exch,w,work,lwork,rwork,info)
      If (info.eq.0) Then
         Call zcopy_(exch*exch,htot,1,Z,1)
      Else
         Write(6,'(A,i10)') 'DIAG:  non-zero Return: INFO=', info
      End If

!----------------------------------------------------------------------!
! deallocate memory:
      If( exch >= 0 ) Then
         Call mma_deallocate(HTOT)
         Call mma_deallocate(WORK)
         If( lmax >= 0 ) Then
            Call mma_deallocate(ibas)
         End If
         Call mma_deallocate(rwork)
      End If

      If( lmax >= 0 ) Then
         Call mma_deallocate(nind)
         Call mma_deallocate(intc)
         Call mma_deallocate(icoord)
      End If
      Return
      End
