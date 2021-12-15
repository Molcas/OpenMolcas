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
      Subroutine pa_prham( exch, npair, i_pair, nneq, neq, nexch,
     &                     nmax, lmax, eso, HLIN1, HLIN3, HLIN9,
     &                     HDIP, HKEX, HDMO, HITO, Dipol, AnisoLines1,
     &                     AnisoLines3, AnisoLines9, KE, DM_Exchange,
     &                     JITO_exchange )
c this function prints the exchange Hamiltonian
c it does not compute any new infromation
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
#include "warnings.fh"
      Integer, intent(in)          :: exch
      Integer, intent(in)          :: npair
      Integer, intent(in)          :: i_pair(npair,2)
      Integer, intent(in)          :: nneq
      Integer, intent(in)          :: neq(nneq)
      Integer, intent(in)          :: nexch(nneq)
      Integer, intent(in)          :: nmax
      Integer, intent(in)          :: lmax
      Real(kind=8), intent(in)    :: eso(nneq,nmax)
      Complex(kind=8), intent(in) :: HLIN1(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HLIN3(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HLIN9(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HDIP(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HKEX(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HDMO(npair,nmax,nmax,nmax,nmax)
      Complex(kind=8), intent(in) :: HITO(npair,nmax,nmax,nmax,nmax)
      Logical, intent(in)          :: Dipol
      Logical, intent(in)          :: DM_exchange
      Logical, intent(in)          :: KE
      Logical, intent(in)          :: AnisoLines1
      Logical, intent(in)          :: AnisoLines3
      Logical, intent(in)          :: AnisoLines9
      Logical, intent(in)          :: JITO_exchange
c local variables
      Complex(kind=8), allocatable :: HTOT(:), H1(:), H2(:), H3(:),
     &                                 H4(:)

      Integer  :: nind(lmax,2), intc(lmax), ibas(exch,lmax),
     &            icoord(lmax),
     &            nb1, nb2, lb1, lb2, i1, i2, is1, is2,
     &            js1, js2, nb, i, j, l, lp, lb, lpr, ibuf
      Integer  :: CtoB, RtoB, ItoB, mem_local
      Integer  :: norder
      Real(kind=8) :: dznrm2_
      External :: norder, dznrm2_

!=======================================================================
      If(npair==0) Then
          Call WarningMessage(2,'PA_PRHAM: npair = 0')
          Return
      End If
      If(nmax==0)  Then
          Call WarningMessage(2,'PA_PRHAM:  nmax = 0')
          Return
      End If
      If(exch==0)  Then
          Call WarningMessage(2,'PA_PRHAM:  exch = 0')
          Return
      End If
      If(lmax==0)  Then
          Call WarningMessage(2,'PA_PRHAM:  lmax = 0')
          Return
      End If
      ibuf=npair*nmax*nmax*nmax*nmax
      If(ibuf>0) Then
        If( (dznrm2_(ibuf,HLIN1,1)==0.0_wp).and.(AnisoLines1) )
     &      Call WarningMessage(2,'PA_PRHAM:  HLIN1 is empty')
        If( (dznrm2_(ibuf,HLIN3,1)==0.0_wp).and.(AnisoLines3) )
     &      Call WarningMessage(2,'PA_PRHAM:  HLIN3 is empty')
        If( (dznrm2_(ibuf,HLIN9,1)==0.0_wp).and.(AnisoLines9) )
     &      Call WarningMessage(2,'PA_PRHAM:  HLIN9 is empty')
        If( (dznrm2_(ibuf,HDIP,1)==0.0_wp).and.(Dipol) )
     &      Call WarningMessage(2,'PA_PRHAM:  HDIP is empty')
        If( (dznrm2_(ibuf,HKEX,1)==0.0_wp).and.(KE) )
     &      Call WarningMessage(2,'PA_PRHAM:  HKEX is empty')
        If( (dznrm2_(ibuf,HDMO,1)==0.0_wp).and.(DM_exchange) )
     &      Call WarningMessage(2,'PA_PRHAM:  HDMO is empty')
        If( (dznrm2_(ibuf,HITO,1)==0.0_wp).and.(JITO_exchange) )
     &      Call WarningMessage(2,'PA_PRHAM:  HITO is empty')
      End If !ibuf
!=======================================================================
! allocate memory
      ItoB=8
      RtoB=8
      CtoB=16
      mem_local=0
      If(exch>=0) Then
        Call mma_allocate(htot,exch,'htot')
        Call mma_allocate(h1,exch,'h1')
        Call mma_allocate(h2,exch,'h2')
        Call mma_allocate(h3,exch,'h3')
        Call mma_allocate(h4,exch,'h4')
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,htot,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,h1,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,h2,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,h3,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,h4,1)
        mem_local=mem_local+5*exch*CtoB
      End If

      !do lp=1,nPair
      !  Write(6,'(A,i2,A,i3)') 'i_Pair(',lp,',1)=',i_pair(lp,1)
      !  Write(6,'(A,i2,A,i3)') 'i_Pair(',lp,',2)=',i_pair(lp,2)
      !end do
c generate the tables:
      nind(:,:)=0
      lpr=0
      l=0
      Do i=1,nneq
        Do j=1,neq(i)
          l=l+1
          nind(l,1)=i
          nind(l,2)=j
        End Do
      End Do
      intc(:)=0
      ibas(:,:)=0
      intc(1)=1
      If (lmax.gt.1) Then
        Do i=2,lmax
          i1=nind(i-1, 1)
          intc(i)=intc(i-1)*nexch(i1)
        End Do
      End If
      Do nb=1,exch
        nb1=nb-1
        Do i=1,lmax
          ibas(nb, lmax-i+1)= nb1 / intc(lmax-i+1)
          nb1=nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
        End Do
      End Do
c
      Write(6,*)
      Write(6,'(100a)') (('%'),j=1,100)
      Write(6,'(30x,a)') 'Hamiltonian of the Total Magnetic '//
     &                   'Interaction'
      Write(6,'(100a)') (('%'),j=1,100)


      Write(6,'(A)') 'Only matrix elements '//
     & 'with non-zero absolute value are listed below'

      If ((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.Dipol
     &                             .and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &        '|       Lines  Model of Interaction   | ',
     &        '|      Dipolar Magnetic Interaction   | ',
     &        '|      Kinetic Exchange Interaction   | ',
     &        '|       TOTAL  Magnetic Interaction   | '
        lpr=4
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.Dipol
     &                        .and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and..not.Dipol
     &                             .and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.Dipol
     &                             .and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If( (AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &              .not.Dipol.and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=2
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                   Dipol.and..not.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=2
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                   .not.Dipol.and.KE.and..not.JITO_exchange ) Then
        Write(6,'(6A)') '                 ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=2
! JITO is active below:
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                    .not.Dipol.and..not.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=2
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         Dipol.and..not.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         .not.Dipol.and.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                              Dipol.and.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=4
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                    .not.Dipol.and..not.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=3
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         Dipol.and..not.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=4
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         .not.Dipol.and.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=4
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                              Dipol.and.KE.and.JITO_exchange) Then
        Write(6,'(6A)') '                 ',
     &   '|       Lines  Model of Interaction   | ',
     &   '|      Dipolar Magnetic Interaction   | ',
     &   '|      Kinetic Exchange Interaction   | ',
     &   '|          ITO Exchange Interaction   | ',
     &   '|       TOTAL  Magnetic Interaction   | '
        lpr=5

      End If
      Write(6,'(110A)') '-----------------',
     & ('|-----  Real  ---------  Imaginary  --| ',i=1,lpr)

      Do nb1 = 1,exch
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,  H1,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,  H2,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,  H3,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,  H4,1)
        Call zcopy_(exch,[(0.0_wp,0.0_wp)],0,HTOT,1)
        Do lp = 1,npair
          icoord(:)=0
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

              H1(nb2)=H1(nb2)+HLIN1(lp,is1,js1,is2,js2)
     &                       +HLIN3(lp,is1,js1,is2,js2)
     &                       +HLIN9(lp,is1,js1,is2,js2)

              H2(nb2)=H2(nb2)+HDIP(lp,is1,js1,is2,js2)

              H3(nb2)=H3(nb2)+HKEX(lp,is1,js1,is2,js2)

              H4(nb2)=H4(nb2)+HITO(lp,is1,js1,is2,js2)

              HTOT(nb2)=HTOT(nb2)+HLIN1(lp,is1,js1,is2,js2)
     &                           +HLIN3(lp,is1,js1,is2,js2)
     &                           +HLIN9(lp,is1,js1,is2,js2)
     &                            +HDIP(lp,is1,js1,is2,js2)
     &                            +HKEX(lp,is1,js1,is2,js2)
     &                            +HDMO(lp,is1,js1,is2,js2)
     &                            +HITO(lp,is1,js1,is2,js2)

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
                HTOT(nb1)=HTOT(nb1)+
     &                          cmplx(eso(i,ibas(nb1,lb)+1),0.0_wp,wp)
                If((lb+1).le.(lmax)) lb=lb+1
              End If
            End Do !j
          End Do !i
        End If !.not.KE
c all matrices for HEXCH( Nb1, xxx) are known:
c proceed to print
        If ((AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &                                        .and.Dipol.and.KE) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H2(nb2)).gt.0.0_wp).or.
     &          (abs(H3(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                            Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H2(nb2), H3(nb2), HTOT(nb2)
            End If
          End Do

        Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.Dipol
     &                                               .and..not.KE ) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H2(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H2(nb2), HTOT(nb2)
            End If
          End Do

        Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                                         .not.Dipol .and.KE ) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H3(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H3(nb2), HTOT(nb2)
            End If
          End Do

        Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                                               Dipol.and.KE ) Then
          Do nb2=nb1,exch
            If( (abs(H2(nb2)).gt.0.0_wp).or.(abs(H3(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H2(nb2), H3(nb2), HTOT(nb2)
            End If
          End Do

        Else If( (AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                                    .not.Dipol.and..not.KE ) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                            Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), HTOT(nb2)
            End If
          End Do

        Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                                        Dipol.and..not.KE ) Then
          Do nb2=nb1,exch
            If( (abs(H2(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                            Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H2(nb2), HTOT(nb2)
            End If
          End Do

        Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9)
     &                                  .and..not.Dipol.and. KE ) Then
          Do nb2=nb1,exch
            If( (abs(H3(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                            Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H3(nb2), HTOT(nb2)
            End If
          End Do
! JITO is active below:
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                    .not.Dipol.and..not.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H4(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                            Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         Dipol.and..not.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H2(nb2)).gt.0.0_wp).or.(abs(H4(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H2(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         .not.Dipol.and.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H3(nb2)).gt.0.0_wp).or.(abs(H4(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H3(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If(.not.(AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                              Dipol.and.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H2(nb2)).gt.0.0_wp).or.(abs(H3(nb2)).gt.0.0_wp).or.
     &          (abs(H4(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                              Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H2(nb2), H3(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                    .not.Dipol.and..not.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H4(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         Dipol.and..not.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H2(nb2)).gt.0.0_wp).or.
     &          (abs(H4(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                              Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H2(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                         .not.Dipol.and.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H3(nb2)).gt.0.0_wp).or.
     &          (abs(H4(nb2)).gt.0.0_wp).or.(abs(HTOT(nb2)).gt.0.0_wp) )
     &                                                              Then
               Write(6,'(A,i4,A,i4,A,4(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H3(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      Else If((AnisoLines1.or.AnisoLines3.or.AnisoLines9).and.
     &                              Dipol.and.KE.and.JITO_exchange) Then
          Do nb2=nb1,exch
            If( (abs(H1(nb2)).gt.0.0_wp).or.(abs(H2(nb2)).gt.0.0_wp).or.
     &          (abs(H3(nb2)).gt.0.0_wp).or.(abs(H4(nb2)).gt.0.0_wp).or.
     &          (abs(HTOT(nb2)).gt.0.0_wp) ) Then
               Write(6,'(A,i4,A,i4,A,5(2E19.11,2x))')
     &                  '<',nb1,'| H |',nb2,'> =',
     &                  H1(nb2), H2(nb2), H3(nb2), H4(nb2), HTOT(nb2)
            End If
          End Do
      End If

      End Do !nb1

      Write(6,'(10A)') '-----------------',
     & ('|-------------------------------------| ',i=1,lpr)

!=======================================================================
! deallocate memory
      If(exch>=0) Then
        Call mma_deallocate(htot)
        Call mma_deallocate(h1)
        Call mma_deallocate(h2)
        Call mma_deallocate(h3)
        Call mma_deallocate(h4)
      End If


      Return
      End
