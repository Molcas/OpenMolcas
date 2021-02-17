************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2009, Giovanni Ghigo                                   *
************************************************************************
C!-----------------------------------------------------------------------!
C!
C!  InterSystem Crossing rate evaluation: Multidimensional Franck-Condon
C!  Modified copy of  FCval by Giovanni Ghigo.
C!  Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
C!  28-Dec-08 - 06-Jan-09 ; June 2009
C!
      Subroutine ISC_FCval(iPrint,        iMaxYes,  nTabDim,
     &  C1,        W1,        det1,       r01,      C2,
     &  W2,        det2,      r02,                  max_mOrd,
     &  max_nOrd,  max_nOrd2, max_mInc,   max_nInc, max_nInc2,
     &  mMat,      nMat,      mInc,       nInc,     mDec,
     &  nDec,      C,         W,          det0,     r00,
     &  L,         U,         FC00,       Alpha1,   Alpha2,
     &  Beta,      nOsc,      nnsiz,      iMx_nOrd,
     &  nYes,      VibWind2,  FCWind2)
C!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "dims.fh"
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),W1(nOsc,nOsc)
      Real*8 W2(nOsc,nOsc),C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8 L (0:max_mOrd,0:max_nInc2)
      Real*8 U (0:max_nOrd,0:max_nOrd2)
      Real*8 Alpha1(nOsc,nOsc),Alpha2(nOsc,nOsc),Beta(nOsc,nOsc)
      Real*8 r00(nOsc),r01(nOsc),r02(nOsc)
      Real*8 FCWind2(nYes)
      Integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2),
     &                             mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2), nInc(0:iMaxYes,nOsc),
     &                             nDec(0:iMaxYes,nOsc)
      Integer VibWind2(nYes)
#include "WrkSpc.fh"
#include "inout.fh"
c      Write(6,*) 'CGGt[ISC_FCval] Enter '                         ! CGGt
c      Write(6,*) '                iMx_nOrd, nYes = ',iMx_nOrd, nYes ! Ct
c      Write(6,*) '                VibWind2 :',(VibWind2(i),i=1,nYes)  !t
c      Write(6,*) '            L matrix:',max_mOrd,max_nInc2       ! CGGt
c      Write(6,*) '            U matrix:',max_nOrd,max_nOrd2       ! CGGt
c      Do i = 0,iMaxYes
c      Write(6,*) (nInc(i,j),j=1,nOsc)
c      EndDo
c      Call XFlush(6)                                              ! CGGt
C!
C!---- Initialize.
      nMaxMat = max(max_mOrd+1,max_nOrd+1)
c      Write(6,*) '            nMaxMat=',nMaxMat                   ! CGGt
      nTabDim = max(nMaxMat,8)
c      Write(6,*) '            nTabDim=',nTabDim                   ! CGGt
      nOscSqr = nOsc**2
c      Write(6,*) '            nOscSqr=',nOscSqr                   ! CGGt
      Call GetMem('temp','Allo','Real',iptemp,nOscSqr)
      Call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
      Call GetMem('temp2','Allo','Real',iptemp2,nOscSqr)

C!
C!---- Setup sqr table.
      n=nTabDim+1
      Call GetMem('sqr','Allo','Real',ipsqr,n+1)
      Do i = 0,nTabDim+1
        Work(ipsqr+i) = sqrt(dble(i))
      End Do
C!
C!---- Calculate alpha1, alpha2 and alpha.
c      Write(6,*) 'CGGt[FCVal] Calculate alpha(s)'                 ! CGGt
c      Call XFlush(6)                                              ! CGGt
      Call GetMem('Alpha','Allo','Real',ipAlpha,nOscSqr)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C1,nOsc,
     &            C1,nOsc,
     &            0.0d0,Alpha1,nOsc)
      call dscal_(nOscSqr,0.5d0,Alpha1,1)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C2,nOsc,
     &            C2,nOsc,
     &            0.0d0,Alpha2,nOsc)
      call dscal_(nOscSqr,0.5d0,Alpha2,1)
c       temp = alpha1+alpha2
      call dcopy_(nOscSqr,Alpha1,1,Work(iptemp),1)
      Call Daxpy_(nOscSqr,1.0d0,Alpha2,1,Work(iptemp),1)
c       alpha = 0.5d0*temp
      call dcopy_(nOscSqr,[0.0d0],0,Work(ipAlpha),1)
      Call Daxpy_(nOscSqr,0.5d0,Work(iptemp),1,Work(ipAlpha),1)

C! Call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
C! call dscal_(nOscSqr,0.5d0,alpha,1)
C!
C!---- Calculate C using a Cholesky factorization of 2*alpha.
C! Call Cholesky(temp,C)
C!
C!---- Calculate W.
C! call dcopy_(nOscSqr,[0.0d0],0,W,1) ; call dcopy_(nOsc,[1.0d0],0,W,nOsc+1)
C! temp = C ; Call Dool_MULA(temp,W,det0)
C!
C!---- Calculate r00.
      Call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
      Call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)
C!
C!---- Calculate beta.
c      Write(6,*) 'CGGt[FCVal] Calculate beta.'                    ! CGGt
c      Call XFlush(6)                                              ! CGGt
      Do i = 1,nOsc
      Do j = 1,nOsc
      Work(iptemp1+j+nOsc*(i-1)-1) = C1(i,j)
      End Do
c      Write(6,*) 'CGGt C1(',i,',j)=',(C1(i,jj),jj=1,nOsc)         ! CGGt
      End Do
c      Call XFlush(6)                                              ! CGGt
c       temp1 = alpha1
      call dcopy_(nOscSqr,Alpha1,1,Work(iptemp1),1)
c       temp  = 2.0d0*alpha
      call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
      Call Daxpy_(nOscSqr,2.0d0,Work(ipAlpha),1,Work(iptemp),1)

      Call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Alpha2,nOsc,
     &            Work(iptemp1),nOsc,
     &            0.0d0,Beta,nOsc)
C!
C!---- Calculate FC00.
c       r_temp1 = r02-r01
      call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
      Call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,Beta,nOsc,
     &            Work(ipr_temp1),nOsc,
     &            0.0d0,Work(ipr_temp2),nOsc)
      FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
      FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
c      Write(6,*) 'CGGt[FCVal] FC00_exp,FC00=',FC00_exp,FC00       ! CGGt
c      Call XFlush(6)                                              ! CGGt
C!
C!---- Calculate A, B and d matrices.
      Call GetMem('A1','Allo','Real',ipA1,nOscSqr)
      Call GetMem('B1','Allo','Real',ipB1,nOscSqr)
      Call GetMem('A2','Allo','Real',ipA2,nOscSqr)
      Call GetMem('B2','Allo','Real',ipB2,nOscSqr)
      Call GetMem('d1','Allo','Real',ipd1,nOsc)
      Call GetMem('d2','Allo','Real',ipd2,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C1,nOsc,
     &            W,nOsc,
     &            0.0d0,Work(ipA1),nOsc)
      Call DGEMM_('T','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,W1,nOsc,
     &            C,nOsc,
     &            0.0d0,Work(iptemp),nOsc)
c       B1 = A1-temp
      call dcopy_(nOscSqr,Work(ipA1),1,Work(ipB1),1)
      Call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB1),1)


      Call DGEMM_('T','N',
     &            nOsc,1,nOsc,
     &            1.0d0,W1,nOsc,
     &            Work(ipr_temp2),nOsc,
     &            0.0d0,Work(ipd1),nOsc)
      const = Work(ipsqr+8)
      call dscal_(nOsc,const,Work(ipd1),1)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C2,nOsc,
     &            W,nOsc,
     &            0.0d0,Work(ipA2),nOsc)
      Call DGEMM_('T','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,W2,nOsc,
     &            C,nOsc,
     &            0.0d0,Work(iptemp),nOsc)
c       B2 = A2-temp
      call dcopy_(nOscSqr,Work(ipA2),1,Work(ipB2),1)
      Call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB2),1)

      Call DGEMM_('T','N',
     &            nOsc,1,nOsc,
     &            1.0d0,W2,nOsc,
     &            Work(ipr_temp2),nOsc,
     &            0.0d0,Work(ipd2),nOsc)
      const =-Work(ipsqr+8)
      call dscal_(nOsc,const,Work(ipd2),1)
C!
C!---- Calculate A1B1T and A2B2T.
      Call GetMem('A1B1T','Allo','Real',ipA1B1T,nOscSqr)
      Call GetMem('A2B2T','Allo','Real',ipA2B2T,nOscSqr)

      Call DGEMM_('N','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Work(ipA1),nOsc,
     &            Work(ipB1),nOsc,
     &            0.0d0,Work(ipA1B1T),nOsc)
      Call DGEMM_('N','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Work(ipA2),nOsc,
     &            Work(ipB2),nOsc,
     &            0.0d0,Work(ipA2B2T),nOsc)
C!

C!
C!---- Initialize L matrix.
      call dcopy_((max_mord+1)*(max_ninc2+1),[0.0d0],0,L,1)
c       L = 0.0d0
      L(0,0) = 1.0d0
C!

C!---- If max_mOrd.gt.0 then set up L(m,0).
      If ( max_mOrd.gt.0 ) Then
        Write(6,*) '*****************************************'
        Write(6,*) ' Hot initial state not implemented yet !'
        Write(6,*) '*****************************************'
        Write(6,*)
        Call Quit_OnUserError()
      End If
C!
C!---- Initialize U matrix.
c      Write(6,*) 'CGGt[FCVal] Initialize U matrix.'               ! CGGt
c      Call XFlush(6)                                              ! CGGt
      call dcopy_((max_nord+1)*(max_nOrd2+1),[0.0d0],0,U,1)
      U(0,0) = 1.0d0
CGGt -------------------------------------------------------------------
c      Do kOsc = 1,nOsc
c      Write(6,*) 'CGGt[FCVal] Work(ipd2..)=',Work(ipd2+kOsc-1)    ! CGGt
c      Call XFlush(6)                                              ! CGGt
c      End Do
CGGt -------------------------------------------------------------------
C!
C!---- If max_nOrd.gt.0 then set up U(n,0).
c      Write(6,*) 'CGGt[FCVal] max_nOrd.gt.0 then set up U(n,0).'  ! CGGt
c      Call XFlush(6)                                              ! CGGt
      Do kOsc = 1,nOsc
c      Write(6,*) 'CGGt[FCVal] nInc(0,kOsc)=',nInc(0,kOsc)         ! CGGt
c      Call XFlush(6)                                              ! CGGt
c      Write(6,*) 'CGGt[FCVal] Work(ipd2..)=',Work(ipd2+kOsc-1)    ! CGGt
c      Call XFlush(6)                                              ! CGGt
        U(nInc(0,kOsc),0) = Work(ipd2+kOsc-1)
      End Do
c      Write(6,*) 'CGGt[FCVal] max_nInc=',max_nInc                 ! CGGt
c      Write(6,*) 'CGGt[FCVal]  iMaxYes=',iMaxYes                  ! CGGt
c      Call XFlush(6)                                              ! CGGt
      max_nInc = min(max_nInc,iMaxYes)
      If ( max_nInc.gt.0 ) Then
        Do iOrd = 1, max_nInc
c      Write(6,*) '              iOrd=',iOrd                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
          kOsc_start = nOsc
          Do While (( nMat(iOrd,kOsc_start).eq.0 ).and.
     &          ( kOsc_start.gt.1 ))
            kOsc_start = kOsc_start-1
          End Do
          Do kOsc = kOsc_start,nOsc
            Do lOsc = 1,nOsc
c      Write(6,*) 'iOrd,kOsc_start,kOsc,lOsc==',iOrd,kOsc_start,kOsc,lOsc
c      Call XFlush(6)                                              ! CGGt
              If ( nMat(iOrd,lOsc).gt.0 ) Then

                U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+
     &            Work(ipsqr+nMat(iOrd,lOsc))*
     &            Work(ipA2B2T+kOsc+nOsc*(lOsc-1)-1)*
     &            U(nDec(iOrd,lOsc),0)
              End If
            End Do

            U(nInc(iOrd,kOsc),0) =
     &        (U(nInc(iOrd,kOsc),0)+Work(ipd2+kOsc-1)*
     &        U(iOrd,0))/
     &        Work(ipsqr+nMat(nInc(iOrd,kOsc),kOsc))
          End Do
        End Do
      End If
C!
C! ---- Use recursion formula to obtain the rest of U.
c      Write(6,*) '            Use recursion ... rest of U.'       ! CGGt
c      Write(6,*) '            max_nOrd2=',max_nOrd2               ! CGGt
c      Call XFlush(6)                                              ! CGGt
      Do jOrd = 1,max_nOrd2
        lOsc = nOsc
        Do While (( nMat(jOrd,lOsc).eq.0 ).and.( lOsc.gt.1 ))
          lOsc = lOsc-1
        End Do
        Do iOrd = 0,max_nOrd
          Do kOsc = 1,nOsc
            If ( nMat(iOrd,kOsc).gt.0 ) Then
c      Write(6,*) '              ',iOrd,kOsc,nMat(iOrd,kOsc) ! CGGt
              U(iOrd,jOrd) = U(iOrd,jOrd)+
     &          (Work(ipsqr+nMat(iOrd,kOsc))/
     &           Work(ipsqr+nMat(jOrd,lOsc)))*
     &           Work(ipA2+kOsc+nOsc*(lOsc-1)-1)*
     &           U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
            End If
          End Do
        End Do
      End Do
C!
      Call GetMem('A2B2T','Free','Real',ipA2B2T,nOscSqr)
      Call GetMem('A1B1T','Free','Real',ipA1B1T,nOscSqr)
      Call GetMem('d2','Free','Real',ipd2,nOsc)
      Call GetMem('d1','Free','Real',ipd1,nOsc)
      Call GetMem('B2','Free','Real',ipB2,nOscSqr)
      Call GetMem('A2','Free','Real',ipA2,nOscSqr)
      Call GetMem('B1','Free','Real',ipB1,nOscSqr)
      Call GetMem('A1','Free','Real',ipA1,nOscSqr)
C!
C!---- Calculate Franck-Condon factors.
      If (iPrint.GE.3) then
        Write(6,*)' Franck-Condon factors for States in the Window: '
        Write(6,'(a,36a)') '  ',('=',i=1,36)
        Write(6,*)'     #     jOrd   FC factor     jSum '
        Write(6,'(a,36a)') '  ',('-',i=1,36)
      EndIf
      Do ii = 1, nYes
        jOrd = VibWind2(ii)
        dFC = FC00*L(0,0)*U(jOrd,0)
        FCWind2(ii) = dFC
        If (iPrint.GE.3) then
          loc_n_max = 0
          Do j=1,nOsc
            loc_n_max = loc_n_max + nMat(jOrd,j)
          EndDo
          Write(6,'(a2,i5,i9,e15.6,a2,i4)')
     &      ' ',ii,jOrd,FCWind2(ii),' ',loc_n_max
        EndIf
      EndDo
      If (iPrint.GE.3) then
        Write(6,'(a,36a)') '  ',('-',i=1,36)
        Write(6,*) ' FC_00 =',FC00
        Write(6,*)
      EndIf
C!
      If (iPrint.GE.4) then
        Write(6,*)
        Write(6,*) ' Full Franck-Condon factors (FC_00=',FC00,'):'
        Write(6,*) ' =================================================='
        Write(6,*) '    jOrd   FC            level                     '
        Write(6,*) ' --------------------------------------------------'
        Do jOrd = 0,max_nOrd
          loc_n_max = 0
          Do j=1,nOsc
            loc_n_max = loc_n_max + nMat(jOrd,j)
          EndDo
          dFC = FC00*L(0,0)*U(jOrd,0)
          Write(6,'(a,i8,e15.6,a2,i4,a2,24i3)')
     &    ' ',jOrd,dFC,' ',loc_n_max,' ',(nMat(jOrd,j),j=1,nOsc)
        End Do
        Write(6,*) ' --------------------------------------------------'
      EndIf
CGGt -------------------------------------------------------------------
C!
      Call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
      Call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
      Call GetMem('Alpha','Free','Real',ipAlpha,nOscSqr)
      Call GetMem('sqr','Free','Real',ipsqr,n+1)
      Call GetMem('temp2','Free','Real',iptemp2,nOscSqr)
      Call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
      Call GetMem('temp','Free','Real',iptemp,nOscSqr)
C!
c      Write(6,*) 'CGGt[FCVal] Exit'
c      Call XFlush(6)                                              ! CGGt
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(max_mInc)
         Call Unused_integer_array(mMat)
         Call Unused_integer_array(mInc)
         Call Unused_integer_array(mDec)
         Call Unused_real_array(r00)
         Call Unused_integer(nnsiz)
         Call Unused_integer(iMx_nOrd)
      End If
      End
