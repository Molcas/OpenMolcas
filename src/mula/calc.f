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
      Subroutine Calc_r00(C1,C2,W1,W2,C,W,alpha1,alpha2,r00,r01,r02,
     &                             det0,det1,det2,FC00,nOsc)
C!
C!  Purpose:
C!    Calculate geometry of the intermediate oscillator.
C!
      Implicit Real*8 ( a-h,o-z )
      Real*8 C1( nOsc,nOsc ), C2( nOsc,nOsc ), C( nOsc,nOsc )
      Real*8 W1( nOsc,nOsc ), W2( nOsc,nOsc ), W( nOsc,nOsc )
      Real*8 alpha1 ( nOsc,nOsc ),alpha2 ( nOsc,nOsc )
      Real*8 r00( nOsc ), r01( nOsc ), r02( nOsc )
#include "WrkSpc.fh"
C!
C!---- Initialize.
      my1=nOsc
      my2=nOsc
      nOscSqr = nOsc**2

      Call GetMem('temp','Allo','Real',iptemp,nOscSqr)
C!
C!---- Calculate alpha1, alpha2 and alpha.
      Call GetMem('alpha','Allo','Real',ipalpha,nOscSqr)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C1,nOsc,
     &            C1,nOsc,
     &            0.0d0,alpha1,nOsc)
      call dscal_(nOscSqr,0.5d0,alpha1,1)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C2,nOsc,
     &            C2,nOsc,
     &            0.0d0,alpha2,nOsc)
      call dscal_(nOscSqr,0.5d0,alpha2,1)

c       temp = alpha1+alpha2
      call dcopy_(nOscSqr,alpha1,1,Work(iptemp),1)
      Call Daxpy_(nOscSqr,1.0d0,alpha2,1,Work(iptemp),1)

c       alpha = 0.5d0*temp
      call dcopy_(nOscSqr,[0.0d0],0,Work(ipalpha),1)
      Call Daxpy_(nOscSqr,0.5d0,Work(iptemp),1,Work(ipalpha),1)
C!
C!---- Calculate C using a Cholesky factorization of 2*alpha.
      Call Cholesky(Work(iptemp),C,nOsc)

C!
C!---- Calculate W.
      call dcopy_(nOscSqr,[0.0d0],0,W,1)
      call dcopy_(nOsc,[1.0d0],0,W,nOsc+1)
      call dcopy_(nOscSqr,C,1,Work(iptemp),1)
c       temp = C
      Call Dool(Work(iptemp),nOsc,nOsc,W,my1,my2,det0)
      det0=abs(det0)
C!
C!---- Calculate r00.
      Call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
      Call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)
      Call GetMem('r_temp','Allo','Real',ipr_temp,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,alpha1,nOsc,
     &            r01,nOsc,
     &            0.0d0,Work(ipr_temp1),nOsc)
      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,alpha2,nOsc,
     &            r02,nOsc,
     &            0.0d0,Work(ipr_temp2),nOsc)
c       r_temp(:,1) = r_temp1+r_temp2
      do i=1,nOsc
      Work(ipr_temp+i-1)=Work(ipr_temp1+i-1)+Work(ipr_temp2+i-1)
      enddo
c       temp = 2.0*alpha
      call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
      Call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

      Call Dool(Work(iptemp),nOsc,nOsc,Work(ipr_temp),nOsc,1,det)
      Call GetMem('beta','Allo','Real',ipbeta,nOscSqr)
c       r00 = r_temp(:,1)
      do i=1,nOsc
      r00(i)=Work(ipr_temp+i-1)
      enddo
C!
C!---- Calculate beta.
      Call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
c       temp1 = alpha1
      call dcopy_(nOscSqr,alpha1,1,Work(iptemp1),1)

c       temp  = 2.0d0*alpha
      call dcopy_(nOscSqr,[0.0d0],0,Work(iptemp),1)
      Call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

      Call Dool(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,alpha2,nOsc,
     &            Work(iptemp1),nOsc,
     &            0.0d0,Work(ipbeta),nOsc)
C!
C!---- Calculate FC00.
c       r_temp1 = r01-r02
      call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
      Call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,Work(ipbeta),nOsc,
     &            Work(ipr_temp1),nOsc,
     &            0.0d0,Work(ipr_temp2),nOsc)
      FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
      FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)

      Call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
      Call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
      Call GetMem('r_temp','Free','Real',ipr_temp,nOsc)
      Call GetMem('alpha','Free','Real',ipalpha,nOscSqr)
      Call GetMem('beta','Free','Real',ipbeta,nOscSqr)
      Call GetMem('temp','Free','Real',iptemp,nOscSqr)
      Call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(W1)
         Call Unused_real_array(W2)
      End If
      End
C!
C!-----------------------------------------------------------------------!
