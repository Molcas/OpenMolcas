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
      FUNCTION qRINT (N,A,C,EXPA)
      Implicit Real*8(A-H,O-Z)
      Real*8 qRint
#include "real.fh"
#include "welcom.fh"
      Real*8 F(kmax+1)
*
      Call qEnter('qRINT')
      qRINT=Zero
      NN=N/2+1
*     Write (*,*) ' N,NN=',n,nn
      BP=Half*C
      START=Sqrt(Pi)
      PRSUM=EXP(BP*BP*A+EXPA)
      SUM=Zero
      ALF=sqrt(A)
      ARG=(BP*ALF)**2
      nT = 1
      call dcopy_(kMax+1,[Zero],0,F,1)
      Call Auxil([Arg],nT,F,nn-1)
*     Call RecPrt(' In qRint:Fm',' ',F,nt,nn)
      GINT=Zero
      Dac = -One
      Do 10 I=0,N
         TAL=(-BP)**(N-I)*Binom(n,i)
         J=(I/2)
*
         IF(J*2.EQ.I) THEN
            Dac = Dac * (Two*DBLE(J)-One)/Two
            FACT=ALF**(-I-1)*START*DAC
            FACT2=BP**(I+1)*F(J+1)
            GINT=(FACT-FACT2)*TAL+GINT
         ELSE
            GAL=One
            HINT=Zero
            Do 101 K=I-1,0,-2
               HINT=HINT+Half/A*BP**K*EXP(-ARG)*GAL
101            GAL=Half*DBLE(K)/A*GAL
            GINT=GINT+TAL*HINT
         EndIF
*
         qRINT=GINT*PRSUM
 10   Continue
*
      Call qExit('qRINT')
      Return
      End
