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
* Copyright (C) 2013, Victor P. Vysotskiy                              *
************************************************************************
      Subroutine NIdiag_New(H,U,n,nv,iOpt)
************************************************************************
*                                                                      *
* This routine is a wrapper that calls appropriate LAPACK routines to  *
* perform diagonalization of symmetric matrices stored in lower        *
* triangular form.                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Victor P. Vysotskiy                                         *
*          Lund university, Sweden                                     *
* Written  2013                                                        *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
* n    - Dimension of matrix                                           *
* nv   - Length of eigenvectors nv>=n                                  *
* H    - Matrix to be diagonalized                                     *
* U    - Eigenvectors                                                  *
* iOpt - Option flag, for future improvements.                         *
*----------------------------------------------------------------------*
#include "WrkSpc.fh"
      External OrbPhase
      Integer n,nv,iOpt
      Real*8  H(*),U(nv,n)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer ipDIA,ipOFF,ipTAU,ipIWRK,ipWork,ipEVL,ipIPSZ,ipHDUP
      Integer lrwrk,liwrk,lh,info,I,M
      Real*8  abstol,dlamch_, Tmp, OrbPhase
      External dlamch_
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      If (n.eq.0) Return
#ifdef _DEBUGPRINT_
      Write(6,*) "New nidiag"
#endif
      Call FZero(U,nv*n)

      lh=n*(n+1)/2
      liwrk=10*n
      lrwrk=20*n

      Call GetMem('DIA','ALLO','REAL',ipDIA,n)
      Call GetMem('EVL','ALLO','REAL',ipEVL,n)
      Call GetMem('OFF','ALLO','REAL',ipOFF,(n-1))
      Call GetMem('TAU','ALLO','REAL',ipTAU,(n-1))
      Call GetMem('IPSZ','ALLO','INTE',ipIPSZ,2*n)
      Call GetMem('IWRK','ALLO','INTE',ipIWRK,liwrk)
      Call GetMem('RWRK','ALLO','REAL',ipWORK,lrwrk)
      Call GetMem('HDUP','ALLO','REAL',ipHDUP,lh)

      call dcopy_(lh,H,1,Work(ipHDUP),1)

      info=0
      call dsptrd_('U',n,Work(ipHDUP),Work(ipDIA),Work(ipOFF),
     &Work(ipTAU),info)

      If(info.ne.0) Then
#ifdef _DEBUGPRINT_
         Write(6,'(A,I4)')"Failed to tridiagonalize matrix",
     &   info
#endif
         Go To 10
      End If
#if defined(_ACML_) && defined(__PGI)
      CALL ILAENVSET(10,'X','X',0,0,0,0,1,INFO)
      CALL ILAENVSET(11,'X','X',0,0,0,0,1,INFO)
#endif
      abstol=dlamch_('Safe minimum')
      info=0
      call dstevr_('V','A',n,Work(ipDIA),Work(ipOFF),
     &Work(ip_Dummy),Work(ip_Dummy),iWork(ip_iDummy),
     &iWork(ip_iDummy),abstol,M,Work(ipEVL),U,nv,
     &iWork(ipIPSZ),Work(ipWORK),lrwrk,iWork(ipIWRK),
     &liwrk,info)

      If(info.ne.0) Then
#ifdef _DEBUGPRINT_
         Write(6,'(A,I4)') "Failed to diagonalize matrix",
     &   info
#endif
         Go To 10
      End If

      call dopmtr_('Left','U','N',N,N,Work(ipHDUP),Work(ipTAU),
     &U,nv,Work(ipWork),info)

      If(info.ne.0) Then
#ifdef _DEBUGPRINT_
        Write(6,'(A,I4)') "Failed to back transform vectors",
     &  info
#endif
        Go To 10
      End If

      call dcopy_(lh,Work(ipHDUP),1,H,1)

      Do I=1,N
         H((I*(I+1))/2)=Work(ipEVL+I-1)
      End Do

10    Continue
      Call GetMem('DIA','FREE','REAL',ipDIA,n)
      Call GetMem('EVL','FREE','REAL',ipEVL,n)
      Call GetMem('OFF','FREE','REAL',ipOFF,(n-1))
      Call GetMem('TAU','FREE','REAL',ipTAU,(n-1))
      Call GetMem('IPSZ','FREE','INTE',ipIPSZ,2*n)
      Call GetMem('RWRK','FREE','REAL',ipWORK,lrwrk)
      Call GetMem('IWRK','FREE','INTE',ipIWRK,liwrk)
      Call GetMem('HDUP','FREE','REAL',ipHDUP,n*(n+1)/2)

      if(info.ne.0) Then
#ifdef _DEBUGPRINT_
         Write(6,'(A)')
     &   "Using the old Givens rot. based routine"
#endif
         Call NIdiag(H,U,n,nv,iOpt)
      End If
*
      Do i = 1, n
         Tmp = OrbPhase(U(1,i),nv)
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_real(Tmp)
#endif
      End
