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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_getIntegralsSelectedByPSD(PrintLevel,
     &                                         Mode,tau,Tol,
     &                                         AB,CD,
     &                                         l_xInt,xInt,
     &                                         IntegralID)
C
C     Thomas Bondo Pedersen, January 2012.
C
C     Purpose: return (AB|CD) integrals in xInt.
C              If the LDF integrals are positive semidefinite:
C                 IntegralID='LDF  '
C                 return LDF integrals
C              else
C                 IntegralID='exact'
C                 return exact integrals
C
C              Positivity is established by means of Cholesky
C              decomposition and Tol defines the tolerance for
C              negative diagonal elements (and hence negativity).
C
C              LDF integrals are computed according to Mode,
C              which is simply passed on to the LDF integral routines
C              along with tau, the integral prescreening threshold.
C
C              Verbosity:
C                PrintLevel<=0: silent
C                PrintLevel=1:  only print stuff in case of error
C                PrintLevel>=2: verbose
C
      Implicit None
      Integer PrintLevel
      Integer Mode
      Real*8  tau
      Real*8  Tol
      Integer AB, CD
      Integer l_xInt
      Real*8  xInt(l_xInt)
      Character*5 IntegralID
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_getIntegralsSelectedByPSD')

      Real*8 Dummy
      Integer iDummy
      Parameter (iDummy=-1234567)

      Logical doDiagonalization
      Parameter (doDiagonalization=.False.)

      Logical Debug
#if defined (_DEBUG_)
      Parameter (Debug=.True.)
#else
      Parameter (Debug=.true.)
#endif

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical  isSymmetric
      External isSymmetric

      real*8 ddot_
      external ddot_

      Integer iPrint
      Integer nAB, nCD, nTot
      Integer ip_LDFInt, l_LDFInt
      Integer ip_Scr, l_Scr
      Integer l, m, ip1, ip2
      Integer ip_EigVal, l_Eigval
      Integer ip_EigVec, l_Eigvec
      Integer ip_Aux, l_Aux
      Integer nFound, iErr

      Real*8 ThrNeg, ThrFail
      Real*8 x
      Real*8 Norm_2(3)

      Integer i, j
      Integer AP_Atoms
      Integer iTri
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      iPrint=PrintLevel
      ThrNeg=-min(1.0d-13,abs(Tol)*1.0d-2)
      If (sign(1.0d0,Tol).lt.0.0d0) Then
         ThrFail=Tol
      Else
         ThrFail=-Tol
      End If

      ! Check dimension
      nAB=LDF_nBas_Atom(AP_Atoms(1,AB))*LDF_nBas_Atom(AP_Atoms(2,AB))
      nCD=LDF_nBas_Atom(AP_Atoms(1,CD))*LDF_nBas_Atom(AP_Atoms(2,CD))
      If (nAB*nCD .gt. l_xInt) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Print header
      If (iPrint.gt.0) Then
         Write(6,'(A,A,I6,A,I6)') SecNam,': AB=',AB,' CD=',CD
         Write(6,'(80A)') ('=',l=1,len(SecNam)+21)
         Write(6,'(3X,A,I2)') 'Mode=',Mode
         Write(6,'(3X,A,1P,D12.4)') 'Prescreening=',tau
         If (.not.doDiagonalization) Then
            Write(6,'(3X,A,1P,D12.4)') 'Negative diagonal zeroing=',
     &                                 ThrNeg
         End If
         Write(6,'(3X,A,1P,D12.4)') 'Tolerance for negative diagonal=',
     &                              ThrFail
      End If

      ! Compute LDF integrals, full AB,CD block
      ! [i.e. including (AB|AB) and (CD|CD)]
      ! Save (AB|CD) in xInt; it may be returned later
      If (AB.eq.CD) Then
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)') 'Computing (AB|AB) LDF integrals'
         End If
         nTot=nAB
         l_LDFInt=nTot**2
         Call GetMem('LDFInt','Allo','Real',ip_LDFInt,l_LDFInt)
         Call LDF_ComputeValenceIntegralsFromC(Mode,tau,
     &                                         AB,AB,l_LDFInt,
     &                                         Work(ip_LDFInt))
         Call dCopy_(l_LDFInt,Work(ip_LDFInt),1,xInt(1),1) !save (AB|AB)
         If (Debug) Then
            Norm_2(1)=dDot_(l_LDFInt,Work(ip_LDFInt),1,
     &                              Work(ip_LDFInt),1)
            Norm_2(2)=0.0d0
            Norm_2(3)=0.0d0
         End If
      Else
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)') 'Computing (AB|AB) LDF integrals'
         End If
         nTot=nAB+nCD
         l_LDFInt=nTot**2
         Call GetMem('LDFInt','Allo','Real',ip_LDFInt,l_LDFInt)
         l_Scr=max(nAB,nCD)**2
         Call GetMem('IntScr','Allo','Real',ip_Scr,l_Scr)
         Call LDF_ComputeValenceIntegralsFromC(Mode,tau,
     &                                         AB,AB,l_Scr,
     &                                         Work(ip_Scr))
         If (Debug) Then
            Norm_2(1)=dDot_(nAB**2,Work(ip_Scr),1,Work(ip_Scr),1)
         End If
         ip1=ip_Scr
         ip2=ip_LDFInt
         Do l=1,nAB
            Call dCopy_(nAB,Work(ip1),1,Work(ip2),1)
            ip1=ip1+nAB
            ip2=ip2+nTot
         End Do
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)') 'Computing (AB|CD) LDF integrals'
         End If
         Call LDF_ComputeValenceIntegralsFromC(Mode,tau,
     &                                         AB,CD,l_Scr,
     &                                         Work(ip_Scr))
         If (Debug) Then
            Norm_2(2)=dDot_(nAB*nCD,Work(ip_Scr),1,Work(ip_Scr),1)
         End If
         Call dCopy_(nAB*nCD,Work(ip_Scr),1,xInt(1),1) !save (AB|CD)
         ip1=ip_Scr
         ip2=ip_LDFInt+nTot*nAB
         Do l=1,nCD
            Call dCopy_(nAB,Work(ip1),1,Work(ip2),1)
            ip1=ip1+nAB
            ip2=ip2+nTot
         End Do
         ip1=ip_Scr
         ip2=ip_LDFInt+nAB
         Do l=1,nAB
            Call dCopy_(nCD,Work(ip1),nAB,Work(ip2),1)
            ip1=ip1+1
            ip2=ip2+nTot
         End Do
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)') 'Computing (CD|CD) LDF integrals'
         End If
         Call LDF_ComputeValenceIntegralsFromC(Mode,tau,
     &                                         CD,CD,l_Scr,
     &                                         Work(ip_Scr))
         If (Debug) Then
            Norm_2(3)=dDot_(nCD**2,Work(ip_Scr),1,Work(ip_Scr),1)
         End If
         ip1=ip_Scr
         ip2=ip_LDFInt+nTot*nAB+nAB
         Do l=1,nCD
            Call dCopy_(nCD,Work(ip1),1,Work(ip2),1)
            ip1=ip1+nCD
            ip2=ip2+nTot
         End Do
         Call GetMem('IntScr','Free','Real',ip_Scr,l_Scr)
      End If

      If (Debug) Then
         If (isSymmetric(Work(ip_LDFInt),nTot,1.0d-12)) Then
            If (iPrint.gt.1) Then
               If (AB.eq.CD) Then
                  Write(6,'(3X,A)') '(AB|AB)  is symmetric'
               Else
                  Write(6,'(3X,A)') '(AB|AB) | (AB|CD)'
                  Write(6,'(3X,A)') '-----------------  is symmetric'
                  Write(6,'(3X,A)') '(CD|AB) | (CD|CD)'
               End If
            End If
         Else
            If (iPrint.gt.1) Then
               Call Cho_Head('Block: (AB|AB)','-',80,6)
               Call Cho_Output(Work(ip_LDFInt),1,nAB,1,nAB,nTot,nTot,
     &                         1,6)
               If (AB.ne.CD) Then
                  Call Cho_Head('Block: (AB|CD)','-',80,6)
                  Call Cho_Output(Work(ip_LDFInt),1,nAB,nAB+1,nTot,
     &                            nTot,nTot,1,6)
                  Call Cho_Head('Block: (CD|AB)','-',80,6)
                  Call Cho_Output(Work(ip_LDFInt),nAB+1,nTot,1,nAB,
     &                            nTot,nTot,1,6)
                  Call Cho_Head('Block: (CD|CD)','-',80,6)
                  Call Cho_Output(Work(ip_LDFInt),nAB+1,nTot,nAB+1,nTot,
     &                            nTot,nTot,1,6)
               End If
            End If
            Call WarningMessage(2,
     &                      SecNam//': integral block is not symmetric')
            Call LDF_Quit(1)
         End If
         iErr=0
         x=0.0d0
         Do i=1,nTot
            If (sign(1.0d0,Work(ip_LDFInt-1+nTot*(i-1)+i)).lt.0.0d0)
     &      Then
               If (Work(ip_LDFInt-1+nTot*(i-1)+i).lt.ThrNeg) Then
                  iErr=iErr+1
                  If (iPrint.gt.0) Then
                     Write(6,'(3X,A,1P,D20.10)')
     &               'Negative diagonal element:',
     &               Work(ip_LDFInt-1+nTot*(i-1)+i)
                  End If
                  If (Work(ip_LDFInt-1+nTot*(i-1)+i).lt.x) Then
                     x=Work(ip_LDFInt-1+nTot*(i-1)+i)
                  End If
               End If
            End If
         End Do
         If (iErr.eq.0) Then
            If (iPrint.gt.1) Then
               Write(6,'(3X,A)') 'Diagonal positive'
            End If
         Else
            If (iPrint.gt.0) Then
               Write(6,'(3X,A,1P,D20.10)')
     &         'Diagonal NOT positive, smallest=',x
            End If
            If (x.lt.ThrFail) Then
               Call WarningMessage(2,SecNam//': negative diagonal')
               Call LDF_Quit(1)
            End If
         End If
         x=0.0d0
         Do j=1,nAB
            x=x+dDot_(nAB,Work(ip_LDFInt+nTot*(j-1)),1,
     &                   Work(ip_LDFInt+nTot*(j-1)),1)
         End Do
         If (abs(x-Norm_2(1)).gt.1.0d-10) Then
            Call WarningMessage(2,SecNam//': (AB|AB) norm error')
            Call LDF_Quit(1)
         End If
         If (CD.ne.AB) Then
            x=0.0d0
            Do j=1,nCD
               x=x+dDot_(nAB,Work(ip_LDFInt+nTot*(nAB+j-1)),1,
     &                      Work(ip_LDFInt+nTot*(nAB+j-1)),1)
            End Do
            If (abs(x-Norm_2(2)).gt.1.0d-10) Then
               Call WarningMessage(2,SecNam//': (AB|CD) norm error')
               Call LDF_Quit(1)
            End If
            x=0.0d0
            Do j=1,nAB
               x=x+dDot_(nCD,Work(ip_LDFInt+nTot*(j-1)+nAB),1,
     &                      Work(ip_LDFInt+nTot*(j-1)+nAB),1)
            End Do
            If (abs(x-Norm_2(2)).gt.1.0d-10) Then
               Call WarningMessage(2,SecNam//': (CD|AB) norm error')
               Call LDF_Quit(1)
            End If
            x=0.0d0
            Do j=1,nCD
               x=x+dDot_(nCD,Work(ip_LDFInt+nTot*(nAB+j-1)+nAB),1,
     &                      Work(ip_LDFInt+nTot*(nAB+j-1)+nAB),1)
            End Do
            If (abs(x-Norm_2(3)).gt.1.0d-10) Then
               Call WarningMessage(2,SecNam//': (CD|CD) norm error')
               Call LDF_Quit(1)
            End If
         End If
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)')
     &      'Block norms OK within 1.0d-10'
         End If
         x=dDot_(nTot**2,Work(ip_LDFInt),1,Work(ip_LDFInt),1)
         If (abs(x-(Norm_2(1)+2.0d0*Norm_2(2)+Norm_2(3))).gt.1.0d-10)
     &   Then
            Call WarningMessage(2,SecNam//': integral norm error')
            Call LDF_Quit(1)
         End If
      End If

      If (doDiagonalization) Then
         ! Check for positivity by diagonalization
         If (iPrint.gt.1) Then
            If (AB.eq.CD) Then
               Write(6,'(3X,A)') 'Triangularization of (AB|AB)'
            Else
               Write(6,'(3X,A)')
     &                        '                     (AB|AB) | (AB|CD)'
               Write(6,'(3X,A)')
     &                        'Triangularization of -----------------'
               Write(6,'(3X,A)')
     &                        '                     (CD|AB) | (CD|CD)'
            End If
         End If
         l_Aux=nTot**2
         l_EigVal=nTot
         l_EigVec=nTot**2
         Call GetMem('Aux','Allo','Real',ip_Aux,l_Aux)
         Call GetMem('EigVal','Allo','Real',ip_EigVal,l_EigVal)
         Call GetMem('EigVec','Allo','Real',ip_EigVec,l_EigVec)
        If (Debug) Call dCopy_(nTot**2,Work(ip_LDFInt),1,Work(ip_Aux),1)
         Call Local_Triang(nTot,Work(ip_LDFInt)) !triangularization
         If (Debug) Then
            iErr=0
            Do m=1,nTot
               Do l=1,nTot
                  If (abs(Work(ip_LDFInt-1+iTri(l,m))
     &                   -Work(ip_Aux-1+nTot*(m-1)+l)).gt.1.0d-12) Then
                     iErr=iErr+1
                     If (iPrint.gt.1) Then
                        Write(6,'(3X,A,I12,A,2I7,1P,2(1X,D15.6))')
     &                  'Triangularization error',iErr,': i,j,sq,tri=',
     &                  l,m,Work(ip_Aux-1+nTot*(m-1)+l),
     &                  Work(ip_Aux-1+iTri(l,m))
                     End If
                  End If
               End Do
            End Do
            If (iErr.ne.0) Then
               If (iPrint.gt.0) Then
                  Write(6,'(3X,A,I12)')
     &            'Triangularization errors:',iErr
               End If
               Call WarningMessage(2,
     &             SecNam//': integral matrix triangularization failed')
               Call LDF_Quit(1)
            Else
               If (iPrint.gt.0) Then
                  Write(6,'(3X,A)')
     &            'Triangularization succeeded'
               End If
            End If
         End If
         If (iPrint.gt.1) Then
            If (AB.eq.CD) Then
               Write(6,'(3X,A)') 'Diagonalizing (AB|AB)'
            Else
               Write(6,'(3X,A)') '              (AB|AB) | (AB|CD)'
               Write(6,'(3X,A)') 'Diagonalizing -----------------'
               Write(6,'(3X,A)') '              (CD|AB) | (CD|CD)'
            End If
         End If
         nFound=0
         iErr=0
         Call Diag_Driver('V','A','L',nTot,Work(ip_LDFInt),Work(ip_Aux),
     &                    nTot,Dummy,Dummy,iDummy,iDummy,
     &                    Work(ip_EigVal),Work(ip_EigVec),nTot,1,
     &                    1,'A',nFound,iErr)
         iErr=0
         l=0
         Do While (l.lt.nTot)
            If (sign(1.0d0,Work(ip_EigVal+l)).lt.0.0d0) Then
               If (Work(ip_EigVal+l).lt.ThrNeg .and. iPrint.gt.1) Then
                  Write(6,'(2X,A,1P,D20.10)')
     &            'Negative eigenvalue=',Work(ip_EigVal+l)
               End If
               If (Work(ip_EigVal+l).lt.ThrFail) Then
                  If (iPrint.gt.0) Then
                     Write(6,'(3X,A,1P,D20.10)')
     &               'NON-PSD! Negative eigenvalue=',Work(ip_EigVal+l)
                  End If
                  iErr=iErr+1
                  l=l+1
               Else
                  l=nTot
               End If
            Else
               l=nTot ! break loop
            End If
         End Do
         Call GetMem('EigVec','Free','Real',ip_Eigvec,l_EigVec)
         Call GetMem('EigVal','Free','Real',ip_EigVal,l_EigVal)
         Call GetMem('Aux','Free','Real',ip_Aux,l_Aux)
      Else
         ! Check for positivity by Cholesky decomposition
         If (iPrint.gt.1) Then
            If (AB.eq.CD) Then
               Write(6,'(3X,A)') 'Cholesky dec. (AB|AB)'
            Else
               Write(6,'(3X,A)') '              (AB|AB) | (AB|CD)'
               Write(6,'(3X,A)') 'Cholesky dec. -----------------'
               Write(6,'(3X,A)') '              (CD|AB) | (CD|CD)'
            End If
         End If
         l_Aux=nTot**2
         Call GetMem('CDVec','Allo','Real',ip_Aux,l_Aux)
         Call CD_InCore_1(Work(ip_LDFInt),nTot,Work(ip_Aux),nTot,l,
     &                    1.0d-12,ThrNeg,ThrFail,iErr)
         If (iErr.ne.0) Then
            If (iPrint.gt.0) Then
               Write(6,'(3X,A,I6)')
     &         'Cholesky decomposition failed with error code',iErr
            End If
            If (iErr.lt.0) Then
               Call WarningMessage(2,SecNam//': error in CD_InCore')
               Call LDF_Quit(1)
            End If
            If (iErr.eq.101 .and. iPrint.gt.1) Then
               x=Work(ip_LDFInt)
               Do i=2,nTot
                  x=min(x,Work(ip_LDFInt-1+nTot*(i-1)+i))
               End Do
               Write(6,'(3X,A,1P,D20.10)')
     &         'Smallest diagonal element=',x
            End If
         Else
            If (iPrint.gt.1) Then
               Write(6,'(3X,A,2I6)')
     &         'Cholesky decomposition succeeded: dimension and rank=',
     &         nTot,l
               If (l.lt.nTot) Then
                  Write(6,'(3X,A)')
     &            'Matrix is positive semidefinite'
               Else If (l.eq.nTot) Then
                  Write(6,'(3X,A)')
     &            'Matrix is positive definite'
               Else
                  Call WarningMessage(2,SecNam//': rank>dimension')
                  Call LDF_Quit(1)
               End If
            End If
         End If
         Call GetMem('CDVec','Free','Real',ip_Aux,l_Aux)
      End If
      Call GetMem('LDFInt','Free','Real',ip_LDFInt,l_LDFInt)

      ! If negative eigenvalues, compute and return conventional
      ! integrals, else return LDF integrals (already in xInt)
      If (iErr.ne.0) Then
         If (iPrint.gt.1) Then
            Write(6,'(3X,A)') 'Computing exact integrals (AB|CD)'
         End If
         IntegralID='exact'
         Call LDF_ComputeValenceIntegrals(AB,CD,l_xInt,xInt)
      Else
         IntegralID='LDF  '
      End If
      If (iPrint.gt.0) Then
         If (AB.eq.CD) Then
            Write(6,'(3X,A,A,A)')
     &      'Returning ',IntegralID,' integrals (AB|AB)'
         Else
            Write(6,'(3X,A,A,A)')
     &      'Returning ',IntegralID,' integrals (AB|CD)'
         End If
      End If

      End
