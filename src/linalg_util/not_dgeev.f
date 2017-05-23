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
      Subroutine not_DGeEV(iOpt,a,lda,w,z,ldz,select,n,aux,naux)
      Implicit Real*8 (a-h,o-z)
      Real*8 a(lda,n),w(2,n), z(2*ldz*n), aux(naux)
      parameter(nw1=200)
      Real*8 w1(nw1)
      Logical select(n)
*
      If (iOpt.eq.2) Then
         Write (6,*) 'not_DGeEV: iOpt=2 is not implemented yet!'
         Call Abend()
      End If
      If (ldz.ne.n) Then
         Write (6,*) 'not_DGeEV: ldz=/=n is not implemented yet!'
         Call Abend()
      End If
      If (iOpt.eq.0) Then
         Write (6,*) 'not_DGeEV: iOpt=0 is not implemented yet!'
         Call Abend()
      End If
      iOff = n+1
      If (nAux.lt.2*n) Then
         Write (6,*) 'not_DGeEV: nAux is too small (naux<2*n)!'
         Call Abend()
      End If
      If (n.gt.nw1) Then
         Write (6,*) 'not_DGeEV: nw1 is too small (nw1<n)!'
         Call Abend()
      End If
      iErr=0
      Call XEIGEN(iOpt,lda,n,a,w,w1,z,aux,aux(iOff),iErr)
      If (iErr.ne.0) Then
         Write (6,*) ' not_DGeEV: iErr=/= 0!'
         Call Abend()
      End If
*     Call RecPrt('w',' ',w,n,1)
*     Call RecPrt('w1',' ',w1,n,1)
*     Call RecPrt('z',' ',z,n,n)
*
*-----Order eigenvalues to ESSL standard
*
      call dcopy_(n,w,1,aux,1)
      Do i = 1, n
         w(1,i)=aux(i)    ! Real
         w(2,i)=w1(i)     ! Imaginary
      End Do
*
*-----Order eigenvector to ESSL standard
*
*
      i = N
      Do While (i.ge.1)
        If ( w(2,i).ne.0.0d0 ) then
           i=i-1
           iOff=(i-1)*n
           call dcopy_(2*N,Z(iOff+1),1,Aux,1)
           iOff=(i-1)*2*n
           call dcopy_(N,Aux(1  ),1,Z(iOff+1),2)
           call dcopy_(N,Aux(1+N),1,Z(iOff+2),2)
           iOff=i*2*n
           call dcopy_(N,Aux(1  ),1,Z(iOff+1),2)
           call dcopy_(N,Aux(1+N),1,Z(iOff+2),2)
           Call DScal_(N,-1.0D0,Z(iOff+2),2)
        Else
           iOff=(i-1)*n
           call dcopy_(N,Z(iOff+1),1,Aux,1)
           iOff=(i-1)*2*n
           call dcopy_(N,Aux(1  ),1,Z(iOff+1),2)
           call dcopy_(N,0.0D0,   0,Z(iOff+2),2)
        End If
        i=i-1
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical_array(select)
      End
