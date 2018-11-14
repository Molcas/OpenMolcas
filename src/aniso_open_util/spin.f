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
      Complex*16 function spin(l,dim,m1,m2)
!  Returns the value of the < m1 | S_l | m2 >
!

      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in):: m1,m2,dim,l
      Integer            :: ipar
      Real(kind=wp)      :: S, MM1, MM2, R, D, F
      Logical            :: dbg

      dbg=.false.

      spin=0.0_wp
      ipar=mod(dim,2)
      S=dble(dim-1)/2.0_wp
      R=0.0_wp

      If(dbg) Write(6,'(A,4I3)') 'l,dim,m1,m2=', l, dim, m1, m2
      If(dbg) Write(6,'(A,I3,F8.3)') ' ipar,  S  =', ipar, S

c set up the true values of the MM1 and MM2
c i.e. the projections of the spin momentum S
      If ( ipar==0 ) Then
         If ( m1<0 ) Then
            MM1=dble(m1) + 0.5_wp
         Else
            MM1=dble(m1) - 0.5_wp
         End If

         If ( m2<0 ) Then
            MM2=dble(m2) + 0.5_wp
         Else
            MM2=dble(m2) - 0.5_wp
         End If
      Else
            MM1=dble(m1)
            MM2=dble(m2)
      End If
      If(dbg) Write(6,'(A,2F8.3)') 'mm1,mm2=',mm1, mm2
      Call xFlush(6)
c compute the value of the matrix element, for each possible cartesian component.
c l=1  => X
c l=2  => Y
c l=3  => Z

      If ( l==1 ) Then

         If (      (MM1-1.0_wp)==MM2 ) Then
            D=S+MM1
            F=S-MM1+1.0_wp
            R=0.5_wp*sqrt(D*F)

            spin=cmplx( R, 0.0_wp, wp)

         Else If ( (MM1+1.0_wp)==MM2 ) Then

            D=S-MM1
            F=S+MM1+1.0_wp
            R=0.5_wp*sqrt(D*F)
            spin=cmplx( R, 0.0_wp, wp)

         Else
            spin=cmplx(0.0_wp,0.0_wp,wp)
         End If
         If(dbg) Write(6,'(A,2F8.3)') 'X, spin =',spin
         Call xFlush(6)

      Else If ( l==2 ) Then

         If (     (MM1-1.0_wp)==MM2 ) Then

            D=S+MM1
            F=S-MM1+1.0_wp
            R=-0.5_wp*sqrt(D*F)
            spin=cmplx(0.0_wp, R, wp)

         Else If( (MM1+1.0_wp)==MM2 ) Then

            D=S-MM1
            F=S+MM1+1.0_wp
            R=0.5_wp*sqrt(D*F)
            spin=cmplx(0.0_wp, R, wp)

         Else
            spin=cmplx(0.0_wp,0.0_wp,wp)
         End If
         If(dbg) Write(6,'(A,2F8.3)') 'Y, spin =',spin
         Call xFlush(6)

      Else If ( l==3 ) Then

         If ( MM1 .ne. MM2 ) Then
            spin=cmplx(0.0_wp,0.0_wp,wp)
         Else
            spin=cmplx(MM1,0.0_wp,wp)
         End If
         If(dbg) Write(6,'(A,2F8.3)') 'Z, spin =',spin
         Call xFlush(6)

      Else
         Write(6,'(A)') 'The spin function gives a wrong number'
         Return
      End If

      If(dbg) Write(6,*) 'Upon Return:  spin =',spin
      Call xFlush(6)
      Return
      End function spin
