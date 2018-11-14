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
      Subroutine decomp(A,Jiso,Jsym,Jantisym,dbg)
      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Real(kind=wp), intent(in) :: A(3,3)
      Real(kind=wp), intent(out):: Jiso, Jsym(3,3), Jantisym(3,3)
      Real(kind=wp), external   :: real_1_trace2
      Logical, intent(in)       :: dbg

      Integer :: i, j
      Real(kind=wp) :: tmp
      Real(kind=wp) :: Dtmp(3,3)

      tmp=0.0_wp
      Jiso=0.0_wp
      Jsym=0.0_wp
      Jantisym=0.0_wp
      !-------------------------------------
      do i=1,3
        tmp=tmp+A(i,i)
      end do
      Jiso=tmp/3.0_wp
      !-------------------------------------
      Do i=1,3
        Jsym(i,i)=A(i,i)-Jiso
      End Do

      ! find the symmetric matrix:
      Do i=1,3
        Do j=1,3
          if(i==j) cycle
          Jsym(i,j)=(A(i,j)+A(j,i))/2.0_wp
        End Do
      End Do
      ! find the anti-symmetric matrix:
      Do i=1,3
        Do j=1,3
          if(i==j) cycle
          Jantisym(i,j)=(A(i,j)-A(j,i))/2.0_wp
        End Do
      End Do

      If (dbg) Then
         Dtmp=0.0_wp
         Do i=1,3
           Dtmp(i,i)=Jiso+Jsym(i,i)+Jantisym(i,i)
           Do j=1,3
             if(i==j) cycle
               Dtmp(i,j)=Jsym(i,j)+Jantisym(i,j)
           End Do
         End Do

         Write(6,*)
         Write(6,*) 'J recovered = '
         Do i=1,3
            Write(6,'(3F24.14)') (Dtmp(i,j), j=1,3)
         End Do
      End If

      End subroutine decomp
