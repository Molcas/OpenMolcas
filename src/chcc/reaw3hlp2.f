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
        subroutine ReaW3hlp2 (Ww,Wx,dima,dimb,no,LunName,LunAux)
c
c        this routine do:
c        reconstruct  Ww(a",be',b,i)  for aSGrp=beSGrp
c        from (a">=be"|b,i) records in V3 file LunName
c
        implicit none
        integer dima,dimb,no,LunAux
        character*8 LunName
        real*8 Ww(1:dima,1:dima,1:dimb,1:no)
        real*8 Wx(*)

c        help variables
c
        integer i,a,be,abebi,b,length

c        read block (a">=be"|b"_i)
        length=(no*dima*(dima+1)*dimb)/2
*       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
        call rea_chcc (LunAux,length,Wx(1))
cmp        call mv0zero (length,length,Wx(1))
        close (LunAux)
c
c          Expand and Set Ww(a",be",b",i) <- Wx(a">=be"|b",i)
        abebi=0
c
        do i=1,no
          do b=1,dimb
          do a=1,dima
          do be=1,a
            abebi=abebi+1
            Ww(a,be,b,i)=Wx(abebi)
            Ww(be,a,b,i)=Wx(abebi)
          end do
          end do
          end do
        end do
c
c
        return
        end
