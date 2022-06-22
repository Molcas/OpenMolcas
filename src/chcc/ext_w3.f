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
        subroutine Ext_W3 (V3,M2,nc,no,dimc,dimcpp,addcpp)
c
c        this routine do:
c         Extract M2(m,c",i) <- V3(m,c',i)
c
        implicit none
        integer nc,no,dimc,dimcpp,addcpp
        real*8 V3(1:nc,1:dimc,1:no)
        real*8 M2(1:nc,1:dimcpp,1:no)
c
c        help variables
        integer m,i,c,cpp
c
        do i=1,no
          c=addcpp
          do cpp=1,dimcpp
          c=c+1
            do m=1,nc
              M2(m,cpp,i)=V3(m,c,i)
            end do
          end do
        end do
c
        return
        end
