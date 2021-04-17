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
        subroutine InsL (Ll,Lg,ncLoc,nc,ncOff,dim1)
c
c        this routine do:
c        insert Llocal(ml,dim1) into Lglobal(m,dim1)
c        on a corresponding place
c
        implicit none
        integer ncLoc,nc,ncOff,dim1
        real*8 Ll(1:ncLoc,1:dim1)
        real*8 Lg(1:nc,1:dim1)
c
c        help variables
        integer m,p
c
        do p=1,dim1
          do m=1,ncLoc
            Lg(ncOff+m,p)=Ll(m,p)
          end do
        end do
c
        return
        end
