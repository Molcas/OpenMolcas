!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
        subroutine NoBlanks(out,n,in)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
        character out*(*), in*(n)
        integer flag
        flag=-1
        j=0
        do 1 i=1,len(in)
         if(flag.eq.-1.and.in(i:i).eq.' ') goto 1
         flag=0
         if(i.le.len(in)-1) then
           if(in(i:i+1).eq.'  ') goto 1
         endif
         j=j+1
         out(j:j)=in(i:i)
1       continue
        out(j+1:)=' '
        return
        end
