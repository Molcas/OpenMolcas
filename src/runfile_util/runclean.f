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
       subroutine init_run_use
#include "run_use_common.fh"
       do i=1,nTocCA
         i_run_CA_used(i)=0
       enddo
       do i=1,nTocDA
         i_run_DA_used(i)=0
       enddo
       do i=1,nTocDS
         i_run_DS_used(i)=0
       enddo
       do i=1,nTocIA
         i_run_IA_used(i)=0
       enddo
       do i=1,nTocIS
         i_run_IS_used(i)=0
       enddo
       return
       end

       subroutine fin_run_use
#include "run_use_common.fh"
       Logical isFalcon
       Common /lFalcon/ isFalcon
       parameter (MakeWarn=40)
c       parameter (MakeErr=100)
       Character*16 Label
       Character*60 Line
c       need_abend=0
       do i=1,nTocCA
         if(i_run_CA_used(i).gt.MakeWarn.and..not.isFalcon) then
c           write(6,*) 'cArray label ',i,' used ',
c     &                       i_run_CA_used(i),' times'
           call lookup_label(i,'cArray labels',Label)
           write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,
     &       ';was used ',i_run_CA_used(i),' times'
           Call WarningMessage(1, Line)
         endif
       enddo
       do i=1,nTocDA
         if(i_run_DA_used(i).gt.MakeWarn.and..not.isFalcon) then
c           write(6,*) 'dArray label ',i,' used ',
c     &                      i_run_DA_used(i),' times'
           call lookup_label(i,'dArray labels',Label)
           write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,
     &       ';was used ',i_run_DA_used(i),' times'
           Call WarningMessage(1,Line)
         endif
       enddo
       do i=1,nTocDS
         if(i_run_DS_used(i).gt.MakeWarn.and..not.isFalcon) then
c           write(6,*) 'dScalar label ',i,' used ',
c     &                      i_run_DS_used(i),' times'
           call lookup_label(i,'dScalar labels',Label)
           write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,
     &       ';was used ',i_run_DS_used(i),' times'
           Call WarningMessage(1,Line)
         endif
       enddo
       do i=1,nTocIA
         if(i_run_IA_used(i).gt.MakeWarn.and..not.isFalcon) then
c           write(6,*) 'iArray label ',i,' used ',
c     &                       i_run_IA_used(i),' times'
           call lookup_label(i,'iArray labels',Label)
           write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,
     &       ';was used ',i_run_IA_used(i),' times'
           Call WarningMessage(1,Line)
         endif
       enddo
       do i=1,nTocIS
         if(i_run_IS_used(i).gt.MakeWarn.and..not.isFalcon) then
c           write(6,*) 'iScalar label ',i,' used ',
c     &                       i_run_IS_used(i), ' times'
           call lookup_label(i,'iScalar labels',Label)
           write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,
     &       ';was used ',i_run_IS_used(i),' times'
           Call WarningMessage(1,Line)
         endif
       enddo
c       if(need_abend.eq.1) call abend()
       return
       end
