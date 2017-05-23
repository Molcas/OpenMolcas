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
* Copyright (C) 2001-2005, Valera Veryazov                             *
************************************************************************
      Subroutine molcas_open(Lu,Name)
      Integer Lu,f_recl,f_iostat
      Character*10 f_access,f_form,f_status
      Logical is_recl,is_error
      Character*(*) Name

      f_recl=1
      f_iostat=100
      f_access='SEQUENTIAL'
      f_form='FORMATTED'
      f_status='UNKNOWN'
      is_recl=.false.
      is_error=.false.

      Call molcas_open_ext2(Lu,Name,f_access,f_form,f_iostat,is_recl,
     &                      f_recl,f_status,is_error)
      If(f_iostat.ne.0) Then
         Write(6,*)
         Write(6,'(3a)')   'molcas_open: Error opening file "',Name,'"'
         Write(6,'(a,i9)') '   iostat is',f_iostat
         Write(6,'(a)')    '   Aborting'
         Write(6,*)
         Call Abend()
      End If
      Return
      End
