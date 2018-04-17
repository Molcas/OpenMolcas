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
* Copyright (C) 2008, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine get scalar double data from the peek/poke buffer.       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: May 2008                                                    *
*                                                                      *
************************************************************************
*  Peek_dScalar
*
*> @brief
*>   Get scalar double data from the peek/poke buffer
*> @author Per-Olof Widmark
*>
*> @details
*> This routine is used to get scalar data of type ``Real*8``
*> from the peek/poke buffer. The data items are identified
*> by a text label.
*>
*> @param[in]  Label Name of field
*> @param[out] Data  Data to get from runfile
*>
*> @copyright All rights reserved by Lund University
************************************************************************
      Subroutine Peek_dScalar(Label,Data)
      Implicit None
#include "pp_ds_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Real*8        Data
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Logical Found
      Integer indx
      Integer i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
      Found=.false.
*     Write(6,'(2a)') 'peek_dscalar: Label is ',Label
*     Write(6,'(a,i8)') 'peek_dscalar: ds_no is ',ds_no
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      indx=-1
      Do i=1,ds_no
         If(ds_label(i).eq.Label) indx=i
      End Do
*     Write(6,'(a,i8)') 'peek_dscalar: indx is ',indx
*----------------------------------------------------------------------*
* Get data from buffer.                                                *
*----------------------------------------------------------------------*
      If(indx.eq.-1) Then
         If(ds_no.ge.nTabDS) Then
            Call SysAbendMsg('Peek_dScalar',
     &         'Too many fields',
     &         'Increase nTabDS and recompile')
         End If
         ds_no=ds_no+1
         indx=ds_no
         Call Qpg_dScalar(Label,Found)
         If(Found) Then
            Call Get_dScalar(Label,Data)
         Else
            Call SysAbendMsg('Peek_dScalar',
     &         'Field not found',
     &         Label)
         End If
         ds_label(indx)=Label
         ds_value(indx)=Data
      Else
         Data=ds_value(indx)
      End If
*     Write(6,'(a,e20.8)') 'peek_dscalar: Data is ',Data
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
