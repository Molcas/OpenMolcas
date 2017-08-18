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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine query the existence of scalar data on runfile.          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Qpg_iScalar
*
*> @brief
*>   Check if a field is available on runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine queries the existence of scalar data on runfile.
*>
*> @param[in]  Label Name of field
*> @param[out] Found Was the field found
*>
*> @see ::Put_iScalar
************************************************************************
      Subroutine Qpg_iScalar(Label,Found)
      Implicit None
#include "pg_is_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Logical       Found
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Integer      RecVal(nTocIS)
      Character*16 RecLab(nTocIS)
      Integer      RecIdx(nTocIS)
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      item
      Integer      nTmp,iTmp
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Read info from runfile.                                              *
*----------------------------------------------------------------------*
      Call ffRun('iScalar labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Found=.False.
         Return
      End If
      Call cRdRun('iScalar labels',RecLab,16*nTocIS)
      Call iRdRun('iScalar values',RecVal,nTocIS)
      Call iRdRun('iScalar indices',RecIdx,nTocIS)
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocIS
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
*
* Do we read an old temporary field?
*
      If(item.ne.-1) Then
         If(RecIdx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, querying temporary iScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
         End If
      End If
*----------------------------------------------------------------------*
* Did we manage to find it?                                            *
*----------------------------------------------------------------------*
      Found=.true.
      If(item.eq.-1) Found=.false.
      If(item.ne.-1.and.RecIdx(item).eq.0) Found=.false.
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
