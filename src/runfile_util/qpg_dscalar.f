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
*  Qpg_dScalar
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
*> @see ::Put_dScalar
************************************************************************
      Subroutine Qpg_dScalar(Label,Found)
      Implicit None
#include "pg_ds_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Logical       Found
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Real*8       RecVal(nTocDS)
      Character*16 RecLab(nTocDS)
      Integer      RecIdx(nTocDS)
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
      Call ffRun('dScalar labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Found=.False.
         Return
      End If
      Call cRdRun('dScalar labels',RecLab,16*nTocDS)
      Call dRdRun('dScalar values',RecVal,nTocDS)
      Call iRdRun('dScalar indices',RecIdx,nTocDS)
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocDS
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
            Write(6,*) '*** Warning, querying temporary dScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
#ifdef _BIGOT_
            Call AbEnd()
#endif
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
