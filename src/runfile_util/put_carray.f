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
* This routine put array character data to the runfile.                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Put_cArray
*
*> @brief
*>   Add/update array data in runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine is used to put array data of type
*> ``Character`` into the runfile. The data items are
*> identified by the \p label. Below is a list of the
*> data items that are recognized. The labels are
*> case insensitive and significant to 16 characters.
*>
*> For development purposes you can use an unsupported
*> label. Whenever such a field is accessed a warning
*> message is printed in the output, to remind the
*> developer to update this routine.
*>
*> List of known labels:
*>
*> '``DFT functional``'     Name of the functional used for the KS-DFT calculation.
*> '``Irreps``'             Names of the irreducible representations.
*> '``Relax Method``'       Name of the method used for geometry optimizations.
*> '``Seward Title``'       The title of the calculation as specified in module SEWARD.
*> '``Slapaf Info 3``'      Misc. information for module SLAPAF.
*> '``Unique Atom Names``'  List of the names of the symmetry unique atoms.
*> '``Unique Basis Names``' List of the basis function names.
*> '``MkNemo.lMole``'       The labels of molecules as specified in mknemo module.
*> '``MkNemo.lCluster``'    The labels of clusters as specified in mknemo module.
*> '``MkNemo.lEnergy``'     The labels of energies as specified in mknemo module.
*>
*> @param[in] Label Name of field
*> @param[in] Data  Data to put on runfile
*> @param[in] nData Length of array
************************************************************************
      Subroutine Put_cArray(Label,Data,nData)
      Implicit None
#include "pg_ca_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Character*16  myLabel
      Integer       nData
      Character*(*) Data
cvv      Character*(*) Data(nData)
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Character*16 RecLab(nTocCA)
      Integer      RecIdx(nTocCA)
      Integer      RecLen(nTocCA)
      Save         RecLab
      Save         RecIdx
      Save         RecLen
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      nTmp
      Integer      item
      Integer      iTmp
      Integer      i, ilen
*----------------------------------------------------------------------*
* Do setup if this is the first call.                                  *
*----------------------------------------------------------------------*
      myLabel=' '
      ilen=len(Label)
      ilen=min(ilen,16)
      myLabel=Label(1:ilen)
      Call ffRun('cArray labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Do i=1,nTocCA
            RecLab(i)=' '
            RecIdx(i)=sNotUsed
            RecLen(i)=0
         End Do
*
*        Observe that label is at most 16 characters!
*
*                     1234567890123456
         RecLab(  1)='DFT functional  '
         RecLab(  2)='Irreps          '
         RecLab(  3)='Relax Method    '
         RecLab(  4)='Seward Title    '
         RecLab(  5)='Slapaf Info 3   '
         RecLab(  6)='Unique Atom Name' !s
         RecLab(  7)='Unique Basis Nam' !es
         RecLab(  8)='LP_L            '
         RecLab(  9)='MkNemo.lMole    '
         RecLab( 10)='MkNemo.lCluster '
         RecLab( 11)='MkNemo.lEnergy  '
         RecLab( 12)='Symbol ZMAT     '
         RecLab( 13)='Tinker Name     '
         RecLab( 14)='ESPF Filename   '
         RecLab( 15)='ChDisp          '
         RecLab( 16)='cmass           '
         RecLab( 17)='BirthCertificate'
         RecLab( 18)='LastEnergyMethod'
         RecLab( 19)='MMO Labels      '
         RecLab( 20)='MCLR Root       '
         RecLab( 21)='Frag_Type       ' ! EFP fragment labels
         RecLab( 22)='ABC             ' ! EFP atom labels
         RecLab( 23)='Un_cen Names    '
*                     1234567890123456
         Call cWrRun('cArray labels',RecLab,16*nTocCA)
         Call iWrRun('cArray indices',RecIdx,nTocCA)
         Call iWrRun('cArray lengths',RecLen,nTocCA)
      Else
         Call cRdRun('cArray labels',RecLab,16*nTocCA)
         Call iRdRun('cArray indices',RecIdx,nTocCA)
         Call iRdRun('cArray lengths',RecLen,nTocCA)
      End If
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=myLabel
      Call UpCase(CmpLab1)
      Do i=1,nTocCA
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
*
* Do we create a new temporary field?
*
      If(item.eq.-1) Then
         Do i=1,nTocCA
            If(RecLab(i).eq.' ') item=i
         End Do
         If(item.ne.-1) Then
            RecLab(item)=myLabel
            RecIdx(item)=sSpecialField
            Call cWrRun('cArray labels',RecLab,16*nTocCA)
            Call iWrRun('cArray indices',RecIdx,nTocCA)
         End If
      End If
*
* Is this a temporary field?
*
      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, writing temporary cArray field'
            Write(6,*) '***   Field: ',myLabel
            Write(6,*) '***'
#ifdef _BIGOT_
            Call AbEnd()
#endif
         End If
      End If
*----------------------------------------------------------------------*
* Write data to disk.                                                  *
*----------------------------------------------------------------------*
      If(item.eq.-1) Then
         Call SysAbendMsg('put_cArray','Could not locate',myLabel)
      End If
      Call cWrRun(RecLab(item),Data,nData)
      If(RecIdx(item).eq.0) Then
         RecIdx(item)=sRegularField
         Call iWrRun('cArray indices',RecIdx,nTocCA)
      End If
      If(RecLen(item).ne.nData) Then
         RecLen(item)=nData
         Call iWrRun('cArray lengths',RecLen,nTocCA)
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
