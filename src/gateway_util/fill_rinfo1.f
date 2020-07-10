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
      Subroutine Fill_rInfo1()
      use Basis_Info
#include "itmax.fh"
#include "info.fh"
#include "rinfo.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Generate stuff for rinfo.fh
*
      krExp=0
      krCof=0
      krBas=0
      krCnt=0

*     Loop over basis sets
      Do iCnttp = 1, nCnttp
*        Loop over distinct centers
         Do icnt = 1, dbsc(iCnttp)%nCntr
            krCnt=krCnt+1
            nAngr(krCnt)=nVal_Shells(iCnttp)-1
*
*           Start with s type shells
            jSh = ipVal(iCnttp)
            Do iAng = 0, nVal_Shells(iCnttp)-1
*
               krBas=krBas+1
               If (krBas.gt.MxAO) Then
                  Call WarningMessage(2,'Too many shells')
                  Write(6,*) 'MORE THAN ',MxAO,' SHELLS'
                  Write(6,*) 'Increase MxAO in info.fh and',
     &                      ' recompile the code!'
                  Call Abend()
               End If
               nPrimr(krBas)=nExp(jSh)
               nBasisr(krBas)=nBasis_Cntrct(jSh)
*
               If (krExp+nExp(jSh).gt.MxPrim) then
                  Call WarningMessage(2,'Too many primitives')
                  write(6,*) 'MORE THAN ',MxPrim,' PRIMITIVES'
                  write(6,*) 'Increase MxPrim in rinfo.fh and',
     &                       'recompile the code!'
                  Call Abend()
               End If
               Do  kExp=1,nExp(jSh)
                  krExp=krExp+1
                  rExp(krExp)=Shells(jSh)%Exp(kExp)
               End Do
*
*
*              Pointer to the untouched contraction matrix as after
*              input.
*
               If (krCof+nExp(jSh)*nBasis(jSh).gt.MxrCof) Then
                  Call WarningMessage(2,
     &                     'Too many contraction coefficients')
                  Write(6,*) 'MORE THAN ',MxrCof,
     &                      ' CONTRACTION COEFFICIENTS'
                  Write(6,*) 'Increase MxrCof in rinfo.fh and',
     &                       'recompile the code!'
                  Call Abend()
               End If
               Do kCof=1,nBasis_Cntrct(jSh)
                  Do  kExp=1,nExp(jSh)
                        krCof=krCof+1
                        rCof(krCof)=Shells(jSh)%Cff_c(kExp,kCof,2)
                  End Do
               End Do
*
               jSh = jSh + 1
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
