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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
***************************************************************
*
* Author :   F. Aquilante
*
*  Get_Int :  driver for the integral generator from Cholesky
*             vectors
***************************************************************
      SUBROUTINE Get_Int(rc,iOpt,iSymp,iSymq,iSymr,iSyms,Xint,lBuf,nMat)

      Implicit Real*8 (a-h,o-z)
      Integer  rc,iOpt
      Integer  iSymp,iSymq,iSymr,iSyms,Npq,Nrs
      Real*8   Xint(*)
      Integer  lBuf,nMat

#include "RdOrd.fh"
#include "TwoRc.fh"

      Character*4 BaseNm
      Character*6 Fname
      Parameter (BaseNm = 'CHFV')

      MulD2h(i,j) = iEor(i-1,j-1)+1

C Check input parameters

      rc = 0
      If (iOpt.ne.1 .and. iOpt.ne.2) Then
         rc = rcRD06
         Write(6,*) 'Get_Int: Invalid option'
         Write(6,*) 'iOpt= ',iOpt
         Call QTrace()
         Call Abend()
      End If
      If (iSymp.lt.iSymq .or. iSymr.lt.iSyms) Then
         rc = rcRD02
         Write(6,*) 'Get_Int: invalid order of symmetry labels'
         Call Qtrace()
         Call Abend()
      End If
      If (MulD2h(iSymp,iSymq) .ne. MulD2h(iSymr,iSyms)) Then
         rc = rcRD01
         Write(6,*) 'Get_Int: wrong symmetry labels, direct product',
     &              ' is not total symmetric'
         Call Qtrace()
         Call Abend()
      End If
      If (lBuf.lt.1) Then
         rc = rcRD04
         Write(6,*) 'Get_Int: invalid buffer size'
         Write(6,*) 'lBuf=',lBuf
         Call Qtrace()
         Call Abend()
      End If

C Open files.
      LuCVec(1) = 7
      Write(Fname,'(A4,I1,I1)') BaseNm,iSymp,iSymq
      Call DANAME_MF_WA (LuCVec(1),Fname)
      if (iSymp.ne.iSymr) then
          LuCVec(2) = 7
          Write(Fname,'(A4,I1,I1)') BaseNm,iSymr,iSyms
          Call DANAME_MF_WA (LuCVec(2),Fname)
      Else
          LuCVec(2) = -1
      end if


         If (iSymp .eq. iSymq) Then
             Npq = nBas(iSymp)*(nBas(iSymp) + 1)/2
         Else
             Npq = nBas(iSymp)*nBas(iSymq)
         End If
         If (iSymr .eq. iSyms) Then
             Nrs = nBas(iSymr)*(nBas(iSyms) + 1)/2
         Else
             Nrs = nBas(iSymr)*nBas(iSyms)
         End If

C --- For debug
C      lBufs=lBuf
C      lBuf=  min(Nrs*min(Npq,2)+1,lBufs)

      If (iOpt.eq.1) then
       pq1  = 1
       nMat = min(Npq,(lBuf - 1)/Nrs)
      else
       If (pq1.lt.1 .or. pq1.gt.Npq) Then
          rc = 999999
          Write(6,*) 'pq1 out of bounds: ',pq1
          Call Qtrace()
          Call Abend()
          nMat = 99999999
       Else
          nMat = min((Npq-pq1+1),(lBuf - 1)/Nrs)
       End If
      end if

       If (nMat.lt.1) Then
          return !no more integrals to compute
       Else
          CALL GEN_INT(rc,iSymp,iSymq,iSymr,iSyms,pq1,nMat,Xint)
       End If
       pq1 = pq1 + nMat

C Close files.
      Do i=1,2
       if (LuCVec(i).ne.-1) then
          Call DACLOS(LuCVec(i))
          LuCVec(i) = -1
       end if
      End do

C      LBuf = lBufs

      Return
      End
