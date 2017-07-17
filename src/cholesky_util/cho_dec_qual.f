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
*  Cho_Dec_Qual
*
*> @author F. Aquilante
*>
*> @details
*> The symmetric positive (semi-)definite matrix
*>
*> \f[ Q_{ab,cd} = Q_{ab,cd} - \sum_J L_{ab,J} L_{cd,J} \f]
*>
*> is evaluated and the diagonal elements are returned in \p QDiag.
*> Subsequently, \p Q is Cholesky decomposed to yield \p NumV vectors
*> \f$ K_{ab,I} \f$ as
*>
*> \f[ Q_{ab,cd} = \sum_I K_{ab,I} K_{cd,I} \f]
*>
*> @note
*> \p Q is destroyed and can not be used after this routine.
*>
*> @param[in]     Diag  Array of diagonal integrals stored as in the first reduced set
*> @param[in]     Lab   Array of Cholesky vectors (symmetry blocked) for the qualified shell pairs
*> @param[in,out] Q     Array (symmetry blocked) containing the matrix of the qualified integrals
*> @param[out]    Kab   Array symmetry blocked of Cholesky vectors resulting from the decomposition
*>                      of the qualified integral matrix
*> @param[out]    iD    Index array. On exit \p iD(k,jSym) contains the index of the qualified diagonals
*>                      from which the \p k -th vector of \p jSym was generated
*> @param[out]    NumV  Array of the number of vectors in each symmetry from the decomposition of the
*>                      qualified integral matrix
*> @param[out]    QDiag Array (symmetry blocked) containing the updated qualified diagonals
************************************************************************
      SUBROUTINE Cho_Dec_Qual(Diag,Lab,Q,Kab,iD,NumV,QDiag)
      Implicit Real*8 (a-h,o-z)

      Real*8   Diag(*),Lab(*),Q(*), QDiag(*)
      Real*8   Kab(*)
      Integer  iD(*),NumV(*)
#include "cholesky.fh"

      Real*8  Dmax(8)
      Integer nVecG(8)
      Logical Sync
      Character*12 SecNam
      Parameter (SecNam = 'Cho_Dec_Qual')
      Parameter (xone = -1.0d0, one = 1.0d0)

      irc=0

C --- Determine the max diagonal (qualified excluded)
C     For 1-center decomposition, this is done later.
      If (Cho_1Center) Then
         Call fZero(Dmax,nSym)
      Else
         Sync = .False.
         Call Cho_P_MaxDX(Diag,Sync,Dmax)
      End If

C --- Set the number of previous Cholesky vectors (global)
      Call Cho_P_GetGV(nVecG,nSym)

      ipD = 1
      ipQ = 1
      ipSQL = 1
      ipSQK = 1
      ipQD = 0
      Do jSym=1,nSym

         mQ = Max(nQual(jSym),1)

C ---  Do the subtraction  Q({ab}|{cd}) -= sum_J  L({ab},J) * L({cd},J)
C ---------------------------------------------------------------------
         Call DGEMM_('N','T',nQual(jSym),nQual(jSym),nVecG(jSym),
     &                      xone,Lab(ipSQL),mQ,
     &                           Lab(ipSQL),mQ,
     &                       one,Q(ipQ),mQ)

C --- Extract diagonal of updated Q({ab}|{cd})
C --------------------------------------------
         Do iQ = 1,nQual(jSym)
            QDiag(ipQD+iQ) = Q(ipQ-1+nQual(jSym)*(iQ-1)+iQ)
         End Do

C --- For 1-Center decomposition, find max. diagonal among qualified
C ------------------------------------------------------------------
         If (Cho_1Center) Then
            Do iQ = 1,nQual(jSym)
               Dmax(jSym) = max(Dmax(jSym),QDiag(ipQD+iQ))
            End Do
         End If

C --- Cholesky decompose Q({ab}|{cd})
C     Note: Q is destroyed in CD_InCore_p
C ---------------------------------------
         Thr = max(ThrCom,Dmax(jSym)*Span)
         Call CD_InCore_p(Q(ipQ),nQual(jSym),Kab(ipSQK),
     &                    nQual(jSym),iD(ipD),NumV(jSym),Thr,irc)

         If (irc.ne.0) Then
            write(6,*) SecNam,' non-zero rc on exit from CD_InCore_p: ',
     &                 irc
            Call Cho_Quit('Decomposition error in '//SecNam,104)
         EndIf

         ipSQL = ipSQL + nQual(jSym)*nVecG(jSym)
         ipSQK = ipSQK + nQual(jSym)**2
         ipQ = ipQ + nQual(jSym)**2
         ipD = ipD + nQual(jSym)
         ipQD = ipQD + nQual(jSym)

      End Do

      Return
      End
