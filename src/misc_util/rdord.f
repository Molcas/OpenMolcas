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
******************************************************************
*
* This subroutine substitutes the previous version of the RdOrd
* to allow the generation of two-elwctron integrals from a set
* of Cholesky vectors
*
*     Data declarations:
*
*        First  : logical variable for identifying the 1st call
*                 to this subroutine
*
*   F. Aquilante
*
******************************************************************

      Subroutine RdOrd(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
******************************************************************
*
*   <Dummy DOC>
*     <Name>RdOrd</Name>
*     <Syntax>Call RdOrd(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)</Syntax>
*     <Arguments>
*       \Argument{rc}{}{}{out}
*       \Argument{iOpt}{}{}{inout}
*       \Argument{iSym}{}{}{in}
*       \Argument{jSym}{}{}{in}
*       \Argument{kSym}{}{}{in}
*       \Argument{lSym}{}{}{in}
*       \Argument{Buf}{}{}{inout}
*       \Argument{lBuf}{}{}{in}
*       \Argument{nMat}{}{}{out}
*     </Arguments>
*     <Purpose>Driver for the actual RdOrd_ routines and the
*              Cholesky vector integral generator</Purpose>
*     <Dependencies></Dependencies>
*     <Author>F. Aquilante</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     </Description>
*    </DOC>
*
*****************************************************************

*----------------------------------------------------------------------*
*                                                                     *
*                                                                     *
*    Calling parameters:                                              *
*    iSym   : irred. representation of "first" symmetry label         *
*    jSym   : irred. representation of "second" symmetry label        *
*    kSym   : irred. representation of "third" symmetry label         *
*    lSym   : irred. representation of "fourth" symmetry label        *
*    Buf    : containts on output the integrals                       *
*    lBuf   : length of the integral buffer                           *
*    nMat   : number of submatrices read in                           *
*    iOpt   : option code (iOpt=1:start reading at first shell        *
*                      distribution in the given product symmetry)    *
*                         (iOpt=2:continue reading)                   *
*    rc     : return code                                             *
*                                                                     *
*    The integrals are returned in the order (lk|ji)                  *
*                                                                     *
*                                                                     *
*----------------------------------------------------------------------*


      Implicit Real*8 (a-h,o-z)

      Integer rc

      Logical First,DoCholesky
      Common / DoCho / DoCholesky

      Data First /.true./
      Save First

#include "RdOrd.fh"

      If (First) Then
        CALL DecideOnCholesky(DoCholesky)
        If (DoCholesky) Then
c   initialize information
          CALL INIT_GETINT(rc)
        End If
        First=.false.
      End If

      If (DoCholesky) Then
       CALL Get_Int(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
      Else
       CALL RdOrd_(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
      end if

c --- Debug printing
c
c      call get_iscalar('nSym',nsym)
c      call get_iarray('nbas',nbas,nsym)
c      if (kSym .eq. lSym) then
c         Nkl = nBas(kSym)*(nBas(kSym)+1)/2
c      else
c         Nkl = nBas(kSym)*nBas(lSym)
c      end if
c      write(6,*)
c      write(6,*)'Integrals: symm(ij|kl), nMat ',iSym,jSym,kSym,lSym,nMat
c      write(6,*)'Nbask,Nbasl,Nkl= ',nBas(kSym),nBas(lSym),Nkl
c      xnrm = sqrt(ddot_(Nkl*nMat,Buf,1,Buf,1))
c      write(6,*)'norm= ',xnrm
c      call cho_output(Buf,1,nMat,1,Nkl,nMat,Nkl,1,6)
c      call xflush(6)
c ---

      Return
      END
