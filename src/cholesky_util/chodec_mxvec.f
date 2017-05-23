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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine ChoDec_MxVec(CD_Col,CD_Vec,MxVec,
     &                  Restart,Thr,Span,MxQual,
     &                  Diag,Qual,Buf,
     &                  iPivot,iQual,
     &                  nDim,lBuf,
     &                  ErrStat,NumCho,
     &                  irc)
************************************************************
*
*   <DOC>
*     <Name>ChoDec\_MxVec</Name>
*     <Syntax>Call ChoDec\_MxVec(CD\_Col,CD\_Vec,MxVec,
*                         Restart,Thr,Span,MxQual,Diag,Qual,
*                         Buf,iPivot,iQual,nDim,lBuf,
*                         ErrStat,NumCho,irc)
*     </Syntax>
*     <Arguments>
*       \Argument{CD\_Col}{External subroutine for retrieving
*                           matrix columns}{}{in}
*       \Argument{CD\_Vec}{External subroutine for filling and
*                           emptying Cholesky vector buffer}{}{in}
*       \Argument{MxVec}{Maximum number of vectors to be generated}
*                         {Integer}{in}
*       \Argument{Restart}{Flag for restarting decomposition}
*                         {Logical}{in}
*       \Argument{Thr}{Decomposition threshold (precision)}
*                     {Real*8}{inout}
*       \Argument{Span}{Span factor}{Real*8}{inout}
*       \Argument{MxQual}{Max. number of qualified columns
*                         per decomposition pass}{Integer}{in}
*       \Argument{Diag}{Array containing the diagonal of the matrix
*                       to be decomposed, dimension is Diag(nDim)}
*                      {Real*8}{in}
*       \Argument{Qual}{Storage array for matrix columns,
*                       dimension is Qual(nDim,0:MxQual)}
*                      {Real*8}{in}
*       \Argument{Buf}{Buffer for storing Cholesky vectors,
*                      dimension is Buf(lBuf)}
*                     {Real*8}{in}
*       \Argument{iPivot}{Pivoting array, dimension is iPivot(nDim)}
*                {Integer}{in}
*       \Argument{iQual}{Array for storing qualified indices,
*                        dimension is iQual(MxQual)}
*                       {Integer}{in}
*       \Argument{nDim}{Linear dimension of matrix}{Integer}{in}
*       \Argument{lBuf}{Buffer size}{Integer}{in}
*       \Argument{ErrStat}{Min, max, and RMS errors, respectively,
*                          for decomposition as judged from the residual
*                          diagonal elements, dimension is ErrStat(3)}
*                         {Real*8}{out}
*       \Argument{NumCho}{Number of Cholesky vectors}{Integer}{inout}
*       \Argument{irc}{Return code}{Integer}{out}
*     </Arguments>
*     <Purpose>Decompose a symmetric positive
*              (semi-) definite matrix</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        ChoDec\_MxVec will Cholesky decompose a positive (semi-)definite
*        matrix without knowledge of all its elements (i.e., using a
*        matrix-direct approach).
*
*        Given the diagonal, the idea is that ChoDec\_MxVec asks the
*        caller for specified columns of the matrix when needed during
*        decomposition. This is done by invoking the external routine
*        passed as CD\_Col. Similarly, the caller must provide an
*        external routine, passed as CD\_Vec, that handles the Cholesky
*        vector buffer - that is, empties or fills the buffer (CD\_Vec
*        could f.ex. do Cholesky vector I/O).
*
*        Unlike ChoDec, however, this routine will stop when either of
*        the following conditions are met:
*
*           1. max. diagonal <= Thr
*           2. number of vectors = MxVec
*
*        The external subroutine CD\_Col must have the following call
*        pattern:
*
*           Call CD\_Col(Col,nDim,iCol,nCol,Buf,lBuf)
*
*        where Col(nDim,nCol) contains the columns, iCol(nCol) specifies
*        the columns requested [i.e., column Col(*,i) should correspond
*        to column iCol(i) in the full matrix], and Buf(lBuf) can be
*        used as (perhaps additional) work space.
*
*        The external subroutine CD\_Vec must have the following call
*        pattern:
*
*           Call CD\_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)
*
*        where iVec1 is the index of the first vector, nVec the number
*        of vectors, and nDim the vector dimension. The buffer Buf has
*        total dimension lBuf. For iOpt=1 ("empty buffer option"),
*        the first nDim*nVec elements of Buf contain the vectors on
*        input - i.e. CD\_Vec is expected to move the vectors to another
*        location so that all of Buf is again available for use.
*        For iOpt=2 ("fill buffer option"), the process is reversed and
*        ChoDec expects CD\_Vec to return vectors
*        iVec1,iVec1+1,...,iVec1+nVec-1, in the first nDim*nVec elements
*        of Buf.
*
*        The error code irc is 0 if success, while a non-zero number
*        is returned if some failure (or if the matrix is not positive
*        (semi-) definite) occurs. In the latter cases, the
*        decomposition is ill-defined and the data returned is junk!
*
*        Note that Thr and Span may be reset in this routine (if they
*        are not properly defined on input).
*     </Description>
*    </DOC>
*
************************************************************

      Implicit Real*8 (a-h,o-z)
      External  CD_Col    ! external routine for matrix columns
      External  CD_Vec    ! external routine for Cholesky vectors
      Logical   Restart
      Dimension Diag(nDim), Qual(nDim,0:MxQual), Buf(lBuf)
      Integer   iPivot(nDim), iQual(MxQual)
      Dimension ErrStat(3)

      Logical   Converged

      Character*12 SecNam
      Parameter (SecNam = 'ChoDec_MxVec')

      Parameter (DefThr = 1.0d-6, DefSpan = 1.0d-2) ! defaults
      Parameter (ThrNeg = -1.0d-13, ThrFail = -1.0d-8)

      Call qEnter(SecNam)

C     Initialize variables.
C     ---------------------

      irc = 0
      ErrStat(1) =  9.876543210d15
      ErrStat(2) = -9.876543210d15
      ErrStat(3) = -9.876543210d15
      If (.not. Restart) NumCho = 0
      Converged = .false.

C     Check dimensions.
C     -----------------

      If (nDim .lt. 1) Then
         irc = 0
         Go To 1  ! exit (nothing to do)
      End If
      If (MxQual .lt. 1) Then
         irc = -1
         Go To 1  ! exit (qualification not possible)
      End If
      mQual  = min(MxQual,nDim)
      MinBuf = nDim + mQual
      If (lBuf .lt. MinBuf) Then
         irc = -2
         Go To 1 ! exit (buffer too small)
      End If
      If (MxVec .lt. 1) Then
         irc = -3
         Go To 1 ! exit (MxVec too small)
      End If
      MxNumCho = min(MxVec,nDim)

C     Check and possibly reset configuration.
C     ---------------------------------------

      If (Thr .lt. 0.0d0) Then
         Thr = DefThr  ! reset threshold
      End If
      If (Span.lt.0.0d0 .or. Span.gt.1.0d0) Then
         Span = DefSpan ! reset span factor
      End If

C     Set up a (possibly updated) copy of the diagonal.
C     Test diagonal for negative elements.
C     -------------------------------------------------

      Call CD_Diag(CD_Vec,
     &             Restart,Converged,Thr,
     &             ThrNeg,ThrFail,
     &             Diag,Qual(1,0),Buf,
     &             nDim,lBuf,
     &             ErrStat,NumCho,
     &             irc)
      If (irc .ne. 0) Go To 1 ! exit (initial diagonal error)

C     Only do decomposition if not converged!
C     ---------------------------------------

      If (.not. Converged .and. NumCho.lt.MxNumCho) Then

C        Cholesky decompose.
C        -------------------

         Call CD_Decomposer(CD_Col,CD_Vec,MxNumCho,
     &                      Thr,Span,mQual,
     &                      ThrNeg,ThrFail,
     &                      Qual(1,0),Qual(1,1),Buf,
     &                      iPivot,iQual,
     &                      nDim,lBuf,
     &                      NumCho,
     &                      irc)
         If (irc .ne. 0) Go To 1  ! exit (decomposition error)

C        Check diagonal and set up error statistics.
C        -------------------------------------------

         Call CD_Diag(CD_Vec,
     &                .true.,Converged,Thr,
     &                ThrNeg,ThrFail,
     &                Diag,Qual(1,0),Buf,
     &                nDim,lBuf,
     &                ErrStat,NumCho,
     &                irc)
         If (irc .ne. 0) Then
            irc = irc + 200 ! makes the error code 400 + n
            Go To 1 ! exit (check error)
         Else If (.not. Converged) Then
            If (NumCho .lt. MxNumCho) Then
               irc = 1
               Go To 1 ! exit (decomposition failure)
            Else If (NumCho .gt. MxNumCho) Then
               Call qTrace
               Call SysAbendMsg(SecNam,'Logical error!',' ')
            End If
         End If

      End If

C     That's it!
C     ----------

    1 Continue
      Call qExit(SecNam)
      End
