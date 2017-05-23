/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
/* Here we exploit the capability of C to define a
 * pointer to a procedure that is selected just once
 * but called several times.
 * Hence we avoid repeated IF tests and floating flags that
 * F90 would demand.  */

#include <stdio.h>
#include "molcastype.h"

#ifdef _CAPITALS_
#define fmm_store_w_contractor FMM_STORE_W_CONTRACTOR
#define fmm_selected_w_contractor FMM_SELECTED_W_CONTRACTOR
#define fmm_store_w_buffer FMM_STORE_W_BUFFER
#define fmm_selected_w_buffer FMM_SELECTED_W_BUFFER
#define fmm_store_test FMM_STORE_TEST
#define fmm_included_pair FMM_INCLUDED_PAIR
#define fmm_stored_t_pair_mould FMM_STORED_T_PAIR_MOULD
#define fmm_store_t_pair_mould1 FMM_STORE_T_PAIR_MOULD1
#define fmm_store_t_pair_mould2 FMM_STORE_T_PAIR_MOULD2
#define fmm_store_t_pair_mould3 FMM_STORE_T_PAIR_MOULD3
#define fmm_store_t_pair_mould4 FMM_STORE_T_PAIR_MOULD4
#define fmm_store_t_contractor FMM_STORE_T_CONTRACTOR
#define fmm_selected_t_contractor FMM_SELECTED_T_CONTRACTOR
#define fmm_store_t_buffer FMM_STORE_T_BUFFER
#define fmm_selected_t_buffer FMM_SELECTED_T_BUFFER
#else
#ifndef ADD_
#define fmm_store_w_contractor fmm_store_w_contractor_
#define fmm_selected_w_contractor fmm_selected_w_contractor_
#define fmm_store_w_buffer fmm_store_w_buffer_
#define fmm_selected_w_buffer fmm_selected_w_buffer_
#define fmm_store_test fmm_store_test_
#define fmm_included_pair fmm_included_pair_
#define fmm_stored_t_pair_mould fmm_stored_t_pair_mould_
#define fmm_store_t_pair_mould1 fmm_store_t_pair_mould1_
#define fmm_store_t_pair_mould2 fmm_store_t_pair_mould2_
#define fmm_store_t_pair_mould3 fmm_store_t_pair_mould3_
#define fmm_store_t_pair_mould4 fmm_store_t_pair_mould4_
#define fmm_store_t_contractor fmm_store_t_contractor_
#define fmm_selected_t_contractor fmm_selected_t_contractor_
#define fmm_store_t_buffer fmm_store_t_buffer_
#define fmm_selected_t_buffer fmm_selected_t_buffer_
#endif
#endif

typedef void (*arg1)(void*);
typedef void (*arg2)(void*, void*);
typedef void (*arg3)(void*, void*, void*);
typedef void (*arg4)(void*, void*, void*, void*);
typedef void (*arg5)(void*, void*, void*, void*, void*);

static arg1 T_contractor;
static arg1 W_contractor;
static arg2 T_buffer;
static arg2 W_buffer;
static arg3 stored_test_pq;
static arg4 T_mould1;
static arg3 T_mould2;
static arg3 T_mould3;
static arg5 T_mould4;

/*---------------------------------------------------------------------------*/

void fmm_store_t_contractor_(arg1 T_contractor_in)
{
  T_contractor = T_contractor_in;
}
void fmm_selected_t_contractor_(void*a)
{
  T_contractor(a);
}

void fmm_store_t_buffer_(arg2 T_buffer_in)
{
  T_buffer = T_buffer_in;
}
void fmm_selected_t_buffer_(void*a, void*b)
{
  T_buffer(a,b);
}

/*---------------------------------------------------------------------------*/

void fmm_store_w_contractor_(arg1 W_contractor_in)
{
  W_contractor = W_contractor_in;
}
void fmm_selected_w_contractor_(void*a)
{
  W_contractor(a);
}

void fmm_store_w_buffer_(arg2 W_buffer_in)
{
  W_buffer = W_buffer_in;
}
void fmm_selected_w_buffer_(void*a, void*b)
{
  W_buffer(a,b);
}

/*---------------------------------------------------------------------------*/

void fmm_store_test_(arg3 test_pq_in)
{
  stored_test_pq = test_pq_in;
}
void fmm_included_pair_(void*a, void*b, void*c)
{
  stored_test_pq(a,b,c);
}

/*---------------------------------------------------------------------------*/

void fmm_store_t_pair_mould1_(arg4 T_mould_in)
{
  T_mould1 = T_mould_in;
}
void fmm_store_t_pair_mould2_(arg3 T_mould_in)
{
  T_mould2 = T_mould_in;
}
void fmm_store_t_pair_mould3_(arg3 T_mould_in)
{
  T_mould3 = T_mould_in;
}
void fmm_store_t_pair_mould4_(arg5 T_mould_in)
{
  T_mould4 = T_mould_in;
}
void fmm_stored_t_pair_mould_(void*LHS,void*RHS, void*id,void*wt,void*T_pair)
{
  T_mould1(LHS,RHS,id,T_pair);
  T_mould2(LHS,id,T_pair);
  T_mould3(RHS,id,T_pair);
  T_mould4(LHS,RHS,id,wt,T_pair);
}

/*---------------------------------------------------------------------------*/
