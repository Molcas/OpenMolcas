/***********************************************************************
 * This file is part of OpenMolcas.                                     *
 *                                                                      *
 * OpenMolcas is free software; you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License, v. 2.1. *
 * OpenMolcas is distributed in the hope that it will be useful, but it *
 * is provided "as is" and without any express or implied warranties.   *
 * For more details see the full text of the license in the file        *
 * LICENSE or in <http://www.gnu.org/licenses/>.                        *
 *                                                                      *
 * Copyright (C) 2015, Marcus Johansson                                 *
 ***********************************************************************/

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>

#include "molcastype.h"

#include "msym.h"

#ifdef _CAPITALS_
#define cmsym_create_context CMSYM_CREATE_CONTEXT
#define cmsym_release_context CMSYM_RELEASE_CONTEXT
#define cmsym_set_elements CMSYM_SET_ELEMENTS
#define cmsym_find_symmetry CMSYM_FIND_SYMMETRY
#define cmsym_symmetrize_molecule CMSYM_SYMMETRIZE_MOLECULE
#define cmsym_symmetrize_orbitals CMSYM_SYMMETRIZE_ORBITALS
#define cmsym_generate_orbital_subspaces CMSYM_GENERATE_ORBITAL_SUBSPACES

#else
#ifndef ADD_
#define cmsym_create_context cmsym_create_context_
#define cmsym_release_context cmsym_release_context_
#define cmsym_set_elements cmsym_set_elements_
#define cmsym_find_symmetry cmsym_find_symmetry_
#define cmsym_symmetrize_molecule cmsym_symmetrize_molecule_
#define cmsym_symmetrize_orbitals cmsym_symmetrize_orbitals_
#define cmsym_generate_orbital_subspaces cmsym_generate_orbital_subspaces_
#endif
#endif

#define MAX_FILE_NAME_LENGTH 255
#define BOHR_TO_ANGSTROM 0.529177249


INT cmsym_create_context(msym_context *pctx, INT *err){
	const char *errstr;
	if(NULL == (*pctx = msymCreateContext())) goto err;
	*err = 0;
	return 0;
err:
	errstr = msymGetErrorDetails();
	fprintf(stderr,"%s\n",errstr);
	*err = 1;
	return 1;
}

INT cmsym_release_context(msym_context *pctx, int*err){
	msym_error_t ret = MSYM_SUCCESS;
	ret = msymReleaseContext(*pctx);
	*pctx = NULL;
	*err = ret;
	return ret;
}

INT cmsym_set_elements(msym_context *pctx, INT *pel, INT *puel, char uelement[][*puel], double xyz[][3], INT *paol, INT basis_ids[][4], INT *err){
	msym_context ctx = *pctx;
	msym_error_t ret = MSYM_SUCCESS;
	const char *errstr;
	int el = *pel, aol = *paol, uel = *puel, mel = 0;
	msym_element_t *elements = calloc(el,sizeof(msym_element_t));
	int *ai = malloc(sizeof(int[el]));

	for(int i = 0; i < el;i++){
		char buf[uel+1];
		memcpy(buf,uelement[i],uel);
		buf[uel] = 0;
		sscanf(buf,"%[A-Za-z]%d",elements[i].name,&ai[i]);
		elements[i].v[0] = BOHR_TO_ANGSTROM*xyz[i][0];
		elements[i].v[1] = BOHR_TO_ANGSTROM*xyz[i][1];
		elements[i].v[2] = BOHR_TO_ANGSTROM*xyz[i][2];
	}

	msym_basis_function_t *bfs = calloc(aol,sizeof(msym_basis_function_t));

	if(MSYM_SUCCESS != (ret = msymSetElements(ctx, el, elements))) goto err;

	for(int i = 0;i < aol;i++){
		bfs[i].element = &elements[basis_ids[i][0]-1];
		bfs[i].type = MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC;
		bfs[i].f.rsh.l = basis_ids[i][2];
		bfs[i].f.rsh.n = basis_ids[i][1] + bfs[i].f.rsh.l;
		bfs[i].f.rsh.m = basis_ids[i][3];
	}

	if(MSYM_SUCCESS != (ret = msymSetBasisFunctions(ctx, aol, bfs))) goto err;

	free(bfs);
	free(ai);
	free(elements);

	*err = ret;
	return ret;

err:
	errstr = msymErrorString(ret);
	fprintf(stderr,"%s",errstr);
	errstr = msymGetErrorDetails();
	fprintf(stderr,": %s\n",errstr);
	free(bfs);
	free(ai);
	free(elements);
	*err = ret;
	return ret;
}


INT cmsym_find_symmetry(msym_context *pctx, char pgname[6], INT *err){
	msym_context ctx = *pctx;
	msym_error_t ret = MSYM_SUCCESS;
	const char *errstr;
	char buf[6];

	if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
	if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), buf))) goto err;

	snprintf(pgname,6,"%s      ",buf);
	pgname[5] = ' ';

	*err = ret;
	return ret;
err:
	errstr = msymErrorString(ret);
	fprintf(stderr,"%s",errstr);
	errstr = msymGetErrorDetails();
	fprintf(stderr,": %s\n",errstr);
	*err = ret;
	return ret;
}

INT cmsym_symmetrize_molecule(msym_context *pctx, char *outfile, INT *err){
	msym_context ctx = *pctx;
	msym_error_t ret = MSYM_SUCCESS;
	const char *errstr = NULL;
	double serr = 0.0;
	msym_element_t *elements = NULL;
	int el = 0;

	if(MSYM_SUCCESS != (ret = msymSymmetrizeElements(ctx,&serr))) goto err;
	if(MSYM_SUCCESS != (ret = msymGetElements(ctx,&el,&elements))) goto err;

	FILE *out = fopen(outfile,"w");
	if(out == NULL){
		fprintf(stderr,"Can't open file %s",outfile);
	} else {
		char buf[14];
		fprintf(out,"%d\n\n",el);

		for(msym_element_t *a = elements; a < (elements + el); a++){
			fprintf(out, "%s ",a->name);
			snprintf(buf,14,"%#+13.10lf",a->v[0]);
			fprintf(out,"%s ", buf);
			snprintf(buf,14,"%#+13.10lf",a->v[1]);
			fprintf(out,"%s ", buf);
			snprintf(buf,14,"%#+13.10lf",a->v[2]);
			fprintf(out,"%s\n", buf);
		}
		fclose(out);
	}

	*err = ret;
	return ret;

err:
	errstr = msymErrorString(ret);
	fprintf(stderr,"%s",errstr);
	errstr = msymGetErrorDetails();
	fprintf(stderr,": %s\n",errstr);
	*err = ret;
	return ret;
}

INT cmsym_generate_orbital_subspaces(msym_context *pctx, INT *l, double c[*l][*l], INT irrep_ids[*l], INT irrep_ind[*l], INT* nirreps, char lbl[*l][8], INT *err){
	msym_context ctx = *pctx;
	msym_error_t ret = MSYM_SUCCESS;
	const char *errstr;
	int mbfsl = 0, mssl = 0;
	msym_basis_function_t *mbfs = NULL;
	const msym_character_table_t *mct = NULL;
	const msym_subrepresentation_space_t *mss = NULL;

	msym_point_group_type_t type = MSYM_POINT_GROUP_TYPE_Kh;
	int n = 0;

	if(MSYM_SUCCESS != (ret = msymGetBasisFunctions(ctx, &mbfsl, &mbfs))) goto err;

	/*///////////////////// TESTING

	  if(MSYM_SUCCESS != (ret = msymGetPointGroupType(ctx, &type, &n))) goto err;
	  if(n == 0 && (MSYM_POINT_GROUP_TYPE_Dnh == type || MSYM_POINT_GROUP_TYPE_Cnv == type)){
	  int lmax = 0;
	  for(int i = 0;i < mbfsl;i++) lmax = mbfs[i].f.sh.l > lmax ? mbfs[i].f.sh.l : lmax;
	  if(MSYM_SUCCESS != (ret = msymSetPointGroupByType(ctx, type, 2*lmax))) goto err;
	  if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
	  }

	////////////////////// END TESTING */

	if(MSYM_SUCCESS != (ret = msymGetSubrepresentationSpaces(ctx, &mssl, &mss))) goto err;
	if(MSYM_SUCCESS != (ret = msymGetCharacterTable(ctx, &mct))) goto err;

	memset(c,0,sizeof(double[*l][*l]));
	memset(lbl,' ', sizeof(char[*l][8]));

	int row = 0, sym = 0;
	for(int i = 0;i < mssl;i++){
		for(int j = 0;j < mss[i].salcl;j++){
			for(int k = 0;k < mss[i].salc[j].fl;k++){
				int index = (int)(mss[i].salc[j].f[k] - mbfs);
				double (*pf)[mss[i].salc[j].fl] = mss[i].salc[j].pf;
				for(int l = 0;l < mss[i].salc[j].d;l++){
					c[row+l][index] = pf[l][k];
					irrep_ids[row+l] = mss[i].s;
					irrep_ind[row+l] = sym + l;
					if(mct->s[mss[i].s].d > 1) snprintf(lbl[sym+l], 8, "%d%s",l,mct->s[mss[i].s].name);
					else snprintf(lbl[sym+l], 8, "%s",mct->s[mss[i].s].name);
					lbl[sym+l][strlen(lbl[sym+l])] = ' ';
				}
			}
			row += mss[i].salc[j].d;
		}
		if(mss[i].salcl) sym += mct->s[mss[i].s].d;
	}
	*nirreps = sym;

	*err = ret;
	return ret;
err:
	errstr = msymErrorString(ret);
	fprintf(stderr,"%s",errstr);
	errstr = msymGetErrorDetails();
	fprintf(stderr,": %s\n",errstr);
	*err = ret;
	return ret;
}

INT cmsym_symmetrize_orbitals(msym_context *pctx, INT *l, double c[*l][*l], INT *err){
	msym_context ctx = *pctx;
	msym_error_t ret = MSYM_SUCCESS;
	const char *errstr;
	if(MSYM_SUCCESS != (ret = msymSymmetrizeWavefunctions(ctx, *l, c, NULL, NULL))) goto err;
	*err = ret;
	return ret;
err:
	errstr = msymErrorString(ret);
	fprintf(stderr,"%s",errstr);
	errstr = msymGetErrorDetails();
	fprintf(stderr,": %s\n",errstr);
	*err = ret;
	return ret;
}
