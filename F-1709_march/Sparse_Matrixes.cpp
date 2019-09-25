#include "Sparse_Matrixes.h"

compressed_matrix::compressed_matrix ()
{
	use_CSR = false;
	size = 0;
	n_entries = 0;

	ig = jg = NULL;

	ggl = ggu = di = di_inv = NULL;
	L = U = Ld = NULL;
	f = x = x0 = NULL;

	ia = ja = NULL;
	a = NULL;
}

compressed_matrix::compressed_matrix (int N, int N_L)
{
	use_CSR = false;
	size = N;
	n_entries = N_L;

	ig = new int[size + 1];
	ggl = new MathVector (n_entries);
	ggu = new MathVector (n_entries);
	di = new MathVector (size);
	di_inv = new MathVector (size);
	f = new MathVector (size);

	jg = new int[n_entries];

	L = new MathVector (n_entries);
	U = new MathVector (n_entries);
	Ld = new MathVector (size);

	x = new MathVector (size);
	x0 = new MathVector (size);

	ia = ja = NULL;
	a = NULL;
}

compressed_matrix::compressed_matrix (const compressed_matrix & crm)
{
	use_CSR = crm.use_CSR;
	if (size != crm.size)
	{
		if (size != 0)
		{
			delete ggl;
			delete ggu;
			delete di;
			delete di_inv;

			delete L;
			delete U;
			delete Ld;

			delete f;
			delete x;
			delete x0;

			delete[] jg;
			delete[] ig;
		}
		size = crm.size;
		n_entries = crm.n_entries;
		ig = new int[size + 1];
		jg = new int[n_entries];

		ggl = new MathVector (n_entries);
		ggu = new MathVector (n_entries);
		di = new MathVector (size);
		di_inv = new MathVector (size);
		
		L = new MathVector (n_entries);
		U = new MathVector (n_entries);
		Ld = new MathVector (size);

		f = new MathVector (size);
		x = new MathVector (size);
		x0 = new MathVector (size);
	}

	ggl->Copy (*(crm.ggl));
	ggu->Copy (*(crm.ggu));
	di->Copy (*(crm.di));

	L->Copy (*(crm.L));
	U->Copy (*(crm.U));
	Ld->Copy (*(crm.Ld));

	f->Copy (*(crm.f));
	x->Copy (*(crm.x));
	x0->Copy (*(crm.x0));

	for (int i = 0; i < size + 1; i++)
	{
		ig[i] = crm.ig[i];
	}
	for (int i = 0; i < n_entries; i++)
	{
		jg[i] = crm.jg[i];
	}
}

compressed_matrix::~compressed_matrix ()
{
	delete ggl;
	delete ggu;
	delete di;
	delete di_inv;

	delete L;
	delete U;
	delete Ld;

	delete f;
	delete x;
	delete x0;

	delete[] jg;
	delete[] ig;

	if (ia != NULL)
		delete ia;
	if (ja != NULL)
		delete ja;
	if (a != NULL)
		delete a;
}

compressed_matrix & compressed_matrix::operator=(compressed_matrix & crm)
{
	if (this != &crm)
	{
		if (size != crm.size)
		{
			if (size != 0)
			{
				delete ggl;
				delete ggu;
				delete di;

				delete L;
				delete U;
				delete Ld;

				delete f;
				delete x;
				delete x0;

				delete[] jg;
				delete[] ig;
			}
			size = crm.size;
			n_entries = crm.n_entries;

			ig = new int[size + 1];
			ggl = new MathVector (n_entries);
			ggu = new MathVector (n_entries);
			di = new MathVector (size);
			f = new MathVector (size);

			jg = new int[n_entries];

			L = new MathVector (n_entries);
			U = new MathVector (n_entries);
			Ld = new MathVector (size);

			x = new MathVector (size);
			x0 = new MathVector (size);
		}

		ggl->Copy (*(crm.ggl));
		ggu->Copy (*(crm.ggu));
		di->Copy (*(crm.di));

		L->Copy (*(crm.L));
		U->Copy (*(crm.U));
		Ld->Copy (*(crm.Ld));

		f->Copy (*(crm.f));
		x->Copy (*(crm.x));
		x0->Copy (*(crm.x0));

		for (int i = 0; i < size + 1; i++)
		{
			ig[i] = crm.ig[i];
		}
		for (int i = 0; i < n_entries; i++)
		{
			jg[i] = crm.jg[i];
		}
	}
	return *this;
}

void compressed_matrix::set_matrix_size (int Size, int N_entries)
{
	if (size != 0)
	{
		delete ggl;
		delete ggu;
		delete di;

		delete L;
		delete U;
		delete Ld;

		delete f;
		delete x;
		delete x0;

		delete[] jg;
		delete[] ig;
	}
	size = Size;
	n_entries = N_entries;

	ig = new int[size + 1];
	ggl = new MathVector (n_entries);
	ggu = new MathVector (n_entries);
	di = new MathVector (size);
	f = new MathVector (size);

	jg = new int[n_entries];

	L = new MathVector (n_entries);
	U = new MathVector (n_entries);
	Ld = new MathVector (size);

	x = new MathVector (size);
	x0 = new MathVector (size);
}

int compressed_matrix::Size ()
{
	return size;
}

void compressed_matrix::factorize () // don't touch it if it works 
{
	int i, j, start, end, row, kU, kL, Uc, Lc, Ue, Le;
	double sL, sU, sD;

	Inner_copy ();

	for (i = 0; i<size; i++)
	{
		//проход по всем элементам L
		start = ig[i];
		end = ig[i + 1];
		sD = 0;
		for (j = start; j<end; j++)
		{
			sL = sU = 0;
			row = jg[j];
			
			//начало элементов строки L и столбца U
			kU = ig[row];
			kL = start;

			//окончание элементов строки L и столбца U
			Ue = ig[row + 1];
			Le = j;

			//столбец элемента L и строка элемента U
			Uc = jg[kU];
			Lc = jg[kL];

			//пока не пройдены все нужные элементы по L и U
			while (kL < Le && kU < Ue)
			{
				//если номер строки текущего элемента U меньше, чем номер столбца элемента L, перейти на след. элемент U
				if (Lc > Uc)
				{
					kU++;
					Uc = jg[kU];
				}

				//если номер строки текущего элемента U больше, чем номер столбца элемента L, перейти на след. элемент L
				if (Uc > Lc)
				{
					kL++;
					Lc = jg[kL];
				}

				if (Uc == Lc && kL < Le && kU < Ue)
				{
					sL += U->getElem (kU) * L->getElem (kL);
					sU += U->getElem (kL) * L->getElem (kU);
					kU++;
					kL++;
					Uc = jg[kU];
					Lc = jg[kL];
				}
			}
			L->setElem (j, L->getElem (j) - sL);
			U->setElem (j, (U->getElem (j) - sU) / Ld->getElem (row));
			sD += U->getElem (j) * L->getElem (j);
		}
		Ld->setElem (i, (Ld->getElem (i) - sD));
	}
}

void compressed_matrix::test_full (int Size, bool symm)
{
	size = Size;
	int elem_counter = 0;
	{
		ig = new int[size + 1];
		ig[0] = 0;
		for (int i = 1; i < size + 1; i++)
		{
			ig[i] = elem_counter;
			elem_counter += i;
		}
	}
	{
		int counter = 0;
		jg = new int[ig[size]];

		for (int i = 1; i < size; i++)
		{
			for (int j = 0; j < i; j++)
			{
				jg[counter] = j;
				counter++;
			}
		}
	}
	{
		double h = 0.1;
		L = new MathVector (ig[size]);
		U = new MathVector (ig[size]);
		ggl = new MathVector (ig[size]);
		ggu = new MathVector (ig[size]);
		for (int i = 0; i < ig[size]; i++)
		{
			ggl->setElem (i, h);
			if (symm)
				ggu->setElem (i, h);
			else
				ggu->setElem (i, sqrt (10.0 * h));

			h += 0.05;
			if (h > 1.5)
				h = 0.1;
		}
	}
	{
		di = new MathVector (size);
		Ld = new MathVector (size);
		for (int i = 0; i < size; i++)
		{
			di->setElem (i, pow ((double)(size), 0.5) + pow ((double)i, 0.1));
		}
	}
	{
		x = new MathVector (size);
		x0 = new MathVector (size);
		for (int i = 0; i < size; i++)
			x->setElem (i, (double)i);
	}
	{
		f = new MathVector (size);
		mult_A_v (*x, f);

		x0->Zero ();
		x->Zero ();
	}
}

void compressed_matrix::test_sparse (int Size, bool symm)
{
	size = Size;
	int elem_counter = 0;
	{
		ig = new int[size + 1];
		ig[0] = 0;
		for (int i = 1; i < size + 1; i++)
		{
			ig[i] = elem_counter;
			elem_counter += (i + 1) / 2;
		}
	}
	{
		int counter = 0;
		jg = new int[ig[size]];

		for (int i = 1; i < size; i++)
		{
			for (int j = (i + 1) % 2; j < i; j += 2)
			{
				jg[counter] = j;
				counter++;
			}
		}
	}
	{
		double h = 0.1;
		double hu = 0.2;
		L = new MathVector (ig[size]);
		U = new MathVector (ig[size]);
		ggl = new MathVector (ig[size]);
		ggu = new MathVector (ig[size]);
		for (int i = 0; i < ig[size]; i++)
		{
			ggl->setElem (i, h);
			if (symm)
				ggu->setElem (i, h);
			else
				ggu->setElem (i, hu);

			h += 0.2;
			if (h > 3.0)
				h = 0.1;

			hu += 0.2;
			if (hu > 4.0)
				hu = 0.2;
		}
	}
	{
		di = new MathVector (size);
		Ld = new MathVector (size);
		for (int i = 0; i < size; i++)
		{
			di->setElem (i, (double)(i + 1) * 10.0);
		}
	}
	{
		x = new MathVector (size);
		x0 = new MathVector (size);
		for (int i = 0; i < size; i++)
			x->setElem (i, (double)i);
	}
	{
		f = new MathVector (size);
		mult_A_v (*x, f);

		x0->Zero ();
		x->Zero ();
	}
}

int compressed_matrix::solve (int * param)
{
	/* param: 
	0 - method
	1 - decomp_type
	2 - depth (for gmres)
	3 - number of threads
	4 - mkl functions */

	if (param[4] == SOLVER_MKL_YES)
	{
		if (!use_CSR)
			convert_to_CSR ();
		refresh_CSR ();
	}

	int solver_iterations = 0;
	switch (param[0])
	{
	case (SOLVER_METHOD_LOS):
	{
		if (param[4] == SOLVER_MKL_YES)
			solver_iterations = solve_LOS_MKL (param[1], param[3]);
		else
			solver_iterations = solve_LOS (param[1], param[3]);
		break;
	}
	case (SOLVER_METHOD_GMRES):
	{
		solver_iterations = solve_GMRES (param[1], param[2], param[3]); // no mkl version
		break;
	}
	case (SOLVER_METHOD_CGM_SYMM):
	{
		solver_iterations = solve_CGM_symm (param[1], param[3]); // no mkl version
		break;
	}
	case (SOLVER_METHOD_PARDISO):
	{
		if (!use_CSR)
			convert_to_CSR ();
		solve_pardiso (param[3]);
		break;
	}
	case (SOLVER_METHOD_LU):
	{
		solve_LU ();
		break;
	}
	case (SOLVER_METHOD_LLT):
	{
		solve_LL ();
		break;
	}
	}
	return solver_iterations;
}

void compressed_matrix::factorize_Cholesky ()
{
	int i, j, start, end, row, kU, kL, Uc, Lc, Ue, Le;
	double sL, sD;

	Inner_copy ();

	for (i = 0; i < size; i++)
	{
		//проход по строке L
		start = ig[i];
		end = ig[i + 1];
		sD = 0;
		for (j = start; j < end; j++)
		{
			sL = 0;
			row = jg[j]; // столбец текущего элемента

			// начало элементов строки L и столбца U
			kU = ig[row];
			kL = start;

			// окончание элементов строки L и столбца U
			Ue = ig[row + 1];
			Le = j;

			//столбец элемента L и строка элемента U
			Uc = jg[kU];
			Lc = jg[kL];

			//пока не пройдены все нужные элементы по L и U
			while (kL < Le && kU < Ue)
			{
				//если номер строки текущего элемента U меньше, чем номер столбца элемента L, перейти на след. элемент U
				if (Lc > Uc)
				{
					kU++;
					Uc = jg[kU];
				}

				//если номер строки текущего элемента U больше, чем номер столбца элемента L, перейти на след. элемент L
				if (Uc > Lc)
				{
					kL++;
					Lc = jg[kL];
				}

				if (Uc == Lc && kL < Le && kU < Ue)
				{
					sL += L->getElem (kU) * L->getElem (kL);
					kU++;
					kL++;
					Uc = jg[kU];
					Lc = jg[kL];
				}
			}
			L->setElem (j, (L->getElem (j) - sL) / Ld->getElem (row));
			U->setElem (j, (U->getElem (j) - sL) / Ld->getElem (row));
			sD += L->getElem (j) * L->getElem (j);
		}
		double ld = Ld->getElem (i);
		if (Ld->getElem (i) - sD < 0)
			printf ("ERROR: sparse_matrix factorize_Cholesky\n");
		Ld->setElem (i, sqrt (Ld->getElem (i) - sD));
	}
}

void compressed_matrix::refresh_CSR ()
{
	int start, end, row;
	int n = ia [size];

	// fill c
	int * c = new int[n];
	for (int i = 0; i < size; i++)
		c[i] = 0;

	// fill a
	for (int i = 0; i < size; i++)
	{
		// lower triangle
		start = ig[i];
		end = ig[i + 1];
		for (int j = start; j < end; j++)
		{
			// lower triangle
			a[ia[i] + c[i]] = ggl->getElem (j);
			c[i]++;
		}
		// diagonal
		a[ia[i] + c[i]] = di->getElem (i);
		c[i]++;
		for (int j = start; j < end; j++)
		{
			// upper triangle
			row = jg[j];
			a[ia[row] + c[row]] = ggu->getElem (j);
			c[row]++;
		}
	}
	delete[] c;
}

void compressed_matrix::mult_U_v (const MathVector & v, MathVector * res)
{
	int i, j;
	int start, end, row;

	for (i = 0; i < size; i++)
	{
		res->setElem(i, v.getElem(i));

		start = ig[i];
		end = ig[i + 1];
		for (j = start; j < end; j++)
		{
			row = jg[j];
			//умножение на U
			res->addToElem(row,  U->getElem(j) * v.getElem(i));
		}
	}
}

void compressed_matrix::inverse_di ()
{
	if (di_inv == NULL)
		di_inv = new MathVector (size);
	for (int i = 0; i < size; i++)
		di_inv->setElem (i, 1.0 / di->getElem (i));
}

void compressed_matrix::mult_A_v (const MathVector & v, MathVector * res)
{
	int i, j;
	int start, end, row;

	for (i = 0; i < size; i++)
	{
		res->setElem (i, di->getElem (i) * v.getElem (i));

		start = ig[i];
		end = ig[i + 1];
		for (j = start; j < end; j++)
		{
			row = jg[j];
			//умножение на матрицу А
			res->addToElem (i, ggl->getElem (j) * v.getElem (row));
			res->addToElem (row, ggu->getElem (j) * v.getElem (i));
		}
	}
}

void compressed_matrix::Inner_copy ()
{
	Ld->Copy (*di);
	L->Copy (*ggl);
	U->Copy (*ggu);
}

void compressed_matrix::Clear ()
{
	ggl->Zero ();
	ggu->Zero ();
	di->Zero ();
	f->Zero ();
}

void compressed_matrix::calc_LvR (const MathVector & r, MathVector * v)
{
	int i, j, start, end, row;
	double s;

	for (i = 0; i < size; i++)
	{
		s = 0;
		start = ig[i];
		end = ig[i + 1];

		for (j = start; j < end; j++)
		{
			row = jg[j];
			s += L->getElem (j) * v->getElem (row);
		}

		v->setElem (i, (r.getElem (i)- s) / Ld->getElem(i));
	}
}

void compressed_matrix::calc_URV (const MathVector & v, MathVector * r)
{
	int i, k, start, end, row;
	double s;

	MathVector * vcopy = new MathVector (size);
	vcopy->Copy (v);

	for (i = size - 1; i >= 0; i--)
	{
		r->setElem (i, vcopy->getElem (i));
		start = ig[i];
		end = ig[i + 1];

		for (k = start; k < end; k++)
		{
			row = jg[k];
			s = U->getElem(k) * r->getElem (i);
			vcopy->addToElem (row, -s);
		}
	}
	delete vcopy;
}

void compressed_matrix::calc_LTRV (const MathVector & v, MathVector * r)
{
	int i, k, start, end, row;
	double s;

	MathVector * vcopy = new MathVector (size);
	vcopy->Copy (v);

	for (i = size - 1; i >= 0; i--)
	{
		r->setElem (i, vcopy->getElem (i) / Ld->getElem (i));
		start = ig[i];
		end = ig[i + 1];

		for (k = start; k < end; k++)
		{
			row = jg[k];
			s = L->getElem (k) * r->getElem (i);
			vcopy->addToElem (row, -s);
		}
	}
	delete vcopy;
}

void compressed_matrix::calc_DvR (const MathVector & r, MathVector * v)
{
	for (int i = 0; i < size; i++)
	{
		v->setElem (i, r.getElem (i) * di_inv->getElem (i));
	}
}

void compressed_matrix::mult_A_v (const MathVector & v, MathVector * res, double * res_omp, int par_num_threads)
{
	int i, j;
	int start, end, row;
	int rank;

	if (size > PARALLEL_MIN_SIZE && par_num_threads > 1)
	{
		P_BLAS::mult_A_v (size, ig, jg, di->elem, ggl->elem, ggu->elem, v.elem, res->elem, res_omp, par_num_threads);
	}
	else
	{
		mult_A_v (v, res);
	}
}

void compressed_matrix::calc_DvR (const MathVector & r, MathVector * v, int par_num_threads)
{
	if (par_num_threads > 1 && size > PARALLEL_MIN_SIZE)
	{
		P_BLAS::calc_DvR (size, di_inv->elem, r.elem, v->elem, par_num_threads);
	}
	else
	{
		calc_DvR (r, v);
	}
}

void compressed_matrix::mult_A_v_MKL (const sparse_matrix_t SMA, const MathVector & v, MathVector * res)
{
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	mkl_sparse_d_mv (SPARSE_OPERATION_NON_TRANSPOSE, 1.0, SMA, descr, v.elem, 0.0, res->elem);
}

void compressed_matrix::calc_DvR_MKL (const MathVector & r, MathVector * v)
{
	vdmul (&size, di_inv->elem, r.elem, v->elem);
}

void compressed_matrix::calc_LvR_MKL (const sparse_matrix_t SMA, const MathVector & r, MathVector * v)
{
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1.0, SMA, descr, r.elem, v->elem);
}

void compressed_matrix::calc_URV_MKL (const sparse_matrix_t SMA, const MathVector & v, MathVector * r)
{
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_UNIT;
	mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1.0, SMA, descr, v.elem, r->elem);
}

int compressed_matrix::solve_GMRES (int d_type, int max_depth)
{
	int k, l, j, i;

	MathVector * v2 = new MathVector (size);
	MathVector * v3 = new MathVector (size);
	MathVector * r = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * w = new MathVector (size);
	MathVector * d, *h, *z;
	double b_norm = f->Norm ();

	switch (d_type)
	{
	case 0:
		inverse_di ();
		x->Copy (*x0);
		mult_A_v (*x0, mr);
		v2->Linear_Combination (*f, 1.0, *mr, -1.0);
		calc_DvR (*v2, r);
		break;
	default:
		factorize ();
		mult_U_v (*x0, x);
		mult_A_v (*x0, mr);
		v2->Linear_Combination (*f, 1.0, *mr, -1.0);
		calc_LvR (*v2, r);
	}


	Matrix * V = new Matrix (size, max_depth + 1);

	Matrix * V_fin;
	Matrix * H = new Matrix (max_depth + 1, max_depth);
	H->IdentityMatrix ();
	Matrix * H_fin;
	MathVector * V_cur = new MathVector (size);
	Matrix * R;
	Matrix * H_H, *H_H_H;

	double c, s, coef, h1, h2;
	bool got_H, found;
	found = false;
	double cur_disc, pred_disc;
	pred_disc = 1.0;
	l = size - 1;
	for (k = 0; k < ITERMAX && !found; k++)
	{
		H->Zero ();
		V->Zero ();
		V_cur->Copy (*r);
		V_cur->MultiplyByNumber (1.0L / r->Norm ());
		V->Insert_Column (0, V_cur);

		got_H = false;
		for (l = 0; l < H->Size1 () && !got_H; l++)
		{
			switch (d_type)
			{
			case 0:
				mult_A_v (*V_cur, v3);
				calc_DvR (*v3, w);
				break;
			default:
				calc_URV (*V_cur, v2);
				mult_A_v (*v2, v3);
				calc_LvR (*v3, w);
			}
			mr->Copy (*w);

			for (j = 0; j <= l; j++)
			{
				V->Get_Column (j, v2);
				H->setElem (j, l, w->Scalar_Product (*v2));
				 
				mr->Linear_Combination (1.0, *v2, -H->Elem (j, l));
			}

			H->setElem (l + 1, l, mr->Norm ());
			if (fabs (H->Elem (l + 1, l)) < LIL)
			{
				got_H = true;
			}
			else
			{
				V_cur->Linear_Combination (0.0, *mr, 1.0L / H->Elem (l + 1, l));
				V->Insert_Column (l + 1, V_cur);
			}
		}

		if (!got_H)
			H_fin = new Matrix (*H);
		else
		{
			H_fin = new Matrix (l + 1, l);
			H_fin->Copy (*H);
		}

		V_fin = new Matrix (V->Size0 (), H->Size1 ());
		V_fin->Copy (*V);

		d = new MathVector (H_fin->Size0 ());
		h = new MathVector (H_fin->Size0 ());
		z = new MathVector (H_fin->Size0 () - 1);
		d->setElem (0, r->Norm ());

		R = new Matrix (H_fin->Size0 ());
		H_H = new Matrix (*H_fin);
		H_H_H = new Matrix (*H_fin);

		for (i = 0; i < H_fin->Size0 () - 1; i++)
		{
			h1 = H_H->Elem (i, i);
			h2 = H_fin->Elem (i + 1, i);

			coef = sqrt (h1 * h1 + h2 * h2);
			c = h1 / coef;
			s = h2 / coef;
			R->Rotation_Matrix (c, s, i);

			R->MultiplyByMatrix (*H_H, H_H_H);
			H_H->Copy (*H_H_H);

			R->MultiplyMatrixByVector (*d, h);
			d->Copy (*h);
		}

		H_H->Cut_Last_Row ();
		d->Cut_Last ();
		H_H->calc_URV (*d, z);

		V_fin->MultiplyMatrixByVector (*z, v2);

		x->Add (*v2);

		// r = S-1 (f - AQ-1x)
		switch (d_type)
		{
		case 0:
			mult_A_v (*x, v3);
			v2->Linear_Combination (*f, 1.0, *v3, -1.0);
			calc_DvR (*v2, r);
			break;
		default:
			calc_URV (*x, v2);
			mult_A_v (*v2, v3);
			v2->Linear_Combination (*f, 1.0, *v3, -1.0);
			calc_LvR (*v2, r);
		}
		cur_disc = r->Norm () / b_norm;
		if (cur_disc < SM_PRECIS || fabs (cur_disc / pred_disc - 1.0L) < 1e-05)
			found = true;
		pred_disc = cur_disc;
		//printf ("%i\t%e\n", k, cur_disc);
		delete H_fin;
		delete V_fin;
		delete R;
		delete H_H;
		delete H_H_H;
		delete d;
		delete h;
		delete z;
	}

	switch (d_type)
	{
	case 0:
		break;
	default:
		v2->Copy (*x);
		calc_URV (*v2, x);
	}
	mult_A_v (*x, v2);
	v2->Substract (*f);
	printf ("GMRES solver: iter %i, result %e\n", k, v2->Norm () / b_norm);

	delete v2;
	delete v3;
	delete r;
	delete mr;
	delete w;
	delete H;
	delete V;
	delete V_cur;

	return k;
}

int compressed_matrix::solve_LOS (int d_type)
{
	int k;
	double alpha, beta, s;
	bool found = false;
		
	MathVector * r = new MathVector (size);
	MathVector * r_smooth = new MathVector (size);
	MathVector * z = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * v2 = new MathVector (size);
	MathVector * p = new MathVector (size);
	MathVector * x_smooth = new MathVector (size);
		
	double b_norm = f->Norm ();

	switch (d_type)
	{
	case 0:
		inverse_di ();
		mult_A_v (*x0, mr);
		mr->Linear_Combination (-1.0, *f, 1.0);
		calc_DvR (*mr, r);
		z->Copy (*r);
		mult_A_v (*z, v2);
		calc_DvR (*v2, p);
		break;
	case 1:
		factorize ();
		// r = S-1 (b - Ax0)
		mult_A_v (*x0, mr);
		mr->Linear_Combination (-1.0, *f, 1.0);
		calc_LvR (*mr, r);
		// z = Q-1r
		calc_URV (*r, z);
		// p = S-1Az
		mult_A_v (*z, v2);
		calc_LvR (*v2, p);
		break;
	default:
		// r = f - Ax
		mult_A_v (*x0, mr);
		r->Linear_Combination (*mr, -1.0, *f, 1.0);
		// z = r
		z->Copy (*r);
		// p = Az
		mult_A_v (*z, p);
	}
	double rn = r->Norm ();
	if (rn < LIL)
		return 0;

 
	x->Copy (*x0);
	x_smooth->Copy (*x0);
	r_smooth->Copy (*r);

	double etta;
	for (k = 0; k < ITERMAX && !found; k++)
	{
		s = p->Scalar_Product (*p);	
		alpha = p->Scalar_Product (*r);
		alpha /= s;

		// x += alpha * z
		x->Linear_Combination (1.0, *z, alpha);

		// r -= alpha * p
		r->Linear_Combination (1.0, *p, -alpha);

		switch (d_type)
		{
		case 0:
			// mr = S-1Ar
			mult_A_v (*r, v2);
			calc_DvR (*v2, mr);
			// v2 = r
			v2->Copy (*r);
			break;
		case 1:
			// mr = S-1AQ-1r
			calc_URV (*r, mr);
			mult_A_v (*mr, v2);
			calc_LvR (*v2, mr);
			// v2 = Q-1r
			calc_URV (*r, v2);
			break;
		default:
			// mr = Ar
			mult_A_v (*r, mr);
		}

		beta = mr->Scalar_Product (*p);
		beta = -beta / s;

		z->Linear_Combination (beta, *v2, 1.0);

		p->Linear_Combination (beta, *mr, 1.0);

		{
			v2->Linear_Combination (*r, 1.0, *r_smooth, -1.0);
			etta = -r_smooth->Scalar_Product (*v2) / v2->Scalar_Product (*v2);
			if (etta > 1.0 - 1e-17)
			{
				x_smooth->Copy (*x);
				r_smooth->Copy (*r);
			}
			else
			{
				if (etta > 1e-17)
				{
					x_smooth->Linear_Combination (1.0 - etta, *x, etta);
					r_smooth->Linear_Combination (1.0 - etta, *r, etta);
				}
			}

		}

		if(r_smooth->Norm () < SM_PRECIS)
			found = true;
	}

	x->Copy (*x_smooth);
	mult_A_v (*x, z);
	z->Substract (*f);
	printf ("LOS solver: iter %i, result %e\n", k, z->Norm () / b_norm);

	delete r;
	delete r_smooth;
	delete x_smooth;
	delete z;
	delete mr;
	delete v2;
	delete p;

	return k;
}

int compressed_matrix::solve_CGM_symm (int d_type)
{
	MathVector * r = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * z = new MathVector (size);
	MathVector * v = new MathVector (size);
	MathVector * v2 = new MathVector (size);
	int k = 0;
	MathVector * x_smooth = new MathVector (size);
	MathVector * r_smooth = new MathVector (size);

	x->Copy (*x0);
	x_smooth->Copy (*x0);
	double etta;

	double b_norm = f->Norm ();

	double zero_norm = f->Norm ();
	if (zero_norm > 1e-30)
	{
		Inner_copy ();
		mult_A_v (*x0, mr);
		r->Linear_Combination (*f, 1.0, *mr, -1.0);
		r_smooth->Copy (*r);

		// z = M-1r
		switch (d_type)
		{
		case 0:
			inverse_di ();
			calc_DvR (*r, z);
			break;
		case 1:
			factorize_Cholesky ();
			calc_LvR (*r, v);
			calc_LTRV (*v, z);
			break;
		default:
			// z = r
			z->Copy (*r);
		}
		mr->Copy (*z);

		double alpha, beta, s;
		bool found = false;
		double norm;

		for (k; k < ITERMAX && !found; k++)
		{			
			s = r->Scalar_Product (*mr);

			mult_A_v (*z, v);
			alpha = s / v->Scalar_Product (*z);

			x->Linear_Combination (1.0, *z, alpha);
			r->Linear_Combination (1.0, *v, -alpha);

			switch (d_type)
			{
			case 0:
				calc_DvR (*r, mr);
				break;
			case 1:
				calc_LvR (*r, v2);
				calc_LTRV (*v2, mr);
				break;
			case 2:
				mr->Copy (*r);
				break;
			}	

			beta = r->Scalar_Product (*mr) / s;

			z->Linear_Combination (beta, *mr, 1.0);
			norm = r->Norm ();

			{
				v2->Linear_Combination (*r, 1.0, *r_smooth, -1.0);
				etta = -r_smooth->Scalar_Product (*v2) / v2->Scalar_Product (*v2);
				if (etta > 1.0 - 1e-17)
				{
					x_smooth->Copy (*x);
					r_smooth->Copy (*r);
				}
				else
				{
					if (etta > 1e-17)
					{
						x_smooth->Linear_Combination (1.0 - etta, *x, etta);
						r_smooth->Linear_Combination (1.0 - etta, *r, etta);
					}
				}

			}

			if (r_smooth->Norm () / b_norm < SM_PRECIS)
				found = true;
		}
	}

	x->Copy (*x_smooth);
	mult_A_v (*x, z);
	z->Substract (*f);
	printf ("CGM LLT solver: iter %i, result %e\n", k, z->Norm () / b_norm);
	// cleanup
	{
		delete r;
		delete mr;
		delete z;
		delete v;
		delete v2;
		delete r_smooth;
		delete x_smooth;
	}
	return k;
}

void compressed_matrix::solve_LU ()
{
	MathVector * z = new MathVector (size);
	factorize ();
	calc_LvR (*f, z);
	calc_URV (*z, x);
	delete z;
}

void compressed_matrix::solve_LL ()
{
	MathVector * z = new MathVector (size);
	factorize_Cholesky ();
	calc_LvR (*f, z);
	calc_LTRV (*z, x);
	delete z;
}

void compressed_matrix::solve_pardiso ()
{
	int start, end, row;
	// для PARDISO
	MKL_INT n = size;
	MKL_INT mtype = 11; // real and non-symmetric 
	MKL_INT nrhs = 1;
	void *pt[64];
	MKL_INT maxfct = 1;
	MKL_INT	mnum = 1;
	MKL_INT msglvl = 0;
	MKL_INT phase = 13;
	MKL_INT *perm = NULL;
	MKL_INT iparm[64];
	MKL_INT info;
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	iparm[0] = 0; //iparm(2) - iparm(64) are filled with default values.
	for (int i = 1; i < 64; i++)
	{
		iparm[i] = 0;
	}

	// solve
	perm = new MKL_INT[size];
	n = size;
	
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, f->elem, x->elem, &info);
	
	MathVector * z = new MathVector (size);
	mult_A_v (*x, z);
	z->Substract (*f);
	printf ("pardiso solver: info %i, result %e\n", info, z->Norm () / f->Norm ());

	delete[] perm;
	delete z;
}

void compressed_matrix::solve_pardiso (int par_num_threads)
{
	int start, end, row;
	// для PARDISO
	MKL_INT n = size;
	MKL_INT mtype = 11; // real and non-symmetric 
	MKL_INT nrhs = 1;
	void *pt[64];
	MKL_INT maxfct = 1;
	MKL_INT	mnum = 1;
	MKL_INT msglvl = 0;
	MKL_INT phase = 13;
	MKL_INT *perm = NULL;
	MKL_INT iparm[64];
	MKL_INT info;
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	iparm[0] = 1;
	for (int i = 1; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[1] = 2; 
	iparm[9] = 13;
	iparm[10] = 0;
	iparm[12] = 0;
	iparm[23] = 1;
	iparm[34] = 1;

	// solve
	refresh_CSR ();
	perm = new MKL_INT[size];
	n = size;

	mkl_set_num_threads (par_num_threads);
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, f->elem, x->elem, &info);

	MathVector * z = new MathVector (size);
	mult_A_v (*x, z);
	z->Substract (*f);
	double f_norm = f->Norm ();
	if (f_norm > 1e-15)
		printf ("pardiso solver: info %i, result %e\n", info, z->Norm () / f->Norm ());
	else
		printf ("pardiso solver: info %i, result %e\n", info, z->Norm ());

	phase = -1;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, f->elem, x->elem, &info);

	{
	//	FILE * pfile = fopen ("Test//par.txt", "w");

	//	fprintf (pfile, "ia:\n");
	//	for (int i = 0; i < size + 1; i++)
	//	{
	//		fprintf (pfile, "%i\n", ia[i]);
	//	}
	//	fprintf (pfile, "\n\n");

	//	fprintf (pfile, "ja:\n");
	//	for (int i = 0; i < size; i++)
	//	{
	//		start = ia[i];
	//		end = ia[i + 1];
	//		for (int j = start; j < end; j++)
	//		{
	//			fprintf (pfile, "%i\n", ja[j]);
	//		}
	//	}
	//	fprintf (pfile, "\n\n");

	//	fprintf (pfile, "a:\n");
	//	for (int i = 0; i < size; i++)
	//	{
	//		start = ia[i];
	//		end = ia[i + 1];
	//		for (int j = start; j < end; j++)
	//		{
	//			fprintf (pfile, "%.13lf\n", a[j]);
	//		}
	//	}
	//	fprintf (pfile, "\n\n");

	//	fprintf (pfile, "f:\n");
	//	for (int i = 0; i < size; i++)
	//	{
	//		fprintf (pfile, "%.13lf\n", f->getElem (i));
	//	}
	//	fprintf (pfile, "\n\n");

	//	fprintf (pfile, "x:\n");
	//	for (int i = 0; i < size; i++)
	//	{
	//		fprintf (pfile, "%.13lf\n", x->getElem (i));
	//	}
	//	fprintf (pfile, "\n\n");

	//	fclose (pfile);
	}

	delete[] perm;
	delete z;
}

int compressed_matrix::solve_LOS (int d_type, int par_num_threads)
{
	int k;
	double alpha, beta, s;
	bool found = false;

	MathVector * r = new MathVector (size);
	MathVector * r_smooth = new MathVector (size);
	MathVector * z = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * v2 = new MathVector (size);
	MathVector * p = new MathVector (size);
	MathVector * x_smooth = new MathVector (size);

	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[size * (par_num_threads - 1)];

	double b_norm = f->Norm (par_num_threads);

	switch (d_type)
	{
	case 0:
		// inverse D
		inverse_di ();
		mult_A_v (*x0, mr, res_omp, par_num_threads);
		mr->Linear_Combination (-1.0, *f, 1.0, par_num_threads);
		calc_DvR (*mr, r, par_num_threads);
		z->Copy (*r);
		mult_A_v (*z, v2, res_omp, par_num_threads);
		calc_DvR (*v2, p, par_num_threads);
		break;
	default:
		factorize ();
		// r = S-1 (b - Ax0)
		mult_A_v (*x0, mr, res_omp, par_num_threads);
		mr->Linear_Combination (-1.0, *f, 1.0, par_num_threads);
		calc_LvR (*mr, r);
		// z = Q-1r
		calc_URV (*r, z);
		// p = S-1Az
		mult_A_v (*z, v2, res_omp, par_num_threads);
		calc_LvR (*v2, p);
	}
	double rn = r->Norm (par_num_threads);
	if (rn < LIL)
		return 0;


	x->Copy (*x0);
	x_smooth->Copy (*x0);
	r_smooth->Copy (*r);

	double etta;
	for (k = 0; k < ITERMAX && !found; k++)
	{
		s = p->Scalar_Product (*p, par_num_threads);
		alpha = p->Scalar_Product (*r, par_num_threads);
		alpha /= s;

		// x += alpha * z
		x->Linear_Combination (1.0, *z, alpha, par_num_threads);
		// r -= alpha * p
		r->Linear_Combination (1.0, *p, -alpha, par_num_threads);

		switch (d_type)
		{
		case 0:
			// mr = S-1Ar
			mult_A_v (*r, v2, res_omp, par_num_threads);
			calc_DvR (*v2, mr, par_num_threads);
			// v2 = r
			v2->Copy (*r);
			break;
		default:
			// mr = S-1AQ-1r
			calc_URV (*r, mr);
			mult_A_v (*mr, v2, res_omp, par_num_threads);
			calc_LvR (*v2, mr);
			// v2 = Q-1r
			calc_URV (*r, v2);
		}

		beta = mr->Scalar_Product (*p, par_num_threads);
		beta = -beta / s;

		z->Linear_Combination (beta, *v2, 1.0, par_num_threads);
		p->Linear_Combination (beta, *mr, 1.0, par_num_threads);

		{
			v2->Linear_Combination (*r, 1.0, *r_smooth, -1.0);
			etta = -r_smooth->Scalar_Product (*v2, par_num_threads) / v2->Scalar_Product (*v2, par_num_threads);
			if (etta > 1.0 - 1e-17)
			{
				x_smooth->Copy (*x);
				r_smooth->Copy (*r);
			}
			else
			{
				if (etta > 1e-17)
				{
					x_smooth->Linear_Combination (1.0 - etta, *x, etta);
					r_smooth->Linear_Combination (1.0 - etta, *r, etta);
				}
			}

		}

		if (r_smooth->Norm (par_num_threads) < SM_PRECIS)
			found = true;
	}

	x->Copy (*x_smooth);
	mult_A_v (*x, z, res_omp, par_num_threads);
	z->Substract (*f);
	printf ("LOS solver: iter %i, result %e\n", k, z->Norm (par_num_threads) / b_norm);

	delete r;
	delete r_smooth;
	delete x_smooth;
	delete z;
	delete mr;
	delete v2;
	delete p;

	if (res_omp != NULL)
		delete[]  res_omp;

	return k;
}

int compressed_matrix::solve_GMRES (int d_type, int max_depth, int par_num_threads)
{
	int k, l, j, i;

	MathVector * v2 = new MathVector (size);
	MathVector * v3 = new MathVector (size);
	MathVector * r = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * w = new MathVector (size);
	MathVector * d, *h, *z;
	double b_norm = f->Norm (par_num_threads);

	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[size * (par_num_threads - 1)];

	switch (d_type)
	{
	case 0:
		inverse_di ();
		x->Copy (*x0);
		mult_A_v (*x0, mr, res_omp, par_num_threads);
		v2->Linear_Combination (*f, 1.0, *mr, -1.0, par_num_threads);
		calc_DvR (*v2, r, par_num_threads);
		break;
	default:
		factorize ();
		mult_U_v (*x0, x);
		mult_A_v (*x0, mr, res_omp, par_num_threads);
		v2->Linear_Combination (*f, 1.0, *mr, -1.0, par_num_threads);
		calc_LvR (*v2, r);
	}


	Matrix * V = new Matrix (size, max_depth + 1);

	Matrix * V_fin;
	Matrix * H = new Matrix (max_depth + 1, max_depth);
	H->IdentityMatrix ();
	Matrix * H_fin;
	MathVector * V_cur = new MathVector (size);
	Matrix * R;
	Matrix * H_H, *H_H_H;

	double c, s, coef, h1, h2;
	bool got_H, found;
	found = false;
	double cur_disc, pred_disc;
	pred_disc = 1.0;
	l = size - 1;
	for (k = 0; k < ITERMAX && !found; k++)
	{
		H->Zero ();
		V->Zero ();
		V_cur->Copy (*r);
		V_cur->MultiplyByNumber (1.0L / r->Norm ());
		V->Insert_Column (0, V_cur);

		got_H = false;
		for (l = 0; l < H->Size1 () && !got_H; l++)
		{
			switch (d_type)
			{
			case 0:
				mult_A_v (*V_cur, v3, res_omp, par_num_threads);
				calc_DvR (*v3, w, par_num_threads);
				break;
			default:
				calc_URV (*V_cur, v2);
				mult_A_v (*v2, v3, res_omp, par_num_threads);
				calc_LvR (*v3, w);
			}
			mr->Copy (*w);

			for (j = 0; j <= l; j++)
			{
				V->Get_Column (j, v2);
				H->setElem (j, l, w->Scalar_Product (*v2));

				mr->Linear_Combination (1.0, *v2, -H->Elem (j, l), par_num_threads);
				//v2->MultiplyByNumber (H->Elem (j, l));
				//mr->Substract (*v2);
			}

			H->setElem (l + 1, l, mr->Norm ());
			if (fabs (H->Elem (l + 1, l)) < LIL)
			{
				got_H = true;
			}
			else
			{
				V_cur->Linear_Combination (0.0, *mr, 1.0L / H->Elem (l + 1, l));
				//mr->MultiplyByNumber (1.0L / H->Elem (l + 1, l));
				//V_cur->Copy (*mr);
				V->Insert_Column (l + 1, V_cur);
			}
		}

		if (!got_H)
			H_fin = new Matrix (*H);
		else
		{
			H_fin = new Matrix (l + 1, l);
			H_fin->Copy (*H);
		}

		V_fin = new Matrix (V->Size0 (), H->Size1 ());
		V_fin->Copy (*V);

		d = new MathVector (H_fin->Size0 ());
		h = new MathVector (H_fin->Size0 ());
		z = new MathVector (H_fin->Size0 () - 1);
		d->setElem (0, r->Norm ());

		R = new Matrix (H_fin->Size0 ());
		H_H = new Matrix (*H_fin);
		H_H_H = new Matrix (*H_fin);

		for (i = 0; i < H_fin->Size0 () - 1; i++)
		{
			h1 = H_H->Elem (i, i);
			h2 = H_fin->Elem (i + 1, i);

			coef = sqrt (h1 * h1 + h2 * h2);
			c = h1 / coef;
			s = h2 / coef;
			R->Rotation_Matrix (c, s, i);

			R->MultiplyByMatrix (*H_H, H_H_H);
			H_H->Copy (*H_H_H);

			R->MultiplyMatrixByVector (*d, h);
			d->Copy (*h);
		}

		H_H->Cut_Last_Row ();
		d->Cut_Last ();
		H_H->calc_URV (*d, z);

		V_fin->MultiplyMatrixByVector (*z, v2);

		x->Add (*v2);

		// r = S-1 (f - AQ-1x)
		switch (d_type)
		{
		case 0:
			mult_A_v (*x, v3, res_omp, par_num_threads);
			v2->Linear_Combination (*f, 1.0, *v3, -1.0, par_num_threads);
			//v2->Copy (*f);
			//v2->Substract (*v3);
			calc_DvR (*v2, r, par_num_threads);
			break;
		default:
			calc_URV (*x, v2);
			mult_A_v (*v2, v3, res_omp, par_num_threads);
			v2->Linear_Combination (*f, 1.0, *v3, -1.0, par_num_threads);
			//v2->Copy (*f);
			//v2->Substract (*v3);
			calc_LvR (*v2, r);
		}
		cur_disc = r->Norm (par_num_threads) / b_norm;
		//printf ("%i %e\n", k, cur_disc);
		if (cur_disc < SM_PRECIS || fabs (cur_disc / pred_disc - 1.0L) < 1e-05)
			found = true;
		pred_disc = cur_disc;
		//printf ("%i\t%e\n", k, cur_disc);
		delete H_fin;
		delete V_fin;
		delete R;
		delete H_H;
		delete H_H_H;
		delete d;
		delete h;
		delete z;
	}

	switch (d_type)
	{
	case 0:
		break;
	default:
		v2->Copy (*x);
		calc_URV (*v2, x);
	}
	mult_A_v (*x, v2, res_omp, par_num_threads);
	v2->Substract (*f);
	printf ("GMRES solver: iter %i, result %e\n", k, v2->Norm (par_num_threads) / b_norm);

	delete v2;
	delete v3;
	delete r;
	delete mr;
	delete w;
	delete H;
	delete V;
	delete V_cur;

	if (res_omp != NULL)
		delete[]  res_omp;

	return k;
}

int compressed_matrix::solve_CGM_symm (int d_type, int par_num_threads)
{
	MathVector * r = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * z = new MathVector (size);
	MathVector * v = new MathVector (size);
	MathVector * v2 = new MathVector (size);
	int k = 0;
	MathVector * x_smooth = new MathVector (size);
	MathVector * r_smooth = new MathVector (size);

	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[size * (par_num_threads - 1)];

	x->Copy (*x0);
	x_smooth->Copy (*x0);
	double etta;

	double b_norm = f->Norm (par_num_threads);

	double zero_norm = f->Norm (par_num_threads);
	if (zero_norm > 1e-30)
	{
		Inner_copy ();
		mult_A_v (*x0, mr, res_omp, par_num_threads);
		r->Linear_Combination (*f, 1.0, *mr, -1.0, par_num_threads);
		r_smooth->Copy (*r);

		// z = M-1r
		switch (d_type)
		{
		case 0:
			inverse_di ();
			calc_DvR (*r, z, par_num_threads);
			break;
		default:
			factorize_Cholesky ();
			calc_LvR (*r, v);
			calc_LTRV (*v, z);

		}
		mr->Copy (*z);

		double alpha, beta, s;
		bool found = false;
		double norm;

		for (k; k < ITERMAX && !found; k++)
		{
			s = r->Scalar_Product (*mr, par_num_threads);

			mult_A_v (*z, v, res_omp, par_num_threads);
			alpha = s / v->Scalar_Product (*z, par_num_threads);

			x->Linear_Combination (1.0, *z, alpha, par_num_threads);
			r->Linear_Combination (1.0, *v, -alpha, par_num_threads);

			switch (d_type)
			{
			case 0:
				calc_DvR (*r, mr, par_num_threads);
				break;
			default:
				calc_LvR (*r, v2);
				calc_LTRV (*v2, mr);
			}

			beta = r->Scalar_Product (*mr, par_num_threads) / s;

			z->Linear_Combination (beta, *mr, 1.0, par_num_threads);
			norm = r->Norm (par_num_threads);

			{
				v2->Linear_Combination (*r, 1.0, *r_smooth, -1.0, par_num_threads);
				etta = -r_smooth->Scalar_Product (*v2, par_num_threads) / v2->Scalar_Product (*v2, par_num_threads);
				if (etta > 1.0 - 1e-17)
				{
					x_smooth->Copy (*x);
					r_smooth->Copy (*r);
				}
				else
				{
					if (etta > 1e-17)
					{
						x_smooth->Linear_Combination (1.0 - etta, *x, etta, par_num_threads);
						r_smooth->Linear_Combination (1.0 - etta, *r, etta, par_num_threads);
					}
				}

			}

			if (r_smooth->Norm (par_num_threads) / b_norm < SM_PRECIS)
				found = true;
		}
	}

	x->Copy (*x_smooth);
	mult_A_v (*x, z, res_omp, par_num_threads);
	z->Substract (*f);
	printf ("CGM LLT solver: iter %i, result %e\n", k, z->Norm (par_num_threads) / b_norm);
	// cleanup
	{
		delete r;
		delete mr;
		delete z;
		delete v;
		delete v2;
		delete r_smooth;
		delete x_smooth;
		delete[] res_omp;
	}
	return k;
}

void compressed_matrix::get_from_files ()
{
	double e;

	FILE * file;

	fopen_s (&file, "matrix//param.txt", "r");
	fscanf_s (file, "%d", &size);
	fclose (file);

	fopen_s (&file, "matrix//ig.txt", "r");
	ig = new int[size + 1];
	for (int i = 0; i < size + 1; i++)
	{
		fscanf_s (file, "%d", &ig[i]);
	}
	fclose (file);

	n_entries = ig[size];

	fopen_s (&file, "matrix//ggl.txt", "r");
	ggl = new MathVector (n_entries);
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%lf", &e);
		ggl->setElem (i, e);
	}
	fclose (file);

	fopen_s (&file, "matrix//ggu.txt", "r");
	ggu = new MathVector (n_entries);
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%lf", &e);
		ggu->setElem (i, e);
	}
	fclose (file);

	fopen_s (&file, "matrix//di.txt", "r");
	di = new MathVector (size);
	for (int i = 0; i < di->getSize (); i++)
	{
		fscanf_s (file, "%lf", &e);
		di->setElem (i, e);
	}
	fclose (file);


	fopen_s (&file, "matrix//jg.txt", "r");
	jg = new int[n_entries];
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%d", &jg[i]);
	}
	fclose (file);

	fopen_s (&file, "matrix//f.txt", "r");
	f = new MathVector (size);
	for (int i = 0; i < size; i++)
	{
		fscanf_s (file, "%lf", &e);
		f->setElem (i, e);
	}
	fclose (file);

	L = new MathVector (n_entries);
	U = new MathVector (n_entries);
	Ld = new MathVector (size);

	x = new MathVector (size);
	x0 = new MathVector (size);
}

void compressed_matrix::get_from_files (char * folder)
{
	char file_name[256];
	double e;

	FILE * file;

	strcpy (file_name, folder);
	strcat (file_name, "//param.txt");
	
	fopen_s (&file, file_name, "r");
	fscanf_s (file, "%d", &size);
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//ig.txt");

	fopen_s (&file, file_name, "r");
	ig = new int[size + 1];
	for (int i = 0; i < size + 1; i++)
	{
		fscanf_s (file, "%d", &ig[i]);
	}
	fclose (file);

	n_entries = ig[size];

	strcpy (file_name, folder);
	strcat (file_name, "//ggl.txt");

	fopen_s (&file, file_name, "r");
	ggl = new MathVector (n_entries);
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%lf", &e);
		ggl->setElem (i, e);
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//ggu.txt");

	fopen_s (&file, file_name, "r");
	ggu = new MathVector (n_entries);
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%lf", &e);
		ggu->setElem (i, e);
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//di.txt");

	fopen_s (&file, file_name, "r");
	di = new MathVector (size);
	for (int i = 0; i < di->getSize (); i++)
	{
		fscanf_s (file, "%lf", &e);
		di->setElem (i, e);
	}
	fclose (file);


	strcpy (file_name, folder);
	strcat (file_name, "//jg.txt");

	fopen_s (&file, file_name, "r");
	jg = new int[n_entries];
	for (int i = 0; i < n_entries; i++)
	{
		fscanf_s (file, "%d", &jg[i]);
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//f.txt"); 

	fopen_s (&file, file_name, "r");
	f = new MathVector (size);
	for (int i = 0; i < size; i++)
	{
		fscanf_s (file, "%lf", &e);
		f->setElem (i, e);
	}
	fclose (file);

	L = new MathVector (n_entries);
	U = new MathVector (n_entries);
	Ld = new MathVector (size);

	x = new MathVector (size);
	x0 = new MathVector (size);
}

void compressed_matrix::set_x0 (double value)
{
	x0->Initialize_const (value);
}

bool compressed_matrix::set_ig_jg (int * Ig, int * Jg)
{
	for (int i = 0; i < size + 1; i++)
		ig[i] = Ig[i];
	for (int i = 0; i < n_entries; i++)
		jg[i] = Jg[i];
	if (n_entries != Ig[size])
		return false;
	return true;
}

void compressed_matrix::convert_to_CSR ()
{
	use_CSR = true;
	int start, end, row;

	int * addr = new int[size];
	for (int i = 0; i < size; i++)
	{
		addr[i] = 0;
	}
	for (int i = 0; i < size; i++)
	{
		addr[i] += ig[i + 1] - ig[i]; // lower triangle
		start = ig[i];
		end = ig[i + 1];
		for (int j = start; j < end; j++)
		{
			row = jg[j];
			addr[row]++; // upper triangle
		}
		addr[i]++; // diagonal
	}
	// amount of all entries
	int n = 0;
	for (int i = 0; i < size; i++)
	{
		n += addr[i];
	}

	if (ia != NULL)
		delete[] ia;
	ia = new MKL_INT[size + 1];
	if (ja != NULL)
		delete[] ja;
	ja = new MKL_INT[n];
	if (a != NULL)
		delete[] a;
	a = new double[n];
	int * c = new int[n];

	// fill ia
	ia[0] = 0;
	for (int i = 1; i < size + 1; i++)
		ia[i] = ia[i - 1] + addr[i - 1];

	// fill c
	for (int i = 0; i < size; i++)
		c[i] = 0;

	// fill ja and a
	for (int i = 0; i < size; i++)
	{
		// lower triangle
		start = ig[i];
		end = ig[i + 1];
		for (int j = start; j < end; j++)
		{
			// lower triangle
			a[ia[i] + c[i]] = ggl->getElem (j);
			ja[ia[i] + c[i]] = jg[j];
			c[i]++;
		}
		// diagonal
		a[ia[i] + c[i]] = di->getElem (i);
		ja[ia[i] + c[i]] = i;
		c[i]++;
		for (int j = start; j < end; j++)
		{
			// upper triangle
			row = jg[j];
			a[ia[row] + c[row]] = ggu->getElem (j);
			ja[ia[row] + c[row]] = i;
			c[row]++;
		}
	}
	delete[] c;
	delete[] addr;
}

void compressed_matrix::make_lut (int *Lia, int *Lja, double * Lt, int *Uia, int *Uja, double * Ut)
{
	int Lcounter = 0;
	int Ucounter = 0;
	for (int i = 0; i < size; i++)
	{
		Lia[i] = Lcounter;
		Uia[i] = Ucounter;
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			Lt[Lcounter] = L->getElem (j);
			Lja[Lcounter] = jg[j];
			Lcounter++;

			Ut[Ucounter] = U->getElem (j);
			Uja[Ucounter] = jg[j];
			Ucounter++;
		}
		Lt[Lcounter] = Ld->getElem (i);
		Lja[Lcounter] = i;
		Lcounter++;
	}
	Lia[size] = Lcounter;
	Uia[size] = Ucounter;
}

int compressed_matrix::solve_LOS_MKL (int d_type, int par_num_threads)
{
	int k;
	double alpha, beta, s;
	bool found = false;

	MathVector * r = new MathVector (size);
	MathVector * r_smooth = new MathVector (size);
	MathVector * z = new MathVector (size);
	MathVector * mr = new MathVector (size);
	MathVector * v2 = new MathVector (size);
	MathVector * p = new MathVector (size);
	MathVector * x_smooth = new MathVector (size);

	double b_norm = f->Norm ();
	double *Lt, *Ut;
	int * Lia, *Lja, *Uia, *Uja;
	Lt = Ut = NULL;
	Lia = Lja = Uia = Uja = NULL;

	mkl_set_num_threads (par_num_threads);
	sparse_matrix_t SMA;
	sparse_status_t stat = mkl_sparse_d_create_csr (&SMA, SPARSE_INDEX_BASE_ZERO, size, size, ia, &ia[1], ja, a);
	sparse_matrix_t SML = NULL;
	sparse_matrix_t SMU = NULL;

	switch (d_type)
	{
	case 0:
		inverse_di ();
		mult_A_v_MKL (SMA, *x0, mr);
		mr->Linear_Combination_MKL (*f, -1.0, 1.0);
		calc_DvR_MKL (*mr, r);
		z->Copy (*r);
		mult_A_v_MKL (SMA, *z, v2);
		calc_DvR_MKL (*v2, p);
		break;
	default:
		factorize ();
		// u is stored as CRC
		Lt = new double[ggl->getSize () + size];
		Lia = new int[ggl->getSize () + size + 1];
		Lja = new int[ggl->getSize () + size];
		Ut = new double[ggu->getSize ()];
		Uia = new int[ggu->getSize () + 1];
		Uja = new int[ggu->getSize ()];

		make_lut (Lia, Lja, Lt, Uia, Uja, Ut);
		mkl_sparse_d_create_csr (&SML, SPARSE_INDEX_BASE_ZERO, size, size, Lia, &Lia[1], Lja, Lt);
		mkl_sparse_d_create_csc (&SMU, SPARSE_INDEX_BASE_ZERO, size, size, Uia, &Uia[1], Uja, Ut);

		// r = S-1 (b - Ax0)
		mult_A_v_MKL (SMA, *x0, mr);
		mr->Linear_Combination_MKL (*f, -1.0, 1.0);
		//calc_LvR (*mr, r);
		calc_LvR_MKL (SML, *mr, r);
		// z = Q-1r
		calc_URV_MKL (SMU, *r, z);
		//calc_URV (*r, z);
		// p = S-1Az
		mult_A_v_MKL (SMA, *z, v2);
		calc_LvR_MKL (SML, *v2, p);
		//calc_LvR (*v2, p);
	}
	double rn = r->Norm ();
	if (rn < LIL)
		return 0;

	x->Copy (*x0);
	x_smooth->Copy (*x0);
	r_smooth->Copy (*r);

	double etta;
	for (k = 0; k < ITERMAX && !found; k++)
	{
		s = p->Scalar_Product_MKL (*p);
		alpha = p->Scalar_Product_MKL (*r);
		alpha /= s;

		// x += alpha * z
		x->Linear_Combination_MKL (*z, 1.0, alpha);;

		// r -= alpha * p
		r->Linear_Combination_MKL (*p, 1.0, -alpha);

		switch (d_type)
		{
		case 0:
			// mr = S-1Ar
			mult_A_v_MKL (SMA, *r, v2);
			calc_DvR_MKL (*v2, mr);
			// v2 = r
			v2->Copy (*r);
			break;
		default:
			// mr = S-1AQ-1r
			calc_URV_MKL (SMU, *r, mr);
			//calc_URV (*r, mr);
			mult_A_v_MKL (SMA, *mr, v2);
			calc_LvR_MKL (SML, *v2, mr);
			//calc_LvR (*v2, mr);
			// v2 = Q-1r
			//calc_URV (*r, v2);
			calc_URV_MKL (SMU, *r, v2);
		}

		beta = mr->Scalar_Product_MKL (*p);
		beta = -beta / s;

		z->Linear_Combination_MKL (*v2, beta, 1.0);

		p->Linear_Combination_MKL (*mr, beta, 1.0);

		{
			v2->Copy (*r);
			v2->Substract (*r_smooth);
			etta = -r_smooth->Scalar_Product_MKL (*v2) / v2->Scalar_Product_MKL (*v2);
			if (etta > 1.0 - 1e-17)
			{
				x_smooth->Copy (*x);
				r_smooth->Copy (*r);
			}
			else
			{
				if (etta > 1e-17)
				{
					x_smooth->Linear_Combination (1.0 - etta, *x, etta);
					r_smooth->Linear_Combination (1.0 - etta, *r, etta);
				}
			}

		}

		if (r_smooth->Norm_MKL () < SM_PRECIS)
			found = true;
	}

	x->Copy (*x_smooth);
	mult_A_v_MKL (SMA, *x, z);
	z->Substract (*f);
	printf ("LOS solver: iter %i, result %e\n", k, z->Norm () / b_norm);
	mkl_sparse_destroy (SMA);
	mkl_sparse_destroy (SML);
	mkl_sparse_destroy (SMU);

	delete r;
	delete r_smooth;
	delete x_smooth;
	delete z;
	delete mr;
	delete v2;
	delete p;
	
	if (Lt != NULL)
		delete[] Lt;
	if (Ut != NULL)
		delete[] Ut;
	if (Lia != NULL)
		delete[] Lia;
	if (Uia != NULL)
		delete[] Uia;
	if (Lja != NULL)
		delete[] Lja;
	if (Uja != NULL)
		delete[] Uja;

	return k;
}

void compressed_matrix::print ()
{
	printf ("di:\n");
	for (int i = 0; i < di->getSize (); i++)
		printf ("%lf\n", di->getElem (i));
	
	printf ("\n");
	printf ("L\tU:\n");
	for (int i = 0; i < ggl->getSize (); i++)
		printf ("%lf\t%lf\n", ggl->getElem (i), ggu->getElem (i));

	printf ("x:\n");
	for (int i = 0; i < x->getSize (); i++)
		printf ("%lf ", x->getElem (i));

	printf ("ig:\n");
	for (int i = 0; i < size + 1; i++)
		printf ("%i\n", ig[i]);

	printf ("jg:\n");
	for (int i = 0; i < ig[size]; i++)
		printf ("%i\n", jg[i]);
}

void compressed_matrix::fprint ()
{
	//printf ("di:\n");
	//for (int i = 0; i < di->getSize (); i++)
	//	printf ("%lf\n", di->getElem (i));
	//
	//printf ("\n");
	//printf ("L\tU:\n");
	//for (int i = 0; i < ggl->getSize (); i++)
	//	printf ("%lf\t%lf\n", ggl->getElem (i), ggu->getElem (i));

	//printf ("x:\n");
	//for (int i = 0; i < x->getSize (); i++)
	//	printf ("%lf\n", x->getElem (i));

	FILE * file = fopen ("matrix.txt", "w");
	fprintf (file, "ig:\n");
	for (int i = 0; i < size + 1; i++)
		fprintf (file, "%i\n", ig[i]);

	fprintf (file, "\n\n");

	fprintf (file, "jg:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%i\n", jg[i]);

	fprintf (file, "di:\n");
	for (int i = 0; i < size; i++)
		fprintf (file, "%.16lf\n", Ld->getElem(i));
	fprintf (file, "ggl:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%.16lf\n", L->getElem (i));
	fprintf (file, "ggu:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%.16lf\n", U->getElem (i));
	fprintf (file, "f:\n");
	for (int i = 0; i < size; i++)
		fprintf (file, "%.16lf\n", f->getElem (i));
	fclose (file);
}

void compressed_matrix::fprint (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	fprintf (file, "ig:\n");
	for (int i = 0; i < size + 1; i++)
		fprintf (file, "%i\n", ig[i]);

	fprintf (file, "\n\n");

	fprintf (file, "jg:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%i\n", jg[i]);

	fprintf (file, "di:\n");
	for (int i = 0; i < size; i++)
		fprintf (file, "%.16lf\n", di->getElem (i));
	fprintf (file, "ggl:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%.16lf\n", ggl->getElem (i));
	fprintf (file, "ggu:\n");
	for (int i = 0; i < ig[size]; i++)
		fprintf (file, "%.16lf\n", ggu->getElem (i));
	fprintf (file, "f:\n");
	for (int i = 0; i < size; i++)
		fprintf (file, "%.16lf\n", f->getElem (i));
	fclose (file);
}

void compressed_matrix::fprint_dif_files (char * folder)
{
	char file_name[256];
	double e;

	FILE * file;

	strcpy (file_name, folder);
	strcat (file_name, "//param.txt");

	fopen_s (&file, file_name, "w");
	fprintf (file, "%d", size);
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//ig.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < size + 1; i++)
	{
		fprintf (file, "%i\n", ig[i]);
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//ggl.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < ig[size]; i++)
	{
		fprintf (file, "%.16lf\n", ggl->getElem (i));
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//ggu.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < ig[size]; i++)
	{
		fprintf (file, "%.16lf\n", ggu->getElem (i));
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//di.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < di->getSize (); i++)
	{
		fprintf (file, "%.16lf\n", di->getElem (i));
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//jg.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < ig[size]; i++)
	{
		fprintf (file, "%i\n", jg[i]);
	}
	fclose (file);

	strcpy (file_name, folder);
	strcat (file_name, "//f.txt");

	fopen_s (&file, file_name, "w");
	for (int i = 0; i < size; i++)
	{
		fprintf (file, "%.16lf\n", f->getElem (i));
	}
	fclose (file);
}

void compressed_matrix::add_to_entry (int I, int J, double value)
{
	CSRC_add_to_entry (I, J, value);
	//if (use_CSR)
	//	CSR_add_to_entry (I, J, value);
}

void compressed_matrix::set_entry (int I, int J, double value)
{
	CSRC_set_entry (I, J, value);
	//if (use_CSR)
	//	CSR_set_entry (I, J, value);
}

void compressed_matrix::CSRC_add_to_entry (int I, int J, double value)
{
	int start, end; // section of entries on I-row / J-column
	int place; // number of entry in sparse value
	bool found;

	if (I == J) // diagonal entry
	{
		Ld->addToElem (I, value);
		di->addToElem (I, value);
	}
	else
	{
		found = false;
		place = -1;

		if (I > J) // L-entry
		{
			start = ig[I];
			end = ig[I + 1];
			
			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == J)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				L->addToElem (place, value);
				ggl->addToElem (place, value);
			}
		}
		else // U-entry
		{
			start = ig[J];
			end = ig[J + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == I)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				U->addToElem (place, value);
				ggu->addToElem (place, value);
			}
		}
	}
}

void compressed_matrix::CSRC_set_entry (int I, int J, double value)
{
	int start, end; // section of entries on I-row / J-column
	int place; // number of entry in sparse value
	bool found;

	if (I == J) // diagonal entry
	{
		di->setElem (I, value);
	}
	else
	{
		found = false;
		place = -1;

		if (I > J) // L-entry
		{
			start = ig[I];
			end = ig[I + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == J)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				ggl->setElem (place, value);
			}
		}
		else // U-entry
		{
			start = ig[J];
			end = ig[J + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == I)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				ggu->setElem (place, value);
			}
		}
	}
}

void compressed_matrix::add_to_f_entry (int I, double value)
{
	f->addToElem (I, value);
}

void compressed_matrix::set_f_entry (int I, double value)
{
	f->setElem (I, value);
}

void compressed_matrix::multiply_entry (int I, int J, double value)
{
	CSRC_multiply_entry (I, J, value);
	if (use_CSR)
		CSR_multiply_entry (I, J, value);
}

double compressed_matrix::get_entry (int I, int J)
{
	if (use_CSR)
		return CSR_get_entry (I, J);
	else
		return CSRC_get_entry (I, J);
}

void compressed_matrix::CSRC_multiply_entry (int I, int J, double value)
{
	int start, end; // section of entries on I-row / J-column
	int place; // number of entry in sparse value
	bool found;

	if (I == J) // diagonal entry
	{
		Ld->multiplyElem (I, value);
		di->multiplyElem (I, value);
	}
	else
	{
		found = false;
		place = -1;

		if (I > J) // L-entry
		{
			start = ig[I];
			end = ig[I + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == J)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				L->multiplyElem (place, value);
				ggl->multiplyElem (place, value);
			}
		}
		else // U-entry
		{
			start = ig[J];
			end = ig[J + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == I)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				U->multiplyElem (place, value);
				ggu->multiplyElem (place, value);
			}
		}
	}
}

double compressed_matrix::CSRC_get_entry (int I, int J)
{
	int start, end; // section of entries on I-row / J-column
	int place; // number of entry in sparse value
	bool found;
	double value = 0.0;

	if (I == J) // diagonal entry
	{
		value = Ld->getElem (I);
		value = di->getElem (I);
	}
	else
	{
		found = false;
		place = -1;

		if (I > J) // L-entry
		{
			start = ig[I];
			end = ig[I + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == J)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				value = L->getElem (place);
				value = ggl->getElem (place);
			}
		}
		else // U-entry
		{
			start = ig[J];
			end = ig[J + 1];

			for (int i = start; i < end && !found; i++)
			{
				if (jg[i] == I)
				{
					found = true;
					place = i;
				}
			}
			if (found)
			{
				value = U->getElem (place);
				value = ggu->getElem (place);
			}
		}
	}
	return value;
}

void compressed_matrix::CSR_add_to_entry (int I, int J, double value)
{
	int start, end;
	start = ia[I];
	end = ia[I + 1];

	bool found = false;
	int place;
	for (int i = start; i < end && !found; i++)
	{
		if (ja[i] == J)
		{
			found = true;
			place = i;
		}
	}
	if (found)
		a[place] += value;
}

void compressed_matrix::CSR_set_entry (int I, int J, double value)
{
	int start, end;
	start = ia[I];
	end = ia[I + 1];

	bool found = false;
	int place;
	for (int i = start; i < end && !found; i++)
	{
		if (ja[i] == J)
		{
			found = true;
			place = i;
		}
	}
	if (found)
		a[place] = value;
}

void compressed_matrix::CSR_multiply_entry (int I, int J, double value)
{
	int start, end;
	start = ia[I];
	end = ia[I + 1];

	bool found = false;
	int place;
	for (int i = start; i < end && !found; i++)
	{
		if (ja[i] == J)
		{
			found = true;
			place = i;
		}
	}
	if (found)
		a[place] *= value;
}

double compressed_matrix::CSR_get_entry (int I, int J)
{
	int start, end;
	start = ia[I];
	end = ia[I + 1];

	bool found = false;
	int place;
	for (int i = start; i < end && !found; i++)
	{
		if (ja[i] == J)
		{
			found = true;
			place = i;
		}
	}
	if (found)
		return a[place];

	return 0.0;
}

void compressed_matrix::get_solution (MathVector * solution)
{
	solution->Copy (*x);
	//for (int i = 0, e = (x->getSize () > solution->getSize () ? solution->getSize () : x->getSize ()); i < e; i++)
	//{
	//	solution->setElem (i, x->getElem(i));
	//}
}

void compressed_matrix::set_starting_point (MathVector * start_x)
{
	x0->Copy (*start_x);
	//for (int i = 0, e = (x0->getSize () > start_x->getSize () ? start_x->getSize () : x0->getSize ()); i < e; i++)
	//{
	//	x0->setElem (i, start_x->getElem (i));
	//}
}

void compressed_matrix::clear_row (int iF)
{
	// set ggl entries from ig[iF] to ig[iF + 1] to 0
	for (int j = ig[iF]; j < ig[iF + 1]; j++)
	{
		ggl->setElem(j, 0.0);
	}
	// go through all entires 
	for (int k = 0; k < n_entries; k++)
	{
		// if their row equals iF, set ggu to 0
		if (jg[k] == iF)
			ggu->setElem (k, 0.0);
	}
}