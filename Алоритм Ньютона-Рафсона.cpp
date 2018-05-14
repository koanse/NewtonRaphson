#include <stdio.h>
#include <conio.h>
#include <math.h>

/*float Z[50] =
{
	88, 86, 85, 83, 80, 75, 76, 89, 86, 85, 86, 83, 75, 70, 84, 91, 74, 86, 83, 81, 75,
	90, 92, 91, 92, 85, 83, 81, 91, 95, 96, 93, 85, 87, 70, 94, 95, 99, 89, 93, 85, 73,
	86, 84, 95, 94, 97, 91, 81, 75
};*/

float Z[50] =
{
	35, 26, 21,  3, 49, 17, 81, 66, 49, 53, 60, 64, 82, 88, 97, 97, 97, 93, 90, 86, 93, 92, 81, 75, 70,
	67, 68, 34, 17, 21, 18, 15, 15, 11, 29, 58, 24, 15, 14, 80, 78, 60, 40, 17, 31, 77, 97, 31, 84, 85
};



float* cov(float *Z, int n, int m);				// Расчет ковариаций для лагов 0,...,m
												// по ряду Z длиной n
float* R(float *C, int m);						// Расчет коэффициентов корреляции по
												// коэффициентам ковариации для лагов 0,...,m
void Fi(float *r, int m, float **Fi, int *p);	// Расчет параметров Ф.. для АК(p)
float* dif(float *Z, int n);					// Расчет ряда первых разностей
float* func(float *W, int n, float *Fi, int p);	// Расчет ряда для нахождения начальных оценок
												// скользящего среднего (длина ряда - результата n - p)
float* Cm(float *C, int q, float *Fi, int p);	// Автоковариации для нахождения начальных оценок СС
float* quadpr(float *Cm, int p, int q, float epsilon);	// Алгоритм Ньютона-Рафсона
main()
{
	int p, q = 10, n = 50;
	float *W = NULL, *Wm = NULL;
    float *C1 = NULL, *C2 = NULL, *C3 = NULL;
	float *r = NULL;
	float *fi = NULL;

	float C[11] = {95, 87, 84, 76, 61, 55, 42, 31, 21, 17, 8};
	
	C1 = cov(Z, n, q);
	r = R(C1, q);
	Fi(r, q, &fi, &p);

	//W = dif(Z, n);
	//W = dif(W, n - 1);
	//W = dif(W, n - 2);
	//W = dif(W, n - 3);
	//W = dif(W, n - 4);
	//W = dif(W, n - 5);
	//Wm = func(W, n - 6, fi, p);

	//printf("Wm:\n");
	//for(int i = 0; i < n - p - 6; i++)
	//	printf("%d\t%f\n", i, Wm[i]);

	//C2 = cov(Wm, n - p - 6, q);
	//C3 = Cm(C2, q, fi, p);
	//quadpr(C3, p, q, 0.01);

	C1 = Cm(C, q, fi, p);
	quadpr(C1, p, q, 0.01);

	return 0;
}

float* cov(float *Z, int n, int m)
{
	int i, j;
	double av1, av2;
	float *C = new float[m + 1];

	for(i = 0; i <= m; i++)
	{
		av1 = 0;
		for(j = i; j < n; j++)
			av1 += Z[j];
		av1 /= n - i;

		av2 = 0;
		for(j = 0; j < n - i; j++)
			av2 += Z[j];
		av2 /= n - i;

		C[i] = 0;
		for(j = i; j < n; j++)
			C[i] += (Z[j] - av1) * (Z[j - i] - av2);
		C[i] /= n - i;
	}
	
	printf("Covariations:\n");
	for(i = 0; i <= m; i++)
		printf("%d\t%f\n", i, C[i]);

	return C;
}

float* R(float *C, int m)
{
	int i;
	float *R = new float[m + 1];
	for(i = 0; i <= m; i++)
		R[i] = C[i] / C[0];
	
	return R;
}

void Fi(float *r, int m, float **Fi, int *p)
{
	int i, j;
	float sum1, sum2, e, max;

	float **MFi = new float*[m];
	for(i = 0; i < m; i++)
		MFi[i] = new float[i];
    
	for(i = 1; i <= m; i++)
	{
		sum1 = 0;
		sum2 = 0;
		for(j = 1; j < i; j++)
		{
			sum1 += MFi[i - 1 - 1][j - 1] * r[i - j];
			sum2 += MFi[i - 1 - 1][j - 1] * r[j];
		}

		MFi[i - 1][i - 1] = (r[i] - sum1) / (1 - sum2);

		for(j = 1; j < i; j++)
			MFi[i - 1][j - 1] = MFi[i - 1 - 1][j - 1] - MFi[i - 1][i - 1] *
				MFi[i - 1 - 1][i - j - 1];
	}

	e = 1 / powf(m, 0.5);
	max = 0;
    for(i = 0; i < m; i++)
		if(MFi[i][i] > max)
		{
			max = MFi[i][i];
			j = i;
		}

	if(max < e)
	{
		printf("Error\n");
		getch();
	}
		
	*p = j + 1;
	*Fi = new float[j + 1];
	for(i = 0; i <= j; i++)
		(*Fi)[i] = MFi[j][i];
	
	for(i = 0; i < m; i++)
	{
		printf("%.4f", MFi[i][i]);
		for(j = 0; j <= i; j++)
			printf("\t%.4f", MFi[i][j]);
		printf("\n");
	}
	getch();

	delete MFi;
	return;
}
float* dif(float *Z, int n)
{
	int i;
	float *W = new float[n - 1];
	for(i = 0; i < n - 1; i++)
		W[i] = Z[i] - Z[i + 1];

	return W;
}
float* func(float *W, int n, float *Fi, int p)
{
    int i, j;
	float sum;
	float *res = new float[n - p];
	for(i = p; i < n; i++)
	{
		sum = 0;
		for(j = 1; j <= p; j++)
			sum += Fi[j - 1] * W[i - j];

		res[i - p] = W[i] - sum;
	}

	return res;
}

float* Cm(float *C, int q, float *Fi, int p)
{
	int i, j, k;
	float sum;
	float *Fim = new float[p + 1];
	float *Cm = new float[q + 1];
	Fim[0] = -1;
	for(i = 1; i <= p; i++)
		Fim[i] = Fi[i - 1];

	for(i = 0; i <= q; i++)
	{
		sum = 0;
		for(j = 0; j <= p; j++)
			sum += Fim[j] * Fim[j];
		Cm[i] = sum * C[i];

		for(j = 0; j <= p; j++)
		{
			sum = 0;
			for(k = 0; k + j <= p; k++)
				sum += Fim[k] * Fim[k + j];
			sum *= C[i + j] + C[abs(i - j)];

			Cm[i] += sum;
		}
	}
	delete Fim;
	return Cm;
}
float* quadpr(float *Cm, int p, int q, float epsilon)
{
	int i, j, k;
	int step = 0;
	float tmp;
	float *tau = new float[q + 1];
	float *tau_next = new float[q + 1];
	float **T = new float*[q + 1], **T_big = new float*[q + 1];
	for(i = 0; i <= q; i++)
	{
		T[i] = new float[q + 1];
		T_big[i] = new float[2 * (q + 1)];
	}
	float *f = new float[q + 1];

	// Начальные значения компонент вектора тау
	tau[0] = powf(Cm[0], 0.5f);
	for(i = 1; i <= q; i++)
		tau[i] = 0;

	for(;;)
	{
        step++;

		// Получение матрицы T
		for(i = 0; i <= q; i++)
			for(j = 0; j <= q; j++)
			{
				if(i + j <= q) T[i][j] = tau[i + j];
				else T[i][j] = 0;
	
				if(j - i >= 0) T[i][j] += tau[j - i];
				else T[i][j] += 0;
			}
	
		// Создание расширенной матрицы
		for(i = 0; i <= q; i++)
			for(j = 0; j <= 2 * q + 1; j++)
			{
				if(j <= q) T_big[i][j] = T[i][j];
				if(j - q - 1 >= 0) T_big[i][j] = 0;
				if(j - q - 1 == i) T_big[i][j] = 1;
			}
	
		
		// Нахождение обратной матрицы
		for(i = 0; i <= q; i++)
		{
			tmp = T_big[i][i];
			for(j = 0; j <= 2 * q + 1; j++)
				T_big[i][j] /= tmp;	

			for(j = 0; j <= q; j++)
			{
				if(i == j) continue;
				tmp = T_big[j][i];
				for(k = 0; k <= 2 * q + 1; k++)
					T_big[j][k] -= T_big[i][k] * tmp;
			}
		}

		// Вычисление вектора f
		for(i = 0; i <= q; i++)
		{
			tmp = 0;
			for(j = 0; j <= q - i; j++)
				tmp += tau[j] * tau[i + j];

			f[i] = tmp - Cm[i];
		}

		// Вычисление нового вектора тау
		for(i = 0; i <= q; i++)
		{
			tmp = 0;
			for(j = 0; j <= q; j++)
				tmp += T_big[i][q + 1 + j] * f[j];

			tau_next[i] = tau[i] - tmp;
		}

		for(i = 0; i <= q; i++)
			tau[i] = tau_next[i];


		// Проверка на окончание работы
		for(i = 0; i <= q; i++)
			if(f[i] > epsilon) break;
		if(i == q + 1) break;
		
		for(i = 0; i <= q; i++)
			printf(" %f", f[i]);
		printf("\n\n");
		getch();
	}

	printf("%d\n", step);
	for(i = 0; i <= q; i++)
		printf("%f ", tau[i]);
	getch();

    return tau;
}
