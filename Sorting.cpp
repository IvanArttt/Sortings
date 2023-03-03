#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stack>


void PrintArr(int* arr, int size)
{
	for (size_t i = 0; i < size; i++)
	{
		printf("%4d", arr[i]);
	}
	printf("\n");
}

void CopyArr(int* arr1, int* arr2, int size)
{
	for (size_t i = 0; i < size; i++)
	{
		arr2[i] = arr1[i];
	}
}

void Check(int* arr, int* arrRand, int* arrMaxMin, int* arrMinMax, int* arrRandSame, int size, void (*Sorting)(int*, int))
{
	CopyArr(arrRand, arr, size);

	clock_t t = clock();

	Sorting(arr, size);

	t = clock() - t;

	printf("Time for sorting Random Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	//PrintArr(arr, size);
	CopyArr(arrMaxMin, arr, size);

	t = clock();

	Sorting(arr, size);

	t = clock() - t;

	printf("Time for sorting MaxMin Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);

	CopyArr(arrMinMax, arr, size);

	t = clock();

	Sorting(arr, size);

	t = clock() - t;

	printf("Time for sorting MinMax Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);

	CopyArr(arrRandSame, arr, size);

	t = clock();

	Sorting(arr, size);

	t = clock() - t;

	printf("Time for sorting RandomSame Numbers is: %d ticks, %3.3f seconds \n\n\n\n", t, ((float)t) / CLOCKS_PER_SEC);
	//PrintArr(arr, size);
}


void CheckShell(int* arr, int* arrRand, int* arrMaxMin, int* arrMinMax, int* arrRandSame, 
	int size, 
	void (*Sorting)(int* arr, int size, int* h, int t),
	int* h1, int* h2, int* h3, int* h4,
	int t1, int t2, int t3, int t4)
{
	printf("\th1:\n");
	CopyArr(arrRand, arr, size);
	clock_t t = clock();
	Sorting(arr, size, h1, t1);
	t = clock() - t;
	printf("\t\tTime for sorting Random Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMaxMin, arr, size);
	t = clock();
	Sorting(arr, size, h1, t1);
	t = clock() - t;
	printf("\t\tTime for sorting MaxMin Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMinMax, arr, size);
	t = clock();
	Sorting(arr, size, h1, t1);
	t = clock() - t;
	printf("\t\tTime for sorting MinMax Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrRandSame, arr, size);
	t = clock();
	Sorting(arr, size, h1, t1);
	t = clock() - t;
	printf("\t\tTime for sorting RandomSame Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);


	printf("\th2:\n");
	CopyArr(arrRand, arr, size);
	t = clock();
	Sorting(arr, size, h2, t2);
	t = clock() - t;
	printf("\t\tTime for sorting Random Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMaxMin, arr, size);
	t = clock();
	Sorting(arr, size, h2, t2);
	t = clock() - t;
	printf("\t\tTime for sorting MaxMin Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMinMax, arr, size);
	t = clock();
	Sorting(arr, size, h2, t2);
	t = clock() - t;
	printf("\t\tTime for sorting MinMax Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrRandSame, arr, size);
	t = clock();
	Sorting(arr, size, h2, t2);
	t = clock() - t;
	printf("\t\tTime for sorting RandomSame Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);

	printf("\th3:\n");
	CopyArr(arrRand, arr, size);
	t = clock();
	Sorting(arr, size, h3, t3);
	t = clock() - t;
	printf("\t\tTime for sorting Random Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMaxMin, arr, size);
	t = clock();
	Sorting(arr, size, h3, t3);
	t = clock() - t;
	printf("\t\tTime for sorting MaxMin Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMinMax, arr, size);
	t = clock();
	Sorting(arr, size, h3, t3);
	t = clock() - t;
	printf("\t\tTime for sorting MinMax Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrRandSame, arr, size);
	t = clock();
	Sorting(arr, size, h3, t3);
	t = clock() - t;
	printf("\t\tTime for sorting RandomSame Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);

	printf("\th4:\n");
	CopyArr(arrRand, arr, size);
	t = clock();
	Sorting(arr, size, h4, t4);
	t = clock() - t;
	printf("\t\tTime for sorting Random Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMaxMin, arr, size);
	t = clock();
	Sorting(arr, size, h4, t4);
	t = clock() - t;
	printf("\t\tTime for sorting MaxMin Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrMinMax, arr, size);
	t = clock();
	Sorting(arr, size, h4, t4);
	t = clock() - t;
	printf("\t\tTime for sorting MinMax Numbers is: %d ticks, %3.3f seconds \n", t, ((float)t) / CLOCKS_PER_SEC);
	CopyArr(arrRandSame, arr, size);
	t = clock();
	Sorting(arr, size, h4, t4);
	t = clock() - t;
	printf("\t\tTime for sorting RandomSame Numbers is: %d ticks, %3.3f seconds \n\n\n", t, ((float)t) / CLOCKS_PER_SEC);
}


/////////////////////////////////////////////////////////////////
void BubbleSort(int* arr, int size)
{
	for (size_t i = 0; i < size - 1; i++)
	{
		for (size_t j = 0; j < size - i - 1; j++)
		{
			if (arr[j + 1] < arr[j])
			{
				int c = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = c;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
void InsertSort(int* arr, int size)
{
	int temp;
	int i, j;

	for (i = 1; i < size; i++)
	{
		temp = arr[i];
		for (j = i - 1; j > -1 && arr[j] > temp; --j)
		{
			arr[j + 1] = arr[j];
		}
		arr[j + 1] = temp;
	}
}


/////////////////////////////////////////////////////////////////
void ShellSort(int* arr, int size, int* h, int t)
{
	int temp;
	int i, j, r;

	for (r = t - 1; r >= 0; --r)
	{
		for (i = h[r]; i < size; ++i)
		{
			temp = arr[i];
			j = i - h[r];

			while (j > -1 && arr[j] > temp)
			{
				arr[j + h[r]] = arr[j];

				j -= h[r];
			}
			arr[j + h[r]] = temp;
		}
	}
}

void Generate(int* h1, int* h2, int* h3, int* h4, int size, int& t1, int& t2, int& t3, int& t4)
{
	int i = 1;
	while (h1[i - 1] < size)
	{
		h1[i] = 3 * h1[i - 1] + 1;
		++i;
		t1++;
	}

	i = 1;
	while (h2[i - 1] < size)
	{
		h2[i] = 2 * h2[i - 1] + 1;
		++i;
		t2++;
	}

	for (int i = 0; (i < size) && (h3[i] < size); i++)
	{
		if (i % 2 == 0) h3[i] = (int)(9 * pow(2, i) - 9 * pow(2, i / 2) + 1);
		else h3[i] = (int)(8 * pow(2, i) - 6 * pow(2, (i + 1) / 2) + 1);
		t3++;
	}

	int temp1 = 0;
	int temp2 = 0;
	int temp = 0;
	int r = 0;
	for (int i = 0; (i < size) && (temp >= 0); i++)
	{
		for (int j = 0; (j < size) && (temp >= 0); j++)
		{
			temp1 = (int)pow(2, i);
			temp2 = (int)pow(3, j);
			temp = temp1 * temp2;
			if (temp < size)
			{
				h4[r] = temp;
				r++;
			}
		}
	}
	t4 = r - 1;
	InsertSort(h4, r-1);
}

/////////////////////////////////////////////////////////////////


template <typename T>
void TurnamentSort(T* k, int size)
{
	int tmp = size;
	int d = 0;
	while (size > pow(2, d))
	{
		d++;
	}
	size = pow(2, d);
	T* arr = new T[size];
	for (int i = 0; i < size; i++)
	{
		if (i < tmp)
			arr[i] = k[i];
		else
			arr[i] = INT_MIN;
	}
	int* index = new int[size - 1];
	for (int i = pow(2, d - 1) - 1; i >= 0; i--)
	{
		if (arr[2 * i + 1] >= arr[2 * i])
			index[(size + 2 * i - 2) / 2] = 2 * i + 1;
		else
			index[(size + 2 * i - 2) / 2] = 2 * i;
	}
	for (int i = size / 2 - 2; i >= 0; i--)
	{
		if (arr[index[2 * i + 2]] >= arr[index[2 * i + 1]])
			index[i] = index[2 * i + 2];
		else
			index[i] = index[2 * i + 1];
	}
	int s = tmp - 1;
	while (s >= 0)
	{
		k[s] = arr[index[0]];
		arr[index[0]] = INT_MIN;
		if (index[0] % 2 == 0)
			index[(size + index[0] - 2) / 2] = index[0] + 1;
		else
			index[(size + index[0] - 3) / 2] = index[0] - 1;
		for (int i = size / 2 - 2; i >= 0; i--)
		{
			if (arr[index[2 * i + 2]] >= arr[index[2 * i + 1]])
				index[i] = index[2 * i + 2];
			else
				index[i] = index[2 * i + 1];
		}
		s--;
	}
	delete[]arr, index;
}

////////////////////////////////////////////////////////////////////////////////////



template <typename T>
void HeapDown(T* k, int n, int i)
{
	int j;
	while (2*i + 1 < n)
	{
		j = 2 * i + 1;
		if ((j + 1) < n && k[j] < k[j + 1]) ++j;
		if (k[i] < k[j])
		{
			std::swap(k[i], k[j]);
			i = j;
		}
		else return;
	}
	return;
}

template <typename T>
void HeapSort(T* k, int size)
{
	int j, n = size;
	T temp, * pq = k;

	for (j = (n - 1)/2; j >= 0; --j) HeapDown(pq, n, j);
	while (n > 1)
	{
		temp = pq[0];
		pq[0] = pq[n - 1];
		pq[--n] = temp;
		HeapDown(pq, n, 0);
	}
}

/////////////////////////////////////////////////////////////



template <typename T>
void QuickSort(T* k, int left, int right)
{
	if (right <= left) return;
	int i = left + 1, p = right;
	T temp, x = k[(left + right) / 2];
	while (true)
	{
		while (i < p && k[i] <= x) ++i;
		while (p >= i && k[p] > x) --p;
		if (i < p)
		{
			temp = k[i];
			k[i++] = k[p];
			k[p--] = temp;
		}
		else break;
	}
	k[left] = k[p];
	k[p] = x;

	QuickSort(k, left, p - 1);
	QuickSort(k, p+1, right);
}

template <typename T>
void QSort(T* k, int size)
{
	QuickSort(k, 0, size - 1);
}


template <typename T>
int QuickSortLeft(T* k, int left, int right)
{
	int i = left + 1, p = right;
	T temp;
	T x = k[left];
	while (true)
	{
		while (i < p && k[i] <= x)
			++i;
		while (p >= i && k[p] > x)
			--p;
		if (i < p)
		{
			temp = k[i];
			k[i++] = k[p];
			k[p--] = temp;
		}
		else
			break;
	}
	k[left] = k[p];
	k[p] = x;
	return p;
}


template <typename T>
void QSortNotRecursive(T* k, int left, int right) //быстрая сортировка нерекурсивная
{
	std::stack<int> A;
	int s = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	A.push(left);
	A.push(right);
	while (!A.empty())
	{
		int r = A.top();
		A.pop();
		int l = A.top();
		A.pop();
		int index = QuickSortLeft(k, l, r);
		if ((index - 1) > l)
		{
			A.push(l);
			A.push(index - 1);
		}
		if ((index + 1) < r)
		{
			A.push(index + 1);
			A.push(r);
		}
	}
}


template <typename T>
void HolandFlag(T* arr, int size)
{
	int left = 0, right = size - 1;
	while (left < right)
	{
		while ((arr[left] == 0) && (left < right)) left++;
		while ((arr[right] == 0) && (left < right)) right--;
		if (left < right)
		{
			std::swap(arr[left], arr[right]);
			left++;
			right--;
		}
	}
}

/////////////////////////////////////////////////////////////////////////


template <typename T>
void MergeSortNoRecursive(T* k, int pos, int m, int p)
{
	int i, j, r;
	T* work = new T[m + p];
	for (i = pos; i < pos + m; ++i)
	{
		work[i] = k[i];
	}
	for (j = 0; j < p; ++j)
	{
		work[pos + p + m - j - 1] = k[pos + m + j];
	}
	i = 0; j = pos + p + m - 1;
	for (r = pos; r < pos + m + p; ++r)
	{
		if (work[j] < work[i])
			k[r] = work[j--];
		else
			k[r] = work[i++];
	}
	delete[] work;
}


template <typename T>
void Merge(T* k, int left, int mid, int right, int s)
{
	int h, i, j;
	T* arr = new T[s];
	h = left;
	i = left;
	j = mid + 1;
	while ((h <= mid) && (j <= right))
	{
		if (k[h] <= k[j])
		{
			arr[i] = k[h];
			h++;
		}
		else
		{
			arr[i] = k[j];
			j++;
		}
		i++;
	}
	if (h > mid)
	{
		for (int m = j; m <= right; m++)
		{
			arr[i] = k[m];
			i++;
		}
	}
	else
	{
		for (int m = h; m <= mid; m++)
		{
			arr[i] = k[m];
			i++;
		}
	}
	for (int m = left; m <= right; m++)
	{
		k[m] = arr[m];
	}
	delete[] arr;
}

template <typename T>
void MergeSort(T* k, int left, int right, int s)
{
	int mid;
	if (left < right) {
		mid = (left + right) / 2;
		MergeSort(k, left, mid, s);
		MergeSort(k, mid + 1, right, s);
		Merge(k, left, mid, right, s);
	}
}


////////////////////////////////////////////////////////////////////


int digit(int k, int d, int p, int s)
{
	for (int i = 1; i < s - p + 1; ++i)
	{
		k /= d;
	}
	return k % d;
}


template <typename T>
void RadixMSD(T* k, int left, int right, int p, int s, int d) //поразрядная сортировка
{
	int* count = new int[d + 1];
	int* count1 = new int[d + 1];
	if (p > s)
		return;
	int i, j, n;
	n = right - left + 1;
	if (n <= 0)
		return;
	T* work = new T[n];
	for (j = 0; j <= d; ++j)
	{
		count[j] = 0;
		count1[j] = 0;
	}
	for (i = left; i <= right; ++i)
	{
		++count[digit(k[i], d, p, s) + 1];
	}
	for (j = 1; j <= d; ++j)
	{
		count[j] += count[j - 1];
	}
	for (j = 1; j <= d; ++j)
	{
		count1[j] = count[j];
	}
	for (i = left; i <= right; ++i)
	{
		work[count1[digit(k[i], d, p, s)]++] = k[i];
	}

	for (i = left; i <= right; ++i)
	{
		k[i] = work[i - left];
	}
	for (j = 0; j < d; ++j)
	{
		SortRazryad(k, left + count[j], left + count[j + 1] - 1, p + 1, s, d);
	}
	delete[] work;
	delete[] count;
	delete[] count1;
}

//////////////////////////////////////////////////////////////////////////////////////////



template <typename T>
void CountSort(T* k, int left, int right, int diapazon)
{
	int s = right + left + 1;
	int* count = new int[diapazon + 1];
	for (int i = 0; i < diapazon + 1; i++)
	{
		count[i] = 0;
	}
	for (int i = left; i < s; i++)
	{
		count[k[i]]++;
	}
	int f = 0;
	for (int i = 0; i < diapazon + 1; i++)
	{
		for (int j = 0; j < count[i]; j++)
		{
			k[f] = i;
			f++;
		}
	}
	delete[] count;
}


int main()
{
	srand(time(NULL));
	const int size = 30000;

	int* arrRand = new int[size];
	int* arrMaxMin = new int[size];
	int* arrMinMax = new int[size];
	int* arrRandSame = new int[size];
	int* arrBubble = new int[size];
	int* arrInsert = new int[size];
	int* arrShell = new int[size];
	int* arrTournir = new int[size];
	int* arrHeap = new int[size];
	int* arrQ = new int[size];
	int* arrMerge = new int[size];

	for (size_t i = 0; i < size; i++)
	{
		arrRand[i] = rand() % (10 * size);
	}

	for (size_t i = 0; i < size; i++)
	{
		arrRandSame[i] = rand() % 10;
	}

	for (size_t i = 0; i < size; i++)
	{
		arrMinMax[i] = i;
	}

	for (size_t i = 0; i < size; i++)
	{
		arrMaxMin[i] = size - i;
	}

	printf("BubbleSort:\n");
	Check(arrBubble, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, BubbleSort);

	printf("InsertSort:\n");
	Check(arrInsert, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, InsertSort);

	int* h1, * h2, * h3, * h4;
	int t1 = 0;
	int t2 = 0;
	int t3 = 0;
	int t4 = 0;
	h1 = new int[(int)(log(size) / log(3)) + 1];
	h2 = new int[(int)(log(size) / log(2)) + 1];
	h3 = new int[(int)(log(size) / log(2)) + 1];
	h4 = new int[size];
	h1[0] = h2[0] = h3[0] = h4[0] = 1;
	Generate(h1, h2, h3, h4, size, t1, t2, t3, t4);
	printf("ShellSort:\n");
	CheckShell(arrShell, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, ShellSort, h1, h2, h3, h4, t1, t2, t3, t4);

	printf("TournirSort:\n");
	Check(arrTournir, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, TurnamentSort);

	printf("HeapSort:\n");
	Check(arrHeap, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, HeapSort);


	printf("QSort:\n");
	Check(arrQ, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, QSort);
	
	//printf("MergeSort:\n");
	//Check(arrMerge, arrRand, arrMaxMin, arrMinMax, arrRandSame, size, MergeSort);
	
	/*for (size_t i = 0; i < size; i++)
	{
		printf("%d\t", arrTournir[i]);
	}*/


	return 0;
}