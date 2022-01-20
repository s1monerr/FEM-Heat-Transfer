#ifndef sort
#define sort

void swap(int& a, int& b) {
	int buf = a;
	a = b;
	b = buf;
}

void heapify(int n, int* tab, int heapsize) {
	int left = 2 * n;
	int right = 2 * n + 1;
	int large;
	if (left <= heapsize && tab[left] > tab[n])
		large = left;
	else
		large = n;

	if (right <= heapsize && tab[right] > tab[large])
		large = right;
	if (large != n)
	{
		swap(tab[n], tab[large]);
		heapify(large, tab, heapsize);
	}
}


void heap_all(int* tab, int heapsize) {
	int n = heapsize / 2;
	while (n > 0)
	{
		heapify(n, tab, heapsize);
		n--;
	}
}


void heapsort(int* tab, int n){ // wersja dla tablicy n - elementow(wszystkie indexy wykorzystane)
	--tab; // wtedy indeks tab[1] wskazuje na tab[0]
	heap_all(tab, n);
	while (n > 1)
	{
		swap(tab[1], tab[n]);
		heapify(1, tab, --n);
	}
}
#endif