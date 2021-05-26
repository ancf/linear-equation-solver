#include "lib.h"
#include <stdlib.h>



Macierz::Macierz(int szer_, int wys_, double  * tab_, double * results_tab_) {
	t_szer = szer_;  // szerokosc macierzy (liczba zmiennych) ukladu rownan
	t_wys = wys_; // wysokosc macierzy (liczba rownan)
	tab = new double[t_szer * t_wys];
	results_tab = new double[t_wys];

	for (int i = 0;i < (t_szer * t_wys);i++) {
		tab[i] = tab_[i];
	}

	for (int i = 0;i < t_wys;i++) {
		results_tab[i] = results_tab_[i];
	}
}

Macierz::~Macierz()
{
	delete[] tab;
	delete[] results_tab;
}



double * Macierz::seidel_solve(int iteration_limit) {
	int n = t_szer;
	double * solutions = new double[n];
	double * guesses = new double[n];
	

	for (int i = 0;i < n;i++) {
		guesses[i] = 0.0;
	//	cout << guesses[i] << endl;
	}
	
	while (iteration_limit > 0)
	{
		
		for (int j = 0; j < n; j++)
		{
			
			solutions[j] = results_tab[j] / tab[j*t_szer + j]; // wartosc prawej strony rownania / wartosc na przekatnej (dominujaca)
			
			for (int k = 0; k < n; k++)
			{
				if (k != j) {
					solutions[j] = solutions[j] - (tab[j*t_szer + k] / tab[j*t_szer + j]) * guesses[k];
					guesses[j] = solutions[j];
				}
			}

			

		//	cout<< solutions[j] << endl;
		}
	//	cout << endl;


	

		iteration_limit--;
	}
	
	delete[] guesses;
	return solutions; 
	
}

bool Macierz::seidel_jacobi_test_compatibility() {
	bool compatibile = true;
	if (t_szer != t_wys) {
	//	cout << "a" << endl;
		compatibile = false;
	}
	else {
		int n = t_wys;
		int * rows = new int[n]; // indeks - nr wiersza, wartoœæ - nr kolumny
		bool * checked_columns = new bool[n];
		for (int i = 0;i < n;i++) {
			checked_columns[i] = false;
		}

		for (int j = 0;j < n;j++) {
			rows[j] = 0;
			for (int k = 0;k < n;k++) {
				//	cout << tab[j*n + rows[j]] << " " << tab[j*n + k] << endl;
				//	cout << j*n + rows[j] << " " << j*n + k << endl;
				//	cout << endl;
				if (fabs(tab[j*n + rows[j]]) < fabs(tab[j*n + k])) { 
					/* je¿eli dla danego wiersza wartosc bezwzglêdna w kolumnie k 
					jest wiêksza ni¿ w kolumnie zapisanej w tablicy to zapisz (nadpisz) wartoœæ w tablicy */

					//	cout << tab[j*n + rows[j]] << " " << tab[j*n + k] << endl;
					//	cout << j*n + rows[j] << " " << j*n + k << endl;
					//	cout << endl;


				//	cout << "j: " << j << " k: " << k << endl;
					rows[j] = k;

				}
			}
			if (checked_columns[rows[j]] == false) {
				checked_columns[rows[j]] = true;
			}
			else {
				//	cout << "b" << endl;
				compatibile = false;
				return compatibile;
			} /* wartoœci dominuj¹ce nie mog¹ wyst¹piæ wiêcej ni¿ 1 raz w danej kolumnie
			  */
		}

		double * tmp = new double[n*n];
		for (int b = 0;b < n*n;b++) {
			tmp[b] = tab[b];
		}
		double * results_tmp = new double[n];
		for (int b = 0;b < n;b++) {
			results_tmp[b] = results_tab[b];
		}



		for (int p = 0; p < n; p++) {



			results_tab[rows[p]] = results_tmp[p];

			for (int q = 0; q < n; q++) {
				tab[rows[p] * n + q] = tmp[p *n + q];

			}
			//	print();
		}



	//	cout << endl;

		/*sprawdzamy czy wartoœæ w danej kolumnie jest nie tylko najwiêksza, ale i dominuj¹ca w wierszu*/

		for (int d = 0; d < n; d++) {
			double row_sum = 0;
			for (int e = 0; e < n; e++) {
				if (e != d) {
				//	cout << fabs(tab[d*n + e]) << endl;
					row_sum += fabs(tab[d*n + e]);
				}

			}
		//	cout << endl;
		//	cout << row_sum <<endl;
			if (row_sum >= fabs(tab[d*n + d])) { 
				compatibile = false;
				return compatibile;
			}

		}
		delete[] rows;
		delete[] checked_columns;
		delete[] tmp;
		delete[] results_tmp;

	}


	
//	cout << "c" << endl;
	return compatibile;
}
double * Macierz::jacobi_solve(double allowed_error) {
	int n = t_wys;
	double * solutions = new double[n];
	double * guesses = new double[n];
	
	
	for (int i = 0; i < n;i++) {
		guesses[i] = 0;
	}
	
	int k = 0;
	while (k != n)
	{
		k = 0;
		for (int j = 0;j<n;j++)
		{
			solutions[j] = (1 / tab[j*n+j])*results_tab[j]; 
			for (int p = 0;p<n;p++)
			{
				if (p != j) {
					solutions[j] = solutions[j] - (1 / tab[j*n + j])*(tab[j*n + p] * guesses[p]);
				}
					
			}
		}
		for (int q = 0;q<n;q++)
		{
			double y = fabs(solutions[q] - guesses[q]);

			
			if (y <= allowed_error)
			{
				k++;
			}
		}
		for (int r = 0;r < n;r++) {
			guesses[r] = solutions[r];
		}
			
	}
	
	delete[] guesses;

	
	return solutions;
}

int Macierz::gaussian_elimination_arg_max(int wier, int t_wys, int t_szer, double * tab, int kol) { // znajdz numer wiersza z najwieksza wartoscia absolutna w danej kolumnie
	int tmp = wier;
	for (int i = wier; i < t_wys; i++) {
		//	cout << "porownaj:"<< abs(tab[i*t_szer + kol]) << " "<< abs(tab[tmp*t_szer + kol]) << endl;
		if (fabs(tab[i*t_szer + kol]) > fabs(tab[tmp*t_szer + kol])) {
			tmp = i;
		}
	}
	// cout <<"maks: " << tmp <<endl;
	return tmp;
}

bool Macierz::gaussian_elimination_rearrange() {
	bool solutions_exist = true;
	if (t_wys < t_szer) { // mniej ukladow rownan niz zmiennych
		solutions_exist = false;
	}
	else {
		int wier = 0;
		int kol = 0;

		//	cout << solutions[0] << endl;

		while (wier < t_wys || kol < t_szer) {
			//	print();
			int maks = gaussian_elimination_arg_max(wier, t_wys, t_szer, tab, kol);
			if (tab[maks*t_szer + kol] == 0) {
				kol++;
				//	print();
			}
			else {

				double tmp;
				for (int k = 0; k < t_szer; k++) {
					if (!(tab[(maks + 1)*t_szer - 1] == 0)) {
						//		cout << "zamien " << tab[maks * t_szer + k] << " " << tab[wier * t_szer + k] << endl;
						//		cout << "adresy" << maks * t_szer + k << " " << wier * t_szer + k << endl;
						tmp = tab[maks * t_szer + k];
						tab[maks * t_szer + k] = tab[wier * t_szer + k];
						tab[wier * t_szer + k] = tmp;
						swap(results_tab[maks], results_tab[wier]);
					}


				}

				//	print();

				for (int i = wier + 1; i < t_wys; i++) {
					if (!(tab[(maks + 1)*t_szer - 1] == 0)) {
						//	cout << "f= " << tab[i*t_szer + kol] << " / " << tab[wier*t_szer + kol] << endl;
						double f = tab[i*t_szer + kol] / tab[wier*t_szer + kol];
						//	cout << f << endl;
						//	cout << "wyzeruj " << i*t_szer + kol << endl;
						tab[i*t_szer + kol] = 0;

						for (int j = kol + 1; j < t_szer; j++) {
							//		cout << tab[i*t_szer + j] << "-" << tab[wier*t_szer + j] << "* f" << endl;
							tab[i*t_szer + j] = tab[i*t_szer + j] - (tab[wier * t_szer + j] * f);

						}
						//	cout << results_tab[i] << "-" << results_tab[wier] << "* f" << endl;
						results_tab[i] = results_tab[i] - results_tab[wier] * f;
						//	print();
					}
				}
				if (wier < t_wys) {
					wier++;
				}
				if (kol < t_szer) {
					kol++;
				}

			}
			//	cout << solutions[0] << endl;
		}
		//	cout << solutions[0] << endl;



		if (t_wys > t_szer) {
			for (int i = t_szer; i < t_wys;i++) {
				if (results_tab[i] != 0) {
					solutions_exist = false;
				}
			}
		}

		
	}
	return solutions_exist;
}

double * Macierz::gaussian_elimination_solve() {
	int n = t_szer;

	double * solutions = new double[n];
	for (int p = 0; p < n;p++) {
		solutions[p] = 0.0;
	}


	

	
		for (int i = n- 1; i >= 0; i--) {
			double sum = 0;

			if (i < n - 1)
				for (int j = i + 1; j < n; j++) {
					//	cout << tab[i*t_szer + j] << "*" <<solutions[j]<< endl;
					sum += tab[i*n + j] * solutions[j]; 
					/* sumowanie iloczynow rozwiazan i wspolczynnikow na prawo od kolumny dla ktorej szukamy rozwiazania */
				}

			solutions[i] = (results_tab[i] - sum) / tab[i*n + i];
		//	cout << solutions[i] << endl;
		}
	
	


	return solutions;

}

void Macierz::print() {
	for (int i = 0; i < t_wys; i++) {
		for (int j = 0; j < t_szer;j++) {
			cout << tab[i*t_szer + j] << " ";
		}
		cout << " = " << results_tab[i] << endl;

	}
}