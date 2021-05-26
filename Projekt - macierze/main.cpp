#include "lib.h"
#include <iostream>



/* Napisa� bibliotek� numeryczn� do rozwi�zywania uk�ad�w r�wna� liniowych.
Biblioteka ma operowa� na obiektach klasy macierz kt�re tworzone s� dynamicznie.
Algorytmy(solvery) licz�ce uklad r�wna� powinny by� podzielone na 2 grupy:
iteracyne(Jackobiego i Gaussa - Seidl'a) oraz og�lne � metoda Gaussa */

int main() {
	double matrixA_left[12] = {6,7,1,5,3,4,2,2,1,4,4,2}; 
	double matrixA_right[4] = {7,5,3,6}; 
	/* tablice reprezentujace uklad:
	6a + 7b + 1c = 7
	5a + 3b + 4c = 5
	2a + 2b + 1c = 3
	4a + 4b + 2c = 6

	rozwiazanie:
	a = -3
	b = 3.2
	c = 2.6

	zrodlo: 
	https://matrixcalc.org/en/slu.html#solve-using-Gaussian-elimination%28%7B%7B6,7,1,0,7%7D,%7B5,3,4,0,5%7D,%7B2,2,1,0,3%7D,%7B4,4,2,0,6%7D%7D%29

	*/

	double matrixB_left[16] = {1, 8, 2, 2, 3, -2, -14, 1, 0, -1, 4, 6, 10, 2, 0 ,-1};
	double matrixB_right[4] = { 5,3,1, 7};
	/*
	
	rozwiazania:
	a = 2279/3592 (ok. 0.63446547884)
	b = 887/1796 (ok. 0.49387527839)
	c = (-899)/7184 (ok. -0.12513919821)
	d = 597/1796 (ok. 0.33240534521)


	
	
	
	*/
	Macierz * A = new Macierz(3, 4, matrixA_left, matrixA_right);
	Macierz * B = new Macierz( 4, 4, matrixB_left, matrixB_right);
	cout << "Macierz A: " << endl;
	A->print();
	
	bool isSolvable = A->gaussian_elimination_rearrange();
	if (isSolvable) {
		cout << "Wynik dla metody Gaussa: " << endl;
		double * solutions = A->gaussian_elimination_solve();

		for (int i = 0;i < A->szer();i++) {
			cout << solutions[i] << endl;
		}

		delete[] solutions;
	}
	else {
		cout << "Brak rozwiazan" << endl;
	}
	
	cout << "Macierz B: " << endl;
	B->print();
	
	bool isGood = B->seidel_jacobi_test_compatibility();
	cout << "Zmodyfikowana macierz B: " << endl;
	B->print();
	if (isGood) {
		cout << "Wynik dla metody Gaussa-Seidla (30sta iteracja): " << endl;
		double * sol1 = B->seidel_solve(30);
		for (int i = 0;i < B->szer();i++) {
			cout << sol1[i] << endl;
		}
		cout << endl;
		cout << "Wynik dla metody Jacobiego (dozwolony blad: 0,0001):" << endl;
		double * sol2 = B->jacobi_solve(0.0001);
		for (int i = 0;i < B->szer();i++) {
			cout << sol2[i] << endl;
		}

		delete[] sol1;
		delete[] sol2;
	}
	else {
		cout << "Uklad musi byc sprowadzalny do macierzy przekatniowo dominujacej!" << endl;
	}

	
		
		delete A;
		delete B;
		
		
	system("PAUSE");
}