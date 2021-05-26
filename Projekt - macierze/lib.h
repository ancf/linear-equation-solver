#pragma once

#include <iostream>
using namespace std;

class Macierz
{
private:
	double * tab;
	double * results_tab;
	size_t t_szer;
	size_t t_wys;
public:
	Macierz(int t_szer, int t_wys, double * tab_, double * results_tab); 
	/* tworzy reprezentacje ukladu rownan:
	tablica tab reprezentuje macierz za� results_tab reprezentuje praw� stron� uk�adu r�wna�
	t_szer i t_wys m�wi� o szeroko�ci (liczba zmiennych) i wysoko�ci (liczba r�wna� i ich rozwi�za�)
	jako �e u�ywamy tablicy jednowymiarowej [t_wys * t_szer] 
	*/
	~Macierz();

	double operator()(int x, int y) const
	{
		return tab[x * t_szer + y];
	}
	double& operator()(int x, int y)
	{
		return tab[x * t_szer + y];
	}

	int szer() const
	{
		return t_szer;
	}

	int wys() const
	{
		return t_wys;
	}

	Macierz(Macierz & m_) :
		tab(m_.tab),
		results_tab(m_.results_tab),
		t_szer(m_.t_szer),
		t_wys(m_.t_wys)
	{
	}

	Macierz& operator=(Macierz const& m_); 

	double * jacobi_solve(double allowed_error); 
	/*
		znajduje rozwiazanie macierzy kwadratowej przek�tniowo dominuj�cej o podanej dok�adno�ci metod� Jacobiego (zwraca tablic� rozwi�za�)
	*/

	bool seidel_jacobi_test_compatibility(); 
	/*
		sprawdza czy macierz jest kwadratowa i przek�tniowo dominujaca, 
		je�li to mo�liwe zmienia kolejno�� wierszy by przekszta�ci� j� na tak�	*/

	double * seidel_solve(int iteration_limit);
	/*
	znajduje rozwiazanie macierzy kwadratowej przek�tniowo dominuj�cej dla podanej iteracji metod� Gaussa-Seidla (zwraca tablic� rozwi�za�)
	*/

	bool gaussian_elimination_rearrange();
	/*
		je�li to mo�liwe przekszta�ca macierz do postaci schodkowej, je�li jest to niemo�liwe lub uk�ad jest sprzeczny zwraca false; 
	*/

	double * gaussian_elimination_solve();
	/*
		zwraca tablic� rozwi�za� macierzy schodkowej
	*/

	int Macierz::gaussian_elimination_arg_max(int wier, int t_wys, int t_szer, double * tab, int kol);
	


	void print();
	/*
		wypisuje reprezentacje ukladu rownan
	*/
	
};
