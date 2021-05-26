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
	tablica tab reprezentuje macierz zaœ results_tab reprezentuje praw¹ stronê uk³adu równañ
	t_szer i t_wys mówi¹ o szerokoœci (liczba zmiennych) i wysokoœci (liczba równañ i ich rozwi¹zañ)
	jako ¿e u¿ywamy tablicy jednowymiarowej [t_wys * t_szer] 
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
		znajduje rozwiazanie macierzy kwadratowej przek¹tniowo dominuj¹cej o podanej dok³adnoœci metod¹ Jacobiego (zwraca tablicê rozwi¹zañ)
	*/

	bool seidel_jacobi_test_compatibility(); 
	/*
		sprawdza czy macierz jest kwadratowa i przek¹tniowo dominujaca, 
		jeœli to mo¿liwe zmienia kolejnoœæ wierszy by przekszta³ciæ j¹ na tak¹	*/

	double * seidel_solve(int iteration_limit);
	/*
	znajduje rozwiazanie macierzy kwadratowej przek¹tniowo dominuj¹cej dla podanej iteracji metod¹ Gaussa-Seidla (zwraca tablicê rozwi¹zañ)
	*/

	bool gaussian_elimination_rearrange();
	/*
		jeœli to mo¿liwe przekszta³ca macierz do postaci schodkowej, jeœli jest to niemo¿liwe lub uk³ad jest sprzeczny zwraca false; 
	*/

	double * gaussian_elimination_solve();
	/*
		zwraca tablicê rozwi¹zañ macierzy schodkowej
	*/

	int Macierz::gaussian_elimination_arg_max(int wier, int t_wys, int t_szer, double * tab, int kol);
	


	void print();
	/*
		wypisuje reprezentacje ukladu rownan
	*/
	
};
