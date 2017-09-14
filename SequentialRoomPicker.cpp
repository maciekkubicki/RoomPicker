/**
 * Przydział pokoi metodą simulated annealing
 * wersja sekwencyjna
 * make seq - kompilacja
 * make runseq - uruchomienie
 */
#include <iostream>
#include <stdlib.h>    
#include <time.h> 
#include <vector>
#include <ctime>
#include <math.h>
#include <climits>
#include <string>
#include <fstream>

#define TRYS 4//ile razy wykonany jest algorytm
#define MAXIT 1000//warunek stopu - liczba iteracji bez zmian + bez spełnienia warunku prawdopodobieństwa


/**
 *  \brief Fukcja oblicza niezadowolenie
 *  
 *  \param [in] antipation Tabela 'nielubienia'
 *  \param [in] rooms Aktualne przyporzadkowanie do pokoi
 *  \param [in] number Liczba studentow
 *  \return Współczynnik niezadowolenia
 *  
 *  
 */
double antipationCount(double **antipation, int *rooms, int number)
{
	double sum = 0.;
	for (int i = 0; i < number; ++i)
		for (int j = 0; j < number; ++j)
			if (rooms[i] == rooms[j])
				sum += antipation[i][j];
	return sum / 2.;
}


int main()
{
	int numOfStudents, numOfRooms; // liczba studentów, liczba pokoi
	
	
	srand(time(NULL));//ustawienie generatora
	bool boolNew = false;//jesli false to musi byc plik antipation.dat i liczba wprowadzonych studentów musi byc zgodna z plikiem
	
	std::vector<time_t> v;//wektor ziaren
	for (int i = 0; i < TRYS; ++i)
	{
		v.push_back(time(NULL)*rand()+i);
		std::cout << v[i] << std::endl;
	}
	
	double duration;//czas wykonania aplikacji
	std::clock_t start;
	

	//wprowadzenie liczby studentów
	std::cout << "Podaj liczbe studentow: ";
	std::cin >> numOfStudents;
	std::cout << std::endl << "Wprowadono: " << numOfStudents << std::endl;
	numOfRooms = numOfStudents / 2;
	
	start = std::clock();//pomiar czasu

	
	//alokacja pamięci
	int *studRooms = new int[numOfStudents];
	double **antipation = new double *[numOfStudents];
	for (int i = 0; i < numOfStudents; ++i)
		antipation[i] = new double[numOfStudents];
	
	if (boolNew) //tworzenie nowej tabeli antypatii
	{
		for (int i = 0; i < numOfStudents; ++i)
			for (int j = i + 1; j < numOfStudents; ++j)
			{
				antipation[i][j] = (double)(rand() / ((double)(RAND_MAX / 10.)));
				antipation[j][i] = antipation[i][j];
				antipation[i][i] = 0.;
			}
		antipation[numOfStudents - 1][numOfStudents - 1] = 0.;
	}
	else//użycie tabeli z pliku
	{
		std::ifstream filea("antipation.dat");
		if (!filea)
		{
			std::cout << "There was an error opening the file.\n";
			return 0;
		}
		for (int i = 0; i < numOfStudents; ++i) {

			for (int j = 0; j < numOfStudents; ++j)
			{
				if (!(filea >> antipation[i][j]))
					std::cout << "Smth went wrong" << std::endl;
			}
		}


	}
	
	std::vector<double> ants;//wektor z wyliczonymi wartosciami antypatii
	for (int trys = 0; trys < TRYS; ++trys)
	{
		srand(v[trys]);//ustawienie generatora
		
		//generowanie początkowego przyporzdkowania do pokoi
		for (int i = 0; i < numOfStudents; ++i)
			studRooms[i] = INT_MAX;

		for (int i = 0; i < numOfRooms; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int student = (int)(rand()) / ((int)(RAND_MAX / (numOfStudents - 1)));
				//std::cout << student << std::endl;
				studRooms[student] == INT_MAX ? studRooms[student] = i : --j;
			}
		}
		
		std::cout << "Początkowa antypacja " << antipationCount(antipation, studRooms, numOfStudents) << " w iteracji " << trys <<std::endl;

		int it = 0;//liczba iteracji bez zmian
		double T = 1.;//temperatura
		while (it < MAXIT)
		{
			double sNow = antipationCount(antipation, studRooms, numOfStudents);//wyliczenie antypacji przed zamianą pokoi
			int r1 = 0, r2 = 0;
			while (studRooms[r1] == studRooms[r2] || r1 == r2)//wybór studentów którzy się zamienią pokojamia
			{
				r1 = (int)(rand()) / ((int)(RAND_MAX / (numOfStudents - 1)));
				r2 = (int)(rand()) / ((int)(RAND_MAX / (numOfStudents - 1)));
			}

			//zamiana
			int t = studRooms[r1];
			studRooms[r1] = studRooms[r2];
			studRooms[r2] = t;
			
			
			double sNew = antipationCount(antipation, studRooms, numOfStudents);//antypacja po zamianie
			double p = (double)(rand()) / (double)(RAND_MAX);//prawdopodobieństwo
			double dt = sNow - sNew;//różnica antypacji przed i po zmianie
			
			if (sNew < sNow || p < exp(dt / T))//jeśli współczynnik niezadowolenia zmniejszył się lub spełniony jest warunek prawdopodobieństwa to zmiana jest akceptowana
				it = 0;//bo zmiana nastąpiła	
			else//jeśli nie to powrót do poprzedniego układu
			{
				t = studRooms[r1];
				studRooms[r1] = studRooms[r2];
				studRooms[r2] = t;
				++it;

			}

			T *= 0.999;//temperatura

		}

		std::cout << "Final " << antipationCount(antipation, studRooms, numOfStudents) << std::endl;
		ants.push_back(antipationCount(antipation, studRooms, numOfStudents));//dodanie wyliczonej antypatii
	}
	//wybór najbardziej optymalnego rozwiązania
	size_t iMin = 0;
	for (size_t i = 1; i < ants.size(); ++i)
		if (ants[iMin] > ants[i])
			iMin = i;

	std::cout << "Najlepsze rozwiązanie dla proby: " << iMin << " ant: " << ants[iMin] << std::endl;

	//sprzątanie
	for (int i = 0; i < numOfStudents; ++i)
		delete[] antipation[i];
	delete[] antipation;
	delete[] studRooms;
	
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;//koniec pomiaru czasu
	std::cout << "Czas: " << duration << std::endl;
	
	return 0;
}
