/**
 * Przydział pokoi metodą simulated annealing
 * wersja z MPI
 * make - kompilacja
 * make run - uruchomienie
 */

#if __GNUC__ > 3
 #include <string.h>
 #include <iostream>
#else
 #include <iostream.h>
#endif
#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 
#include <vector>
#include <ctime>
#include <math.h>
#include <climits>
#include <string>
#include <fstream>

#include "sprng_cpp.h" 
#include "mpi.h" 



#define MAXIT 1000//warunek stopu - liczba iteracji bez zmian + bez spełnienia warunku prawdopodobieństwa
#define MASTER 0//numer wezła głównego


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

/**
 *  Struktura przechowujaca wartość współczynnika niezadowolenia w danym węźle rank.
 */
struct antStruct
{
		double ant;
		int rank;
};


int main(int argc, char* argv[])
{
	
	int myid, numprocs; // id wezła, liczba procesów
	int numOfStudents, numOfRooms; // liczba studentów, liczba pokoi
	double duration; // czas wykonania programu
	
	std::clock_t start;
	std::ofstream fs;

	Sprng *stream;	//generator
	int seed, gtype;	//ziarno i typ generatora

	antStruct minimum; // struktura do której przypiasne zostanie najlepsze rozwiazanie

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	
	if (numprocs < 2)//jeśli liczba wezłów jest mniejsza niż 2 to koniec programu
	{
		MPI_Finalize();
		return 0;
	}

	double **antipation;//tabela antypacji
	double *temp;
	
	int	*studRooms;//tablica przechowujaca przypiasnie studentow do pokoi
	int *finalRooms;//finalne rozwiazanie
	
	
	if (myid == MASTER)//wczytanie liczby studentow i typu generatora
	{

		std::cout << "Num of procs " << numprocs << std::endl;
		std::cout << "Podaj liczbe studentow: ";
		std::cin >> numOfStudents;
		if (numOfStudents % 2 == 1)
			++numOfStudents;
		std::cout << std::endl << "Przyjeta liczba studentow: " << numOfStudents << std::endl;
		
		#include "gen_types_menu.h"
        std::cout << "Type in a generator type (integers: 0,1,2,3,4,5):  ";
		std::cin >> gtype;
   
		start = std::clock();//start pomiaru czasu


	}
	
	MPI_Bcast(&numOfStudents, 1, MPI_INT, 0, MPI_COMM_WORLD);//rozesłanie liczby studentów
	MPI_Bcast(&gtype,1,MPI_INT,0,MPI_COMM_WORLD );//rozesłanie typu generatora
	
	seed = make_sprng_seed();//generowanie ziarna
	
	//inicjalizacja generatora
	stream = SelectType(gtype);
	stream->init_sprng(myid,numprocs,seed,SPRNG_DEFAULT);
	
	//alokacja pamieci
	antipation = new double *[numOfStudents];
	temp = new double[numOfStudents*numOfStudents];
	for (int i = 0; i < numOfStudents; ++i)
		antipation[i] = &temp[i*numOfStudents];

	if (myid == MASTER)//generowanie tabeli antypatii i zapis do pliku
	{

		for (int i = 0; i < numOfStudents; ++i) {
			for (int j = i + 1; j < numOfStudents; ++j) {
				antipation[i][j] = (double)(stream->sprng() * 10.);
				antipation[j][i] = antipation[i][j];
				antipation[i][i] = 0.;
			}
		}
		antipation[numOfStudents - 1][numOfStudents - 1] = 0.;
		
		fs.open("antipation.dat", std::fstream::out | std::ios::trunc);
		for (int i = 0; i < numOfStudents; ++i)
		{
			for (int j = 0; j < numOfStudents; ++j)
			{
				fs << antipation[i][j] << " ";
			}
			fs << std::endl;
		}
		fs.close();
	}

	MPI_Bcast(&antipation[0][0], numOfStudents * numOfStudents, MPI_DOUBLE, 0, MPI_COMM_WORLD);//rozesłanie tabeli antypatii
	

	//std::cout << "Test wysylania " << myid << " " << numOfStudents << " MY seed " << seed << " " << numOfStudents << std::endl;
	
	numOfRooms = numOfStudents / 2; //pokoje są dwuosobowe
	studRooms = new int[numOfStudents]; //alokacja pamięci na przydział do pokoi
	
		//generowanie początkowego przydziału do pokoi
		for (int i = 0; i < numOfStudents; ++i)
			studRooms[i] = INT_MAX;

		for (int i = 0; i < numOfRooms; ++i)
			for (int j = 0; j < 2; ++j)
			{
				int student = (int)(stream->sprng() * (numOfStudents));
				//std::cout << student << std::endl;
				studRooms[student] == INT_MAX ? studRooms[student] = i : --j;
			}
	
	
	
	std::cout << "Proces: " << myid << " antypacja poczatkowa: " << antipationCount(antipation, studRooms, numOfStudents) << std::endl;
	
	int it = 0;//liczba iteracji bez zmian
	double T = 1.;//temperatura
	
	while (it < MAXIT)
	{
		
		double sNow = antipationCount(antipation, studRooms, numOfStudents);//wyliczenie antypatii
		int r1 = 0, r2 = 0;
		while (studRooms[r1] == studRooms[r2] || r1 == r2)//wybór studentów, którzy zamienią się pokojami
		{
			r1 = (int)(stream->sprng() * (numOfStudents));
			r2 = (int)(stream->sprng() * (numOfStudents));
		}
		
		//zamiana
		int t = studRooms[r1];
		studRooms[r1] = studRooms[r2];
		studRooms[r2] = t;
		
		
		double sNew = antipationCount(antipation, studRooms, numOfStudents);//wyliczenie współczynnika antypatii po zmianie
		double p = (double)(stream->sprng());//prawdopodobeństow
		double dt = sNow - sNew;//różnica wpsółczynników antypatii przed i po zmianie
		if (sNew < sNow || p < exp(dt / T))//jeśli współczynnik niezadowolenia zmniejszył się lub spełniony jest warunek prawdopodobieństwa to zmiana jest akceptowana
			it = 0;//bo zmiana nastąpiła	
		else//jeśli nie to powrót do poprzedniego układu
		{
			t = studRooms[r1];
			studRooms[r1] = studRooms[r2];
			studRooms[r2] = t;
			++it;

		}

		T *= 0.999;//spadek temperatury
	 

	}

	std::cout << "Proces: " << myid << " antypacja koncowa: " << antipationCount(antipation, studRooms, numOfStudents) << std::endl;
	
	antStruct local;//struktura ze współczynnikiem antypatii na danym węźle
	local.rank = myid;
	local.ant = antipationCount(antipation, studRooms, numOfStudents);
	
	MPI_Allreduce(&local, &minimum, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);//redukcja każdy do każdego minimum
	//MPI_Reduce(&local, &minimum, 1, MPI_DOUBLE_INT, MPI_MINLOC, MASTER, MPI_COMM_WORLD);
	
	if (myid == MASTER)
		std::cout << "Minimum w procesie nr: " << minimum.rank << " wartosc: " << minimum.ant << std::endl;
	
	if (minimum.rank == MASTER)//jeśli minimumjest w węźle MASTER to je zapisz do pliku
	{
		fs.open("result.dat",  std::fstream::out | std::ios::trunc);
		fs << "Student " << "Room" << std::endl;
		finalRooms = new int[numOfStudents];
		for (int i = 0; i < numOfStudents; ++i)
		{
			finalRooms[i] = studRooms[i];
			//std::cout << "Student id: " << i << " Pokoj: " << finalRooms[i] << std::endl;
			fs << i << "         " << finalRooms[i] << std::endl;
		}
		fs.close();
		delete[] finalRooms;
	}
	else if (myid == minimum.rank)//jeśli minimum jest na innym weźle niż MASTER to je tam wyślij
		MPI_Send(studRooms, numOfStudents, MPI_INT, MASTER, minimum.rank, MPI_COMM_WORLD);//w sumie to wysylanie jest nie potrzebne, bo moznaby bezposrednie z tego węzła zapisać do pliku
	
	if (myid == MASTER && minimum.rank != MASTER)//odebranie na MASTER rozwiązania i zapis do pliku
	{
		fs.open("result.dat", std::fstream::out | std::ios::trunc);
		fs << "Student " << "Room" << std::endl;
		
		finalRooms = new int[numOfStudents];
		MPI_Recv(finalRooms, numOfStudents, MPI_INT, minimum.rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		for (int i = 0; i < numOfStudents; ++i)
		{
			//std::cout << "Student id: " << i << " Pokoj: " << finalRooms[i] << std::endl;
			fs << i << "        " << finalRooms[i] << std::endl;
		}
		fs.close();
		delete[] finalRooms;
	}
	
	

	//sprzątanie
	delete[] antipation[0];
	delete[] antipation;
	delete[] studRooms;
	
	MPI_Barrier(MPI_COMM_WORLD);//zaczekanie na pozostałe procesy
	if (myid == MASTER)
	{
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;//koniec pomiaru czasu
		std::cout << "Czas: " << duration << std::endl;
	}

	MPI_Finalize();

	return 0;
}
