#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

struct Node // struktura wezlow
{
	double x;
	double y;
	bool BC = false;
};

struct Element // struktura elementow
{
	int ID[4];
};

struct Grid // siatka MES
{
	int nN;
	int nE;
	vector <Node> ND;
	vector <Element> EL;

	void readGrid()
	{
		ifstream file("Test3_31_31_kwadrat.txt");

		string tekst;
		int wartosc, war1, war2, war3, war4;
		double wartosc1, wartosc2;

		for (int i = 0; i < 8; i++)
		{
			file >> tekst >> wartosc;
		}

		file >> tekst >> tekst >> wartosc;
		nN = wartosc;

		file >> tekst >> tekst >> wartosc;
		nE = wartosc;

		file >> tekst;

		for (int i = 0; i < nN; i++) // wartosci x i y node
		{
			int number;
			file >> number >> tekst >> wartosc1 >> tekst >> wartosc2;
			Node node;
			node.x = wartosc1;
			node.y = wartosc2;
			ND.push_back(node);
		}

		file >> tekst >> tekst;

		for (int i = 0; i < nE; i++) // ID elementow
		{
			int number;
			file >> number >> tekst >> war1 >> tekst >> war2 >> tekst >> war3 >> tekst >> war4;
			Element* element = new Element;
			element->ID[0] = war1;
			element->ID[1] = war2;
			element->ID[2] = war3;
			element->ID[3] = war4;
			EL.push_back(*element);
		}

		file >> tekst;
		int wartosc5;

		while (!file.eof()) // nody z flaga BC = 1 
		{
			file >> wartosc5 >> tekst;
			ND[wartosc5 - 1].BC = true;
		}

		file.close();
	}
};

struct GlobalData
{
	int simulationTime;
	int simulationStepTime;
	int conductivity;
	int alfa;
	int tot;
	int initialTemp;
	int density;
	int specificHeat;

	void readGlobalData()
	{
		ifstream file("Test3_31_31_kwadrat.txt");

		string tekst;
		int wartosc;

		for (int i = 0; i < 8; i++)
		{
			file >> tekst >> wartosc;

			if (tekst == "SimulationTime")
			{
				simulationTime = wartosc;
			}
			else if (tekst == "SimulationStepTime")
			{
				simulationStepTime = wartosc;
			}
			else if (tekst == "Conductivity")
			{
				conductivity = wartosc;
			}
			else if (tekst == "Alfa")
			{
				alfa = wartosc;
			}
			else if (tekst == "Tot")
			{
				tot = wartosc;
			}
			else if (tekst == "InitialTemp")
			{
				initialTemp = wartosc;
			}
			else if (tekst == "Density")
			{
				density = wartosc;
			}
			else if (tekst == "SpecificHeat")
			{
				specificHeat = wartosc;
			}
		}
		file.close();
	}
};

struct Element4 // dwupunktowy schemat calkowania
{
	int size = 4;
	double eta[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };
	double ksi[4] = { -1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) };
	double tabksi[4][4];
	double tabeta[4][4];
	double N[4][4];

	Element4()
	{
		for (int i = 0; i < size; i++)
		{
			tabksi[i][0] = -0.25 * (1 - eta[i]);
			tabksi[i][1] = 0.25 * (1 - eta[i]);
			tabksi[i][2] = 0.25 * (1 + eta[i]);
			tabksi[i][3] = -0.25 * (1 + eta[i]);
		}

		for (int i = 0; i < size; i++)
		{
			tabeta[i][0] = -0.25 * (1 - ksi[i]);
			tabeta[i][1] = -0.25 * (1 + ksi[i]);
			tabeta[i][2] = 0.25 * (1 + ksi[i]);
			tabeta[i][3] = 0.25 * (1 - ksi[i]);
		}

		for (int i = 0; i < size; i++)   // funkcje ksztaltu potrzebne do calkowania macierzyC
		{
			N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
			N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
			N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
			N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
		}
	}
};

struct Element9 // trzypunktowy schemat calkowania
{
	int size1 = 4;
	int size2 = 9;
	double eta[9] = { -sqrt((double)3 / 5), -sqrt((double)3 / 5), -sqrt((double)3 / 5), 0, 0, 0, sqrt((double)3 / 5), sqrt((double)3 / 5), sqrt((double)3 / 5) };
	double ksi[9] = { -sqrt((double)3 / 5), 0, sqrt((double)3 / 5), -sqrt((double)3 / 5), 0, sqrt((double)3 / 5), -sqrt((double)3 / 5), 0, sqrt((double)3 / 5) };
	double tabksi[9][4];
	double tabeta[9][4];
	double N[9][4];
	double waga1[9] = { 0.55555556, 0.88888889, 0.555555556, 0.55555556, 0.88888889, 0.555555556 , 0.55555556, 0.88888889, 0.555555556 }; // odpowiednie wagi do macierzyH
	double waga2[9] = { 0.55555556, 0.55555556, 0.55555556, 0.88888889, 0.88888889, 0.88888889, 0.55555556 , 0.55555556 , 0.55555556 };

	Element9()
	{
		for (int i = 0; i < size2; i++)
		{
			tabksi[i][0] = -0.25 * (1 - eta[i]);
			tabksi[i][1] = 0.25 * (1 - eta[i]);
			tabksi[i][2] = 0.25 * (1 + eta[i]);
			tabksi[i][3] = -0.25 * (1 + eta[i]);
		}

		for (int i = 0; i < size2; i++)
		{
			tabeta[i][0] = -0.25 * (1 - ksi[i]);
			tabeta[i][1] = -0.25 * (1 + ksi[i]);
			tabeta[i][2] = 0.25 * (1 + ksi[i]);
			tabeta[i][3] = 0.25 * (1 - ksi[i]);
		}

		for (int i = 0; i < size2; i++)     // // funkcje ksztaltu potrzebne do calkowania macierzyC
		{
			N[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
			N[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
			N[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
			N[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
		}
	}
};

struct Boki4 // struktura bokow dla dwupunktowego schematu calkowania
{
	double** N;

	Boki4()
	{
		double** punktCalkowania = new double* [8];

		for (int i = 0; i < 8; i++)
		{
			punktCalkowania[i] = new double[2];
		}

		punktCalkowania[0][0] = (-1 / sqrt(3));
		punktCalkowania[0][1] = -1;
		punktCalkowania[1][0] = 1 / sqrt(3);
		punktCalkowania[1][1] = -1;

		punktCalkowania[2][0] = 1;
		punktCalkowania[2][1] = (-1 / sqrt(3));
		punktCalkowania[3][0] = 1;
		punktCalkowania[3][1] = 1 / sqrt(3);

		punktCalkowania[4][0] = (-1 / sqrt(3));
		punktCalkowania[4][1] = 1;
		punktCalkowania[5][0] = 1 / sqrt(3);
		punktCalkowania[5][1] = 1;

		punktCalkowania[6][0] = -1;
		punktCalkowania[6][1] = (-1 / sqrt(3));
		punktCalkowania[7][0] = -1;
		punktCalkowania[7][1] = (1 / sqrt(3));

		N = new double* [8];

		for (int i = 0; i < 8; i++)
		{
			N[i] = new double[4];
		}

		for (int i = 0; i < 8; i++)    // funkcje ksztaltu do macierzyHBC
		{
			N[i][0] = 0.25 * (1 - punktCalkowania[i][0]) * (1 - punktCalkowania[i][1]);
			N[i][1] = 0.25 * (1 + punktCalkowania[i][0]) * (1 - punktCalkowania[i][1]);
			N[i][2] = 0.25 * (1 + punktCalkowania[i][0]) * (1 + punktCalkowania[i][1]);
			N[i][3] = 0.25 * (1 - punktCalkowania[i][0]) * (1 + punktCalkowania[i][1]);
		}
	}
};

struct Boki9 // struktura bokow dla trzypunktowego schematu calkowania
{
	double waga[12] = { 0.55555556, 0.88888889, 0.555555556,  0.55555556, 0.88888889, 0.555555556,  0.55555556, 0.88888889, 0.555555556,  0.55555556, 0.88888889, 0.555555556 }; // odpowiednie wagi do macierzyHBC
	double** N;

	Boki9()
	{
		double** punktCalkowania = new double* [12];

		for (int i = 0; i < 12; i++)
		{
			punktCalkowania[i] = new double[2];
		}

		punktCalkowania[0][0] = -sqrt((double)3 / 5);
		punktCalkowania[0][1] = -1;
		punktCalkowania[1][0] = 0;
		punktCalkowania[1][1] = -1;
		punktCalkowania[2][0] = sqrt((double)3 / 5);
		punktCalkowania[2][1] = -1;

		punktCalkowania[3][0] = 1;
		punktCalkowania[3][1] = -sqrt((double)3 / 5);
		punktCalkowania[4][0] = 1;
		punktCalkowania[4][1] = 0;
		punktCalkowania[5][0] = 1;
		punktCalkowania[5][1] = sqrt((double)3 / 5);

		punktCalkowania[6][0] = sqrt((double)3 / 5);
		punktCalkowania[6][1] = 1;
		punktCalkowania[7][0] = 0;
		punktCalkowania[7][1] = 1;
		punktCalkowania[8][0] = -sqrt((double)3 / 5);
		punktCalkowania[8][1] = 1;

		punktCalkowania[9][0] = -1;
		punktCalkowania[9][1] = sqrt((double)3 / 5);
		punktCalkowania[10][0] = -1;
		punktCalkowania[10][1] = 0;
		punktCalkowania[11][0] = -1;
		punktCalkowania[11][1] = -sqrt((double)3 / 5);

		N = new double* [12];

		for (int i = 0; i < 12; i++)
		{
			N[i] = new double[4];
		}

		for (int i = 0; i < 12; i++)   // funkcje ksztlaltu do macierzyHBC
		{
			N[i][0] = 0.25 * (1 - punktCalkowania[i][0]) * (1 - punktCalkowania[i][1]);
			N[i][1] = 0.25 * (1 + punktCalkowania[i][0]) * (1 - punktCalkowania[i][1]);
			N[i][2] = 0.25 * (1 + punktCalkowania[i][0]) * (1 + punktCalkowania[i][1]);
			N[i][3] = 0.25 * (1 - punktCalkowania[i][0]) * (1 + punktCalkowania[i][1]);
		}
	}
};

double** wyliczanieMacierzyH(Element4 elem4, Element9 elem9, Boki4 boki4, Boki9 boki9, int pc, double x[4], double y[4], bool bc[4], GlobalData globaldata)
{
	int size = pc * pc;
	double** macierzJakobiego = new double* [size];

	for (int i = 0; i < size; i++)
	{
		macierzJakobiego[i] = new double[4];
	}

	for (int i = 0; i < size; i++)
	{
		double dxdksi = 0;
		double dxdeta = 0;
		double dydksi = 0;
		double dydeta = 0;

		for (int j = 0; j < 4; j++)
		{
			if (pc == 2)
			{
				dxdksi = dxdksi + elem4.tabksi[i][j] * x[j];
				dxdeta = dxdeta + elem4.tabeta[i][j] * x[j];
				dydksi = dydksi + elem4.tabksi[i][j] * y[j];
				dydeta = dydeta + elem4.tabeta[i][j] * y[j];
			}
			else if (pc == 3)
			{
				dxdksi = dxdksi + elem9.tabksi[i][j] * x[j];
				dxdeta = dxdeta + elem9.tabeta[i][j] * x[j];
				dydksi = dydksi + elem9.tabksi[i][j] * y[j];
				dydeta = dydeta + elem9.tabeta[i][j] * y[j];
			}
		}

		macierzJakobiego[i][0] = dxdksi;
		macierzJakobiego[i][1] = dydksi;
		macierzJakobiego[i][2] = dxdeta;
		macierzJakobiego[i][3] = dydeta;
	}

	double* detJ = new double[size];
	for (int i = 0; i < size; i++)
	{
		detJ[i] = macierzJakobiego[i][0] * macierzJakobiego[i][3] - macierzJakobiego[i][1] * macierzJakobiego[i][2];
	}

	for (int i = 0; i < size; i++)
	{
		double temp = macierzJakobiego[i][0];
		macierzJakobiego[i][0] = macierzJakobiego[i][3];
		macierzJakobiego[i][1] = macierzJakobiego[i][1] * (-1);
		macierzJakobiego[i][2] = macierzJakobiego[i][2] * (-1);
		macierzJakobiego[i][3] = temp;
	}

	double* odwroconydetJ = new double[size];
	for (int i = 0; i < size; i++)
	{
		odwroconydetJ[i] = 1 / detJ[i];
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzJakobiego[i][j] = macierzJakobiego[i][j] * odwroconydetJ[i];
		}
	}

	double** pochodneX = new double* [size];
	double** pochodneY = new double* [size];

	for (int i = 0; i < size; i++)
	{
		pochodneX[i] = new double[4];
		pochodneY[i] = new double[4];
	}

	if (pc == 2)
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				pochodneX[i][j] = macierzJakobiego[i][0] * elem4.tabksi[i][j] + macierzJakobiego[i][1] * elem4.tabeta[i][j];
			}
		}
	}

	else if (pc == 3)
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				pochodneX[i][j] = macierzJakobiego[i][0] * elem9.tabksi[i][j] + macierzJakobiego[i][1] * elem9.tabeta[i][j];
			}
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (pc == 2)
			{
				pochodneY[i][j] = macierzJakobiego[i][2] * elem4.tabksi[i][j] + macierzJakobiego[i][3] * elem4.tabeta[i][j];
			}
			else if (pc == 3)
			{
				pochodneY[i][j] = macierzJakobiego[i][2] * elem9.tabksi[i][j] + macierzJakobiego[i][3] * elem9.tabeta[i][j];
			}
		}
	}

	// MACIERZ H 
	double** macierzH = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		macierzH[i] = new double[4];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzH[i][j] = 0;
		}
	}

	if (pc == 2)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[0][i] * pochodneX[0][j] + pochodneY[0][i] * pochodneY[0][j]) * detJ[0]);
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[1][i] * pochodneX[1][j] + pochodneY[1][i] * pochodneY[1][j]) * detJ[1]);
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[2][i] * pochodneX[2][j] + pochodneY[2][i] * pochodneY[2][j]) * detJ[2]);
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[3][i] * pochodneX[3][j] + pochodneY[3][i] * pochodneY[3][j]) * detJ[3]);
			}
		}

	}
	else if (pc == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[0][i] * pochodneX[0][j] + pochodneY[0][i] * pochodneY[0][j]) * detJ[0]) * elem9.waga1[0] * elem9.waga2[0];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[1][i] * pochodneX[1][j] + pochodneY[1][i] * pochodneY[1][j]) * detJ[1]) * elem9.waga1[1] * elem9.waga2[1];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[2][i] * pochodneX[2][j] + pochodneY[2][i] * pochodneY[2][j]) * detJ[2]) * elem9.waga1[2] * elem9.waga2[2];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[3][i] * pochodneX[3][j] + pochodneY[3][i] * pochodneY[3][j]) * detJ[3]) * elem9.waga1[3] * elem9.waga2[3];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[4][i] * pochodneX[4][j] + pochodneY[4][i] * pochodneY[4][j]) * detJ[4]) * elem9.waga1[4] * elem9.waga2[4];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[5][i] * pochodneX[5][j] + pochodneY[5][i] * pochodneY[5][j]) * detJ[5]) * elem9.waga1[5] * elem9.waga2[5];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[6][i] * pochodneX[6][j] + pochodneY[6][i] * pochodneY[6][j]) * detJ[6]) * elem9.waga1[6] * elem9.waga2[6];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[7][i] * pochodneX[7][j] + pochodneY[7][i] * pochodneY[7][j]) * detJ[7]) * elem9.waga1[7] * elem9.waga2[7];
				macierzH[i][j] += (globaldata.conductivity * (pochodneX[8][i] * pochodneX[8][j] + pochodneY[8][i] * pochodneY[8][j]) * detJ[8]) * elem9.waga1[8] * elem9.waga2[8];

			}
		}
	}

	// MACIERZ HBC
	double detJ2[4];
	for (int i = 0; i < 4; i++)
	{
		if (i == 3)
		{
			detJ2[i] = sqrt(pow(x[i] - x[0], 2) + pow(y[i] - y[0], 2)) / 2;
		}
		else
		{
			detJ2[i] = sqrt(pow(x[i] - x[i + 1], 2) + pow(y[i] - y[i + 1], 2)) / 2;
		}
	}

	double** macierzHBC = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		macierzHBC[i] = new double[4];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzHBC[i][j] = 0;
		}
	}

	if (pc == 2)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (bc[0] == true && bc[1] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * (boki4.N[0][i] * boki4.N[0][j] + boki4.N[1][i] * boki4.N[1][j]) * detJ2[0];
				}
				if (bc[1] == true && bc[2] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * (boki4.N[2][i] * boki4.N[2][j] + boki4.N[3][i] * boki4.N[3][j]) * detJ2[1];
				}
				if (bc[2] == true && bc[3] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * (boki4.N[4][i] * boki4.N[4][j] + boki4.N[5][i] * boki4.N[5][j]) * detJ2[2];
				}
				if (bc[3] == true && bc[0] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * (boki4.N[6][i] * boki4.N[6][j] + boki4.N[7][i] * boki4.N[7][j]) * detJ2[3];
				}
			}
		}
	}

	else if (pc == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (bc[0] == true && bc[1] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * ((boki9.N[0][i] * boki9.N[0][j]) * boki9.waga[0] + (boki9.N[1][i] * boki9.N[1][j]) * boki9.waga[1] + (boki9.N[2][i] * boki9.N[2][j]) * boki9.waga[2]) * detJ2[0];
				}
				if (bc[1] == true && bc[2] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * ((boki9.N[3][i] * boki9.N[3][j]) * boki9.waga[3] + (boki9.N[4][i] * boki9.N[4][j]) * boki9.waga[4] + (boki9.N[5][i] * boki9.N[5][j]) * boki9.waga[5]) * detJ2[1];
				}
				if (bc[2] == true && bc[3] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * ((boki9.N[6][i] * boki9.N[6][j]) * boki9.waga[6] + (boki9.N[7][i] * boki9.N[7][j]) * boki9.waga[7] + (boki9.N[8][i] * boki9.N[8][j]) * boki9.waga[8]) * detJ2[2];
				}
				if (bc[3] == true && bc[0] == true)
				{
					macierzHBC[i][j] += globaldata.alfa * ((boki9.N[9][i] * boki9.N[9][j]) * boki9.waga[9] + (boki9.N[10][i] * boki9.N[10][j]) * boki9.waga[10] + (boki9.N[11][i] * boki9.N[11][j]) * boki9.waga[11]) * detJ2[3];
				}
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzH[i][j] += macierzHBC[i][j];
		}
	}

	return macierzH;
}

double* liczenieP(Boki4 boki4, Boki9 boki9, int pc, double x[4], double y[4], bool bc[4], GlobalData globaldata)
{
	double* P = new double[4];

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
	}

	double detJ2[4];
	for (int i = 0; i < 4; i++)
	{
		if (i == 3)
		{
			detJ2[i] = sqrt(pow(x[i] - x[0], 2) + pow(y[i] - y[0], 2)) / 2;
		}
		else
		{
			detJ2[i] = sqrt(pow(x[i] - x[i + 1], 2) + pow(y[i] - y[i + 1], 2)) / 2;
		}
	}

	if (pc == 2)
	{
		for (int i = 0; i < 4; i++)
		{
			if (bc[0] == true && bc[1] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki4.N[0][i] + boki4.N[1][i]) * detJ2[0];
			}
			if (bc[1] == true && bc[2] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki4.N[2][i] + boki4.N[3][i]) * detJ2[1];
			}
			if (bc[2] == true && bc[3] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki4.N[4][i] + boki4.N[5][i]) * detJ2[2];
			}
			if (bc[3] == true && bc[0] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki4.N[6][i] + boki4.N[7][i]) * detJ2[3];
			}
		}
	}

	else if (pc == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			if (bc[0] == true && bc[1] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki9.N[0][i] * boki9.waga[0] + boki9.N[1][i] * boki9.waga[1] + boki9.N[2][i] * boki9.waga[2]) * detJ2[0];
			}
			if (bc[1] == true && bc[2] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki9.N[3][i] * boki9.waga[3] + boki9.N[4][i] * boki9.waga[4] + boki9.N[5][i] * boki9.waga[5]) * detJ2[1];
			}
			if (bc[2] == true && bc[3] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki9.N[6][i] * boki9.waga[6] + boki9.N[7][i] * boki9.waga[7] + boki9.N[8][i] * boki9.waga[8]) * detJ2[2];
			}
			if (bc[3] == true && bc[0] == true)
			{
				P[i] += globaldata.alfa * globaldata.tot * (boki9.N[9][i] * boki9.waga[9] + boki9.N[10][i] * boki9.waga[10] + boki9.N[11][i] * boki9.waga[11]) * detJ2[3];
			}
		}
	}

	return P;
}

double** liczenieC(Element4 elem4, Element9 elem9, int pc, double x[4], double y[4], GlobalData globaldata)
{
	int size = pc * pc;
	double** macierzJakobiego = new double* [size];

	for (int i = 0; i < size; i++)
	{
		macierzJakobiego[i] = new double[4];
	}

	for (int i = 0; i < size; i++)
	{
		double dxdksi = 0;
		double dxdeta = 0;
		double dydksi = 0;
		double dydeta = 0;

		for (int j = 0; j < 4; j++)
		{
			if (pc == 2)
			{
				dxdksi = dxdksi + elem4.tabksi[i][j] * x[j];
				dxdeta = dxdeta + elem4.tabeta[i][j] * x[j];
				dydksi = dydksi + elem4.tabksi[i][j] * y[j];
				dydeta = dydeta + elem4.tabeta[i][j] * y[j];
			}
			else if (pc == 3)
			{
				dxdksi = dxdksi + elem9.tabksi[i][j] * x[j];
				dxdeta = dxdeta + elem9.tabeta[i][j] * x[j];
				dydksi = dydksi + elem9.tabksi[i][j] * y[j];
				dydeta = dydeta + elem9.tabeta[i][j] * y[j];
			}
		}

		macierzJakobiego[i][0] = dxdksi;
		macierzJakobiego[i][1] = dydksi;
		macierzJakobiego[i][2] = dxdeta;
		macierzJakobiego[i][3] = dydeta;
	}

	double* detJ = new double[size];

	for (int i = 0; i < size; i++)
	{
		detJ[i] = macierzJakobiego[i][0] * macierzJakobiego[i][3] - macierzJakobiego[i][1] * macierzJakobiego[i][2];
	}

	double** macierzC = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		macierzC[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			macierzC[i][j] = 0;
		}
	}

	if (pc == 2)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem4.N[0][i] * elem4.N[0][j]) * detJ[0];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem4.N[1][i] * elem4.N[1][j]) * detJ[1];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem4.N[2][i] * elem4.N[2][j]) * detJ[2];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem4.N[3][i] * elem4.N[3][j]) * detJ[3];
			}
		}
	}

	else if (pc == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[0][i] * elem9.N[0][j]) * elem9.waga1[0] * elem9.waga2[0] * detJ[0];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[1][i] * elem9.N[1][j]) * elem9.waga1[1] * elem9.waga2[1] * detJ[1];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[2][i] * elem9.N[2][j]) * elem9.waga1[2] * elem9.waga2[2] * detJ[2];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[3][i] * elem9.N[3][j]) * elem9.waga1[3] * elem9.waga2[3] * detJ[3];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[4][i] * elem9.N[4][j]) * elem9.waga1[4] * elem9.waga2[4] * detJ[4];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[5][i] * elem9.N[5][j]) * elem9.waga1[5] * elem9.waga2[5] * detJ[5];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[6][i] * elem9.N[6][j]) * elem9.waga1[6] * elem9.waga2[6] * detJ[6];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[7][i] * elem9.N[7][j]) * elem9.waga1[7] * elem9.waga2[7] * detJ[7];
				macierzC[i][j] += globaldata.specificHeat * globaldata.density * (elem9.N[8][i] * elem9.N[8][j]) * elem9.waga1[8] * elem9.waga2[8] * detJ[8];
			}
		}
	}

	return macierzC;
}

void agregacja(double** macierzGlobalnaH, double** macierzLokalnaH, int id[4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzGlobalnaH[id[i] - 1][id[j] - 1] += macierzLokalnaH[i][j];
		}
	}
}

void agregacjaWektora(double* wektorGlobalny, double* wektorLokalny, int id[4])
{
	for (int i = 0; i < 4; i++)
	{
		wektorGlobalny[id[i] - 1] += wektorLokalny[i];
	}
}

double* ukladRownan(double** macierz, double* wektor, int n)
{
	int i;
	int j;
	int k;

	double* x = new double[n];
	double** a = new double* [n];

	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n + 1];
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i][j] = macierz[i][j];
		}
	}

	for (i = 0; i < n; i++)
	{
		a[i][n] = wektor[i];
	}

	for (i = 0; i < n; i++)
	{
		for (k = i + 1; k < n; k++)
		{
			if (abs(a[i][i]) < abs(a[k][i]))
			{
				for (j = 0; j <= n; j++)
				{
					double temp = a[i][j];
					a[i][j] = a[k][j];
					a[k][j] = temp;
				}
			}
		}
	}

	for (i = 0; i < n - 1; i++)
	{
		for (k = i + 1; k < n; k++)
		{
			double t = a[k][i] / a[i][i];
			for (j = 0; j <= n; j++)
			{
				a[k][j] = a[k][j] - t * a[i][j];
			}
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		x[i] = a[i][n];
		for (j = i + 1; j < n; j++)
		{
			if (j != i)
			{
				x[i] = x[i] - a[i][j] * x[j];
			}
		}
		x[i] = x[i] / a[i][i];
	}

	// Dodatkowo znajdowanie minimum i maksimum
	double MIN;
	double MAX;
	int p = 2;

	if (x[0] > x[1])
	{
		MIN = x[1];
		MAX = x[0];
	}
	else
	{
		MIN = x[0];
		MAX = x[1];
	}

	while (p + 2 <= n)
	{
		if (x[p] > x[p + 1])
		{
			if (x[p] > MAX)
				MAX = x[i];
			if (x[p + 1] < MIN)
				MIN = x[p + 1];
		}
		else
		{
			if (x[p + 1] > MAX)
				MAX = x[p + 1];
			if (x[p] < MIN)
				MIN = x[p];
		}
		p += 2;
	}
	if (n % 2 == 1)
	{
		if (x[p] > MAX) MAX = x[p];
		if (x[p] < MIN) MIN = x[p];
	}

	cout << endl;
	cout << "MIN:" << " " << MIN << endl;;
	cout << "MAX: " << MAX << endl;;

	return x;
}

int main()
{
	Grid grid;
	grid.readGrid();

	GlobalData globalData;
	globalData.readGlobalData();

	Element4 el1;
	Element9 el2;
	Boki4 boki4;
	Boki9 boki9;

	int rozmiar = grid.nN;

	// tworzenie macierzy
	double** macierzH = new double* [rozmiar];
	double** macierzC = new double* [rozmiar];
	double** macierzHmacierzC = new double* [rozmiar];

	for (int i = 0; i < rozmiar; i++)
	{
		macierzH[i] = new double[rozmiar];
		macierzC[i] = new double[rozmiar];
		macierzHmacierzC[i] = new double[rozmiar];
	}

	for (int i = 0; i < rozmiar; i++)
	{
		for (int j = 0; j < rozmiar; j++)
		{
			macierzH[i][j] = 0;
			macierzC[i][j] = 0;
			macierzHmacierzC[i][j] = 0;
		}
	}

	// tworzenie wektorow
	double* wektorP = new double[rozmiar];
	double* finalwektorP = new double[rozmiar];
	double* wektort0 = new double[rozmiar];

	for (int i = 0; i < rozmiar; i++)
	{
		wektorP[i] = 0;
		finalwektorP[i] = 0;
		wektort0[i] = globalData.initialTemp;
	}

	int id[4];
	double x[4];
	double y[4];
	bool bc[4];

	// przejscie po elementach siatki
	for (int i = 0; i < grid.nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			id[j] = grid.EL[i].ID[j];
			x[j] = grid.ND[id[j] - 1].x;
			y[j] = grid.ND[id[j] - 1].y;
			bc[j] = grid.ND[id[j] - 1].BC;
		}

		// agregacja macierzH
		double** lokalnaMacierzH = wyliczanieMacierzyH(el1, el2, boki4, boki9, 3, x, y, bc, globalData);
		agregacja(macierzH, lokalnaMacierzH, id);

		// agregacja macierzC
		double** lokalnaMacierzC = liczenieC(el1, el2, 3, x, y, globalData);
		agregacja(macierzC, lokalnaMacierzC, id);

		// agregacja wektorP
		double* lokalnyWektorP = liczenieP(boki4, boki9, 3, x, y, bc, globalData);
		agregacjaWektora(wektorP, lokalnyWektorP, id);
	}

	for (int i = 0; i < rozmiar; i++)
	{
		for (int j = 0; j < rozmiar; j++)
		{
			macierzHmacierzC[i][j] = macierzH[i][j] + (macierzC[i][j] / globalData.simulationStepTime);
		}
	}

	for (int t = 0; t < 10; t++)
	{
		for (int i = 0; i < rozmiar; i++)
		{
			for (int j = 0; j < rozmiar; j++)
			{
				finalwektorP[i] += (macierzC[i][j] / globalData.simulationStepTime) * wektort0[j];
			}
		}

		for (int i = 0; i < rozmiar; i++)
		{
			finalwektorP[i] += wektorP[i];
		}

		// OSTATECZNE ROZWIAZANIE, UZYSKANE TEMPERATURY 
		cout << "Czas: " << (t + 1) * globalData.simulationStepTime;

		wektort0 = ukladRownan(macierzHmacierzC, finalwektorP, grid.nN);
		cout << endl;

		for (int i = 0; i < grid.nN; i++)
		{
			finalwektorP[i] = 0;
		}
	}

	return 0;
}