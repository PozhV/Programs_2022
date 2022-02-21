/*
 * DSMC.h
 *
 *  Created on: 6 апр. 2021 г.
 *      Author: User
 */
#pragma once
#include "utils.h"
#include "vector3d.h"
#include"Ascillator.h"
#include "Randoms.h"
#include<fstream>
#include <iostream>
#include <limits>
#include <string>
enum { NO_COLLISION, LOWER_WALL, UPPER_WALL, LEFT_WALL, RIGHT_WALL, 
	ASCILLATOR_UP, ASCILLATOR_DOWN, ASCILLATOR_LEFT, ASCILLATOR_RIGHT};


struct Molecule
{
  Vector3d r;
  Vector3d v;
};

struct Cell
{
  int molcount;
  std::vector<Molecule> molecules;
  double cvrMax;
  double collisionRemainder;
};
struct Parameters
{
	const double L;
	const double H_L;
	const double l_L;
	const double A_H;
	const double h_H;
	const double Kn;
	const double u_c;

};
class DSMC
{
public:
	DSMC(Parameters P, const double T_, int y_cellcount, double dtMult, int numPartPerCell, int CellsAverNum);

  void Step();

  int CalculateAveragedVtk();
  int ExportToVtk();

  int CalculateTotalMolCount();
  std::pair<int, double> FindNearestCollisionSurface(const Molecule& molecule, double dt, double t) const;
  
private:
  void MoveParticles();
  void UpdateMoleculesInCells();
  void CollideParticles();
  void InitParticles();
  void Collision(Molecule& mol1, Molecule& mol2);
  double MolcountToNumDensity(int N);
  double   CalculateTranslationalTemperature(const Cell& cell);
  Vector3d CalculateMeanVelocity(const Cell& cell);

private:
  // Geometric parameters
  double cellSize;
  double cellVolume;
  double cellSquare;
  Vector3d domainSize;
  double T;
  double numDensity;
  double p;
  // Number of cells
  int xcellcount;
  int ycellcount;
  int zcellcount;
  double CellsCalcNum;
  int xcellcountforVtk;
  int ycellcountforVtk;
  double betta;
  double gammaFunc;
public:
  double dt;
  Ascillator as;
  double Period;
  //Params for stats
private:
  double Time = 0.;
  int counter = 0;
  int ExportNum = 0;
  int VTKcount = 0;
  int SubcellCounts[MaxSubcell][MaxSubcell];
  int SubCellMol[MaxSubcell][MaxSubcell][100];
  //Termodinamics params
  double MN;
  std::ofstream file;

  //Microparams
  std::vector< std::vector< std::vector<Vector3d> > > meanV;
  std::vector< std::vector< std::vector<double> > >   Temp;
  std::vector< std::vector< std::vector<double> > >   pressure;
  std::vector< std::vector< std::vector<int> > >      molPerCell;
  std::vector< std::vector< std::vector<Cell> > >   cells;
  double press = 0;
};
