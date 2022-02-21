/*
 * DSMC.cpp
 *
 *  Created on: 13 февр. 2021 г.
 *      Author: Asus
 */
#include "DSMC.h"
#include<algorithm>
using namespace std;
DSMC::DSMC(Parameters P, const double T_, int y_cellcount, double dtMult, int NPC, int CellsAverNum) :
    cellSize(P.L* P.H_L / y_cellcount),
    cellSquare(pow(cellSize, 2)),
    cellVolume(pow(cellSize, 3)),
    domainSize({ P.L, P.L * P.H_L, cellSize }),
    T(T_),
    numDensity(1. / (P.Kn * domainSize.y * sqrt(2) * SIGMA_REF)),
    p(numDensity* BOLTZMANN_K* T),
    xcellcount((int)(domainSize.x / cellSize)),
    ycellcount((int)(domainSize.y / cellSize)),
    zcellcount((int)(domainSize.z / cellSize)),
    CellsCalcNum(CellsAverNum),
    xcellcountforVtk(int(xcellcount / CellsCalcNum)),
    ycellcountforVtk(int(ycellcount / CellsCalcNum)),
    betta(1. / sqrt(3. * R_GAS_CONSTANT * T)),
    gammaFunc(Gamma(2.5 - OMEGA)),
    dt(dtMult* (cellSize* betta)),
    MN(numDensity* cellVolume / NPC),
    as(M_PI* P.u_c / (betta * P.A_H * domainSize.y), 0., domainSize.y* P.A_H, P.l_L* domainSize.x, P.h_H* domainSize.y,
        (domainSize.x - P.l_L * domainSize.x) / 2., (domainSize.x + P.l_L * domainSize.x) / 2.),
    Period(2 * M_PI / as.Get_om()),
    file("pressure_" + to_string(P.l_L) + "_" + to_string(P.Kn) + ".txt"),
    meanV(xcellcountforVtk, vector<vector<Vector3d> >(ycellcountforVtk, vector<Vector3d>(zcellcount, Vector3d{ 0, 0, 0 }))),
    Temp(xcellcountforVtk, vector<vector<double> >(ycellcountforVtk, vector<double>(zcellcount, 0))),
    pressure(xcellcountforVtk, vector<vector<double> >(ycellcountforVtk, vector<double>(zcellcount, 0))),
    molPerCell(xcellcountforVtk, vector<vector<int> >(ycellcountforVtk, vector<int>(zcellcount, 0))),
    cells(xcellcount, vector<vector<Cell> >(ycellcount, vector<Cell>(zcellcount, Cell())))
{
  cout<<1./betta<<endl;
  cout << cellSize << endl;
  cout << "Dens: " << numDensity << endl;
  cout<<"dt: " << dt << endl;
  cout << "Period: " << Period << endl;
  cout << as.Get_om() << endl;
  for (int i = 0; i < xcellcount; i++)
    for (int j = 0; j < ycellcount; j++)
      for (int k = 0; k < zcellcount; k++)
      {
        cells[i][j][k].molcount = 0;
        cells[i][j][k].cvrMax   = (SIGMA_REF / betta) * gammaFunc;
      }
InitParticles();
}
void DSMC::InitParticles()
{
    double ed1 = as.Get_edges().first;
    double ed2 = as.Get_edges().second;
    double y0 = as.Calc_pos(0);
	double n_reserv   = p / (BOLTZMANN_K *T);
	int    N       = int(n_reserv * cellVolume / MN);
	  for (int i = 0; i < xcellcount; i++)
	  {
	    for (int j = 0; j < ycellcount; j++)
	    {
	      for (int k = 0; k < zcellcount; k++)
	      {
              if ((i >= floor(ed1 / cellSize) && i <= ceil(ed2 / cellSize)) && j >= floor((y0 - as.h) / cellSize) && j <= ceil((y0 + as.h)/cellSize))
                  continue;
	        Cell& cell = cells[i][j][k];
	        cell.molecules.resize(N);
	        cell.molcount = N;
	        for (int w = 0; w < cell.molcount; w++)
	        {
	          cell.molecules[w].r    = cellSize * Vector3d{ i + RandN(), j + RandN(), k + RandN() };
	          cell.molecules[w].v    = { RandG(betta), RandG(betta), RandG(betta) };
            }
	      }
	    }
	  }
}


void DSMC::Step()
{
    press = 0;
    MoveParticles();
    UpdateMoleculesInCells();
    CollideParticles();
    press *= MOLECULAR_MASS * MN / dt;
    file <<(Time / Period)<<" "<< press<<endl;
    Time += dt;
    counter++;
}

std::pair<int, double> DSMC::FindNearestCollisionSurface(const Molecule& molecule, const double _dt, const double t) const
{
  const double x = molecule.r.x;
  const double y = molecule.r.y;
  const double vx = molecule.v.x;
  const double vy = molecule.v.y;
  const double ed1 = as.Get_edges().first;
  const double ed2 = as.Get_edges().second;
  int bound_indx;
  int bound_indy;
  double Tx, Ty;
  double t1;
  if (vx > 0)
  {
      t1 = (ed1 - x) / vx;
      Tx = abs((domainSize.x - x) / vx) - DBL_MIN;
      bound_indx = 4;
      if (t1 < 0)
      {
          t1 = 0;
      }
      else if (abs(as.F(t + t1, y + vy * t1)) <= as.h)
      {
          if (t1 <= _dt)
              return { 7, t1 };
          else
              return { 0, _dt };
      }

  }
  else
  {
      t1 = (ed2 - x) / vx;
      Tx = abs(x / vx) - DBL_MIN;
      bound_indx = 3;
      if (t1 < 0)
      {
          t1 = 0;
      }
      else if (abs(as.F(t + t1, y + vy * t1)) <= as.h)
      {
          if (t1 <= _dt)
              return { 8, t1 };
          else
              return { 0, _dt };
      }
  }
  if (vy > 0)
  {
      bound_indy = 2;
      Ty = abs((y - domainSize.y) / vy) - DBL_MIN;
  }
  else
  {
      bound_indy = 1;
      Ty = abs(y / vy) - DBL_MIN;
  }
  const double Tmin = min(min(Tx, Ty), _dt);
  if((vx > 0 && x <= ed2 || vx < 0 && x >= ed1) && Tmin > t1)
  {
      t1 += t;
      double t2 = t1;
      int indef = 0;
      int i = 0;
      do
      {
          double t3 = t2;
          t2 = as.extrem(t2, vy, t + _dt - t2, i);
          bool f1 = as.F_up(t1, y + vy * (t1 - t));
          bool f2 = as.F_up(t2, y + vy * (t2 - t));
          bool f3 = as.F_down(t1, y + vy * (t1 - t));
          bool f4 = as.F_down(t2, y + vy * (t2 - t));
          
          if ((f1^f2) && (f3^f4))
          {
              if (as.dif_F((t1 + t2) / 2., vy) > 0)
              {
                  indef = 1;
                  break;
              }
              else
              {
                  indef = 2;
                  break;
              }
          }
          if (f1^f2)
          {
              indef = 1;
              break;
          }
          if (f3^f4)
          {
              indef = 2;
              break;
          }
          t1 = t2;
          i++;
          if (i > 1000)
          {
              cout << t2 - t - _dt << "  " << t2 - t3 << endl;
          }
      } while (t2 < t + _dt - DBL_MIN);
      if (indef > 0)
      {
          bool (Ascillator::*Func)(double t, double y) const = &Ascillator::F_up;
          double d = as.h;
          if (indef == 2)
          {
              d = -as.h;
              Func = &Ascillator::F_down;
          }
          double eps = 1e-5 * as.maxV * dt;
          double med = (t1 + t2) / 2.0;
          double dif = abs(as.F(med, y + vy * (med - t)) + d);
          int i = 0;
          while (dif > eps)
          {
              bool f1 = (as.*Func)(t1, y + vy * (t1 - t));
              bool f2 = (as.*Func)(t2, y + vy * (t2 - t));
              bool f_med = (as.*Func)(med, y + vy * (med - t));

              if (f1 ^ f_med)
              {
                  t2 = med;
              }
              else if (f_med ^ f2)
              {
                  t1 = med;
              }
              
              med = (t1 + t2) / 2.;
              dif = abs(as.F(med, y + vy * (med - t)) + d);
              i++;
              if (i > 1000)
                  cout << dif << endl;
          }
          
          if ((x + vx * (med - t) <= ed2) && (x + vx * (med - t) >= ed1))
          {
              if (indef == 1)
                  return { 5, med - t};
              else
                  return { 6, med - t};
          }
      }
  }
  if (Tmin < _dt)
  {
      if (Tx < Ty)
      {
          return{ bound_indx, Tx };
      }
      else 
      {
          return{ bound_indy, Ty };
      }
  }
  else
  {
      return { 0, _dt };
  }

}


void DSMC::MoveParticles()
{
  for (int i = 0; i < xcellcount; i++)
    for (int j = 0; j < ycellcount; j++)
      for (int k = 0; k < zcellcount; k++)
        for (int n = 0; n < cells[i][j][k].molcount; n++)
        {
          Molecule& molecule = cells[i][j][k].molecules[n]; // Remembering this specific molecule inside cell (this is just for code-reading simplicity)

          double dtTmp = 0;
          while (dtTmp < dt) 
          {
            std::pair <int, double> res = FindNearestCollisionSurface(molecule, dt - dtTmp, Time + dtTmp);
            int    boundaryType = res.first;
            double dt1          = res.second; // time until collision with the nearest boundary
            double prev_x = molecule.r.x;
            molecule.r += molecule.v * dt1; // Advance particle position until collision with boundary

            if (boundaryType == NO_COLLISION)
            {
              // Applying periodic conditions for Z
                if (molecule.r.z > domainSize.z)
                {
                    molecule.r.z -= floor( molecule.r.z / domainSize.z ) * domainSize.z;
                }
                else if (molecule.r.z < 0)
                {
                    molecule.r.z += (floor( abs(molecule.r.z) / domainSize.z ) + 1) * domainSize.z;
                }
            }
            else if (boundaryType == LOWER_WALL)
            {
              molecule.v    = { RandG(betta), RandGFL(betta), RandG(betta) }; // Diffuse reflection:
            }
            else if (boundaryType == UPPER_WALL)
            {
            	molecule.v    = { RandG(betta), (-1)*RandGFL(betta), RandG(betta) };
            }
            else if (boundaryType == LEFT_WALL)
            {
                molecule.v = { RandGFL(betta), RandG(betta), RandG(betta) };
            }
            else if (boundaryType == RIGHT_WALL)
            {
                molecule.v = { (- 1) * RandGFL(betta), RandG(betta), RandG(betta)};
            }
            else if (boundaryType == ASCILLATOR_UP)
            {
                press -= as.dif_F(Time + dtTmp, molecule.v.y);
                molecule.v = { RandG(betta), RandGFL(betta) + as.CurVelocity(Time + dtTmp), RandG(betta)};
                press += as.dif_F(Time + dtTmp, molecule.v.y);
            }
            else if (boundaryType == ASCILLATOR_DOWN)
            {
                press -= as.dif_F(Time + dtTmp, molecule.v.y);
                molecule.v = { RandG(betta), (-1) * RandGFL(betta) + as.CurVelocity(Time + dtTmp), RandG(betta) };
                press += as.dif_F(Time + dtTmp, molecule.v.y);
            }
            else if (boundaryType == ASCILLATOR_LEFT)
            {
                molecule.v = { (- 1) * RandGFL(betta), RandG(betta) + as.CurVelocity(Time + dtTmp), RandG(betta)};
            }
            else
            {
                molecule.v = { RandGFL(betta), RandG(betta) + as.CurVelocity(Time + dtTmp), RandG(betta) };
            }
            dtTmp += dt1;
          }
        }

  return;
}

void DSMC::CollideParticles()
{
  for (int i = 0; i < xcellcount; i++)
    for (int j = 0; j < ycellcount; j++)
      for (int k = 0; k < zcellcount; k++)
      {
        Cell& cell = cells[i][j][k];

        if (cell.molcount < MinMolCount) // Only if we have at least two molecules inside cell
          continue;

        double lambda = 1. / (sqrt(2) * SIGMA_HS_REF * MolcountToNumDensity(cell.molcount)); // mean free-path
        double lambda1 = lambda / 3;

        int    NumCollisionsTotal = 0;   // we need this only to calculate quality
        double meanDist           = 0;   // we need this only to calculate quality

        if ((cellSize > 2 * lambda1) && USE_ADAPTIVE_CELLS) // Condition to divide cell into subCells
        {

        }
        else
        {
        	double NumCollisionsDbl = (1/2.*cell.molcount * cell.molcount * cell.cvrMax * dt / cellVolume * MN) + cell.collisionRemainder;
          int    NumCollisionsInt = (int)(NumCollisionsDbl);
          cell.collisionRemainder = NumCollisionsDbl - NumCollisionsInt;
          NumCollisionsTotal = NumCollisionsInt;

          while (NumCollisionsInt > 0)
          {

            int n1 = (int)(RandU(0, cell.molcount));
            int n2 = (int)(RandU(0, cell.molcount));
            if (n1 == n2)
              continue;
            meanDist  += (cell.molecules[n1].r - cell.molecules[n2].r).length(); // Distance between molecules
            double cr  = (cell.molecules[n1].v - cell.molecules[n2].v).length();
            double cvr = cr * SIGMA_REF * (pow((VEL_REF / cr), (2.0 * OMEGA - 1)) / gammaFunc);
            if (cvr > RandN() * cell.cvrMax) // Condition to accept collision
              Collision(cell.molecules[n1], cell.molecules[n2]);
            cell.cvrMax = std::max(cvr, cell.cvrMax);
            NumCollisionsInt--;
          }
       }
      }

  return;
}
Vector3d DSMC::CalculateMeanVelocity(const Cell& cell)
{
  Vector3d meanV = { 0.0, 0.0, 0.0 };
  for (int n = 0; n < cell.molcount; n++)
    meanV += cell.molecules[n].v;

  meanV = meanV / std::max(cell.molcount, 1);

  return meanV;
}

double DSMC::CalculateTranslationalTemperature(const Cell& cell)
{
  Vector3d meanV = CalculateMeanVelocity(cell);

  double sumVV = 0;
  for (int n = 0; n < cell.molcount; n++)
    sumVV += (cell.molecules[n].v - meanV).norm();

  sumVV = sumVV / std::max(cell.molcount, 1);

  return sumVV / (3.0 * R_GAS_CONSTANT);
}

int DSMC::CalculateTotalMolCount()
{
  int totalMolCount = 0;

  for (int i = 0; i < xcellcount; i++)
    for (int j = 0; j < ycellcount; j++)
      for (int k = 0; k < zcellcount; k++)
        totalMolCount += cells[i][j][k].molcount;

  return totalMolCount;
}


double DSMC::MolcountToNumDensity(int N)
{
  return N * MN / cellVolume;
}


void DSMC::UpdateMoleculesInCells()
{
  for (int i = 0; i < xcellcount; i++)
    for (int j = 0; j < ycellcount; j++)
      for (int k = 0; k < zcellcount; k++)
      {
        Cell& curCell = cells[i][j][k];

        for (int n = 0; n < curCell.molcount; n++)
        {
          Molecule mol = curCell.molecules[n];

          int ci = (int)(floor(mol.r.x / cellSize));
          int cj = (int)(floor(mol.r.y / cellSize));
          int ck = (int)(floor(mol.r.z / cellSize));

          if (ci != i || cj != j || ck != k)
          {
            curCell.molecules[n] = curCell.molecules[curCell.molcount - 1]; // sets n-th molecule ot be equal the last
            curCell.molcount--;
            n--;

            if ((ci >= 0 && ci < xcellcount) &&
                (cj >= 0 && cj < ycellcount) &&
                (ck >= 0 && ck < zcellcount)) // In other case the molecule left the area through open boundary
            {
              Cell& newCell = cells[ci][cj][ck];

              if (newCell.molcount < int(newCell.molecules.size()))
                newCell.molecules[newCell.molcount] = mol;
              else
                newCell.molecules.push_back(mol);

              newCell.molcount++;
            }
          }
        }
      }

  return;
}

void DSMC::Collision(Molecule& mol1, Molecule& mol2)
{
  Vector3d vRel;
  vRel= mol1.v - mol2.v;
  double cr = vRel.length();

  double cosTheta = RandU(-1., 1.);
  double sinTheta = sqrt(1. - cosTheta * cosTheta);
  double fi       = RandU(0, 2 * M_PI); //

  Vector3d vr;
  vr.x = cr * cosTheta;
  vr.y = cr * sinTheta * cos(fi);
  vr.z = cr * sinTheta * sin(fi);

  // center-of-mass velocity
  Vector3d cm = (mol1.v + mol2.v) / 2.0;

  // Velocities after collision
  mol1.v = cm + vr / 2.0;
  mol2.v = cm - vr / 2.0;

  if ((mol1.v.length() > 1e+5) || (mol2.v.length() > 1e+5))
    printf("bad velocity\n");

  return;
}


int DSMC::CalculateAveragedVtk()
{
  VTKcount++;
  Vector3d tmpV(0,0,0);
  double   temp=0;
  double   numDensity=0;
  int aver_molcount=0, i1, j1;
  for (int i = 0; i < xcellcountforVtk; i++)
  {
	  i1=int(i*CellsCalcNum);
    for (int j = 0; j < ycellcountforVtk; j++)
    {
    	 j1=int(j*CellsCalcNum);
      for (int k = 0; k < zcellcount; k++)
      {
    	for(int k1=0; k1<CellsCalcNum; k1++)
    	{
        	for(int k2=0; k2<CellsCalcNum; k2++)
        	{
        		Cell& curCell = cells[i1+k1][j1+k2][k];
        		tmpV       += CalculateMeanVelocity(curCell);
        		temp       += CalculateTranslationalTemperature(curCell);
        		aver_molcount+=curCell.molcount;
        		numDensity+=MolcountToNumDensity(curCell.molcount);
        	}
    	}
		tmpV       /=pow(CellsCalcNum, 2);
		temp       /=pow(CellsCalcNum, 2);
		aver_molcount=int(aver_molcount/pow(CellsCalcNum, 2));
		numDensity/=pow(CellsCalcNum, 2);
        molPerCell[i][j][k]  += aver_molcount;
        meanV[i][j][k]       += tmpV;
        Temp[i][j][k]           += temp;
        pressure[i][j][k]    += numDensity * BOLTZMANN_K * temp; // pressure (in Pascal)
        aver_molcount=0;
        temp=0;
        numDensity=0;
        tmpV.x=0; tmpV.y=0; tmpV.z=0;
      }
    }
  }
  return 0;
}

int DSMC::ExportToVtk()
{
  ExportNum++;
  int i1, j1;
  const string fname = "./field_averaged_" + to_string(ExportNum) + ".vtk";
  ofstream vtkfile(fname);
  vtkfile<<"# vtk DataFile Version 2.0\n";
  vtkfile<<"DSMC Results\n";
  vtkfile<<"ASCII\n";
  vtkfile<<"DATASET STRUCTURED_POINTS\n";
  vtkfile<< "DIMENSIONS " << xcellcount << " " << ycellcount << " " << zcellcount << "\n";
  vtkfile<<"ORIGIN 0 0 0\n";
  vtkfile<<"SPACING 1 1 1\n";
  vtkfile << "POINT_DATA " << xcellcount * ycellcount * zcellcount << "\n";
  int k = 0;
  vtkfile << "\nSCALARS MolPerCell float\n";
  vtkfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j < ycellcount; j++)
  {
	  j1=int(j/CellsCalcNum);
	  for (int i = 0; i < xcellcount; i++)
	  {
		  i1=int(i/CellsCalcNum);
      	  vtkfile<<int(molPerCell[i1][j1][k] / VTKcount)<<" ";
	  }
      vtkfile << "\n";
  }
  vtkfile << "\nSCALARS Temperature_Translational float\n";
  vtkfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j < ycellcount; j++)
  {
		j1=int(j/CellsCalcNum);
		for (int i = 0; i < xcellcount; i++)
		{
			i1=int(i/CellsCalcNum);
			vtkfile<<Temp[i1][j1][k] / VTKcount<<" ";
		}
		vtkfile << "\n";
  }

  vtkfile << "\nSCALARS Vx float\n";
  vtkfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j < ycellcount; j++)
  {
	  j1=int(j/CellsCalcNum);
	  for (int i = 0; i < xcellcount; i++)
	  {
		  i1=int(i/CellsCalcNum);
		  vtkfile << meanV[i1][j1][k].x / VTKcount<<" ";
	  }
	  vtkfile << "\n";
  }

  vtkfile << "\nSCALARS Vy float\n";
  vtkfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j < ycellcount; j++)
  {
	  j1=int(j/CellsCalcNum);
	  for (int i = 0; i < xcellcount; i++)
	  {
		  i1=int(i/CellsCalcNum);
		  vtkfile << meanV[i1][j1][k].y / VTKcount<<" ";
	  }
    vtkfile << "\n";
  }
  vtkfile << "\nSCALARS Pressure float\n";
  vtkfile << "LOOKUP_TABLE default\n"; 
  for (int j = 0; j < ycellcount; j++)
  {
	  j1=int(j/CellsCalcNum);
	  for (int i = 0; i < xcellcount; i++)
	  {
		  i1=int(i/CellsCalcNum);
		  vtkfile << pressure[i1][j1][k] / VTKcount<<" ";
	  }
	vtkfile << "\n";
  }
  vtkfile.close();
  meanV.      assign(xcellcountforVtk, std::vector< std::vector<Vector3d> >(ycellcountforVtk, std::vector<Vector3d>(zcellcount, Vector3d{ 0, 0, 0 })));
  Temp.          assign(xcellcountforVtk, std::vector< std::vector<double> >(ycellcountforVtk, std::vector<double>(zcellcount, 0)));
  pressure.   assign(xcellcountforVtk, std::vector< std::vector<double> >(ycellcountforVtk, std::vector<double>(zcellcount, 0)));
  molPerCell. assign(xcellcountforVtk, std::vector< std::vector<int> >(ycellcountforVtk, std::vector<int>(zcellcount, 0)));

  VTKcount = 0;
  return 0;
}







