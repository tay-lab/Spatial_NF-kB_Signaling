/**********************************
* Class header file:              *
* Definition of the "Cell" class  *
**********************************/
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define Nvar 12
#define Xcells 16  // Must be even for boundary conditions to be correct.
#define Ycells 130 // Must be even for boundary conditions to be correct.

class Cell
{
 private:
  double kN, B, delta, Ki; 
  double ta, gma, Kn1, p;
  double tla, Aa, ga, tA, gAm, tlA, gA;
  double Ktot, mu, beta, sigma, KTNF, T_prod, scale; 
  int Hill;
  
  double k1[Nvar], k2[Nvar], k3[Nvar], k4[Nvar];
  double yval[Nvar];

  // Pointers til vectorer, der initialiseres i constructor:
  std::vector<double> *y;
  std::vector<double> *dy;
  std::vector<double> *ynew; 

  double sqr(double a);
  double power(double a, int Hill);

public:
  //constructor
  Cell();
  //destructor
  ~Cell();
  //method (functions):
  void AssignParameters ();
  void Initialize (double initNF, double initIm, double initIa, double initA20m, double initA20, double initIKK, double initIKKi);
  void InitializeT (double initvalues);
  void InitializeGenes (double initvalues);
  void Der ();
  void Timestep (double dt);
  void ReassignT (double NewTNF);
  void AddTNFout (double Add);
  double GetTNFoutValue ();
  void ReAssignIM (double NewIM);
  void ReAssignN (double NewN);
  void ReAssignKTNF (double NewKTNF);
  void ReassignT_prod (double NewT_prod);
  
  double NOut();
  double TOut();
  double G1Out();
  double G2Out();
  double G3Out();
};
