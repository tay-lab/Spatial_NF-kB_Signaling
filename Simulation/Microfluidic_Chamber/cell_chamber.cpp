/************************************
* Class implementation file. 		*
* Functions used in the Cell class. *
************************************/

# include "cell_chamber.h"	// The header file containing the Cell class.
using namespace std;

// Constructor, defines what should be initiated when the Cell class is called.
Cell::Cell(){
	y = new std::vector<double>(Nvar,0.0);		// Nvar doubles with value 0.0.
	dy = new std::vector<double>(Nvar,0.0);		// The derivatives.	
	ynew = new std::vector<double>(Nvar,0.0);	// Used when performing time step.
};

// Destructor:
Cell::~Cell(){
	delete y;
	delete dy;
};


// Functions:
void Cell::AssignParameters ()
{
  // Units hours and muM
  kN = 11.3; B = 1.09; delta = 0.029; Ki = 0.035;	// NFkB
  p = 0.0; ta = 59.5; gma = 2.; Kn1 = 0.5;			// Ima
  Aa = 63.; ga = 0.4; tla = 14.4;					// Ia

  tA = 5.; gAm = 1.;								// A20m
  tlA = 15.; gA = 0.1;								// A20
  Ktot = 1.; mu = 10; beta = 2.; sigma = 0.25; KTNF = 1.; scale = 1; // IKK/IKKi
  T_prod = 0; Hill = 3;								// TNF
};

//******************************************************
void Cell::Initialize(double initNF, double initIm, double initIa, double initA20m, double initA20, double initIKK, double initIKKi){
  (*y)[0] = initNF; 
  (*y)[1] = initIm; 
  (*y)[2] = initIa;
  (*y)[3] = initA20m;
  (*y)[4] = initA20;
  (*y)[5] = initIKK;
  (*y)[6] = initIKKi;
}

void Cell::InitializeT (double initvalues) {
  (*y)[7] = initvalues; 
};

void Cell::InitializeGenes (double initvalues) {
  (*y)[8] = initvalues;  
  (*y)[9] = initvalues; 
  (*y)[10] = initvalues; 
};

//**************************************************
// Calculating derivatives:
//**************************************************
void Cell::Der () {
  
  double dummy = (1 - (*y).at(0)) / (Ki + ((*y).at(2)));
  // NFkB
  (*dy).at(0) = kN*dummy - B * (*y).at(2) *(*y).at(0)/((*y).at(0)+delta);
  // Ima
  (*dy).at(1) = p + ta*sqr((*y).at(0))/(Kn1+sqr((*y).at(0))) - gma*(*y).at(1);
  // Ia
  (*dy).at(2) = tla*(*y).at(1) - Aa*(*y).at(5)*dummy*(*y).at(2) - ga*(*y).at(2);
  // A20m
  (*dy).at(3) = tA*sqr((*y).at(0)) - gAm*(*y).at(3);
  // A20
  (*dy).at(4) = tlA*(*y).at(3) - gA*(*y).at(4);
  // IKK
  //(*dy).at(5) = KTNF/10.*(*y).at(7)*(Ktot-(*y).at(5)-(*y).at(6)) - mu*sqr((*y).at(5));
  (*dy).at(5) = scale*power( (*y).at(7), Hill)/( power( (*y).at(7), Hill)  + power(KTNF, Hill) )*(Ktot-(*y).at(5)-(*y).at(6)) - mu*sqr((*y).at(5));
  // IKKi
  (*dy).at(6) = mu*sqr((*y).at(5)) - beta*(*y).at(6)/(sigma*sqr((*y).at(4))+1);
  // TNF
  (*dy).at(7) = T_prod;
  //Gene1:
  (*dy).at(8) = 4.*(*y).at(0) - (*y).at(9)/(1./4.);
  //Gene2:
  (*dy).at(9) = 2.*(*y).at(0) - (*y).at(10)/2.;
  //Gene3:
  (*dy).at(10) = 1/2.*(*y).at(0) - (*y).at(11)/64.;
};

//******************************************************
// The following function calculates the new concentration values after a single time step
void Cell::Timestep(double dt){
// An implementation of the Runga-Kutta method for solving differential equations:

	for(int i = 0;i < Nvar; i++){
		yval[i] = (*y).at(i);}

	Der();	// Calculate derivatives.
	for(int i = 0; i < Nvar; i++){
		k1[i] = dt * (*dy).at(i);}	// Replace the y-values with the modified form.
	for (int i = 0; i < Nvar; i++){
		(*y).at(i) = (*y).at(i) + 1./2.*k1[i];}

	// Redo for step 2, 3 and 4:
	Der();	// Calculate derivatives.
	for(int i = 0; i < Nvar; i++){
		k2[i] = dt * (*dy).at(i);	// Replace the y-values with the modified form.
		(*y).at(i) = yval[i];}		// Reset to original form before calculating new value.
	for (int i = 0; i < Nvar; i++){
		(*y).at(i) = (*y).at(i) + 1./2.*k2[i];
	}

	Der();	// Calculate derivatives.
	for(int i = 0; i < Nvar; i++){
		k3[i] = dt * (*dy).at(i);	// Replace the y-values with the modified form.
		(*y).at(i) = yval[i];}		// Reset to original form before calculating new value.
	for (int i = 0; i < Nvar; i++){
		(*y).at(i) = (*y).at(i) + k3[i];
	}

	Der();	// Calculate derivatives.
	for(int i = 0; i < Nvar; i++){
		k4[i] = dt * (*dy).at(i);	// Replace the y-values with the modified form.
		(*y).at(i) = yval[i];}		// Reset to original form before calculating new value.

	for (int i=0; i<Nvar; i++){		// Calculate the final function value.
		(*ynew).at(i) = (*y).at(i) + 1./6.*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);
		if((*ynew).at(i)<0.0){
			(*ynew).at(i) = 0;};
	};

	for (int i=0;i<Nvar;i++){
		y = ynew;	// Replace old y with the newly calculated.
	}
}


// Current concentration of variables:
double Cell::NOut(void){
	return((*y).at(0));}
double Cell::TOut(){
	return((*y).at(7));}
double Cell::G1Out(){
	return((*y).at(8));}
double Cell::G2Out(){
	return((*y).at(9));}
double Cell::G3Out(){
	return((*y).at(10));}	
	
	
void Cell::ReassignT_prod(double newT_prod){
	T_prod = newT_prod;
}
void Cell::ReassignT(double newT){
	(*y).at(7) = newT;}
	
void Cell::ReAssignKTNF (double NewKTNF){
  //Reassigns KTNF:
  KTNF = NewKTNF;
};

//**************************************************
// Small functions:
//**************************************************
double Cell::sqr(double a) {
  return a*a;
  };

double Cell::power(double a, int Hill) {
  double powera = 1;
  for (int i = 0; i < Hill; i++) powera = powera*a;
  return powera;
  };
	