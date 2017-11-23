

#ifndef HELPFUNCTIONS_H_
#define HELPFUNCTIONS_H_

#include <iostream>
#include <vector>
#include <math.h> 
using namespace std;
struct Max {
   //
   int value;
   int index;
};
struct Min {
   //
   int value;
   int index;
};
struct MaxMatrix {
    //
    std::vector<double> value;
    std::vector<int>  index;
};

//
double checkPositive(std::vector<int> a);
//
std::vector<int> unionV(std::vector<int>A, std::vector<int>B);
//
MaxMatrix maxMatrix(std::vector<std::vector<double> > a,int dim);
//
std::vector<std::vector<double> > log2M(std::vector<std::vector<double> >m);
std::vector<double> log2M(std::vector<double> m);
//
std::vector<std::vector<double> > logicalMatrix(std::vector<std::vector<double> >a,double b);
std::vector<double>  logicalMatrix(std::vector<double> a,double b);
//
std::vector<double> multEl(std::vector<double> a,std::vector<double> b);
std::vector<std::vector<double> > multEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b);
std::vector<std::vector<std::vector<double> > > multEl(std::vector<std::vector<std::vector<double> > >a,double  b);
std::vector<double> multEl(std::vector<double> a,double b);
std::vector<std::vector<std::vector<double> > > addEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > >b);

//
std::vector<double> divEl(std::vector<double> a,double b);
std::vector<std::vector<std::vector<double> > > divEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > > b);
std::vector<std::vector<double> > divEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b);
std::vector<std::vector<double> > divEl(std::vector<std::vector<double> > a,double b);
//
std::vector<double> addEl(std::vector<double> a, std::vector<double> b);
std::vector<std::vector<double> > addEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b);
std::vector<std::vector<std::vector<double> > > addEl(std::vector<std::vector<std::vector<double> > >a,double c);
std::vector<std::vector<std::vector<double> > > addEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > >b);
//
std::vector<int> subEl(std::vector<int> a,int b);
std::vector<std::vector<int> > subEl(std::vector<std::vector<int> > a , int b);
std::vector<std::vector<double> > subEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b);
//
double sum(std::vector<double> a);
std::vector<double> sum(std::vector<std::vector<double> > a,int dim);
//
std::vector<std::vector<double> > matrixMultiplication(std::vector<double> a,std::vector<double> b);
std::vector<std::vector<double> > matrixMultiplication(std::vector<std::vector<double> > a,std::vector<std::vector<double> >b);
//
std::vector<double> zeros(int len);
std::vector<int> zerosINT(int len);
std::vector<std::vector<double> >  zeros(int row,int columns);
std::vector<double>  ones(int len);
std::vector<int>  createSequence(int len,int number);
std::vector<double>  createSequenceDouble(int len,double number);
std::vector<int>  createSequence(int len,int start,int end);
//
std::vector< std::vector<double> > transpose( std::vector< std::vector<double> > a);
std::vector< std::vector<double> > transpose( std::vector<double>  a);
//
std::vector<double> arrangeItems(std::vector<double>a,std::vector<int>b);
std::vector<std::vector<double> > arrangeItems(std::vector<std::vector<double> >a,std::vector<int>b);
std::vector<int> arrangeItems(std::vector<int>a,std::vector<int>b);
//
std::vector<int>  setdiff(std::vector<int> a,int v);
int intersect(std::vector<int> a,std::vector<int> b);
int intersect(std::vector<int> a,int b);
std::vector<double> tabulate(std::vector<int> a);
//
std::vector<std::vector<double> >  repmat(std::vector<double> v, int rep,int dim);
std::vector<double>  repmat(double v, int rep,int dim);
//
std::vector<double> matrixToVector(std::vector<std::vector<double> > a);
std::vector<double> getColumn(std::vector<std::vector<double> > a, int c);
std::vector<int> getColumn(std::vector<std::vector<int> > a, int c);
std::vector<std::vector<double> >replaceColumn(std::vector<std::vector<double> > a,std::vector<double> b, int c);
//
std::vector<int> findIndex(std::vector<std::vector<int> > a,int dim,int v);
std::vector<double> findIndex(std::vector<double> a,int v);
int countIndex(std::vector<std::vector<int> > a,int dim,int v);
int countIndex(std::vector<int> a,int v);
std::vector<double> getValuesIndex(std::vector<double> a,std::vector<int> index);

std::vector<int> arrayToVector(int a[]);
std::vector<double> arrayToVector(double a[]);
//
double** vectorToArray(std::vector<std::vector<double> > a);
double* vectorToArray1D(std::vector<double> a);
double** vectorToArray(std::vector<double> a);
//
Max max(std::vector<double> a);
Max max(std::vector<int> a);
Min min(std::vector<int> a);
std::vector<int> add(std::vector<int> a,int v);
std::vector<std::vector<double> > initialise(int n);
//
void print(int v);
void printArray(int a[],int lengthA);
void printArray(double a[],int lengthA);
void printVector(std::vector<int> a);
void printVector(std::vector<std::vector<double> > a);
void printVector(std::vector<std::vector<int> > a);
void printVector(std::vector<double> a);
void printFactorGraph(struct factorGraph fg);
//

#endif /* HELPFUNCTIONS_H_ */
