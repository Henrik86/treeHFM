/*
 * helpFunctions.cpp
 *
 *  Created on: 19.12.2012
 *      Author: henrikfailmezger
 */
#include <iostream>
#include <vector>
#include <math.h>
#include "helpFunctions.h"
using namespace std;


double checkPositive(std::vector<int> a){ //check for any positive element
    bool pos=0;
    for (int i=0;i<a.size();i++){
        if(a.at(i)>0){
            pos=1;
        }
    }
    return pos;
}

double sum(std::vector<double> a){
	double s=0.0;
	for (int i=0;i<a.size();i++){
		s=s+a.at(i);
	}
	return s;
}
std::vector<double> sum(std::vector<std::vector<double> > a,int dim){
	std::vector<double> c;
	if(dim==1){
		for (int i=0;i<a.size();i++){
			double s=0.0;
			for (int j=0;j<a.at(i).size();j++){
				s=s+a.at(i).at(j);
			}
			c.push_back(s);
		}
	}else{
		for (int i=0;i<a.at(0).size();i++){
					double s=0.0;
					for (int j=0;j<a.size();j++){
						s=s+a.at(j).at(i);
					}
					c.push_back(s);
		}
	}
	return c;
}
std::vector<double> zeros(int len){
    std::vector<double> a(len,0.0);
    //for (int i=0;i<len;i++){
    //	a.push_back(0.0);
    //}
    return a;
}

std::vector<int> zerosINT(int len){
    std::vector<int> a(len,0);
    //for (int i=0;i<len;i++){
    //	a.push_back(0.0);
    //}
    return a;
}
std::vector<std::vector<double> >  zeros(int row,int columns){
    std::vector<std::vector<double> > a(row, std::vector<double>(columns,0.0));
	//for (int i=0;i<len;i++){
	//	a.push_back(0.0);
	//}
	return a;
}
std::vector<double>  ones(int len){
	std::vector<double> a(len,1.0);
	//for (int i=0;i<len;i++){
	//	a.push_back(1.0);
	//}
	return a;
}
std::vector<int>  createSequence(int len,int number){
    std::vector<int> a(len,number);
    //for (int i=0;i<len;i++){
    //	a.push_back(1.0);
    //}
    return a;
}
std::vector<int>  createSequence(int len,int start,int end){
    int counter=start;
    std::vector<int> a;
    for (int i=0;i<len;i++){
    	a.push_back(counter++);
    }
    return a;
}
std::vector<double>  createSequenceDouble(int len,double number){
    std::vector<double> a(len,number);
    //for (int i=0;i<len;i++){
    //	a.push_back(1.0);
    //}
    return a;
}
std::vector<std::vector<double> > log2M(std::vector<std::vector<double> >m){ // counts how often a number is in the array
    for (int i = 0; i<m.size(); i++){
        for (int j = 0; j<m.at(0).size(); j++){
                m.at(i).at(j)=log2(m.at(i).at(j));
        }
    }
    return m;
}
std::vector<double>  log2M(std::vector<double> m){ // counts how often a number is in the array
    for (int i = 0; i<m.size(); i++){
            m.at(i)=log2(m.at(i));
    }
    return m;
}
std::vector<std::vector<double> >  repmat(std::vector<double> v, int rep,int dim){
	std::vector<std::vector<double> > a;
	if(dim==1){
		for (int i=0;i<rep;i++){
			a.push_back(v);
		}
    }else{
        for (int i=0;i<v.size();i++){
            std::vector<double> rowVector;
            for (int j=0;j<rep;j++){
                rowVector.push_back(v.at(i));
            }
            a.push_back(rowVector);
        }
    }
	return a;
}
std::vector<double>   repmat(double v, int rep,int dim){
	std::vector<double>  a;
	if(dim==1){
		for (int i=0;i<rep;i++){
			a.push_back(v);
		}
	}
	return a;
}
std::vector<int> setdiff(std::vector<int> a,int v){ //difference between vector, value
	std::vector<int> f;
	for (int i = 0; i<a.size(); i++){
		if(a.at(i) !=v){
			f.push_back(a.at(i));
		}
	}
	return f;
}

std::vector<int> findIndex(std::vector<std::vector<int> > a,int dim,int v){
	std::vector<int> f;
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
               if(a[i][dim]==v){
		    	  f.push_back(i);
		      }
	}

	return f;
}
std::vector<double> findIndex(std::vector<double> a,int v){
    std::vector<double> f;
    for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
        if(a[i]==v){
            f.push_back(i);
        }
    }
    
    return f;
}
int countIndex(std::vector<std::vector<int> > a,int dim,int v){ // counts how often a number is in the array
	int c=0;
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
				if(a[i][dim]==v){
		    	  c++;
		      }
	}
	return c;
}
int countIndex(std::vector<int> a,int v){ // counts how often a number is in the array
    int c=0;
    for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
        if(a[i]==v){
            c++;
        }
    }
    return c;
}
std::vector<double> getValuesIndex(std::vector<double> a,std::vector<int> index){ // counts how often a number is in the array
    std::vector<double> f;
    for (int i = 0; i<index.size(); i++){      // durchl�uft alle Werte von 0 bis 5
        f.push_back(a.at(index.at(i)));
    }
    return f;
}
std::vector<std::vector<double> > logicalMatrix(std::vector<std::vector<double> >a,double b){ // counts how often a number is in the array
    std::vector<std::vector<double> > m(a.size(), std::vector<double>(a.at(0).size()));
    for (int i = 0; i<a.size(); i++){
        std::vector<double> c;
        for (int j = 0; j<a.at(0).size(); j++){
            if(a.at(i).at(j)==b){// durchl�uft alle Werte von 0 bis 5
                m.at(i).at(j)=1;
            }else{
                m.at(i).at(j)=0;
            }
        }
    }
    return m;
}
std::vector<double>  logicalMatrix(std::vector<double> a,double b){ // counts how often a number is in the array
    std::vector<double> m(a.size());
    for (int i = 0; i<a.size(); i++){
        if(a.at(i)==b){// durchl�uft alle Werte von 0 bis 5
                m.at(i)=1.0;
        }else{
                m.at(i)=0.0;
        }
        
    }
    return m;
}
std::vector<double> arrangeItems(std::vector<double>a,std::vector<int>b){ // counts how often a number is in the array
	std::vector<double> c;
	for (int i = 0; i<b.size(); i++){      // durchl�uft alle Werte von 0 bis 5
				c.push_back(a.at(b.at(i)));
	}
	return c;
}
std::vector<int> arrangeItems(std::vector<int>a,std::vector<int>b){ // counts how often a number is in the array
    std::vector<int> c;
    for (int i = 0; i<b.size(); i++){      // durchl�uft alle Werte von 0 bis 5
        c.push_back(a.at(b.at(i)));
    }
    return c;
}
std::vector<std::vector<double> > arrangeItems(std::vector<std::vector<double> >a,std::vector<int>b){ // counts how often a number is in the array
	std::vector<std::vector<double> > d;
	for (int i = 0; i<a.size(); i++){
		std::vector<double> c;
		for (int j = 0; j<b.size(); j++){      // durchl�uft alle Werte von 0 bis 5
				c.push_back(a.at(i).at(b.at(j)));
		}
		d.push_back(c);
	}
	return d;
}
//
std::vector< std::vector<double> > transpose( std::vector< std::vector<double> > a){
	
    std::vector<double> v(a.size());
	std::vector<vector<double> > b(a.at(0).size(),v);
	
	for (int i = 0; i < a.size(); i++){
		for (int j = 0; j < a.at(0).size(); j++){
			b.at(j).at(i) = a.at(i).at(j);
		}
	}
	return b;
}
std::vector< std::vector<double> > transpose( std::vector<double>  a){ //row to column vector
    
    std::vector< std::vector<double> > v(a.size(),std::vector<double>(1));
    
    for (int i = 0; i < a.size(); i++){
            v.at(i).at(0) = a.at(i);
    }
    return v;
}
//
std::vector<double> multEl(std::vector<double> a,std::vector<double> b){ // element-by-element product
	if(a.size()==b.size()){
		for (int i = 0; i<a.size(); i++){
			a.at(i)=a.at(i)*b.at(i);
		}
	}
	return a;
}
std::vector<std::vector<double> > multEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b){ // element-by-element product
	if(a.size()==b.size()){
		for (int i = 0; i<a.size(); i++){
			for (int j = 0; j<a.at(i).size(); j++){
				a.at(i).at(j)=a.at(i).at(j)*b.at(i).at(j);
			}
		}
	}
	return a;
}
std::vector<std::vector<std::vector<double> > > multEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > > b){ // element-by-element product
    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                for (int k = 0; k<a.at(i).at(j).size(); k++){
                    a.at(i).at(j).at(k)=a.at(i).at(j).at(k)*b.at(i).at(j).at(k);
                }
            }
        }
    }
    return a;
}
std::vector<std::vector<std::vector<double> > > multEl(std::vector<std::vector<std::vector<double> > >a,double  b){ // element-by-element product
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                for (int k = 0; k<a.at(i).at(j).size(); k++){
                    a.at(i).at(j).at(k)=a.at(i).at(j).at(k)*b;
                }
            }
        }
    
    return a;
}
std::vector<double> multEl(std::vector<double> a,double b){ // element-by-element product
		for (int i = 0; i<a.size(); i++){
			a.at(i)=a.at(i)*b;
		}

	return a;
}
std::vector<double> divEl(std::vector<double> a,double b){ // element-by-element product
    for (int i = 0; i<a.size(); i++){
        a.at(i)=a.at(i)/b;
    }
    
    return a;
}
std::vector<std::vector<std::vector<double> > > divEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > > b){ // element-by-element product
    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                for (int k = 0; k<a.at(i).at(j).size(); k++){
                    a.at(i).at(j).at(k)=a.at(i).at(j).at(k)/b.at(i).at(j).at(k);
                }
            }
        }
    }
    return a;
}
std::vector<std::vector<double> > divEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b){ // element-by-element product
    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                a.at(i).at(j)=a.at(i).at(j)/b.at(i).at(j);
            }
        }
    }
    return a;
}
std::vector<std::vector<double> > divEl(std::vector<std::vector<double> > a,double b){ // element-by-element product
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                a.at(i).at(j)=a.at(i).at(j)/b;
            }
        }
    
    return a;
}
std::vector<std::vector<double> > matrixMultiplication(std::vector<double> a,std::vector<double> b){ // element-by-element product
    std::vector<std::vector<double> > result(a.size(), std::vector<double>(b.size()));
    
    int n=a.size();
    int p=b.size();
    
    for(int i=0; i<n; i++) {
        for(int j=0; j<p; j++) {
            result.at(i).at(j)=a.at(i)*b.at(j);
        }
    }

    return result;
}
std::vector<std::vector<double> > matrixMultiplication(std::vector<std::vector<double> > a,std::vector<std::vector<double> >b){ // element-by-element product
    
    
    int n =a.size();
    int m = a[0].size();
    int p = b[0].size();
    //
    std::vector<std::vector<double> > c(n, std::vector<double>(p));
    //
    for(int i=0; i<n; i++) {
        for(int j=0; j<p; j++) {
            double sum =0;
            for (int k=0;k<m;k++){
                sum=sum+a.at(i).at(k)*b.at(k).at(j);
            }
            c.at(i).at(j)=sum;
        }
    }
    
    return c;
}



std::vector<double> addEl(std::vector<double> a,std::vector<double> b){ // element-by-element product
	if(a.size()==b.size()){
		for (int i = 0; i<a.size(); i++){
			a.at(i)=a.at(i)+b.at(i);
		}
	}
	return a;
}
std::vector<std::vector<double> > addEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b){ // element-by-element product
    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                a.at(i).at(j)=a.at(i).at(j)+b.at(i).at(j);
            }
        }
    }
    return a;
}
std::vector<std::vector<std::vector<double> > > addEl(std::vector<std::vector<std::vector<double> > >a,double c){ // element-by-element product
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                for (int k = 0; k<a.at(i).size(); k++){
                    a.at(i).at(j).at(k)=a.at(i).at(j).at(k)+c;
                }
            }
        }
    
    return a;
}
std::vector<std::vector<std::vector<double> > > addEl(std::vector<std::vector<std::vector<double> > >a,std::vector<std::vector<std::vector<double> > > b){ // element-by-element product
    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                for (int k = 0; k<a.at(i).at(j).size(); k++){
                    a.at(i).at(j).at(k)=a.at(i).at(j).at(k)+b.at(i).at(j).at(k);
                }
            }
        }
    }
    return a;
}
std::vector<int> subEl(std::vector<int> a,int b){ // subtract element from vector
	for (int i = 0; i<a.size(); i++){
			a.at(i)=a.at(i)-b;
	}
	return a;
}
std::vector<std::vector<int> > subEl(std::vector<std::vector<int> > a , int b){ // subtract element from vector
	for (int i = 0; i<a.size(); i++){
		for (int j = 0; j<a.at(i).size(); j++){
			a.at(i).at(j)=a.at(i).at(j)-b;
		}
	}
	return a;
}
std::vector<std::vector<double> > subEl(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b){ // element-by-element product

    if(a.size()==b.size()){
        for (int i = 0; i<a.size(); i++){
            for (int j = 0; j<a.at(i).size(); j++){
                a.at(i).at(j)=a.at(i).at(j)-b.at(i).at(j);
            }
        }
    }
    return a;
}

//
int intersect(std::vector<int> a, std::vector<int> b){ //boolean, 1 if any element interset
	int intersect =0;
	for(int i=0;i<a.size(); i++){
		for(int j=0;j<b.size(); j++){
			if(a.at(i)==b.at(j)){
				intersect =1;
				return intersect;
			}
		}
	}
	return intersect;
}
int intersect(std::vector<int> a, int b){ //boolean, 1 if any element interset
	int intersect =0;
	for(int i=0;i<a.size(); i++){
			if(a.at(i)==b){
				intersect =1;
				return intersect;
			}
		}
	return intersect;
}
//
std::vector<std::vector<double> >replaceColumn(std::vector<std::vector<double> > a,std::vector<double> b, int c){
    for(int i=0;i<a.size(); i++){
            a.at(i).at(c)=b.at(i);
    }
    return a;
}
std::vector<double> getColumn(std::vector<std::vector<double> > a, int c){
    std::vector<double> b; //puts everything below each other, jumps over empty entries
    for(int i=0;i<a.size(); i++){
        if(a.at(i).size()>c){
            b.push_back(a.at(i).at(c));
        }
    
    }
    return b;
}
std::vector<int> getColumn(std::vector<std::vector<int> > a, int c){
    std::vector<int> b; //puts everything below each other, jumps over empty entries
    for(int i=0;i<a.size(); i++){
        if(a.at(i).size()>=c){
            b.push_back(a.at(i).at(c));
        }
        
    }
    return b;
}
Max max(std::vector<double> a){
	int m=-INFINITY;
	int mIndex=0;
	Max mv;
	for (int i=0; i<a.size();i++){
		if(a[i]>m){
			m=a[i];
			mIndex=i;
		}
	}
	mv.value=m;
	mv.index=mIndex;
	return mv;
}
Max max(std::vector<int> a){
        int m=-INFINITY;
		int mIndex=0;
		Max mv;
		for (int i=0; i<a.size();i++){
			if(a[i]>m){
				m=a[i];
				mIndex=i;
			}
		}
		mv.value=m;
		mv.index=mIndex;
		return mv;
}
std::vector<double> matrixToVector(std::vector<std::vector<double> > a){
    std::vector<double> b;
    for(int i=0;i<a.size(); i++){
        for(int j=0;j<a.at(0).size(); j++){
            b.push_back(a.at(i).at(j));
        }
    }
    return b;
}

MaxMatrix maxMatrix(std::vector<std::vector<double> > a,int dim){ //finds maximum row/column values in a matrix
        //cout<<"maxMatrix"<<endl;
        std::vector<double> mValues;
        std::vector<int> mIndices;
        double m=-INFINITY;
		int mIndex=0;
		MaxMatrix mv;
        if(dim==1)
            for (int i=0; i<a.size();i++){
                m=-INFINITY;
                mIndex=0;
                for (int j=0; j<a.at(0).size();j++){
                    if(a.at(i).at(j)>m){
                        m=a.at(i).at(j);
                        mIndex=j;
                    }
                }
                mValues.push_back(m);
                mIndices.push_back(mIndex);
            }
        else{
            for (int i=0; i<a.at(0).size();i++){
                m=-INFINITY;
                mIndex=0;
                for (int j=0; j<a.size();j++){
                   
                    if(a.at(j).at(i)>m){
                        m=a.at(j).at(i);
                        mIndex=j;
                        
                    }
                }
                
                mValues.push_back(m);
                mIndices.push_back(mIndex);
            }
        }
		mv.value=mValues;
		mv.index=mIndices;
		return mv;
}
 

Min min(std::vector<int> a){
	int m=a[0];
	int mIndex=0;
	for (int i=0; i<a.size();i++){
		if(a[i]<=m){
			m=a[i];
			mIndex=i;
		}
	}
	Min mv;
	mv.index=mIndex;
	mv.value=m;
	return mv;
}
std::vector<int> add(std::vector<int> a,int v){
	std::vector<int> s;
	for (int i = 0; i < a.size(); i++){
		s.push_back(a[i]+v);
	}
	return s;
}
std::vector<std::vector<double> > initialise(int n){
	std::vector<std::vector<double> > s;
	for (int i = 0; i < n; i++){
		std::vector<double> b;
		s.push_back(b);
	}
	return s;
}
std::vector<double> tabulate(std::vector<int> a){
    //
    Max maxValue;
    maxValue=max(a);
    std::vector<double> tVector=zeros(maxValue.value+1);
    //
    for(int i=0;i<a.size(); i++){
        tVector.at(a.at(i))=tVector.at(a.at(i))+1;
    }
    return tVector;
}
//
double* vectorToArray1D(std::vector<double> a){
    double* temp;
    temp = new double[a.size()];
    for(unsigned i=0; (i < a.size()); i++)
    {
        temp[i] = a[i];
    }
    return temp;
}
//
double** vectorToArray(std::vector<std::vector<double> > a){
    double** temp;
    temp = new double*[a.size()];
    for(unsigned i=0; (i < a.size()); i++)
    {
        temp[i] = new double[a.at(0).size()];
        for(unsigned j=0; (j < a.at(0).size()); j++)
        {
            temp[i][j] = a[i][j];
        }
    }
    return temp;
}
//
double** vectorToArray(std::vector<double> a){
    double** temp;
    temp = new double*[a.size()];
    for(unsigned i=0; (i < a.size()); i++)
    {       temp[i] = new double[1];
            temp[i][0] = a[i];
    }
    return temp;
}
void print(int v){
	for (int i = 0; i<v; i++){      // durchl�uft alle Werte von 0 bis 5
		      cout << i <<"," << endl;
		}
}

void printArray(int a[],int lengthA) {
	for (int i = 0; i<lengthA; i++){      // durchl�uft alle Werte von 0 bis 5
	      cout << a[i] <<",";
	}
}
void printArray(double a[],int lengthA) {
	for (int i = 0; i<lengthA; i++){      // durchl�uft alle Werte von 0 bis 5
	      cout << a[i] <<",";
	}
}
void printVector(std::vector<int> a) {
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
	      cout << a[i] <<",";
	}
	cout<<endl;
}
void printVector(std::vector<std::vector<double> > a) {
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
		if(a.at(i).size()>0){
			for (int j = 0; j<a.at(i).size(); j++){
				cout << a.at(i).at(j) <<",";
			}
			cout<<endl;
		}
	}
	cout<<endl;
}
void printVector(std::vector<std::vector<int> > a) {
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
		for (int j = 0; j<a[i].size(); j++){
			cout << a[i][j] <<",";
		}
		cout<<endl;
	}
	cout<<endl;
}
void printVector(std::vector<double> a) {
	for (int i = 0; i<a.size(); i++){      // durchl�uft alle Werte von 0 bis 5
	      cout << a[i] <<",";
	}
	cout<<endl;
}




